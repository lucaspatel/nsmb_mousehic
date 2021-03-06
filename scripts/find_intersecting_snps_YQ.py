import sys, pysam, gzip, pdb, argparse, array
#from pympler import asizeof

#Edited by DG. Changes:
#       -Now handles soft clipping. [see lines 333-335]
#       -Throws out reads where allele present at a given SNP does not match reference or alternative [see line 351]
#       -Counts totals of reads eliminated and prints out to stderr
#       -Includes matched sequence legnth in output read name in remap.fq so this can be checked in next step of pipeline [see lines 247-252]
#       -sequence in new fastq matches direction and complement of sequence in original fastq   [see lines 255-259]

dgcounter=0             #DG EDIT HERE   [counter for reads with "0 potential seqs", which are also eliminated. These correspond to reads that have non-M/S in the cigar, have a non-ref/alt allele at the SNP position, overlap > 10 SNPs, or other potential reasons]
dgcounter2=0            #DG EDIT HERE   [counter for reads that overlap >10 SNPs, which are give "0 potential seqs" and eliminated]
dgcounter3=0            #DG EDIT HERE   [counter for reads with >10 potential seqs (overlapping 4-10 SNPs), which are eliminated]
dgcounter4=0            #DG EDIT HERE   [counter for cases in which the allele present at a given SNP position in a read does not match reference or alternative]

#### Class to hold the data for a single SNP
class SNP:
    def __init__(self,snp_line):
        snp_split=snp_line.strip().split()
        self.pos=int(snp_split[0])-1
        self.alleles=[snp_split[1],snp_split[2]]
        self.ptype="snp"
        self.max_len=0
        for i in range(len(self.alleles)):
            if self.alleles[i]=="-":
                self.alleles[i]=""
                self.ptype="indel"
            elif len(self.alleles[i])>self.max_len:
                self.max_len=len(self.alleles[i])
        if self.max_len>1:
            self.ptype="indel"

        haplo_chars=snp_split[5:]
        haplos=array.array('B',[0])*len(haplo_chars)

    def add_allele(self,new_alleles):
        for new_allele in new_alleles:
            if new_allele=="-":
                self.ptype="indel"
                new_allele=""
            if not (new_allele in self.alleles):
                self.alleles.append(new_allele)
                if len(new_allele)>self.max_len:
                    self.max_len=len(new_allele)
        if self.max_len>1:
            self.ptype="indel"

    def shift_indel(self):
        self.pos+=1
        self.max_len-=1
        i=0
        while i < len(self.alleles):
            if len(self.alleles)<=1:
                self.alleles.pop(i)
            else:
                self.alleles[i]=self.alleles[i][1:]
                i+=1
        self.alleles.append("")

#### Class to keep track of all the information read in from the bamfile/snpfile
class Bam_scanner:
    # Constructor: opens files, creates initial table
    def __init__(self,is_paired_end,max_window,file_name,keep_file_name,remap_name,remap_num_name,fastq_names,snp_dir):
        self.is_paired_end=is_paired_end

        ### Read in all input files and create output files
        self.snp_dir=snp_dir
        sys.stderr.write("{}\n".format(file_name))
        self.bamfile=pysam.Samfile(file_name, "rb")
        self.keep_bam=pysam.Samfile(keep_file_name, "wb", template=self.bamfile)
        self.remap_bam=pysam.Samfile(remap_name, "wb", template=self.bamfile)
        self.remap_num_file=gzip.open(remap_num_name,"w")
        self.fastqs=[gzip.open(fqn,"w") for fqn in fastq_names]
        try:
            self.cur_read=self.bamfile.next()
        except:
            sys.stderr.write("No lines available for input")
            return()
        self.end_of_file=False

        self.remap_num=1
        self.ref_match=0
        self.alt_match=0
        self.no_match=0
        self.toss=0
        self.nosnp=0
        self.remap=0
        self.tot=0
        self.printstats = True

        self.num_reads=0

        ### Initialize the indel tracking dictionary

        self.pos=self.cur_read.pos
        self.chr_num=self.cur_read.tid
        self.chr_name=self.bamfile.getrname(self.cur_read.tid)
        self.max_window=max_window

        self.num_reads=0
        ### Initialize the read tracking tables
        self.read_table=[[] for x in range(self.max_window)]

        ### Initialize the SNP and indel tracking tables
        self.switch_chr()
        self.fill_table()


    # fills the table of reads starting from the current position and extending for the next <max_window> base pairs
    def fill_table(self):
        if self.end_of_file:
            return()
        if self.num_reads==0:
            self.pos=self.cur_read.pos
            self.init_snp_table()
            #self.num_reads+=1000
        while self.cur_read.pos<self.pos+self.max_window:
            self.num_reads+=1
            self.read_table[self.cur_read.pos % self.max_window].append(self.cur_read)
            try:
                self.cur_read=self.bamfile.next()
            except:
                self.empty_table()
                self.end_of_file=True
                return()
            if self.cur_read.tid != self.chr_num:
                self.empty_table()
                self.chr_num=self.cur_read.tid
                try:
                    self.chr_name=self.bamfile.getrname(self.chr_num)
                except:
                    sys.stderr.write("Problem with tid: "+str(self.chr_num))
                    self.skip_chr()
                self.pos=self.cur_read.pos
                self.switch_chr()
                self.fill_table()

    # Switches to looking for SNPs on the next chromosome
    def switch_chr(self):
        chr_match=False
        while not chr_match and not self.end_of_file:
            try:
                #print "%s/%s.snps.txt.gz"%(self.snp_dir,self.chr_name)
                self.snpfile = gzip.open("%s/%s.snps.txt.gz"%(self.snp_dir,self.chr_name))
                sys.stderr.write("Starting on chromosome "+self.chr_name+"\n")
                chr_match=True
            except:
                sys.stderr.write("SNP file for chromosome "+self.chr_name+" is not found. Skipping these reads.\n")
                self.skip_chr()

        self.end_of_snp_file=False
        self.get_next_snp()

    # Initializes the SNP table
    def init_snp_table(self):
        # create an empty SNP table
        self.num_snps=0
        self.indel_dict={}
        self.snp_table=[0 for x in range(self.max_window)]
        self.indel_table=[[] for x in range(self.max_window)]
        # skips SNPs that are upstream of the current read
        while not self.end_of_snp_file and self.cur_snp.pos<self.pos:
            self.get_next_snp()

        # adds SNPs downstream of the current read and within the current window
        while not self.end_of_snp_file and self.cur_snp.pos<self.pos+self.max_window:
            if self.cur_snp.ptype=="snp":
                self.add_snp()
            else:
                self.add_indel()
            self.get_next_snp()

        #sys.stderr.write(str(self.num_snps)+"\n")

    def add_snp(self):
        cur_pos=self.cur_snp.pos % self.max_window
        if self.snp_table[cur_pos]==0:
            self.num_snps+=1
            self.snp_table[cur_pos]=self.cur_snp
        elif isinstance(self.snp_table[cur_pos],SNP):
            self.snp_table[cur_pos].add_allele(self.cur_snp.alleles)

    def add_indel(self):
        position=self.cur_snp.pos
        if self.indel_dict.has_key(position):
            start=self.indel_dict[position].max_len
            self.indel_dict[position].add_allele(self.cur_snp.alleles)
        else:
            self.indel_dict[position]=self.cur_snp
            start=0
        end=self.indel_dict[position].max_len
        i=start
        while i<end and self.cur_snp.pos+i<self.pos+self.max_window:
            self.indel_table[(self.cur_snp.pos+i)%self.max_window].append(position)
            i+=1

    #to read in next SNP or signals end of file
    def get_next_snp(self):
        snp_line=self.snpfile.readline()
        if snp_line:
            self.cur_snp=SNP(snp_line)
        else:
            self.end_of_snp_file=True

    # Skips all of the reads coming from this chromosome and moves on to the next
    # Used if the SNP file can't be located
    def skip_chr(self):
        while self.cur_read.tid == self.chr_num:
            try:
                self.cur_read=self.bamfile.next()
            except:
                self.empty_table()
                self.end_of_file=True
                return()

        self.chr_num=self.cur_read.tid
        try:
            self.chr_name=self.bamfile.getrname(self.chr_num)
        except:
            sys.stderr.write("Problem with tid: "+str(self.chr_num))
            self.skip_chr()

    # Processes all reads that map to the current position and removes them from the read table
    # Treats reads as single-end
    def empty_slot_single(self):
	global dgcounter                                        #DG EDIT HERE========================================================================
	global dgcounter3                                       #DG EDIT HERE========================================================================
        cur_slot=self.pos % self.max_window
        while len(self.read_table[cur_slot])>0:
            self.tot+=1
            read=self.read_table[cur_slot].pop()
            self.num_reads-=1
            seqs=self.check_for_snps(read,0)
            num_seqs=len(seqs)
            if num_seqs==0:
                dgcounter3+=1                                   #DG EDIT HERE=========================================================================
                continue
	    if num_seqs>10:
		dgcounter+=1                                    #DG EDIT HERE=========================================================================
		continue
            if num_seqs==1:
                self.keep_bam.write(read)
            else:
                self.remap_num_file.write("%i\n"%(num_seqs-1))
                self.remap_num_file.flush()
                self.remap_bam.write(read)
                for seq in seqs[1:]:
                    dg_match=0                                  #DG EDIT HERE=========================================================================
                    for cigar in read.cigar:                    #DG EDIT HERE=========================================================================
                        if cigar[0] == 0:                       #DG EDIT HERE=========================================================================
                            dg_match=cigar[1]                   #DG EDIT HERE=========================================================================
                    loc_line="%i:%s:%i:%i:%i" % (self.remap_num,self.chr_name,read.pos,num_seqs-1,dg_match)     #DG EDIT HERE=========================
                    #loc_line="%i:%s:%i:%i" % (self.remap_num,self.chr_name,read.pos,num_seqs-1)                #DG EDIT HERE=========================

                    #self.fastqs[0].write("@%s\n%s\n+%s\n%s\n"%(loc_line,seq,loc_line,read.qual))                               #DG EDIT HERE=========
                    if 16 & read.flag == 16:                                                                                    #DG EDIT HERE=========
                        self.fastqs[0].write("@%s\n%s\n+%s\n%s\n"%(loc_line,self.reverse_complement(seq),loc_line,read.qual))   #DG EDIT HERE=========
                    else:                                                                                                       #DG EDIT HERE=========
                        self.fastqs[0].write("@%s\n%s\n+%s\n%s\n"%(loc_line,seq,loc_line,read.qual))                            #DG EDIT HERE=========

                self.remap_num+=1
        #if self.printstats:
        #    sys.stderr.write(str(self.tot)+" "+str(self.nosnp)+" "+str(self.remap)+" "+str(self.toss)+"\n")
        #    self.printstats = False
        #sys.stderr.write(str(self.ref_match)+" "+str(self.alt_match)+" "+str(self.no_match)+"\n")
        self.pos+=1
        self.shift_SNP_table()


    # Processes all reads that map to the current position and removes them from the read table
    # Treats reads as paired-end
    def empty_slot_paired(self):
	cur_slot=self.pos % self.max_window
        while len(self.read_table[cur_slot])>0:  #While there are reads in this slot
            read=self.read_table[self.pos % self.max_window].pop() #Pop the first read in the slot
            self.num_reads-=1
            pair_chr_num=read.rnext #Figure out the matching read position
            pair_pos=read.mpos
            if pair_chr_num != self.chr_num or pair_pos-self.pos > self.max_window:
                continue
            pair_slot=pair_pos % self.max_window # Find the slot the matching read in
            for indx in range(len(self.read_table[pair_slot])):
                #if self.read_table[pair_slot][indx].qname==read.qname:
                if self.read_table[pair_slot][indx].qname.split(":")[-1]==read.qname.split(":")[-1]: #for testing purposes
                    pair_read=self.read_table[pair_slot].pop(indx)
                    self.num_reads-=1
                    seq1s=self.check_for_snps(read,0)
                    seq2s=self.check_for_snps(pair_read,read.mpos-read.pos)
                    num_seqs=len(seq1s)*len(seq2s)
                    if num_seqs==0 or num_seqs>32:
                        break
                    if num_seqs==1:
                        self.keep_bam.write(read)
                        self.keep_bam.write(pair_read)
                    else:
                        self.remap_bam.write(read)
                        self.remap_bam.write(pair_read)
                        self.remap_num_file.write("%i\n"% (2*(num_seqs-1)))
                        first=True
                        for seq1 in seq1s:
                            for seq2 in seq2s:
                                if not first:
                                    left_pos=min(read.pos,pair_read.pos)
                                    right_pos=max(read.pos,pair_read.pos)
                                    loc_line="%i:%s:%i:%i:%i" % (self.remap_num,self.chr_name,left_pos,right_pos,num_seqs-1)
                                    self.fastqs[0].write("@%s\n%s\n+%s\n%s\n"%(loc_line,seq1,loc_line,read.qual))
                                    self.fastqs[1].write("@%s\n%s\n+%s\n%s\n"%(loc_line,self.reverse_complement(seq2),loc_line,pair_read.qual))
                                first=False
                        self.remap_num+=1
                    break # stop searching for the pair since it was found

        #sys.stderr.write(str(self.ref_match)+" "+str(self.alt_match)+" "+str(self.no_match)+" "+str(self.toss)+"\n")
        self.pos+=1
        self.shift_SNP_table()


    # Checks a single aligned read for overlapping SNPs and created alternative sequences for remapping
    def check_for_snps(self,read,start_dist):
	global dgcounter2                                                               #DG EDIT HERE=============================================================================
        global dgcounter4                                                               #DG EDIT HERE=============================================================================
        indx=read.pos % self.max_window
        p=0
        num_snps=0
        seg_len=start_dist
        seqs=[read.seq]
        if start_dist>0:
            has_junc=False
        for cigar in read.cigar:
            seg_len+=cigar[1]
            if seg_len>self.max_window:
                sys.stderr.write("Segment distance (from read pair and junction separation) is too large. A read has been thrown out. Consider increasing the max window size.\n")
                return([])
            if cigar[0] == 4:                                                           #DG EDIT HERE=============================================================================
                p=p+cigar[1]                                                            #DG EDIT HERE=============================================================================
            elif cigar[0] == 0:                                                         #DG EDIT HERE=============================================================================
            #if cigar[0] in (0,4): #if CIGAR indicates a match alignment to the reference genome or a soft clipping
                for i in range(cigar[1]):
                    if len(self.indel_table[indx])==0:
                        snp=self.snp_table[indx]
                        if snp!=0:
                            num_snps+=1
                            if num_snps>10: #If there are more than 10 snps overlapping, throw out the read to prevent memory blow-up
                                dgcounter2+=1                                           #DG EDIT HERE=============================================================================
				return([])
                            init_seqs=list(seqs)
                            for seq in init_seqs:
                                matches=0
                                if seq[p] not in snp.alleles:
                                    dgcounter4+=1
                                    #sys.stderr.write(str(read.pos)+" "+seq[p]+"  "+str(snp.alleles)+"\n")              #DG EDIT HERE=============================================
                                    return([])                                                                          #DG EDIT HERE=============================================
                                    #sys.stderr.write(str(start_dist)+" "+seq[p]+"  "+str(snp.alleles)+"\n")            #DG EDIT HERE=============================================
                                for geno in snp.alleles:
                                    if seq[p]==geno:
                                        matches+=1
                                        for alt_geno in snp.alleles:
                                            if not alt_geno == geno:
                                                new_seq=seq[:p]+alt_geno+seq[p+1:]
                                                seqs.append(new_seq)
                                if matches==0:
                                    self.no_match+=1
                                else:
                                    self.ref_match+=1
                    else:  #it's an indel, throw it out
                        self.toss+=1
                        return([])
                    indx=(indx+1) % self.max_window
                    p+=1
            elif cigar[0]==3:   # if it is skipped in the reference genome (splice junction)
                indx=(indx+cigar[1]) % self.max_window
                has_junc=True
            else: #if there is a non-N/M in the read CIGAR, throw out the read
                self.toss+=1
                return([])
        #sys.stderr.write(str(len(seqs))+" "+str(p)+"\n")
        if len(seqs)==1:
            self.nosnp+=1
        else:
            self.remap+=1
        return seqs

    # Shifts the SNP table over one position and makes sure that indels are not lost
    def shift_SNP_table(self):
        ### Current slot to fill is the position + max_window - 1
        cur_slot=(self.pos-1)%self.max_window

        ### Delete indels that are no longer used (if they ended at the previous position)
        for indel_pos in self.indel_table[cur_slot]:
            #sys.stderr.write(str(indel_pos+self.indel_dict[indel_pos].max_len-1)+"\t"+str(self.pos-1)+"\t"+str(self.indel_dict[indel_pos].max_len)+"\n")
            if indel_pos+self.indel_dict[indel_pos].max_len-1==self.pos-1:
                del self.indel_dict[indel_pos]

        self.indel_table[cur_slot]=[]
        ### Carry over indels from the previous slot
        for indel_pos in self.indel_table[cur_slot-1]:
            if indel_pos+self.indel_dict[indel_pos].max_len-1>=self.pos+self.max_window-1:
                self.indel_table[cur_slot].append(indel_pos)

        if self.snp_table[cur_slot]!=0:
            self.num_snps-=1
        self.snp_table[cur_slot]=0

        #### See if there is a SNP overlapping the current spot
        while not self.end_of_snp_file and self.pos+self.max_window-1 > self.cur_snp.pos:
            sys.stderr.write(str(self.num_snps)+" "+str(self.pos)+" "+str(self.cur_snp.pos)+" !!!\n")
            sys.stderr.write("SNP out of order has been skipped\n")
            self.get_next_snp()

        while not self.end_of_snp_file and self.cur_snp.pos==self.pos+self.max_window-1:
            if self.cur_snp.ptype=="snp":
                self.add_snp()
            else:
                self.add_indel()
                if not self.cur_snp.pos in self.indel_table[cur_slot]:
                    self.indel_table[cur_slot].append(cur_snp.pos)
            self.get_next_snp()

    # Completely empties the read_table by repeatedly calling empty_slot function
    def empty_table(self):
        end_pos=self.pos+self.max_window
        while self.pos < end_pos:
            if self.is_paired_end:
                self.empty_slot_paired()
            else:
                self.empty_slot_single()

    def complement(self,letter):
        if letter=='A':
            return('T')
        elif letter=='T':
            return('A')
        elif letter=='C':
            return('G')
        elif letter=='G':
            return('C')
        else:
            return(letter)

    def reverse_complement(self,read):
        reverse=""
        for letter in read:
            reverse=self.complement(letter)+reverse
        return reverse

def main():
    global dgcounter                                                    #DG EDIT HERE=============================================================================
    global dgcounter2                                                   #DG EDIT HERE=============================================================================
    global dgcounter3                                                   #DG EDIT HERE=============================================================================
    global dgcounter4                                                   #DG EDIT HERE=============================================================================
    parser=argparse.ArgumentParser()
    parser.add_argument("-p", action='store_true', dest='is_paired_end', default=False)
    parser.add_argument("-m", action='store', dest='max_window', type=int, default=100000)
    parser.add_argument("infile", action='store')
    parser.add_argument("snp_dir", action='store')
    parser.add_argument("pref", action='store')

    options=parser.parse_args()
    infile=options.infile
    snp_dir=options.snp_dir
    name_split=infile.split(".")

    ####YQ EDIT HERE####
    #if len(name_split)>1:
    #    pref=".".join(name_split[:-1])
    #else:
    #    pref=name_split[0]
    pref =  options.pref
    ####################

    pysam.sort("-o", pref+".sort.bam", infile)

    sort_file_name=pref+".sort.bam"
    keep_file_name=pref+".keep.bam"
    remap_name=pref+".to.remap.bam"
    remap_num_name=pref+".to.remap.num.gz"

    if options.is_paired_end:
        fastq_names=[pref+".remap.fq1.gz",pref+".remap.fq2.gz"]
    else:
        fastq_names=[pref+".remap.fq.gz"]

    bam_data=Bam_scanner(options.is_paired_end,options.max_window,sort_file_name,keep_file_name,remap_name,remap_num_name,fastq_names,snp_dir)
    bam_data.fill_table()
    #i=0
    while not bam_data.end_of_file:
        #i+=1
        #if i>50000:
            #sys.stderr.write(str(asizeof.asizeof(bam_data))+"\t"+str(asizeof.asizeof(bam_data.snp_table))+"\t"+str(asizeof.asizeof(bam_data.read_table))+"\t"+str(bam_data.num_reads)+"\t"+str(bam_data.num_snps)+"\n")
            #sys.stderr.write(str(asizeof.asizeof(bam_data))+"\t"+str(bam_data.num_reads)+"\t"+str(bam_data.num_snps)+"\t"+str(len(bam_data.indel_dict))+"\n")
            #i=0
        if options.is_paired_end:
            bam_data.empty_slot_paired()
        else:
            bam_data.empty_slot_single()
        bam_data.fill_table()
    temp_dg = str(dgcounter)                                                            #DG EDIT HERE=============================================================================
    temp_dg2 = str(dgcounter2)                                                          #DG EDIT HERE=============================================================================
    temp_dg3 = str(dgcounter3)                                                          #DG EDIT HERE=============================================================================
    temp_dg4 = str(dgcounter4)                                                          #DG EDIT HERE=============================================================================
    sys.stderr.write("\nReads with >10 computed sequences (eliminted):\t")              #DG EDIT HERE=============================================================================
    print >> sys.stderr, temp_dg                                                        #DG EDIT HERE=============================================================================
    sys.stderr.write("\nReads with 0 computed sequences (eliminated):\t")               #DG EDIT HERE=============================================================================
    print >> sys.stderr, temp_dg3                                                       #DG EDIT HERE=============================================================================
    sys.stderr.write("\nReads overlapping >10 SNPs (given 0 computed seqs):\t")         #DG EDIT HERE=============================================================================
    print >> sys.stderr, temp_dg2                                                       #DG EDIT HERE=============================================================================
    sys.stderr.write("\nSNP positions that don't match given alleles:\t")               #DG EDIT HERE=============================================================================
    print >> sys.stderr, temp_dg4                                                       #DG EDIT HERE=============================================================================
    #sys.stderr.write(str(dgcounter))                                                   #DG EDIT HERE=============================================================================
    sys.stderr.write("\nFinished!\n")                                                   #DG EDIT HERE=============================================================================

main()
