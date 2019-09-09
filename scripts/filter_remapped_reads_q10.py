import sys, pysam, gzip, pdb, argparse, pdb

parser = argparse.ArgumentParser()
parser.add_argument("-p", action='store_true', dest='is_paired_end', default=False)
parser.add_argument("orig_bam")
parser.add_argument("remap_bam")
parser.add_argument("keep_bam")
parser.add_argument("orig_num_file")

options= parser.parse_args()

orig_bam=pysam.Samfile(options.orig_bam,"rb")
remap_bam=pysam.Samfile(options.remap_bam,"rb")
keep_bam=pysam.Samfile(options.keep_bam,"wb",template=orig_bam)
orig_num_file=gzip.open(options.orig_num_file)

correct_maps=[]
end_of_file=False

# Get a list of reads that remapped correctly
try:
    remap_read=remap_bam.next()
except:
    end_of_file=True

while not end_of_file:
    chrm=remap_read.qname.strip().split(":")[1]                         #DG_EDIT==========================================remove handling of reverse read here. This issue was dealt with by modifying the find_intersecting_snps.py script to account for stand at an earlier stage.
    #if remap_read.is_reverse:                                          #DG_EDIT==========================================
    #    pos=int(remap_read.qname.strip().split(":")[3])                #DG_EDIT==========================================
    #else:                                                              #DG_EDIT==========================================
    #    pos=int(remap_read.qname.strip().split(":")[2])                #DG_EDIT==========================================
    pos=int(remap_read.qname.strip().split(":")[2])                     #DG_EDIT==========================================

    read_num=int(remap_read.qname.strip().split(":")[0])

    m_len=-1                                                            #DG_EDIT==========================================
    for dg_cigar in remap_read.cigar:                                   #DG_EDIT==========================================
         if dg_cigar[0] == 0:                                           #DG_EDIT==========================================
            m_len=dg_cigar[1]                                           #DG_EDIT==========================================
    orig_m_len=int(remap_read.qname.strip().split(":")[4])              #DG_EDIT==========================================

    #if remap_read.tid != -1 and remap_read.pos==pos and remap_bam.getrname(remap_read.tid)==chrm:                              #DG_EDIT==========================================
    if remap_read.tid != -1 and remap_read.mapping_quality >= 10 and remap_read.pos==pos and remap_bam.getrname(remap_read.tid)==chrm and m_len==orig_m_len:         #DG_EDIT==========================================
        dels=0   #Throw out the remapped read if it remapped with a deletion...for now
        for cig in remap_read.cigar:
            #if not cig[0] in (0,3,4):
            if not cig[0] in (0,4):                                     #DG_EDIT==========================================
                dels+=1
        if dels==0:
            correct_maps.append(read_num)
    try:
        remap_read=remap_bam.next()
    except:
        end_of_file=True

# Sort this list
correct_maps.sort()

#pdb.set_trace()
sys.stderr.write(str(len(correct_maps))+" reads remapped to the correct position\n")

# Pull out original aligned reads if all of the alternatives mapped correctly
map_indx=0
correct=0
end_of_file=False
line_num=1

try:
    orig_read=orig_bam.next()
except:
    end_of_file=True
try:
    orig_num=int(orig_num_file.readline().strip())
except ValueError:
    orig_num=-1

while not end_of_file and map_indx< len(correct_maps) and line_num <= correct_maps[-1]:
    if line_num < correct_maps[map_indx]:
        if orig_num==correct:
            keep_bam.write(orig_read)
        if options.is_paired_end:
            try:
                orig_read=orig_bam.next()
            except:
                sys.stderr.write("File ended unexpectedly (no pair found)")
                exit()
            if orig_num==correct:
                keep_bam.write(orig_read)

        line_num+=1
        correct=0
        try:
            orig_read=orig_bam.next()
            orig_num=int(orig_num_file.readline().strip())
        except:
            end_of_file=True
    elif line_num == correct_maps[map_indx]:
        correct+=1
        map_indx+=1
    else:
        sys.stderr.write("There was a problem with the index sorting\n")
        exit()

if orig_num==correct:
    keep_bam.write(orig_read)
