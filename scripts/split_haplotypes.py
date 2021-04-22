
# %load "/oasis/tscc/scratch/lpatel/02_aligned/pipeline_base/scripts/split_haplotypes.py"
# split haplotypes
import sys, os, argparse, pysam
from collections import defaultdict

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

def cmdline_args():
    # Make parser object
    p = argparse.ArgumentParser(description=
        """
        This is a test of the command line argument parser in Python.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    p.add_argument("input_bam", help="The input BAM/SAM file to be split")
    p.add_argument("-ref", default="b6", help="Label for reference strain (default: b6)")
    p.add_argument("-alt", default="cast", help="Label for alternate strain (default: cast)")
    p.add_argument("-v", "--verbosity", type=int, choices=[0,1,2], default=0, help="increase output verbosity")
    p.add_argument("--sam", default=False, action="store_true", help="Set to True if using SAM")

    args = p.parse_args()
    print(args.input_bam)
    
    if not is_valid_file(p, args.input_bam):
        p.error("The file %s does not exist!" % args.input_bam)
    elif args.input_bam[-4:] != ".bam" and args.input_bam[-4:] != ".sam":
        p.error("Strange input extension: %s" % args.input_bam)
    elif args.input_bam[-4:] == ".sam" and not args.sam:
        p.error("Missing '--sam' flag, (set to True if using SAM)")
        
    return(args)

def split_haplotypes(args):
    
    BINARY_FLAG = ""
    if not args.sam:
        BINARY_FLAG = "b"
    
    reference_label = args.ref
    alternate_label = args.alt
    homologous_label = "homo"
    ambiguous_label = "amb"
    unassignable_label = "un"
    
    split_input = args.input_bam.split(os.extsep, 1)
    out_prefix, out_suffix = split_input[0], os.path.splitext(split_input[1])[0]
    print(out_prefix, reference_label, out_suffix, ("sam" if args.sam else "bam"))
   
    ref_path = f"{out_prefix}.{reference_label}.{out_suffix}.{'sam' if args.sam else 'bam'}"
    alt_path = f"{out_prefix}.{alternate_label}.{out_suffix}.{'sam' if args.sam else 'bam'}"
    homo_path = f"{out_prefix}.{homologous_label}.{out_suffix}.{'sam' if args.sam else 'bam'}"
    amb_path = f"{out_prefix}.{ambiguous_label}.{out_suffix}.{'sam' if args.sam else 'bam'}"
    un_path = f"{out_prefix}.{unassignable_label}.{out_suffix}.{'sam' if args.sam else 'bam'}"
    
    bam_file = pysam.AlignmentFile(args.input_bam, "r%s" % BINARY_FLAG)
    ref_file = pysam.AlignmentFile(ref_path, "w%s" % BINARY_FLAG, template=bam_file)
    alt_file = pysam.AlignmentFile(alt_path, "w%s" % BINARY_FLAG, template=bam_file)
    homo_file = pysam.AlignmentFile(homo_path, "w%s" % BINARY_FLAG, template=bam_file)
    amb_file = pysam.AlignmentFile(amb_path, "w%s" % BINARY_FLAG, template=bam_file)
    un_file = pysam.AlignmentFile(un_path, "w%s" % BINARY_FLAG, template=bam_file) 
    
    r1 = None
    r2 = None
    read_pairs = 0
    for read in bam_file:
        # from https://www.biostars.org/p/306041/
        if not read.is_paired or read.mate_is_unmapped or read.is_duplicate:
              continue
        if read.is_read2:
            r2 = read
        else:
            r1 = read
            r2 = None
            continue
        if not r1 is None and not r2 is None and r1.query_name == r2.query_name:
            # read pair found
            
            r1_refs = int(r1.get_tag("YB"))
            r1_alts = int(r1.get_tag("YC"))
            r2_refs = int(r2.get_tag("YB"))
            r2_alts = int(r2.get_tag("YC"))
            
            # ref group criteria:
            # read1 has reference alleles, read2 has reference alleles
            # read1 has reference alleles, read2 has no alleles
            # read1 has no alleles, read2 has reference alleles
            if ((r1_refs > 0 and r1_alts == 0 and r2_refs > 0 and r2_alts == 0) or \
                (r1_refs > 0 and r1_alts == 0 and r2_refs == 0 and r2_alts == 0) or \
                (r1_refs == 0 and r1_alts == 0 and r2_refs > 0 and r2_alts == 0)):
                ref_file.write(r1)
                ref_file.write(r2)
            
            # alt group criteria:
            # read1 has alternate alleles, read2 has alternate alleles
            # read1 has alternate alleles, read2 has no alleles
            # read1 has no alleles, read2 has alternate alleles
            if ((r1_refs == 0 and r1_alts > 0 and r2_refs == 0 and r2_alts > 0) or \
                (r1_refs == 0 and r1_alts > 0 and r2_refs == 0 and r2_alts == 0) or \
                (r1_refs == 0 and r1_alts == 0 and r2_refs == 0 and r2_alts > 0)):
                alt_file.write(r1)
                alt_file.write(r2)
            
            # homologous group criteria:
            # read1 has reference alleles, read2 has alternate alleles
            # read1 has alternate alleles, read2 has reference alleles
            if ((r1_refs > 0 and r1_alts == 0 and r2_refs == 0 and r2_alts > 0) or \
                (r1_refs == 0 and r1_alts > 0 and r2_refs > 0 and r2_alts == 0)):
                homo_file.write(r1)
                homo_file.write(r2)
            
            # ambiguous group criteria:
            # read1 has no alleles, read2 has no alleles
            if ((r1_refs == 0 and r1_alts == 0 and r2_refs == 0 and r2_alts == 0) or \
                (r1_refs == 0 and r1_alts > 0 and r2_refs == 0 and r2_alts == 0) or \
                (r1_refs == 0 and r1_alts == 0 and r2_refs == 0 and r2_alts > 0)):
                amb_file.write(r1)
                amb_file.write(r2)
            
            # unassignable group criteria:
            # read1 has both reference and alternate alleles, read2 has both reference and alternate alleles
            # read1 has both reference and alternate alleles, read2 has reference alleles
            # read1 has both reference and alternate alleles, read2 has alternate alleles
            # read1 has both reference and alternate alleles, read2 has no alleles
            # read1 has reference alleles, read2 has both reference and alternate alleles
            # read1 has alternate alleles, read2 has both reference and alternate alleles
            # read1 has no alleles, read2 has both reference and alternate alleles
            if ((r1_refs > 0 and r1_alts > 0 and r2_refs > 0 and r2_alts > 0) or \
                (r1_refs > 0 and r1_alts > 0 and r2_refs > 0 and r2_alts == 0) or \
                (r1_refs > 0 and r1_alts > 0 and r2_refs == 0 and r2_alts > 0) or \
                (r1_refs > 0 and r1_alts > 0 and r2_refs == 0 and r2_alts == 0) or \
                (r1_refs > 0 and r1_alts == 0 and r2_refs > 0 and r2_alts > 0) or \
                (r1_refs == 0 and r1_alts > 0 and r2_refs > 0 and r2_alts > 0) or \
                (r1_refs == 0 and r1_alts == 0 and r2_refs > 0 and r2_alts > 0)):
                un_file.write(r1)
                un_file.write(r2)

if __name__ == '__main__':
    
    if sys.version_info<(3,0,0):
        sys.stderr.write("You need python 3.0 or later to run this script\n")
        sys.exit(1)
        
    try:
        args = cmdline_args()
    except Exception as e:
        print(e)
        print('Try $python <script_name> <bam_file>')
        
    split_haplotypes(args)
