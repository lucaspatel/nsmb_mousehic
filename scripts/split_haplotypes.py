# split haplotypes
import sys
import argparse

def cmdline_args():
    # Make parser object
    p = argparse.ArgumentParser(description=
        """
        This is a test of the command line argument parser in Python.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    p.add_argument("required_positional_arg",
                   help="input bam")
    p.add_argument("-v", "--verbosity", type=int, choices=[0,1,2], default=0,
                   help="increase output verbosity")
                   
    group1 = p.add_mutually_exclusive_group(required=True)
    group1.add_argument('--enable',action="store_true")
    group1.add_argument('--disable',action="store_false")

    return(p.parse_args())


    
# Try running with these args
#
# "Hello" 123 --enable
if __name__ == '__main__':
    
    if sys.version_info<(3,0,0):
        sys.stderr.write("You need python 3.0 or later to run this script\n")
        sys.exit(1)
        
    try:
        args = cmdline_args()
        print(args)
    except:
        print('Try $python <script_name> "Hello" 123 --enable')

    print()
