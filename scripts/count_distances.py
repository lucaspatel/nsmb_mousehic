
# author: lucas
# note, 11/08 changing index values to make this work, not sure what broke...
# okay I must have changed it without noting so for Yanxiao's pipeline output

import sys, numpy as np, argparse, os

def cmdline_args():
    # Make parser object
    p = argparse.ArgumentParser(description=
        """
        This is a test of the command line argument parser in Python.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    p.add_argument("input_pairs", help="The input pairs file to be counted")
    p.add_argument("--all", default=False, action="store_true", help="Set to True to count file globally")

    args = p.parse_args()
        
    return(args)

def count(args):
    pre_dedup=os.path.splitext(args.input_pairs)[0]
    pre=os.path.splitext(pre_dedup)[0]
    chro=1 if not args.all else "all"
    seen=[]
    breaks = [(10000*1.12**i) for i in range(0,100)]
    counts = [0 for i in breaks]

    with open(args.input_pairs, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            bits = line.split('\t')

            # check not same chr
            if bits[1] != bits[3]: 
                continue
            
             
            # check new chr
            if bits[1][3:] != chro and not args.all:
                
                # check unsorted
                if bits[1][3:] in seen:
                    raise ValueError("Unexpected chromosome order; input pairs likely not sorted...")
                    
                with open("%s_%s_intracounts.tsv" % (pre,chro), 'w+') as outfile:
                    outfile.write("breaks\tcounts\n")
                    for b, c in zip(breaks,counts):
                        outfile.write("%i\t%i\n" % (b,c))
                       
                seen.append(chro)
                chro = bits[1][3:]
                counts = [0 for i in breaks]
                
            diff = abs(int(bits[4])-int(bits[2])) 
            counts[np.digitize([diff],breaks)[0]] += 1
            
        with open("%s_%s_intracounts.tsv" % (pre,chro), 'w+') as outfile:
            outfile.write("breaks\tcounts\n")
            for b, c in zip(breaks,counts):
                outfile.write("%i\t%i\n" % (b,c))

    breaks = [0] + breaks[:-1]
    if len(breaks) != len(counts):
        raise ValueError("Something has gone wrong; list size mismatch...")

if __name__ == "__main__":
    # usage: python count.py input
    try:
        args = cmdline_args()
        print(args)
    except Exception as e:
        print(e)
        print('Try $python <pairs_file>')
        
    count(args)
