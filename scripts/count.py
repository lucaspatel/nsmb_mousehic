# author: lucas
# note, 11/08 changing index values to make this work, not sure what broke...
# okay I must have changed it without noting so for Yanxiao's pipeline output

import sys, numpy as np

def main(input, pre, chro):
  breaks = [(10000*1.12**i) for i in range(0,100)]
  counts = [0 for i in breaks]

  with open(input, 'r') as f:
    for line in f:
      if line.startswith('#'):
        continue
      bits = line.split('\t')
      
      #not same chr
      if bits[1] != bits[3] or \
      (bits[1] != "chr"+str(chro) or bits[3] != "chr"+str(chro)): # or bits[1] == 'chrX' or bits[3] == 'chrX':
        continue
      diff = abs(int(bits[4])-int(bits[2])) 
      print bits[2],bits[4],diff
      counts[np.digitize([diff],breaks)[0]] += 1


  breaks = [0] + breaks[:-1]

  #counts = counts[1:] + [0]
  #print counts
  #print len(breaks),len(counts)

  if len(breaks) != len(counts):
    raise ValueError("Something has gone wrong; list size mismatch...")

  with open("%s_%s_intracounts.tsv" % (pre,chro), 'w+') as f:
    f.write("breaks\tcounts\n")
    for b, c in zip(breaks,counts):
      f.write("%i\t%i\n" % (b,c))

if __name__ == "__main__":
  #usage python count.py input output chr
  main(sys.argv[1], sys.argv[2], sys.argv[3])
