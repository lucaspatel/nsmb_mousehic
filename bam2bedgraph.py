# author: lucas
# This script takes a bedgraph file, generally representing a ChIP-Seq run, as well as an annotated compartment bedgraph file with an A/B call in the fifth column. Should be used in conjunction with output from 'analyze_hotspots.py'.

# TODO: This script could easily be combined with 'analyze_hotspots.py' to make a useful tool for general compartment calls to A/B assignments, agnosting of data type. Would be a great general tool.
# usage: python analyze_bedgraph.py <path to bedgraph> <path to annotated bedgraph-like compartments> <path to output>

import sys
from collections import defaultdict

comps = defaultdict(dict)

def main(bedgraph, compartments, out_path):
  binsize = 0

  #infer binsize
  with open(compartments, 'r') as f:
    while "NaN" in next(f):
      next(f)
    l = next(f)
    bits = l.split('\t')
    binsize = int(bits[2])-int(bits[1])

  #build A/B index
  with open(compartments, 'r') as f:
    for line in f:
      bits = line.split('\t')
      comps[str(bits[0])][int(bits[1])] = bits[4]
 
  #read bedgraph and write output
  out = open(out_path, 'w+')

  #bedgraph
  with open(bedgraph, 'r') as f:
    for line in f:
      bits = line.split('\t')
      position = (int(bits[2])+int(bits[1]))/2
      position_val = (position/binsize)*binsize

      if bits[0] not in comps:
	continue
      out.write(line.strip('\n') + '\t' + comps[str(bits[0])][position_val])

  #sam
  with open(bedgraph, 'r') as f:
    pass
    for line in f:
      bits = line.split('\t')
      position = (int(bits[2])+int(bits[1]))/2
      position_val = (position/binsize)*binsize

if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2], sys.argv[3])
