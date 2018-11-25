# author: lucas
# This script takes compartment data in .wig format and DSB hotspots to quantify the A/B compartment correlation with hotspot signal strength.
# usage: python analyze_hotspots.py <path to hotspots> <path to compartments>

import sys

pos_score, neg_score = 0,0
pos_count, neg_count, nanz_count = 0,0,0

def get_comp(val):
  global pos_score, neg_score, pos_count, neg_count, nanz_count

  if pos_score > abs(neg_score):
    #then probably positive is A
    return "A" if val > 0 else "B"
  else:
    #then probably negative is A
    return "A" if val < 0 else "B"

def main(hotspots, compartments, out_path):
  global pos_score, neg_score, pos_count, neg_count, nanz_count

  chro = "chr" + compartments.split('_')[1]
  binsize = -1
  comp_index = {}

  #infer chromosome and binsize
  with open(compartments, 'r') as f:
    l1 = next(f)
    l2 = next(f)
    if l1.split()[0] == "variableStep":
      l1 = next(f)
    binsize = int(l2.split()[0])-int(l1.split()[0])
    if binsize == 0 or binsize % 2 != 0:
      raise ValueError("Something wrong with binsize inferrence.")
  
  #build compartments index
  with open(compartments, 'r') as f:
    for line in f:
      if line.split()[0] == "variableStep":
	continue
      bits = line.split('\t')
      comp_index[int(bits[0])] = float(bits[1].strip())

  #print comp_index

  with open(hotspots, 'r') as f:
    for peak in f:
      bits = peak.split('\t')
      if str(bits[0]) != chro:
        continue
      position = int(bits[1])
      position_val = (position/binsize)*binsize

      #assignment
      if comp_index[position_val] > 0:
        pos_score += comp_index[position_val]
        pos_count += 1
      elif comp_index[position_val] < 0:
        neg_score += comp_index[position_val]
        neg_count += 1
      else:
        nanz_count += 1
  
  #confirm bin assignment
  out_args = [out_path+":", pos_score, neg_score, pos_count, neg_count, nanz_count]
  print('\t'.join(str(x) for x in out_args))
  if pos_score > abs(neg_score) and pos_count < neg_count \
    or pos_score < abs(neg_score) and pos_count > neg_count:
    print("Ambiguous bin assignment is inconsistent for %s." % out_path)
  #check with gene transcription start sites
  genes = open("../../../../00_ref/genes_ucsc2.gtf", 'r')
  gene_A_score, gene_B_score = 0,0
  for gene in genes:
    #comment and nonchro skipping
    bits = gene.split('\t')
    if gene.startswith('#') or bits[1] != chro:
      continue
    
    #strand determination
    start = 0
    if bits[2] == '+':
      start = int(bits[3])
    elif bits[2] == '-':
      start = int(bits[4])
    else:
      raise ValueError("Strange strand assignment.")
    
    #check index
    start = (start/binsize)*binsize

    if get_comp(comp_index[start]) == "A":
      gene_A_score += 1
    else:
      gene_B_score += 1

  #report
  if gene_A_score > gene_B_score:
    print("Transcriptional start sites indicate proper A/B assignment for %s (A: %i, B: %i)" % (chro,gene_A_score,gene_B_score))
  else:
    print("Transcriptional start sites indicate INCORRECT A/B assignment for %s (A: %i, B:%i)" % (chro,gene_A_score,gene_B_score))

  #write output file
  a_count,b_count = 0,0
  out = open(out_path, 'w+')
  with open(hotspots, 'r') as f:
    for peak in f:
      bits = peak.split('\t')
      if str(bits[0]) != chro:
        continue
      position = int(bits[1])
      position_val = (position/binsize)*binsize
      if get_comp(comp_index[position_val]) == "A":
	a_count += 1
      else:
	b_count += 1

      #assignment
      if comp_index[position_val] != 0 and comp_index[position_val] != "NaN":
        out.write(peak.strip()+'\t'+str(comp_index[position_val])+'\t'+get_comp(comp_index[position_val])+'\n')
      else:
	out.write(peak.strip()+'\t'+str(comp_index[position_val])+'\n')

  out2 = open(out_path+"_whole", 'w+')
  with open("control_500000_comp.bedgraph", 'r') as f:
    for line in f:
      bits = line.split('\t')
      if bits[0] != chro:
	  continue
      out2.write(line.strip('\n') + '\t' + get_comp(comp_index[int(bits[1])])+'\n')

  print("%s has %i A hotspots and %i B hotspots. (%f active)" % (out_path,a_count,b_count,(float(a_count)/(a_count+b_count))))

if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2], sys.argv[3])
