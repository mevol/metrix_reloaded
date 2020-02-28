#!/bin/env python3

import os
from Bio import SeqIO, pairwise2, SeqRecord
from Bio.SeqUtils import molecular_weight



class SeqData(object):
  '''Store items of interest from a FASTA sequence file; check provided sequence file
     and make sure it exists and contains protein sequences; pick the sequence which
     is longest and is protein; do a pairwise alignment between the provided target
     sequence and the model sequence which was isolated from PDB file; calculate
     similarity as percentage and pass on to Phaser for MR; count number of methionenes
     M to get number of sites to search for when running EP on Se-SAD data''' 
  
  def __init__(self, filename):
    
    #check file exists and is a sequence file
    if filename is None:
  	  raise RuntimeError('Need to specify seqin filename')
    elif not os.path.exists(filename):
      raise RuntimeError('%s does not exist' % filename)
      
    for record_1 in SeqIO.parse(filename, "fasta"):
#      print(record_1)
      if ":A|" in record_1.id:
        first_record = record_1.seq
   
    # only accept files with a single sequence
      if 'X' in first_record:
        first_record_noX = first_record.replace("X", "")
        self.seq = first_record_noX
      else:
        self.seq = first_record
      return
    
  def get_seq_identity(self, search):
    aligns = pairwise2.align.globalxx(self.seq, search)
    scores = [a[2] for a in aligns]
    return 100 * max(scores) / len(self.seq)

  def mol_weight(self):
    return molecular_weight(str(self.seq), 'protein')
    
    
  #count the number of methionins in a sequence file; to be used for Se-SAD  
  def num_methionin(self):
    return self.seq.count('M')
