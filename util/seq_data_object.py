#!/usr/bin/env ccp4-python
# -*- coding: utf-8 -*-

from __future__ import division



class SeqData(object):
  '''Store items of interest from a FASTA sequence file; check provided sequence file
     and make sure it exists and contains protein sequences; pick the sequence which
     is longest and is protein; do a pairwise alignment between the provided target
     sequence and the model sequence which was isolated from PDB file; calculate
     similarity as percentage and pass on to Phaser for MR; count number of methionenes
     M to get number of sites to search for when running EP on Se-SAD data''' 
  
  def __init__(self, filename):
    import os.path
    from Bio import SeqIO, pairwise2, SeqRecord
    
    #check file exists and is a sequence file
    if filename is None:
  	  raise RuntimeError('Need to specify seqin filename')
    elif not os.path.exists(filename):
      raise RuntimeError('%s does not exist' % filename)

    #self.seq_file = list(SeqIO.parse(filename, "fasta"))
    
    first_record = SeqIO.parse(filename, "fasta").next()
   
    # only accept files with a single sequence
    #records = [e for e in self.seq_file]
    #assert len(records) == 1
    if 'X' in first_record:
      first_record_noX = first_record.replace("X", "")
      self.seq = str(first_record_noX.seq)
    else:
      self.seq = str(first_record.seq)
   	
    #self.seq = str(records[0].seq)
    #self.seq = str(first_record.seq) #this would be the original line without
				      #replacing "X"
   	
    return
    
  def get_seq_identity(self, search):
    from Bio import pairwise2

    aligns = pairwise2.align.globalxx(self.seq, search)
    scores = [a[2] for a in aligns]
  
  #from Bio.pairwise2 import format_alignment
  #for a in aligns:
  #  print format_alignment(*a)
    return 100 * max(scores) / len(self.seq)

  def mol_weight(self):
    from Bio.SeqUtils import molecular_weight
    return molecular_weight(self.seq, 'protein')
    
    
  #count the number of methionins in a sequence file; to be used for Se-SAD  
  def num_methionin(self):
    return self.seq.count('M')
    
    
    
    
    