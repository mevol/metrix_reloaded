#!/usr/bin/env ccp4-python

from __future__ import division



class PdbData(object):
  '''Store items of interest from a PDB search model; this does not work with
     complexes yet where there are possibly several protein chains of various lengths
     and level of completeness; the current state assumes a single protein species;
     needs fixing; does not include nucleic acids at the moment'''
  
  def __init__(self, filename):
    import os.path
    from iotbx import pdb
    
    #check file exists and is a PDB file
    if filename is None:
  	  raise RuntimeError('Need to specify xyzin filename')
    elif not os.path.exists(filename):
      raise RuntimeError('%s does not exist' % filename)
    
    #describe hierarchy of a PDB file which defines it as model-chain-conformer
    #of aa; get sequence of it and check length of sequence and find chain which
    #is longest and isolate it as search model
    self.pdb_obj = pdb.input(file_name=filename)
    hierarchy = self.pdb_obj.construct_hierarchy()
    longest = 0
    self.fasta = None
    self.seq = None
    for model in hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          try:
	        seqlen = len(conformer.as_padded_sequence())
          except IndexError:
            continue
          if seqlen > longest:
            longest = seqlen
            self.fasta = conformer.format_fasta()
            self.seq = conformer.as_padded_sequence()
            self.chainid = chain.id   
    return