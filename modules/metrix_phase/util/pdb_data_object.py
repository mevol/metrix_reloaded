#!/usr/bin/env ccp4-python

import os
import gemmi

class PdbData(object):
  '''Store items of interest from a PDB search model; this does not work with
     complexes yet where there are possibly several protein chains of various lengths
     and level of completeness; the current state assumes a single protein species;
     needs fixing; does not include nucleic acids at the moment'''
  
  def __init__(self, filename):
    #check file exists and is a PDB file
    if filename is None:
  	  raise RuntimeError('Need to specify xyzin filename')
    elif not os.path.exists(filename):
      raise RuntimeError('%s does not exist' % filename)
    
    #describe hierarchy of a PDB file which defines it as model-chain-conformer
    #of aa; get sequence of it and check length of sequence and find chain which
    #is longest and isolate it as search model

    longest = 0
    self.fasta = None
    self.seq = None

    structure = gemmi.read_structure(filename)
    structure.setup_entities()
    for model in structure:
      for chain in model:
        name = chain.name
        polymer = model[name].get_polymer()
        seqlen = len(polymer)
        if seqlen > longest:
          self.seq = polymer.make_one_letter_sequence()
    
    return
