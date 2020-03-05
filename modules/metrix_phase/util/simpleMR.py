#!/bin/env python3
'''run a very simple Phaser job for MR'''

import argparse

from mtz_data_object import MtzData
from pdb_data_object import PdbData
from seq_data_object import SeqData
from matth_coeff_function_object import MattCoeff
from phaser_object import PhaserSearch
   
   
#define main function
def simpleMR(xyzin, hklin, seqin, mode):
  
  #load MTZ file and look-up header
  mtzdata = MtzData(hklin) 
  
  #load PDB file for search model
  pdbdata = PdbData(xyzin)
  
  #load sequence of target molecule
  seqdata = SeqData(seqin)
  print(seqdata.seq)
  with open("temp.seq", "w") as temp:
    temp.write(str(seqdata.seq))
  seq = "temp.seq"
  print(seq)
  
  #determine molecular weight of target from sequence
  mw = seqdata.mol_weight()
  
  #run Matthews coefficient to get number of molecules in asu
  matt_obj = MattCoeff(mw,
                       mtzdata.cell,
                       mtzdata.sg_num,
                       hklin)
                       
  nmol = matt_obj.num_molecules()

  #determine sequence identity between target and search model
  seq_identity = seqdata.get_seq_identity(pdbdata.seq)
  
  #run Phaser in space group search with number of molecules and sequence identity
  #determined above; use provided search model
  phaser_obj = PhaserSearch(xyzin,
                            hklin,
                            seq,
                            mtzdata.fcols,
                            mtzdata.qcols,
                            nmol,
                            seq_identity,
                            mode)
  
  
  return


if __name__ == '__main__':
  import os.path
  parser = argparse.ArgumentParser()
  parser.add_argument('--xyzin', help='name of input PDB file')
  parser.add_argument('--hklin', help='name of input MTZ file')
  parser.add_argument('--seqin', help='name of input sequence file')
  parser.add_argument('--mode', help='Phaser mode to use')
  args = parser.parse_args()
  

  if args.xyzin is None:
  	raise RuntimeError('Need to specify xyzin filename')
  else:
  	if not os.path.exists(args.xyzin):
  		raise RuntimeError('%s does not exist' % args.xyzin)
  
  if args.hklin is None:
  	raise RuntimeError('Need to specify hklin filename')
  else:
  	if not os.path.exists(args.hklin):
  		raise RuntimeError('%s does not exist' % args.hklin)
  		
  if args.seqin is None:
  	raise RuntimeError('Need to specify seqin filename')
  else:
  	if not os.path.exists(args.seqin):
  		raise RuntimeError('%s does not exist' % args.seqin)
  
  simpleMR(args.xyzin, args.hklin, args.seqin, args.mode)

