#!/bin/env python3

###############################################################################
#
#  library import
#
###############################################################################

import os
import datetime
import argparse
import procrunner
import logging
from multiprocessing import Pool
from pathlib import Path

###############################################################################
#
#  define command line arguments
#
###############################################################################

def parse_command_line():
  '''defining the command line input to make it runable'''
  parser = argparse.ArgumentParser(
           description='Integration of diffraction images with xia2 and dials')
  
  parser.add_argument(
    '--mtz_dir', 
    type=str, 
    dest="mtz_dir",
    default="",
    help='Give toplevel directory containing MTZ files from integration')

  parser.add_argument(
    '--seq_dir',
    type=str,
    dest='seq_dir',
    default='',
    help='Specify directory containing FASTA sequence(s)')

  parser.add_argument(
    '--model_dir',
    type=str,
    dest='model_dir',
    default='',
    help='Specify directory containing PDB search model(s)')

  parser.add_argument(
    '--lookup_lst',
    type=str,
    dest='lookup_lst',
    default='',
    help='Specify lookup table for target-search PDB pairs')

  parser.add_argument(
    '--mode',
    type=str,
    dest='mode',
    default='AUTO',
    help='Specify Phaser run mode; AUTO for full space group search; MR_ELLG for eLLG calculation')
    
  parser.add_argument(
    '--output_dir',
    type=str,
    dest='output_dir',
    default='',
    help='Specify output directory')

  parser.add_argument(
    '--processes',
    type=str,
    dest='processes',
    default='',
    help='Number of processes to use in multiprocessing')

  parser.add_argument(
    '--chuncksize',
    type=str,
    dest='chuncksize',
    default='',
    help='Number of datasets to process in parallel')

  args = parser.parse_args()
  if args.mtz_dir == '':
    parser.print_help()
    exit(0)
  return args

###############################################################################
#
#  create output directory
#
###############################################################################

def make_output_folder(output_dir):
  date = datetime.datetime.today().strftime('%Y%m%d')
  date_dir = os.path.join(output_dir, "MR_phasing_"+date)
  os.makedirs(date_dir, exist_ok = True)
  return date_dir

###############################################################################
#
#  define class
#
###############################################################################

class MolecularReplacement(object):
  def __init__(self, date_dir, mtz_dir, seq_dir, model_dir, lookup_lst, mode):
   self.date_dir = date_dir
   self.mtz_dir = mtz_dir
   self.seq_dir = seq_dir
   self.model_dir = model_dir
   self.lookup_lst = lookup_lst
   self.mode = mode

  def __call__(self, pdb_id):
    print('Running molecular replacement for %s \n' %(pdb_id))
    logging.info(f"Running molecular replacement for {pdb_id}")

    os.mkdir(os.path.join(self.date_dir, pdb_id))
    os.chdir(os.path.join(self.date_dir, pdb_id)) 

    target_dir = os.path.join(self.mtz_dir, pdb_id)
    data_dir = os.path.join(target_dir, 'DataFiles')
    mtzdata = os.path.join(data_dir, "AUTOMATIC_DEFAULT_free.mtz")
    
    seq = os.path.join(self.seq_dir, pdb_id+".fasta")
    
    lookup = dict(line.strip().split(',') for line in open(self.lookup_lst).readlines())

    search_id = lookup[pdb_id]
    search_model = os.path.join(self.model_dir, search_id+".pdb")
    
    __location__ = os.path.realpath(os.path.join(os.getcwd(),
                                    os.path.dirname(__file__)))
                                    
    simplemr = os.path.join(__location__, "util/simpleMR.py")

    command = [simplemr,
                 '--xyzin=%s' %(search_model),
                 '--hklin=%s' %(mtzdata),
                 '--seqin=%s' %(seq),
                 '--mode=%s' %(self.mode)]

    print(' '.join(command))
    logging.info(f"{command}")
    result = procrunner.run(command)

###############################################################################
#
#  run the class and functions
#
###############################################################################

def run():
  logging.basicConfig(filename="MR_phasing.log", level=logging.INFO)

  args = parse_command_line()
  
  pdb_id_list = []
  for d in os.listdir(args.mtz_dir):
    pdb_id_list.append(d)
  
  num_jobs = len(pdb_id_list)
  print("Number of MR jobs: ", num_jobs)
  logging.info(f"Number of MR jobs:{num_jobs}") 
  
  date_dir = make_output_folder(args.output_dir) 
  os.chdir(date_dir)
  
  iterable = [line.strip() for line in pdb_id_list]

  with Pool(processes = int(args.processes)) as pool:
    pool.map(MolecularReplacement(date_dir,
                                  args.mtz_dir,
                                  args.seq_dir,
                                  args.model_dir,
                                  args.lookup_lst,
                                  args.mode),
                                  iterable,
                                  int(args.chuncksize))

  os.chdir(date_dir)
