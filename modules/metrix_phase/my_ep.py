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
import pickle
import traceback
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
    '--output_dir',
    type=str,
    dest='output_dir',
    default='',
    help='Specify output directory')

#  parser.add_argument(
#    '--type',
#    type=str,
#    dest='type',
#    default='',
#    help='SHELXE type to run')

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
  date_dir = os.path.join(output_dir, "EP_phasing_"+date)
  os.makedirs(date_dir, exist_ok = True)
  return date_dir

###############################################################################
#
#  define class
#
###############################################################################

class ExperimentalPhasing(object):
  def __init__(self, date_dir, mtz_dir, seq_dir):
   self.date_dir = date_dir
   self.mtz_dir = mtz_dir
   self.seq_dir = seq_dir

  def __call__(self, pdb_id):
    print('Processing anomalous data for %s \n' %(pdb_id))
    logging.info(f"Processing anomalous data for {pdb_id}")

    os.mkdir(os.path.join(self.date_dir, pdb_id))
    os.chdir(os.path.join(self.date_dir, pdb_id)) 

    pdb_dir = os.path.join(self.mtz_dir, pdb_id)
    mtz_dir = os.path.join(pdb_dir, 'DataFiles')
    
    seq = os.path.join(self.seq_dir, pdb_id+".fasta")
    
    __location__ = os.path.realpath(os.path.join(os.getcwd(),
                                    os.path.dirname(__file__)))
                                    
    shelxcde = os.path.join(__location__, "util/simple_xia2_to_shelxcde_new.py")

    command = [shelxcde,
                 '--xia2dir=%s' %(pdb_dir),
                 '--name=%s' %(pdb_id),
                 '--atom=Se',
                 '--seqin=%s' %(seq),
                 '--type=trace'
                 ]

    print(' '.join(command))
    logging.info(f"{command}")
    result = procrunner.run(command)

    if result is not None:
        with open('simple_xia2_to_shelxcde.log', 'w') as f:
          f.writelines(result.stdout_lines)
          f.write('\n')
#    except Exception:
#      pickle.dump((e, traceback.format_exc()), open("%s_error.pickle" % pdb_id, "w"))



###############################################################################
#
#  run the class and functions
#
###############################################################################

def run():
  logging.basicConfig(filename="integration.log", level=logging.INFO)

  args = parse_command_line()
  
  pdb_id_list = []
  for d in os.listdir(args.mtz_dir):
    pdb_id_list.append(d)
  
  num_jobs = len(pdb_id_list)
  print("Number of Shelx jobs: ", num_jobs)
  logging.info(f"Number of Shelx jobs:{num_jobs}")
  
  date_dir = make_output_folder(args.output_dir) 
  os.chdir(date_dir)
  
  iterable = [line.strip() for line in pdb_id_list]

  with Pool(processes = int(args.processes)) as pool:
    pool.map(ExperimentalPhasing(date_dir, args.mtz_dir, args.seq_dir), iterable, int(args.chuncksize))

  os.chdir(date_dir)
