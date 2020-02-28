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
    '--image_dir', 
    type=str, 
    dest="image_dir",
    default="",
    help='Give toplevel directory holding a collection diffraction images')
    
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
  if args.image_dir == '':
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
  date_dir = os.path.join(output_dir, "integration_"+date)
  os.makedirs(date_dir, exist_ok = True)
  return date_dir

###############################################################################
#
#  define class
#
###############################################################################

class IntegrateData(object):
  def __init__(self, date_dir, image_dir):
   self.date_dir = date_dir
   self.image_dir = image_dir

  def __call__(self, pdb_id):
    print('Processing anomalous data for %s \n' %(pdb_id))
    logging.info(f"Processing anomalous data for {pdb_id}")

    os.mkdir(os.path.join(self.date_dir, pdb_id))
    os.chdir(os.path.join(self.date_dir, pdb_id)) 

    pdb_dir = os.path.join(self.image_dir, pdb_id)
    
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__))
    )
    xia2 = os.path.join(__location__, "util/shell_scripts/xia2.sh")

    command = [xia2,
               'pipeline=dials',
               'atom=Se',
               pdb_dir,
               'dials.find_spots.global_threshold=Auto']

    print(' '.join(command))
    logging.info(f"{command}")
    result = procrunner.run(command)

###############################################################################
#
#  run the class and functions
#
###############################################################################

def run():
  logging.basicConfig(filename="integration.log", level=logging.INFO)

  args = parse_command_line()
  
  image_dir = args.image_dir
  pdb_id_list = []
  for d in os.listdir(args.image_dir):
    pdb_id_list.append(d)
  
  num_jobs = len(pdb_id_list)
  print("Number of xia2 jobs: ", num_jobs)
  logging.info(f"Number of xia2 jobs:{num_jobs}")
  
  date_dir = make_output_folder(args.output_dir) 
  os.chdir(date_dir)
  
  iterable = [line.strip() for line in pdb_id_list]

  with Pool(processes = int(args.processes)) as pool:
    pool.map(IntegrateData(date_dir, image_dir), iterable, int(args.chuncksize))

  os.chdir(date_dir)
