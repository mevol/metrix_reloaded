#!/bin/env python3

###############################################################################
#
#  library import
#
###############################################################################

import os
import subprocess
import datetime
import argparse
import procrunner
from multiprocessing import Pool, Process
#from libtbx import easy_run, easy_mp

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
  date_dir = os.path.join(output_dir, date)
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

    os.mkdir(os.path.join(self.date_dir, pdb_id))
    os.chdir(os.path.join(self.date_dir, pdb_id)) 

    pdb_dir = os.path.join(self.image_dir, pdb_id)

    command = ['xia2',
               'pipeline=dials',
               'atom=Se',
               pdb_dir,
               'dials.find_spots.global_threshold=Auto']

    print(' '.join(command))
    
    result = procrunner.run(command)

###############################################################################
#
#  run the class and functions
#
###############################################################################

def run():
  args = parse_command_line()
  
  image_dir = args.image_dir
  pdb_id_list = []
  for d in os.listdir(args.image_dir):
    pdb_id_list.append(d)
  
  num_jobs = len(pdb_id_list)
  print("Number of xia2 jobs: ", num_jobs)
  
  date_dir = make_output_folder(args.output_dir) 
  os.chdir(date_dir)
  
  iterable = [line.strip() for line in pdb_id_list]

  agents = 5
  chunksize = 3 
  with Pool(processes = agents) as pool:
    pool.map(IntegrateData(date_dir, image_dir), iterable, chunksize)

  os.chdir(date_dir)
