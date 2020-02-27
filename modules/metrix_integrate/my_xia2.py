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
#def make_output_folder(output_dir, pdb_id, date_dir):
  date = datetime.datetime.today().strftime('%Y%m%d')
  date_dir = os.path.join(output_dir, date)
  os.makedirs(date_dir, exist_ok = True)
  #os.chdir(date_dir)
  #pdb_out = os.path.join(date_dir, pdb_id)
  #os.makedirs(pdb_out, exist_ok=True)
  #return pdb_out
  return date_dir

###############################################################################
#
#  define class
#
###############################################################################

class IntegrateData(object):
  #def __init__(self, pdb_out, pdb_id, pdb_dir, date_dir):
  def __init__(self, date_dir, image_dir):
   self.date_dir = date_dir
   self.image_dir = image_dir
#   self.env = env

#    self.pdb_out = pdb_out
#    print(self.pdb_out)
#    self.pdb_id = pdb_id
#    print(self.pdb_id)
#    self.pdb_dir = pdb_dir
#    print(self.pdb_dir)
#    self.date_dir = date_dir
#    print(self.date_dir)
#    self.run_anomalous_xia2()
    #self.env = env

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
    
    
#    try:
#      easy_run.fully_buffered(command).raise_if_errors()
#    except Exception as e:
#      print
#    fn_txt = os.path.join(self.date_dir, pdb_id, 'xia2.txt')
#    xia2_file = open(fn_txt)
#    success = False
#    for line in xia2_file:
#      if line.startswith('Status: normal termination'):
#        print('Finished anomalous processing with xia2 for %s' %pdb_id)
#    success = True
    
    


################################################################################
#  def run_anomalous_xia2(self):
#    print("I am here", self.pdb_out)
#    os.chdir(self.pdb_out)
#    command = ["xia2",
#               "pipeline=dials",
#               "atom=Se",
#               self.pdb_dir,
#               "dials.find_spots.global_threshold=Auto"]
#
#    print(' '.join(command))
#    
#    result = procrunner.run(command,
#                            #stdin=b_keywords,
#                            print_stdout=False,
#                            timeout=5)
#    
#    assert not result['exitcode']
#    assert result['exitcode'] == 0 # alternatively
#    
#    assert not result['stderr']
#    assert result['stderr'] == b'' # alternatively
#    
#    os.chdir(self.date_dir)
#
#
#
###########################################################################    
#    try:
#      easy_run.fully_buffered(command).raise_if_errors()
#    except Exception, e:
#      print e
#    fn_txt = os.path.join(self.pdb_out, 'xia2.txt')
#    xia2_file = open(fn_txt)
#    success = False
#    for line in xia2_file:
#      if line.startswith('Status: normal termination'):
#	      print('Finished anomalous processing with xia2 for %s' %self.pdb_id)
#	      success = True
	      
	      
	      
	      
	      
#	with open(os.path.join(self.output_dir, 'anomalous_processing.txt'), 'a') as text_file:
#	  text_file.write('Finished anomalous processing with xia2 for %s \n' %pdb_id)
#	text_file.close()
#	success = True
#    if not success:
#      with open(os.path.join(self.output_dir, 'anomalous_processing.txt'), 'a') as text_file:
#	text_file.write('Anomalous processing with xia2 for %s FAILED\n' %pdb_id)
#      text_file.close()
#      print 'Anomalous processing with xia2 for %s failed' %pdb_id
    




#    result = procrunner.run(command,
#                            stdin=b_keywords,
#                            print_stdout=False,
#                            timeout=5)

#  def __call__(self):
#    print('Processing anomalous data for %s \n' %(self.pdb_id))
  
    #os.mkdir(os.path.join(self.output_dir,pdb_id))
#    os.chdir(self.pdb_out)
    #os.environ.update(self.env)
      
#    print("I am here", self.pdb_out)
    
#      command = ['xia2',
#               'pipeline=dials',
#               'atom=Se',
#               '/dls/metrix/metrix/Newcastle_data/%s' %pdb_id,
#               'dials.find_spots.global_threshold=Auto']
    
#      print(' '.join(command))

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

#  dir_list = []
#  for pdb_id in pdb_id_list:
#    pdb_dir = os.path.join(args.image_dir, pdb_id)
#    dir_list.append(pdb_dir)
  
  date_dir = make_output_folder(args.output_dir) 
  os.chdir(date_dir)
  
  iterable = [line.strip() for line in pdb_id_list]
  print(iterable)


#  for i in iterable:
#    p = Process(target=IntegrateData, args=(date_dir, image_dir))
#    p.start()


  agents = 5
  chunksize = 3 
  with Pool(processes = agents) as pool:
    pool.map(IntegrateData(date_dir, image_dir), iterable, chunksize)
  
#  print("Result: " + str(result))



#  easy_mp.parallel_map(IntegrateData(date_dir, image_dir),
#  iterable, processes=20, method='sge', qsub_command="qsub -q test-medium.q -pe smp 12")

  os.chdir(date_dir)
  

#  integrate_data = IntegrateData(date_dir, image_dir)  
  

  
  
  
  
#  date = datetime.datetime.today().strftime('%Y%m%d')
#  date_dir = os.path.join(args.output_dir, date)
#  #os.makedirs(date_dir, exist_ok = True)
  
#  dir_list = []  
  
  # Loop through the pdbs and add each entry to the database
#  for pdb_id in pdb_id_list:
#    pdb_dir = os.path.join(args.image_dir, pdb_id)
#    pdb_out = make_output_folder(args.output_dir, pdb_id, date_dir)
#    integrate_data = IntegrateData(pdb_out, pdb_id, pdb_dir, date_dir)


#    dir_list.append(pdb_dir)
#  print(dir_list)
  
#  agents = 5
#  chunksize = 3
#  with Pool(processes = agents) as pool:
#    results = pool.map(IntegrateData, dir_list, chunksize)
#  
#  print("Result: " + str(result))
  
  
  

