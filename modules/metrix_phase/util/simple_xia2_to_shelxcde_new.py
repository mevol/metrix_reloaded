#!/bin/env python3

import argparse
import fileinput
import sys
import shutil
import os
import procrunner
import gemmi
import re
from xia2_json_reader import xia2_json_reader
from mtz_data_object import MtzData
from seq_data_object import SeqData
from matth_coeff_function_object import MattCoeff, matt_coeff_factory
from generate_possible_spacegroups import generate  
from multiprocessing import Pool, Process

def simpleSHELXC(name, cell, wavelengths, sg, find, ntry=1000):
  print("SHELXC")
  print("======")
  cell_round = tuple(map(lambda x: isinstance(x, float) and round(x, 2) or x, cell))
  
  for w in wavelengths:
    label = w["name"]
    sca = os.path.relpath(w['sca'])
  
    keywords_shelxc = """{0} {1}
CELL {2} {3} {4} {5} {6} {7}
SPAG {8}
FIND {9}
NTRY {10}
"""
    fmt = (label,) + (sca,) + cell_round + (sg,) + (find,) + (ntry,)
    keywords = keywords_shelxc.format(*fmt)
#    print(keywords)

    result = procrunner.run(["shelxc", name],
                            stdin=keywords.encode("utf-8"),
                            print_stdout=True)

  return result
  
def simpleSHELXD(name):

  print("SHELXD")
  print("======")
  # check required file exists
  fa = name + '_fa'
  if not os.path.exists(fa + '.ins'):
    raise RuntimeError('Could not find {0}'.format(fa + '.ins'))
  else:
    with open('%s.ins' %fa) as f:
      s = f.read()
    s = s.replace('SEED 1', 'SEED 42')
    with open('%s.ins' %fa, 'w') as f:
      f.write(s)
    f.close()
  
  result = procrunner.run(["shelxd", fa],
                          print_stdout=False)

  return result
  
def simpleSHELXE(name, find, solvent_frac=0.5, inverse_hand=False):
  fa = name + '_fa'

  print(name)
  print(fa)
  print(solvent_frac)
  print(find)
  solvent = "-s{0}".format(str(solvent_frac))
  print(solvent)
  h_value = "-h{0}".format(find)
  print(h_value)
  z_value = "-z{0}".format(find)
  print(z_value)
  
  if not inverse_hand:  
    result = procrunner.run(["shelxe", name, fa, solvent, h_value, z_value, "-a5", "-q"],
                          print_stdout=True)  
    
  if inverse_hand:
    result = procrunner.run(["shelxe", name, fa, solvent, h_value, z_value, "-a5", "-q", "-i"],
                          print_stdout=True)  
  
  msg = "SHELXE - {0} hand"
  if inverse_hand:
    msg = msg.format("inverse")
  else:
    msg = msg.format("original")
  print(msg)
  print("=" * len(msg))

  # check required files exist
  if not os.path.exists(fa + '.ins'):
    raise RuntimeError('Could not find {0}'.format(fa + '.ins'))
    
  return result

  #fix to use newer shelxe; use line below when CCP4 has been updated
#  cmd = "/dls_sw/apps/shelx/64/2017-1/shelxe {0} {1} -s{2} -m200 -h{3} -z{3} -e1".format(name, fa, solvent_frac, find)
  #below is original line which is not used for this run
  #cmd = "/dls_sw/apps/shelx/64/2017-1/shelxe {0} {1} -s{2} -m -h{3} -z{3} -a5 -q".format(name, fa, solvent_frac, find)
  #cmd = "shelxe {0} {1} -s{2} -m0 -h{3} -z{3} ".format(name, fa, solvent_frac, find)#NO solvent and initial phases
  #cmd = "shelxe {0} {1} -s{2} -h{3} -z{3} ".format(name, fa, solvent_frac, find)#solvent and initial phases
  #cmd = "shelxe {0} {1} -s{2} -m20 -h{3} -z{3} ".format(name, fa, solvent_frac, find)#solvent flattening
  #cmd = "shelxe {0} {1} -s{2} -m200 -h{3} -z{3} ".format(name, fa, solvent_frac, find)#more solvent flattening
  #cmd = "shelxe {0} {1} -s{2} -m200 -h{3} -z{3} -e1".format(name, fa, solvent_frac, find)#solvent flattening and free-lunch algorithm
  
  
  #keywords_shelxe = '''{0} {1} -s{2} -m -h{3} -z{3} -a5 -q'''
  
#  fa = name + '_fa'
  
  #keywords = keywords_shelxe.format(name, fa, solvent_frac, find)
  
  #cmd = 
  
#  cmd = "shelxe {0} {1} -s{2} -m -h{3} -z{3} -a5 -q".format(name, fa, solvent_frac, find)
#  if inverse_hand: cmd += " -i"
#
#  msg = "SHELXE - {0} hand"
#  if inverse_hand:
#    msg = msg.format("inverse")
#  else:
#    msg = msg.format("original")
#  print(msg)
#  print("=" * len(msg))
#
#  # check required files exist
#  if not os.path.exists(fa + '.ins'):
#    raise RuntimeError('Could not find {0}'.format(fa + '.ins'))
#    
#  result = procrunner.run(["shelxe"], stdin=cmd.encode("utf-8"), print_stdout=False)  
#  return result

def copy_sca_locally(wavelengths):
  '''Copy .sca files locally to work around problem at DLS where SHELX fails
  to find the files'''
  for w in wavelengths:
    f = os.path.basename(w['sca'])
    shutil.copy(w['sca'], '.')
    w['sca'] = os.path.abspath(f)
  return wavelengths

if __name__ == '__main__':
########################################################################
###  receive command line arguments
########################################################################

  parser = argparse.ArgumentParser()
  parser.add_argument('--name', help='filename stem for SHELX')
  parser.add_argument('--xia2dir', help='path to a xia2 processing directory')
  parser.add_argument('--seqin', help='name of input sequence file')
  parser.add_argument('--atom', help="heavy atom (default to 'Se')")
  parser.add_argument('--ntry', help="override ntry for SHELX (default 1000)")
  parser.add_argument('--type', help="override SHELXE type (default 'trace')")
  args = parser.parse_args()

  if args.name is None:
    raise RuntimeError('Need to specify a filename stem for SHELX')

  if args.xia2dir is None:
    raise RuntimeError('Need to specify the path to the xia2 processing directory')
  else:
    if not os.path.exists(args.xia2dir):
      raise RuntimeError('%s does not exist' % args.xia2dir)

  if args.atom is None:
    print("Defaulting to --atom Se")
    args.atom = "Se"

  if args.ntry is None:
    print("Defaulting to --ntry 1000")
    args.ntry = 1000

########################################################################
###  find xia2.json and extract information
########################################################################

  # Find xia2.json
  xia2json = os.path.join(args.xia2dir, 'xia2.json')
  if not os.path.exists(xia2json):
    raise RuntimeError('%s does not exist' % xia2json)
  
  # Extract data from xia2.json
  xia2_dat = xia2_json_reader(xia2json)

########################################################################
###  get xia2 point group and create list of possible space groups
########################################################################

  # determine spacegroups for pointgroup
  space_groups = generate(gemmi.SpaceGroup(xia2_dat.sg_name))
  
########################################################################
###  copy SCA file locally
########################################################################

  # copy .sca files locally and update the wavelengths dictionary
  print("xia2 wave", xia2_dat.wavelengths)
  wl = copy_sca_locally(xia2_dat.wavelengths)
  print(99999999, wl)


########################################################################
###  identify experimental phasing wavelengths
########################################################################
#
#  # Try to identify the wavelengths given the chosen scatterer
#  xia2_dat.identify_wavelengths(args.atom)
#  print(22222222, xia2_dat.identify_wavelengths(args.atom))


## FIX ME; returns NONE atm

########################################################################
###  find scaled MTZ file
########################################################################

  # Find scaled mtz
  scaled_mtz = xia2_dat.scaled_mtz

########################################################################
###  use Matthews coefficient and sequence to get number of sites
###  to look for
########################################################################

  # If we have the sequence, make SeqData and MattCoeff objects to calculate various
  # quantities
  find = 10
  if args.seqin is None:
    seqdata = None
    matt_coeff = None
    #num_molecules = None
    print("No sequence supplied. The number of sites will default to 10")
  else:
    seqdata = SeqData(args.seqin)
    num_met = seqdata.num_methionin()
    print("Num Met: ", num_met)
    matt_coeff = matt_coeff_factory(scaled_mtz, args.seqin)

    # Set 'find' if SeMet
    if args.atom.upper() == "SE":
      find = matt_coeff.num_molecules() * seqdata.num_methionin()
      #num_molecules = matt_coeff.num_molecules()
      #print("#MET per asu = ", find)
      #print 'Num_molecules in ASU: %s \n' %num_molecules
    else:
      print("Sequence supplied, but the scatterer is not Se, so the number of "
            "sites will default to 10 anyway")

########################################################################
###  get working directory
########################################################################

  cwd = os.path.normpath(os.getcwd())
  
  print(88888888, cwd)

########################################################################
###  create sub-directories for each space group and create lists to
###  monitor phasing results
########################################################################

  cfom = []
  c_output = []

  for sg in space_groups:
    sg_str = str(sg)
    sg_str = str(sg_str).strip('<gemmi.SpaceGroup("')
    sg_str = str(sg_str).strip('")>')
    sg_str = "".join(sg_str.split())
    print("Trying {0} \n".format(sg_str))
    os.chdir(cwd)
    if "(" in sg_str:
      pass
    else:
      if not os.path.exists(sg_str):
        os.mkdir(sg_str)

      print("********", os.getcwd())
    
      os.chdir(sg_str)
    
      print("new working dir", os.getcwd())

########################################################################
###  run shelx C and D functions
########################################################################

      # Run SHELXC
      c_output.append(simpleSHELXC(args.name,
                                   xia2_dat.cell,
                                   wl,
                                   sg_str,
                                   find,
                                   args.ntry))
      
      #print(c_output)

      # Run SHELXD
      d_output = simpleSHELXD(args.name)
    
      # Get CFOM from .res
      print("********", os.getcwd())
      try:
        with open(args.name + '_fa.res', "r") as f:
          for line in f:
            if "SHELXD" in line:
              print(11111111, line)
              print(22222222, line.split()[-1])
              cfom.append(float(line.split()[-1]))
#      print(4444444, cfom)
      except IOError:
        cfom.append(0)
#    print(cfom)
#    print(c_output)
    os.chdir(cwd)

  best_idx = cfom.index(max(cfom))
  print("Best index", best_idx)
  best_sg = str(space_groups[best_idx]).strip('<gemmi.SpaceGroup("')
  best_sg = str(best_sg).strip('")>')
  best_sg = "".join(best_sg.split())
  print("Best sg", best_sg)
  print("Best space group: %s with CFOM=%s" %(str(best_sg),
                                              str(cfom[best_idx])))
  print(c_output[best_idx])
  c_output = c_output[best_idx]
  for line in c_output:
    if line.startswith(' Resl'): print(line)
    if line.startswith(' <d"/sig>'): print(line)
    if line.startswith(' CC(1/2)'): print(line)
  
  os.chdir(best_sg)
  print("********", os.getcwd())

########################################################################
###  run shelx E functions for best_sg
########################################################################

  # Run SHELXE
  solvent_frac = 0.5 if matt_coeff is None else matt_coeff.solvent_fraction(
    matt_coeff.num_molecules())
    
  # original hand  
  e_output_ori = simpleSHELXE(args.name,
                              find,
                              solvent_frac)
                              
  try:
    with open(args.name + '.lst') as f:
      for line in f.readlines():
        if line.startswith('Best trace'): print(line)
  except IOError:
    pass
  
  # inverse hand
  e_output_inv = simpleSHELXE(args.name,
                              find,
                              solvent_frac,
                              inverse_hand=True)
                              
  try:
    with open(args.name + '_i.lst') as f:
      for line in f.readlines():
        if line.startswith('Best trace'): print(line)
  except IOError:
    pass

