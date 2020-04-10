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
    print(keywords)

    result = procrunner.run(["shelxc", name],
                            stdin=keywords.encode("utf-8"),
                            print_stdout=True)

  return result
  
def simpleSHELXD(name):

  print("SHELXD")
  print("======")
  print("&&&&&&&&", name)
  # check required file exists
  fa = name + '_fa'
  print(55555555, fa)
  if not os.path.exists(fa + '.ins'):
    raise RuntimeError('Could not find {0}'.format(fa + '.ins'))
  else:
    with open('%s.ins' %fa) as f:
      s = f.read()
    s = s.replace('SEED 1', 'SEED 42')
    with open('%s.ins' %fa, 'w') as f:
      f.write(s)
    f.close()
    print(fa)
  
  keywords = fa  
  result = procrunner.run(["shelxd", fa],
                          #stdin=keywords.encode("utf-8"),
                          print_stdout=False)

  return result
  
def simpleSHELXE(name, solvent_frac=0.5, inverse_hand=False):
  return

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
  print(33333333, xia2_dat.sg_name)
  #gemmi_sg = gemmi.SpaceGroup(x2_dat.sg_name)
  #print(gemmi_sg)
  space_groups = generate(gemmi.SpaceGroup(xia2_dat.sg_name))
  
  print(77777777, space_groups)

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
  print(11111111, scaled_mtz)

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
    if not os.path.exists(sg_str): os.mkdir(sg_str)

    print("********", os.getcwd())
    
    os.chdir(sg_str)
    
    print("new working dir", os.getcwd())

########################################################################
###  run shelx functions
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
    try:
      with open(args.name + '_fa.res') as f:
        cfom.append(float(f.readline().split()[-1]))
    except IOError:
      cfom.append(0)

  os.chdir(cwd)
