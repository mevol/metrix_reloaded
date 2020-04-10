#!/bin/env python3

'''Attempt to solve a structure with SHELXC/D/E given the path to a xia2.json
file and access to the processed data.'''

# This is the shelxc MAD example from http://strucbio.biologie.uni-konstanz.de/ccp4wiki/index.php/SHELX_C/D/E:
#
#   shelxc jia <<EOF
#   NAT jia_nat.hkl
#   HREM jia_hrem.sca
#   PEAK jia_peak.sca
#   INFL jia_infl.sca
#   LREM jia_lrem.sca
#   CELL 96.00 120.00 166.13 90 90 90
#   SPAG C2221
#   FIND 8
#   NTRY 10
#   EOF
#   shelxd jia_fa
#   shelxe jia jia_fa -s0.6 -m20
#   shelxe jia jia_fa -s0.6 -m20 -i

import argparse
import fileinput
import sys
import shutil
import os
import procrunner
import gemmi
import re
#from multiprocessing import Pool, Process
from xia2_json_reader import xia2_json_reader
from mtz_data_object import MtzData
from seq_data_object import SeqData
from matth_coeff_function_object import MattCoeff, matt_coeff_factory
from generate_possible_spacegroups import generate  

# Function stolen from fast_ep...
#def sanitize_spacegroup(spacegroup_name):
#  spacegroup_name.replace(' ','')
#  if not ':' in spacegroup_name:
#    return spacegroup_name
#  assert(spacegroup_name.startswith('R'))
#  return 'H%s' % spacegroup_name.split(':')[0][1:]

def simpleSHELXC(name, cell, wavelengths, sg, find, ntry=1000):
  print("SHELXC")
  print("======")
  # wavelengths should be be a dictionary with keys 'label' and 'sca'
  # prepare keywords
  #keywords = []

#  for w in wavelengths:
#    label = w['label']
#    sca = os.path.relpath(w['sca'])
#    keywords.append("{0} {1}".format(label, sca))
#    
#  cell_n = tuple(re.findall(r'\d+(?:.\d+)?', str(cell)))
#  print(cell_n)
#
#  cell_n_float = (float(cell_n[0]), float(cell_n[1]), float(cell_n[2]), float(cell_n[3]), float(cell_n[4]), float(cell_n[5]),)
#  print(cell_n_float)
#  
#  cell_n_round = tuple(map(lambda x: isinstance(x, float) and round(x, 2) or x, cell_n_float))
#  print(cell_n_round)  
#
#  keywords.append('cell {}'.format(cell))
#  keywords.append('spag {}'.format(sg))
#  keywords.append('find {}'.format(find))
#  keywords.append('ntry {}'.format(ntry))
#
#  print(11111111, keywords)
#
#  result = procrunner.run(["shelxc", name],
#           stdin=keywords.encode("utf-8"),
#           print_stdout=True)


  for w in wavelengths:
    label = w['label'].upper()
    sca = os.path.relpath(w['sca'])
    
    print(label, sca)

    cell_n = re.findall(r'\d+(?:.\d+)?', str(cell))
    print(cell_n)
    
    cell_n_round = [round(float(num), 2) for num in cell_n]
    
    print(1111111, cell_n_round)
    
    #cell_n_round = ''.join(str(cell_n_round))
    
    #print(cell_n_round)
    
#    for param in cell_n:
#      param_round = round(float(param), 2)
#      print(param_round)
    
#    cell_n_float = (float(cell_n[0]), float(cell_n[1]), float(cell_n[2]), float(cell_n[3]), float(cell_n[4]), float(cell_n[5]),)
#    print(cell_n_float)
    #cell_n_round = str(map(lambda x: isinstance(x, float) and round(x, 2) or x, cell_n_float))
    #print(cell_n_round)
    #fmt = label + sca + cell_n_round + sg + find + ntry
    
    #print(fmt)

    keywords = """
               {0} {1}
               CELL {2} {3} {4} {5} {6} {7}
               SPAG {8}
               FIND {9}
               NTRY {10}
               """

    keywords = keywords.format(label,
                              sca,
                              cell_n_round[0],
                              cell_n_round[1],
                              cell_n_round[2],
                              cell_n_round[3],
                              cell_n_round[4],
                              cell_n_round[5],
                              sg,
                              find,
                              ntry)
  print(keywords)

  result = procrunner.run(["shelxc", name],
           stdin=keywords.encode("utf-8"),
           print_stdout=True)

    
#  for w in wavelengths:
#    label = w['label'].upper()
#    sca = os.path.relpath(w['sca']) 
#    
#    keywords = '''
#               {0} {1}
#               CELL {2} {3} {4} {5} {6} {7}
#               SPAG {8}
#               FIND {9}
#               NTRY {10}
#               SFAC SE
#               '''
#    
#    print(11111111, keywords)
#    
#    cell_n = tuple(re.findall(r'\d+(?:.\d+)?', str(cell)))
#    
#    print(cell_n)
#    
#    cell_n_float = (float(cell_n[0]), float(cell_n[1]), float(cell_n[2]), float(cell_n[3]), float(cell_n[4]), float(cell_n[5]),)
#    
#    print(cell_n_float)
#    
#    cell_n_round = tuple(map(lambda x: isinstance(x, float) and round(x, 2) or x, cell_n_float))
#    
#    print(cell_n_round)
#    
#    fmt = (label,) + (sca,) + cell_n_round + (sg,) + (find,) + (ntry,)
#    
#    keywords = keywords.format(*fmt)
#  
#    print(keywords)
  
  #with open("shelxc.inp", "w") as output:
  #  output.write(keywords) 
  #
  #__location__ = os.path.realpath(os.path.join(os.getcwd(),
  #                                os.path.dirname(__file__)))
  #shelx = os.path.join(__location__,
  #                             "shell_scripts/shelx.sh")
                               
  #print(os.getcwd())
#
#  script = os.path.join(os.getcwd(), "shelxc.inp")
#  print(script)
  
  #with open("shelxc.inp", "r") as shelxc_input:
  #  
  # 
#    result = procrunner.run(["shelxc", name],
#             stdin=keywords.encode("utf-8"),
#             print_stdout=True)
  
  #shelx_name = name

  #result = procrunner.run(["shelxc", shelx_name],
  #      stdin=keywords.encode("utf-8"),
  #      print_stdout=True)

           
  print(result)
  
  out = result["stdout"].decode("utf-8")
  
  print(result["stdout"].decode("utf-8"))
  
  with open("shelxc.log", "w") as output:
    output.write(out)                          

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
  
  keywords = fa  
  result = procrunner.run(["shelxd"],
                          stdin=keywords.encode("utf-8"),
                          print_stdout=False)
  return result

def simpleSHELXE(name, find, solvent_frac=0.5, inverse_hand=False):
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
  
  fa = name + '_fa'
  
  #keywords = keywords_shelxe.format(name, fa, solvent_frac, find)
  
  #cmd = 
  
  cmd = "shelxe {0} {1} -s{2} -m -h{3} -z{3} -a5 -q".format(name, fa, solvent_frac, find)
  if inverse_hand: cmd += " -i"

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
    
  result = procrunner.run(["shelxe"], stdin=cmd.encode("utf-8"), print_stdout=False)  
  return result



def copy_sca_locally(wavelengths):
  '''Copy .sca files locally to work around problem at DLS where SHELX fails
  to find the files'''
  for w in wavelengths:
    f = os.path.basename(w['sca'])
    shutil.copy(w['sca'], '.')
    w['sca'] = os.path.abspath(f)
  return wavelengths

if __name__ == '__main__':
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

  # Find xia2.json
  x2j = os.path.join(args.xia2dir, 'xia2.json')
  if not os.path.exists(x2j):
    raise RuntimeError('%s does not exist' % x2j)
  
  # Extract data from xia2.json
  x2_dat = xia2_json_reader(x2j)

  # Find scaled mtz
  scaled_mtz = x2_dat.scaled_mtz
  #print(3333, scaled_mtz)

  # Try to identify the wavelengths given the chosen scatterer
  x2_dat.identify_wavelengths(args.atom)
  print(x2_dat.identify_wavelengths(args.atom))

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
#    print("Num Met: ", num_met)
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

#  # copy .sca files locally and update the wavelengths dictionary
#  wl = copy_sca_locally(x2_dat.wavelengths)
#  print(wl)
  
  # determine spacegroups for pointgroup
  print(6666666666, x2_dat.sg_name)
  #gemmi_sg = gemmi.SpaceGroup(x2_dat.sg_name)
  #print(gemmi_sg)
  space_groups = generate(gemmi.SpaceGroup(x2_dat.sg_name))
  
  print(7777777, space_groups)

  cwd = os.path.normpath(os.getcwd())
  
  print(88888888, cwd)
  
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

#    print("********", os.getcwd())

    # Run SHELXC
 #   c_output.append(simpleSHELXC(args.name, x2_dat.cell, wl,
 #     sg_str, find, args.ntry))
      
 #   print(c_output)
    
    
    
    os.chdir(sg_str)
    
    print("********", os.getcwd())

    # copy .sca files locally and update the wavelengths dictionary
    wl = copy_sca_locally(x2_dat.wavelengths)
    print(wl)



    # Run SHELXC
    c_output.append(simpleSHELXC(args.name, x2_dat.cell, wl,
      sg_str, find, args.ntry))
      
    print(c_output)

    # Run SHELXD
    d_output = simpleSHELXD(args.name)
    
    # Get CFOM from .res
    try:
      with open(args.name + '_fa.res') as f:
        cfom.append(float(f.readline().split()[-1]))
    except IOError:
      cfom.append(0)

  os.chdir(cwd)

  best_idx = cfom.index(max(cfom))
  best_sg = space_groups[best_idx]
  print("Best space group: {0} with CFOM={1}").format(best_sg, cfom[best_idx])
  c_output = c_output[best_idx]
  for line in c_output:
    if line.startswith(' Resl'): print(line)
    if line.startswith(' <d"/sig>'): print(line)
    if line.startswith(' CC(1/2)'): print(line)
  
  os.chdir(best_sg)

  # Run SHELXE
  solvent_frac = 0.5 if matt_coeff is None else matt_coeff.solvent_fraction(
    matt_coeff.num_molecules())
  e_output_ori = simpleSHELXE(args.name, find, solvent_frac)
  #with Pool(processes = 2) as pool:
  #  e_output_ori = pool.map(simpleSHELXE(args.name, find, solvent_frac))

  #use print to write things to standard out and then it will be dumped to log file
  #print 'Solvent fraction: %s' %solvent_frac
  #cell_volume = matt_coeff.cell_volume()
  #print 'Unit cell volume: %s' %cell_volume
  #for line in e_output_ori:
  #  if line.startswith('Best trace'): print line
  try:
    with open(args.name + '.lst') as f:
      for line in f.readlines():
        if line.startswith('Best trace'): print(line)
  except IOError:
    pass
  e_output_inv = simpleSHELXE(args.name, find, solvent_frac, inverse_hand=True)
  #with Pool(processes = 2) as pool:
  #  e_output_inv = pool.map(simpleSHELXE(args.name, find, solvent_frac, inverse_hand=True))
  #for line in e_output_inv:
  #  if line.startswith('Best trace'): print line
  try:
    with open(args.name + '_i.lst') as f:
      for line in f.readlines():
        if line.startswith('Best trace'): print(line)
  except IOError:
    pass
