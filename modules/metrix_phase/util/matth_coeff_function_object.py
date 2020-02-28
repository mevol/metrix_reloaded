#!/bin/env python3

import sys
import procrunner
import re
import os
from mtz_data_object import MtzData
from seq_data_object import SeqData

class MattCoeff(object):
  '''get nmol in asu from cell dimensions and sequence; this is to get more information
     about the unknown target for which MR should be run'''

  def __init__(self, mw, cell, sg, hklin):    
  #get cell and symmetry from HKL file; get MW from sequence file
  #Prepare keywords for Matthew's coefficient as is required by CCP4
    
    keywords_matth = '''
                     MOLWEIGHT {0}
                     CELL {1} {2} {3} {4} {5} {6}
                     SYMMETRY {7}
                     AUTO
                     '''

    cell_n = tuple(re.findall(r'\d+(?:.\d+)?', str(cell)))
    sg_n = tuple(re.findall(r'"(.*?)"', str(sg)))
    fmt = (mw,) + cell_n + sg_n
    keywords = keywords_matth.format(*fmt)
    
    print(keywords)

    __location__ = os.path.realpath(os.path.join(os.getcwd(),
                                    os.path.dirname(__file__)))
                                    
    matt_coeff = os.path.join(__location__, "util/shell_scripts/matt_coeff.sh")


    self.matt_log = procrunner.run(matt_coeff,
                                   stdin=keywords.encode('utf-8'),
                                   print_stdout=False)

    self.extract_table()

  def extract_table(self):
    '''Extract the table of values from the log text'''

    in_table = False
    table = []
    for line in self.matt_log:
      if line.strip() == "_" * 48:
        in_table = not in_table
      elif in_table:
        table.append(line)

    self._table_lines = []
    for line in table:
      vals = [float(v) for v in line.split()]
      self._table_lines.append({
        'nmol_per_asu':vals[0],
        'coefficient':vals[1],
        'percent_solvent':vals[2],
        'probability':vals[3],
        })
    return

  def num_molecules(self):
    '''Extract most probable nmol from matthews_coef output; recognise pattern in
       output table to find line giving most likely number of molecules to search
       for during MR'''
    nmol = 0
    best_p = 0
    for line in self._table_lines:
      p = line['probability']
      if p > best_p:
        best_p = p
        nmol = int(line['nmol_per_asu'])
    return nmol

  def solvent_fraction(self, nmol):
    '''Extract estimated solvent fraction from matthews_coef output'''   
    for line in self._table_lines:
      if int(line['nmol_per_asu']) == nmol:
        return line['percent_solvent'] / 100    
    return None
  
def matt_coeff_factory(hklin, seqin):
  """Create a MattCoeff object given an MTZ file and a sequence file"""

  # get information from the sequence
  seqdata = SeqData(seqin)

  # load MTZ file and look-up header
  mtzdata = MtzData(hklin)
  
  # get number of molecules in asu
  matt_obj = MattCoeff(seqdata.mol_weight(),
                       mtzdata.cell,
                       mtzdata.sg_num,
                       hklin)
  
  return matt_obj

if __name__ == "__main__":
  matt_coeff_factory(sys.argv[1], sys.argv[2])
