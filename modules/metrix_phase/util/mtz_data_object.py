#!/bin/env python3

import os.path
import gemmi
import re

class MtzData(object):
  '''Store items of interest from an MTZ file; this makes unit cell, space group
     and certain column types accessible to other programs'''
  
  def __init__(self, filename):
    #check if an MTZ file has ben provided; if not raise error and request file
    if filename is None:
  	  raise RuntimeError('Need to specify hklin filename')
    elif not os.path.exists(filename):
      raise RuntimeError('%s does not exist' % filename)
    self.mtz_file = gemmi.read_mtz_file(filename)
    
    self.sg_num = self.mtz_file.spacegroup.number
    
    f = self.mtz_file.columns_with_type("F")
    f = re.findall(r'[A-Z]+', str(f))
    if "F" and "WAVE" in f:
      f = "F_WAVE1"
      self.fcols = f
    else:
      f = "F"
      self.fcols = f

    q = self.mtz_file.columns_with_type("Q")
    q = re.findall(r'[A-Z]+', str(q))
    if "SIGF" and "WAVE" in q:
      q = "SIGF_WAVE1"
      self.qcols = q
    else:
      q = "SIGF"
      self.qcols = q

    self.cell = self.mtz_file.cell
    
    return

