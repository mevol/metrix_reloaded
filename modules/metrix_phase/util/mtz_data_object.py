#!/bin/env python3

import os.path
import gemmi

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
    
    self.sg_num = self.mtz_file.spacegroup
    fcols = self.mtz_file.columns_with_type("F")
    qcols = self.mtz_file.columns_with_type("Q")[0]
    self.cell = self.mtz_file.cell
    
    return

