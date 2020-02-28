#!/usr/bin/env ccp4-python
# -*- coding: utf-8 -*-

from __future__ import division


class MtzData(object):
  '''Store items of interest from an MTZ file; this makes unit cell, space group
     and certain column types accessible to other programs'''
  
  def __init__(self, filename):
    from iotbx import mtz
    import os.path
    #check if an MTZ file has ben provided; if not raise error and request file
    if filename is None:
  	  raise RuntimeError('Need to specify hklin filename')
    elif not os.path.exists(filename):
      raise RuntimeError('%s does not exist' % filename)
    self.mtz_file = mtz.object(filename)
    
    self.sg_num = self.mtz_file.space_group_number()
    
    cols = self.mtz_file.column_labels()
    col_types = self.mtz_file.column_types()
    
    #pull out first column of type F; if there are multiple F and Q columns then
    #subsequent ones will be ignored; this may need amending for different phasing methods
    #getting column with structure factors'F'
    findex = col_types.index('F')
    self.F = cols[findex]
    #getting column with errors of structure factors 'Q'
    qindex = findex + 1
    assert col_types[qindex] == 'Q'
    self.Q = cols[qindex]
    
    self.cell = self._get_cell()
    
    return
    
  def _get_cell(self):
    '''Extract unit cell from the MTZ; compare different unit cells if multiple are
       present in one MTZ file and insure they are isomorphous; if they are not
       isomorphous raise an error'''
    
    from libtbx.test_utils import approx_equal
    xls = self.mtz_file.crystals()
    ucs = [e.unit_cell() for e in xls]
    
    cell0 = ucs[0]
    tst = all([cell0.is_similar_to(e) for e in ucs])
    
    if not tst:
      print "Multiple unit cells found! Only the first will be used:"
      for cell in ucs: print cell.parameters()
      
      #need to find a solution here so one can continue by picking a particular cell
    
    return cell0