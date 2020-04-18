#!/bin/env python3

import json
import sys
from math import ceil
import procrunner
import gemmi
import os
import xraydb

class xia2_json_reader(object):
  def __init__(self, json_path):
    try:
      with open(json_path) as f:
        self.json_string = f.read()
    except IOError:
      print('Please provide a valid xia2.json')
      raise

    self.json = json.loads(self.json_string)

    # try to extract a crystal
    if len(self.json['_crystals']) != 1:
      raise TypeError('Expected 1 crystal, found more.'
                      ' Please check xia2.json')
    
    
    self.crystals = self.json['_crystals']
    
    
    for name in self.crystals.keys():
      # extract scaler cell
      self.cell = self.crystals[name]['_scaler']['_scalr_cell']
      # extract space group info
      self.sg_name = self.crystals[name]['_scaler']['_scalr_likely_spacegroups'][0]
      # extract name and path of scaled MTZ file
      self.scaled_mtz = str(self.crystals[name]['_scaler']['_scalr_scaled_reflection_files']['mtz'])

      # extract wavelength details
      self.wavelengths = [{'name':k,
                         'wavelength':v['_wavelength']} \
        for k, v in self.crystals[name]['_wavelengths'].items()]
      self.wavelengths = sorted(self.wavelengths, key=lambda d: d['wavelength'])

      # include paths to unmerged sca files
      sca = self.crystals[name]['_scaler']['_scalr_scaled_reflection_files']['sca']
      for w in self.wavelengths:
        w['sca'] = sca[w['name']]

      return

  def identify_wavelengths(self, scatterer="Se"):
    """Try to identify peak, infl, hrem and lrem"""

    # get label for single wavelength --> SAD
    if len(self.wavelengths) == 1:
      self.wavelengths[0]['label'] = 'SAD'
      return

    # ignore more than 4 wavelengths MAD
    if len(self.wavelengths) >= 4:
      raise RuntimeError("Can't handle > 4 wavelength MAD today :-(")

    # get f' and f'' for 3-wavelength MAD and asign labels    
    min_wl = self.wavelengths[0]['wavelength']
    max_wl = self.wavelengths[-1]['wavelength']
    wl_range =  max_wl - min_wl
    npoints = int(ceil(wl_range / 0.000001))
    lambdas = [min_wl + 0.000001 * e for e in range(npoints + 1)]

    element = gemmi.Element(scatterer).atomic_number
    fp = [(gemmi.cromer_libermann(z=int(element), energy=12398.4197386209/float(l)))[0] \
      for l in lambdas]
    fdp = [(gemmi.cromer_libermann(z=int(element), energy=12398.4197386209/float(l)))[1] \
      for l in lambdas]
    infl_idx = fp.index(min(fp))
    peak_idx = fdp.index(max(fdp))

    theor_peak_wl = lambdas[peak_idx]
    infl_wl = lambdas[infl_idx]
    
    close_to_peak = [d for d in self.wavelengths \
      if abs(d['wavelength'] - theor_peak_wl) < 0.001]
    if len(close_to_peak) == 1:
      close_to_peak[0]['label'] = 'peak'
      peak_wl = close_to_peak[0]['wavelength']
      #print("single peak", peak_wl)
    elif len(close_to_peak) == 2:
      close_to_peak[1]['label'] = 'infl'
      close_to_peak[0]['label'] = 'peak'
      peak_wl = close_to_peak[0]['wavelength']
      #print("double peak", peak_wl)
    else:
      raise RuntimeError('Could not distinguish peak and inflection '
                         'point wavelengths')

    for wl in self.wavelengths:
      if 'label' in self.wavelengths:
        continue
      if wl['wavelength'] < peak_wl - 0.008:
        wl['label'] = 'hrem'
      elif wl['wavelength'] > peak_wl + 0.008:
        wl['label'] = 'lrem'

    return self.wavelengths

if __name__ == '__main__':
  if len(sys.argv) < 2:
    print('Please provide a xia2.json')
    sys.exit()

  xia2_json = xia2_json_reader(sys.argv[1])
  xia2_json.identify_wavelengths()
