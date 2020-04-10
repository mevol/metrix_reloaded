#!/bin/env python3

import json
import sys
from math import ceil
import procrunner
import gemmi
import os
#from cctbx.eltbx import henke as table_module
#from iotbx import mtz

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
 #     print(self.cell)
      # extract space group info
      self.sg_name = self.crystals[name]['_scaler']['_scalr_likely_spacegroups'][0]
#      print(self.sg_name)
      

      self.scaled_mtz = str(self.crystals[name]['_scaler']['_scalr_scaled_reflection_files']['mtz'])
      #print("My MTZ", self.scaled_mtz)
#    
#      mtz = gemmi.read_mtz_file(self.scaled_mtz)    
#    mtz_obj = mtz.object(self.scaled_mtz)
#    self.sg_name = mtz_obj.space_group_name()
#    mtz = gemmi.read_mtz_file()

      # extract wavelength details
      self.wavelengths = [{'name':k,
                         'wavelength':v['_wavelength']} \
        for k, v in self.crystals[name]['_wavelengths'].items()]
      self.wavelengths = sorted(self.wavelengths, key=lambda d: d['wavelength'])
      print("xia2 wavelengths", self.wavelengths)

      # include paths to unmerged sca files
      sca = self.crystals[name]['_scaler']['_scalr_scaled_reflection_files']['sca']
      for w in self.wavelengths:
        w['sca'] = sca[w['name']]

      return

  def identify_wavelengths(self, scatterer="Se"):
    """Try to identify peak, infl, hrem and lrem"""

    #t=table_module.table(scatterer)

    if len(self.wavelengths) == 1:
      self.wavelengths[0]['label'] = 'SAD'
      return

    if len(self.wavelengths) > 4:
      raise RuntimeError("Can't handle > 4 wavelength MAD today :-(")

    # create range of lambdas spanning all datasets, incrementing
    # by 0.00001 A
    min_wl = self.wavelengths[0]['wavelength']
    max_wl = self.wavelengths[-1]['wavelength']
    wl_range =  max_wl - min_wl
    npoints = int(ceil(wl_range / 0.000001))
    lambdas = [min_wl + 0.000001 * e for e in range(npoints + 1)]

    # get theoretical f' and f''
    
    for l in lambdas:
      keywords_crossec = "\n".join([
                                    f"ATOM {scatterer}",
                                    f"NWAV 1 {l}",
                                    f"END",
                                   ])
    
      __location__ = os.path.realpath(os.path.join(os.getcwd(),
                                      os.path.dirname(__file__)))
      crossec_shell = os.path.join(__location__,
                                   "shell_scripts/crossec.sh")

      print(0000000000, crossec_shell)

      result = procrunner.run(["crossec"],
                              stdin=keywords_crossec.encode('utf-8'),
                              print_stdout=False)

      print(result)
      #keywords_matth = '''
      #                 ATOM {0}
      #                 NWAV {1} {2}
      #                 END
      #                 '''
      #fmt = (scatterer,) + lambdas + (l,)
      
    
    fp = [t.at_angstrom(l).fp() for l in lambdas]
    fdp = [t.at_angstrom(l).fdp() for l in lambdas]
    infl_idx = fp.index(min(fp))
    peak_idx = fdp.index(max(fdp))

    theor_peak_wl = lambdas[peak_idx]
    infl_wl = lambdas[infl_idx]

    # look for wavelengths within 0.001 of the theoretical peak
    close_to_peak = [d for d in self.wavelengths \
      if abs(d['wavelength'] - theor_peak_wl) < 0.001]
    if len(close_to_peak) == 1:
      close_to_peak[0]['label'] = 'peak'
      peak_wl = close_to_peak[0]['wavelength']
    elif len(close_to_peak) == 2:
      close_to_peak[1]['label'] = 'infl'
      close_to_peak[0]['label'] = 'peak'
      peak_wl = close_to_peak[0]['wavelength']
    else:
      raise RuntimeError('Could not distinguish peak and inflection '
                         'point wavelengths')

    for wl in self.wavelengths:
      if wl.has_key('label'): continue
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
