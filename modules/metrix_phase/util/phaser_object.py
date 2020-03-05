#!/bin/env python3

import procrunner
import os

class PhaserSearch(object):
  '''simple space group search with single search model; using MTZ file for data and
     therein structure factor columns; use PDB file as search model; give sequence
     identity; give number of molecules to look for; give sequence of target structure'''

  def __init__(self, xyzin, hklin, seqin, f, sigf,
                            num, seq_identity, mode):
    '''Produce keywords for Phaser'''
    
    keywords_phaser = '''
      MODE {0}
      SGALTERNATIVE SELECT ALL
      HKLIN {1}
      LABIN F = {2} SIGF = {3}
      ENSEMBLE ens PDB {4} IDENTITY {7} 
      COMPOSITION PROTEIN SEQUENCE {5} NUM {6}
      SEARCH ENSEMBLE ens
      '''
          
    keywords = keywords_phaser.format(mode, hklin, f, sigf, xyzin, seqin,
                                num, seq_identity)
                                
    __location__ = os.path.realpath(os.path.join(os.getcwd(),
                                    os.path.dirname(__file__)))
                                    
    phaser = os.path.join(__location__, "shell_scripts/phaser.sh")
    
    self.phaser_log = procrunner.run(["/bin/bash", phaser],
                                     stdin=keywords.encode("utf-8"))
                                     
    out = self.phaser_log["stdout"].decode("utf-8")
    
    with open("phaser.log", "w") as output:
      output.write(out)                          
  

