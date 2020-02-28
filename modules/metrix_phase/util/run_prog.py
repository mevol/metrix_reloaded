#!/usr/bin/env ccp4-python
# -*- coding: utf-8 -*-

from __future__ import division



def run_prog(prog_name, cmd="", keywords="", show_stdout=False):
  try:
    from CCP4Dispatchers import dispatcher_builder
  except ImportError:
    import os, sys
    sys.path.insert(0,
      os.path.join(os.environ["CCP4"],"share","python"))
  from CCP4Dispatchers import dispatcher_builder  
  d = dispatcher_builder(prog_name, cmd, keywords)
  print "Running: {0}".format(prog_name + " " + cmd)
  d.call()
  if d.call_err:
    print "Error calling {0}:".format(prog_name), d.call_err
    
  if d.call_val !=0:
    print "Unusual return value from {0}:".format(prog_name), d.call_val
    
  if d.stderr_data:
    print "{0} error stream output".format(prog_name)
    for line in d.stderr_data: print line
    
  if d.stdout_data and show_stdout:
    print "{0} output".format(prog_name)
    for line in d.stdout_data: print line
  
  
  return d.stdout_data