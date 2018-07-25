#!/usr/bin/env python3
''' Simple utilities for interacting with autogen runs.'''

import pickle as pkl
import sys
import argparse

def get_info(pickle):
  print("Info about %s..."%pickle)
  man=pkl.load(open(pickle,'rb'))

  print("  Queue id: {}".format(man.runner.queueid))
  try:
    print("  Properties queue id: {}".format(man.prunner.queueid))
  except AttributeError:
    pass

def set_attribute(pickle,attr,val):
  ''' Set an attribute in a pickled manager and save the updated manager back to that pickle.
  Args:
    pickle (str): path to pickled manager.
    attr (str): name of attribute to change. 
    val: value of attribute to set.
  '''
  print("Setting %s in %s..."%(attr,pickle))
  man=pkl.load(open(pickle,'rb'))
  man.__dict__[attr]=val
  pkl.dump(man,open(pickle,'wb'))

if __name__=='__main__':

  parser=argparse.ArgumentParser("Autogen untilities.")
  parser.add_argument('manager',type=str,help='Pickle file to look at.')
  # Can add more options as needed.

  args=parser.parse_args()
  get_info(args.manager)
  
