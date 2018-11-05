import subprocess as sub
import importlib
import os
import shutil
import sys
import time

#####################################################################################
class LocalSubmitter:
  """Abstract submission class. Defines interaction with queuing system (for
  nonqueuing defines how to execute jobs). Child classes must define:
  __init__: 
     can have any parameters, but should set up any variables like queue time and number of
     processor cores
  _submit_job(self,
              inpfns : list of input filenames (list of strings)
              outfn  : where stdout should be sent (string)
              jobname : job name for the queue (string)
              loc :  directory where job should be run (string)
            )
            returns a list of queue ids (list of strings)
  """
#-------------------------------------------------------

def check_BW_stati(queueids):
  """Utility function to determine the status of a set Blue Waters job.
  Args: 
    queueids (list): list of queueids as string representation of int, e.g. ['4819103','4819104'].
  """
  try:
    qstat = sub.check_output(
        "qstat ", stderr=sub.STDOUT, shell=True
      ).decode()
  except sub.CalledProcessError:
    return "unknown"
  qstat=qstat.split('\n')
  for qid in queueids:
    for line in qstat:
      spl=line.split()
      if qid in line and len(spl) > 4:
        stat=line.split()[-2]
        if stat == "R" or stat == "Q":
          return "running"
  return 'unknown'

#-------------------------------------------------------
def check_PBS_stati(queueids):
  """Utility function to determine the status of a set PBS job.
  Args: 
    queueids (list): list of queueids as string representation of int, e.g. ['4819103','4819104'].
  """
  try:
    qstat = sub.check_output(
        "qstat ", stderr=sub.STDOUT, shell=True
      ).decode()
  except sub.CalledProcessError:
    return "unknown"
  qstat=qstat.split('\n')
  for qid in queueids:
    for line in qstat:
      spl=line.split()
      if qid in line and len(spl) > 4:
        stat=line.split()[4]
        if stat == "R" or stat == "Q":
          return "running"
  return 'unknown'

