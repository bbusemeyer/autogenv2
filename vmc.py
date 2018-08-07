from __future__ import print_function
import os
import average_tools as avg
from trialfunc import export_qwalk_trialfunc

####################################################
class VMCWriter:
  def __init__(self,
        system,
        trialfunc,
        errtol=0.1,
        minblocks=10,
        nblock=100,
        optimize_det=False,
        averages=''
      ):
    ''' Object for producing input into a VMC QWalk run. 
    Args:
      options (dict): editable options are as follows.
        system (System or str): system object or system section.
        trialfunc (TrialFunc or str): TrialFunc object or trial wavefunction section.
        errtol (float): tolerance for the estimated energy error. 
        averages (str): averages section in qwalk.
        optimize_det (bool): treat coefficients of slater determinant as parameters (for derivatives).
        minblocks (int): minimum number of VMC steps to take, considering equillibration time.
    '''
    self.system=system
    self.trialfunc=trialfunc
    self.errtol=errtol
    self.minblocks=minblocks
    self.nblock=nblock
    self.averages=averages
    self.optimize_det=optimize_det

    self.qmc_abr='VMC'
    self.completed=False

  #-----------------------------------------------
  def qwalk_input(self,infile):
    if type(self.system) is not str:
      system=self.system.export_qwalk_sys()
    else:
      system=self.system
    if type(self.trialfunc) is not str:
      trialfunc=export_qwalk_trialfunc(self.trialfunc,optimize_det=self.optimize_det)
    else:
      trialfunc=self.trialfunc

    outlines=[
        "method { VMC",
        "  nblock %i"%(self.nblock)
      ]+['  '+line for line in self.averages.split('\n')]+[
        "}"
      ]
    outlines.append(system)
    outlines.append(trialfunc)

    with open(infile,'w') as f:
      f.write('\n'.join(outlines))

    self.completed=True

     
####################################################
import subprocess as sub
import json
class VMCReader:
  ''' Reads results from a VMC calculation. 

  Attributes:
    output (dict): results of calculation. 
    completed (bool): whether the run has converged to a final answer.
  '''
  def __init__(self,errtol=0.01,minblocks=15):
    self.output={}
    self.completed=False

    self.errtol=errtol
    self.minblocks=minblocks
    self.gosling="gosling"

  def read_outputfile(self,outfile):
    ''' Read output file results.

    Args:
      outfile (str): output to read.
    '''
    return json.loads(sub.check_output([self.gosling,"-json",outfile.replace('.o','.log')]).decode())

  def check_complete(self):
    ''' Check if a VMC run is complete.
    Returns:
      bool: If self.results are within error tolerances.
    '''
    completed=True
    if len(self.output)==0:
      return False # No results yet.
    if self.output['properties']['total_energy']['error'][0] > self.errtol:
      print("VMC incomplete: (%f) does not meet tolerance (%f)"%\
          (self.output['properties']['total_energy']['error'][0],self.errtol))
      completed=False
    if self.output['total blocks']-self.output['warmup blocks'] < self.minblocks:
      print("VMC incomplete: Run completed %d blocks, but requires %d."%\
          (self.output['total blocks']-self.output['warmup blocks'],self.minblocks))
      completed=False
    return completed
          
  #------------------------------------------------
  def collect(self,outfile):
    ''' Collect results for an output file and resolve if the run needs to be resumed. 

    Args: 
      outfiles (list): list of output file names to open and read.
    Returns:
      str: status of run = {'ok','restart'}
    '''
    # Gather output from files.
    self.completed=True
    status='unknown'
    if os.path.exists(outfile):
      self.output=self.read_outputfile(outfile)
      self.output['file']=outfile

    # Check files.
    self.completed=self.check_complete()
    if not self.completed:
      status='restart'
    else:
      status='ok'
    return status
      
  #------------------------------------------------
  def write_summary(self):
    ''' Print out all the items in output. '''
    print("#### Variational Monte Carlo")
    for f,out in self.output.items():
      print(f,out)
