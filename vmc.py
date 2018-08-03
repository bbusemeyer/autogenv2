from __future__ import print_function
import os
import average_tools as avg
####################################################
class VMCWriter:
  def __init__(self,options={}):
    ''' Object for producing input into a VMC QWalk run. 
    Args:
      options (dict): editable options are as follows.
        trialfunc (str): system and trial wavefunction section.
        errtol (float): tolerance for the estimated energy error. 
        extra_observables (list): see `average_tools.py` for how to use this.
        minblocks (int): minimum number of VMC steps to take, considering equillibration time.
        iterations (int): number of VMC steps to attempt.
    '''
    self.trialfunc=''
    self.errtol=0.1
    self.minblocks=10 
    self.nblock=100
    self.extra_observables=[]

    self.qmc_abr='VMC'
    self.completed=False
    self.set_options(options)

  #-----------------------------------------------
    
  def set_options(self, d):
    ''' Save setting of options.
    Args: 
      d (dict): attributes to update.
    '''
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        raise ValueError("Error:",k,"not a keyword for VMCWriter")
      selfdict[k]=d[k]

    # Check completeness of average generator options.
    for avg_generator in self.extra_observables:
      avg.check_opts(avg_generator)

  #-----------------------------------------------
  def qwalk_input(self,infile):
    if self.trialfunc=='':
      print(self.__class__.__name__,": Trial function not ready. Postponing input file generation.")
      self.completed=False
    else:
      outlines=[
          "method { VMC timestep %g nblock %i"%(self.timestep,self.nblock)
        ]
      for avg_opts in self.extra_observables:
        outlines+=avg.average_section(avg_opts)
      outlines+=["}"]
      outlines+=self.trialfunc.split('\n')

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
