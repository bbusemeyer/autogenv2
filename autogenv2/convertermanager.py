import autogenv2
from autogenv2.autogen_tools import resolve_status, update_attributes
from autogenv2.autopaths import paths
import qwalk_objects as obj
import os
import pickle as pkl
import shutil as sh

class ConverterManager:
  """ Internal class managing process of running a DFT job though crystal.
  Has authority over file names associated with this task.""" 
  def __init__(self,realonly=True,maxbands=(None,None),name='converter',path=None):
    ''' CrystalManager manages the writing of a Crystal input file, it's running, and keeping track of the results.
      Args:
        name (str): identifier for this job. This names the files associated with run.
        path (str): where to operate the jobs in.
    '''
    # Where to save self.
    self.name=name
    self.pickle="%s.pkl"%self.name

    # Ensure path is set up correctly.
    if path is None:
      path=os.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    self.logname="%s@%s"%(self.__class__.__name__,self.path+self.name)

    # Internal variables.
    self.system = None
    self.orbitals = None

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      #print(self.logname,": rebooting old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      self.recover(old)

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #------------------------------------------------
  def recover(self,other):
    ''' Recover old class by copying over data. Retain variables from old that may change final answer.'''
    # Practically speaking, the run will preserve old `take_keys` and allow new changes to `skip_keys`.
    # This is because you are taking the attributes from the older instance, and copying into the new instance.

    update_attributes(copyto=self,copyfrom=other,
        skip_keys=['path','logname','name'],
        take_keys=['restarts','completed','qwalk_orbs','qwalk_sys','bundle_ready','scriptfile'])

  #----------------------------------------
  def nextstep(self):
    ''' Determine and perform the next step in the calculation.'''
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    print(self.logname,": next step.")
    cwd=os.getcwd()
    os.chdir(self.path)

    if self.system is None or self.orbitals is None:
      print(self.logname,": converting solutions to qwalk.")
      self.system, self.orbitals = obj.pack_objects(maxbands=self.maxbands,realonly=self.realonly)

    # Update the file.
    with open(self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    os.chdir(cwd)

  #----------------------------------------
  def status(self):
    if self.completed:
      return 'ok'
    else:
      return 'not_finished'

