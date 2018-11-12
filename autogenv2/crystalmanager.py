import numpy as np
from autogenv2.manager import resolve_status, update_attributes, Manager
from autogenv2.autorunner import RunnerPBS
from autogenv2.autopaths import paths
from qwalk_objects.crystal import CrystalReader
from qwalk_objects.propertiesreader import PropertiesReader
import os
import pickle as pkl
import shutil as sh

class CrystalManager(Manager):
  """ Internal class managing process of running a DFT job though crystal.
  Has authority over file names associated with this task.""" 
  def __init__(self,writer,runner,creader=None,name='crystal_run',path=None, preader=None,prunner=None,
      bundle=False,max_restarts=2):
    ''' CrystalManager manages the writing of a Crystal input file, it's running, and keeping track of the results.
    Args:
      writer (PySCFWriter): writer for input.
      reader (PySCFReader): to read PySCF output.
      runner (runner object): to run job.
      creader (CrystalReader): Reads the crystal results, (None implies use default reader).
      preader (PropertiesReader): Reads properties results, if any (None implies use default reader).
      prunner (runner object): run properties if needed (None implies use same runner as crystal).
      name (str): identifier for this job. This names the files associated with run.
      bundle (bool): Whether you'll use a bundling tool to run these jobs.
      max_restarts (int): maximum number of times you'll allow restarting before giving up (and manually intervening).
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

    #print(self.logname,": initializing")

    # Handle reader and runner defaults.
    self.writer=writer
    if creader is None: self.creader=CrystalReader()
    else: self.creader=creader
    if preader is None: self.preader=PropertiesReader()
    else: self.preader=preader
    if prunner is None: self.prunner=runner
    else: self.prunner=prunner
    if runner is None: self.runner=RunnerPBS()
    else: self.runner=runner

    # Internal.
    self.crysinpfn=self.name
    self.propinpfn=self.name+'.prop'
    self.crysoutfn=self.crysinpfn+'.o'
    self.propoutfn=self.propinpfn+'.o'
    self.restarts=0
    self.completed=False
    self.bundle=bundle

    # Smart error detection.
    self.max_restarts=max_restarts
    self.savebroy=[]

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      print(self.logname,": rebooting old manager.")
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
        skip_keys=['writer','runner','creader','preader','prunner','savebroy',
                   'path','logname','name',
                   'max_restarts','bundle'],
        take_keys=['restarts','completed','qwalk_orbs','qwalk_sys','bundle_ready','scriptfile'])

    # Update queue settings, but save queue information.
    update_attributes(copyto=self.runner,copyfrom=other.runner,
        skip_keys=['queue','walltime','np','nn','jobname','mode','account','prefix','postfix'],
        take_keys=['queueid'])
    update_attributes(copyto=self.prunner,copyfrom=other.prunner,
        skip_keys=['queue','walltime','np','nn','jobname','mode','account','prefix','postfix'],
        take_keys=['queueid'])

    update_attributes(copyto=self.creader,copyfrom=other.creader,
        skip_keys=[],
        take_keys=['completed','output'])

    update_attributes(copyto=self.preader,copyfrom=other.preader,
        skip_keys=[],
        take_keys=['completed','output'])

    updated=update_attributes(copyto=self.writer,copyfrom=other.writer,
        skip_keys=['maxcycle','edifftol'],
        take_keys=['completed','modisymm','restart','guess_fort','guess_fort13','_elements'])
    if updated:
      self.writer.completed=False

  #----------------------------------------
  def nextstep(self):
    ''' Determine and perform the next step in the calculation.'''
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    print(self.logname,": next step.")
    cwd=os.getcwd()
    os.chdir(self.path)

    # Generate input files.
    if not self.writer.completed:
      if self.writer.guess_fort is not None:
        sh.copy(self.writer.guess_fort,'fort.20')
      if self.writer.guess_fort13 is not None:
        sh.copy(self.writer.guess_fort,'in.fort.13') # save copy in case it's overwritten.
        sh.copy(self.writer.guess_fort,'fort.13')
      with open(self.crysinpfn,'w') as f:
        self.writer.write_crys_input(self.crysinpfn)
      with open(self.propinpfn,'w') as f:
        self.writer.write_prop_input(self.propinpfn)

    # Check on the CRYSTAL run
    status=resolve_status(self.runner,self.creader,self.crysoutfn)
    print(self.logname,": status= %s"%(status))

    if status=="not_started":
      self.runner.add_command("cp %s INPUT"%self.crysinpfn)
      self.runner.add_task("%s &> %s"%(paths['Pcrystal'],self.crysoutfn))

    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.creader.collect(self.crysoutfn)
      print(self.logname,": status %s"%status)
      if status=='killed':
        if self.restarts >= self.max_restarts:
          print(self.logname,": restarts exhausted (%d previous restarts). Human intervention required."%self.restarts)
        else:
          print(self.logname,": attempting restart (%d previous restarts)."%self.restarts)
          self.writer.restart=True
          sh.copy(self.crysinpfn,"%d.%s"%(self.restarts,self.crysinpfn))
          sh.copy(self.crysoutfn,"%d.%s"%(self.restarts,self.crysoutfn))
          sh.copy('fort.79',"%d.fort.79"%(self.restarts))
          self.writer.guess_fort='./fort.79'
          sh.copy(self.writer.guess_fort,'fort.20')
          self.writer.write_crys_input(self.crysinpfn)
          self.runner.add_command("cp %s INPUT"%self.crysinpfn)
          self.runner.add_task("%s &> %s"%(paths['Pcrystal'],self.crysoutfn))
          self.restarts+=1

    # Ready for bundler or else just submit the jobs as needed.
    if not self.bundle:
      qsubfile=self.runner.submit()

    self.completed=self.creader.completed

    # Update the file.
    with open(self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    os.chdir(cwd)

  #----------------------------------------
  def collect(self):
    ''' Call the collect routine for readers.'''
    print(self.logname,": collecting results.")
    self.creader.collect(self.path+self.crysoutfn)

    self.update_pickle()

  #------------------------------------------------
  def ready_properties(self):
    ''' Run properties for WF exporting.
    Returns:
      bool: whether it was successful.'''
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    ready=False
    self.nextstep()

    if not self.completed:
      return False

    cwd=os.getcwd()
    os.chdir(self.path)

    # Check on the properties run
    status=resolve_status(self.prunner,self.preader,self.propoutfn)
    print(self.logname,": properties status= %s"%(status))
    if status=='not_started':
      ready=False
      self.prunner.add_command("cp %s INPUT"%self.propinpfn)
      self.prunner.add_task("%s &> %s"%(paths['Pproperties'],self.propoutfn))

      if not self.bundle:
        qsubfile=self.prunner.submit()
    elif status=='ready_for_analysis':
      self.preader.collect(self.propoutfn)

    if self.preader.completed:
      ready=True
      print(self.logname,": properties completed successfully.")
    else:
      ready=False
      print(self.logname,": properties run incomplete.")

    os.chdir(cwd)
    self.update_pickle()

    return ready

  #----------------------------------------
  def export_record(self):
    ''' Combine input and results into convenient dict.'''
    res = {}
    spins_consistent = True
    if self.writer.spin_polarized:
      coarse_moments = coarsen_moments(self.creader.output['mag_moments'])
      if (coarse_moments != np.array(self.writer.initial_spins)).any():
        print("Spin moments changed from initial setting.")
        print("Inital setting: {}".format(np.array(self.writer.initial_spins)))
        print("Current spins:  {}".format(coarse_moments))
        spins_consistent = False
    res['spins_consistent'] = (spins_consistent)

    res['manager'] = self.__class__.__name__
    res['path'] = (self.path)
    res['name'] = (self.name)
    res['completed'] = (self.creader.completed)
    res['initial_spins'] = (self.writer.initial_spins)
    res['majority_guess'] = (self.writer.majority_guess)
    res['minority_guess'] = (self.writer.minority_guess)
    res['supercell'] = (self.writer.supercell)
    res['functional'] = (self.writer.functional)
    res['diis'] = (self.writer.diis)
    res['diis_opts'] = (self.writer.diis_opts)
    res['kmesh'] = (self.writer.kmesh)
    res['tolinteg'] = (self.writer.tolinteg)
    res['xml'] = (self.writer.xml_name)
    for prop in ['total_energy','mag_moments','atomic_charges']:
      if prop in self.creader.output:
        res[prop] = (self.creader.output[prop])
        #print(self.creader.output[prop])
        #print("  Found %s."%prop)
      else:
        res[prop] = (None)
        #print("  Didn't find %s."%prop)
    return res

def coarsen_moments(moments,cutoff_to_zero=0.5):
  moments = np.array(moments)
  coarse = np.zeros(moments.shape,dtype=int)
  coarse[moments > cutoff_to_zero] = 1
  coarse[moments < -cutoff_to_zero] = -1
  return coarse
