import autogenv2
from autogenv2.autogen_tools import resolve_status, update_attributes
from autogenv2.autorunner import RunnerPBS
from autogenv2.autopaths import paths
import qwalk_objects as obj
import os
import pickle as pkl
import shutil as sh

class CrystalManager:
  """ Internal class managing process of running a DFT job though crystal.
  Has authority over file names associated with this task.""" 
  def __init__(self,writer,runner,creader=None,name='crystal_run',path=None, preader=None,prunner=None,
      trylev=False,bundle=False,max_restarts=2):
    ''' CrystalManager manages the writing of a Crystal input file, it's running, and keeping track of the results.
    Args:
      writer (PySCFWriter): writer for input.
      reader (PySCFReader): to read PySCF output.
      runner (runner object): to run job.
      creader (CrystalReader): Reads the crystal results, (None implies use default reader).
      preader (PropertiesReader): Reads properties results, if any (None implies use default reader).
      prunner (runner object): run properties if needed (None implies use same runner as crystal).
      name (str): identifier for this job. This names the files associated with run.
      trylev (bool): When restarting use LEVSHIFT option to encourage convergence, then do a rerun without LEVSHIFT.
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
    if creader is None: self.creader=obj.crystal.CrystalReader()
    else: self.creader=creader
    if preader is None: self.preader=obj.propertiesreader.PropertiesReader()
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
    self.trylev=trylev
    self.max_restarts=max_restarts
    self.savebroy=[]
    self.lev=False

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
        skip_keys=['writer','runner','creader','preader','prunner','lev','savebroy',
                   'path','logname','name',
                   'trylev','max_restarts','bundle'],
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
          if self.trylev:
            print(self.logname,": trying LEVSHIFT.")
            self.writer.levshift=[10,1] # No mercy.
            self.savebroy=deepcopy(self.writer.broyden)
            self.writer.broyden=[]
            self.lev=True
          sh.copy(self.crysinpfn,"%d.%s"%(self.restarts,self.crysinpfn))
          sh.copy(self.crysoutfn,"%d.%s"%(self.restarts,self.crysoutfn))
          sh.copy('fort.79',"%d.fort.79"%(self.restarts))
          self.writer.guess_fort='./fort.79'
          sh.copy(self.writer.guess_fort,'fort.20')
          self.writer.write_crys_input(self.crysinpfn)
          self.runner.add_command("cp %s INPUT"%self.crysinpfn)
          self.runner.add_task("%s &> %s"%(paths['Pcrystal'],self.crysoutfn))
          self.restarts+=1
    elif status=='done' and self.lev:
      # We used levshift to converge. Now let's restart to be sure.
      print("Recovering from LEVSHIFTer.")
      self.writer.restart=True
      self.writer.levshift=[]
      self.creader.completed=False
      self.lev=False
      sh.copy(self.crysinpfn,"%d.%s"%(self.restarts,self.crysinpfn))
      sh.copy(self.crysoutfn,"%d.%s"%(self.restarts,self.crysoutfn))
      sh.copy('fort.79',"%d.fort.79"%(self.restarts))
      self.writer.guess_fort='./fort.79'
      sh.copy(self.writer.guess_fort,'fort.20')
      self.writer.write_crys_input(self.crysinpfn)
      sh.copy(self.crysinpfn,'INPUT')
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

  #----------------------------------------
  def submit(self):
    ''' Submit any work and update the manager.'''
    qsubfile=self.runner.submit()

    self.update_pickle()

    return qsubfile

  #----------------------------------------
  def release_commands(self):
    ''' Release the runner of any commands it was tasked with and update the manager.'''
    commands=self.runner.release_commands()
    self.update_pickle()

    return commands

  #------------------------------------------------
  def update_queueid(self,qid):
    ''' If a bundler handles the submission, it can update the queue info with this.
    Args:
      qid (str): new queue id from submitting a job. The Manager will check if this is running.
    '''
    self.runner.queueid.append(qid)
    self.update_pickle()

  #------------------------------------------------
  def update_pickle(self):
    ''' If you make direct changes to the internals of the pickle, you need to call this to insure they are saved.'''
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

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
  def status(self):
    if self.completed:
      return 'ok'
    else:
      return 'not_finished'
