from autogenv2.manager import resolve_status, update_attributes, Manager
from autogenv2.autorunner import RunnerPBS
from autogenv2.autopaths import paths
from qwalk_objects.trialfunc import export_qwalk_trialfunc,separate_jastrow,Jastrow
import os
import pickle as pkl

#######################################################################
class QWalkManager(Manager):
  def __init__(self,writer,reader,runner=None,trialfunc=None,
      name='qw_run',path=None,bundle=False):
    ''' QWalkManager managers the writing of a QWalk input files, it's running, and keeping track of the results.
    Args:
      writer (qwalk writer): writer for input.
      reader (qwalk reader): to read job.
      runner (Runner object): to run job.
      trialfunc (TrialFunction): TrialFunction object for generating trail function input. 
        Note: This is only used if write.trailfunc arguement==''. 
      name (str): identifier for this job. This names the files associated with run.
      path (str): directory where this manager is free to store information.
      bundle (bool): False - submit jobs. True - dump job commands into a script for a bundler to run.
      qwalk (str): absolute path to qwalk executible.
    '''
    self.name=name
    self.pickle="%s.pkl"%(self.name)

    # Ensure path is set up correctly.
    if path is None:
      path=os.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    self.logname="%s@%s"%(self.__class__.__name__,self.path+self.name)

    #print(self.logname,": initializing")

    self.writer=writer
    self.reader=reader
    self.trialfunc=trialfunc
    if runner is not None: self.runner=runner
    else: self.runner=RunnerPBS()
    self.bundle=bundle

    self.completed=False
    self.infile=name
    self.outfile="%s.o"%self.infile
    self.stdout="%s.out"%self.infile

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

    #TODO this forbids all changes to trialfunc's managers even their runners (for instance). Should allows safe changes.
    update_attributes(copyto=self,copyfrom=other,
        skip_keys=['writer','runner','reader','path','logname','name','bundle'],
        take_keys=['restarts','completed'])

    # Update queue settings, but save queue information.
    update_attributes(copyto=self.runner,copyfrom=other.runner,
        skip_keys=['queue','walltime','np','nn','jobname'],
        take_keys=['queueid'])

    update_attributes(copyto=self.reader,copyfrom=other.reader,
        skip_keys=['errtol','minblocks'],
        take_keys=['completed','output'])

    updated=update_attributes(copyto=self.writer,copyfrom=other.writer,
        skip_keys=['maxcycle','nblock','savetrace'],
        take_keys=['completed','tmoves','extra_observables','timestep','sys','trialfunc'])
    if updated:
      self.writer.completed=False

  #------------------------------------------------
  def nextstep(self):
    ''' Perform next step in calculation. trialfunc managers are updated if they aren't completed yet.'''
    # Recover old data.
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    print(self.logname,": next step.")

    # Check dependency is completed first.
    if self.writer.trialfunc=='':
      print(self.logname,": checking trial function.")
      self.writer.trialfunc = export_qwalk_trialfunc(self.trialfunc)

    # Work on this job.
    cwd=os.getcwd()
    os.chdir(self.path)

    # Write the input file.
    if not self.writer.completed:
      self.writer.qwalk_input(self.infile)
    
    status=resolve_status(self.runner,self.reader,self.outfile)
    print(self.logname,": %s status= %s"%(self.name,status))
    if status=="not_started" and self.writer.completed:
      exestr="%s %s &> %s"%(paths['qwalk'],self.infile,self.stdout)
      self.runner.add_task(exestr)
      print(self.logname,": %s status= submitted"%(self.name))
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.reader.collect(self.outfile)
      if status=='ok':
        print(self.logname,": %s status= %s, task complete."%(self.name,status))
        self.completed=True
      else:
        print(self.logname,": %s status= %s, attempting rerun."%(self.name,status))
        exestr="%s %s &> %s"%(paths['qwalk'],self.infile,self.stdout)
        self.runner.add_task(exestr)
    elif status=='done':
      self.completed=True

    # Ready for bundler or else just submit the jobs as needed.
    if not self.bundle:
      qsubfile=self.runner.submit()

    # Update the file.
    with open(self.pickle,'wb') as outf:
      pkl.dump(self,outf)

    os.chdir(cwd)

  #----------------------------------------
  def collect(self):
    ''' Call the collect routine for readers.'''
    print(self.logname,": collecting results.")
    self.reader.collect(self.path+self.outfile)

    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #----------------------------------------
  def export_jastrow(self,optimizebasis=True,freezeall=False):
    ''' Make a Jastrow function from any resulting trial function optimization.
    Returns:
      bool: Whether it was successful.'''
    # Theoretically more than just Jastrow can be provided, but practically that's the only type of wavefunction we tend to export.

    # Recover old data.
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    wfout = self.path+self.outfile.replace('.o','.wfout')
    if not self.completed or not os.path.exists(wfout):
      raise AssertionError("QWalk run not ready, or doesn't make a Jastrow.")

    return Jastrow(separate_jastrow(wfout,optimizebasis=optimizebasis,freezeall=freezeall))
  
  def export_record(self,obdmfunc=None,tbdmfunc=None,obdmerrfunc=None,tbdmerrfunc=None):
    ''' Combine input and output into convenient run record.'''
    if obdmfunc is None: obdmfunc       = lambda x: x
    if obdmerrfunc is None: obdmerrfunc = lambda x: x
    if tbdmfunc is None: tbdmfunc       = lambda x: x
    if tbdmerrfunc is None: tbdmerrfunc = lambda x: x

    res = {}
    res['path'] = self.path
    res['name'] = self.name
    res['manager'] = self.__class__.__name__
    res['completed'] = self.completed

    # This presupposes slater, but maybe that's ok for now.
    res['states'] = self.writer.trialfunc.slater.states
    res['coefs'] = self.writer.trialfunc.slater.weights
    res['kweight'] = self.writer.trialfunc.slater.orbitals.kweight
    res['kpoint'] = self.writer.trialfunc.slater.orbitals.kpoint

    # Currently there are different formats for different QMC runs.
    # TODO unify QMC output formats.
    if 'total_energy' in self.reader.output:
      res['total_energy'] = self.reader.output['total_energy']
      res['total_energy_err'] = self.reader.output['total_energy_err']

    if 'properties' in self.reader.output:
      if 'total_energy' in self.reader.output['properties']:
        res['total_energy'] = self.reader.output['properties']['total_energy']['value'][0]
        res['total_energy_err'] = self.reader.output['properties']['total_energy']['error'][0]
      if 'tbdm_basis' in self.reader.output['properties']:
        res['basis'] = self.reader.output['properties']['tbdm_basis']['states']
        if 'tbdm' in self.reader.output['properties']['tbdm_basis']:
          res['tbdm'] = tbdmfunc([[self.reader.output['properties']['tbdm_basis']['tbdm'][spini+spinj] for spinj in ('up','down')] for spini in ('up','down')])
          res['tbdm_err'] = tbdmerrfunc([[self.reader.output['properties']['tbdm_basis']['tbdm'][spini+spinj+'_err'] for spinj in ('up','down')] for spini in ('up','down')])
        if 'obdm' in self.reader.output['properties']['tbdm_basis']:
          res['obdm'] = obdmfunc([self.reader.output['properties']['tbdm_basis']['obdm'][spin] for spin in ('up','down')])
          res['obdm_err'] = obdmerrfunc([self.reader.output['properties']['tbdm_basis']['obdm'][spin+'_err'] for spin in ('up','down')])

    return res
