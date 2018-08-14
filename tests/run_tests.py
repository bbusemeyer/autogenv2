''' 
Runs a set of tests to check new updates are not breaking any features.
This is also a good demo of advanced features.
'''

from crystal import CrystalWriter
from crystalmanager import CrystalManager
from qwalkmanager import QWalkManager
from variance import VarianceWriter,VarianceReader
from linear import LinearWriter,LinearReader
from dmc import DMCWriter,DMCReader
from trialfunc import SlaterJastrow
from autorunner import RunnerPBS
<<<<<<< HEAD
import numpy as np

###################################################################################################################
# Individual tests.
=======
from pickle import load
import numpy as np

# TODO 
# Dmc jobs aren't saving their queueid and are getting run multiple times, I think.

###################################################################################################################
# Individual test definitions.
>>>>>>> master
def h2_crystal_equil_test():
  jobs=[]

  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'spin_polarized':False
    })
  cwriter.set_struct_fromxyz(open('structures/h2_equil.xyz','r').read())

  cman=CrystalManager(
      name='crys',
      path='test_h2eq_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=16,
          walltime='0:01:00'
        )
    )
  jobs.append(cman)

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:05:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  jobs.append(var)

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs.append(lin)

  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='0:10:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)

  return jobs

def h2_crystal_stretch_test():
  jobs=[]

  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'spin_polarized':True,
      'initial_spins':[1,-1]
    })
  cwriter.set_struct_fromxyz(open('structures/h2_stretch.xyz','r').read())

  cman=CrystalManager(
      name='crys',
      path='test_h2st_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=16,
          walltime='0:01:00'
        )
    )
  jobs.append(cman)

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:05:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  jobs.append(var)

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs.append(lin)

  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='0:10:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)

  return jobs

def si_crystal_test():
  ''' Simple tests that check PBC is working Crystal, and that QMC can be performed on the result.'''
  jobs=[]

  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'kmesh':(3,3,3),
      'spin_polarized':False
    })
  cwriter.set_struct_fromcif(open('structures/si.cif','r').read(),primitive=True)

  cman=CrystalManager(
      name='crys',
      path='test_si_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=16,
          walltime='0:10:00'
        )
    )
  jobs.append(cman)

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  jobs.append(var)

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='1:00:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs.append(lin)

  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='1:00:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)

  return jobs

def si_pyscf_test():
  ''' Simple tests that check PBC is working Crystal, and that QMC can be performed on the result.'''
  jobs=[]

  pwriter=PySCFPBCWriter({
      'gmesh':[16,16,16],
      'cif':open('si.cif','r').read()
    })
  pman=PySCFManager(
      name='scf',
      path='test_si_pyscf',
      writer=pwriter,
      runner=PySCFRunnerPBS(
          queue='secondary',
          nn=1,
          np=16,
          walltime='4:00:00'
        )
    )
  jobs.append(pman)

  var=QWalkManager(
      name='var',
      path=pman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(pman,kpoint=0)
    )
  jobs.append(var)

  lin=QWalkManager(
      name='lin',
      path=pman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='0:30:00'
        ),
      trialfunc=SlaterJastrow(slatman=pman,jastman=var,kpoint=0)
    )
  jobs.append(lin)

  return jobs

def mno_test():
  ''' Test spinful calculations in PBC and real kpoint sampling. Also test trialfunc dependency. '''
  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'kmesh':(3,3,3),
      'initial_spins':[1,0],
      'total_spin':5
    })
  cwriter.set_struct_fromcif(open('structures/mno.cif','r').read(),primitive=True)
  cwriter.set_options({})

  cman=CrystalManager(
      name='crys',
      path='test_mno_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=20,
          walltime='2:00:00'
        )
    )

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=20,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=20,queue='secondary',walltime='1:00:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs=[cman,var,lin]
  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='1:00:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)
  return jobs

###################################################################################################################
# Test suite definitions.
def quick_crystal():
  ''' Fast-running tests of main crystal features. '''
  jobs=[]

  # Spinless molecule
  jobs+=h2_crystal_equil_test()

  # Spinless crystal
  jobs+=si_crystal_test()

  # Spinful molecule
  jobs+=h2_crystal_stretch_test()

  # Spinful crystal
  jobs+=mno_test()

  return jobs

def quick_pyscf():
  ''' Fast-running tests of main pyscf features. '''
  # Spinless molecule

  # Spinless crystal

  # Spinful molecule

  # Spinful crystal
  return [] #TODO

def thorough_crystal():
  ''' Longer tests that look over additional crystal features, compared to the quick tests.
  You should probably run the quick tests before running these.'''
  return [] # TODO

def thorough_pyscf():
  ''' Longer tests that look over additional pyscf features, compared to the quick tests.
  You should probably run the quick tests before running these.'''
  return [] # TODO

###################################################################################################################
# Test operations
def run_tests():
  ''' Generate data using current generation of autogen.'''
  jobs=[]

  # Tests that are run.
  jobs+=quick_crystal()

  for job in jobs:
    job.nextstep()

def check_tests(reffn='refs.json'):
  ''' Check the results of run_tests with reference data. '''
  from json import load
  jobs=[]
  report=[]

  # Tests that are compared to reference data.
  jobs+=quick_crystal()

  with open('refs.json','r') as inpf:
    refs=load(inpf)

  for job in jobs:
    jobtype=job.__class__.__name__
    if jobtype=='CrystalManager':
      issame=compare_crystal(job,refs[job.path+job.name])
    elif jobtype=='QWalkManager':
      issame=compare_qwalk(job,refs[job.path+job.name])
    else:
      raise NotImplementedError("No routine for checking %s results yet"%jobtype)

    if not issame:
      report.append("%s: jobs differ more than tolerance, or else aren't run yet."%(job.path+job.name))

  print("#######################################")
  print("### Results of tests ##################" )
  print("%d jobs failed to match reference."%len(report))
  print('\n'.join(report))

def compare_crystal(job,ref):
  ''' Make sure two crystal jobs have the same results within machine precision.
  Args:
    job (Manager): What you'd like to check.
    ref (dict): Reference data.
  '''
  issame=False

  job.collect()

  issame=abs(job.creader.output['total_energy']-ref['total_energy'])<1e-15

  return issame

def compare_qwalk(job,ref,nsigma=3):
  ''' Make sure two qwalk jobs have the same results within error.
  Args:
    job (Manager): What you'd like to check. 
    ref (dict): Reference data.
    nsigma (float): Number of standard devations to allow before test declairs the results are different.
  '''
  issame=False

  job.collect()

  try:
    issame=abs(job.reader.output['total_energy'] - ref['total_energy']) < \
        nsigma*(job.reader.output['total_energy_err']**2 + ref['total_energy_err']**2)**0.5
  except KeyError:
    pass # probably not the right QMC type.

  try:
    issame=abs(job.reader.output['properties']['total_energy']['value'][0] - ref['total_energy']) < \
        nsigma*(job.reader.output['properties']['total_energy']['error'][0]**2 + ref['total_energy_err']**2)**0.5
  except KeyError:
    pass # probably not the right QMC type.

  try:
    issame=abs(job.reader.output['sigma'] - ref['sigma']) < \
        ref['sigma'] # Not quite right, but will catch really bad runs.
  except KeyError:
    pass # probably not the right QMC type.

  return issame

def update_refs(reffn='refs.json'):
  ''' Update the references file.
  Args:
    reffn: Reference file name that will be produced.
  '''
  from json import dump
  jobs=[]
  refs={}
  jobs+=quick_crystal()

  for job in jobs:
    refs[job.path+job.name]=grab_job_ref(job)

  with open(reffn,'w') as outf:
    dump(refs,outf)

def grab_job_ref(job):
  ref={}
  job.collect()

  try:
    reader=job.reader
  except AttributeError:
    reader=job.creader

  try:
    ref['total_energy'] = reader.output['total_energy']
  except KeyError:
    try:
      ref['total_energy'] =reader.output['properties']['total_energy']['value'][0]
    except KeyError:
      pass

  try:
    ref['total_energy_err']=reader.output['total_energy_err']
  except KeyError:
    try:
      ref['total_energy_err']=reader.output['properties']['total_energy']['error'][0]
    except KeyError:
      pass

  try:
    ref['sigma']=reader.output['sigma']
  except KeyError:
    try:
      ref['sigma']=reader.output['properties']['total_energy']['sigma'][0]
    except KeyError:
      pass

  return ref

if __name__=='__main__':
  run_tests()
  check_tests()
  #update_refs()
