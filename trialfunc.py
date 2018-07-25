# You can easily define new wave functions here. 
# The only requirement is to define the export() method, which defines how to generate the QWalk input. 
import os
from numpy import array
from manager_tools import separate_jastrow

#######################################################################
def export_trialfunc(wf):
  ''' Convenience function for wrapping the wave function output with trialfunc. 
  Args:
    wf (wave function object): wave function you want as your trial function. Should have an export routine.
  Returns:
    str: trialfunc section.
  '''
  outlines=['trialfunc { ']+\
      ['  '+line for line in wf.export().split('\n')]+\
      ['}']
  return '\n'.join(outlines)+'\n'

#######################################################################
class Slater:
  def __init__(self,orbitals,orbfile,states,detweights=(1.0,)):
    ''' Slater wave function object.

    Args: 
      orbitals (Orbitals): Orbitals from which you'll select the orbitals. 0-based indexing.
      orbfile (str): file name where orbitals has written to. See orbitals.write_qwalk_orb().
        Same string as QWalk orbfile.
      states (array-like): states[spin channel][determinant][orbital] select orbitals for determinants.
      detweights (array-like): Weights of determinants for multideterminant expansion. First should be 1.0.
    '''
    self.orbitals=orbitals
    self.orbfile=orbfile
    self.detweights=array(detweights)
    self.states=array(states)
    assert (len(self.states.shape)==3) and (self.states.shape[0]==2) and (self.states.shape[1]==self.detweights.shape[0]),\
        "States array should be nspin by ndeterminant by nelectrons. One detweight per determinant."
   
  #------------------------------------------------
  def export(self):
    ''' Export a Slater section for use in a trialfunction.
    Returns:
      str: Slater wave function section for QWalk.
    '''

    states = array(self.states)+1

    if any([(e.imag!=0.0).any() for e in self.orbitals.eigvecs]):
      orbstr = "corbitals"
    else:
      orbstr = "orbitals"

    upstatelines = [' '.join(det) for det in states[0].astype(str)]
    downstatelines = [' '.join(det) for det in states[1].astype(str)]

    outlines = [
        "slater",
        "{0} {{".format(orbstr),
        "  magnify 1",
        "  nmo {0}".format(states.max()),
        "  orbfile {0}".format(self.orbfile),
        self.orbitals.export_qwalk_basis(),
        "  centers { useglobal }",
        "}",
        "detwt {{ {} }}".format(' '.join(self.detweights.astype(str))),
        "states {"
      ]
    for detweight,upline,downline in zip(self.detweights,upstatelines,downstatelines):
      outlines+=[
          "  # Spin up orbitals detweight {}.".format(detweight), 
          "  " + upline,
          "  # Spin down orbitals detweight {}.".format(detweight), 
          "  " + downline
        ]
    outlines+=['}']
    return "\n".join(outlines)

  #------------------------------------------------
  export_trialfunc=export_trialfunc

#######################################################################
# Quick and dirty version for now.
class Jastrow:
  def __init__(self,jastfn=None):
    ''' Jastrow wave function object.
    Args:
      jastfn (str): file name to read with Jastrow section..
    '''
    self.text=''
    if jastfn is not None:
      self.text=separate_jastrow(jastfn)

  #------------------------------------------------
  def export(self):
    ''' Export a Slater section for use in a trialfunction.
    Returns:
      str: Jastrow wave function section for QWalk.
    '''

    return self.text
  #------------------------------------------------
  export_trialfunc=export_trialfunc

#######################################################################
class SlaterJastrow:
  ''' Wave function that multiplies a slater/multislater wave function by a jastrow.'''
  def __init__(self,slater,jastrow):
    self.slater=slater
    self.jastrow=jastrow

  #------------------------------------------------
  def export(self):
    ''' Export a Slater section for use in a trialfunction.
    Args: 
      orbfn: File name to write the orbitals to. Defaults to "slater_%d_%d_%d"%self.orbitals.kpoint.
    Returns:
      str: Slater wave function section for QWalk.
    '''

    outlines=[
        'slater-jastrow',
        '  wf1 {'
        ]+['    '+line for line in self.slater.export().split('\n')]+[
        '  }',
        '  wf2 {'
        ]+['    '+line for line in self.jastrow.export().split('\n')]+[
        '  }'
      ]
    
    return '\n'.join(outlines)
  #------------------------------------------------
  export_trialfunc=export_trialfunc
