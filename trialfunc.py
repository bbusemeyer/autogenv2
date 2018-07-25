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
  def __init__(self,orbitals,uporbs,downorbs,detweights=(1.0,)):
    ''' Slater wave function object.

    Args: 
      orbitals (Orbitals): Orbitals from which you'll select the orbitals. 0-based indexing.
      uporbs (array-like): which orbitals used for up electrons.
      downorbs (array-like): which orbitals used for down electrons.
      kpoint (tuple): kpoint in same form as in Orbitals.
    '''
    uporbs=array(uporbs)
    downorbs=array(downorbs)

    # A lot can go wrong with this input. Raise hell if something is wrong.
    if len(detweights)>1:
      assert len(uporbs.shape)==2 or len(downorbs.shape)==2,\
          "Multiple deterimants require specifying weights."
      assert len(detweights)==uporbs.shape[0] or len(detweights)==downorbs.shape[0],\
          "Number of determinant weights should be same as number of deteriminants."
    else:
      uporbs=array([uporbs])
      downorbs=array([downorbs])
      assert (len(uporbs.shape)==1 or uporbs.shape[0]==1),\
          "Multiple deterimants require specifying weights."

    self.orbitals=orbitals
    self.detweights=array(detweights)
    self.uporbs=array(uporbs)
    self.downorbs=array(downorbs)
   
  #------------------------------------------------
  def export(self):
    ''' Export a Slater section for use in a trialfunction.
    Returns:
      str: Slater wave function section for QWalk.
    '''

    uporbs = array(self.uporbs)+1
    downorbs = array(self.downorbs)+1

    if any([(e.imag!=0.0).any() for e in self.orbitals.eigvecs]):
      orbstr = "corbitals"
    else:
      orbstr = "orbitals"

    print(uporbs.astype(str))
    uporblines = [' '.join(det) for det in uporbs.astype(str)]
    downorblines = [' '.join(det) for det in downorbs.astype(str)]

    outlines = [
        "slater",
        "{0} {{".format(orbstr),
        "  magnify 1",
        "  nmo {0}".format(max(uporbs.max(),downorbs.max())),
        "  orbfile {0}".format(self.orbitals.last_orbfile),
        self.orbitals.export_qwalk_basis(),
        "  centers { useglobal }",
        "}",
        "detwt {{ {} }}".format(' '.join(self.detweights.astype(str))),
        "states {"
      ]
    for detweight,upline,downline in zip(self.detweights,uporblines,downorblines):
      outlines+=[
          "  # Spin up orbitals detweight {}.".format(detweight), 
          "  " + upline,
          "  # Spin down orbitals detweight {}.".format(detweight), 
          "  " + downline,
          "}"
        ]
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
