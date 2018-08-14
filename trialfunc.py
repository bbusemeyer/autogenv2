''' Trial function objects facilitate combining orbitals and other options to generating trial function input.'''
from numpy import array

#######################################################################
def export_qwalk_trialfunc(wf,**kwargs):
  ''' Wrapper for making making a wave function section into a trialfunc.
  Args:
    wf (TrialFunc): wave function for QWalk.
  Returns:
    str: trialfunc section.
  '''
  outlines=['trialfunc { ']+\
      ['  '+line for line in wf.export_qwalk_wf(**kwargs).split('\n')]+\
      ['}']
  return '\n'.join(outlines)

#######################################################################
class TrialFunc:
  ''' Abstract wave function class defined minimal interface.'''
  def export_qwalk_wf(self,**kwargs):
    ''' This must be implmented. '''
    raise NotImplementedError('Trial function not completely implemented')

#################################################################################################
# Simple implementation for now. 
# In the future this should actually store the data for the wave function. 
class Jastrow(TrialFunc):
  ''' Class representing a Jastrow wave function. '''
  def __init__(self,jaststr):
    self.jaststr=jaststr

  def export_qwalk_wf(self,**kwargs):
    return self.jaststr

#################################################################################################
class Slater:
  ''' Class representing a slater determinant wave function. '''
  def __init__(self,weights,states,orbfile,orbitals,shift_downorb=False):
    '''
    Args: 
      weights (array-like): Weights of determinants for multideterminant expansion. First should be 1.0.
      states (array-like): states[determinant][spin channel][orbital] select orbitals for determinants.
        Indicies should reference whats written in the orbfile.
      orbfile (str): where orbitals are stored on disk (see write_qwalk_orb).
      orbitals (Orbitals): Something that can export_qwalk_orbitals(orbfile).
      shift_downorb (bool): Shift states[1] by number of up orbitals. 
        Useful for unrestricted calculations.
    '''
    self.weights=weights
    self.states=states
    self.orbfn=orbfile
    self.orbitals=orbitals
    self.shift_downorb=shift_downorb*(orbitals.eigvecs[0].shape[0])

  def export_qwalk_wf(self,optimize_det=False):
    ''' Write out a qwalk wave function section .
    Returns:
      str: wave function section.
    '''

    # Convert python indexing to QWalk indexing.
    states=array(self.states)+1
    states[:,1,:]+=self.shift_downorb
    weights=array(self.weights)

    # Check input validity.
    assert (len(states.shape)==3) and (states.shape[1]==2) and (states.shape[0]==weights.shape[0]),\
        "States array should be nweights({}) by nspin(2) by nelectrons. One detweight per determinant. Its shape is {}".format(weights.shape[0],states.shape)

    upstatelines = [' '.join(det[0]) for det in states.astype(str)]
    downstatelines = [' '.join(det[1]) for det in states.astype(str)]

    if optimize_det: optimize_det_lines=['optimize_det']
    else:            optimize_det_lines=[]

    outlines = [
        "slater",
      ]+optimize_det_lines+[
        self.orbitals.export_qwalk_orbitals(self.orbfn),
        "detwt {{ {} }}".format(' '.join(weights.astype(str))),
        "states {"
      ]
    for detweight,upline,downline in zip(weights,upstatelines,downstatelines):
      outlines+=[
          "  # Spin up orbitals detweight {}.".format(detweight), 
          "  " + upline,
          "  # Spin down orbitals detweight {}.".format(detweight), 
          "  " + downline
        ]
    outlines+=['}']
    return "\n".join(outlines)

#################################################################################################
class SlaterJastrow:
  ''' Class rpresenting a slater determinant wave function. '''
  def __init__(self,slater,jastrow):
    '''
    Args: 
      slater (Slater): Slater wave function to multiply.
      jastrow (Jastrow): Jastrow wave function to multiply.
    '''
    self.slater=slater
    self.jastrow=jastrow
  def export_qwalk_wf(self,**kwargs):
    ''' Write out a qwalk wave function section .
    Returns:
      str: wave function section.
    '''
    outlines=[
        'slater-jastrow',
        '  wf1 {'
        ]+['    '+line for line in self.slater.export_qwalk_wf(**kwargs).split('\n')]+[
        '  }',
        '  wf2 {'
        ]+['    '+line for line in self.jastrow.export_qwalk_wf(**kwargs).split('\n')]+[
        '  }'
      ]
    return '\n'.join(outlines)
