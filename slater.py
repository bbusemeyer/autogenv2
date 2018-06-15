class Slater:
  ''' Class defining a Slater or multi-Slater determinant wave function. 
  - Determinant weights, occupations
  - Orbitals.
  - CSF 
  '''
  def __init__():
    self.occupation=[{'occ':(,),'weight':1.0}]
    self.orbitals=None

  def export_qwalk_wf():
    ''' Export a QWalk wave function section. '''
    raise NotImplementedError

  def export_pyscf_wf():
    ''' Export a PySCF wave function section. '''
    raise NotImplementedError
