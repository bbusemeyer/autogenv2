class Orbitals:
  ''' Class containing information defining a set of single-particle orbitals.
  - Basis set for the orbitals.
  - Coefficients of the orbitals in the basis.
  '''

  def __init__():
    self.basis={}
    self.coefs=[]

  def export_qwalk_orbs(occpy=(,),deoccpy=(,),nspin=(1,1)):
    ''' Export an orb file for QWalk.

    Args:
      occpy (iterable): Modify the ground state by occupying these orbitals.
      deoccupy (iterable): Modify the ground state by deoccupying these orbitals.
    '''
    raise NotImplementedError

  def export_qwalk_basis():
    ''' Export the basis set for these orbitals.'''
    raise NotImplementedError

  def export_qwalk_orbs(occpy=(,),deoccpy=(,),nspin=(1,1)):
    ''' Export an SCF object with these orbitals.'''
    raise NotImplementedError
