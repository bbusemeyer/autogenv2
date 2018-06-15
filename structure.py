class Structure:
  ''' Class containing information:
    - atom placement,
    - atom types,
    - pseudopotentials (if any),
    - spin,
    - k-points (if any),
    - lattice vectors (if periodic).

  Can also export to various file outputs.'''

  def __init__():
    self.positions={}
    self.psuedo={}
    self.nspin=(1,1)
    self.kpoint=(0,0,0)
    self.lattice=((1,0,0),(0,1,0),(0,0,1))
    self.groupnumber=1
    self.latparm={}

  def import_xyz():
    ''' Generate a molecule's Structure from xyz file.'''
    raise NotImplementedError

  def import_cif():
    ''' Generate a crystal's Structure from xyz file.'''
    raise NotImplementedError

  def export_crystal_geom():
    ''' Generate the geometry section of a Crystal input.'''
    raise NotImplementedError

  def export_qwalk_sys():
    ''' Generate a system section for QWalk.'''
    raise NotImplementedError

