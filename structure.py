from pymatgen.io.xyz import XYZ
periodic_table = [
  "h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p","s","cl","ar",
  "k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br",
  "kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te",
  "i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm",
  "yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn",
  "fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr",
  "rf","db","sg","bh","hs","mt","ds","rg","cp","uut","uuq","uup","uuh","uus","uuo"
]

###########################################################################################
def scrub_err(numstr):
  ''' Remove error bar notation.

  Example 
  > scrub_err('123.190(2)')
  > 123.190
  '''
  if '(' in numstr:
    return float(numstr[:numstr.find('(')])
  else:
    return float(numstr)

###########################################################################################
def space_group_format(group_number):  
  ''' Generate the format of the symmetry group input format for crystal.  
  Returns: 
    str: .format()-able string.  
  '''  
  format_string="" 

  if 1<=group_number<3: 
    #print("Case triclinic") 
    format_string="{a} {b} {c} {alpha} {beta} {gamma}" 

  elif 3<=group_number<16:  
    #print("Case monoclinic")  
    format_string="{a} {b} {c} {beta}" 

  elif 16<=group_number<75: 
    #print("Case orthorhombic")  
    format_string="{a} {b} {c}"  
     
  elif 75<=group_number<143:  
    #print("Case tetragonal")  
    format_string="{a} {c}"  

  elif 143<=group_number<168: 
    #print("Case trigonal")  
    format_string="{a} {c}"  

  elif 168<=group_number<195: 
    #print("Case trigonal")  
    format_string="{a} {c}"  

  elif 167<=group_number<231: 
    #print("Case cubic") 
    format_string="{a}"  

  else:  
    raise AssertionError("Invalid group_number") 

  return format_string

###########################################################################################
class Structure:
  ''' Class containing information:
    - atom placement,
    - atom types,
    - pseudopotentials (if any),
    - spin,
    - k-points (if any),
    - lattice vectors (if periodic).

  Can also export to various file outputs.'''

  # ----------------------------------------------------------------------------------------
  def __init__(self):
    self.positions=[]
    self.pseudo={}
    #self.nspin=(1,1) # Not easy to extract. Is it really a property of structure?
    self.group_number=1
    self.supercell=None
    self.latparm={}

  # ----------------------------------------------------------------------------------------
  def import_xyz(self):
    ''' Generate a molecule's Structure from xyz file.'''
    struct=XYZ.from_string(xyzstr).molecule.as_dict()

  # ----------------------------------------------------------------------------------------
  def import_cif_pymatgen(self,cifstr,primitive=True):
    ''' Import the positions and lattice parameters from CIF file using pymatgen.
    This importer generates all symmetry inequivilent positions here, but cannot handle symmetry.'''
    from pymatgen.io.cif import CifParser
    pystruct=CifParser.from_string(cifstr).get_structures(primitive=primitive)[0].as_dict()

    for site in pystruct['sites']:
      assert len(site['species'])==1,\
          'Multiple site not tested. Check this works ok and remove this assertion.'
      element=site['species'][0]['element']
      self.positions.append({
            'species':element,
            'abc':[float(c) for c in site['abc']],
            'xyz':[float(c) for c in site['xyz']]
          })
    for key in 'alpha','beta','gamma','a','b','c':
      self.latparm[key]=pystruct['lattice'][key]

  # This is the safer import method for now.
  import_cif=import_cif_pymatgen

  # ----------------------------------------------------------------------------------------
  def import_cif_raw(self,ciffn):
    ''' Import the positions and lattice parameters from CIF file.
    This importer relies on the code to produce symmetry inequivilent points, but can handle symmetry.'''
    from CifFile import CifFile,ReadCif

    struct=ReadCif(ciffn)
    assert len(struct)==1,\
        'Can only handle CIFs containing one structure for now.'
    struct=struct[struct.keys()[0]]
    for label in struct['_atom_site_label']:
      self.positions.append({
        'species':''.join([c for c in label if c.isalpha()])
        })
    for idx,abc in enumerate(zip(struct['_atom_site_fract_x'],struct['_atom_site_fract_y'],struct['_atom_site_fract_z'])):
      self.positions[idx]['abc']=[scrub_err(num) for num in abc]

    for key in 'alpha','beta','gamma':
      self.latparm[key]=scrub_err(struct['_cell_angle_%s'%key])
    for key in 'a','b','c':
      self.latparm[key]=scrub_err(struct['_cell_length_%s'%key])

    self.group_number=int(struct['_symmetry_Int_Tables_number'])

  # ----------------------------------------------------------------------------------------
  def lookup_pseudopotential(self,xml_name='BFD_Library.xml',species_list=None):
    from xml.etree.ElementTree import ElementTree
    ''' Lookup pseudopotentials for the atoms in the structure. 

    Args:
      xml (str): path to xml for lookup.
      species (list): List of species to have pseudopotentials for. Default: all atoms currently in structure.
    '''
    if species_list is None:
      species_list=set([s['species'] for s in self.positions])
    for species in species_list:
      pseudo={}
      tree = ElementTree()
      tree.parse(xml_name)
      element = tree.find('./Pseudopotential[@symbol="{}"]'.format(species))
      eff_core_charge = element.find('./Effective_core_charge').text
      local_path = './Gaussian_expansion/Local_component'
      non_local_path = './Gaussian_expansion/Non-local_component'
      local_list = element.findall(local_path)
      non_local_list = element.findall(non_local_path)
      proj_list = element.findall('./Gaussian_expansion/Non-local_component/Proj')
      pseudo['local']=[{
        'exp':float(lc.find('./Exp').text),
        'coef':float(lc.find('./Coeff').text),
        'r_to_n':int(lc.find('./r_to_n').text)
        } for lc in local_list]
      pseudo['nonlocal']=[{
        'angular':int(pl.text),
        'exp':float(nlc.find('./Exp').text),
        'coef':float(nlc.find('./Coeff').text),
        'r_to_n':int(nlc.find('./r_to_n').text)
        } for pl,nlc in zip(proj_list,non_local_list)]
      self.pseudo[species]=pseudo

  # ----------------------------------------------------------------------------------------
  def export_crystal_geom(self):
    ''' 
    Returns:
      list: List of the lines (str) making up the geometry section.
    '''

    geomlines=[
        "CRYSTAL","0 0 0",
        str(self.group_number),
        space_group_format(self.group_number).format(**self.latparm)
      ]

    geomlines+=["%i"%len(self.positions)]
    for site in self.positions:
      elemz=periodic_table.index(site['species'].lower())+1
      # TODO assumes psuedopotential.
      geomlines+=[str(elemz+200)+" %g %g %g"%tuple(site['abc'])]

    if self.supercell is not None:
      geomlines+=["SUPERCELL"]
      for row in self.supercell:
        geomlines+=[' '.join(map(str,row))]

    return geomlines

  # ----------------------------------------------------------------------------------------
  def export_qwalk_sys(self,nspin):
    ''' Generate a system section for QWalk.
    Args:
      nspin (tuple): Number of up and down electrons in the system.
    Returns: 
      list: List of the lines (str) making up the qwalk system section.
    '''
