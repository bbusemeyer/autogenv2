import numpy as np
periodic_table = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
    'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po',
    'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
    'Mt', 'Ds', 'Rg', 'Cp', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo']
bohr=1.88972598858 # in Angstroms

###########################################################################################
class System:
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
    self.latparm={}

  # ----------------------------------------------------------------------------------------
  def import_xyz(self):
    ''' Generate a molecule's Structure from xyz file.'''
    from pymatgen.io.xyz import XYZ
    struct=XYZ.from_string(xyzstr).molecule.as_dict()

  # ----------------------------------------------------------------------------------------
  def import_cif_pymatgen(self,cifstr,primitive=True):
    ''' Import the positions and lattice parameters from CIF file using pymatgen.
    This importer generates all symmetry inequivilent positions here, but cannot handle symmetry.
    Args:
      cifstr (str): *contents* (not filename) of a CIF file.
      primitive (bool): pymatgen will attempt to find smallest possible cell.
    '''
    from pymatgen.io.cif import CifParser
    from pymatgen.core.periodic_table import Element
    pydict=CifParser.from_string(cifstr).get_structures(primitive=primitive)[0].as_dict()

    for site in pydict['sites']:
      assert len(site['species'])==1,\
          'Multiple site not tested. Check this works ok and remove this assertion.'
      element=site['species'][0]['element']
      self.positions.append({
            'species':element,
            'abc':[float(c) for c in site['abc']], # TODO store only one.
            'xyz':[float(c) for c in site['xyz']]
          })
    for key in 'alpha','beta','gamma','a','b','c':
      self.latparm[key]=pydict['lattice'][key]
    self.latparm['latvecs']=np.array(pydict['lattice']['matrix'])*bohr

  # This is the safer import method for now.
  import_cif=import_cif_pymatgen

  # ----------------------------------------------------------------------------------------
  def import_cif_raw(self,ciffn):
    ''' Import the positions and lattice parameters from CIF file.
    This importer relies on the code to produce symmetry inequivilent points, but can handle symmetry.
    This also doesn't provide lattice vectors and the like, so you'll need a different routine to get those.
    pymatgen can generate them based off the lengths and angles.

    Nonetheless you can use this to generate a crystal output with symmetry.
    '''
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
      pseudo['core_charge']=eff_core_charge
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
  def export_crystal_geom(self,supercell=None):
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
      elemz=periodic_table.index(site['species'])+1
      # TODO assumes psuedopotential.
      geomlines+=[str(elemz+200)+" %g %g %g"%tuple(site['abc'])]

    if supercell is not None:
      geomlines+=["SUPERCELL"]
      for row in supercell:
        geomlines+=[' '.join(map(str,row))]

    return geomlines

  # ----------------------------------------------------------------------------------------
  def export_qwalk_sys(self,cutoff_divider,nspin,kpoint=(0.0,0.0,0.0)):
    ''' Generate a system section for QWalk.
    Args:
      cutoff divider (float): calculate using Basis object. 7.5 is a good default for finite systems.
      kpoint (tuple): boundary condition.
    Returns: 
      list: List of the lines (str) making up the qwalk system section.
    '''
    outlines = []

    # Assumes 0-d or 3-d here.
    if self.latparm != {}:
      outlines += [
          "system { periodic",
          "  nspin {{ {} {} }}".format(*nspin),
          "  latticevec {",
        ]
      for i in range(3):
        outlines.append("    {:< 15} {:< 15} {:< 15}".format(*[lc for lc in self.latparm['latvecs'][i]]))
      outlines += [
          "  }",
          "  origin { 0 0 0 }",
          "  cutoff_divider {0}".format(cutoff_divider),
          "  kpoint {{ {:4}   {:4}   {:4} }}".format(*kpoint)
        ]
    else: # is molecule.
      outlines += [
          "system { molecule",
          "  nspin {{ {} {} }}".format(*nspin)
        ]
    for position in self.positions:
      if position['species'] in self.pseudo:
        charge=self.pseudo[position['species']]['core_charge']
      else:
        assert 0, "Check part and delete this assertion!"
        charge=periodic_table.find(position['species'])+1
      outlines.append(
        "  atom {{ {0} {1} coor {2} }}".format(
          position['species'],
          charge,
          "{:< 15} {:< 15} {:< 15}".format(*[p for p in position['xyz']])
        )
      )
    outlines.append("}")
    done = []
    for species in self.pseudo:
      nonlocal_counts=[0,0,0,0,0,0] # Extra long to be safe.
      for term in self.pseudo[species]['nonlocal']:
        nonlocal_counts[term['angular']]+=1
      for idx in range(len(nonlocal_counts)-1):
        assert nonlocal_counts[idx]!=0 or nonlocal_counts[idx+1]==0,\
          'Not sure if QWalk can handle a pseudopotential like this,'
      nonlocal_counts=[c for c in nonlocal_counts if c>0]
      num_types=1 + len(nonlocal_counts)

      if num_types > 2: aip = 12
      else:             aip =  6

      outlines += [
          "pseudo {",
          "  {}".format(species),
          "  aip {:d}".format(aip),
          "  basis {{ {}".format(species),
          "    rgaussian",
          "    oldqmc {",
          "      0.0 {:d}".format(num_types),
          "      "+' '.join(["{}" for i in range(num_types)])\
              .format(*(nonlocal_counts+[len(self.pseudo[species]['local'])]))
        ]
      last=-1
      for gaussian in self.pseudo[species]['nonlocal']:
        assert gaussian['angular']>last,\
            "Wait, I thought the order is always ascending! Well this part needs some more code."
        last=gaussian['angular']
        outlines.append("      {:d}   {:<12} {:< 12}".format(
          gaussian['r_to_n']+2,
          gaussian['exp'],
          gaussian['coef']
        ))
      for gaussian in self.pseudo[species]['local']:
        outlines.append("      {:d}   {:<12} {:< 12}".format(
          gaussian['r_to_n']+2,
          gaussian['exp'],
          gaussian['coef']
        ))
      outlines += ["    }","  }","}"]
    return '\n'.join(outlines)+'\n'

  # ----------------------------------------------------------------------------------------
  def export_jast(self,threebody=False):
    ''' Makes a unoptimized Jastrow wave function section which should be optimizized before using in QMC.
    Returns:
      str: Jastrow section.
    '''
    if threebody: raise NotImplementedError("Should be simple to add three-body, but haven't bothered yet.")
      
    basis_cutoff = find_basis_cutoff(self.latparm['latvecs'])
    atom_types = set([pos['species'] for pos in self.positions])
    outlines = [
        "jastrow2",
        "group {",
        "  optimizebasis",
        "  eebasis {",
        "    ee",
        "    cutoff_cusp",
        "    gamma 24.0",
        "    cusp 1.0",
        "    cutoff {0}".format(basis_cutoff),
        "  }",
        "  eebasis {",
        "    ee",
        "    cutoff_cusp",
        "    gamma 24.0",
        "    cusp 1.0",
        "    cutoff {0}".format(basis_cutoff),
        "  }",
        "  twobody_spin {",
        "    freeze",
        "    like_coefficients { 0.25 0.0 }",
        "    unlike_coefficients { 0.0 0.5 }",
        "  }",
        "}",
        "group {",
        "  optimize_basis",
      ]
    for atom_type in atom_types:
      outlines += [
        "  eibasis {",
        "    {0}".format(atom_type),
        "    polypade",
        "    beta0 0.2",
        "    nfunc 3",
        "    rcut {0}".format(basis_cutoff),
        "  }"
      ]
    outlines += [
        "  onebody {",
      ]
    for atom_type in atom_types:
      outlines += [
        "    coefficients {{ {0} 0.0 0.0 0.0}}".format(atom_type),
      ]
    outlines += [
        "  }",
        "  eebasis {",
        "    ee",
        "    polypade",
        "    beta0 0.5",
        "    nfunc 3",
        "    rcut {0}".format(basis_cutoff),
        "  }",
        "  twobody {",
        "    coefficients { 0.0 0.0 0.0 }",
        "  }",
        "}"
      ]

    return "\n".join(outlines)

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

######################################################################################
def find_cutoff_divider(latvecs,min_exp):
  basis_cutoff=find_basis_cutoff(latvecs)
  cutoff_length=(-np.log(1e-8)/min_exp)**.5

  return basis_cutoff*2.0 / cutoff_length

######################################################################################
def find_basis_cutoff(latvecs):
  cutoff_divider = 2.000001
  cross01 = np.cross(latvecs[0], latvecs[1])
  cross12 = np.cross(latvecs[1], latvecs[2])
  cross02 = np.cross(latvecs[0], latvecs[2])

  heights = [0,0,0]
  heights[0]=abs(np.dot(latvecs[0], cross12)/np.dot(cross12,cross12)**.5)
  heights[1]=abs(np.dot(latvecs[1], cross02)/np.dot(cross02,cross02)**.5)
  heights[2]=abs(np.dot(latvecs[2], cross01)/np.dot(cross01,cross01)**.5)
  return min(heights)/cutoff_divider
