import crystal2qmc
from numpy import array
from autogen_tools import normalize_eigvec

#######################################################################
def export_trialfunc(wfsection):
  ''' Wrapper for making making a wave function section into a trialfunc.
  Args:
    wfsection (str): wave function section for QWalk.
  Returns:
    str: trialfunc section.
  '''
  outlines=['trialfunc { ']+\
      ['  '+line for line in wfsection.split('\n')]+\
      ['}']
  return '\n'.join(outlines)+'\n'

#######################################################################
def make_slaterjastrow(slatersection,jastrowsection):
  ''' Wrapper for joining Slater and Jastrow in a sacred union. 

  Args:
    slatersection (str): slater wave function section for QWalk.
    jastrowsection (str): slater wave function section for QWalk.
  Returns:
    str: Slater-Jastrow wave function section.
  '''
  outlines=[
      'slater-jastrow',
      '  wf1 {'
      ]+['    '+line for line in slatersection.split('\n')]+[
      '  }',
      '  wf2 {'
      ]+['    '+line for line in jastrowsection.split('\n')]+[
      '  }'
    ]
  return '\n'.join(outlines)

#################################################################################################
class Orbitals:
  ''' Class containing information defining a set of single-particle orbitals.
  - Basis set for the orbitals.
  - Coefficients of the orbitals in the basis.
  '''

  #----------------------------------------------------------------------------------------------
  def __init__(self):
    self.basis={}
    self.eigvecs=[]
    self.eigvals=[]
    self.atom_order=()
    self.kpt_weight=0.0
    self.kpoint=(0.0,0.0,0.0)
    self.last_orbfile=None # Last path orbfile was written to.

  #----------------------------------------------------------------------------------------------
  def write_qwalk_orb(self,outfn,nperspin):
    ''' Generate a orb file for QWalk. 
    This just writes to the file because orbfile are necessarily separate in QWalk.

    Args:
      nperspin (tuple): (Nup,Ndown) number of electrons in each spin channel.
      outfn (str): file to write to.
    '''
    nspin=len(self.eigvecs)
    outf=open(outfn,'w')

    nao_atom = count_naos(self.basis)
    totnmo = (max(nperspin))*nspin

    # Do the printing.
    coef_cnt = 0
    for moidx in range(totnmo):
      for atidx,atom in enumerate(self.atom_order):
        for aoidx in range(nao_atom[atom]):
          outf.write(" {:5d} {:5d} {:5d} {:5d}\n"\
              .format(moidx+1,aoidx+1,atidx+1,coef_cnt+1))
          coef_cnt += 1
    eigvec_flat = [normalize_eigvec(self.eigvecs[s].copy(),basis).ravel() for s in range(nspin)]
    print_cnt = 0
    outf.write("COEFFICIENTS\n")
    if any([(e.imag!=0.0).any() for e in self.eigvecs]):
      for eigv in eigvec_flat: 
        for r,i in zip(eigv.real,eigv.imag):
          outf.write("({:<.12e},{:<.12e}) "\
              .format(r,i))
          print_cnt+=1
          if print_cnt%5==0: outf.write("\n")
    else:                        # Real coefficients
      for eigr in eigvec_flat:
        for r in eigr:
          outf.write("{:< 15.12e} ".format(r))
          print_cnt+=1
          if print_cnt%5==0: outf.write("\n")

  #----------------------------------------------------------------------------------------------
  def export_pyscf_basis(self):
    from pyscf.gto.basis import parse
    angmap={
          'S':'s',
          'P':'p',
          '5D':'d',
          '7F_crystal':'f',
          'G':'g',
          'H':'h'
        }

    pyscfbasis={}
    for species in self.basis:
      outlines=[]
      for element in self.basis[species]:
        outlines+=['%s %s'%(species,angmap[element['angular']])]
        for exp,coef in zip(element['exponents'],element['coefs']):
          outlines+=['  %.16f %.16f'%(exp,coef)]
      pyscfbasis[species]=parse('\n'.join(outlines))

    return pyscfbasis

  #----------------------------------------------------------------------------------------------
  def export_qwalk_basis(self):
    ''' Generate a basis section for QWalk.
    Returns:
      list: lines pertaining to basis section.
    '''
    outlines=[]

    for species in self.basis:
      outlines+=[
          'basis { ',
          '  %s'%species.capitalize(),
          '  aospline',
          '  normtype CRYSTAL',
          '  gamess {'
        ]
      for element in self.basis[species]:
        numprim=element['coefs'].shape[0]
        outlines+=['    %s %d'%(element['angular'],numprim)]
        outlines+=['    %d %.16f %.16f'%(idx+1,exp,coef)
            for idx,exp,coef in zip(range(numprim),element['exponents'],element['coefs'])
          ]
      outlines+=['  }','}']
    return '\n'.join(outlines)

  #------------------------------------------------
  def export_qwalk_slater(self,weights,states,orbfile,write_orbfile=True):
    ''' Export a Slater section for use in a trialfunction.
    It is far more efficient to write an orbfile and set the write_orbfile to False.
    In that case, pass the orbfile from that write.

    Args: 
      weights (array-like): Weights of determinants for multideterminant expansion. First should be 1.0.
      states (array-like): states[determinant][spin channel][orbital] select orbitals for determinants.
      orbfile (str): where orbitals are stored on disk (see write_qwalk_orb).
      write_orbfile (bool): Make a new orb file with these specifications. 
        If False, will assume such an orbfile already is set up.
    Returns:
      str: Slater wave function section for QWalk.
    '''
    # Convert python indexing to QWalk indexing.
    states=array(states)+1
    weights=array(weights)

    # Check input validity.
    assert (len(states.shape)==3) and (states.shape[1]==2) and (states.shape[0]==weights.shape[0]),\
        "States array should be nspin by ndeterminant by nelectrons. One detweight per determinant."

    if write_orbfile:
      self.write_qwalk_orb(orbfile)

    if any([(e.imag!=0.0).any() for e in self.eigvecs]):
      orbstr = "corbitals"
    else:
      orbstr = "orbitals"

    upstatelines = [' '.join(det[0]) for det in states.astype(str)]
    downstatelines = [' '.join(det[1]) for det in states.astype(str)]

    outlines = [
        "slater",
        "{0} {{".format(orbstr),
        "  magnify 1",
        "  nmo {0}".format(states.max()),
        "  orbfile {0}".format(orbfile),
        self.export_qwalk_basis(),
        "  centers { useglobal }",
        "}",
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

###############################################################################
def find_min_exp(basis):
  ''' Find minimum exponent in basis. '''
  allexponents=[]
  for species in basis:
    for basisel in basis[species]:
      allexponents+=list(basisel['exponents'])
  return min(allexponents)

###############################################################################
def count_naos(basis):
  ''' How many AOs are there in each atom?
  Args:
    basis (dict): basis dictionary like in Orbitals.
  Returns:
    int: number of AOs for each atomic species in the basis.
  '''
  countmap={'S':1,'P':3,'5D':5,'7F_crystal':7,'G':9,'H':11}
  results={}

  for atom in basis:
    results[atom]=0
    for basis_element in basis[atom]:
      results[atom]+=countmap[basis_element['angular']]
  return results
