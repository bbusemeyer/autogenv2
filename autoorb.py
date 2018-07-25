import crystal2qmc

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
    self.last_orbfile=''

  #------------------------------------------------------------------------------------------
  def find_min_exp(self):
    ''' Find minimum exponent in basis. '''
    allexponents=[]
    for species in self.basis:
      for basisel in self.basis[species]:
        allexponents+=list(basisel['exponents'])
    return min(allexponents)

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

  #----------------------------------------------------------------------------------------------
  def write_qwalk_orb(self,nperspin,outfn):
    ''' Generate a orb file for QWalk. 
    This just writes to the file because orbfile are necessarily separate in QWalk.

    Args:
      outfn (str): file to write to.
      nperspin (tuple): (Nup,Ndown) number of electrons in each spin channel.
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
    eigvec_flat = [self.eigvecs[s].ravel() for s in range(nspin)]
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

    self.last_orbfile=outfn.split('/')[-1]

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

