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
    return outlines

  #----------------------------------------------------------------------------------------------
  def write_qwalk_orb(self,outfn,nperspin,nvirtual=0,spin_restricted=False):
    ''' Generate a orb file for QWalk. 
    This just writes to the file because orbfile are necessarily separate in QWalk.

    Args:
      outfn (str): file to write to.
      nperspin (tuple): (Nup,Ndown) number of electrons in each spin channel.
      nvirtual (int): Number of virtual (unoccupied in ground state) orbitals to print for each spin channel.
      spin_restricted (bool): Up and down orbitals are restricted to be the same.
    '''
    # Figure out sizes of things.
    if spin_restricted: nspin=1
    else:               nspin=2

    outf=open(outfn,'w')

    nao_atom = count_naos(self.basis)
    totnmo = (max(nperspin)+nvirtual)*nspin

    # Do the printing.
    coef_cnt = 0
    for moidx in range(totnmo):
      for atidx,atom in enumerate(self.atom_order):
        for aoidx in range(nao_atom[atom]):
          outf.write(" {:5d} {:5d} {:5d} {:5d}\n"\
              .format(moidx+1,aoidx+1,atidx+1,coef_cnt+1))
          coef_cnt += 1
    eigvec_flat = [self.eigvecs[s][0:nperspin[s]+nvirtual].flatten() for s in range(nspin)] # TODO not memory efficient.
    print_cnt = 0
    outf.write("COEFFICIENTS\n")
    if (eigvec_flat[0].imag!=0.0).any(): # Complex coefficients
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


