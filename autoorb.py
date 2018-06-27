
#################################################################################################
class Orbitals:
  ''' Class containing information defining a set of single-particle orbitals.
  - Basis set for the orbitals.
  - Coefficients of the orbitals in the basis.
  '''

  #----------------------------------------------------------------------------------------------
  def __init__(self):
    self.basis={}
    self.eigvecs={}
    self.eigvals={}
    self.is_complex={}
    self.atom_order=()
    self.kpt_weights={}

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
  def write_orb(self,kpt,outfn):#,maxmo_spin=-1):
    ''' Generate a orb file for QWalk. 
    This just writes to the file because the orbital sections are often quite large.

    Args:
      kpt (tuple): kpoint coordinate of the orbitals. See self.eigsys for options.
      outfn (str): file to write to.
      maxmo_spin (int): TODO Number of orbitals to print for each spin channel.
        You should set this to the maximum over the two channels.
    '''
    outf=open(outfn,'w')

    num_angular={'S':1,'P':3,'5D':5,'7F_crystal':7,'G':9,'H':11}

    nmo=len(self.eigvecs[(0,0,0)]['real'])*self.eigvecs[(0,0,0)]['real'][0].shape[0]

    coef_cnt=0
    for moidx in range(nmo):
      for atidx,atom in enumerate(self.atom_order):
        nao=sum(( num_angular[element['angular']] for element in self.basis[atom] ))
        for aoidx in range(nao):
          coef_cnt += 1
          outf.write(" {:5d} {:5d} {:5d} {:5d}\n"\
              .format(moidx+1,aoidx+1,atidx+1,coef_cnt))
    eigreal_flat = [e[0:,:].flatten() for e in self.eigvecs[kpt]['real']]
    eigimag_flat = [e[0:,:].flatten() for e in self.eigvecs[kpt]['imag']]
    print_cnt = 0
    outf.write("COEFFICIENTS\n")
    if self.is_complex[kpt]: #complex coefficients
      for eigr,eigi in zip(eigreal_flat,eigimag_flat):
        for r,i in zip(eigr,eigi):
          outf.write("({:<.12e},{:<.12e}) "\
              .format(r,i))
          print_cnt+=1
          if print_cnt%5==0: outf.write("\n")
    else: #Real coefficients
      for eigr in eigreal_flat:
        for r in eigr:
          outf.write("{:< 15.12e} ".format(r))
          print_cnt+=1
          if print_cnt%5==0: outf.write("\n")
    outf.close()
