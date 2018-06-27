
#################################################################################################
class Orbitals:
  ''' Class containing information defining a set of single-particle orbitals.
  - Basis set for the orbitals.
  - Coefficients of the orbitals in the basis.
  '''

  #----------------------------------------------------------------------------------------------
  def __init__(self):
    self.basis={}
    self.eigsys={}

  #----------------------------------------------------------------------------------------------
  def export_qwalk_basis(self):
    ''' Generate a basis section for QWalk.

    Args:
    Returns:
      list: lines pertaining to system section.
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
  def write_orb(self,kpt,outfn,maxmo_spin=-1):
    ''' Generate a orb file for QWalk. 
    This just writes to the file because the orbital sections are often quite large.

    Args:
      kpt (tuple): kpoint coordinate of the orbitals. See self.eigsys for options.
      outfn (str): file to write to.
      maxmo_spin (int): Number of orbitals to print for each both spin channels.
        You should set this to the maximum over the two channels.
    '''
    outf=open(outfn,'w')
    if maxmo_spin < 0:
      maxmo_spin=basis['nmo']

    eigvecs_real = self.eigsys['eigvecs'][kpt]['real']
    eigvecs_imag = self.eigsys['eigvecs'][kpt]['imag']
    atidxs = np.unique(basis['atom_shell'])-1
    nao_atom = np.zeros(atidxs.size,dtype=int)
    for shidx in range(len(basis['nao_shell'])):
      nao_atom[basis['atom_shell'][shidx]-1] += basis['nao_shell'][shidx]
    #nao_atom = int(round(sum(basis['nao_shell']) / len(ions['positions'])))
    coef_cnt = 1
    totnmo = maxmo_spin*self.eigsys['nspin'] #basis['nmo'] * eigsys['nspin']
    for moidx in np.arange(totnmo)+1:
      for atidx in atidxs+1:
        for aoidx in np.arange(nao_atom[atidx-1])+1:
          outf.write(" {:5d} {:5d} {:5d} {:5d}\n"\
              .format(moidx,aoidx,atidx,coef_cnt))
          coef_cnt += 1
    eigreal_flat = [e[0:maxmo_spin,:].flatten() for e in eigvecs_real]
    eigimag_flat = [e[0:maxmo_spin,:].flatten() for e in eigvecs_imag]
    print_cnt = 0
    outf.write("COEFFICIENTS\n")
    if self.eigsys['ikpt_iscmpx'][kpt]: #complex coefficients
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
