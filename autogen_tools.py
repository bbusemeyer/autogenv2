import numpy as np
import os 

def resolve_status(runner,reader,outfile):
  #Check if the reader is done
  if reader.completed:
    return 'done'

  #Check if the job is in the queue or running. If so, we just return that.
  currstat=runner.check_status()
  if currstat=='running':
    return currstat
  
  #Now we are in a state where either there was an error,
  #the job hasn't been run, or we haven't collected the results
  if not os.path.exists(outfile):
    return 'not_started'

  #We are in an error state or we haven't collected the results. 
  return "ready_for_analysis"

######################################################################
def deep_compare(d1,d2):
  '''I have to redo dict comparison because numpy will return a bool array when comparing.'''
  if type(d1)!=type(d2):
    return False
  if type(d1)==dict:
    if d1.keys()!=d2.keys():
      return False
    allsame=True
    for key in d1.keys():
      allsame=allsame and deep_compare(d1[key],d2[key])
    return allsame
  else:
    try:
      return np.array_equal(d1,d2)
    except TypeError:
      return d1==d2

######################################################################
def update_attributes(copyto,copyfrom,skip_keys=[],take_keys=[]):
  ''' Save update of class attributes. If copyfrom has additional attributes, they are ignored.

  Args:
    copyto (obj): class who's attributes are being updated.
    copyfrom (obj): class who's attributes will be copied from.
    skip_keys (list): list of attributes (str) not to update. 
    take_keys (list): list of attributes (str) that are ok to update. Others will raise warning and be skipped.
  Returns:
    bool: Whether any changes were made.
  '''
  updated=False
  for key in copyfrom.__dict__.keys():
    if key in skip_keys: 
      #print("Skipping key (%s)"%key)
      pass
    elif key not in copyto.__dict__.keys():
      print("Warning: Object update. An attribute (%s) was skipped because it doesn't exist in both objects."%key)
    elif not deep_compare(copyto.__dict__[key],copyfrom.__dict__[key]):
      if key not in take_keys:
        print("Warning: update to attribute (%s) cancelled, because it requires job to be rerun."%key)
      else:
        #print("Copy",key)
        copyto.__dict__[key]=copyfrom.__dict__[key]
        updated=True
    else:
      #print("Keys match (%s)"%key)
      pass
  return updated

######################################################################
def separate_jastrow(wffile,optimizebasis=True,freezeall=False):
  ''' Seperate the jastrow section of a QWalk wave function file.
  Args:
    optimizebasis (bool): whether to leave optimizebasis keyword in the section. 
    freezeall (bool): whether to add freeze to all sections with parameters. 
  Returns:
    str: Jastrow section for qwalk.
  '''
  wff=open(wffile,'r')
  tokens=wff.read().split('\n')
  in_jastrow=False
  nopen=0
  nclose=0
  jastlines=[]
  for line in tokens:

    if 'jastrow2' in line.lower() or 'jastrow3' in line.lower():
      in_jastrow=True

    if in_jastrow:

      if not optimizebasis and 'optimizebasis' in line.lower():
        line=' '.join([word for word in line.split() if word.lower() != 'optimizebasis'])
        if line=='': continue

      if 'BODY' in line and freezeall:
        line=line+' FREEZE'

      nopen+=line.count("{")
      nclose+=line.count("}")

      if nopen<nclose:
        in_jastrow=False
      else:
        jastlines.append(line)

  return '\n'.join(jastlines)

###############################################################################
# f orbital normalizations are from 
# <http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html>
def normalize_eigvec(eigvec,basis,atom_order):
  ''' Changes crystal normalization to qwalk normalization.
  eigvec will also be changed in place (i.e. does not copy).
  Args:
    eigvec (array): eigenvectors indexed by vector then AO.
    basis (dict): basis information. 
  Returns:
    array: normalized view of eigvec.
  '''
  snorm = 1./(4.*np.pi)**0.5
  pnorm = snorm*(3.)**.5
  dnorms = [
      .5*(5./(4*np.pi))**.5,
      (15./(4*np.pi))**.5,
      (15./(4*np.pi))**.5,
      .5*(15./(4.*np.pi))**.5,
      (15./(4*np.pi))**.5
    ]
  fnorms = [
      ( 7./(16.*np.pi))**.5,
      (21./(32.*np.pi))**.5,
      (21./(32.*np.pi))**.5,
      (105./(16.*np.pi))**.5, # xyz
      (105./(4.*np.pi))**.5,
      (35./(32.*np.pi))**.5,
      (35./(32.*np.pi))**.5
    ]

  # Find order of angular momentum channels.
  angular_order = []
  for species in atom_order:
    for element in basis[species]:
      angular_order.append(element['angular'])

  # This is to align properly with the d-components of eigvecs.
  dnorms *= angular_order.count('5D')
  fnorms *= angular_order.count('7F_crystal')

  # Duplicate according to number of Lz states to get AO labels for eigvecs.
  nangular = {'S':1,'P':3,'5D':5,'7F_crystal':7,'G':9,'H':11}
  ao_type=[]
  for btype in angular_order:
    ao_type+=nangular[btype]*[btype]

  ao_type=np.array(ao_type)

  eigvec[:,ao_type=='S'] *= snorm
  eigvec[:,ao_type=='P'] *= pnorm
  eigvec[:,ao_type=='5D'] *= dnorms
  eigvec[:,ao_type=='7F_crystal'] *= fnorms

  return eigvec
