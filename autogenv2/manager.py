import numpy as np
import os 
import pickle as pkl

######################################################################
class Manager:
  ''' Skeleton for managers.'''
  def __init__(name='AGmanager',path='./'):
    # Where to save self.
    self.name=name
    self.pickle="%s.pkl"%self.name

    # Ensure path is set up correctly.
    if path is None:
      path=os.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    self.logname="%s@%s"%(self.__class__.__name__,self.path+self.name)

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      print(self.logname,": rebooting old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      self.recover(old)

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #----------------------------------------
  def nextstep(self):
    ''' Redefine to have the manager do something.'''
    pass

  #----------------------------------------
  def collect(self):
    ''' Redefine to have manager farm data from its reader.'''
    pass

  #----------------------------------------
  def status(self):
    if self.completed:
      return 'ok'
    else:
      return 'not_finished'

  #------------------------------------------------
  def update_pickle(self):
    ''' If you make direct changes to the internals of the pickle, you need to call this to insure they are saved.'''
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #----------------------------------------
  def submit(self):
    ''' Submit any work and update the manager.'''
    qsubfile=self.runner.submit()

    self.update_pickle()

    return qsubfile

  #----------------------------------------
  def release_commands(self):
    ''' Release the runner of any commands it was tasked with and update the manager.'''
    commands=self.runner.release_commands()
    self.update_pickle()

    return commands

  #------------------------------------------------
  def update_queueid(self,qid):
    ''' If a bundler handles the submission, it can update the queue info with this.
    Args:
      qid (str): new queue id from submitting a job. The Manager will check if this is running.
    '''
    self.runner.queueid.append(qid)
    self.update_pickle()

  #----------------------------------------
  def export_results(self):
    return {
        'manager':self.__class__.__name__,
        'name':self.name,
        'path':self.path
      }

######################################################################
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
