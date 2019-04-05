import numpy as np
import subprocess as sub
import os

class Bundler:
  ''' Class for handling the bundling of several jobs of approximately the same 
  length, but possibly in different locations. ''' 
  def __init__(self,queue='normal',
                    walltime='48:00:00',
                    jobname='AGBundler',
                    npb=16,ppn=32,
                    mode='xe',
                    account='batr',
                    prefix=None,
                    postfix=None
                    ):
    ''' npb is the number of nodes desired per bundle. '''
    self.npb=npb
    self.ppn=ppn
    self.jobname=jobname
    self.mode=mode
    self.account=account
    self.queue=queue
    self.walltime=walltime
    if prefix is None: self.prefix=[]
    else:              self.prefix=prefix
    if postfix is None: self.postfix=[]
    else:               self.postfix=postfix
    self.queueid=[]

  def submit(self,mgrs,jobname=None):
    ''' Submit a list of managers in bundles.
    Args:
      mgrs (list): list of managers to submit.
      jobname (str): what will appear in qstat.
    '''
    print(self.__class__.__name__,"Submitting bundles of jobs.")
    if jobname is None: jobname=self.jobname

    assign=np.cumsum([mgr.runner.nn for mgr in mgrs])
    assign=((assign-0.1)//self.npb).astype(int)

    print(assign)

    for bidx in range(assign[-1]+1):
      self._submit_bundle(np.array(mgrs)[assign==bidx],"%s_%d"%(jobname,bidx))

  def _submit_bundle(self,mgrs,jobname=None,nn=None):
    ''' Submit a set of runners that require the correct number of nodes.
    This is usually called by submit, after it determines the break-up of jobs.
    Args: 
      mgrs (list): list of managers ready for submission. 
      jobname (str): what appears in qstat.
      nn (int): number of nodes to be used for all jobs (default:sum of nn in each manager).
    '''
    if nn is None:      nn=sum([mgr.runner.nn for mgr in mgrs])
    if jobname is None: jobname=self.jobname
    cwd=os.getcwd()

    qsublines=[
        "#PBS -q %s"%self.queue,
        "#PBS -l nodes=%i:ppn=%i:%s"%(nn,self.ppn,self.mode),
        "#PBS -l walltime=%s"%self.walltime,
        "#PBS -j oe ",
        "#PBS -A %s"%self.account,
        "#PBS -N %s "%jobname,
        "#PBS -o %s.out "%jobname,
        "cd %s"%cwd
      ] + self.prefix
    for mgr in mgrs:
      # This might be better without an error-out.
      lines=mgr.release_commands()
      if len(lines)>0:
        qsublines+=["cd %s"%mgr.path]+lines+["cd %s"%cwd]
    qsublines+=["wait"]+self.postfix

    qsubfile=jobname+".qsub"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsublines))
    try:
      result=sub.check_output("qsub %s"%(qsubfile),shell=True)
      queueid=result.decode().split()[0].split('.')[0]
      print("Submitted as %s"%queueid)
    except sub.CalledProcessError:
      print("Error submitting job. Check queue settings.")

    for mgr in mgrs:
      mgr.update_queueid(queueid)


