'''



'''
import os, subprocess, time, sys
from sched import scheduler
from comMethods import *

class clustersjob:
    """ The largest numbers of clusters (up to 900*5) take 64G of memory and a single processor.  
    We need to run uncle 10 (build clusters) with a separate job on a large memory job
    that builds the clusters"""

    def __init__(self):
        """  """
        self.currJobIds = []

    def clustBuild(self):
        '''Runs the process of building clusters'''
        self.makeRunClusters()
        self.waitclustersjob()
        if not os.path.exists('clustersjob.out'):
            sys.exit('Cluster job never started.  Stopping program')
        else:
            lines = readfile('clustersjob.out')
            if 'done' not in lines[-1]:
                sys.exit('Cluster job failed.  Stopping program')            
            
    def waitclustersjob(self):
        """Waits until the cluster job is done """
        s = scheduler(time.time, time.sleep)    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 10 #check every x seconds
            s.enterabs(event_time, 1, self.reportFinished, (self.currJobIds))
            s.run()
            finished = self.reportFinished(self.currJobIds)
           
    def makeRunClusters(self): #bch all
        """Creates cluster jobfile and starts the run """     
        jobFile = open('clustersjob','w')
        jobFile.write("#!/bin/bash\n\n")
        jobFile.write("#SBATCH --time=00:35:00\n")
        jobFile.write("#SBATCH --ntasks=1\n")
        jobFile.write("#SBATCH --mem-per-cpu=64G\n")
        jobFile.write("#SBATCH --mail-user=joseph-lawrence@hotmail.com\n")              
        jobFile.write("#SBATCH --mail-type=FAIL\n")
        #jobFile.write("#SBATCH --mail-type=END\n") 
        jobFile.write("#SBATCH --job-name=clusters\n" )           
        jobFile.write("uncle 10 > clustersjob.out\n") 
        jobFile.close()
        proc = subprocess.Popen(['sbatch','clustersjob'], stdout=subprocess.PIPE)
        jobid = proc.communicate()[0].split()[3]
        subprocess.call(['echo', 'Submitted cluster job ' + jobid])
        self.currJobIds.append(jobid)
    
    def reportFinished(self, jobIds):
        """ Returns true if all the jobs with IDs given in 'jobIds' have finished VASP
            calculations, false otherwise. """
        devnull = open(os.devnull, 'w')
        for jobid in jobIds:
            proc = subprocess.Popen(['squeue', '--job', jobid], stdout=subprocess.PIPE, stderr=devnull)
            output = proc.communicate()[0].split()
            if len(output) != 8 and len(output) != 0:   # It will list either all the headers or
                return False                            # sometimes an error and outputs nothing.
                                                        # The error in this case is an "invalid
                                                        # job id" error because the job is no
        return True                                     # longer on the supercomputer.    
