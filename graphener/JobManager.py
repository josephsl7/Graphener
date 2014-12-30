'''
Created on Aug 26, 2014

@author: eswens13
'''
import os, subprocess, time
from sched import scheduler
import RunVasp
from comMethods import *


class JobManager:
    """ This class is responsible for scheduling and watching the VASP calculations.  It starts
        the low-precision relaxation, waits for all structures to either timeout or complete, 
        reports how many structures converged and how many did not, then starts the 
        normal-precision relaxation for all structures that converged.  It repeats this process
        for the normal-precision relaxation and the DOS runs.  An instance of the RunVasp class
        is used to prepare all the directories, retrieve the needed files, and submit the jobs to
        the supercomputer. """

    def __init__(self, atoms,ediffg):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.vaspRunner = RunVasp.RunVasp(self.atoms,ediffg)
    
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
    
    def reportLowStats(self, vstructsToStart):
        """ Displays the percentage of structures that converged during the low-precision VASP
            calculations.  Also displays the number of structures that converged and the number
            of structures that didn't. """
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for structure in vstructsToStart[i]:
                    structureDir = atomDir + '/' + structure
                    if os.path.isdir(structureDir):
#                        print structure ,'checks',finishCheck(structureDir), convergeCheck(structureDir, getNSW(structureDir))
                        if finishCheck(structureDir) and convergeCheck(structureDir, getNSW(structureDir)):
                            total += 1
                            converged += 1
                        else:
                            total += 1
            
            if total == 0:
                subprocess.call(['echo','\tNo structures were submitted for ' + self.atoms[i]])
            else:
                percent = float(float(converged) / float(total)) * 100.0
                notConverged = total - converged
                subprocess.call(['echo','\t%5.2f %% of the structures converged during low-precision relaxation.' % percent])
                subprocess.call(['echo','\t%d structures converged.' % converged])
                subprocess.call(['echo','\t%d structures did not converge.' % notConverged])
                            
    def reportDOSStats(self, vstructsToStart):
        """ Reports the percentage of structures that converged during the Density of States VASP
            calculations.  Also reports the number of structures that converged and the number
            that didn't. """
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for struct in vstructsToStart[i]:
                    structDir = atomDir + '/' + struct
                    if os.path.isdir(structDir):
                        dosDir = structDir + '/DOS'
                        if os.path.isdir(dosDir):
                            if finishCheck(dosDir) and convergeCheck(dosDir, 2):
                                total += 1
                                converged += 1
                            else:
                                total += 1            
            if total == 0:
                subprocess.call(['echo','\tNo structures were submitted for ' + self.atoms[i]])
            else:
                percent = float(float(converged) / float(total)) * 100.0
                notConverged = total - converged
                subprocess.call(['echo','\t%5.2f %% of the structures converged during the DOS run.' % percent])
                subprocess.call(['echo','\t%d structures converged.' % converged])
                subprocess.call(['echo','\t%d structures did not converge.' % notConverged]) 


    def reportNormalStats(self, vstructsToStart):
        """ Reports the percentage of structures that converged during normal-precision VASP
            calculations.  Also reports the number of structures that converged and the number
            that didn't. """
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for struct in vstructsToStart[i]:
                    structDir = atomDir + '/' + struct
                    if os.path.isdir(structDir):
                        normalDir = structDir + '/normal'
                        if os.path.isdir(normalDir):
                            if finishCheck(normalDir) and convergeCheck(normalDir, getNSW(normalDir)):
                                total += 1
                                converged += 1
                            else:
                                total += 1
            if total == 0:
                subprocess.call(['echo','\tNo structures were submitted for ' + self.atoms[i]])
            else:
                percent = float(float(converged) / float(total)) * 100.0
                notConverged = total - converged
                subprocess.call(['echo','\t%5.2f %% of the structures converged during normal-precision relaxation.' % percent])
                subprocess.call(['echo','\t%d structures converged.' % converged])
                subprocess.call(['echo','\t%d structures did not converge.' % notConverged])
    
    def runDOSJobs(self, vstructsToStart, vstructsRestart):
        """ Starts the Density of States VASP calculations for all of the structures in 
            'vstructsToStart' and waits for all of the jobs to finish. It checks on the jobs every ten 
            minutes. """
        subprocess.call(['echo','\nStarting DOS run. . .\n'])
        vstructsToRun = joinLists([vstructsRestart,vstructsToStart])
        self.vaspRunner.run(3, vstructsToStart, vstructsToRun)   
        s = scheduler(time.time, time.sleep)    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 60 #seconds
            s.enterabs(event_time, 1, self.reportFinished, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinished(self.vaspRunner.getCurrentJobIds())   
        self.reportDOSStats(vstructsToRun)

    def runHexMono(self):    
        subprocess.call(['echo','\nPreparing and running hexagonal metal monolayers directories\n']) #bch   
        self.vaspRunner.makeRunHexMono()  
            
    def runLowJobs(self, vstructsToStart,vstructsRestart):
        """ Starts the low-precision VASP calculations for all of the structures in 'vstructsToStart'
            and waits for all of the jobs to finish. It checks on the jobs every ten minutes. """
        subprocess.call(['echo','\nPreparing directories for VASP. . .\n'])
        vstructsToRun = joinLists([vstructsRestart,vstructsToStart])
        self.vaspRunner.prepareForVasp(vstructsToStart)
        self.vaspRunner.prepareRestarts(vstructsRestart) #contcar->poscar, incar and job file
    
        s = scheduler(time.time, time.sleep)
        
        if len(flat(vstructsToRun)) != 0:
            subprocess.call(['echo','\nStarting low-precision ionic relaxation. . .\n'])
            self.vaspRunner.run(1,vstructsToStart,vstructsToRun)
        
            finished = False
            start_time = time.time()
            event_time = start_time
            while not finished:
                event_time += 5
                s.enterabs(event_time, 1, self.reportFinished, ([self.vaspRunner.getCurrentJobIds()]))
                s.run()
                finished = self.reportFinished(self.vaspRunner.getCurrentJobIds())
    
#        self.reportLowStats(vstructsToRun)
    
    def runNormalJobs(self, vstructsToStart,vstructsRestart):
        """ Starts the normal-precision VASP calculations for all of the structures in 'vstructsToStart'
            and waits for all of the jobs to finish. It checks on the jobs every ten minutes. """
        subprocess.call(['echo','\nStarting normal-precision ionic relaxation. . .\n'])
        self.vaspRunner.run(2, vstructsToStart,vstructsToRun)
        vstructsToRun = joinLists([vstructsRestart,vstructsToStart])
        
        s = scheduler(time.time, time.sleep)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished or len(self.vaspRunner.getCurrentJobIds())==0:
            event_time += 60
            s.enterabs(event_time, 1,self.reportFinished(self.vaspRunner.getCurrentJobIds())) #bch: was self.reportFinished, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinished(self.vaspRunner.getCurrentJobIds())
    
        self.reportNormalStats(vstructsToRun)

    def runSingleAtoms(self):  
        subprocess.call(['echo','\nPreparing and running single atom directories\n']) #bch   
        self.vaspRunner.makeRunSingleDirectories()
#    
#    def runPureSys(self):  
#        subprocess.call(['echo','\nPreparing and running pure systems\n']) #bch   
#        self.vaspRunner.makeRunSingleDirectories()    
