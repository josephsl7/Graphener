'''
Created on Aug 26, 2014

@author: eswens13
'''
import os, subprocess, time
from sched import scheduler
import RunVasp

class JobManager:
    """ This class is responsible for scheduling and watching the VASP calculations.  It starts
        the low-precision relaxation, waits for all structures to either timeout or complete, 
        reports how many structures converged and how many did not, then starts the 
        normal-precision relaxation for all structures that converged.  It repeats this process
        for the normal-precision relaxation and the DOS runs.  An instance of the RunVasp class
        is used to prepare all the directories, retrieve the needed files, and submit the jobs to
        the supercomputer. """

    def __init__(self, atoms):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.vaspRunner = RunVasp.RunVasp(self.atoms)
    
    def reportFinshed(self, jobIds):
        for jobid in jobIds:
            proc = subprocess.Popen(['squeue', '--job', jobid], stdout=subprocess.PIPE)
            output = proc.communicate()[0].split()
            if len(output) != 8 and len(output) != 0:   # It will list either all the headers or
                return False                            # sometimes it outputs an error and outputs
                                                        # nothing.
        return True
    
    def reportLowStats(self):
        for atom in self.atoms:
            print "\nFor atom " + atom + ":"
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                contents = os.listdir(atomDir)
                for item in contents:
                    itemDir = atomDir + '/' + item
                    if os.path.isdir(itemDir):
                        if self.FinishCheck(itemDir):
                            total += 1
                            converged += 1
                        else:
                            total += 1
            
            percent = float(float(converged) / float(total)) * 100.0
            notConverged = total - converged
            print "\t%5.2f %% of the structures converged during low-precision relaxation." % percent
            print "\t%d structures converged." % converged
            print "\t%d structures did not converge." % notConverged
                            
    def FinishCheck(self, folder):
        """ Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR.  The input
            parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """
            
        lastfolder = os.getcwd()
        os.chdir(folder)
        
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)   
         
        return newstring[0].find('Voluntary') > -1 #True/False

    def reportNormalStats(self):
        for atom in self.atoms:
            print "\nFor atom " + atom + ":"
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                contents = os.listdir(atomDir)
                for item in contents:
                    itemDir = atomDir + '/' + item
                    if os.path.isdir(itemDir):
                        normalDir = itemDir + '/normal'
                        if os.path.isdir(normalDir):
                            if self.FinishCheck(normalDir):
                                total += 1
                                converged += 1
                            else:
                                total += 1
            
            percent = float(float(converged) / float(total)) * 100.0
            notConverged = total - converged
            print "\t%5.2f %% of the structures converged during normal-precision relaxation." % percent
            print "\t%d structures converged." % converged
            print "\t%d structures did not converge." % notConverged
    
    def reportDOSStats(self):
        for atom in self.atoms:
            print "\nFor atom " + atom + ":"
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                contents = os.listdir(atomDir)
                for item in contents:
                    itemDir = atomDir + '/' + item
                    if os.path.isdir(itemDir):
                        dosDir = itemDir + '/DOS'
                        if os.path.isdir(dosDir):
                            if self.FinishCheck(dosDir):
                                total += 1
                                converged += 1
                            else:
                                total += 1
            
            percent = float(float(converged) / float(total)) * 100.0
            notConverged = total - converged
            print "\t%5.2f %% of the structures converged during the DOS run." % percent
            print "\t%d structures converged." % converged
            print "\t%d structures did not converge." % notConverged
    
    def runLowJobs(self):
        print "\nPreparing directories for VASP. . .\n"
        self.vaspRunner.prepareForVasp()
    
        s = scheduler(time.time, time.sleep)
    
        print "\nStarting low-precision ionic relaxation. . .\n"
        self.vaspRunner.run(1)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 600
            s.enterabs(event_time, 1, self.reportFinshed, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinshed(self.vaspRunner.getCurrentJobIds())
    
        self.reportLowStats()
    
    def runNormalJobs(self):
        print "\nStarting normal-precision ionic relaxation. . .\n"
        self.vaspRunner.run(2)
        
        s = scheduler(time.time, time.sleep)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 600
            s.enterabs(event_time, 1, self.reportFinshed, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinshed(self.vaspRunner.getCurrentJobIds())
    
        self.reportNormalStats()

    def runDOSJobs(self):
        print "\nStarting DOS run. . .\n"
        self.vaspRunner.run(3)
    
        s = scheduler(time.time, time.sleep)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 600
            s.enterabs(event_time, 1, self.reportFinshed, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinshed(self.vaspRunner.getCurrentJobIds())
    
        self.reportDOSStats()









