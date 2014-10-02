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
        devnull = open(os.devnull, 'w')
        for jobid in jobIds:
            proc = subprocess.Popen(['squeue', '--job', jobid], stdout=subprocess.PIPE, stderr=devnull)
            output = proc.communicate()[0].split()
            if len(output) != 8 and len(output) != 0:   # It will list either all the headers or
                return False                            # sometimes an error and outputs nothing.
                                                        # The error in this case is an "invalid
                                                        # job id" error because the job is no
        return True                                     # longer on the supercomputer.    
    
    def reportLowStats(self, structList):
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for structure in structList[i]:
                    structureDir = atomDir + '/' + structure
                    if os.path.isdir(structureDir):
                        if self.FinishCheck(structureDir) and self.convergeCheck(structureDir, 400):
                            total += 1
                            converged += 1
                        else:
                            total += 1
            
            percent = float(float(converged) / float(total)) * 100.0
            notConverged = total - converged
            subprocess.call(['echo','\t%5.2f %% of the structures converged during low-precision relaxation.' % percent])
            subprocess.call(['echo','\t%d structures converged.' % converged])
            subprocess.call(['echo','\t%d structures did not converge.' % notConverged])
                            
    def FinishCheck(self, folder):
        # From Dr. Hess
        """ Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR.  The input
            parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """
            
        lastfolder = os.getcwd()
        os.chdir(folder)
        
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)   
         
        return newstring[0].find('Voluntary') > -1 #True/False

    def convergeCheck(self, folder, NSW):
        """Tests whether force convergence is done by whether the last line of Oszicar is less than NSW."""
        try:
            value = self.getSteps(folder)
            return value < NSW #True/False
        except:
            return False #True/False

    def getSteps(self, folder):
        '''number of steps in relaxation, as an integer'''
        lastfolder = os.getcwd()
        os.chdir(folder)
        if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
            os.chdir(lastfolder) 
            return -9999
        oszicar = open('OSZICAR','r')
        laststep = oszicar.readlines()[-1].split()[0]
        oszicar.close()
        os.chdir(lastfolder)  
        try:
            value = int(laststep)
            return value
        except:
            return 9999

    def reportNormalStats(self, structList):
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for struct in structList[i]:
                    structDir = atomDir + '/' + struct
                    if os.path.isdir(structDir):
                        normalDir = structDir + '/normal'
                        if os.path.isdir(normalDir):
                            if self.FinishCheck(normalDir) and self.convergeCheck(normalDir, 400):
                                total += 1
                                converged += 1
                            else:
                                total += 1
            
            percent = float(float(converged) / float(total)) * 100.0
            notConverged = total - converged
            subprocess.call(['echo','\t%5.2f %% of the structures converged during normal-precision relaxation.' % percent])
            subprocess.call(['echo','\t%d structures converged.' % converged])
            subprocess.call(['echo','\t%d structures did not converge.' % notConverged])
    
    def reportDOSStats(self, structList):
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for struct in structList[i]:
                    structDir = atomDir + '/' + struct
                    if os.path.isdir(structDir):
                        dosDir = structDir + '/DOS'
                        if os.path.isdir(dosDir):
                            if self.FinishCheck(dosDir) and self.convergeCheck(dosDir, 2):
                                total += 1
                                converged += 1
                            else:
                                total += 1
            
            percent = float(float(converged) / float(total)) * 100.0
            notConverged = total - converged
            subprocess.call(['echo','\t%5.2f %% of the structures converged during the DOS run.' % percent])
            subprocess.call(['echo','\t%d structures converged.' % converged])
            subprocess.call(['echo','\t%d structures did not converge.' % notConverged])
    
    def runLowJobs(self, structList):
        subprocess.call(['echo','\nPreparing directories for VASP. . .\n'])
        self.vaspRunner.prepareForVasp(structList)
    
        s = scheduler(time.time, time.sleep)
    
        subprocess.call(['echo','\nStarting low-precision ionic relaxation. . .\n'])
        self.vaspRunner.run(1, structList)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 600
            s.enterabs(event_time, 1, self.reportFinshed, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinshed(self.vaspRunner.getCurrentJobIds())
    
        self.reportLowStats(structList)
    
    def runNormalJobs(self, structList):
        subprocess.call(['echo','\nStarting normal-precision ionic relaxation. . .\n'])
        self.vaspRunner.run(2, structList)
        
        s = scheduler(time.time, time.sleep)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 600
            s.enterabs(event_time, 1, self.reportFinshed, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinshed(self.vaspRunner.getCurrentJobIds())
    
        self.reportNormalStats(structList)

    def runDOSJobs(self, structList):
        subprocess.call(['echo','\nStarting DOS run. . .\n'])
        self.vaspRunner.run(3, structList)
    
        s = scheduler(time.time, time.sleep)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 600
            s.enterabs(event_time, 1, self.reportFinshed, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinshed(self.vaspRunner.getCurrentJobIds())
    
        self.reportDOSStats(structList)









