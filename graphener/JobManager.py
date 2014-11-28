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
    
    def reportLowStats(self, vstructsCurrent):
        """ Displays the percentage of structures that converged during the low-precision VASP
            calculations.  Also displays the number of structures that converged and the number
            of structures that didn't. """
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for structure in vstructsCurrent[i]:
                    structureDir = atomDir + '/' + structure
                    if os.path.isdir(structureDir):
                        if self.FinishCheck(structureDir) and self.convergeCheck(structureDir, 400):
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
                            
    def FinishCheck(self, folder):
        """ Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR.  The input
            parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """
            
        lastfolder = os.getcwd()
        os.chdir(folder)
        
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)   
         
        return newstring[0].find('Voluntary') > -1 #True/False

    def convergeCheck(self, folder, NSW):
        """ Tests whether force convergence is done by whether the last line of OSZICAR (the last
            ionic relaxation step) is less than NSW."""
        try:
            value = self.getSteps(folder)
            return value < NSW #True/False
        except:
            return False #True/False

    def getSteps(self, folder):
        """Returns the number of steps of ionic relaxation, as an integer. """
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

    def reportDOSStats(self, vstructsCurrent):
        """ Reports the percentage of structures that converged during the Density of States VASP
            calculations.  Also reports the number of structures that converged and the number
            that didn't. """
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for struct in vstructsCurrent[i]:
                    structDir = atomDir + '/' + struct
                    if os.path.isdir(structDir):
                        dosDir = structDir + '/DOS'
                        if os.path.isdir(dosDir):
                            if self.FinishCheck(dosDir) and self.convergeCheck(dosDir, 2):
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


    def reportNormalStats(self, vstructsCurrent):
        """ Reports the percentage of structures that converged during normal-precision VASP
            calculations.  Also reports the number of structures that converged and the number
            that didn't. """
        for i in xrange(len(self.atoms)):
            subprocess.call(['echo','\nFor atom ' + self.atoms[i] + ':'])
            total = 0
            converged = 0
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                for struct in vstructsCurrent[i]:
                    structDir = atomDir + '/' + struct
                    if os.path.isdir(structDir):
                        normalDir = structDir + '/normal'
                        if os.path.isdir(normalDir):
                            if self.FinishCheck(normalDir) and self.convergeCheck(normalDir, 400):
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
    
    def runDOSJobs(self, vstructsCurrent):
        """ Starts the Density of States VASP calculations for all of the structures in 
            'vstructsCurrent' and waits for all of the jobs to finish. It checks on the jobs every ten 
            minutes. """
        subprocess.call(['echo','\nStarting DOS run. . .\n'])
        self.vaspRunner.run(3, vstructsCurrent)   
        s = scheduler(time.time, time.sleep)    
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 60 #seconds
            s.enterabs(event_time, 1, self.reportFinished, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinished(self.vaspRunner.getCurrentJobIds())   
        self.reportDOSStats(vstructsCurrent)

    def runHexMono(self):    
        subprocess.call(['echo','\nPreparing and running hexagonal metal monolayers directories\n']) #bch   
        self.vaspRunner.makeRunHexMono()  
            
    def runLowJobs(self, vstructsCurrent):
        """ Starts the low-precision VASP calculations for all of the structures in 'vstructsCurrent'
            and waits for all of the jobs to finish. It checks on the jobs every ten minutes. """
        subprocess.call(['echo','\nPreparing directories for VASP. . .\n'])
        self.vaspRunner.prepareForVasp(vstructsCurrent)
    
        s = scheduler(time.time, time.sleep)
        
        subprocess.call(['echo','\nStarting low-precision ionic relaxation. . .\n'])
        self.vaspRunner.run(1, vstructsCurrent)
        
        finished = False
        start_time = time.time()
        event_time = start_time
        while not finished:
            event_time += 20
            s.enterabs(event_time, 1, self.reportFinished, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinished(self.vaspRunner.getCurrentJobIds())
    
        self.reportLowStats(vstructsCurrent)
    
    def runNormalJobs(self, vstructsCurrent):
        """ Starts the normal-precision VASP calculations for all of the structures in 'vstructsCurrent'
            and waits for all of the jobs to finish. It checks on the jobs every ten minutes. """
        subprocess.call(['echo','\nStarting normal-precision ionic relaxation. . .\n'])
        print 'vstructsCurrent', vstructsCurrent
        self.vaspRunner.run(2, vstructsCurrent)
        
        s = scheduler(time.time, time.sleep)
    
        finished = False
        start_time = time.time()
        event_time = start_time
        print 'job ids'
        print self.vaspRunner.getCurrentJobIds()
        while not finished or len(self.vaspRunner.getCurrentJobIds())==0:
            event_time += 60
            s.enterabs(event_time, 1,self.reportFinished(self.vaspRunner.getCurrentJobIds())) #bch: was self.reportFinished, ([self.vaspRunner.getCurrentJobIds()]))
            s.run()
            finished = self.reportFinished(self.vaspRunner.getCurrentJobIds())
    
        self.reportNormalStats(vstructsCurrent)
        print 'done with runNormalJobs'

    def runSingleAtoms(self):  
        subprocess.call(['echo','\nPreparing and running single atom directories\n']) #bch   
        self.vaspRunner.makeRunSingleDirectories()  
