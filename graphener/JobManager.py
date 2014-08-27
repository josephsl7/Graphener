'''
Created on Aug 26, 2014

@author: eswens13
'''
import os
import subprocess


class JobManager:


    def __init__(self, atoms):
        """ CONSTRUCTOR """
        self.atoms = atoms
    
    def reportFinshed(self, jobIds):
        for jobid in jobIds:
            proc = subprocess.Popen(['squeue', '--job', jobid], stdout=subprocess.PIPE)
            output = proc.communicate()[0].split()
            if len(output) != 8:
                return False
        
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
        
        