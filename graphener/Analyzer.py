'''
Created on Aug 27, 2014

@author: eswens13
'''
import os
import re
import subprocess


class Analyzer:


    def __init__(self, atoms):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.structList = []
        
        self.setStructList()
    
    def makeAnalysisDir(self):
        subprocess.call(['mkdir','analysis'])
    
    def setStructList(self):
        atomDir = os.getcwd() + '/' + self.atoms[0]
        contents = os.listdir(atomDir)
        for item in contents:
            itemDir = atomDir + '/' + item
            if os.path.isdir(itemDir):
                self.structList.append(item)
        
    def analyze(self):
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            for struct in self.structList:
                structDir = atomDir + '/' + struct + '/'
                self.writeAnalysis(structDir)
                if os.path.isdir(structDir + 'normal'):
                    self.writeAnalysis(structDir + 'normal/')
                if os.path.isdir(structDir + 'DOS'):
                    self.writeAnalysis(structDir + 'DOS/')

    def writeAnalysis(self, folder):
        # folder must end with a '/'
        outfile = open("results", "w")
        outfile.write("********** INCAR SETTINGS **********\n\n")
        incarLines = self.readfile(folder + "INCAR")
        for line in incarLines:
            outfile.write(line)
    
        outfile.write("\n********** FINISH CHECK **********\n")
        finished = self.FinishCheck(folder)
        if finished:
            outfile.write("\nFinished\n")
        else:
            outfile.write("\nDid not finish\n")

        if finished:
            outfile.write("\n********** CONVERGE CHECK **********\n")

            theLine = incarLines[0]
            for line in incarLines:
                match = re.search('.*NSW.*', line)
                if match:
                    theLine = line
                    break
    
            match = re.match('(\s*NSW\s*)=\s*([0-9]*)', theLine)
            if match:
                NSW = int(match.group(2))
                if self.convergeCheck(folder, NSW):
                    outfile.write("\nConverged in " + str(self.getSteps(folder)) + " steps.\n")
                else:
                    outfile.write("\nDid not converge.\n")
            else:
                outfile.write("\nBad INCAR format\n")
        
            outfile.write("\n********** CPU TIME **********\n")

            CPUTime = self.cpuTime(folder)
            CPUTime = CPUTime / 60.0
            hours = 0;
            minutes = 0;
            seconds = 0;
        
            if CPUTime > 1:
                if CPUTime > 60.0:
                    hours = int(CPUTime / 60.0)
                    hoursString = str(hours).zfill(2)
                
                    minutes = int(CPUTime) - (hours * 60)
                    minutesString = str(minutes).zfill(2)
                
                    seconds = int((CPUTime % 1) * 60)
                    secondsString = str(seconds).zfill(2)
            
                    outfile.write('\n' + hoursString + ':' + minutesString + ':' + secondsString + '\n')
            
                else:
                    minutes = int(CPUTime)
                    minutesString = str(minutes).zfill(2)
                
                    seconds = int((CPUTime % 1) * 60)
                    secondsString = str(seconds).zfill(2)
            
                    outfile.write('\n00:'+ minutesString + ':' + secondsString + '\n')
            else:
                outfile.write("\nLess than one minute.\n")
        
            outfile.write("\n********** ENERGY **********\n")
        
            outfile.write("\n" + str(self.writeEnergiesOszicar([folder])) + " eV\n")
    
        outfile.close()
    
        subprocess.call(['cp', 'results', folder + 'end_of_run_info'])
        subprocess.call(['rm', 'results'])

    def readfile(self, filepath):
        file1 = open(filepath,'r')
        lines = file1.readlines()
        file1.close()
        return lines

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
    
    def cpuTime(self, folder):
        lastfolder = os.getcwd()
        os.chdir(folder)
        proc = subprocess.Popen(['tail OUTCAR | grep User'],shell=True,stdout=subprocess.PIPE)
        result =  proc.communicate()
        os.chdir(lastfolder) 
        try:
            time =  float(result[0][-10:-3].strip()) #last few characters
        except:
            time = 0
        
        return time
    
    def writeEnergiesOszicar(self, alist):           
        for i in alist:
            try:
                oszicar = open(i+'/OSZICAR','r')
                energy = oszicar.readlines()[-1].split()[2]
                oszicar.close()
            except:
                energy = '0'
        
            return energy