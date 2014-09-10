'''
Created on Aug 27, 2014

@author: eswens13
'''
import os, subprocess

class Extractor:
    """ This class is responsible for using UNCLE to build the clusters necessary to generate a
        set of 'training' structures and creating pseudo-POSCAR files for each structure in the
        set. The Structs2Poscar class will take this set of pseudo-POSCARs and prepare them for 
        VASP calculations. """

    def __init__(self, clusterNums, trainStructNum):
        """ CONSTRUCTOR """
        
        self.clusterNums = clusterNums
        self.trainStructNum = trainStructNum
        
        self.uncleExec = os.path.abspath('needed_files/uncle.x')
        self.extractExec = os.path.abspath('needed_files/makestr.x')
        self.enumFile = 'enum/struct_enum.out'
        
        self.structList = []
    
    def changeEnumFile(self):
        subprocess.call(['mv',self.enumFile, self.enumFile + '_OLD'])
        
        oldfile = open(self.enumFile + '_OLD','r')
        oldlines = [line for line in oldfile]
        oldfile.close()
        
        newfile = open(self.enumFile, 'w')
        for i in xrange(len(oldlines)):
            if i == 1:
                newfile.write('bulk\n')
            else:
                newfile.write(oldlines[i])
        
        newfile.close()
    
    def buildClusters(self):
        oldLatFile = 'needed_files/lat.in'
        oldFile = open(oldLatFile, 'r')
        oldLines = [line for line in oldFile]
        oldFile.close()
        
        newFile = open('enum/lat.in','w')
        for i in xrange(len(oldLines)):
            if i == 38:
                for num in self.clusterNums:
                    newFile.write(str(num) + " ")
                newFile.write("\n")
            else:
                newFile.write(oldLines[i])
        newFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        
        proc = subprocess.Popen([self.uncleExec, '10'])
        print "\npid = " + str(proc.pid) + "\n"
        [procid, exitstatus] = os.waitpid(proc.pid, 0)
        print "Exit status from os.waitpid() = " + str(exitstatus)
        
        print "\nContinued? " + str(os.WIFCONTINUED(exitstatus))
        print "\nExited? " + str(os.WIFEXITED(exitstatus))
        print "\n\tExit Status was " + str(os.WEXITSTATUS(exitstatus))
        print "\nSignaled? " + str(os.WIFSIGNALED(exitstatus))
        print "\n\tSignal was " + str(os.WSTOPSIG(exitstatus))
        print "\nStopped? " + str(os.WIFSTOPPED(exitstatus))
        
        os.chdir(lastDir)

    def chooseTrainingStructures(self):
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        
        subprocess.call([self.uncleExec, '42', str(self.trainStructNum)])
        
        os.chdir(lastDir)

    def setStructList(self):
        trainFile = open('enum/training_set_structures.dat','r')
        for line in trainFile:
            self.structList.append(line.strip().split()[1])
        
        trainFile.close()
          
    def extract(self):
        #******************************************************************************************
        #  This is the issue.  The call to self.chooseTrainingStructures() is getting called before
        #  UNCLE is done building the clusters.
        #******************************************************************************************
        self.changeEnumFile()
        print "\nGenerating clusters. . .\n"
        self.buildClusters()
        print "\nChoosing i.i.d. structures. . .\n"
        self.chooseTrainingStructures()
        self.setStructList()
        
        print "\nExtracting structures from struct_enum.out"
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        for struct in self.structList:
            subprocess.call([self.extractExec, 'struct_enum.out', struct])
        
        os.chdir(lastDir)
        

        