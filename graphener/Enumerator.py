'''
Created on Aug 26, 2014

@author: eswens13
'''
import os, subprocess


class Enumerator:
    """ Enumerates symmetrically unique structures using UNCLE.  After this class finishes its 
        work, the Extractor class will take over and extract pseudo-POSCAR files from the
        struct_enum.out file. """
  
    def __init__(self, atoms, volRange, clusterNums, trainStructNum, uncleOutput):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.volRange = volRange
        
        self.clusterNums = clusterNums
        self.trainStructNum = trainStructNum
        
        self.uncleExec = os.path.abspath('needed_files/uncle.x')
        self.enumFile = 'enum/struct_enum.out'
        self.enumExec = os.path.abspath('needed_files/enum.x')
        self.uncleOut = uncleOutput

    def changeEnumFile(self):
        """ In order to build the clusters that will be used in the cluster expansion correctly, 
            we have to change the 'surf' setting in struct_enum.out (from UNCLE enumeration) to 
            'bulk'.  This method takes care of that. """
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
        """ Uses UNCLE to build the number of each n-body clusters specified in the settings.in
            file. """
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
        
        subprocess.call([self.uncleExec, '10'], stdout=self.uncleOut)
        
        os.chdir(lastDir)

    def chooseTrainingStructures(self):
        """ Chooses a list of i.i.d. structures from struct_enum.out. The length of the list 
            is determined by the TRAINING_STRUCTS setting in settings.in. """
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        
        subprocess.call([self.uncleExec, '42', str(self.trainStructNum)], stdout=self.uncleOut)
        
        os.chdir(lastDir)
    
    def enumerate(self):
        """ Runs through the whole process of enumeration, cluster building, and choosing an
            i.i.d. set of training structures. """
        subprocess.call(['mkdir','enum'])
        infile = open('needed_files/struct_enum.in','r')
        inlines = []
        for line in infile:
            inlines.append(line)
        infile.close()
        
        structFile = open('enum/struct_enum.in','w')
        for i in xrange(len(inlines)):
            if i == 9:
                structFile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + " ")
                structFile.write("# Starting and ending cell sizes for search\n")
            else:
                structFile.write(inlines[i])
        structFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        subprocess.call([self.enumExec,'struct_enum.in'], stdout=self.uncleOut)
        
        os.chdir(lastDir)
        
        self.changeEnumFile()
        subprocess.call(['echo','\nGenerating clusters. . .\n'])
        self.buildClusters()
        subprocess.call(['echo','\nChoosing i.i.d. structures. . .\n'])
        self.chooseTrainingStructures()
        
            
            
        