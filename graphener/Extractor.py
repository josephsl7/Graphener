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

    def __init__(self):
        """ CONSTRUCTOR """
        
        self.extractExec = os.path.abspath('needed_files/makestr.x')
        self.structList = []

    def setTrainingStructs(self):
        self.structList = []
        trainFile = open('enum/training_set_structures.dat','r')
        for line in trainFile:
            self.structList.append(line.strip().split()[1])
        
        trainFile.close()
    
    def setStructList(self, alist):
        self.structList = []
        for item in alist:
            self.structList.append(str(item))
    
    def getStructList(self):
        structs = []
        for struct in self.structList:
            structs.append(struct)
        
        structs.sort()
        return structs
          
    def extract(self):
        # You must call either the setTrainingStructs() or setStructList() functions before calling
        # this method.
        subprocess.call(['echo','\nExtracting structures from struct_enum.out'])
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        for struct in self.structList:
            subprocess.call([self.extractExec, 'struct_enum.out', struct])
        
        os.chdir(lastDir)
        

        