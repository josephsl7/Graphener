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

    def __init__(self, atoms, uncleOutput):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.extractExec = os.path.abspath('needed_files/makestr.x')
        self.uncleOut = uncleOutput
        self.structList = []

    def setTrainingStructs(self):
        self.structList = []
        trainStructs = []
        trainFile = open('enum/training_set_structures.dat','r')
        for line in trainFile:
            trainStructs.append(line.strip().split()[1])
        trainFile.close()
        
        # Make sure the pure structures are in the training structures list.
        if not self.contains('1', trainStructs):
            trainStructs[0] = '1'
            
        if not self.contains('3', trainStructs):
            trainStructs[1] = '3'
        
        for atom in self.atoms:
            self.structList.append(trainStructs)
    
    def contains(self, struct, alist):
        if len(alist) == 0:
            return False
        
        for i in xrange(len(alist)):
            if struct == alist[i]:
                return True
        
        return False
    
    def setStructList(self, alist):
        """ The list being passed to this method (alist) will actually be a list of lists.  It will
            have a structure list for each atom that still has not finished the main convergence
            loop. """
        self.structList = []
        for atomStructs in alist:
            self.structList.append(atomStructs)
    
    def getStructList(self):
        structs = []
        for atomStructs in self.structList:
            structs.append(atomStructs)
        
        return structs
          
    def extract(self):
        # You must call either the setTrainingStructs() or setStructList() functions before calling
        # this method.
        subprocess.call(['echo','\nExtracting structures from struct_enum.out'])
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        uniqueSet = set(self.structList[0])
        
        # Only extract the union of all the sets of structures.  (No duplicates)
        for i in xrange(1,len(self.structList)):
            uniqueSet = uniqueSet.union(self.structList[i])
        
        for struct in uniqueSet:
            subprocess.call([self.extractExec, 'struct_enum.out', struct], stdout=self.uncleOut)
        
        os.chdir(lastDir)
        

        