'''
Created on Aug 27, 2014

@author: eswens13
'''
import os, subprocess

class Extractor:
    """ This class is responsible for creating "pseudo-POSCAR" files for each structure in the
        set of structures that we want to run through VASP calculations.  It does this by 
        "extracting" the information about the structure from struct_enum.out using the makestr.x
        routine from the enumlib library in UNCLE.  The StructsToPoscar class will then take this 
        set of pseudo-POSCARs and prepare them for VASP calculations. """

    def __init__(self, atoms, uncleOutput, startFromExisting):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.extractExec = os.path.abspath('needed_files/makestr.x')
        self.uncleOut = uncleOutput
        self.exStructList = []
        self.startFromExisting = startFromExisting


    def checkPureInCurrent(self, iterNum, vstructsToStart, vstructsFinish):  #OK
        """ This checks that the pure elements are in the list to calculate or in the finished structures
        Only called on the first iteration for the firs """
#        self.exStructList = []
        for iatom in xrange(len(self.atoms)):
            if not self.contains('1', vstructsToStart[iatom]+vstructsFinish[iatom]):
                vstructsToStart[iatom].append ('1')
            if not self.contains('3', vstructsToStart[iatom]+vstructsFinish[iatom]):
                vstructsToStart[iatom].append ('3')
        return vstructsToStart

    def contains(self, struct, alist):
        """ Returns True if 'alist' contains the item 'struct', False otherwise. """
        if len(alist) == 0:
            return False
        
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
        
        return False

    def extract(self,vstructsToStart):
        """ This method uses the makestr.x executable from the enumlib in UNCLE to 
            create the pseudo-POSCAR files for each structure in self.exStructList. These files are 
            generally called something like "vasp.000241" indicating the structure number in 
            struct_enum.out.  We only want to extract the union of all the lists in 
            self.exStructList. """
        subprocess.call(['echo','\nExtracting structures from struct_enum.out\n'])
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        uniqueSet = set()
        
        # Only extract the union of all the sets of structures.  (No duplicates)
        
        for i in xrange(len(vstructsToStart)):
            uniqueSet = uniqueSet.union(vstructsToStart[i])
        for struct in uniqueSet:
            subprocess.call([self.extractExec, 'struct_enum.out', struct], stdout=self.uncleOut)
        
        os.chdir(lastDir)
    
    def setStructsFromGSS(self, alist): # this flattens the list because the pseudo-POSCAR is not unique to an atom
        """ Sets the list of structures that we want to run through VASP calculations.  The list 
            being passed to this method (alist) will actually be a list of lists.  It will have a 
            structure list for each atom that still has not finished the main convergence loop. """
        self.exStructList = []
        for atomStructs in alist:
            self.exStructList.append(atomStructs)
    
#    def getexStructList(self): # this flattens the list because the pseudo-POSCAR is not unique to an atom 
#        """ Returns the list of structures for each atom that has not finished the convergence
#            loop. """
#        structs = []
#        for atomStructs in self.exStructList:
#            structs.append(atomStructs)
#        
#        return structs
          

        

        