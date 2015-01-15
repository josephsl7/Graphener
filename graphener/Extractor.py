'''
Created on Aug 27, 2014


'''
import os, subprocess

class Extractor:
    """ This class is responsible for creating "pseudo-POSCAR" files for each structure in the
        set of structures that we want to run through VASP calculations.  It does this by 
        "extracting" the information about the structure from struct_enum.out using the makestr.x
        routine from the enumlib library in UNCLE.  The StructsToPoscar class will then take this 
        set of pseudo-POSCARs and prepare them for VASP calculations. """

    def __init__(self, atoms, uncleOutput, startMethod,pureMetal):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.extractExec = os.path.abspath('needed_files/makestr.x')
        self.uncleOut = uncleOutput
        self.exStructList = []
        self.startMethod = startMethod
        self.pureMetal = pureMetal


    def checkPureInCurrent(self, iterNum, vstructsToStart, vstructsFinish):  #OK
        """ This checks that the pure elements are in the list to calculate or in the finished structures
        Only called on the first iteration for the firs """
#        self.exStructList = []
        for iatom in xrange(len(self.atoms)):
            if not self.contains('1', vstructsToStart[iatom]+vstructsFinish[iatom]):
                vstructsToStart[iatom].append ('1')
            if not self.contains(self.pureMetal, vstructsToStart[iatom]+vstructsFinish[iatom]):
                vstructsToStart[iatom].append (self.pureMetal)
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
        # Only extract the union of all the sets of structures. (No duplicates)
        for i in xrange(len(vstructsToStart)):
            uniqueSet = uniqueSet.union(vstructsToStart[i])
        for struct in uniqueSet:
            subprocess.call([self.extractExec, 'struct_enum.out', struct], stdout=self.uncleOut) 
        os.chdir(lastDir)
          

        

        