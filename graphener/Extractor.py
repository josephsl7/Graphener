'''
Created on Aug 27, 2014

@author: eswens13
'''
import os, subprocess

class Extractor:
    """ This class is responsible for creating "pseudo-POSCAR" files for each structure in the
        set of structures that we want to run through VASP calculations.  It does this by 
        "extracting" the information about the structure from struct_enum.out using the makestr.x
        routine from the enumlib library in UNCLE.  The Structs2Poscar class will then take this 
        set of pseudo-POSCARs and prepare them for VASP calculations. """

    def __init__(self, atoms, uncleOutput, startFromExisting):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.extractExec = os.path.abspath('needed_files/makestr.x')
        self.uncleOut = uncleOutput
        self.exStructList = []
        self.startFromExisting = startFromExisting

    def setStructsFromTraining(self, iterNum, vstructsAll):
        """ This method sets the structure list for each atom to the list of training structures
            generated by UNCLE i.i.d. selection.  It DOES NOT choose the set of training structures,
            it just sets the member list of this class.  We always include the pure structures in 
            the fit, so it makes sure that they are done at least by the end of the first 
            iteration. Appends the chosen structures to past_structs.dat for each atom. """
        self.exStructList = []
        for i in xrange(len(self.atoms)):
            if self.startFromExisting[i] and iterNum == 1:
                self.exStructList.append([])
            else:
                trainStructs = []
                trainFile = open(self.atoms[i] + '/enumpast/training_set_structures.dat','r')
                for line in trainFile:
                    trainStructs.append(line.strip().split()[1])
                trainFile.close()
            
                # Make sure the pure structures are in the first iteration's training structures list
                # if they are not in the atom's past_structs.dat file.  We want to get the pure
                # structures in the fit on the first iteration because they are used to calculate
                # formation energies for everything else.
                if iterNum == 1:
                    if not self.contains('1', vstructsAll[i]):
                        if not self.contains('1', trainStructs):
                            trainStructs[0] = '1'
                    if not self.contains('3', vstructsAll[i]):
                        if not self.contains('3', trainStructs):
                            trainStructs[1] = '3'
#                # Append the newly chosen structures to past_structs.dat so they don't get chosen again
#                pastFile = open(self.atoms[i] + '/enumpast/past_structs.dat','a')
#                for struct in trainStructs:
#                    pastFile.write(struct + '\n')
#                pastFile.close()
            
                self.exStructList.append(trainStructs)

    def contains(self, struct, alist):
        """ Returns True if 'alist' contains the item 'struct', False otherwise. """
        if len(alist) == 0:
            return False
        
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
        
        return False
    
    def setStructsFromGSS(self, alist):
        """ Sets the list of structures that we want to run through VASP calculations.  The list 
            being passed to this method (alist) will actually be a list of lists.  It will have a 
            structure list for each atom that still has not finished the main convergence loop. """
        self.exStructList = []
        for atomStructs in alist:
            self.exStructList.append(atomStructs)
    
    def getexStructList(self):
        """ Returns the list of structures for each atom that has not finished the convergence
            loop. """
        structs = []
        for atomStructs in self.exStructList:
            structs.append(atomStructs)
        
        return structs
          
    def extract(self):
        """ You must call either the setTrainingStructs() or setexStructList() functions before calling
            this method.  This method uses the makestr.x executable from the enumlib in UNCLE to 
            create the pseudo-POSCAR files for each structure in self.exStructList. These files are 
            generally called something like "vasp.000241" indicating the structure number in 
            struct_enum.out.  We only want to extract the union of all the lists in 
            self.exStructList. """
        subprocess.call(['echo','\nExtracting structures from struct_enum.out\n'])
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        uniqueSet = set(self.exStructList[0])
        
        # Only extract the union of all the sets of structures.  (No duplicates)
        for i in xrange(1,len(self.exStructList)):
            uniqueSet = uniqueSet.union(self.exStructList[i])
        
        for struct in uniqueSet:
            subprocess.call([self.extractExec, 'struct_enum.out', struct], stdout=self.uncleOut)
        
        os.chdir(lastDir)
        

        