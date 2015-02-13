'''
'''
import os, subprocess
import ClustersBuild
from comMethods import *

class Enumerator:
    """ This class enumerates symmetrically unique structures in a given volume range using UNCLE.  
        It then builds the clusters necessary to perform a cluster expansion and chooses a 
        specified number of "training structures" to perform a first fit on.  After this class 
        finishes its work, the Extractor class will take over and extract pseudo-POSCAR files from 
        the struct_enum.out file that is produced. The methods in this class are only needed for 
        the first iteration of the main convergence loop. """
  
    def __init__(self, atoms, volRange, clusterNums, uncleOutput):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.volRange = volRange
        
        self.clusterNums = clusterNums
        
        self.uncleExec = os.path.abspath('needed_files/uncle.x')
        self.enumFile = os.path.abspath('enum/struct_enum.out')
        self.enumExec = os.path.abspath('needed_files/enum.x')
        self.uncleOut = uncleOutput
        self.clusterNums = clusterNums
        self.makeAtomDirectories()

    def buildClusters(self):
        """ Uses UNCLE to build the number of each n-body clusters specified in the settings.in
            file. """
        oldLatFile = 'needed_files/lat.in'
        oldFile = open(oldLatFile, 'r')
        oldLines = [line for line in oldFile]
        oldFile.close()
        
        newFile = open('enum/lat.in','w')
        for i in xrange(len(oldLines)):
            if 'Number pairs' in oldLines[i-1] and i>=1: #bch use label on previous line
                for num in self.clusterNums:
                    newFile.write(str(num) + " ")
                newFile.write("\n")
            else:
                newFile.write(oldLines[i])
        newFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
#        subprocess.call([self.uncleExec, '10'], stdout=self.uncleOut)
        if sum(self.clusterNums)<=1500: #the 1500 rule of thumb assumes you are running Main with 16G. 
            subprocess.call([self.uncleExec, '10'], stdout=self.uncleOut)
        else:
#            subprocess.call(['echo','Warning: BLOCKING CLUSTER JOB to save time'])
            clustersjob = ClustersBuild.clustersjob()
            clustersjob.clustBuild()
        os.chdir(lastDir)

    def changeEnumFile(self):
        """ In order to build the clusters that will be used in the cluster expansion correctly, 
            we have to change the 'surf' setting in struct_enum.out (from UNCLE enumeration) to 
            'bulk'.  It changes the name of the old 'surf' version to 'struct_enum.out_OLD'. """
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
        
    def chooseTrainingStructures(self,iteration, startMethod,nNew,ntot,distribute):
        """ If startMethod is not the same for each atomChooses a list of i.i.d. structures from struct_enum.out for each different metal atom,
         
            The UNCLE option that we run to choose the training structures should look for a file 
            called 'past_structs.dat' so as not to choose structures from that list. (Dr. Hess 
            added that option in the UNCLE executable we're using.) The length of the list of
            training structures for each atom is determined by the TRAINING_STRUCTS setting in 
            settings.in. """
        lastDir = os.getcwd()
        natoms = len(self.atoms)
        iidStructs = [[]]*natoms
        
#        subprocess.call(['echo','Warning: BLOCKING IID selection  to save time'])
#            for iatom,atom in enumerate(self.atoms):
#                if 2 <= nNew[iatom] < ntot:
#                    lines = readfile(atomDir + '/enumpast/training_set_structures.dat')
#                    iidStructs[iatom] = [line.strip().split()[1] for line in lines]           

        if (iteration == 1 and startMethod == 'empty folders') or natoms == 1: #initialize training_set_structures in enumpast/.  Compute iid structures once, and copy to all atom folders that need them
            subprocess.call(['echo','\nChoosing i.i.d. structures for all\n'])                         
            os.chdir('enum')
            if 2 <= nNew[0] < ntot:  #uncle requires at least 2 iid structures.
                subprocess.call([self.uncleExec, '42', str(nNew[0])], stdout=self.uncleOut) 
                lines = readfile('training_set_structures.dat')
            elif nNew[0] == ntot: #asking for all the structures for small enumerations, so just list them
                structlines = ['{}   {}\n'.format(str(i+1),str(i+1)) for i in range(ntot)]
                writefile(structlines,'training_set_structures.dat')
                lines = readfile('training_set_structures.dat')                 
            else: 
                subprocess.call(['echo','\t Number of iid structures requested is less than 2...skipping iids\n']) 
                lines = []  
            for iatom,atom in enumerate(self.atoms):
                atomDir = lastDir + '/' + atom
                iidList = [line.strip().split()[1] for line in lines]                    
                subprocess.call(['echo','\nCopying i.i.d. structures for ' + atom + ' . . .\n'])                         
                epDir = lastDir + '/' + atom + '/enumpast'
                subprocess.call(['cp','training_set_structures.dat',epDir])
                iidStructs[iatom] = iidList 
            os.chdir(lastDir)                                       
        else: # must get separate iid structures for each atom, so parallelize        
            #prep
            for iatom,atom in enumerate(self.atoms):
                if 2 <= nNew[iatom] < ntot:
                    atomDir = lastDir + '/' + atom
    #            try:
                    os.chdir(atomDir + '/enumpast')
                    subprocess.call(['echo','\nChoosing i.i.d. structures for ' + atom + ' . . .\n'])
                    subprocess.call(['ln','-s','../../enum/struct_enum.out'])
                    subprocess.call(['ln','-s','../../enum/lat.in']) 
                    subprocess.call(['ln','-s','../../enum/enum_PI_matrix.out'])
                    subprocess.call(['ln','-s','../../enum/clusters.out'])                          
    #                subprocess.call([self.uncleExec, '42', str(nNew[iatom])], stdout=self.uncleOut)
                    os.chdir(lastDir)
            if distribute and natoms > 1:
                #make job files
                os.chdir(lastDir)
                jobIds = []
                if sum(nNew[1:])>0:
                    mem = '16' #Gb
                    walltime = 8.0 #hrs
                    subdir = 'enumpast'
                    execString = self.uncleExec + ' 42 '
                    atomStrings = [str(n) for n in nNew]
                    parallelJobFiles(self.atoms,subdir,walltime,mem,execString,atomStrings) 
                    #submit jobs for atoms 2 an above
                    jobIds = parallelAtomsSubmit(self.atoms[1:],subdir)
                #use this job to calculate the first atom:
                if 2 <= nNew[0] < ntot:
                    os.chdir(lastDir + '/' + self.atoms[0]  + '/' + subdir)
                    subprocess.call(['echo','\tThis job calculating the first atom: {}. Submitted jobs for the others.\n'.format(self.atoms[0])])
                    subprocess.call([self.uncleExec, '42',str(nNew[0])], stdout=self.uncleOut)             
                    os.chdir(lastDir)      
                #wait for others
                if len(jobIds)>0: parallelAtomsWait(jobIds) 
            else: #run tasks sequentially
                for iatom, atom in enumerate(atoms):
                    if 2 <= nNew[iatom] < ntot:
                        os.chdir(lastDir + '/' + self.atoms[iatom]  + '/' + subdir)
                        subprocess.call(['echo','\tCalculating atom: {}.\n'.format(atom)])
                        subprocess.call([self.uncleExec, '42',str(nNew[iatom])], stdout=self.uncleOut)             
                        os.chdir(lastDir)                        
            #get the iidStructs from training_set_structures.dat for each atom
            for iatom,atom in enumerate(self.atoms):
                if 2 <= nNew[iatom] < ntot:
                    lines = readfile(atomDir + '/enumpast/training_set_structures.dat')
                    iidStructs[iatom] = [line.strip().split()[1] for line in lines]           
#            except:
#                    subprocess.call(['echo','\n~~~~~~~~~~ Could not choose i.i.d. structures for ' + atom + '! ~~~~~~~~~~\n'])
        os.chdir(lastDir)       
        return iidStructs
    
    def enumerate(self):
        """ Runs through the whole process of enumeration, cluster building, and choosing an
            i.i.d. set of training structures. """
        if not os.path.isdir('enum'): subprocess.call(['mkdir','enum'])
        infile = open('needed_files/struct_enum.in','r')
        inlines = []
        for line in infile:
            inlines.append(line)
        infile.close()
        
        structFile = open('enum/struct_enum.in','w')
        npoints = int(inlines[6].split()[0])  
        for i in xrange(len(inlines)):
            if i == 7 + npoints: 
                structFile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + " ")
                structFile.write("# Starting and ending cell sizes for search\n")
            else:
                structFile.write(inlines[i])
        structFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        subprocess.call(['echo','\nEnumerating symmetrically unique structures. . .\n'])
        subprocess.call([self.enumExec,'struct_enum.in'], stdout=self.uncleOut)
               
        self.changeEnumFile() #writes 'bulk' in place of surface
        os.chdir(lastDir)
        subprocess.call(['echo','\nGenerating clusters. . .\n'])
        self.buildClusters()
        os.chdir(lastDir + '/enum')
        #Run the smallest iid job possible to calculate enum_PI_matrix.out.  All the atoms need it.
        subprocess.call(['echo','\nCalculating enum_PI_matrix.out\n']) 
        subprocess.call(['rm', 'enum_PI_matrix.out'])        
        subprocess.call([self.uncleExec, '42', '2'], stdout=self.uncleOut) 
        os.chdir(lastDir)

    def getNtot(self,dir):
        """Gets total number of structures in enumeration"""
        return int(readfile(dir + '/struct_enum.out')[-1].split()[0])
                          
    def makeAtomDirectories(self):
        """ Creates a directory for each atom in the atom list specified in settings.in.  All the 
            VASP and UNCLE files for the atom will be placed in this directory. """
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if not os.path.isdir(atomDir):
                subprocess.call(['mkdir',atomDir])    

                  