'''
Created on Aug 26, 2014

@author: eswens13
'''
import os, subprocess


class Enumerator:
    """ This class enumerates symmetrically unique structures in a given volume range using UNCLE.  
        It then builds the clusters necessary to perform a cluster expansion and chooses a 
        specified number of "training structures" to perform a first fit on.  After this class 
        finishes its work, the Extractor class will take over and extract pseudo-POSCAR files from 
        the struct_enum.out file that is produced. The methods in this class are only needed for 
        the first iteration of the main convergence loop. """
  
    def __init__(self, atoms, volRange, clusterNums, trainStructNum, uncleOutput):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.volRange = volRange
        
        self.clusterNums = clusterNums
        self.trainStructNum = trainStructNum
        
        self.uncleExec = os.path.abspath('needed_files/uncle.x')
        self.enumFile = os.path.abspath('enum/struct_enum.out')
        self.enumExec = os.path.abspath('needed_files/enum.x')
        self.uncleOut = uncleOutput

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
        """ Chooses a list of i.i.d. structures from struct_enum.out for each different metal atom. 
            The UNCLE option that we run to choose the training structures should look for a file 
            called 'past_structs.dat' so as not to choose structures from that list. (Dr. Hess 
            added that option in the UNCLE executable we're using.) The length of the list of
            training structures for each atom is determined by the TRAINING_STRUCTS setting in 
            settings.in. """
        lastDir = os.getcwd()
        
        for atom in self.atoms:
            neededFilesDir = lastDir + '/needed_files'
            atomDir = lastDir + '/' + atom
            try:
                # Look for the past_structs.dat file in needed_files folder.  If it is there, copy
                # it to the atom's enum/ directory. If there's not, make an empty one for that atom.
                pastStructFile = neededFilesDir + '/past_structs.' + atom + '.dat'
                if os.path.exists(pastStructFile):
                    subprocess.call(['cp', pastStructFile, atomDir + '/enum/past_structs.dat'])
                else:
                    emptyFile = open(atomDir + '/enum/past_structs.dat','w')
                    emptyFile.close()
                    
                os.chdir(atomDir + '/enum')
                subprocess.call(['echo','\nChoosing i.i.d. structures for ' + atom + ' . . .\n'])
                subprocess.call([self.uncleExec, '42', str(self.trainStructNum)], stdout=self.uncleOut)
                os.chdir(lastDir)
            except:
                subprocess.call(['echo','\n~~~~~~~~~~ Could not choose i.i.d. structures for ' + atom + '! ~~~~~~~~~~\n'])

    def makeAtomDirectories(self):
        """ Creates a directory for each atom in the atom list specified in settings.in.  All the 
            VASP and UNCLE files for the atom will be placed in this directory. """
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if not os.path.isdir(atomDir):
                subprocess.call(['mkdir',atomDir])
    
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
        subprocess.call(['echo','\nEnumerating symmetrically unique structures. . .\n'])
        subprocess.call([self.enumExec,'struct_enum.in'], stdout=self.uncleOut)
        
        os.chdir(lastDir)
        
        self.changeEnumFile()
        
        subprocess.call(['echo','\nGenerating clusters. . .\n'])
        self.buildClusters()
        
        self.makeAtomDirectories()
        for atom in self.atoms:
            subprocess.call(['cp','-r','enum', atom + '/'])
            if os.path.exists('needed_files/structures.start.' + atom):
                subprocess.call(['cp','needed_files/structures.start.' + atom, atom + '/structures.in.base'])
        
        self.chooseTrainingStructures()
        
            
            
        