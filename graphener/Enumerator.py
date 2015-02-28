'''
'''
import os, subprocess, shutil
import ClustersBuild
from comMethods import *
from numpy import *

class Enumerator:
    """ This class enumerates symmetrically unique structures in a given volume range using UNCLE.  
        It then builds the clusters necessary to perform a cluster expansion and chooses a 
        specified number of "training structures" to perform a first fit on.  After this class 
        finishes its work, the Extractor class will take over and extract pseudo-POSCAR files from 
        the struct_enum.out file that is produced. The methods in this class are only needed for 
        the first iteration of the main convergence loop. """
  
    def __init__(self, atoms, volRange, clusterNums, uncleOutput, distribute, enumVcNum):
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
	self.case = len(atoms[0].split(','))
        self.distribute = distribute
        self.enumVcNum = enumVcNum

    def buildClusters(self):
        """ Uses UNCLE to build the number of each n-body clusters specified in the settings.in
            file. """
        oldLatFile = 'needed_files/lat.in'
        oldFile = open(oldLatFile, 'r')
        oldLines = [line for line in oldFile]
        oldFile.close()
        
        newFile = open('enum/lat.in','w')
        for i in xrange(len(oldLines)):
            if i>=1 and 'Number pairs' in oldLines[i-1]: #bch use label on previous line
                for num in self.clusterNums:
                    newFile.write(str(num) + " ")
                newFile.write("\n")
            elif '#case' in oldLines[i]:
                newFile.write(str(self.case) + "  #case\n")
            else:
                newFile.write(oldLines[i])
        newFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
#        subprocess.call([self.uncleExec, '10'], stdout=self.uncleOut)
        if sum(self.clusterNums)<=1500: #the 1500 rule of thumb assumes you are running Main with 16G. 
            subprocess.call([self.uncleExec, '10'], stdout=self.uncleOut)
        else:
            subprocess.call(['echo','Warning: BLOCKING CLUSTER JOB to save time'])
#            clustersjob = ClustersBuild.clustersjob()
#            clustersjob.clustBuild()
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
            if oldlines[i].strip().split()[0] == 'start':
                startline = i + 1
                break

        vec1 = oldlines[2].strip().split()[:3]
        vec2 = oldlines[3].strip().split()[:3]
        vec3 = oldlines[4].strip().split()[:3]
       
        vec1comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec1] #to handle ******
        vec2comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec2]
        vec3comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec3]
    
        #Parent lattice vectors
        pLV = array([[vec1comps[0],vec2comps[0],vec3comps[0]], [vec1comps[1],vec2comps[1],vec3comps[1]], [vec1comps[2],vec2comps[2],vec3comps[2]]])
        nD = int(oldlines[5].strip().split()[0]) #Number of sites per cell
        pBas = array([[float(comp) for comp in line.strip().split()[:3]] for line in oldlines[6:6 + nD]]) #Base vectors for sites in the primitive cell

        structNum = 1
        newHNF = True
        for i in xrange(startline, len(oldlines)):
            if i != startline:
                if oldlines[i-1].strip().split()[8:26] != oldlines[i].strip().split()[8:26]:
                    newHNF = True
            if newHNF == True:
                S = array([int(x) for x in oldlines[i].strip().split()[8:11]]) #SNF vector
                [a,b,c,d,e,f] = [int(x) for x in oldlines[i].strip().split()[11:17]] #From HNF matrix
                L = array([[int(x) for x in oldlines[i].strip().split()[17:20]], [int(x) for x in oldlines[i].strip().split()[20:23]], [int(x) for x in oldlines[i].strip().split()[23:26]]]) #Left transform
                aBas = array(zeros((3, nD*a*c*f))) #The base vectors for sites in the structure
                gIndx = array(zeros(nD*a*c*f)) #Mapping function showing the site that goes with each structure in the labeling

                ic = -1  #counter
                for iD in range(nD): # Loop over the number at sites/parent cell (the d set)
                    for z1 in range(a):
                        for z2 in range((b*z1)/a, c+(b*z1)/a):
                            for z3 in range(z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c):
                                ic = ic + 1
                                aBas[:,ic]=dot(pLV,array([z1,z2,z3]))+pBas[iD,:] # Atomic basis vector in Cartesian coordinates
                                greal = dot(L,array([z1,z2,z3])) # Map position into the group.
                                g = greal.astype(int) # Convert the g-vector from real to integer
                                g = fmod(g,S) # Bring the g-vector back into the first tile.
                                gIndx[ic] = iD*S[0]*S[1]*S[2]+g[0]*S[1]*S[2]+g[1]*S[2]+g[2]+1

                forbid = []
                for vecNum1 in range(nD*a*c*f - 1):
                    for vecNum2 in range(vecNum1 + 1, nD*a*c*f):
                        if linalg.norm(aBas[:,vecNum2]-aBas[:,vecNum1]) < 1.9:
                            forbid.append([int(gIndx[vecNum1]) - 1, int(gIndx[vecNum2]) -1])
                newHNF = False

            oldline = oldlines[i]
            label = oldline.strip().split()[-1]
            invalidPairs = [label[pair[0]] != str(self.enumVcNum) and label[pair[1]] != str(self.enumVcNum) for pair in forbid]
            if True not in invalidPairs or self.enumVcNum == -1:
                newline = '%11i%10i%8i%9i%9i%12i%4i%6i%4i%3i%3i%5i%3i%3i%3i%3i%3i%7i%5i%5i%5i%5i%5i%5i%5i%5i' % tuple([structNum] + [int(x) for x in oldline.strip().split()[1:-1]]) + '   ' + oldline.strip().split()[-1] + '\n'
                newfile.write(newline)
                structNum += 1
        
        newfile.close()
        
    def chooseTrainingStructures(self,iteration, startMethod,nNew,nTotClusters):
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
#                if 2 <= nNew[iatom] < nTotClusters:
#                    lines = readfile(atomDir + '/enumpast/training_set_structures.dat')
#                    iidStructs[iatom] = [line.strip().split()[1] for line in lines]           

        if (iteration == 1 and startMethod == 'empty folders') or natoms == 1: #initialize training_set_structures in enumpast/.  Compute iid structures once, and copy to all atom folders that need them
            subprocess.call(['echo','\nChoosing i.i.d. structures for all\n'])                         
            os.chdir('enum')
            if 2 <= nNew[0] < nTotClusters:  #uncle requires at least 2 iid structures.
                subprocess.call([self.uncleExec, '42', str(nNew[0])], stdout=self.uncleOut) 
                lines = readfile('training_set_structures.dat')
            elif nNew[0] == nTotClusters: #asking for all the structures for small enumerations, so just list them
                structlines = ['{}   {}\n'.format(str(i+1),str(i+1)) for i in range(nTotClusters)]
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
                if 2 <= nNew[iatom] < nTotClusters:
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
            if self.distribute and natoms > 1:
                #make job files
                os.chdir(lastDir)
                jobIds = []
                if sum(nNew[1:])>0:
                    mem = '16' #Gb
                    walltime = 16.0 #hrs
                    subdir = 'enumpast'
                    execString = self.uncleExec + ' 42 '
                    atomStrings = [str(n) for n in nNew]
                    parallelJobFiles(self.atoms,subdir,walltime,mem,execString,atomStrings) 
                    #submit jobs for atoms 2 an above
                    jobIds = parallelAtomsSubmit(self.atoms[1:],subdir)
                #use this job to calculate the first atom:
                if 2 <= nNew[0] < nTotClusters:
                    os.chdir(lastDir + '/' + self.atoms[0]  + '/' + subdir)
                    subprocess.call(['echo','\tThis job calculating the first atom: {}. Submitted jobs for the others.\n'.format(self.atoms[0])])
                    subprocess.call([self.uncleExec, '42',str(nNew[0])], stdout=self.uncleOut)             
                    os.chdir(lastDir)      
                #wait for others
                if len(jobIds)>0: parallelAtomsWait(jobIds) 
            else: #run tasks sequentially
                for iatom, atom in enumerate(atoms):
                    if 2 <= nNew[iatom] < nTotClusters:
                        os.chdir(lastDir + '/' + self.atoms[iatom]  + '/' + subdir)
                        subprocess.call(['echo','\tCalculating atom: {}.\n'.format(atom)])
                        subprocess.call([self.uncleExec, '42',str(nNew[iatom])], stdout=self.uncleOut)             
                        os.chdir(lastDir)                        
            #get the iidStructs from training_set_structures.dat for each atom
            for iatom,atom in enumerate(self.atoms):
                if 2 <= nNew[iatom] < nTotClusters:
                    lines = readfile(atomDir + '/enumpast/training_set_structures.dat')
                    iidStructs[iatom] = [line.strip().split()[1] for line in lines]           
#            except:
#                    subprocess.call(['echo','\n~~~~~~~~~~ Could not choose i.i.d. structures for ' + atom + '! ~~~~~~~~~~\n'])
        os.chdir(lastDir)       
        return iidStructs
    
    def enumerate(self):
        """ Runs through the whole process of enumeration, cluster building, and choosing an
            i.i.d. set of training structures. """
        if os.path.isdir('enum'):
            shutil.rmtree('enum')

        subprocess.call(['mkdir','enum'])

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
            elif i == 5:
                structFile.write(" " + str(self.case) + " -nary case\n")
            elif i in range(7, 8 + npoints):
                structFile.write(inlines[i].split()[0]+" ")
                structFile.write(inlines[i].split()[1]+" ")
                structFile.write(inlines[i].split()[2]+"    ")
                for num in range(self.case):
                    structFile.write(str(num))
                    if num < self.case-1:
                        structFile.write("/")
                structFile.write("   # d0" + str(i-6) + " d-vector\n")
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
       
        #Run the smallest iid job possible to calculate enum_PI_matrix.out.  All the atoms need it.
        subprocess.call(['echo','\nCalculating enum_PI_matrix.out\n']) 
        subprocess.call(['rm', 'enum_PI_matrix.out'])  
        if self.distribute: #not really distribution over atoms, but using this method to submit a job. still may need larger memory for this  
                mem = '64' #Gb
                walltime = 16.0 #hrs
                subdir = '../enum/'
                execString = self.uncleExec + ' 42 '
                atomStrings = ['2']
                parallelJobFiles([self.atoms[0]],subdir,walltime,mem,execString,atomStrings)             
                jobIds = parallelAtomsSubmit([self.atoms[0]],subdir)  
                subprocess.call(['echo','\tSubmitted job {}\n'.format(jobIds[0])]) 
                if len(jobIds)>0: parallelAtomsWait(jobIds)
        else:
            os.chdir(lastDir + '/enum')
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

                  
