'''
Created on Aug 29, 2014


'''
import os,sys,subprocess,time,shutil
from comMethods import *
from numpy import ceil

class Fitter:
    """ This class performs the UNCLE fits to the VASP data that has been gathered so far.  It also
        keeps track of the fitting errors, prediction errors, and summaries of cluster expansion 
        terms from iteration to iteration. """

    def __init__(self, atoms, M_fitStructures, N_subsets, vstructsFinished, uncleOutput,distribute):
        """ CONSTRUCTOR """  
        self.vstructsFinished = vstructsFinished  
        self.atoms = atoms
        self.M_fitStructures = M_fitStructures
        self.N_subsets = N_subsets
        self.enumFolder = os.getcwd() + '/enum/'
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
        #self.header = "peratom\nnoweights\nposcar\n"
        self.vstructsFinished = vstructsFinished
        self.distribute = distribute

    def filterStructuresIn(self, fitsDir, iteration, maxE):
        '''Remove structs from structures.in, just before fitting, that don't fit criteria here. 
           For now, structures above maxE will be removed'''
        if maxE < 100: # < the absurd default
            inFile = fitsDir + '/structures.in'
            lines = readfile(inFile)
            subprocess.call(['mv',inFile,inFile + '_{}full'.format(iteration)]) 
            outFile = open(fitsDir + '/structures.in','w')
            outFile.write("peratom\nnoweights\nposcar\n")
            for j in xrange(len(lines)):
                if list(lines[j].strip().split()[0])[:2] == ['#','-']:
                    if j != len(lines) - 1:
                        FE = float(lines[j+1].strip().split()[6].strip(','))
                        if FE <= maxE:
                            natomsList = [int(i) for i in lines[j+6].strip().split()]
                            natoms = natomsList[0]+natomsList[1]
                            outFile.writelines(lines[j:j+10+natoms])
                        else:
                            struct = lines[j + 1].strip().split()[3]
                            subprocess.call(['echo','\tNot including struct {} in fit (exceeds maxE)'.format(struct)])
            outFile.close() 
            
    def filterStructuresInFrac(self, fitsDir, iteration, cullFrac):
        '''Remove structs from structures.in, just before fitting, that don't fit criteria here. 
           For now, the top cullFrac of uncle formation energies will be removed'''
        if cullFrac > 0: 
            # Get the structs and FE ordered by energy 
            FEdata = zeros(20000,dtype = [('struct', int32),('FE', float),('keep',bool)])
            inFile = fitsDir + '/structures.in'
            lines = readfile(inFile)
            subprocess.call(['mv',inFile,inFile + '_{}full'.format(iteration)])
            nstructs = 0
            for j in xrange(len(lines)):
                if list(lines[j].strip().split()[0])[:2] == ['#','-']:
                    if j != len(lines) - 1:
                        struct = lines[j + 1].strip().split()[3]
                        FE = float(lines[j+1].strip().split()[6].strip(','))
                        nstructs += 1
                        FEdata[nstructs]['FE'] = FE
                        FEdata[nstructs]['struct'] = struct
            sort(FEdata,order = ['FE','struct']) #from low to high
            #keep only part of them
            nkeep = ceil((1-cullFrac)*nstructs)
            FEdata[0:nkeep]['keep'] = [True]*nkeep
            FEdata[nkeep:nstructs-nkeep+1]['keep'] = [False]*(nstructs-nkeep)
            #scan lines again and write only the ones to keep
            outFile = open(fitsDir + '/structures.in','w')
            outFile.write("peratom\nnoweights\nposcar\n")
            for j in xrange(len(lines)):
                if list(lines[j].strip().split()[0])[:2] == ['#','-']:
                    if j != len(lines) - 1:
                        struct = int(lines[j + 1].strip().split()[3])
                        istruct = next(i for i in FEdata[:nstructs+1]['struct'] if i == struct)
                        if FEdata[istruct]['keep']:
                            natomsList = [int(i) for i in lines[j+6].strip().split()]
                            natoms = natomsList[0]+natomsList[1]
                            outFile.writelines(lines[j:j+10+natoms])
                        else:                       
                            subprocess.call(['echo','\tNot including struct {} in fit (removing top {} in FE))'.format(struct,cullFrac)])
            outFile.close()
                
    def fitVASPData(self, iteration, maxE):
        """ Performs the UNCLE fit to the VASP data. """
        subprocess.call(['echo','\nFitting VASP data . . .\n'])
        natoms = len(self.atoms)
        lastDir = os.getcwd()
        subdir = 'fits'
        #prep for all
        for iatom, atom in enumerate(self.atoms):
            nfinished = len(self.vstructsFinished[iatom])
            if nfinished > 1: #don't try fitting if structures.in is too small
                atomDir = lastDir + '/' + atom
                if nfinished < 100:
                    cullFrac = 0
                else:
                    cullFrac = 0.01
                if os.path.isdir(atomDir):
                    fitsDir = atomDir + '/fits'
                    if os.path.isdir(fitsDir):
                        os.chdir(fitsDir)
                        subprocess.call(['cp', atomDir + '/structures.in', '.' ]) #so we have the latest version here 
#                        self.filterStructuresInFrac(fitsDir,iteration, cullFrac) #remove some structures at the top of the FE list.                   
#                            check = subprocess.check_output([self.uncleExec, '15'])
#                            subprocess.call(['echo','Uncle 15 feedback'+ check])s
        if self.distribute and natoms > 1: #parallelize the atom jobs
            os.chdir(lastDir)
            mem = '16' #Gb
            walltime = 2.0 #hrs
            execString = self.uncleExec + ' 15'
            atomStrings = ['']*natoms
            parallelJobFiles(self.atoms,subdir,walltime,mem,execString,atomStrings) 
            #submit jobs
            jobIds = parallelAtomsSubmit(self.atoms[1:],subdir)
            #use this job to calculate the first atom:
            os.chdir(lastDir + '/' + self.atoms[0]  + '/' + subdir)
            subprocess.call(['echo','\n\tThis job calculating the first atom: {}. Submitted jobs for the others.\n'.format(self.atoms[0])])
            subprocess.call([self.uncleExec, '15'], stdout=self.uncleOut)             
            os.chdir(lastDir)
            #wait
            parallelAtomsWait(jobIds)
        else: #run tasks sequentially
            for iatom, atom in enumerate(self.atoms):
                os.chdir(lastDir + '/' + self.atoms[iatom]  + '/' + subdir)
                subprocess.call(['echo','\tCalculating atom: {}.\n'.format(atom)])
                subprocess.call([self.uncleExec, '15'], stdout=self.uncleOut)             
                os.chdir(lastDir)                  
                                  
        #post calc work for all
        for iatom, atom in enumerate(self.atoms):
            if len(self.vstructsFinished[iatom]) > 1: #don't try fitting if structures.in is too small
                atomDir = lastDir + '/' + atom
                if os.path.isdir(atomDir):
                    fitsDir = atomDir + '/fits'
                    if os.path.isdir(fitsDir):
                        os.chdir(fitsDir) 
                        subprocess.call(['mv','fitting_errors.out','fitting_errors_' + str(iteration) + '.out'])
                        subprocess.call(['mv','prediction_errors.out','prediction_errors_' + str(iteration) + '.out'])
                        subprocess.call(['mv','J.1.summary.out','J.1.summary_' + str(iteration) + '.out'])
                        subprocess.call(['cp','structures.in', 'structures.in_' + str(iteration)]) #also leaves a copy of the file to be appended to
                        subprocess.call(['cp','structures.holdout', 'structures.holdout_' + str(iteration)]) #leave the file in case a
        os.chdir(lastDir)
    
    def makeFitDirectories(self):
        """ Creates the 'fits' directories for each atom and populates the directories with the 
            files that UNCLE needs to perform a fit.  These files are lat.in, CS.in, clusters.out, 
            and the current structures.in and structures.holdout files. """
        for iatom, atom in enumerate(self.atoms):
            atomDir = os.path.abspath(atom)
            fitsDir = atomDir + '/fits'
            if os.path.isdir(fitsDir): #remove it...start clean because must have current files

                try:
                    check = subprocess.check_output(['rm','-r',fitsDir])
                except:
                    subprocess.call(['echo','ERROR in removing /fits for atom {}'.format(atom)])

            subprocess.call(['mkdir',fitsDir])
            subprocess.call(['ln','-s',self.enumFolder + '/struct_enum.out',fitsDir])
            subprocess.call(['ln','-s',self.enumFolder + '/lat.in',fitsDir])
            subprocess.call(['ln','-s',self.enumFolder + '/clusters.out',fitsDir])

            infile = open(self.neededFilesDir + '/CS.in','r')
            inlines = [line for line in infile]
            infile.close()
            # TODO:  This doesn't work right now unless it's a negative number in settings.in
            outfile = open(fitsDir + '/CS.in','w')
            for i in xrange(len(inlines)):
                if i == 60:
                    if (self.M_fitStructures > len(self.vstructsFinished[iatom]) and self.M_fitStructures > 0):
                        outfile.write(str(len(self.vstructsFinished[iatom])) + "\n")
                    else:
                        outfile.write(str(self.M_fitStructures) + "\n")
                elif i == 62:
                    outfile.write(str(self.N_subsets) + "\n")
                else:
                    outfile.write(inlines[i])
            outfile.close()

    def writeHoldout(self, N, structs,vdata):
        '''Writes structures.holdout from a list of struct names for each atom.'''            
        for iatom in xrange(len(self.atoms)):
            nmax = min(N,len(structs[iatom]))
            atomDir = os.path.abspath(self.atoms[iatom])
            if os.path.exists('{}/structures.holdout'.format(atomDir)):
                subprocess.call(['rm', atomDir + '/structures.holdout'])
            structuresWrite(nmax,atomDir,self.vstructsFinished[iatom],\
                            vdata[iatom,:nmax]['FE'],vdata[iatom,:nmax]['conc'],\
                            vdata[iatom,:nmax]['energy'],'.holdout','w')        
            subprocess.call(['cp', atomDir + '/structures.holdout', atomDir + '/fits/structures.holdout'])
        
