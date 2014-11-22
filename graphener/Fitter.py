'''
Created on Aug 29, 2014

@author: eswens13
'''
import os
import subprocess


class Fitter:
    """ This class performs the UNCLE fits to the VASP data that has been gathered so far.  It also
        keeps track of the fitting errors, prediction errors, and summaries of cluster expansion 
        terms from iteration to iteration. """

    def __init__(self, atoms, M_fitStructures, N_subsets, vstructsFinished, uncleOutput):
        """ CONSTRUCTOR """  
        self.vstructsFinished = vstructsFinished  
        self.atoms = atoms
        self.M_fitStructures = M_fitStructures
        self.N_subsets = N_subsets
        
        self.enumFolder = os.getcwd() + '/enum/'
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
        self.header = "peratom\nnoweights\nposcar\n"

        
    def makeFitDirectories(self):
        """ Creates the 'fits' directories for each atom and populates the directories with the 
            files that UNCLE needs to perform a fit.  These files are lat.in, CS.in, clusters.out, 
            and the current structures.in and structures.holdout files. """
        for n in xrange(len(self.atoms)):
            atomDir = os.path.abspath(self.atoms[n])
            fitsDir = atomDir + '/fits'
            if not os.path.isdir(fitsDir): subprocess.call(['mkdir',fitsDir])
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
                    if (self.M_fitStructures > len(self.vstructsFinished[n]) and self.M_fitStructures > 0):
                        outfile.write(str(len(self.vstructsFinished[n])) + "\n")
                    else:
                        outfile.write(str(self.M_fitStructures) + "\n")
                elif i == 62:
                    outfile.write(str(self.N_subsets) + "\n")
                else:
                    outfile.write(inlines[i])
            outfile.close()
               
    def fitVASPData(self, iteration):
        """ Performs the UNCLE fit to the VASP data. After the fit is done, it adds the iteration
            onto the end of the files we want to keep track of from iteration to iteration. """
        lastDir = os.getcwd()
        for atom in self.atoms:
            atomDir = lastDir + '/' + atom
            if os.path.isdir(atomDir):
                subprocess.call(['echo','\nFitting VASP data for ' + atom + '. . .\n'])
                fitsDir = atomDir + '/fits'
                if os.path.isdir(fitsDir):
                    os.chdir(fitsDir)
                    subprocess.call([self.uncleExec, '15'], stdout=self.uncleOut)
                    subprocess.call(['mv','fitting_errors.out','fitting_errors_' + str(iteration) + '.out'])
                    subprocess.call(['mv','prediction_errors.out','prediction_errors_' + str(iteration) + '.out'])
                    subprocess.call(['mv','J.1.summary.out','J.1.summary_' + str(iteration) + '.out'])
                    subprocess.call(['cp','structures.in', 'structures.in_' + str(iteration)]) #leave the file to be appended to
                    subprocess.call(['cp','structures.holdout', 'structures.holdout_' + str(iteration)]) #leave the file in case a
                    os.chdir(lastDir)

    def holdoutFromIn(self, atoms,startFromExisting):
        '''Writes structures.holdout for the first iteration, for a given atom
        If starting from existing struct calculations, takes up to N structs in the top of the 
        structures.in file.  Useful when starting from existing structs.  In this case, 
        they should be the lowest N FE structs, since past_structs.dat should be ordered at first.'''            
     
        nmax = 100
        for iatom, atom in enumerate(atoms):
            if startFromExisting[iatom]:
                atomDir = os.getcwd() + '/' + atom
                
                infile = open(atomDir + '/fits/structures.in', 'r')
                holdoutFile = open(atomDir + '/fits/structures.holdout', 'w')
#                holdoutFile.write(self.header) #bch why are 2 headers being written if I uncomment this?
                count = 0
                for line in infile:
                    if list(line.strip().split()[0])[:2] == ['#', '-']:
                        if count >= nmax:
                            break
                        count += 1    
                    holdoutFile.write(line)
    
            infile.close()
            holdoutFile.close()






