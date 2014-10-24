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

    def __init__(self, atoms, M_fitStructures, N_subsets, structuresInLengths, uncleOutput):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.M_fitStructures = M_fitStructures
        self.N_subsets = N_subsets
        self.structuresInLengths = structuresInLengths
        
        self.enumFolder = os.getcwd() + '/enum/'
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
        
    def makeFitDirectories(self):
        """ Creates the 'fits' directories for each atom and populates the directories with the 
            files that UNCLE needs to perform a fit.  These files are lat.in, CS.in, clusters.out, 
            and the current structures.in and structures.holdout files. """
        for n in xrange(len(self.atoms)):
            atomDir = os.path.abspath(self.atoms[n])
            fitsDir = atomDir + '/fits/'
            subprocess.call(['mkdir',fitsDir])
            subprocess.call(['cp',self.enumFolder + 'lat.in', fitsDir])
            subprocess.call(['cp',self.enumFolder + 'clusters.out', fitsDir])
            subprocess.call(['cp',atomDir + '/structures.in', fitsDir])
            subprocess.call(['cp',atomDir + '/structures.holdout', fitsDir])
            
            infile = open(self.neededFilesDir + 'CS.in','r')
            inlines = [line for line in infile]
            infile.close()
            
            outfile = open(fitsDir + 'CS.in','w')
            for i in xrange(len(inlines)):
                if i == 60:
                    if (self.M_fitStructures > self.structuresInLengths[n] and self.M_fitStructures > 0):
                        outfile.write(str(self.structuresInLengths[n]) + "\n")
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
                    os.chdir(lastDir)








