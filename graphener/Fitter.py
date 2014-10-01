'''
Created on Aug 29, 2014

@author: eswens13
'''
import os
import subprocess


class Fitter:


    def __init__(self, atoms, M_fitStructures, N_subsets):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.M_fitStructures = M_fitStructures
        self.N_subsets = N_subsets
        self.enumFolder = os.getcwd() + '/enum/'
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        
    def makeFitDirectories(self):
        for atom in self.atoms:
            atomDir = os.path.abspath(atom)
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
                    outfile.write(str(self.M_fitStructures) + "\n")
                elif i == 62:
                    outfile.write(str(self.N_subsets) + "\n")
                else:
                    outfile.write(inlines[i])
            outfile.close()
               
    def fitVASPData(self, iteration):
        lastDir = os.getcwd()
        for atom in self.atoms:
            atomDir = lastDir + '/' + atom
            if os.path.isdir(atomDir):
                subprocess.call(['echo','\nFitting VASP data for ' + atom + '. . .\n'])
                fitsDir = atomDir + '/fits'
                if os.path.isdir(fitsDir):
                    os.chdir(fitsDir)
                    subprocess.call([self.uncleExec, '15'])
                    subprocess.call(['mv','fitting_errors.out','fitting_errors_' + str(iteration) + '.out'])
                    subprocess.call(['mv','J.1.summary.out','J.1.summary_' + str(iteration) + '.out'])
                    os.chdir(lastDir)








