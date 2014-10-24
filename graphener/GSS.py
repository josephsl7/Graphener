'''
Created on Aug 29, 2014

@author: eswens13
'''
import os, subprocess

class GSS:
    """ This class performs the UNCLE ground state search for the lowest formation energy
        structures. From this, we get a list of every structure that was enumerated and its
        corresponding formation energy as predicted by UNCLE. We use this list to decide which new
        structures to add into the model. We also keep track of the list and the plots of the list
        from iteration to iteration. """

    def __init__(self, atoms, volRange, plotTitle, xlabel, ylabel, uncleOutput):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.volRange = volRange
        self.plotTitle = plotTitle
        self.xlabel = xlabel
        self.ylabel = ylabel
        
        self.enumFolder = os.getcwd() + '/enum/'
        self.fitsDir = None
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
    
    def makeGSSDirectories(self):
        """ Creates the 'gss' directories for each different metal atom and populates them with 
            the files that UNCLE needs to perform a ground state search.  These files are 
            struct_enum.out, lat.in, J.1.out, groundstatesearch.in, and gss_plot.gp. """
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                self.fitsDir = atomDir + '/fits/'
                gssDir = atomDir + '/gss/'
                subprocess.call(['mkdir',gssDir])
                subprocess.call(['cp',self.enumFolder + 'struct_enum.out', gssDir])
                subprocess.call(['cp',self.enumFolder + 'lat.in', gssDir])
                subprocess.call(['cp',self.fitsDir + 'J.1.out', gssDir])
                
                infile = open(self.neededFilesDir + 'groundstatesearch.in','r')
                inlines = [line for line in infile]
                infile.close()
                
                outfile = open(gssDir + 'groundstatesearch.in','w')
                for i in xrange(len(inlines)):
                    if i == 4:
                        outfile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + "\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
                
                infile = open(self.neededFilesDir + 'gss_plot.gp','r')
                inlines = [line for line in infile]
                infile.close()
                
                outfile = open(gssDir + 'gss_plot.gp','w')
                for i in xrange(len(inlines)):
                    if i == 3:
                        outfile.write("set xlabel \"" + self.xlabel + "\"\n")
                    elif i == 4:
                        outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                    elif i == 5:
                        outfile.write("set title \"" + self.plotTitle + " (" + atom + ")\"\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
    
    def performGroundStateSearch(self, iteration):
        """ Performs the ground state search with the current fit from UNCLE. """
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo','\nPerforming ground state search for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                subprocess.call([self.uncleExec, '21'], stdout=self.uncleOut)
                os.chdir(lastDir)

    def makePlots(self, iteration):
        """ Creates the plots of the predicted energies of all the structures that have been 
            enumerated. Adds the iteration number onto the end of the filenames for the plots and
            the lists. """
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo', '\nMaking plot of ground states for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                subprocess.call(['gnuplot', 'gss_plot.gp'])
                subprocess.call(['mv','gss.out','gss_' + str(iteration) + '.out'])
                subprocess.call(['mv','gss.pdf','gss_' + str(iteration) + '.pdf'])
                os.chdir(lastDir)
    
    def contains(self, struct, alist):
        """ Returns true if 'struct' is found in 'alist', false otherwise. """
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
    
        return False
    
    def getAllGSSStructures(self, iteration, failedStructs):
        """ Returns a list of all the structures sorted by their predicted formation energies.
            It does this for each metal atom that has been specified by the user so this will 
            actually return a list of lists--a list for each atom. """
        allStructs = []
        for n in xrange(len(self.atoms)):
            atomStructs = []
            structsByEnergy = []
            gssFile = os.getcwd() + '/' + self.atoms[n] + '/gss/gss_' + str(iteration) + '.out'
            infile = open(gssFile, 'r')
            
            i = 0
            for line in infile:
                if i >= 2:
                    formEnergy = float(line.strip().split()[7])
                    struct = int(line.strip().split()[0])
                    
                    # Do not include structures that failed VASP calculations.
                    if not self.contains(struct, failedStructs[n]):
                        structsByEnergy.append([formEnergy, struct])
                i += 1
            infile.close()
            
            structsByEnergy.sort()
            
            for struct in structsByEnergy:
                atomStructs.append(str(struct[1]))
            
            allStructs.append(atomStructs)
            
        return allStructs
            
                







        