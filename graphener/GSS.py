'''
Created on Aug 29, 2014

@author: eswens13
'''
import os, subprocess

class GSS:
    

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
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo','\nPerforming ground state search for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                subprocess.call([self.uncleExec, '21'], stdout=self.uncleOut)
                os.chdir(lastDir)

    def makePlots(self, iteration):
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
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
    
        return False
    
    def getAllGSSStructures(self, iteration, failedStructs):
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
            
                







        