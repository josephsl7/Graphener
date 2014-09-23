'''
Created on Aug 29, 2014

@author: eswens13
'''
import os, subprocess

class GSS:
    

    def __init__(self, atoms, volRange, plotTitle, xlabel, ylabel):
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
    
    def performGroundStateSearch(self):
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo','\nPerforming ground state search for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                subprocess.call([self.uncleExec, '21'])
                os.chdir(lastDir)

    def makePlots(self):
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo', '\nMaking plot of ground states for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                subprocess.call(['gnuplot', 'gss_plot.gp'])
                os.chdir(lastDir)
    
    def getAllGSSStructures(self):
        allStructs = []
        for atom in self.atoms:
            atomStructs = []
            structsByEnergy = []
            gssFile = os.getcwd() + '/' + atom + '/gss/gss.out'
            infile = open(gssFile, 'r')
            
            i = 0
            for line in infile:
                if i >= 2:
                    structsByEnergy.append([float(line.strip().split()[7]), int(line.strip().split()[0])])
            infile.close()
            
            structsByEnergy.sort()
            for struct in structsByEnergy:
                atomStructs.append(struct[1])
            
            allStructs.append(atomStructs)
        
        return allStructs
            
                







        