'''
Created on Aug 29, 2014

@author: eswens13
'''
import os, subprocess
from random import random, seed

class MakeUncleFiles:


    def __init__(self, atoms):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
    
        self.structList = []
        self.pureHenergy = 0.0
        self.pureMenergy = 0.0
        
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0
    
    def initialize(self):
        self.structList = []
        
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0      
    
    def FinishCheck(self, folder):
        """ Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR.  The input
            parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """
            
        lastfolder = os.getcwd()
        os.chdir(folder)
        
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)   
         
        return newstring[0].find('Voluntary') > -1 #True/False
    
    def setStructureList(self, atomDir):
        """ Initializes the list of structures to add to the structures.in and structures.holdout
            files by adding only the structures that VASP finished relaxing to the member
            structList. """
        
        print "In setStructureList()"
        
        self.structList = []
        
        pureHdir = os.getcwd() + '/1'
        pureMdir = os.getcwd() + '/3'
        
        self.setAtomCounts(pureHdir)
        self.setEnergy(pureHdir)
        self.pureHenergy = float(self.energy)
        
        self.setAtomCounts(pureMdir)
        self.setEnergy(pureMdir)
        self.pureMenergy = float(self.energy)
        
        print "Got pure energies"
        print "Pure H = " + str(self.pureHenergy)
        print "Pure M = " + str(self.pureMenergy)
        
        conclist = []
        energyList = []
        dirList = os.listdir(os.getcwd())
        for item in dirList:
            fullPath = os.path.abspath(item)
            if os.path.isdir(fullPath):
                if os.path.isdir(fullPath + '/DOS'):
                    if self.FinishCheck(fullPath + '/DOS'):
                        
                        # Check for concentration
                        self.setAtomCounts(fullPath)
                        self.setEnergy(os.getcwd() + '/' + item)
                        structEnergy = float(self.energy)
                        
                        concentration = 0.0
                        if self.atomCounts[0] == 0:
                            concentration = 1.0
                        else:
                            concentration = float(float(self.atomCounts[1]) / float(self.atomCounts[0] + self.atomCounts[1]))
                        
                        formationEnergy = structEnergy - (concentration * self.pureMenergy + (1.0 - concentration) * self.pureHenergy)
                        
                        conclist.append([concentration, fullPath])
                        energyList.append([formationEnergy, fullPath])
        
        conclist.sort()
        energyList.sort()
        
        for i in xrange(len(conclist)):
            self.structList.append(energyList[i][1])    
    
    def getStructureList(self):
        """ Returns the list of usable structures. """
        return self.structList
    
    def setLatticeVectors(self, structFile):
        """ Gets the lattice vectors from the first structure in the structList and sets
            the corresponding member components. """
        vecFile = open(structFile + '/POSCAR','r')
        vecFileLines = vecFile.readlines()
        vecFile.close()
        
        vec1 = vecFileLines[2].strip().split()
        vec2 = vecFileLines[3].strip().split()
        vec3 = vecFileLines[4].strip().split()
        
        vec1comps = [float(comp) for comp in vec1]
        vec2comps = [float(comp) for comp in vec2]
        vec3comps = [float(comp) for comp in vec3]
            
        
        self.vec1x = vec1comps[0]
        self.vec1y = vec1comps[1]
        if vec1comps[2] == 15.0:
            self.vec1z = 1000.0
        else:
            self.vec1z = vec1comps[2]
        
        self.vec2x = vec2comps[0]
        self.vec2y = vec2comps[1]
        if vec2comps[2] == 15.0:
            self.vec2z = 1000.0
        else:
            self.vec2z = vec2comps[2]
        
        self.vec3x = vec3comps[0]
        self.vec3y = vec3comps[1]
        if vec3comps[2] == 15.0:
            self.vec3z = 1000.0
        else:
            self.vec3z = vec3comps[2]
    
    def closeOutFiles(self):
        """ Closes both the structures.in and structures.holdout files. """
        self.infile.close()
        self.holdoutFile.close()

    def setIDString(self, poscarDir):
        """ Sets the first written line of each structure to the form:
                C H (Metal Atom)  Structure:  #(Decimal) (#(Binary)) """
        poscar = open(poscarDir + '/POSCAR', 'r')
        ID = poscar.readlines()[0].strip()
        
        self.idString = ID

    def setAtomPositions(self, poscarDir):
        """ Retrieves the positions of each of the atoms.  Appends the x-coordinate to the xPositions
            list, the y-coordinate to the yPositions list.  For a surface in UNCLE the z-coordinate
            is always zero. """
        poscar = open(poscarDir + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomPositions = poscarLines[7 + sum(self.atomCounts):7 + (2 * sum(self.atomCounts))]
        self.atomPositions = [line.strip().split() for line in self.atomPositions]
           
        for pos in self.atomPositions:
            self.xPositions.append(float(pos[0]))
            self.yPositions.append(float(pos[1]))
            self.zPositions.append(0.0)

    def setAtomCounts(self, poscarDir):
        """ Retrieves the number of H atoms and the number of M atoms from the POSCAR file
            and sets the corresponding members. """
        self.atomCounts = []

        poscar = open(poscarDir + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        counts = poscarLines[5].strip().split()
        
        if len(counts) == 3:
            self.atomCounts.append(int(counts[1]))
            self.atomCounts.append(int(counts[2]))
        elif len(counts) == 2:
            if poscarLines[0].split()[1] == 'H':
                self.atomCounts.append(int(counts[1]))
                self.atomCounts.append(0)
            elif poscarLines[0].split()[1] == 'M':
                self.atomCounts.append(0)
                self.atomCounts.append(int(counts[1]))

    def setEnergy(self, directory):  
        """ Retrieves the energy of the structure from the OSZICAR file and sets the corresponding 
            member. """  
        print "Directory = " + directory       
        try:
            oszicar = open(directory + '/DOS/OSZICAR','r')
            energy = oszicar.readlines()[-1].split()[2]
            oszicar.close()
        except:
            energy = 0
        
        energy = float(energy)
        peratom = energy / sum(self.atomCounts)
        
        self.energy = str(peratom)

    def writeHeader(self):
        """ Writes the headers of the structures.in and structures.holdout files. """
        self.infile.write(self.header)
        self.holdoutFile.write(self.header)
    
    def writeDashedLine(self):
        """ Writes a dashed line in the structures.in/.holdout files as a separator between
            different structures. """
        self.outfile.write("#------------------------------------------------\n")
    
    def writeIDString(self):
        """ Writes the ID string of the current structure to either the structures.in or 
            structures.holdout file. """
        concentration = float(float(self.atomCounts[1]) / float(sum(self.atomCounts)))
        formationEnergy = float(self.energy) - (concentration * self.pureMenergy + (1.0 - concentration) * self.pureHenergy)
        
        self.outfile.write(self.idString + " FE = " + str(formationEnergy) + ", Concentration = " + str(concentration) + "\n")
        
    def writeLatticeVecs(self):
        """ Writes the lattice vectors of the current structure to the structures.in or
            structures.holdout file. """
        self.outfile.write("1.0\n")
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec1z, self.vec1x, self.vec1y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec2z, self.vec2x, self.vec2y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec3z, self.vec3x, self.vec3y))
    
    def writeAtomCounts(self):
        """ Writes the number of H atoms and the number of M atoms to the structures.in or 
            structures.holdout file. """
        if len(self.atomCounts) == 2:
            self.outfile.write(str(self.atomCounts[0]) + " " + str(self.atomCounts[1]) + "\n")
        elif len(self.atomCounts) == 1:
            self.outfile.write(str(self.atomCounts[0]) + "\n")
    
    def writeAtomPositions(self):
        """ Writes the positions of the atoms in the current structure to the structures.in
            or structures.holdout file.  The positions are written in the form:
                z-coord   x-coord   y-coord
            because this is the convention that UNCLE uses. """
        self.outfile.write("Cartesian\n")
        for i in xrange(len(self.atomPositions)):
            self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % 
                               (self.zPositions[i], self.xPositions[i], self.yPositions[i]))
     
    def writeEnergy(self):
        """ Writes the energy of the current structure to the structures.in or structures.holdout
            file. """
        self.outfile.write("#Energy:\n")
        self.outfile.write(str(self.energy) + "\n") 
    
    def writePOSCAR(self, poscarDir):
        """ Calls all the methods needed to write all the needed information about the current
            structure to the structures.in or structures.holdout files.  Puts a maximum of 10%
            of the structures in the structures.holdout file. """
        if self.holdoutCount / float(len(self.structList)) < .10 and random() < .15:
            self.outfile = self.holdoutFile
            self.holdoutCount += 1
        else:
            self.outfile = self.infile
            self.inCount += 1
        
        self.setIDString(poscarDir)
        self.setLatticeVectors(poscarDir)
        self.setAtomCounts(poscarDir)
        self.setAtomPositions(poscarDir)
        self.setEnergy(poscarDir)
        
        self.writeDashedLine()
        self.writeIDString()
        self.writeLatticeVecs()
        self.writeAtomCounts()
        self.writeAtomPositions()
        self.writeEnergy()

    def makeUncleFiles(self):
        lastDir = os.getcwd()
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                os.chdir(atomDir)
                print "\nCreating structures.in and structures.holdout files for " + atom + "\n"
                self.initialize()
                self.infile = open('structures.in','w')
                self.holdoutFile = open('structures.holdout','w')
                self.setStructureList(atomDir)
                self.writeHeader()
                
                for structure in self.structList:
                    self.writePOSCAR(structure)
                
                self.closeOutFiles()
                os.chdir(lastDir)









