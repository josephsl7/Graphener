'''
Created on Aug 13, 2014


'''
from glob import glob
from numpy import sqrt,dot,rint
from numpy.linalg import inv,det
from copy import deepcopy
import os, subprocess
from comMethods import *

class structsToPoscar:
    """ Converts the pseudo-POSCAR files generated by the Extractor class to standard POSCAR files
        that VASP can use.  The Converter class below is the one that does most of the work in 
        conversion.  This class is the administrative class, managing the files and directories.  
        It creates directories for each of the metal atoms specified by the user and, within each 
        atom's directory, creates directories for each of the structures to be run through VASP
        calculations.  It populates these directories with the converted POSCAR file corresponding 
        to that structure. """
        
    def __init__(self, atoms, s2pStructList):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.s2pStructList = deepcopy(s2pStructList)
    
    def makePlots(self, plotDir):
        """ NOT WORKING CURRENTLY.  This method puts a visual representation of each structure in 
            its directory. """
        toPlot = os.listdir(self.getExportDir())
        for struct in toPlot:
            subprocess.call(['python', 'PlotGraphene.py', self.getPlotDir(), struct, '-u'])
    
    def convertOutputsToPoscar(self):
        """ Converts all the pseudo-POSCARs created by the Extractor class into POSCAR files that
            VASP can use to run calculations and puts them in a directory that will contain all
            the information for that structure. """
        for name in glob('enum/vasp.0*'):
            structNum = str(self.retrieveStructNum(name))
            if structNum>0:
                self.changeToPoscar(name)           
                for i in xrange(len(self.s2pStructList)):
                    if self.contains(structNum, self.s2pStructList[i]):
                        structDir = os.getcwd() + '/' + self.atoms[i] + '/' + structNum
                        if os.path.isdir(structDir):
                            subprocess.call('rm -r ' + structDir + '/*', shell=True)
                        else:
                            subprocess.call(['mkdir', structDir])                    
                        
                        subprocess.call(['cp','POSCAR',structDir])
                        self.s2pStructList[i].remove(str(structNum))
                
                subprocess.call(['rm',name])
                subprocess.call(['rm','POSCAR'])

            
    def retrieveStructNum(self, structFile):
        """ Returns the structure number from the name of the pseudo-POSCAR file.  For example,
            the file "vasp.000478" would return 478 as the structure number. """
        fileChars = list(structFile)
        i = len(fileChars) - 1
        while fileChars[i] != '.':
            i -= 1
        try:
            structNum = int(''.join(fileChars[i + 1:]))
        except:
            structNum = 0
        return structNum

    def contains(self, struct, alist):
        """ Retuns true if the list 'alist' contains the structure 'struct', false otherwise. """
        for i in xrange(len(alist)):
            if struct == alist[i]:
                return True
        
        return False

    def changeToPoscar(self, structFile):
        """ The hook for the Converter class.  Uses the Converter class to convert a pseudo-POSCAR
            file into a standard VASP POSCAR file. """
        infile = open(structFile, 'r')
        inLines = infile.readlines()
        infile.close()
        
        converter = Converter(structFile)
        converter.convert()
        
        poscar = open('POSCAR','w')
        
        if converter.isPure():
            if converter.getAtomCounts()[1] == 0:
                poscar.write("PURE M " + inLines[0])
            elif converter.getAtomCounts()[2] == 0:
                poscar.write("PURE H " + inLines[0])
        else:
            poscar.write(inLines[0])
            
        poscar.write('1.0\n')
        
        # These are in the 'scrambled', (z, x, y) form from UNCLE.  We need to unscramble them.
        latticevecs = converter.getLatticeVectors()
        
        for vec in latticevecs:
            poscar.write('%12.8f %12.8f %12.8f\n' % (vec[1], vec[2], vec[0]))
        
        atomCounts = converter.getAtomCounts()
        
        if converter.isPure():
            if atomCounts[1] == 0:
                poscar.write(str(atomCounts[0]) + ' ' + str(atomCounts[2]) + '\n')
            elif atomCounts[2] == 0:
                poscar.write(str(atomCounts[0]) + ' ' + str(atomCounts[1]) + '\n')
        else:
            poscar.write(str(atomCounts[0]) + ' ' + str(atomCounts[1]) + ' ' + str(atomCounts[2]) + '\n')
        
        poscar.write('Cartesian\n')
        
        Cpos = converter.getCPositions()
        for pos in Cpos:
            poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], pos[2]))
        
        Hpos = converter.getHPositions()
        for pos in Hpos:
            poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], pos[2]))
        
        Mpos = converter.getMPositions()
        for pos in Mpos:
            poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], pos[2]))
        
        poscar.close()
        

class Converter:
    """ Converts one of the files made from the 'makestr.x' routine in UNCLE to a POSCAR file
        that is ready to be used by VASP. """
        
    def __init__(self, uncleFile):
        """ CONSTRUCTOR """
        self.uncleFile = uncleFile
        self.atomCounts = []
        
        self.lattVec1 = []
        self.lattVec2 = []
        self.lattVec3 = []
        
        # 2D distance between hexagonal C atoms
        self.distance = 1.42085901
        
        # Distance from the plane when buckled
        self.dC = .22856
        
        # Bond distances from C atoms
        self.dH = 1.1
        self.dM = 2.2
        
        self._2dCpos = []
        self._3dCpos = []
        
        self._2dHpos = []
        self._3dHpos = []
        
        self._2dMpos = []
        self._3dMpos = []
        
        self.pure = False
    
    def convert(self):
        """ Runs through the whole conversion process for a single POSCAR file. """
        uncleFile = open(self.uncleFile, 'r')
        uncleLines = uncleFile.readlines()
        uncleFile.close()
        
        # Extract the lattice vectors.
        vectorLines = [line.strip().split() for line in uncleLines[2:5]]
        for i in xrange(len(vectorLines)):
            for j in xrange(len(vectorLines[i])):
                if list(vectorLines[i][j])[0] == '*':
                    vectorLines[i][j] = '15.0'
                    break
            
        self.lattVec1 = [float(vectorLines[0][0]), float(vectorLines[0][1]), float(vectorLines[0][2])]
        self.lattVec2 = [float(vectorLines[1][0]), float(vectorLines[1][1]), float(vectorLines[1][2])]
        self.lattVec3 = [float(vectorLines[2][0]), float(vectorLines[2][1]), float(vectorLines[2][2])]
        LV = zeros((3,3),dtype =float)
        LV[:,0] = self.lattVec1
        LV[:,1] = self.lattVec2
        LV[:,2] = self.lattVec3 
        cellVol = det(LV)
        PLV = array(  [[2.13128850,  -1.23050000,   0.00000000], 
                       [2.13128850,   1.23050000,   0.00000000], 
                       [0.00000000,   0.00000000,  15.00000000]])
        primCellVol = det(PLV)
        dxC = PLV[0,0]*2*0.333333333333333 #distance between 2 C atoms
        nCatoms = 2 * int(rint(cellVol/primCellVol))
        # Get the number of each type of atom.
        adatomCounts = uncleLines[5].strip().split()
        adatomCounts = [int(count) for count in adatomCounts] 
        nadatoms = sum(adatomCounts)      
        self.atomCounts = [nCatoms, adatomCounts[0], adatomCounts[1]]
        if self.atomCounts[1] == 0 or self.atomCounts[2] == 0:
            self.pure = True

        positionLines = uncleLines[7:]
        
        self.set3dCpositionsFromDirectCoordinates(positionLines)
        self.set3dHpositionsFromDirectCoordinates(positionLines)
        self.set3dMpositionsFromDirectCoordinates(positionLines)
#        print 'test', nCatoms, cellVol, primCellVol,nCatoms/cellVol,2.0/primCellVol
        if not isequal(nadatoms/cellVol,2.0/primCellVol): #need another C atom for each site (only one uncle site per cell)
            self.addCpositions(dxC)
            self.atomCounts = [nCatoms, adatomCounts[0], adatomCounts[1]]

    def addCpositions(self,dxC):
        '''Each primitive cell has 2 carbon atoms.  If uncle has only one site/primitive cell, then we have 
        to add the second carbon , which will not be a site for an adataom, and is below the plane
        if the first is above'''
        Cpos = deepcopy(self._3dCpos)
        for pos in Cpos:
            self._3dCpos.append([pos[0] + dxC, pos[1], -pos[2]])

            
    def get2DDistance(self, atom1, atom2):
        """ Returns the two-dimensional distance between two atoms. """
        xcomp = atom2[0] - atom1[0]
        ycomp = atom2[1] - atom1[1]
        
        return sqrt(pow(xcomp, 2) + pow(ycomp, 2))

    def set3dCpositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the C atoms given the 2D 'surface' positions from UNCLE.
            Another way to think of this is that it introduces the "buckling" into the sheet of
            C atoms. """
        self._3dCpos = []
        for line in positionLines:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            comp3 = [position[2] * self.lattVec3[0], position[2] * self.lattVec3[1], position[2] * self.lattVec3[2]]
            
            old_z = comp1[0] + comp2[0] + comp3[0]
            x = comp1[1] + comp2[1] + comp3[1]
            y = comp1[2] + comp2[2] + comp3[2]
            
            z = 0.0
            if old_z > 0:
                z = self.dC
            else:
                z = -self.dC
            
            self._3dCpos.append([x, y, z])
            
    def set3dHpositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the H atoms in the system. """
        self._3dHpos = []
        for line in positionLines[0:self.atomCounts[1]]:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            comp3 = [position[2] * self.lattVec3[0], position[2] * self.lattVec3[1], position[2] * self.lattVec3[2]]
            
            old_z = comp1[0] + comp2[0] + comp3[0]
            x = comp1[1] + comp2[1] + comp3[1]
            y = comp1[2] + comp2[2] + comp3[2]
            
            z = 0.0
            if old_z > 0:
                z = self.dC + self.dH
            else:
                z = -self.dC - self.dH
            
            self._3dHpos.append([x, y, z])

    def set3dMpositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the metal atoms of the system. """
        self._3dMpos = []
        for line in positionLines[self.atomCounts[1]:]:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            comp3 = [position[2] * self.lattVec3[0], position[2] * self.lattVec3[1], position[2] * self.lattVec3[2]]
            
            old_z = comp1[0] + comp2[0] + comp3[0]
            x = comp1[1] + comp2[1] + comp3[1]
            y = comp1[2] + comp2[2] + comp3[2]
            
            z = 0.0
            if old_z > 0:
                z = self.dC + self.dM
            else:
                z = -self.dC - self.dM
            
            self._3dMpos.append([x, y, z])
        
    def getLatticeVectors(self):
        """ Returns the lattice vectors of the structure. """
        return [self.lattVec1, self.lattVec2, self.lattVec3]

    def getAtomCounts(self):
        """ Returns the number of each atom in the structure. """
        return self.atomCounts

    def getCPositions(self):
        """ Returns the list of C positions in the structure. """
        return self._3dCpos
    
    def getHPositions(self):
        """ Returns the list of H positions in the structure. """
        return self._3dHpos
    
    def getMPositions(self):
        """ Returns the list of M positions in the structure. """
        return self._3dMpos

    def isPure(self):
        """ Returns true if the structure is a pure structure, false otherwise. """
        return self.pure
    
    



    