'''
Created on Aug 13, 2014

@author: eswens13
'''
from glob import glob
from math import sqrt
import os, subprocess


class Structs2Poscar:
    """ Converts the pseudo-POSCAR files that UNCLE gives us after extracting from the enumeration
        to standard POSCAR files that VASP can use.  The Converter class below is the one that 
        does most of the work in conversion.  This class is the administrative class, managing the 
        files and directories. """
        
    def __init__(self, atoms):
        """ CONSTRUCTOR """
        self.atoms = atoms
    
    def makeAtomDirectories(self):
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                subprocess.call('rm -r ' + atomDir + '/*', shell=True)
            else:
                subprocess.call(['mkdir',atomDir])
    
    def makePlots(self, plotDir):
        
        toPlot = os.listdir(self.getExportDir())
        for struct in toPlot:
            subprocess.call(['python', 'PlotGraphene.py', self.getPlotDir(), struct, '-u'])
    
    def convertOutputsToPoscar(self):
        self.makeAtomDirectories()
        
        for name in glob('vasp.0*'):
            structNum = str(self.retrieveStructNum(name))
            self.changeToPoscar(name)
            for atom in self.atoms:
                structDir = os.getcwd() + '/' + atom + '/' + structNum
                if os.path.isdir(structDir):
                    subprocess.call('rm -r ' + structDir + '/*', shell=True)
                else:
                    subprocess.call(['mkdir', structDir])
                
                subprocess.call(['cp','POSCAR',structDir])
            
            subprocess.call(['rm',name])
            subprocess.call(['rm','POSCAR'])
            
    def retrieveStructNum(self, structFile):
        fileChars = list(structFile)
        i = len(fileChars) - 1
        while fileChars[i] != '.':
            i -= 1

        return int(''.join(fileChars[i + 1:]))

    def changeToPoscar(self, structFile):
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
        
        # Get the number of each type of atom.
        nonCatomCounts = uncleLines[5].strip().split()
        nonCatomCounts = [int(count) for count in nonCatomCounts]
        
        self.atomCounts = [sum(nonCatomCounts), nonCatomCounts[0], nonCatomCounts[1]]
        if self.atomCounts[1] == 0 or self.atomCounts[2] == 0:
            self.pure = True
                
        positionLines = uncleLines[7:]
        
        self.set3dCpositionsFromDirectCoordinates(positionLines)
        self.set3dHpositionsFromDirectCoordinates(positionLines)
        self.set3dMpositionsFromDirectCoordinates(positionLines)
              
    def set2dCPositionsFromDirectCoordinates(self, positionLines):
        self._2dCpos = []
        for line in positionLines:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            
            x = comp1[0] + comp2[0]
            y = comp1[1] + comp2[1]
            
            self._2dCpos.append([x,y])
    
    def set2dHPositionsFromDirectCoordinates(self, positionLines):
        self._2dHpos = []
        for line in positionLines[:self.atomCounts[1]]:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            
            x = comp1[0] + comp2[0]
            y = comp1[1] + comp2[1]
            
            self._2dHpos.append([x,y])
        
    def set2dMPositionsFromDirectCoordinates(self, positionLines):
        self._2dMpos = []
        for line in positionLines[self.atomCounts[1]:]:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            
            x = comp1[0] + comp2[0]
            y = comp1[1] + comp2[1]
            
            self._2dMpos.append([x,y])
        
    def get2DDistance(self, atom1, atom2):
        xcomp = atom2[0] - atom1[0]
        ycomp = atom2[1] - atom1[1]
        
        return sqrt(pow(xcomp, 2) + pow(ycomp, 2))

    def set3dCpositionsFromDirectCoordinates(self, positionLines):
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

    def set3dCpositionsOld(self):
        print "In set3dCPositions"
    
        eps = .02
        done = []
        
        allList = []
        for pos in self._2dCpos:
            allList.append(pos)
        
        # Pick a point, make it above the plane, and add it to the list of points that are done.
        firstPos = [allList[0][0], allList[0][1], self.dC]
        done.append(firstPos)
        allList.remove(allList[0])
        
        doneInd = 0
        print "Beginning while loop"
        while len(allList) > 0:
            neighbors = []
            i = 0
            while i < len(allList):
                toCompare = [done[doneInd][0], done[doneInd][1]]
                compDistance = self.get2DDistance(toCompare, allList[i])
                if compDistance > self.distance - eps and compDistance < self.distance + eps:
                    neighbors.append(allList[i])
                    allList.remove(allList[i])
                else:
                    i += 1
            
            sizeChanged = False
            for pos in neighbors:
                if done[doneInd][2] > 0:
                    newPos = [pos[0], pos[1], -self.dC]
                    done.append(newPos)
                    sizeChanged = True
                else:
                    newPos = [pos[0], pos[1], self.dC]
                    done.append(newPos)
                    sizeChanged = True
            
            if sizeChanged:
                doneInd += 1
        
        print "Exiting while loop"
        
        self._3dCpos = []
        for pos in done:
            self._3dCpos.append(pos)
        
        print "Exiting set3dCPositions"
        
    def set3dHPositionsOld(self):
        self._3dHpos = []
        for hpos in self._2dHpos:
            for cpos in self._3dCpos:
                if cpos[0] == hpos[0] and cpos[1] == hpos[1]:
                    if cpos[2] > 0:
                        self._3dHpos.append([cpos[0], cpos[1], cpos[2] + self.dH])
                    else:
                        self._3dHpos.append([cpos[0], cpos[1], cpos[2] - self.dH])
    
    def set3dMPositionsOld(self):
        self._3dMpos = []
        for mpos in self._2dMpos:
            for cpos in self._3dCpos:
                if cpos[0] == mpos[0] and cpos[1] == mpos[1]:
                    if cpos[2] > 0:
                        self._3dMpos.append([cpos[0], cpos[1], cpos[2] + self.dM])
                    else:
                        self._3dMpos.append([cpos[0], cpos[1], cpos[2] - self.dM])

    def getLatticeVectors(self):
        return [self.lattVec1, self.lattVec2, self.lattVec3]

    def getAtomCounts(self):
        return self.atomCounts

    def getCPositions(self):
        return self._3dCpos
    
    def getHPositions(self):
        return self._3dHpos
    
    def getMPositions(self):
        return self._3dMpos

    def isPure(self):
        return self.pure
    
    



    