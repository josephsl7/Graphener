'''
Created on Aug 20, 2014

@author: eswens13
'''

from math import sqrt
from numpy import dot, transpose
from numpy.linalg.linalg import inv
import os, subprocess


class DistanceInfo:

    def __init__(self, atoms):
        """ CONSTRUCTOR """
    
        self.atoms = atoms
    
        self.structList = []
        self.setStructureList()
        
        self.inPlaneDistances = []
        self.normalDistances = []
        self.distances = []
        
        self.structureInfo = []
        self.origVecs = []
        self.relaxedVecs = []
        self.minOrigPositions = []
        self.relaxedPositions = []
        
    def setStructureList(self):
        atomDir = os.getcwd() + '/' + self.atoms[0]
        contents = os.listdir(atomDir)
        unsortedStructs = []
        for item in contents:
            if os.path.isdir(atomDir + '/' + item):
                structDir = atomDir + '/' + item
                poscar = open(structDir + '/POSCAR', 'r')
                poscarLines = [line.strip() for line in poscar]
                poscar.close()

                counts = [int(count) for count in poscarLines[5].split()]
            
                concentration = 0.0
                if poscarLines[0].split()[1] == 'H':
                    concentration = 0.0
                elif poscarLines[0].split()[1] == 'M':
                    concentration = 1.0
                else:
                    concentration = float(float(counts[2]) / float(counts[1] + counts[2]))
                
                unsortedStructs.append([concentration, item])
        
        unsortedStructs.sort()  # Now it's sorted by concentration
        
        for struct in unsortedStructs:
            self.structList.append(struct[1])
    
    def getStructureList(self):
        return self.structList
    
    def getAtomList(self):
        return self.atoms
    
    def getOriginalPositions(self, structureDir):
        # These are going to be in Cartesian coordinates
        
        fullStructPath = os.path.abspath(structureDir)
        poscar = open(fullStructPath + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        stringvecs = poscarLines[2:5]
        stringvecs = [line.strip().split() for line in stringvecs]
        
        self.origVecs = []
        for vec in stringvecs:
            newvec = [float(comp) for comp in vec]
            self.origVecs.append(newvec)
        
        counts = poscarLines[5].strip().split()
        counts = [int(count) for count in counts]
        
        self.structureInfo = []
        self.structureInfo.append(counts[0] / 2)
        self.structureInfo.append(counts)
    
        total = sum(counts)
        positionLines = poscarLines[7:7 + total]
        
        # Convert to Direct coordinates in terms of the !! NEW !! lattice vectors and then back to
        # Cartesian coordinates.
        positions = []
        for line in positionLines:
            newStringPosition = line.strip().split()
            newPosition = [float(comp) for comp in newStringPosition]
            
            directs = []
            r = dot(inv(transpose(self.relaxedVecs)), transpose(newPosition)) # Change to direct coordinates.
            
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        directs.append([r[0] + i, r[1] + j, r[2] + k])
            
            carts = []
            for pos in directs:
                rnew = dot(transpose(self.relaxedVecs), transpose(pos)) # Change back to cartesian coordinates.
                carts.append(rnew)
                
            positions.append(carts)
        
        return positions
        
    def getRelaxedPositions(self, structureDir):
        # These are going to be in Direct coordinates
        self.relaxedPositions = []
        
        fullStructPath = os.path.abspath(structureDir)
        poscar = open(fullStructPath + '/DOS/CONTCAR','r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        stringvecs = poscarLines[2:5]
        stringvecs = [line.strip().split() for line in stringvecs]
    
        self.relaxedVecs = []
        for vec in stringvecs:
            newvec = [float(comp) for comp in vec]
            self.relaxedVecs.append(newvec)
        
        counts = poscarLines[6].strip().split()
        counts = [int(count) for count in counts]
        
        total = sum(counts)
        positionLines = poscarLines[8:8 + total]
        
        newDirectPositions = []
        for line in positionLines:
            newPosition = line.strip().split()
            newDirectPositions.append([float(pos) for pos in newPosition])
        
        # Step 1:  Convert to Cartesian coordinates
        for position in newDirectPositions:
            
            rnew = dot(transpose(self.relaxedVecs), transpose(position))
            self.relaxedPositions.append(rnew)
        
        """# Step 2:  Convert to Direct coordinates in terms of the !! OLD !! lattice vectors.
        oldDirPos = []
        for position in newCartPos:
            rnew = dot(inv(transpose(self.origVecs)), transpose(position))
            oldDirPos.append(rnew)
        
        # Step 3:  Translate all the positions into the parallelepiped defined by the !! OLD !!
        #          lattice vectors.
        for r in oldDirPos:
            eps = 1e-4
            for i in xrange(3):
                while r[i] > 1.0 - eps or r[i] < 0.0 - eps:
                    if r[i] > 1.0 - eps:
                        r[i] = r[i] - 1
                    elif r[i] < 0.0 - eps:
                        r[i] = r[i] + 1
        
        # Step 4:  Convert back to Cartesian coordinates.    
        positions = []
        for position in oldDirPos:
            rnew = dot(transpose(self.origVecs), transpose(position)) # Transform to Cartesian coordinates
            positions.append(rnew)
            
        return positions """
    
    def getDistances(self, origPositions):
        if len(origPositions) != len(self.relaxedPositions):
            subprocess.call(['echo','\nERROR:  There are ' + str(len(origPositions)) + ' original positions and ' + str(len(self.relaxedPositions)) + ' relaxed positions.\n'])
            return None 
        else:
            self.minOrigPositions = []
            self.distances = []
            self.inPlaneDistances = []
            self.normalDistances = []
            for i in xrange(len(origPositions)):
                xr = self.relaxedPositions[i][0]
                yr = self.relaxedPositions[i][1]
                zr = self.relaxedPositions[i][2]
                
                trialDistances = []
                for j in xrange(len(origPositions[i])):
                
                    xo = origPositions[i][j][0]
                    yo = origPositions[i][j][1]
                    zo = origPositions[i][j][2]
                
                    newDistance = sqrt(pow((xr - xo), 2) + pow((yr - yo), 2) + pow((zr - zo), 2))
                    
                    trialDistances.append(newDistance)
                
                position = origPositions[i][self.find(min(trialDistances), trialDistances)]
                self.minOrigPositions.append(position)
                
                inPlaneDistance = sqrt(pow((xr - position[0]), 2) + pow((yr - position[1]), 2))
                normalDistance = sqrt(pow((zr - position[2]), 2))
                
                self.inPlaneDistances.append(inPlaneDistance)
                self.normalDistances.append(normalDistance)
                self.distances.append(min(trialDistances))
    
    def find(self, toFind, alist):
        if len(alist) == 0:
            return -1
        else:
            for i in xrange(len(alist)):
                if alist[i] == toFind:
                    return i
            
            return -1
    
    def getMCBondingInfo(self, structureDir):
        [Cindexes, Mcount] = self.getMCIndexes(structureDir)
        
        if len(Cindexes) != Mcount:
            subprocess.call(['echo','ERROR:  Did not retrieve the same number of M atoms and corresponding C atoms.'])
        else:
            contcar = open(structureDir + '/DOS/CONTCAR','r')
            contcarLines = [line.strip() for line in contcar]
            contcar.close()
            
            counts = [int(count) for count in contcarLines[6].strip().split()]
            allPositions = contcarLines[8:8 + sum(counts)]
            
            Cpos = []
            for i in Cindexes:
                stringPos = allPositions[i].split()
                floatPos = [float(comp) for comp in stringPos]
                cartPos = dot(transpose(self.relaxedVecs), transpose(floatPos))
                Cpos.append(cartPos)
            
            Mpos = []
            Mlines = allPositions[-Mcount:]
            for line in Mlines:
                position = [float(comp) for comp in line.split()]
                Mpos.append(position)
            
            bondLengths = []
            for i in xrange(len(Cpos)):
                trialMpos = []
                
                for x in [-1,0,1]:
                    for y in [-1,0,1]:
                        for z in [-1,0,1]:
                            newPos = [Mpos[i][0] + x, Mpos[i][1] + y, Mpos[i][2] + z]
                            trialMpos.append(dot(transpose(self.relaxedVecs), transpose(newPos)))
                
                cx = Cpos[i][0]
                cy = Cpos[i][1]
                cz = Cpos[i][2]
                
                trialLengths = []
                for trialPos in trialMpos:
                    mx = trialPos[0]
                    my = trialPos[1]
                    mz = trialPos[2]
                
                    length = sqrt(pow((mx - cx), 2) + pow((my - cy), 2) + pow((mz - cz), 2))
                    trialLengths.append(length)
                    
                bondLengths.append(min(trialLengths))
            
            return [min(bondLengths), max(bondLengths)]    
        
    def getMCIndexes(self, structureDir):
        origposcar = open(structureDir + '/POSCAR', 'r')
        poscarLines = [line.strip() for line in origposcar]
        origposcar.close()
        
        counts = poscarLines[5].strip().split()
        counts = [int(count) for count in counts]
        
        Ccount = counts[0]
        Mcount = 0
        if structureDir.split('/')[-1] == '3':
            Mcount = counts[1]
        else:
            Mcount = counts[2]
        
        Clines = poscarLines[7:7 + Ccount]
        Cpos = [line.split() for line in Clines]
        
        Mlines = poscarLines[-Mcount:]
        Mpos = [line.split() for line in Mlines]
        
        Cindexes = []
        for pos in Mpos:
            for i in xrange(len(Cpos)):
                if Cpos[i][0] == pos[0] and Cpos[i][1] == pos[1]:
                    Cindexes.append(i)
                    
        return [Cindexes, Mcount]
    
    def getBucklingInfo(self, structureDir):
        poscar = open(structureDir + '/POSCAR', 'r')
        poscarLines = [line.strip() for line in poscar]
        poscar.close()
        
        contcar = open(structureDir + '/DOS/CONTCAR','r')
        contcarLines = [line.strip() for line in contcar]
        contcar.close()
        
        poscarCounts = [int(count) for count in poscarLines[5].strip().split()]
        
        oldCpos = []
        for line in poscarLines[7:7 + poscarCounts[0]]:
            stringPos = line.split()
            floatPos = [float(comp) for comp in stringPos]
            oldCpos.append(floatPos)
        
        contcarCounts = [int(count) for count in contcarLines[6].strip().split()]
        
        newCpos = []
        for line in self.relaxedPositions[:contcarCounts[0]]:
            newCpos.append(line)
            
        NNpairs = self.getNNPairs(oldCpos)
        
        buckleDistances = []
        for pair in NNpairs:
            ind1 = pair[0]
            ind2 = pair[1]
            
            pos1 = newCpos[ind1]
            pos2 = newCpos[ind2]
            
            zpos1 = pos1[2]
            
            trialPos2 = []
            pos2Dir = dot(inv(transpose(self.relaxedVecs)), transpose(pos2))
            for i in [-1,0,1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        newPos = [pos2Dir[0] + i, pos2Dir[1] + j, pos2Dir[2] + k]
                        newCartPos = dot(transpose(self.relaxedVecs), transpose(newPos))
                        trialPos2.append(newCartPos)
            
            trialDistances = []
            for pos in trialPos2:
                zpos2 = pos[2]
                
                distance = abs(zpos1 - zpos2)
                trialDistances.append(distance)
            
            buckleDistances.append(min(trialDistances))
        
        return sum(buckleDistances) / len(buckleDistances)
    
    def getNNPairs(self, Cpos):
        eps = 1e-4
        pairs = []
        for i in xrange(len(Cpos)):
            for j in xrange(len(Cpos)):
                if i != j:
                    d = self.distance(Cpos[i], Cpos[j])
                    if d >= 1.42085899 - eps and d <= 1.42085899 + eps and not self.NNPairExists([i,j], pairs):
                        pairs.append([i,j])
        
        return pairs         
    
    def NNPairExists(self, pair, pairlist):
        if len(pairlist) == 0:
            return False
        else:
            for x in pairlist:
                if x[0] == pair[0] and x[1] == pair[1]:
                    return True
                elif x[0] == pair[1] and x[1] == pair[0]:
                    return True
                
            return False
        
    def distance(self, point1, point2):
        return sqrt(pow((point1[0] - point2[0]), 2) + pow((point1[1] - point2[1]), 2))
    
    def getMaxInPlaneDistance(self):
        return max(self.inPlaneDistances)
    
    def getMaxNormalDistance(self):
        return max(self.normalDistances)
                
    def getRMSDistance(self):
        n = len(self.distances)
        theSum = 0.0
        
        for d in self.distances:
            theSum = theSum + pow(d, 2)
        
        average = theSum / n
        
        rms = sqrt(average)
        
        return rms
  
    def writeInfoToFile(self, structureDir):
        fullpath = os.path.abspath(structureDir)
        outfile = open(fullpath + '/distance_info','w')
        
        structName = fullpath.split('/')[-1]
        
        self.getRelaxedPositions(structureDir)
        origPositions = self.getOriginalPositions(structureDir)
        
        self.getDistances(origPositions)
        
        outfile.write("*************************************************\n")
        outfile.write("      DISTANCE INFO FOR STRUCTURE " + structName + "\n")
        outfile.write("*************************************************\n\n")
        
        outfile.write("Volume Factor: " + str(self.structureInfo[0]) + "\t")
        
        concentration = 0.0
        if structName == '1':
            concentration = 0.0
        elif structName == '3':
            concentration = 1.0
        else:
            Hnum = self.structureInfo[1][1]
            Mnum = self.structureInfo[1][2]
            concentration = float(float(Mnum) / float(Hnum + Mnum))
        
        outfile.write("Concentration [M / (H + M)]: " + str(concentration) + "\n\n")
        
        rmsDistance = self.getRMSDistance()
        outfile.write("RMS distance moved in relaxation: " + str(rmsDistance) + "\n\n")
        
        outfile.write("Maximum distance moved in the plane: " + str(self.getMaxInPlaneDistance()) + "\n\n")
        
        outfile.write("Maximum distance moved perpendicular to the plane: " + str(self.getMaxNormalDistance()) + "\n\n")
        
        if structName == '1':
            outfile.write("Minimum M-C bond length: ------\n\n")
            outfile.write("Maximum M-C bond length: ------\n\n")
        else:
            [minMC, maxMC] = self.getMCBondingInfo(structureDir)
            outfile.write("Minimum M-C bond length: " + str(minMC) + "\n\n")
            outfile.write("Maximum M-C bond length: " + str(maxMC) + "\n\n")
        
        aveBuckle = self.getBucklingInfo(structureDir)
        outfile.write("Average buckling distance: " + str(aveBuckle) + "\n")
        
        outfile.write("\nOriginal Positions:\t\t\t\tRelaxedPositions:\n")
        for i in xrange(len(self.minOrigPositions)):
            o = self.minOrigPositions[i]
            r = self.relaxedPositions[i]
            outfile.write("%12.8f %12.8f %12.8f\t\t%12.8f %12.8f %12.8f\n" % (o[0], o[1], o[2], r[0], r[1], r[2]))
        
        outfile.close()

    def exportToCSV(self):
        topDir = os.getcwd()
        
        outfile = open('analysis/distance_summary.csv','w')
        
        outfile.write("Structure, vol. factor, M / (M + H), d - rms, Max d - para, Max d - perp, Min M - C, Max M - C, Ave Buckle\n")
        
        for atom in self.atoms:
            outfile.write("\n")
            outfile.write("\t, ,\t,*********************, ********** " + atom + " **********, *********************\n")
            atomDir = topDir + '/' + atom
            if os.path.isdir(atomDir):
                for structure in self.structList:
                    structDir = atomDir + '/' + structure
                    if os.path.isdir(structDir):
                        dfile = structDir + '/distance_info'
                        if os.path.isfile(dfile):
                            infile = open(dfile, 'r')
                            inlines = infile.readlines()
                            infile.close()
                            
                            volFactor = int(inlines[4].strip().split()[2])
                            concentration = float(inlines[4].strip().split()[9])
                            rms = float(inlines[6].strip().split()[5])
                            inplane = float(inlines[8].strip().split()[6])
                            normal = float(inlines[10].strip().split()[7])
                            if inlines[12].strip().split()[4] == '------':
                                minMC = 0.0
                            else:
                                minMC = float(inlines[12].strip().split()[4])
                            if inlines[14].strip().split()[4] == '------':
                                maxMC = 0.0
                            else:
                                maxMC = float(inlines[14].strip().split()[4])
                            buckle = float(inlines[16].strip().split()[3])
                            
                            outfile.write("%s, %d, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n" %
                                          (structure, volFactor, concentration, rms, inplane, normal, minMC, maxMC, buckle))
                            
                        else:
                            subprocess.call(['echo','ERROR:  The file ' + dfile + ' does not exist.'])
                    else:
                        subprocess.call(['echo','ERROR:  There is no directory ' + structDir])
            else:
                subprocess.call(['echo','ERROR:  There is no directory ' + atomDir])
                
        outfile.close()

    def getDistanceInfo(self):
        topDir = os.getcwd()
        for atom in self.getAtomList():
            atomDir = topDir + '/' + atom
            if os.path.isdir(atomDir):
                subprocess.call(['echo','********************'])
                subprocess.call(['echo','    ATOM ' + atom])
                subprocess.call(['echo','********************'])
                os.chdir(atomDir)
            
                for structure in self.structList:
                    structDir = atomDir + '/' + structure
                    if os.path.isdir(structDir):
                        DOSdir = structDir + '/DOS/'
                        if os.path.isdir(DOSdir):
                            subprocess.call(['echo','Working on structure ' + structure])
                            self.writeInfoToFile(structDir)
                        else:
                            subprocess.call(['echo','Structure ' + structure + ' did not converge. Skipping. . .'])
                    else:
                        subprocess.call(['echo','\nERROR: There is no directory ' + structDir + '\n'])
                os.chdir(topDir)
            else:
                subprocess.call(['echo','\nERROR: There is no directory ' + atomDir + '\n'])

        self.exportToCSV()







