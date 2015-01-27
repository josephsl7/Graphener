'''
Created on Aug 20, 2014


'''

from math import sqrt
from numpy import dot, transpose
from numpy.linalg.linalg import inv, norm
import os, subprocess
from comMethods import *
import StructsToPoscar

def correctz(directvec):
    '''Translates direct with z positions farther than 0.5 to the cell edge to closer than that'''
    if directvec[2] > 0.5:
        directvec[2] -= 1
    elif directvec[2] < -0.5:
        directvec[2] += 1
    return directvec
    

class DistanceInfo:

    def __init__(self, atoms,pureMetal,iteration):
        """ CONSTRUCTOR """
    
        self.atoms = atoms
    
        self.disStructList = []
        self.setStructureList()
        
        self.inPlaneDistances = []
        self.normalDistances = []
        self.distances = []
        
        self.structureInfo = []
        self.atomCounts = [] 
        self.origLattVecs = [] #all atoms positions
        self.relaxedLattVecs = [] #all atoms positions
        self.minOrigPositions = []
        self.relaxedPositions = []
        self.pureMetal = pureMetal
        self.iteration = iteration
        
    def setStructureList(self):
        '''Read finished structure names from vaspHFE.out '''
        
#        
#        
#        
#        atomDir = os.getcwd() + '/' + self.atoms[0]
#        contents = os.listdir(atomDir)
#        unsortedStructs = []
#        for item in contents:
#            if os.path.isdir(atomDir + '/' + item):
#                structDir = atomDir + '/' + item
#                poscar = open(structDir + '/POSCAR', 'r')
#                poscarLines = [line.strip() for line in poscar]
#                poscar.close()
#
#                counts = [int(count) for count in poscarLines[5].split()]
#            
#                concentration = 0.0
#                if poscarLines[0].split()[1] == 'H':
#                    concentration = 0.0
#                elif poscarLines[0].split()[1] == 'M':
#                    concentration = 1.0
#                else:
#                    concentration = float(float(counts[2]) / float(counts[1] + counts[2]))
#                
#                unsortedStructs.append([concentration, item])
#        
#        unsortedStructs.sort()  # Now it's sorted by concentration
#        
#        for struct in unsortedStructs:
#            self.disStructList.append(struct[1])
    
    def getStructureList(self):
        return self.disStructList
    
    def getatoms(self):
        return self.atoms
    
    def getOriginalPositions(self, struct):
        '''These are going to be in Cartesian coordinates. Go back to psuedoPOSCAR because POSCAR
        might have been overwritten by CONTCAR'''
        
        subprocess.call(['../needed_files/makestr.x','../enum/struct_enum.out',str(struct)])
        vfile = 'vasp.' + '0'*(6-len(str(struct))) + str(struct)
        toPoscar = StructsToPoscar.structsToPoscar([],[])
        toPoscar.convertOne(vfile) 
        subprocess.call(['rm', vfile])         
        poscarLines = readfile('POSCAR')
        subprocess.call(['rm', 'POSCAR'])       
        lattvecs = poscarLines[2:5]
        lattvecs = [line.strip().split() for line in lattvecs]

        self.origLattVecs = []
        for vec in lattvecs:
            newvec = [float(comp) for comp in vec]
            self.origLattVecs.append(newvec)
        counts = poscarLines[5].strip().split()
        self.atomCounts = [int(count) for count in counts]
        
        self.structureInfo = []
        self.structureInfo.append(self.atomCounts[0] / 2)
        self.structureInfo.append(self.atomCounts)
    
        total = sum(self.atomCounts)
        positionLines = poscarLines[7:7 + total]
        
        # Convert to Direct coordinates in terms of the !! NEW !! lattice vectors and then back to
        # Cartesian coordinates.
        positions = []
        for line in positionLines:
            newStringPosition = line.strip().split()
            newPosition = [float(comp) for comp in newStringPosition]
            
            directs = []
            r = dot(inv(transpose(self.relaxedLattVecs)), transpose(newPosition)) # Change to direct coordinates.
            
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        directs.append([r[0] + i, r[1] + j, r[2] + k])
            
            carts = []
            for pos in directs:
                rnew = dot(transpose(self.relaxedLattVecs), transpose(pos)) # Change back to cartesian coordinates.
                carts.append(rnew)
                
            positions.append(carts)
        
        return positions
        
    def getRelaxedPositions(self, structureDir):
        # These are going to be in Direct coordinates
        self.relaxedPositions = []
        
        fullStructPath = os.path.abspath(structureDir)
#        poscar = open(fullStructPath + '/CONTCAR','r')
        poscar = open(fullStructPath + '/CONTCAR','r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        lattvecs = poscarLines[2:5]
        lattvecs = [line.strip().split() for line in lattvecs]
    
        self.relaxedLattVecs = []
        for vec in lattvecs:
            newvec = [float(comp) for comp in vec]
            self.relaxedLattVecs.append(newvec)

        
        total = sum(self.atomCounts)
        positionLines = poscarLines[8:8 + total]
        
        newDirectPositions = []
        for line in positionLines:
            newPosition = line.strip().split()
            newDirectPositions.append([float(pos) for pos in newPosition])
        
        # Step 1:  Convert to Cartesian coordinates
        for position in newDirectPositions:
            
            rnew = dot(transpose(self.relaxedLattVecs), transpose(position))
            self.relaxedPositions.append(rnew)
        
        """# Step 2:  Convert to Direct coordinates in terms of the !! OLD !! lattice vectors.
        oldDirPos = []
        for position in newCartPos:
            rnew = dot(inv(transpose(self.origLattVecs)), transpose(position))
            oldDirPos.append(rnew)
        
        # Step 3:  Translate all the positions into the parallelepiped defined by the !! OLD !!
        #          lattice vectors.
        for r in oldDirPos:
            eps = 1e-4
            for i in range(3):
                while r[i] > 1.0 - eps or r[i] < 0.0 - eps:
                    if r[i] > 1.0 - eps:
                        r[i] = r[i] - 1
                    elif r[i] < 0.0 - eps:
                        r[i] = r[i] + 1
        
        # Step 4:  Convert back to Cartesian coordinates.    
        positions = []
        for position in oldDirPos:
            rnew = dot(transpose(self.origLattVecs), transpose(position)) # Transform to Cartesian coordinates
            positions.append(rnew)
            
        return positions """
    
    def getDistances(self):
        if len(self.origPositions) != len(self.relaxedPositions):
            subprocess.call(['echo','\nERROR:  There are ' + str(len(self.origPositions)) + ' original positions and ' + str(len(self.relaxedPositions)) + ' relaxed positions.\n'])
            return None 
        else:
            self.minOrigPositions = []
            self.distances = []
            self.inPlaneDistances = []
            self.normalDistances = []
            for i in range(len(self.origPositions)):
                xr = self.relaxedPositions[i][0]
                yr = self.relaxedPositions[i][1]
                zr = self.relaxedPositions[i][2]
                
                trialDistances = []
                for j in range(len(self.origPositions[i])):
                
                    xo = self.origPositions[i][j][0]
                    yo = self.origPositions[i][j][1]
                    zo = self.origPositions[i][j][2]
                
                    newDistance = sqrt(pow((xr - xo), 2) + pow((yr - yo), 2) + pow((zr - zo), 2))
                    
                    trialDistances.append(newDistance)
                
                position = self.origPositions[i][self.find(min(trialDistances), trialDistances)]
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
            for i in range(len(alist)):
                if alist[i] == toFind:
                    return i
            
            return -1
    
    def getMCBondingInfo(self, structureDir):
        '''Find distances to nearest C atom'''
        struct = structureDir.split('/')[-1]      
        contcar = open(structureDir + '/CONTCAR','r')
        contcarLines = [line.strip() for line in contcar]
        contcar.close()
        
        
        nC = self.atomCounts[0]
        if struct == self.pureMetal:
            nM = self.atomCounts[1]
        else:
            nM = self.atomCounts[2]
        allPositions = contcarLines[8:8 + sum(self.atomCounts)]
        
        Cpos = zeros((nC,3),dtype = float)
        for i in range(nC):
            stringPos = allPositions[i].split()
            floatPos = [float(comp) for comp in stringPos]
            floatPos = correctz(floatPos)
            cartPos = dot(transpose(self.relaxedLattVecs), transpose(floatPos))
            Cpos[i,:] = cartPos
        Mpos = zeros((nM,3),dtype = float)
        Mlines = allPositions[-nM:]
        for i,line in enumerate(Mlines):
            position = correctz([float(comp) for comp in line.split()])
            Mpos[i,:] = dot(transpose(self.relaxedLattVecs), transpose(position)) #convert to cartesian
#        for i in range(len(Mpos)):
#            trialMpos = []     D0 WE REALLY NEED COPIES TO TEST???????????????
#            for x in [-1,0,1]: #making copies in neighboring cells
#                for y in [-1,0,1]:
##                    for z in [-1,0,1]:
#                        z=0
#                        newPos = [Mpos[i,0] + x, Mpos[i,1] + y, Mpos[i,2] + z]
#                        trialMpos.append(dot(transpose(self.relaxedLattVecs), transpose(newPos)))        
        allClosest = []
        for Mvec in Mpos:
            trialLengths = []
            for Cvec in Cpos:
                trialLengths.append(norm(Mvec-Cvec))
            allClosest.append(min(trialLengths)) #closest for this M atom   
        return [min(allClosest), max(allClosest)]    
        
#    def getMCIndexes(self, structureDir):
#        '''Returns which C atoms match which M atoms in the original POSCAR (because they are on top sites'''
#        origposcar = open(structureDir + '/POSCAR', 'r')
#        poscarLines = [line.strip() for line in origposcar]
#        origposcar.close()
#        
#        counts = poscarLines[5].strip().split()
#        counts = [int(count) for count in counts]
#        
#        Ccount = counts[0]
#        Mcount = 0
#        if structureDir.split('/')[-1] == self.pureMetal:
#            Mcount = counts[1]
#        else:
#            Mcount = counts[2]
#        
#        Clines = poscarLines[7:7 + Ccount]
#        Cpos = [line.split() for line in Clines]
#        
#        Mlines = poscarLines[-Mcount:]
#        Mpos = [line.split() for line in Mlines]
#        
#        Cindexes = []
#        for pos in Mpos:
#            for i in range(len(Cpos)):
#                if Cpos[i][0] == pos[0] and Cpos[i][1] == pos[1]:
#                    Cindexes.append(i)
#                    
#        return [Cindexes, Mcount]
    
    def getBucklingInfo(self, structureDir):
        poscar = open(structureDir + '/POSCAR', 'r')
        poscarLines = [line.strip() for line in poscar]
        poscar.close()
        
        contcar = open(structureDir + '/CONTCAR','r')
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
        print 'NNpairs';print NNpairs
        print contcarCounts
        print newCpos        
        buckleDistances = []
        for pair in NNpairs:
            ind1 = pair[0]
            ind2 = pair[1]
            
            pos1 = newCpos[ind1]
            pos2 = newCpos[ind2]
            
            zpos1 = pos1[2]
            
            trialPos2 = []
            pos2Dir = dot(inv(transpose(self.relaxedLattVecs)), transpose(pos2))
            for i in [-1,0,1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        newPos = [pos2Dir[0] + i, pos2Dir[1] + j, pos2Dir[2] + k]
                        newCartPos = dot(transpose(self.relaxedLattVecs), transpose(newPos))
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
        for i in range(len(Cpos)):
            for j in range(len(Cpos)):
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
        self.origPositions = self.getOriginalPositions(structName)   
        self.getDistances()
        
        outfile.write("*************************************************\n")
        outfile.write("      DISTANCE INFO FOR STRUCTURE " + structName + "\n")
        outfile.write("*************************************************\n\n")
        
        outfile.write("Volume Factor: " + str(self.structureInfo[0]) + "\t")
        
        concentration = 0.0
        if structName == '1':
            concentration = 0.0
        elif structName == self.pureMetal:
            concentration = 1.0
        else:
            Hnum = self.structureInfo[1][1]
            Mnum = self.structureInfo[1][2]
            concentration = float(float(Mnum) / float(Hnum + Mnum))
        
        outfile.write("Concentration [x = M / (H + M)]: " + str(concentration) + "\n\n")
        
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
        for i in range(len(self.minOrigPositions)):
            o = self.minOrigPositions[i]
            r = self.relaxedPositions[i]
            outfile.write("%12.8f %12.8f %12.8f\t\t%12.8f %12.8f %12.8f\n" % (o[0], o[1], o[2], r[0], r[1], r[2]))
        
        outfile.close()

    def exportToCSV(self):
        topDir = os.getcwd()
        
        outfile = open('analysis/distance_summary.csv','w')
        
        outfile.write("Structure, vol. factor, M / (M + H), moved_rms, Max moved_para, Max moved_perp, Min d_MC, Max d_MC, Ave Buckle dz\n")
        
        for atom in self.atoms:
            outfile.write("\n")
            outfile.write("\t, ,\t,*********************, ********** " + atom + " **********, *********************\n")
            atomDir = topDir + '/' + atom
            if os.path.isdir(atomDir):
                for structure in self.disStructList:
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
    def getStructList(self,atomDir):
        structList = []
        flines = readfile(atomDir + '/gss/vaspFE_{}.out'.format(self.iteration))
        for line in flines:
            structList.append(line.split()[0])
        return structList

    def getDistanceInfo(self):
        topDir = os.getcwd()
        for atom in self.getatoms():
            atomDir = topDir + '/' + atom
            structList = self.getStructList(atomDir)
            if os.path.isdir(atomDir):
                subprocess.call(['echo','********************'])
                subprocess.call(['echo','    ATOM ' + atom])
                subprocess.call(['echo','********************'])
                os.chdir(atomDir)
                for structure in structList:
                    print structure
                    structDir = atomDir + '/' + structure
                    if structure =='447':
                        self.writeInfoToFile(structDir)
                os.chdir(topDir)
            else:
                subprocess.call(['echo','\nERROR: There is no directory ' + atomDir + '\n'])
            
#                for structure in structList:
#                    structDir = atomDir + '/' + structure
#                    if os.path.isdir(structDir):
#                        DOSdir = structDir + '/DOS/'
#                        if os.path.isdir(DOSdir):
#                            subprocess.call(['echo','Working on structure ' + structure])
#                            self.writeInfoToFile(structDir)
#                        else:
#                            subprocess.call(['echo','Structure ' + structure + ' did not converge. Skipping. . .'])
#                    else:
#                        subprocess.call(['echo','\nERROR: There is no directory ' + structDir + '\n'])
#                os.chdir(topDir)
#            else:
#                subprocess.call(['echo','\nERROR: There is no directory ' + atomDir + '\n'])

        self.exportToCSV()        

    def plotDisplPara(self):
        '''Plot the vaspHFE vs energy, but make the marker color reflect the max distance moved parallel to the plane.'''
        
#        gnuplot> set view map
#        gnuplot> set pm3d map
#        gnuplot> splot "test_map.dat" with points palette pt 9







