'''
Created on Aug 20, 2014


'''

from numpy import sqrt, array, mod,dot, transpose, std, size, zeros, where, int32
from numpy.linalg.linalg import inv, norm
import os, subprocess
from comMethods import readfile,trimSmall
import StructsToPoscar
import matplotlib.pyplot as plt

def correctz(directvec):
    '''Translates direct vectors with z positions farther than 0.5 to the cell edge to closer than that,
    because we like to keep z centered around zero'''
    eps = 0.01
#    for i in range(2):
#        if directvec[i] > 1.0 + eps:
#            directvec[i] -= 1
#        elif directvec[i] < 0 - eps:
#            directvec[i] += 1
    if directvec[2] > 0.5 :
        directvec[2] -= 1
    elif directvec[2] < -0.5:
        directvec[2] += 1    
    return directvec
    
class DistanceInfo:
    from comMethods import setAtomCounts
    def __init__(self, atoms,pureMetal,iteration):
        """ CONSTRUCTOR """
    
        self.atoms = atoms
    
        self.disStructList = []
        self.setStructureList()
        
        self.relaxLattVecs = []
        self.inPlaneMove = []
        self.normalMove = []
        self.totalMove = []
        self.displace = None
        self.atomCounts = [] 
        self.nC = 0
        self.nH = 0
        self.nM = 0
        self.ntot = 0
        self.origCPos = []
        self.origMPos = []
        self.origHPos = []
        self.nomapOrigPos = [] #need to easily find NNs
        self.origPos = [] #all atoms' positions of current structure
        self.relaxPos = [] #all atoms' positions of current structure
        self.relaxCPos = []
        self.relaxMPos = []
        self.relaxHPos = []
        self.origLattVecs = [] 
        self.relaxLattVecs = []
        self.hfe = []

        self.pureMetal = pureMetal
        self.iteration = iteration
        self.struct = ''
        self.structDir = ''
        self.buckleAvg = 0.0
        self.conc = 0.0
        self.moveTotal = 0.0
        self.volFactor = 0
        self.moveInPlane = 0.0
        self.moveNormal = 0.0
        self.maxMC = 0.0
        self.minMC = 0.0
               
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
    
#    def getatoms(self):
#        return self.atom
    
    def getOrigPos(self):
        '''These are really the original positions mapped into the relaxed unit cell.  Will be in Cartesian coordinates. Go back to psuedoPOSCAR because POSCAR
        might have been overwritten by CONTCAR.  '''
        
        subprocess.call(['../needed_files/makestr.x','../enum/struct_enum.out',str(self.struct)])
        vfile = 'vasp.' + '0'*(6-len(str(self.struct))) + str(self.struct)
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
        positionLines = poscarLines[7:7 + self.nTot] # from vfile, these are already in direct lattice coordinates
        # Convert to Direct coordinates in terms of the !! NEW !! lattice vectors and then back to
        # Cartesian coordinates.
        directPos = []
        for line in positionLines:
            position = line.strip().split()
            directPos.append([float(pos) for pos in position])
        # Convert to Cartesian coordinates
        origPos = []
        nomapOrigPos = []
        for position in directPos:
            #convert to direct
            directPos = dot(inv(transpose(self.origLattVecs)), transpose(position)) # Change to direct coordinates in the old cell
#            directPos = correctz(directPos)
            nomapOrigPos.append(dot(transpose(self.origLattVecs), transpose(directPos)))
            origPos.append(dot(transpose(self.relaxLattVecs), transpose(directPos)))
 
#            r = dot(inv(transpose(self.origLattVecs)), transpose(newPosition)) # Change to direct coordinates in the old cell
#            r = correctz(r)
#            rnew = dot(transpose(self.relaxLattVecs), transpose(r)) # Map into new cell
#            origPos.append(rnew)
        self.nomapOrigPos = trimSmall(array(nomapOrigPos)) #need these to find NNpairs
        self.origPos = array(origPos)
        self.origCPos = self.origPos[:self.nC]
        self.origHPos = self.origPos[self.nC:self.nC+self.nH]
        self.origMPos = self.origPos[:-self.nM]
        
    def getRelaxPos(self):
        '''Will be in Cartesian coordinates. Also gets atomCounts'''
        self.setAtomCounts(self.structDir)
        self.nC = self.atomCounts[0] 
        self.nH = self.atomCounts[1] 
        self.nM = self.atomCounts[2]
        self.nTot = sum(self.atomCounts) 
        self.volFactor = self.nC/2
        contcarLines = readfile(self.structDir + '/CONTCAR')
        lattvecs = contcarLines[2:5]
        lattvecs = [line.strip().split() for line in lattvecs]    
        self.relaxLattVecs = []
        for vec in lattvecs:
            newvec = [float(comp) for comp in vec]
            self.relaxLattVecs.append(newvec)
        positionLines = contcarLines[8:8 + self.nTot]     
        newDirectPos = []
        for line in positionLines:
            newPosition = line.strip().split()
            newDirectPos.append([float(pos) for pos in newPosition])
        # Already direct from Contcar.  Convert to Cartesian coordinates
        relaxPos = []
        for position in newDirectPos:
            position = correctz(position)
            rnew = dot(transpose(self.relaxLattVecs), transpose(position))
            relaxPos.append(rnew)
        self.relaxPos = trimSmall(array(relaxPos))
        self.relaxCPos = self.relaxPos[:self.nC]
        self.relaxHPos = self.relaxPos[self.nC:self.nC+self.nH]
        self.relaxMPos = self.relaxPos[-self.nM:]

    def getMove(self):
        '''Movements of atoms''' 
        self.displace = self.relaxPos - self.origPos
        self.moveInPlane = zeros(self.nTot,dtype = float)
        self.moveNormal = zeros(self.nTot,dtype = float)
        self.moveTotal = zeros(self.nTot,dtype = float)
        for i in range(len(self.displace)):
            self.moveInPlane[i] = norm(self.displace[i,:2]) #only x,y parts
            self.moveNormal[i] = abs(self.displace[i,2])
            self.moveTotal[i] = norm(self.displace[i,:])
                   
    def find(self, toFind, alist):
        if len(alist) == 0:
            return -1
        else:
            for i in range(len(alist)):
                if alist[i] == toFind:
                    return i
            
            return -1
    
    def getMCBondingInfo(self):
        '''Find distances to nearest C atom''' 
        allClosest = []         
        for Mvec in self.relaxMPos:
            trialLengths = []
            for Cvec in self.relaxCPos:
                trialLengths.append(norm(Mvec-Cvec))
            allClosest.append(min(trialLengths)) #closest for this M atom 
        self.minMC = min(allClosest)
        self.maxMC = max(allClosest)
           
    def getBucklingInfo(self): 
        NNpairs = self.getNNPairs(self.nomapOrigPos[:self.nC])   
        buckleDistances = []
        for pair in NNpairs:
            buckleDistances.append(abs(self.relaxCPos[pair[0],2]-self.relaxCPos[pair[1],2])) #z components 
        self.buckleAvg =  sum(buckleDistances) / len(buckleDistances)
    
    def getNNPairs(self, Cpos):
        eps = 1e-4
        pairs = []
        for i in range(len(Cpos)):
            for j in range(len(Cpos)):
                if i != j:
                    d = norm(Cpos[i,:2] - Cpos[j,:2]) #take off the z components...looking for inplane distances
                    if abs(d - 1.42085899)<= eps and not self.NNPairExists([i,j], pairs):
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
        
#    def distance(self, point1, point2):
#        return sqrt(pow((point1[0] - point2[0]), 2) + pow((point1[1] - point2[1]), 2))
    
#    def getMaxInPlaneDistance(self):
#        return max(self.inPlaneMove)
#    
#    def getMaxNormalDistance(self):
#        return max(self.normalMove)
                
#    def getRMSDistance(self):
#        n = len(self.distances)
#        theSum = 0.0
#        
#        for d in self.distances:
#            theSum = theSum + pow(d, 2)
#        
#        average = theSum / n
#        
#        rms = sqrt(average)
#        
#        return rms
  
    def writeInfoToFile(self):
        outfile = open(self.structDir + '/movement_info','w')       
        outfile.write("*************************************************\n")
        outfile.write("      DISTANCE INFO FOR STRUCTURE " + self.struct + "\n")
        outfile.write("*************************************************\n\n")
        
        outfile.write("Volume Factor: " + str(self.volFactor) + "\t")
        
        self.conc = float(float(self.nM) / float(self.nH + self.nM))
        
        outfile.write("Metal concentration: " + str(self.conc) + "\n\n")
        
        outfile.write("RMS distance moved in relaxation: " + str(std(self.moveTotal)) + "\n\n")
        
        outfile.write("Maximum distance moved in the plane: " + str(max(self.moveInPlane)) + "\n\n")
        
        outfile.write("Maximum distance moved perpendicular to the plane: " + str(max(self.moveNormal)) + "\n\n")
        
        if self.struct == '1':
            outfile.write("Minimum M-C bond length: NONE \n\n")
            outfile.write("Maximum M-C bond length: NONE \n\n")
        else:
            self.getMCBondingInfo()
            outfile.write("Minimum M-C bond length: " + str(self.minMC) + "\n\n")
            outfile.write("Maximum M-C bond length: " + str(self.maxMC) + "\n\n")
        
        
        outfile.write("Average buckling distance: " + str(self.buckleAvg) + "\n")
        
        outfile.write("\nOriginal Pos : \tRelaxed Position: \t\t\t Movement\n")
        outfile.write("(mapped to new cell)\n")
        for i in range(len(self.origPos)):
            o = self.origPos[i]
            r = self.relaxPos[i]
            if 0 <= i < self.nC:
                atomtype = 'C'
            elif self.nC <= i < self.nC + self.nH:
                atomtype = 'H'
            elif self.nC + self.nH <= i < self.nTot:
                atomtype = 'M'
            outfile.write("%s  %7.3f %7.3f %7.3f \t\t %7.3f %7.3f %7.3f \t\t %7.3f\n" \
                          % (atomtype, o[0], o[1], o[2], r[0], r[1], r[2], self.moveTotal[i]))        
        outfile.close()
    
    def writeToCSV(self,csvfile):
        # "Structure, vol, conc, moved_rms, Max moved_para, Max moved_perp, Min d_MC, Max d_MC, Ave Buckle dz\n")
#        HFEindex = where(self.hfe['struct']==self.struct)
        HFEindex = self.hfe['struct'].tolist().index(int(self.struct))
        csvfile.write("%s, %d, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f\n" %
                      (self.struct, self.volFactor , self.conc , self.hfe[HFEindex]['HFE'], std(self.moveTotal), \
                       max(self.moveInPlane), max(self.moveNormal), self.minMC, self.maxMC, self.buckleAvg))

    def getStructList(self,atomDir):
        structList = []
        hlines = readfile(atomDir + '/gss/vaspHFE_{}.out'.format(self.iteration))
        for line in hlines:
            structList.append(line.split()[0])
        self.hfe = zeros((len(structList)),dtype = [('struct', int32),('HFE', float)])
        for i,line in enumerate(hlines):
            self.hfe[i]['struct'] = line.split()[0]
            self.hfe[i]['HFE'] = line.split()[2]
        return structList

    def correctRelaxPos(self):
        '''Contcar moves any position outside of cell into cell.  
         We want to change that, because we can't easily find movements from initial positions
         that are close to the cell edge'''
        for i,relaxpos in enumerate(self.relaxPos):
             #convert to direct
             relaxpos = dot(inv(transpose(self.relaxLattVecs)), transpose(relaxpos)) # Change to direct coordinates in the old cell
             origpos = dot(inv(transpose(self.origLattVecs)), transpose(self.origPos[i]))
             dtest = zeros(4,dtype = float)
             dmin = 100
             for j in range(-2,2):
                 for k in range(-2,2):                 
                     dtest = norm(relaxpos - origpos + array([j,k,0]))
                     if dtest<dmin:
                         dmin = dtest
                         bestshift = array([j,k,0])
             relaxpos = relaxpos + bestshift
             #convert to cartesian
             self.relaxPos[i] = dot(transpose(self.relaxLattVecs), transpose(relaxpos))            

    def getDistanceInfo(self):
        topDir = os.getcwd()
        for atom in self.atoms:
            atomDir = topDir + '/' + atom
            moveDir = atomDir + '/movement'
            if not os.path.exists(moveDir): os.mkdir(moveDir)
            csvfile = open(moveDir +'/movement.csv','w')        
            csvfile.write("Structure, vol, conc, HFE, moved_rms, Max moved_para, Max moved_perp, Min d_MC, Max d_MC, Ave Buckle dz\n")
            structList = self.getStructList(atomDir)
            if os.path.isdir(atomDir):
                subprocess.call(['echo','********************'])
                subprocess.call(['echo','    ATOM ' + atom])
                subprocess.call(['echo','********************'])
                os.chdir(atomDir)
                for istruct, structure in enumerate(structList):
                    if mod(istruct+1,100) == 0 or istruct+1 ==len(structList): subprocess.call(['echo','\tAnalyzed {} of {} structures for {}'.format(istruct+1,len(structList),atom)])
                    print structure
                    self.struct = structure
                    self.structDir = atomDir + '/' + structure
                    #if structure =='447':
                    self.getRelaxPos()
                    self.getOrigPos()
                    self.correctRelaxPos()   
                    self.getMove()
                    self.getBucklingInfo()
                    self.writeInfoToFile()
                    self.writeToCSV(csvfile)
                plotHFEMovePara(atomDir)
                os.chdir(topDir)
                
            else:
                subprocess.call(['echo','\nERROR: There is no directory ' + atomDir + '\n'])
            csvfile.close()       

    def plotHFEMovePara(self,atomDir):
        '''Plot the vaspHFE vs energy, but make the marker color reflect the max distance moved parallel to the plane.'''
        lines = readfile(atomDir + '/movement/movement.csv')
        struct = []
        conc = []
        HFE = []
        move = []
        for line in lines[1:]:
            struct.append(float(line.split(',')[0]))
            conc.append(float(line.split(',')[2]))
            HFE.append(float(line.split(',')[3]))
            move.append(float(line.split(',')[5]))
        #plot
        fig = plt.figure()
        plt.xlabel('{} concentration x'.format(atomDir.split('/')[0]))
        plt.ylabel('Energy (eV)')
        plt.title('Formation energy coded by in-plane movement')
        plt.scatter(conc, HFE, c=move , cmap='autumn')
        plt.scatter(conc, HFE, c=move , cmap='jet')
        plt.colorbar()
        plt.show()
        fig.savefig(atomDir + '/gss/HFE_moveInPlane_{}'.format(self.iteration))
     
        









