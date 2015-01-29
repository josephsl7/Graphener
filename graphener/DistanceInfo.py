'''
Created on Aug 20, 2014


'''

from numpy import sqrt, array, mod,dot, transpose, std, size, zeros, where, int32
from numpy.linalg.linalg import inv, norm
import os, subprocess,sys
from comMethods import readfile,trimSmall
import StructsToPoscar
import matplotlib.pyplot as plt
  
class DistanceInfo:
    from comMethods import setAtomCounts
    def __init__(self, atoms,pureMetal,iteration):
        """ CONSTRUCTOR """
    
        self.atoms = atoms
    
        self.disStructList = []
#        self.setStructureList()
        
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
        self.origLattVecs = None 
        self.relaxLattVecs = None
        self.indexLVz = None
        self.hfe = []
        self.output = []

        self.pureMetal = pureMetal
        self.iteration = iteration
        self.struct = ''
        self.structDir = ''
        self.buckleAvg = 0.0
        self.conc = 0.0
        self.moveTotal = 0.0
        self.volFactor = 0
        self.moveInPlane = 0.0
        self.moveMInPlane = 0.0
        self.moveHInPlane = 0.0
        self.moveNormal = 0.0
        self.maxMC = 0.0
        self.minMC = 0.0
               
    def correctz(self,directvec):
        '''Translates direct vectors with z positions farther than 0.5 to the cell edge to closer than that,
        because we like to keep z centered around zero'''
        eps = 0.01
        if directvec[self.indexLVz] > 0.5 :
            directvec[self.indexLVz] -= 1
        elif directvec[self.indexLVz] < -0.5:
            directvec[self.indexLVz] += 1    
        return directvec
    
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
        self.origLattVecs = zeros((3,3),dtype = float)
        for i,vec in enumerate(lattvecs):
            newvec = [float(comp) for comp in vec]
            self.origLattVecs[i,:] = newvec
        positionLines = poscarLines[7:7 + self.nTot] # from vfile, these are already in direct lattice coordinates
        # Convert to Direct coordinates in terms of the !! NEW !! lattice vectors and then back to
        # Cartesian coordinates.
        positions = []
        for line in positionLines:
            position = line.strip().split()
            positions.append([float(comp) for comp in position])
        self.nomapOrigPos = array(positions) #need these to find NNpairs
        origPos = []
        for position in positions:
            #convert to direct
            directPos = dot(inv(transpose(self.origLattVecs)), transpose(position)) # Change to direct coordinates in the old cell
            origPos.append(dot(transpose(self.relaxLattVecs), transpose(directPos)))
        self.origPos = array(origPos)
        self.origCPos = self.origPos[:self.nC]
        self.origHPos = self.origPos[self.nC:self.nC+self.nH]
        self.origMPos = self.origPos[:-self.nM]
        
    def atomType(self,i):
        if 0 <= i < self.nC:
            atomtype = 'C'
        elif self.nC <= i < self.nC + self.nH:
            atomtype = 'H'
        elif self.nC + self.nH <= i < self.nTot:
            atomtype = 'M'
        return atomtype
        
    def getRelaxPos(self):
        '''Will be in Cartesian coordinates. Also gets atomCounts'''
        self.setAtomCounts(self.structDir)
        self.nC = self.atomCounts[0] 
        self.nH = self.atomCounts[1] 
        self.nM = self.atomCounts[2]
        self.nTot = sum(self.atomCounts) 
        self.volFactor = self.nC/2
        self.relaxLattVecs = zeros((3,3),dtype = float)
        contcarLines = readfile(self.structDir + '/CONTCAR')
        lattvecs = contcarLines[2:5]
        lattvecs = [line.strip().split() for line in lattvecs]    
        for i,vec in enumerate(lattvecs):
            newvec = [float(comp) for comp in vec]
            self.relaxLattVecs[i,:] = newvec
            cart = dot(transpose(self.relaxLattVecs), transpose(newvec))
            if abs(newvec[2]) > max(abs(cart[1]),abs(cart[1])):
                self.indexLVz = i #which LVector is in the z dir
        positionLines = contcarLines[8:8 + self.nTot]     
        relaxPos = []
        for line in positionLines:
            newPosition = line.strip().split()
            newPosition = [float(pos) for pos in newPosition]
            newPosition = self.correctz(newPosition)
            # Already direct from Contcar.  Convert to Cartesian coordinates
            relaxPos.append(dot(transpose(self.relaxLattVecs), transpose(newPosition)))
        self.relaxPos = array(relaxPos)

    def getMove(self):
        '''Movements of atoms'''         
        self.displace = self.relaxPos - self.origPos
        self.moveInPlane = zeros(self.nTot,dtype = float)
        self.moveNormal = zeros(self.nTot,dtype = float)
        self.moveTotal = zeros(self.nTot,dtype = float)
        self.moveMInPlane = zeros(self.nTot,dtype = float)
        self.moveHInPlane = zeros(self.nTot,dtype = float)
        
        for i in range(len(self.displace)):
            self.moveInPlane[i] = norm(self.displace[i,:2]) #only x,y parts
            self.moveNormal[i] = abs(self.displace[i,2])
            self.moveTotal[i] = norm(self.displace[i,:])
            if self.atomType(i) =='M': 
                self.moveMInPlane[i] = self.moveInPlane[i]
            elif self.atomType(i) =='H':
                self.moveHInPlane[i] = self.moveInPlane[i] 
            
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
        
        outfile.write("\nOriginal Pos : \t\t\t\t Relaxed Position: \t\t\t Movement\t\t delta-x \t delta-y \t delta-z\n")
        outfile.write("(mapped to new cell)\n")
        for i in range(len(self.origPos)):
            o = self.origPos[i,:]
            r = self.relaxPos[i,:]
            outfile.write("%s  %7.3f %7.3f %7.3f \t\t %7.3f %7.3f %7.3f \t\t %7.3f \t\t %7.3f \t %7.3f \t %7.3f\n" \
                          % (self.atomType(i), o[0], o[1], o[2], r[0], r[1], r[2],\
                              self.moveTotal[i], self.displace[i,0], self.displace[i,1], self.displace[i,2]))        
        outfile.close()
    
    def writeToCSV(self,csvfile):
        # "Structure, vol, conc, moved_rms, Max moved_para, Max moved_perp, Min d_MC, Max d_MC, Ave Buckle dz\n")
#        HFEindex = where(self.hfe['struct']==self.struct)
        HFEindex = self.hfe['struct'].tolist().index(int(self.struct))
        csvfile.write("%s, %d, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f\n" %
                      (self.struct, self.volFactor , self.conc , self.hfe[HFEindex]['HFE'], std(self.moveTotal), \
                       max(self.moveInPlane), max(self.moveMInPlane), max(self.moveHInPlane), max(self.moveNormal), self.minMC, self.maxMC, self.buckleAvg))

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
         We want to change that, because we can't easily find true small movements from initial positions
         that are close to the cell edge'''
        iplanar = [0,1,2]
        del iplanar[self.indexLVz]
        for i,relaxpos in enumerate(self.relaxPos):
            origpos =  self.origPos[i,:]           
            dmin = 100
            for j in range(-1,2):
                for k in range(-1,2): 
                    shift = j*self.relaxLattVecs[iplanar[0],:] + k*self.relaxLattVecs[iplanar[1],:]              
                    dtest = norm(relaxpos + shift - origpos)
                    if dtest<dmin:
                        dmin = dtest
                        bestshift = shift
            relaxpos = relaxpos + bestshift
            self.relaxPos[i,:] = relaxpos
        self.relaxCPos = self.relaxPos[:self.nC,:]
        self.relaxHPos = self.relaxPos[self.nC:self.nC+self.nH,:]
        self.relaxMPos = self.relaxPos[-self.nM:,:]


    def plotHFEMove(self,atomDir,moveType):
        '''Plot the vaspHFE vs energy, but make the marker color reflect the parameter moveType, 
        which is a string that is part of the movement.csv header, the substrings of which 
        make up self.output'''
        lines = readfile(atomDir + '/movement/movement.csv')
        print lines, atomDir + '/movement/movement.csv'
        print readfile('/panfs/pan.fsl.byu.edu/scr/usr/81/bch/top.tm_row1.v15/Cr_pv/movement/movement.csv')
 
        struct = []
        conc = []
        HFE = []
        move = []
        imove = self.output.index(moveType)
        for line in lines[1:]:
            print line
            struct.append(float(line.split(',')[0]))
            conc.append(float(line.split(',')[2]))
            HFE.append(float(line.split(',')[3]))
            move.append(float(line.split(',')[imove]))
        
        print 'imove',imove
        print move
        print conc
        #plot
        fig = plt.figure()
        plt.xlabel('{} concentration x'.format(atomDir.split('/')[0]))
        plt.ylabel('Energy (eV)')
        plt.title('Formation energy colored by movement')
#        plt.scatter(conc, HFE, c=move , cmap='autumn')
        plt.scatter(conc, HFE, c=move , cmap='jet')
        plt.colorbar()
        plt.show()
        fig.savefig(atomDir + '/gss/HFE_{}_{}'.format(moveType.replace(' ','_',),self.iteration))
        
    def getDistanceInfo(self):
        topDir = os.getcwd()
        for atom in self.atoms:
            atomDir = topDir + '/' + atom
            moveDir = atomDir + '/movement'
            if not os.path.exists(moveDir): os.mkdir(moveDir)
            csvfile = open(moveDir +'/movement.csv','w')
            header =  "Structure, vol, conc, HFE, moved_rms, Max moved_para, Max M moved_para, Max H moved_para, Max moved_perp, Min d_MC, Max d_MC, Ave Buckle dz"       
            csvfile.write(header + "\n")
            self.output = [item.strip() for item in header.split(',')]
            structList = self.getStructList(atomDir)
            if os.path.isdir(atomDir):
                subprocess.call(['echo','********************'])
                subprocess.call(['echo','    ATOM ' + atom])
                subprocess.call(['echo','********************'])
                os.chdir(atomDir)
                for istruct, structure in enumerate(structList[:3]):
                    if mod(istruct+1,100) == 0 or istruct+1 ==len(structList): subprocess.call(['echo','\tAnalyzed {} of {} structures for {}'.format(istruct+1,len(structList),atom)])
                    print structure
                    self.struct = structure
                    self.structDir = atomDir + '/' + structure
                    self.getRelaxPos()
                    self.getOrigPos()
                    self.correctRelaxPos()   
                    self.getMove()
                    self.getBucklingInfo()
                    self.writeInfoToFile()
                    self.writeToCSV(csvfile)
                self.plotHFEMove(atomDir,'Max M moved_para')
                os.chdir(topDir)
                
            else:
                subprocess.call(['echo','\nERROR: There is no directory ' + atomDir + '\n'])
            csvfile.close()







