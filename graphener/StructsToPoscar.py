'''
Created on Aug 13, 2014


'''
from glob import glob
from numpy import sqrt,dot,rint,transpose,sign
from numpy.linalg import inv,det,norm
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
    
    def convertOne(self,name):
        self.changeToPoscar(name)

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
        """ Returns true if the list 'alist' contains the structure 'struct', false otherwise. """
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
        print 'struct',structFile
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

        latticevecs = converter.getLatticeVectors()
        
        for i in range(3):
            poscar.write('%12.8f %12.8f %12.8f\n' % (latticevecs[0,i], latticevecs[1,i], latticevecs[2,i]))
        
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
        
        self.LV = None  #lattice vectors
        self.Mmat = None
        self.PLV = None
        self.dC = .22856 #buckling distance...not essential
        # Bond distances from C atoms
        self.dH = 1.1
        self.dM = 2.2
        self.Cpos = []
        self.Hpos = []
        self.Mpos = []
        self.pure = False
    
    def convert(self):
        """ Runs through the whole conversion process for a single POSCAR file. """
        uncleFile = open(self.uncleFile, 'r')
        uncleLines = uncleFile.readlines()
        uncleFile.close()
        # Extract the lattice vectors.
        vectorLines = [line.strip().split() for line in uncleLines[2:5]]          
        vec1comps = [float(comp) if comp[0]!='*' else 15.0 for comp in vectorLines[0]] #to handle ******
        vec2comps = [float(comp) if comp[0]!='*' else 15.0 for comp in vectorLines[1]]
        vec3comps = [float(comp) if comp[0]!='*' else 15.0 for comp in vectorLines[2]]
        self.LV = zeros((3,3),dtype =float)
        self.LV[:,0] = array([vec1comps[1],vec1comps[2],vec1comps[0]]) #store lattice vectors as columns; switch from zxy to xyz
        self.LV[:,1] = array([vec2comps[1],vec2comps[2],vec2comps[0]])  
        self.LV[:,2] = array([vec3comps[1],vec3comps[2],vec3comps[0]]) 
        cellVol = det(self.LV)
        self.PLV = array(  [[2.13128850,  -1.23050000,   0.00000000],  #parent lattice vectors as rows
                       [2.13128850,   1.23050000,   0.00000000], 
                       [0.00000000,   0.00000000,  15.00000000]])
        self.PLV = transpose(self.PLV) #vectors as columns
        primCellVol = det(self.PLV)
        nCatoms = 2 * int(rint(cellVol/primCellVol))
        #get the integer matrix M that defines the lattice vectors in terms of the primitive ones:
        # LV = PLV * M
        self.Mmat = dot(inv(self.PLV),self.LV)   
        # Get the number of each type of atom.
        adatomCounts = uncleLines[5].strip().split()
        adatomCounts = [int(count) for count in adatomCounts] 
        nadatoms = sum(adatomCounts)      
        self.atomCounts = [nCatoms, adatomCounts[0], adatomCounts[1]]
        if self.atomCounts[1] == 0 or self.atomCounts[2] == 0:
            self.pure = True
        positionLines = uncleLines[7:]
        self.getCpositions()
        self.getHpositionsFromDirectCoordinates(positionLines)
        self.getMpositionsFromDirectCoordinates(positionLines)

    def getCpositions(self):
        '''Finds all the C positions within the unit cell'''
        dxC = self.PLV[0,0]*2*0.333333333333333 #distance between 2 C atoms
        pos1 = array([dxC,0,self.dC])  #Our 2 C atoms are at 1/3 and 2/3 along the x axis
        pos2 = array([2*dxC,0,-self.dC])  
        #find all graphene lattice points in region of unit cell defined by self.LV
        #the columns of M give us the 3 lattice vectors.
        reps = [] #copies of the two atoms in the primitive cell
        for iMcolumn in range(3): #for the 3 LVs        
            m0limit = int(self.Mmat[0,iMcolumn]+1*sign(self.Mmat[0,iMcolumn])) #3rd PLV is in z direction...don't need it.  
            m1limit = int(self.Mmat[1,iMcolumn]+1*sign(self.Mmat[1,iMcolumn]))
            for m0 in range(min(0,m0limit),max(0,m0limit)+1): 
                for m1 in range(min(0,m1limit),max(0,m1limit)+1):
                    reps.append(pos1 + m0*self.PLV[:,0] + m1*self.PLV[:,1])
                    reps.append(pos2 + m0*self.PLV[:,0] + m1*self.PLV[:,1])
        print reps
        positions = []
        for vec in reps:
            vecplanar = vec; vecplanar[2] = 0.0 # ignore z component for testing whether in unit cell
            directPos = dot(inv(self.LV), transpose(vecplanar)) # Change to direct coordinates, 
            if max(directPos)< 1.0 and min(directPos)>=0.0: #then in the unit cell
                exists = False
                for pos in positions:
                    if norm(pos-vec)<1e-6:
                        exists = True
                        break
                if not exists:
                    positions.append(vec)
        if len(positions) != self.atomCounts[0]:
            sys.exit('Error in getCpositions: number of C atoms does not match area of cell. Stopping')
        self.Cpos = positions # a list of numpy vectors
            
    def getHpositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the H atoms in the system. """
        self.Hpos = []
        for line in positionLines[0:self.atomCounts[1]]:
            direct = line.strip().split()
            direct = array([float(comp) for comp in direct])
            cart = dot(self.LV, transpose(direct))
            print cart
            #displace from plane
            onTop = False
            for pos in self.Cpos: #check to see if this is on top of a C atom:
                if norm(pos[:2]-cart[:2])<0.1: #only test planar distance
                    onTop = True
                    break        
            if onTop: 
                cart[2] = sign(cart[2])* (self.dC + self.dH) 
            else:
                cart[2] = sign(cart[2])* self.dH       
            self.Hpos.append(cart)
        print "H", self.Hpos


    def getMpositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the M atoms in the system. """
        self.Mpos = []
        for line in positionLines[self.atomCounts[1]:]:
            direct = line.strip().split()
            direct = array([float(comp) for comp in direct])
            cart = dot(self.LV, transpose(direct))
            print cart
            #displace from plane
            onTop = False
            for pos in self.Cpos: #check to see if this is on top of a C atom:
                if norm(pos[:2]-cart[:2])<0.1: #only test planar distance
                    onTop = True
                    break        
            if onTop: 
                cart[2] = sign(cart[2])* (self.dC + self.dM) 
            else:
                cart[2] = sign(cart[2])* self.dM
            self.Mpos.append(cart)
        print "M", self.Mpos 
        
    def getLatticeVectors(self):
        """ Returns the lattice vectors of the structure. """
        return self.LV

    def getAtomCounts(self):
        """ Returns the number of each atom in the structure. """
        return self.atomCounts

    def getCPositions(self):
        """ Returns the list of C positions in the structure. """
        return self.Cpos
    
    def getHPositions(self):
        """ Returns the list of H positions in the structure. """
        return self.Hpos
    
    def getMPositions(self):
        """ Returns the list of M positions in the structure. """
        return self.Mpos

    def isPure(self):
        """ Returns true if the structure is a pure structure, false otherwise. """
        return self.pure