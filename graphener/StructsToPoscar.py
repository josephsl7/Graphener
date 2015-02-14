'''
Created on Aug 13, 2014


'''
from glob import glob
from numpy import sqrt,dot,rint,transpose,sign
from numpy.linalg import inv,det,norm
from copy import deepcopy
import os, subprocess
from comMethods import *
import shutil

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
        self.dC = .22856 #buckling distance...not essential
        # Bond distances from C atoms
        self.dH = 1.1
        self.dM = 2.2
    
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
                        shutil.rmtree(structDir)
                    subprocess.call(['mkdir', structDir])                    
            
                    infile = open('POSCAR', 'r')
                    inLines = infile.readlines()
                    infile.close()

                    poscar = open(structDir + '/POSCAR','w')

                    line = inLines[0]
                    elements = ['C'] + self.atoms[i].split(',')

                    poscar.write(line[:line.find(':') + 1] + ' ')

                    vacancies = 0

                    for num, element in enumerate(elements):
                        if self.atomCounts[num] != 0 and element.find('Vc') == -1 and num != 0:
                            poscar.write(str(num) + ' ')
                        elif element.find('Vc') != -1:
                            vacancies = self.atomCounts[num]

                    if vacancies == 0:
                        poscar.write(line[line.find('-'):])
                    else:
                        poscar.write(line[line.find('-'):-1].strip() + ' -- ' + 'Vacancies: ' + str(vacancies) + '\n')

                    poscar.write(inLines[1] + inLines[2] + inLines[3] + inLines[4])

                    countLine = ''        

                    for num, count in enumerate(self.atomCounts):
                        if count != 0 and elements[num].find('Vc') == -1:
                            countLine = countLine + str(count) + ' '
                    poscar.write(countLine.strip() + '\n')

                    poscar.write(inLines[6])

                    linenum = 7
                    for num, count in enumerate(self.atomCounts):
                        if count != 0:
                            for line in range(count):
                                if elements[num].find('Vc') == -1:
                                    pos = inLines[linenum].strip().split()
                                    pos = [float(comp) for comp in pos]
                                    z = pos[2]
                                    if elements[num] == 'H':
                                        if z == self.dC or z == -self.dC:
                                            z += sign(z)*self.dH
                                    elif elements[num] != 'C':
                                        if z == self.dC or z == -self.dC:
                                             z += sign(z)*self.dM
                                        elif z == self.dH or z == -self.dH:
                                            z = sign(z)*self.dM
                                    poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], z))
                                linenum += 1
                           
                    poscar.close()
                    subprocess.call(['cp', 'POSCAR', structDir + '/struct'])
            subprocess.call(['rm', name])
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
        self.struct = str(self.retrieveStructNum(structFile))
        converter.convert()

        self.atomCounts = converter.getAtomCounts()
        
        poscar = open('POSCAR','w')
        
        poscar.write("Contains atoms: ")
        for i, count in enumerate(self.atomCounts):
            if count != 0 and i != 0:
                poscar.write(str(i) + " ")
        
        if converter.isPure():
            poscar.write("- (PURE CASE) ")

        poscar.write("- " + inLines[0])
            
        poscar.write('1.0\n')

        latticevecs = converter.getLatticeVectors()
        
        for i in range(3):
            poscar.write('%12.8f %12.8f %12.8f\n' % (latticevecs[0,i], latticevecs[1,i], latticevecs[2,i]))
        
        countLine = ''        

        for count in self.atomCounts:
            if count != 0:
                countLine = countLine + str(count) + ' '
        poscar.write(countLine.strip() + '\n')
        
        poscar.write('Cartesian\n')
        
        Cpos = converter.getCPositions()
        for pos in Cpos:
            poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], pos[2]))
        
        Pos = converter.getPositions()
        for pos in Pos:
            poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], pos[2]))
        
        poscar.close()

class Converter:
    """ Converts one of the files made from the 'makestr.x' routine in UNCLE to a POSCAR file
        that is ready to be used by VASP. """
        
    def __init__(self, uncleFile):
        """ CONSTRUCTOR """
        self.uncleFile = uncleFile
        self.atomCounts = []
        self.struct = None
        self.LV = None  #lattice vectors
        self.Mmat = None
        self.PLV = None
        self.dC = .22856 #buckling distance...not essential
        # Bond distances from C atoms
        self.dH = 1.1
        self.dM = 2.2
        self.Cpos = []
        self.Pos = []
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
    
        self.atomCounts = [nCatoms] + adatomCounts
        if sum([count != 0 for count in adatomCounts]) == 1:
            self.pure = True
        positionLines = uncleLines[7:]
        self.getCpositions()
        self.getPositionsFromDirectCoordinates(positionLines)

    def getCpositions(self):
        '''Finds all the C positions within the unit cell'''
        eps = 1e-6
        dxC = self.PLV[0,0]*2*0.333333333333333 #distance between 2 C atoms
        pos1 = array([dxC,0,self.dC])  #Our 2 C atoms are at 1/3 and 2/3 along the x axis
        pos2 = array([2*dxC,0,-self.dC])  
        #find all graphene lattice points in region of unit cell defined by self.LV
        #the columns of M give us the 3 lattice vectors.
        reps = [] #copies of the two atoms in the primitive cell
        for iMcolumn in range(3): #for the 3 LVs        
            m0limit = int(self.Mmat[0,iMcolumn]+1*sign(self.Mmat[0,iMcolumn])) #3rd PLV is in z direction...don't need it.  
            m1limit = int(self.Mmat[1,iMcolumn]+1*sign(self.Mmat[1,iMcolumn]))
#            m0limit = 100; m1limit = 100
            extra = 2
            for m0 in range(min(0,m0limit)-extra,max(0,m0limit)+extra): 
                for m1 in range(min(0,m1limit)-extra,max(0,m1limit)+extra):
                    reps.append(pos1 + m0*self.PLV[:,0] + m1*self.PLV[:,1])
                    reps.append(pos2 + m0*self.PLV[:,0] + m1*self.PLV[:,1])
#        print reps
        positions = []
        for vec in reps:
            vecplanar = vec; vecplanar[2] = 0.0 # ignore z component for testing whether in unit cell
            directPos = dot(inv(self.LV), transpose(vecplanar)) # Change to direct coordinates, 
            if min(directPos) >= 0.0 - eps and max(directPos) < 1.0 - eps: #then in the unit cell
                exists = False
                for pos in positions:
                    if norm(pos-vec) < eps:
                        exists = True
                        break
                if not exists:
#                    print directPos
                    positions.append(vec)
        if len(positions) != self.atomCounts[0]:
            sys.exit('Error in getCpositions for {}: number of C atoms ({}) does not match area of cell ({}). Stopping'.format(self.struct,len(positions),self.atomCounts[0]))
        self.Cpos = positions # a list of numpy vectors
#        print 'C',self.Cpos
     
    def getPositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the  adatoms in the system. """
        self.Pos = []
        for line in positionLines[0:]:
            direct = line.strip().split()
            direct = array([float(comp) for comp in direct])
            cart = dot(self.LV, transpose(direct))
            #displace from plane
            onTop = False
            for pos in self.Pos: #check to see if this is on top of a C atom:
                if norm(pos[:2]-cart[:2])<0.1: #only test planar distance
                    onTop = True
                    break        
            if onTop: 
                cart[2] = sign(cart[2])* (self.dC) 
            else:
                cart[2] = sign(cart[2])* self.dH       
            self.Pos.append(cart)
#        print "H", self.Pos
        
    def getLatticeVectors(self):
        """ Returns the lattice vectors of the structure. """
        return self.LV

    def getAtomCounts(self):
        """ Returns the number of each atom in the structure. """
        return self.atomCounts

    def getCPositions(self):
        """ Returns the list of C positions in the structure. """
        return self.Cpos
    
    def getPositions(self):
        """ Returns the list of H positions in the structure. """
        return self.Pos
    
    def isPure(self):
        """ Returns true if the structure is a pure structure, false otherwise. """
        return self.pure
    
    



