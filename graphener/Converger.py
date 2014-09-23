'''
Created on Sep 10, 2014

@author: eswens13
'''
import os

class Converger:


    def __init__(self, atoms):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.pairList = []
    
    def contains(self, alist, element):
        for item in alist:
            if item == element:
                return True
        
        return False
    
    def readGSSfile(self, atom):
        atomDir = os.getcwd() + '/' + atom
        
        gssfile = atomDir + '/gss/gss.out'
        self.structList = []
        
        infile = open(gssfile, 'r')
        gssLines = infile.readlines()
        infile.close()
        
        for line in gssLines[2:]:
            formEnergy = float(line.strip().split()[7])
            structNum = int(line.strip().split()[0])
            self.pairList.append([formEnergy, structNum])
        
        self.pairList.sort()
        
    def removeExistingStructs(self, atom):
        infile = open(os.getcwd() + '/' + atom + '/structures.in','r')
        inlines = infile.readlines()
        infile.close()
        
        existingStructs = []
        for i in xrange(len(inlines)):
            if i >= 3:
                if inlines[i].strip() == "#------------------------------------------------":
                    if inlines[i + 1].split.strip()[0] != 'PURE':
                        existingStructs.append(inlines[i + 1].strip().split()[3])
        
        for pair in self.pairList:
            structNum = str(pair[1])
            if self.contains(existingStructs, structNum):
                self.pairList.remove(pair)
        
        self.pairList.sort()
    
    def getLowestEnergyStructures(self, atom):
        self.readGSSfile(atom)
        self.removeExistingStructs(atom)
        
        structList = []
        for pair in self.pairList[:100]:
            structList.append(pair[1])
        
        return structList
            
        
        
            
        
        
        
        
        
        
        