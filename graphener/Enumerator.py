'''
Created on Aug 26, 2014

@author: eswens13
'''
import os, subprocess


class Enumerator:
    """ Enumerates symmetrically unique structures using UNCLE and then extracts the pseudo-POSCAR
        files from the struct_enum.out file produced by UNCLE.  After this class finishes its 
        work, the Structs2Poscar class will take over and convert the pseudo-POSCAR files into
        standard POSCAR files. """
  
    def __init__(self, atoms):
        """ CONSTRUCTOR """
        self.atoms = atoms
        
        self.enumExec = '/fslhome/eswens13/compute/uncleStuff/theuncleApr2014/enumlib/trunk/enum.x'
        self.extractExec = '/fslhome/eswens13/compute/uncleStuff/theuncleApr2014/enumlib/trunk/makestr.x'
    
    def enumerate(self):
        subprocess.call(['cp','/fslhome/eswens13/compute/uncleStuff/theuncleApr2014/enumlib/trunk/training/struct_enum.in.gr8', os.getcwd()])
        subprocess.call([self.enumExec,'struct_enum.in.gr8'])
    
    def extract(self, structList):
        #TODO: Still need to figure out where this list is going to come from.
        
        for struct in structList:
            subprocess.call([self.extractExec, 'struct_enum.out',struct])
        
            
            
        