'''
Created on Aug 26, 2014

@author: eswens13
'''
import os, subprocess


class Enumerator:
    """ Enumerates symmetrically unique structures using UNCLE.  After this class finishes its 
        work, the Extractor class will take over and extract pseudo-POSCAR files from the
        struct_enum.out file. """
  
    def __init__(self, atoms, volRange):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.volRange = volRange
        
        self.enumExec = os.path.abspath('needed_files/enum.x')
    
    def enumerate(self):
        subprocess.call(['mkdir','enum'])
        infile = open('needed_files/struct_enum.in','r')
        inlines = []
        for line in infile:
            inlines.append(line)
        infile.close()
        
        structFile = open('enum/struct_enum.in','w')
        for i in xrange(len(inlines)):
            if i == 9:
                structFile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + " ")
                structFile.write("# Starting and ending cell sizes for search\n")
            else:
                structFile.write(inlines[i])
        structFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        subprocess.call([self.enumExec,'struct_enum.in'])
        
        os.chdir(lastDir)
        
            
            
        