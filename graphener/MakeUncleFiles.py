'''
Created on Aug 29, 2014

@author: eswens13
'''
from numpy import zeros, mod, count_nonzero 
import os, subprocess, sys
from random import random


class MakeUncleFiles:


    def __init__(self, atoms, startFromExisting, iteration,finalDir):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.startFromExisting = startFromExisting
        self.iteration = iteration
        self.finalDir = finalDir   
        
        self.newlyFinished = [[]*len(atoms)]  #list of structures, not paths
        self.newlyFailed = [[]*len(atoms)]
        self.pureHenergy = 0.0
        self.pureMenergy = 0.0
        
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0   
        self.singleE = [] 
        self.hexE = [] 
        self.vdata = []

    def analyzeNewVasp(self,vstructsCurrent):
        """ Initializes the list of structures to add to the structures.in and structures.holdout
            files by adding only the structures that VASP finished to the member
            newlyFinished. Sorts the list by metal concentration (N_M / N_total). Adds the structures that
            failed VASP calculations to the member 'failedList'. """   
        lastDir = os.getcwd()
        self.newlyFinished = [[]*len(self.atoms)]
        self.newlyFailed = [[]*len(self.atoms)]
        
        for iatom, atom in enumerate(self.atoms):
            # If it is the first iteration and we are starting from an existing structures.start
            # file, we just append an empty list to the newlyFinished and the failedList.  Else, 
            # proceed as normal.
            if self.iteration == 1 and self.startFromExisting[iatom]:
                'Do nothing...'
            else:
                atomDir = lastDir + '/' + atom
                os.chdir(atomDir)
#                pureHdir = os.getcwd() + '/1'
#                pureMdir = os.getcwd() + '/3'
#            
#                if os.path.exists(pureHdir):
#                    self.setAtomCounts(pureHdir)
#                    self.setEnergy(pureHdir)
#                    self.pureHenergy = float(self.energy)
#                else:
#                    etest = self.getPureEnergyFromExisting('H')
#                    if etest != 999999:
#                        self.pureHenergy = etest
#                    else:
#                        subprocess.call(['echo', '\nERROR:  The pure H structure is not part of structures.start for ' + atom + '.\n'])
#            
#                if os.path.exists(pureMdir):
#                    self.setAtomCounts(pureMdir)
#                    self.setEnergy(pureMdir)
#                    self.pureMenergy = float(self.energy)
#                else:
#                    etest = self.getPureEnergyFromExisting('M')
#                    if etest != 999999:
#                        self.pureMenergy = etest
#                    else:
#                        subprocess.call(['echo', '\nERROR:  The pure M structure is not part of structures.start for ' + atom + '.\n'])
#            
                conclist = []
                finished = []
                failed = []
                for i, struct in enumerate(vstructsCurrent[iatom]):
                    if mod(i+1,100) == 0: print 'Checking',i+1,'of',len(dirList2), 'structures in', atom  
#                    fullPath = os.path.abspath(struct)
                    if os.path.isdir(struct):
                        if os.path.isdir(self.finalDir):
                            if self.FinishCheck(struct + self.finalDir) and self.convergeCheck(struct + self.finaldir, self.getNSW(struct + self.finaldir)): #finaldir                           
                               # Check for concentration
                                self.setAtomCounts(struct)                            
                                concentration = 0.0
                                if self.atomCounts[0] == 0:
                                    concentration = 1.0
                                else:
                                    concentration = float(float(self.atomCounts[1]) / float(self.atomCounts[0] + self.atomCounts[1]))                           
                                conclist.append([concentration, struct])
                            else:
                                failed.append(struct)
                        else:
                            subprocess.call(['echo', '\nERROR: directory does not exist: ' + struct]) 
                self.newlyFailed[iatom] = failed #for atom i                
                conclist.sort()                
                for i in xrange(len(conclist)):
                    finished.append(conclist[i][1]) 
                self.newlyFinished[iatom]= finished #sorted by concentration, for atomi
                os.chdir(lastDir)

    def contains(self, struct, alist):
        """ Returns True if the list 'alist' contains the structure 'struct', False otherwise. """
        if len(alist) <= 0:
            return False

        for elem in alist:
            if str(struct) == str(elem):
                return True

        return False
   
    def convergeCheck(self, folder, NSW):
        """ Tests whether force convergence is done by whether the last line of OSZICAR (the last
            ionic relaxation step) is less than NSW."""
        try:
            value = self.getSteps(folder)
            return value < NSW and self.energyDropCheck(folder)
        except:
            return False  # True/False

    def elConvergeCheck(self,folder,NELM):  
        """Tests electronic convergence is done by whether the electronic step is less than NELM."""
        try:
            value = self.getElSteps(folder)
            return value < NELM #True/False
        except:
            return False #True/False
                        
    def electronicConvergeFinish(self,dir): 
        '''Test requires electronic convergence AND vasp finishing'''
        #get NELM, the max electronic steps allowed in the run.  Using first directory in dirslist
        proc = subprocess.Popen(['grep','-i','NELM',dir+'/INCAR'],stdout=subprocess.PIPE)
        result =  proc.communicate()[0]
        NELM = int(result.split('=')[1].split()[0])
        return self.elConvergeCheck(dir,NELM) and self.FinishCheck(dir)
            
    def energyDropCheck(self,dir):
        '''tests whether the energies have dropped in OSZICAR...rising energies are unphysical 
        and show a numerical convergence problem. The factor like 0.99 allows for very small energy rises only'''
        lines = self.readfile(dir+'/OSZICAR') 
        energies = []
        for line in lines:
            if 'F=' in line:
                energies.append(float(line.split()[2]))
        return energies[-1] <= 0.99*energies[0] 
    
    def FinishCheck(self, folder):
        """ Tests whether VASP is done by finding "Voluntary" in last line of OUTCAR.  The input
            parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """   
        lastfolder = os.getcwd()
        os.chdir(folder)        
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'], stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)            
        return newstring[0].find('Voluntary') > -1

    def getElSteps(self,folder): 
        '''number of electronic steps for isolated runs'''
        lastfolder = os.getcwd()
        os.chdir(folder)
        try:
            oszicar = open('OSZICAR','r') 
            laststep = oszicar.readlines()[-2].split()[1] # Next to last line, 2nd entry
            oszicar.close()
            os.chdir(lastfolder) 
            value = int(laststep)
            return value         
        except:
            os.chdir(lastfolder)         
            return 9999    
        
    def getEnergy(self,dir): 
        lines = self.readfile(dir+'/OSZICAR')
        if len(lines[-1].split())>1:
            energy = float(lines[-1].split()[2])  #Free energy
        else: 
            energy = 0.0
        return energy

    def getPureEs(self,iatom):
        lastDir = os.getcwd()
        dir = lastDir + '/' + self.atoms[iatom]
        os.chdir(dir)
        pureHdir =  '/1'
        pureMdir =  '/3'
        
        if os.path.exists(pureHdir):
            self.setAtomCounts(pureHdir)
            self.setEnergy(pureHdir)
            self.pureHenergy = float(self.energy)
        else:
            self.pureHenergy = self.getPureEnergyFromExisting('H', dir+'/enumpast/structures.start') 
        if os.path.exists(pureMdir):        
            self.setAtomCounts(pureMdir)
            self.setEnergy(pureMdir)
            self.pureMenergy = float(self.energy)
        else:
            self.pureMenergy = self.getPureEnergyFromExisting('M', dir+'/enumpast/structures.start') 
#                    if etest != 999999:
#                        self.pureMenergy = etest
#                    else:
# 

    def getPureEnergyFromExisting(self, label, path):
        """ This method is used to extract the energy of the pure structures when they are in the
            structures.start file and hence will not be calculated. If the pure structure is not
            in the structures.start file, return 999999. """
        lines = self.readfile(path)       
        startLooking = False
        if label == 'H':
            for i in xrange(len(lines)):
                lineParts = lines[i].strip().split()
                if lineParts[0].lower() == 'pure' and lineParts[1].lower() == 'h':
                    startLooking = True
                
                if startLooking:
                    if lineParts[0] == '#Energy:':
                        return float(lines[i + 1].strip())
                        
        elif label == 'M':
            for i in xrange(len(lines)):
                lineParts = lines[i].strip().split()
                if lineParts[0].lower() == 'pure' and lineParts[1].lower() == 'm':
                    startLooking = True
                
                if startLooking:
                    if lineParts[0] == '#Energy:':
                        return float(lines[i + 1].strip())
        
        return 999999

    def getNSW(self,dir): 
        proc = subprocess.Popen(['grep','-i','NSW',dir+'/INCAR'],stdout=subprocess.PIPE) 
        return int(proc.communicate()[0].split('=')[-1])   

    def writeHoldoutFromIn(self, atomDir):
        '''Writes structures.holdout for the first iteration, for a given atom
        If starting from existing struct calculations, takes up to N structs in the top of the 
        structures.in file.  Useful when starting from existing structs.  In this case, 
        they should be the lowest N FE structs, since past_structs.dat should be ordered at first.'''            
 
        nmax = 100
        infile = open(atomDir + '/fits/structures.in', 'r')
        holdoutFile = open(atomDir + '/fits/structures.holdout', 'w')

        count = 0
        for line in infile:
            if list(line.strip().split()[0])[:2] == ['#', '-']:
                if count >= nmax:
                    break
                count += 1

            holdoutFile.write(line)

        infile.close()
        holdoutFile.close()

    def getSteps(self, folder):
        """ Returns the number of ionic relaxation steps that VASP performed, as an integer. """
        lastfolder = os.getcwd()
        os.chdir(folder)
        if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
            os.chdir(lastfolder) 
            return -9999
        oszicar = open('OSZICAR', 'r')
        laststep = oszicar.readlines()[-1].split()[0]
        oszicar.close()
        os.chdir(lastfolder)  
        try:
            value = int(laststep)
            return value
        except:
            return 9999 
           
    def hexMonolayerEnergies(self,dir1): 
        file = open(dir1 +'/hex_monolayer_refs/hex_energies','w')
        self.hexE = zeros(len(self.atoms),dtype = float) +100  #default to large number so can tell if not read
        subprocess.call(['echo', '\nReading hexagonal monolayer energies\n'])
        for iatom,atom in enumerate(self.atoms):
            dir2 = dir1 + '/hex_monolayer_refs'+'/'+atom
            if self.FinishCheck(dir2) and self.convergeCheck(dir2, self.getNSW(dir2)): #finaldir
                print'{} monolayer (per atom): {:8.4f} '.format(atom,self.getEnergy(dir2))
                file.write('{} monolayer (per atom): {:8.4f} \n'.format(atom,self.getEnergy(dir2)))
                self.hexE[iatom] = self.getEnergy(dir2) 
            else:
                file.write('{} monolayer not converged \n'.format(atom))
        os.chdir(dir1)  

    def makeUncleFiles(self, iteration, holdoutStructs,vstructsCurrent,vdata):
        """ Runs through the whole process of creating structures.in and structures.holdout files
            for each metal atom. """
        self.vdata = vdata
        self.singleAtomsEnergies(os.getcwd())   
        self.hexMonolayerEnergies(os.getcwd())    
        self.analyzeNewVasp(vstructsCurrent)
        for iatom,atom in enumerate(self.atoms):
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                subprocess.call(['echo', '\nCreating structures.in and structures.holdout files for ' + atom + '\n'])
                self.reinitialize()              
                if iteration == 1 and self.startFromExisting[iatom]: #need only structures.holdout
                    startfile = os.getcwd() + '/needed_files/structures.start.' + atom                     
                    subprocess.call(['cp',startfile,atomDir + '/fits/structures.in'])
                    subprocess.call(['cp',startfile,atomDir + '/enumpast/structures.start'])
                    self.getPureEs(iatom)
                    self.writeHoldoutFromIn(atomDir)                      
                else: #need both structures.in and .holdout                   
                    self.getPureEs(iatom)
                    if len(self.newlyFinished[iatom]) > 0: self.vaspToVdata(iatom)
                    structsInPath = atomDir + '/fits/structures.in'
                    if os.path.exists(structsInPath):
                        outfile = open(structsInPath, 'a')
                    else:
                        outfile = open(structsInPath, 'w')                                   
                        outfile.write(self.header)                                                       
                    for structure in self.newlyFinished[iatom]:
                        self.writePOSCAR(structure,outfile)
                    outfile.close                                    
                    
                    #Fix below:  it's possible that none of the structures in holdoutStructs
                    #have folders and POSCARS in the atom directory, because they came from 
                    #structure.start.  So the poscars need to come not from  writePOSCAR(dir,..), but 
                    #to generate the psuedo-POSCAR, convert it and write it to holdout. 
#                    outfile = open(atomDir + '/fits/structures.holdout', 'w')
#                    outfile.write(self.header) 
#                    for structure in holdoutStructs[iatom]:
#                            # Make sure the structure has converged before trying to write it to
#                            # structures.holdout
#                        if self.contains(structure, self.newlyFinished[iatom]): self.writePOSCAR(structure, self.holdoutFile)
#                    outfile.close 
                    # Replace below when fix above is done: Just using a default holdout
                    subprocess.call(['cp','/needed_files/holdount.in',atomDir + '/fits/'])
                self.vFE2PlotFiles(iatom) #record vasp formation/binding energies and write to files for plotting in gss
        return self.newlyFinished, self.newlyFailed, vdata
                    
    def readfile(self,filepath): 
        file1 = open(filepath,'r')
        lines = file1.readlines()
        file1.close()
        return lines

    def reinitialize(self):
        """ Re-initializes the class for a new metal atom. """
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0 

    def setAtomCounts(self, poscarDir):
        """ Retrieves the number of H atoms and the number of M atoms from the POSCAR file and sets 
            the corresponding members. """
        self.atomCounts = []

        poscar = open(poscarDir + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        counts = poscarLines[5].strip().split()
        
        if len(counts) == 3:
            self.atomCounts.append(int(counts[1]))
            self.atomCounts.append(int(counts[2]))
        elif len(counts) == 2:
            if poscarLines[0].split()[1] == 'H':
                self.atomCounts.append(int(counts[1]))
                self.atomCounts.append(0)
            elif poscarLines[0].split()[1] == 'M':
                self.atomCounts.append(0)
                self.atomCounts.append(int(counts[1]))

    def setAtomPositions(self, poscarDir):
        """ Retrieves the positions of each of the atoms.  Appends the x-coordinate to the 
            xPositions list, the y-coordinate to the yPositions list.  For a surface in UNCLE the 
            z-coordinate is always zero. """
        poscar = open(poscarDir + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomPositions = poscarLines[7 + sum(self.atomCounts):7 + (2 * sum(self.atomCounts))]
        self.atomPositions = [line.strip().split() for line in self.atomPositions]
           
        for pos in self.atomPositions:
            self.xPositions.append(float(pos[0]))
            self.yPositions.append(float(pos[1]))
            self.zPositions.append(0.0)

    def setEnergy(self, directory):  
        """ Retrieves the energy of the structure from the OSZICAR file and sets the corresponding 
            member. """   
        try:
            oszicar = open(directory + self.finalDir + '/OSZICAR','r')
            energy = oszicar.readlines()[-1].split()[2]
            oszicar.close()
        except:
            energy = 0
        
        energy = float(energy)
        peratom = energy / sum(self.atomCounts)       
        self.energy = str(peratom)
        
    def setIDString(self, poscarDir):
        """ Sets the first written line of each structure to the form:
                C H (Metal Atom)  Structure:  #(Decimal) (#(Binary)) """
        print 'poscarDir', poscarDir
        poscar = open(poscarDir + '/POSCAR', 'r')
        ID = poscar.readlines()[0].strip()       
        self.idString = ID
        
    def setLatticeVectors(self, structFile):
        """ Gets the lattice vectors from the first structure in the newlyFinished and sets the 
            corresponding member components. """
        vecFile = open(structFile + '/POSCAR', 'r')
        vecFileLines = vecFile.readlines()
        vecFile.close()
        
        vec1 = vecFileLines[2].strip().split()
        vec2 = vecFileLines[3].strip().split()
        vec3 = vecFileLines[4].strip().split()
        
        vec1comps = [float(comp) for comp in vec1]
        vec2comps = [float(comp) for comp in vec2]
        vec3comps = [float(comp) for comp in vec3]
            
        
        self.vec1x = vec1comps[0]
        self.vec1y = vec1comps[1]
        if vec1comps[2] == 15.0:
            self.vec1z = 1000.0
        else:
            self.vec1z = vec1comps[2]
        
        self.vec2x = vec2comps[0]
        self.vec2y = vec2comps[1]
        if vec2comps[2] == 15.0:
            self.vec2z = 1000.0
        else:
            self.vec2z = vec2comps[2]
        
        self.vec3x = vec3comps[0]
        self.vec3y = vec3comps[1]
        if vec3comps[2] == 15.0:
            self.vec3z = 1000.0
        else:
            self.vec3z = vec3comps[2]

    def singleAtomsEnergies(self,dir1): 
        self.singleE = zeros(len(self.atoms),dtype = float) +100  #default to large number so can tell if not read
        subprocess.call(['echo', '\nReading single atom energies\n'])
        file = open(dir1 +'/single_atoms/single_atom_energies','w')
        for iatom,atom in enumerate(self.atoms):
            dir2 = dir1 + '/single_atoms'+'/'+atom
            if self.electronicConvergeFinish(dir2): 
                print 'Energy of {} atom: {:8.4f} \n'.format(atom,self.getEnergy(dir2))
                file.write('{} atom: {:12.8f} \n'.format(atom,self.getEnergy(dir2)))
                self.singleE[iatom] = self.getEnergy(dir2)
        file.close()  
        os.chdir(dir1) 

    def vFE2PlotFiles(self, iatom):
        """ For all finished structs, record the different vasp formation and binding energies to files 
        for plots.  vaspToVdata should be run first"""
        
        # TODO:  We should probably figure out how to sort the structures in existing 
        #        structures.start files together with the structures we have calculated in VASP 
        #        during the loop.
#        lastDir = os.getcwd()
#        os.chdir(lastDir + '/' + self.atoms[iatom])
        eIsolatedH = -1.115 
        eIsolatedC = -1.3179 
        eH2 = -6.7591696/2.0 
        energyGraphene = -18.456521 #for 2 C atoms 
        vaspFEfile = open('vaspFE.out','w')  
        vaspBEfile = open('vaspBE.out','w')  
        vaspHFEfile = open('vaspHFE.out','w')  
        nfinished = count_nonzero(self.vdata[iatom,:]['struct'])
        istruct = 0 #creating for all finished structs
        for i in range(nfinished):
            struct = self.vdata[iatom,istruct]['struct']
            conc = self.vdata[iatom,istruct]['conc']
            natoms = self.vdata[iatom,istruct]['natoms']
            nmetal = int(conc*natoms)
            nH = natoms - nmetal 
            ncarbon = nH + nmetal #all C sites have an adatom           
            structEnergy = self.vdata[iatom,istruct]['energy']                 
            formationEnergy = structEnergy - (conc * self.pureMenergy + (1.0 - conc) * self.pureHenergy)
            self.vdata[iatom,istruct]['FE'] = formationEnergy
            vaspFEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,formationEnergy))            
            bindEnergy = structEnergy - (nH*eIsolatedH + nmetal*self.singleE[iatom] + ncarbon*energyGraphene/2.0)/ float(natoms) #2 atoms in graphene 
            self.vdata[iatom,istruct]['BE'] = bindEnergy
            vaspBEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,bindEnergy))  
            hexFormationEnergy = structEnergy - energyGraphene/2  - (conc * self.hexE[iatom] + (1.0 - conc) * eH2)
            self.vdata[iatom,istruct]['HFE'] = hexFormationEnergy
            vaspHFEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,hexFormationEnergy))    
            istruct += 1                          
        vaspFEfile.close()
        vaspBEfile.close()
        vaspHFEfile.close()                
#        os.chdir(lastDir)
    
    def vaspToVdata(self, iatom):
        """ Record the newly finished structures into vdata, and 
        leaves the newly finished structures sorted by formation energy"""
        # TODO:  We should probably figure out how to sort the structures in existing 
        #        structures.start files together with the structures we have calculated in VASP 
        #        during the loop.
        lastDir = os.getcwd()
        os.chdir(lastDir + '/' + self.atoms[iatom])
#        pureHdir = os.getcwd() + '/1'
#        pureMdir = os.getcwd() + '/3'
#        
#        self.setAtomCounts(pureHdir)
#        self.setEnergy(pureHdir)
#        self.pureHenergy = float(self.energy)
#        
#        self.setAtomCounts(pureMdir)
#        self.setEnergy(pureMdir)
#        self.pureMenergy = float(self.energy)
        
        formEnergyList = []
        sortedStructs = []
        nfinished = count_nonzero(vdata[iatom,:]['struct'])
        istruct = nfinished #starting position for vdata array
        for struct in self.newlyFinished[iatom]:
            self.setAtomCounts(struct) #reads in atom numbers
            self.setEnergy(struct)
            self.vdata[iatom,istruct]['struct'] = struct
            structEnergy = float(self.energy)
            self.vdata[iatom,istruct]['energy'] = structEnergy 
            conc = 0.0
            natoms =  float(self.atomCounts[0] + self.atomCounts[1])
            self.vdata[iatom,istruct]['natoms'] = natoms
            if self.atomCounts[0] == 0:
                conc = 1.0
            else:
                conc = float(float(self.atomCounts[1])/natoms)
            self.vdata[iatom,istruct]['conc'] = conc                       
            formationEnergy = structEnergy - (conc * self.pureMenergy + (1.0 - conc) * self.pureHenergy)
            formEnergyList.append([formationEnergy, structDir])
            istruct += 1                                
        formEnergyList.sort()       
        for pair in formEnergyList:
            sortedStructs.append(pair[1])        
        self.newlyFinished[iatom] = sortedStructs           
        os.chdir(lastDir)
        
    def writeAtomCounts(self):
        """ Writes the number of H atoms and the number of M atoms to the structures.in or 
            structures.holdout file. """
        if len(self.atomCounts) == 2:
            self.outfile.write(str(self.atomCounts[0]) + " " + str(self.atomCounts[1]) + "\n")
        elif len(self.atomCounts) == 1:
            self.outfile.write(str(self.atomCounts[0]) + "\n")
    
    def writeAtomPositions(self):
        """ Writes the positions of the atoms in the current structure to the structures.in
            or structures.holdout file.  The positions are written in the form:
                z-coord   x-coord   y-coord
            because this is the convention that UNCLE uses. """
        self.outfile.write("Cartesian\n")
        for i in xrange(len(self.atomPositions)):
            self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % 
                               (self.zPositions[i], self.xPositions[i], self.yPositions[i]))   

    def writeDashedLine(self):
        """ Writes a dashed line in the structures.in/.holdout files as a separator between
            different structures. """
        self.outfile.write("#------------------------------------------------\n")

    def writeEnergy(self):
        """ Writes the energy of the current structure to the structures.in or structures.holdout
            file. """
        self.outfile.write("#Energy:\n")
        self.outfile.write(str(self.energy) + "\n")    

#    def writeHeader(self,index):
#        """ Writes the headers of the structures.in and structures.holdout files. """
#        if not self.startFromExisting[index]:
#            self.infile.write(self.header)
#        self.holdoutFile.write(self.header)
       
    def writeIDString(self):
        """ Writes the ID string of the current structure to either the structures.in or 
            structures.holdout file. """
        concentration = float(float(self.atomCounts[1]) / float(sum(self.atomCounts)))
        formationEnergy = float(self.energy) - (concentration * self.pureMenergy + (1.0 - concentration) * self.pureHenergy)
        
        self.outfile.write(self.idString + " FE = " + str(formationEnergy) + ", Concentration = " + str(concentration) + "\n")
        
    def writeLatticeVecs(self):
        """ Writes the lattice vectors of the current structure to the structures.in or
            structures.holdout file. """
        self.outfile.write("1.0\n")
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec1z, self.vec1x, self.vec1y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec2z, self.vec2x, self.vec2y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec3z, self.vec3x, self.vec3y))
    
    def writePOSCAR(self, poscarDir, outfile):
        """ Calls all the methods needed to write all the needed information about the current
            structure to the 'structures.in' or 'structures.holdout' files.  Uses an input list
            to decide which structures to put in the 'structures.holdout' file. """
        self.outfile = outfile
        self.setIDString(poscarDir)
        self.setLatticeVectors(poscarDir)
        self.setAtomCounts(poscarDir)
        self.setAtomPositions(poscarDir)
        self.setEnergy(poscarDir)
        
        # Make sure the pure structures go in structures.in
        if self.idString.split()[0] == 'PURE':
            self.outfile = outfile
        
        self.writeDashedLine()
        self.writeIDString()
        self.writeLatticeVecs()
        self.writeAtomCounts()
        self.writeAtomPositions()
        self.writeEnergy()









