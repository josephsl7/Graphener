'''
Created on Aug 29, 2014

@author: eswens13
'''
from numpy import zeros, mod, count_nonzero,sort
import os, subprocess, sys
from random import random
from Main import slurmProblem
from comMethods import *
 
class MakeUncleFiles:
    from comMethods import setAtomCounts
    def __init__(self, atoms, startMethod, iteration,finalDir,restartTimeout):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.startMethod = startMethod
        self.iteration = iteration
        self.finalDir = finalDir   
        self.restartTimeout = restartTimeout
        
        self.newlyFinished = [] #list of structures, not paths
        self.newlyFailed = []
        self.newlyToRestart = []
        self.pureHenergy = 0.0
        self.pureMenergy = 0.0
        
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        #self.header = "peratom\nnoweights\nposcar\n"
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

    def analyzeNewVasp(self,vstructsToRun):
        """ Initializes the list of structures to add to the structures.in and structures.holdout
            files by adding only the structures that VASP finished to the member
            newlyFinished. Sorts the list by metal concentration (N_M / N_total). Adds the structures that
            failed VASP calculations to the member 'failedList'. """   
        lastDir = os.getcwd()
        self.newlyFinished = [[]]*len(self.atoms)
        self.newlyFailed = [[]]*len(self.atoms)
        self.newlyToRestart = [[]]*len(self.atoms)
        
        for iatom, atom in enumerate(self.atoms):
            # If it is the first iteration and we are starting from an existing structures.start
            # file, we just append an empty list to the newlyFinished and the failedList.  Else, 
            # proceed as normal.
            if len(vstructsToRun[iatom])>0:
                atomDir = lastDir + '/' + atom
                os.chdir(atomDir)        
                conclist = []
                finished = []
                failed = []
                restart = []

                for istruct, struct in enumerate(vstructsToRun[iatom]):
                    if mod(istruct+1,100) == 0 or istruct+1 == len(vstructsToRun[iatom]): 
                        subprocess.call(['echo','\tChecked {} of {} structures in {}'.format(istruct+1,len(vstructsToRun[iatom]),atom)])
                    if os.path.isdir(atomDir + '/' + struct):
                        vaspDir = atomDir + '/' + struct + self.finalDir
                        if os.path.isdir(vaspDir):
                            if finishCheck(vaspDir) and convergeCheck(vaspDir, getNSW(vaspDir)) and energyDropCheck(vaspDir): #finalDir                           
                               # Check for concentration
                                self.setAtomCounts(struct)                            
                                concentration = 0.0
                                if self.atomCounts[0] == 0:
                                    concentration = 1.0
                                else:
                                    concentration = float(float(self.atomCounts[1]) / float(self.atomCounts[0] + self.atomCounts[1]))                           
                                conclist.append([concentration, struct])
                            elif self.restartTimeout and not slurmProblem(vaspDir):
                                restart.append(struct)
                            else:
                                failed.append(struct)
                        else:
                            subprocess.call(['echo', '\nERROR: directory does not exist: ' + struct]) 
                self.newlyFailed[iatom] = failed #for atom i                
                conclist.sort()                
                for i in xrange(len(conclist)): finished.append(conclist[i][1]) 
                self.newlyFinished[iatom]= finished
                self.newlyToRestart[iatom]= restart #sorted by concentration, for atomi
                subprocess.call(['echo','\tAtom {} vasp calcs: {} finished, {} to restart, {} failed\n'.format(atom,len(self.newlyFinished[iatom]),len(self.newlyToRestart[iatom]),len(self.newlyFailed[iatom]))])

                os.chdir(lastDir)

    def contains(self, struct, alist):
        """ Returns True if the list 'alist' contains the structure 'struct', False otherwise. """
        if len(alist) <= 0:
            return False

        for elem in alist:
            if str(struct) == str(elem):
                return True

        return False

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
        return self.elConvergeCheck(dir,NELM) and finishCheck(dir)
            
    def energyDropCheck(self,dir):
        '''tests whether the energies have dropped in OSZICAR...rising energies are unphysical 
        and show a numerical convergence problem. The factor like 0.99 allows for very small energy rises only'''
        lines = readfile(dir+'/OSZICAR') 
        energies = []
        for line in lines:
            if 'F=' in line:
                energies.append(float(line.split()[2]))
        return energies[-1] <= 0.99*energies[0] 
    

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
        lines = readfile(dir+'/OSZICAR')
        if len(lines[-1].split())>1:
            energy = float(lines[-1].split()[2])  #Free energy
        else: 
            energy = 0.0
        return energy

    def getPureEs(self,iatom):
        lastDir = os.getcwd()
        dir = lastDir + '/' + self.atoms[iatom]
        os.chdir(dir)
        pureHdir =  dir + '/1'
        pureMdir =  dir + '/3'
        
        if os.path.exists(pureHdir):
            self.setAtomCounts(pureHdir)
            self.setEnergy(pureHdir)
            self.pureHenergy = float(self.energy)
        else:
            subprocess.call(['echo','Missing pure H energy folder'])
        if os.path.exists(pureMdir):        
            self.setAtomCounts(pureMdir)
            self.setEnergy(pureMdir)
            self.pureMenergy = float(self.energy)
        else:
            subprocess.call(['echo','Missing pure M energy folder']) 
#                    if etest != 999999:
#                        self.pureMenergy = etest
#                    else:
        os.chdir(lastDir)
        
#    def getPureEnergyExisting(self, label, path):
#        """ This method is used to extract the energy of the pure structures when they are in the
#            structures.start file and hence will not be calculated. If the pure structure is not
#            in the structures.start file, return 999999. """
#        lines = readfile(path)       
#        startLooking = False
#        if label == 'H':
#            for i in xrange(len(lines)):
#                lineParts = lines[i].strip().split()
#                if lineParts[0].lower() == 'pure' and lineParts[1].lower() == 'h':
#                    startLooking = True
#                
#                if startLooking:
#                    if lineParts[0] == '#Energy:':
#                        return float(lines[i + 1].strip())
#                        
#        elif label == 'M':
#            for i in xrange(len(lines)):
#                lineParts = lines[i].strip().split()
#                if lineParts[0].lower() == 'pure' and lineParts[1].lower() == 'm':
#                    startLooking = True
#                
#                if startLooking:
#                    if lineParts[0] == '#Energy:':
#                        return float(lines[i + 1].strip())
#        
#        return 999999

#    def getNSW(self,dir): 
#        proc = subprocess.Popen(['grep','-i','NSW',dir+'/INCAR'],stdout=subprocess.PIPE) 
#        return int(proc.communicate()[0].split('=')[-1])   

#    def writeHoldoutFromIn(self, atomDir):
#        '''Writes structures.holdout for the first iteration, for a given atom
#        If starting from existing struct calculations, takes up to N structs in the top of the 
#        structures.in file.  Useful when starting from existing structs.  In this case, 
#        they should be the lowest N FE structs, since past_structs.dat should be ordered at first.'''            
# 
#        nmax = 100
#        infile = open(atomDir + '/structures.in', 'r')
#        holdoutFile = open(atomDir + '/fits/structures.holdout', 'w')
#
#        count = 0
#        for line in infile:
#            if list(line.strip().split()[0])[:2] == ['#', '-']:
#                if count >= nmax:
#                    break
#                count += 1
#
#            holdoutFile.write(line)
#
#        infile.close()
#        holdoutFile.close()
           
    def hexMonolayerEnergies(self,dir1,iteration): 
        file = open(dir1 +'/hex_monolayer_refs/hex_energies','w')
        self.hexE = zeros(len(self.atoms),dtype = float) +100  #default to large number so can tell if not read
        if iteration == 1: subprocess.call(['echo', '\nReading hexagonal monolayer energies\n'])
        for iatom,atom in enumerate(self.atoms):
            dir2 = dir1 + '/hex_monolayer_refs'+'/'+atom
            if finishCheck(dir2) and convergeCheck(dir2, getNSW(dir2)) and energyDropCheck(vaspDir): #finalDir
                if iteration == 1: subprocess.call(['echo','{} monolayer (per atom): {:8.4f} '.format(atom,self.getEnergy(dir2))])
                file.write('{} monolayer (per atom): {:8.4f} \n'.format(atom,self.getEnergy(dir2)))
                self.hexE[iatom] = self.getEnergy(dir2) 
            else:
                file.write('{} monolayer not converged \n'.format(atom))
        os.chdir(dir1)  

    def makeUncleFiles(self, iteration, holdoutStructs,vstructsToRun,vdata):
        """ Runs through the whole process of creating structures.in and structures.holdout files
            for each metal atom. """

        self.vdata = vdata
        self.singleAtomsEnergies(os.getcwd(),iteration)   
        self.hexMonolayerEnergies(os.getcwd(),iteration)    
        self.analyzeNewVasp(vstructsToRun)

        for iatom,atom in enumerate(self.atoms):
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                subprocess.call(['echo', '\nCreating structures.in and structures.holdout files for ' + atom + '\n'])
                self.reinitialize()                             
                self.getPureEs(iatom)
                nNewFinished = len(self.newlyFinished[iatom])
                if nNewFinished > 0: 
                    newPos = self.vaspToVdata(iatom) #starting position of new data
                    i1 = newPos; i2 = newPos + nNewFinished
                    structuresInWrite(atomDir,vdata[iatom,i1:i2]['struct'], vdata[iatom,i1:i2]['FE'],\
                         vdata[iatom,i1:i2]['conc'],vdata[iatom,i1:i2]['energy'],'a')                                  
                    print "vdata[iatom,i1:i2]['struct']",vdata[iatom,i1:i2]['struct']
                    print "vdata[iatom,i1:i2]['FE']", vdata[iatom,i1:i2]['FE']
                    print "vdata[iatom,i1:i2]['conc']", vdata[iatom,i1:i2]['conc']
                    print"vdata[iatom,i1:i2]['energy']"
                
                    #Fix below:  it's possible that none of the structures in holdoutStructs
                    #have folders and POSCARS in the atom directory, because they came from 
                    #structure.start.  So the poscars need to come not from  writeUnclePOSCAR(dir,..), but 
                    #to generate the pseudo-POSCAR, convert it and write it to holdout. 
#                    outfile = open(atomDir + '/fits/structures.holdout', 'w')
#                    outfile.write(self.header) 
#                    for structure in holdoutStructs[iatom]:
#                            # Make sure the structure has converged before trying to write it to
#                            # structures.holdout
#                        if self.contains(structure, self.newlyFinished[iatom]): self.writeUnclePOSCAR(structure, self.holdoutFile)
#                    outfile.close 
                    # Replace below when fix above is done: Just using a default holdout
                self.vdataToPlotFiles(iatom) #record vasp formation/binding energies and write to files for plotting in gss
            subprocess.call(['cp','needed_files/structures.holdout',atomDir + '/fits/'])
    
        return self.newlyFinished, self.newlyToRestart, self.newlyFailed, self.vdata

    def reinitialize(self):
        """ Re-initializes the class for a new metal atom. """
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        #self.header = "peratom\nnoweights\nposcar\n"
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

#    def setAtomCounts(self, poscarDir):
#        """ Retrieves the number of H atoms and the number of M atoms from the POSCAR file and sets 
#            the corresponding members. """
#        self.atomCounts = []
#
#        poscar = open(poscarDir + '/POSCAR', 'r')
#        poscarLines = poscar.readlines()
#        poscar.close()
#        
#        counts = poscarLines[5].strip().split()
#        
#        if len(counts) == 3:
#            self.atomCounts.append(int(counts[1]))
#            self.atomCounts.append(int(counts[2]))
#        elif len(counts) == 2:
#            if poscarLines[0].split()[1] == 'H':
#                self.atomCounts.append(int(counts[1]))
#                self.atomCounts.append(0)
#            elif poscarLines[0].split()[1] == 'M':
#                self.atomCounts.append(0)
#                self.atomCounts.append(int(counts[1]))

#    def setAtomPositions(self, poscarDir):
#        """ Retrieves the positions of each of the atoms.  Appends the x-coordinate to the 
#            xPositions list, the y-coordinate to the yPositions list.  For a surface in UNCLE the 
#            z-coordinate is always zero. """
#        poscar = open(poscarDir + '/POSCAR', 'r')
#        poscarLines = poscar.readlines()
#        poscar.close()
#        
#        self.atomPositions = []
#        self.xPositions = []
#        self.yPositions = []
#        self.zPositions = []
#        
#        self.atomPositions = poscarLines[7 + sum(self.atomCounts):7 + (2 * sum(self.atomCounts))]
#        self.atomPositions = [line.strip().split() for line in self.atomPositions]
#           
#        for pos in self.atomPositions:
#            self.xPositions.append(float(pos[0]))
#            self.yPositions.append(float(pos[1]))
#            self.zPositions.append(0.0)


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
        poscar = open(poscarDir + '/POSCAR', 'r')
        ID = poscar.readlines()[0].strip()       
        self.idString = ID
        
#    def setLatticeVectors(self, structDir):
#        """ Gets the lattice vectors from the first structure in the newlyFinished and sets the 
#            corresponding member components. """
#        vecFile = open(structDir + '/POSCAR', 'r')
#        vecFileLines = vecFile.readlines()
#        vecFile.close()
#        
#        vec1 = vecFileLines[2].strip().split()
#        vec2 = vecFileLines[3].strip().split()
#        vec3 = vecFileLines[4].strip().split()
#        
#        vec1comps = [float(comp) for comp in vec1]
#        vec2comps = [float(comp) for comp in vec2]
#        vec3comps = [float(comp) for comp in vec3]
#            
#        
#        self.vec1x = vec1comps[0]
#        self.vec1y = vec1comps[1]
#        if vec1comps[2] == 15.0:
#            self.vec1z = 1000.0
#        else:
#            self.vec1z = vec1comps[2]
#        
#        self.vec2x = vec2comps[0]
#        self.vec2y = vec2comps[1]
#        if vec2comps[2] == 15.0:
#            self.vec2z = 1000.0
#        else:
#            self.vec2z = vec2comps[2]
#        
#        self.vec3x = vec3comps[0]
#        self.vec3y = vec3comps[1]
#        if vec3comps[2] == 15.0:
#            self.vec3z = 1000.0
#        else:
#            self.vec3z = vec3comps[2]

    def singleAtomsEnergies(self,dir1,iteration): 
        self.singleE = zeros(len(self.atoms),dtype = float) +100  #default to large number so can tell if not read
        if iteration == 1: subprocess.call(['echo', '\nReading single atom energies\n'])
        file = open(dir1 +'/single_atoms/single_atom_energies','w')
        for iatom,atom in enumerate(self.atoms):
            dir2 = dir1 + '/single_atoms'+'/'+atom
            if self.electronicConvergeFinish(dir2): 
                if iteration == 1: subprocess.call(['echo', 'Energy of {} atom: {:8.4f} \n'.format(atom,self.getEnergy(dir2))])
                file.write('{} atom: {:12.8f} \n'.format(atom,self.getEnergy(dir2)))
                self.singleE[iatom] = self.getEnergy(dir2)
        file.close()  
        os.chdir(dir1) 

#    def timeoutCheck(self,structDir):
#        '''checks the latest slurm file for the words TIME LIMIT'''
#        slurmlist = []               
#        structfiles = os.listdir(structDir)
#        for file in structfiles:
#            filePath = structDir + '/' + file
#            if 'slurm' in file and os.stat(filePath).st_size > 0:
#                slurmlist.append(file)
#        if len(slurmlist)>0: #only use the last slurm
#            slurmlist.sort()
#            lastslurm = slurmlist[-1]
#            for line in readfile(structDir + '/' + lastslurm):
#                if 'TIME LIMIT' in line: 
#                    return True
#        return False        

    def vdataToPlotFiles(self, iatom):
        """ For all finished structs, record the different vasp formation and binding energies to files 
        for plots.  vaspToVdata should be run first"""
        lastDir = os.getcwd()
        os.chdir(lastDir + '/' + self.atoms[iatom])
        eIsolatedH = -1.1073   
        eIsolatedC = -1.3179 
        eH2 = -6.7591696/2.0 
        energyGraphene = -18.456521 #for 2 C atoms 
        vaspFEfile = open('vaspFE.out','w')  
        vaspBEfile = open('vaspBE.out','w')  
        vaspHFEfile = open('vaspHFE.out','w')  
        nfinished = count_nonzero(self.vdata[iatom,:]['struct'])
        istruct = 0 #creating for all finished structs
        #this does not need to be sorted. 
        for j in range(nfinished):
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
            #note that the hex monolayers have one atom per cell
            hexFormationEnergy = structEnergy - energyGraphene/2  - (conc * self.hexE[iatom] + (1.0 - conc) * eH2)
            self.vdata[iatom,istruct]['HFE'] = hexFormationEnergy
            vaspHFEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,hexFormationEnergy))    
            istruct += 1                          
        vaspFEfile.close()
        vaspBEfile.close()
        vaspHFEfile.close()                
        os.chdir(lastDir)
    
    def vaspToVdata(self, iatom):
        """ Record the newly finished structures into vdata, and 
        leaves the newly finished structures sorted by formation energy"""
        # TODO:  We should probably figure out how to sort the structures in existing 
        #        structures.start files together with the structures we have calculated in VASP 
        #        during the loop.
        lastDir = os.getcwd()
        os.chdir(lastDir + '/' + self.atoms[iatom])
        formEnergyList = []
        sortedStructs = []
        nfinished = count_nonzero(self.vdata[iatom,:]['struct'])
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
            formEnergyList.append([formationEnergy, struct])
            istruct += 1                                
        formEnergyList.sort()       
        for pair in formEnergyList:
            sortedStructs.append(pair[1])        
        self.newlyFinished[iatom] = sortedStructs           
        os.chdir(lastDir)
        return nfinished #the position in vdata before adding the data
            
#    def writeAtomCounts(self):
#        """ Writes the number of H atoms and the number of M atoms to the structures.in or 
#            structures.holdout file. """
#        if len(self.atomCounts) == 2:
#            self.outfile.write(str(self.atomCounts[0]) + " " + str(self.atomCounts[1]) + "\n")
#        elif len(self.atomCounts) == 1:
#            self.outfile.write(str(self.atomCounts[0]) + "\n")
#    
#    def writeAtomPositions(self):
#        """ Writes the positions of the atoms in the current structure to the structures.in
#            or structures.holdout file.  The positions are written in the form:
#                z-coord   x-coord   y-coord
#            because this is the convention that UNCLE uses. """
#        self.outfile.write("Cartesian\n")
#        for i in xrange(len(self.atomPositions)):
#            self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % 
#                               (self.zPositions[i], self.xPositions[i], self.yPositions[i]))   
#
#    def writeDashedLine(self):
#        """ Writes a dashed line in the structures.in/.holdout files as a separator between
#            different structures. """
#        self.outfile.write("#------------------------------------------------\n")
#
#    def writeEnergy(self):
#        """ Writes the energy of the current structure to the structures.in or structures.holdout
#            file. """
#        self.outfile.write("#Energy:\n")
#        self.outfile.write(str(self.energy) + "\n")    
#      
#    def writeIDString(self):
#        """ Writes the ID string of the current structure to either the structures.in or 
#            structures.holdout file. """
#        concentration = float(float(self.atomCounts[1]) / float(sum(self.atomCounts)))
#        formationEnergy = float(self.energy) - (concentration * self.pureMenergy + (1.0 - concentration) * self.pureHenergy)
#        
#        self.outfile.write(self.idString + " FE = " + str(formationEnergy) + ", Concentration = " + str(concentration) + "\n")
        
#    def writeLatticeVecs(self):
#        """ Writes the lattice vectors of the current structure to the structures.in or
#            structures.holdout file. """
#        self.outfile.write("1.0\n")
#        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec1z, self.vec1x, self.vec1y))
#        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec2z, self.vec2x, self.vec2y))
#        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec3z, self.vec3x, self.vec3y))
#    
#    def writeUnclePOSCAR(self, poscarDir, outfile, idString):
#        """ Calls all the methods needed to write all the needed information about the current
#            structure to the 'structures.in' or 'structures.holdout' files.  Uses an input list
#            to decide which structures to put in the 'structures.holdout' file. """
#        self.outfile = outfile
#        if idString == '':
#            self.setIDString(poscarDir)
#        else: 
#            self.setIDString = idString
#        self.setLatticeVectors(poscarDir)
#        self.setAtomCounts(poscarDir)
#        self.setAtomPositions(poscarDir)
#        self.setEnergy(poscarDir)
#        
#        # Make sure the pure structures go in structures.in
#        if self.idString.split()[0] == 'PURE':
#            self.outfile = outfile       
#        self.writeDashedLine(); outfile.flush()       
#        self.writeIDString(); outfile.flush()        
#        self.writeLatticeVecs(); outfile.flush()     
#        self.writeAtomCounts(); outfile.flush()     
#        self.writeAtomPositions(); outfile.flush()     
#        self.writeEnergy(); outfile.flush()     





