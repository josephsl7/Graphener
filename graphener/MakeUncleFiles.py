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
    def __init__(self, atoms, startMethod, iteration,finalDir,restartTimeout,pureMetal):
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
        self.pureMetal = pureMetal
        
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
                                nadatoms =  float(self.atomCounts[1] + self.atomCounts[2])                         
                                concentration = 0.0
                                if self.atomCounts[1] == 0:
                                    concentration = 1.0
                                else:
                                    concentration = float(float(self.atomCounts[2]) / float(nadatoms))                           
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
        pureMdir =  dir + '/' + self.pureMetal
        
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
        


    def hexMonolayerEnergies(self,dir1,iteration): 
        file = open(dir1 +'/hex_monolayer_refs/hex_energies','w')
        self.hexE = zeros(len(self.atoms),dtype = float) +100  #default to large number so can tell if not read
        if iteration == 1: subprocess.call(['echo', '\nReading hexagonal monolayer energies\n'])
        for iatom,atom in enumerate(self.atoms):
            dir2 = dir1 + '/hex_monolayer_refs'+'/'+atom
            if finishCheck(dir2) and convergeCheck(dir2, getNSW(dir2)) and energyDropCheck(dir2): #finalDir
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
                    structuresWrite('all',atomDir,vdata[iatom,i1:i2]['struct'], vdata[iatom,i1:i2]['FE'],\
                                    vdata[iatom,i1:i2]['conc'],vdata[iatom,i1:i2]['energy'],'.in','a')                                                 
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
        peratom = energy / sum(self.atomCounts[1:])       
        self.energy = str(peratom)
        
    def setIDString(self, poscarDir):
        """ Sets the first written line of each structure to the form:
                C H (Metal Atom)  Structure:  #(Decimal) (#(Binary)) """
        poscar = open(poscarDir + '/POSCAR', 'r')
        ID = poscar.readlines()[0].strip()       
        self.idString = ID
        
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
            nadatoms = self.vdata[iatom,istruct]['nadatoms']
            nmetal = int(conc*nadatoms)
            nH = nadatoms - nmetal 
            ncarbon = self.vdata[iatom,istruct]['nCarbon'] 
            #multiply stored energy by nadatoms so we have vasp run energy
            structEnergy = nadatoms * self.vdata[iatom,istruct]['energy'] 
            formationEnergy = (structEnergy - nmetal*self.pureMenergy - nH*self.pureHenergy)/float(nadatoms)
            self.vdata[iatom,istruct]['FE'] = formationEnergy
            vaspFEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,formationEnergy))            
            bindEnergy = (structEnergy - ncarbon*energyGraphene/2.0 - nH*eIsolatedH - nmetal*self.singleE[iatom])/ float(nadatoms) #2 atoms in graphene 
            self.vdata[iatom,istruct]['BE'] = bindEnergy
            vaspBEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,bindEnergy))  
            #note that the hex monolayers have one atom per cell
            hexFormationEnergy = (structEnergy - energyGraphene*ncarbon/2.0  - nmetal *self.hexE[iatom] - nH*eH2)/float(nadatoms)
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
            nCarbon =  atomCounts[0]
            nadatoms =  float(self.atomCounts[1] + self.atomCounts[2])
            self.vdata[iatom,istruct]['nadatoms'] = nadatoms
            if self.atomCounts[1] == 0:
                conc = 1.0
            else:
                conc = float(float(self.atomCounts[2])/nadatoms)
            self.vdata[iatom,istruct]['conc'] = conc                       
            formationEnergy = structEnergy - (conc * self.pureMenergy + (1.0 - conc) * self.pureHenergy)
            self.vdata[iatom,istruct]['FE'] = formationEnergy
            formEnergyList.append([formationEnergy, struct])
            istruct += 1                                
        formEnergyList.sort()       
        for pair in formEnergyList:
            sortedStructs.append(pair[1])        
        self.newlyFinished[iatom] = sortedStructs           
        os.chdir(lastDir)
        return nfinished #the position in vdata before adding the data
