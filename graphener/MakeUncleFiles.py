'''
Created on Aug 29, 2014


'''
from numpy import zeros, mod, count_nonzero,sort
import os, subprocess, sys
from random import random
from Main import slurmProblem
from comMethods import *
 
class MakeUncleFiles:
    from comMethods import setAtomCounts,hexMonolayerEnergies,singleAtomsEnergies,getPureEs,\
                    contains,elConvergeCheck,electronicConvergeFinish,getElSteps,getEnergy,setEnergy
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

        self.case = len(atoms[0].split(','))

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
                    vaspDir = atomDir + '/' + struct + self.finalDir
                    if outcarWarn(vaspDir): 
                        subprocess.call(['echo','\tOUTCAR warning for struct{}: failed'.format(struct)])
                        failed.append(struct)
                    elif not energyDropCheck(vaspDir): 
                        subprocess.call(['echo','\tEnergy rose unphysically for struct {}: failed'.format(struct)])
                        failed.append(struct)                                            
                    else:
                        if finishCheck(vaspDir) and convergeCheck(vaspDir, getNSW(vaspDir)):                           
                           # Check for concentration
                            self.setAtomCounts(struct) 
                            nadatoms =  float(sum(self.atomCounts) - self.atomCounts[0])                       

                            concentration = [float(float(x) / nadatoms) for x in self.atomCounts[1:]]
                          
                            conclist.append([concentration[0], struct])
                        elif self.restartTimeout and not slurmProblem(vaspDir):
                            restart.append(struct)
                        else:
                            failed.append(struct)
                            subprocess.call(['echo','\tStruct {}: failed'.format(struct)])
                self.newlyFailed[iatom] = failed #for atom i                
                conclist.sort()                
                for i in xrange(len(conclist)): finished.append(conclist[i][1]) 
                self.newlyFinished[iatom]= finished
                self.newlyToRestart[iatom]= restart #sorted by concentration, for atomi
                subprocess.call(['echo','\tAtom {} vasp calcs: {} finished, {} to restart, {} failed\n'.format(atom,len(self.newlyFinished[iatom]),len(self.newlyToRestart[iatom]),len(self.newlyFailed[iatom]))])
                os.chdir(lastDir)
        
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
        
    def setIDString(self, poscarDir):
        """ Sets the first written line of each structure to the form:
                C H (Metal Atom)  Structure:  #(Decimal) (#(Binary)) """
        poscar = open(poscarDir + '/POSCAR', 'r')
        ID = poscar.readlines()[0].strip()       
        self.idString = ID
        
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
        elements = self.atoms[iatom].split(',')
        #this does not need to be sorted.
        for j in range(nfinished):
            struct = self.vdata[iatom,istruct]['struct']
            conc = self.vdata[iatom,istruct]['conc']
            nadatoms = self.vdata[iatom,istruct]['nadatoms']
            nmetal = int(conc*nadatoms)
            nH = nadatoms - nmetal 
            ncarbon = self.vdata[iatom,istruct]['nCarbon'] 
            #multiply stored energy by nadatoms so we have vasp run energy
            structEnergy = self.vdata[iatom,istruct]['energy'] 
            self.setAtomCounts(str(struct))

            formationEnergy = structEnergy
            for i, count in enumerate(self.atomCounts[1:]):
                formationEnergy -= float(count) / float(nadatoms) * self.pureEnergies[i]

            self.vdata[iatom,istruct]['FE'] = formationEnergy
            vaspFEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,formationEnergy))            

            bindEnergy = structEnergy - energyGraphene / 2.0 * ncarbon / float(nadatoms)
            for element, count in zip(self.atoms[iatom].split(','), self.atomCounts[1:]):
                if element == 'Vc':
                    pass
                elif element == 'H':
                    bindEnergy -= count / float(nadatoms) * eIsolatedH
                else:
                    try:
                        bindEnergy -= count / float(nadatoms) * self.singleE[element]
                    except:
                        bindEnergy -= count / float(nadatoms) * 100

#            print 'struct',struct, 'Energy',structEnergy,'BE',bindEnergy,'nadatoms',nadatoms
#            print ncarbon*energyGraphene/2.0,nH*eIsolatedH,nmetal*self.singleE[iatom]
            self.vdata[iatom,istruct]['BE'] = bindEnergy
            vaspBEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,bindEnergy))  
            #note that the hex monolayers have one atom per cell

            hexEnergy = structEnergy - energyGraphene / 2.0 * ncarbon / float(nadatoms)
            for element, count in zip(self.atoms[iatom].split(','), self.atomCounts[1:]):
                if element == 'Vc':
                    pass
                elif element == 'H':
                    hexEnergy -= count / float(nadatoms) * eH2
                else:
                    try:
                        hexEnergy -= count / float(nadatoms) * self.hexE[element]
                    except:
                        hexEnergy -= count / float(nadatoms) * 100

            self.vdata[iatom,istruct]['HFE'] = hexEnergy
            vaspHFEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,conc,hexEnergy))    
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

            nCarbon =  self.atomCounts[0]

            nadatoms =  float(sum(self.atomCounts) - nCarbon)                       
            self.vdata[iatom,istruct]['nadatoms'] = nadatoms

            conc = [float(float(x) / nadatoms) for x in self.atomCounts[1:]][0]
            self.vdata[iatom,istruct]['conc'] = conc    
                   
            formationEnergy = structEnergy
            for i, count in enumerate(self.atomCounts[1:]):
                formationEnergy -= float(count) / float(nadatoms) * self.pureEnergies[i]

            self.vdata[iatom,istruct]['FE'] = formationEnergy
            formEnergyList.append([formationEnergy, struct])
            istruct += 1                                
        formEnergyList.sort()       
        for pair in formEnergyList:
            sortedStructs.append(pair[1])        
        self.newlyFinished[iatom] = sortedStructs           
        os.chdir(lastDir)
        return nfinished #the position in vdata before adding the data
