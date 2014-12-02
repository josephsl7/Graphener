'''
Created on Aug 26, 2014

@author: eswens13
'''

import os, subprocess,sys,re
from random import seed
from numpy import zeros,array,sqrt,std,amax,amin,int32,sort,count_nonzero,delete
from copy import deepcopy

import Enumerator, Extractor, StructsToPoscar, JobManager, MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo     
from MainMethods import *


# -------------------------------------------- MAIN -----------------------------------------------
          
if __name__ == '__main__':
    maindir = '/fslhome/bch/cluster_expansion/graphene/testtm2'  
##    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1.continue'
#    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1'
#    maindir = os.getcwd()
    subprocess.call(['echo','Starting in ' + maindir])
    
    os.chdir(maindir)
    pathMax = maindir + '/needed_files/diffMax'
    if os.path.exists(pathMax):
        diffMax = float(readfile(pathMax)[0].strip()) 
        subprocess.call(['echo','\ndiffMax read from file (convergence criterium): '+str(diffMax)])
    else:
        subprocess.call(['echo','\nMissing diffMax file with the convergence criterium diffMax in eV. Stopping'])
        sys.exit('Stop')
    seed()
    [atoms, volRange, clusterNums, runTypes, PriorOrIID, niid, mfitStructs, nfitSubsets, priorNum, plotTitle, xlabel, ylabel,restartTimeout] = readSettingsFile()
    uncleOutput = open('uncle_output.txt','w') # All output from UNCLE will be written to this file.
    natoms = len(atoms)
    vstructsAll = [[]]*natoms #every structure vasp has attempted before this iteration, a list for each atom
    vstructsToStart = [[]]*natoms #the structures vasp must create folders for and start this iteration, a list for each atom   
    newlyFinished = [[]]*natoms
    newlyFailed = [[]]*natoms
    newlyToRestart = [[]]*natoms
    newStructsPrior = [[]]*natoms  #structures from previous iteration with highest priority for vasp calculation
    holdoutStructs = [[]]*natoms

    if not os.path.isdir('single_atoms'):
        manager1 = JobManager.JobManager(atoms)
        manager1.runSingleAtoms()
    if not os.path.isdir('hex_monolayer_refs'):
        manager1 = JobManager.JobManager(atoms)
        manager1.runHexMono()

#    vstructsTemp,vdata = parseStartStructures(atoms,startFromExisting,vstructsFinished,vdata,)
#    vstructsFinished = joinLists(vstructsFinished,vstructsTemp)
#    vstructsAll = joinLists(vstructsAll,vstructsTemp)

    [vstructsFinished, vstructsRestart, vstructsFailed,startFromExisting,vdata] = checkInitialFolders(atoms,restartTimeout) #assigns all struct folders to either finished, restart, or failed''' 
    vstructsAll = joinLists(vstructsFinished,vstructsFailed)
    vstructsAll = joinLists(vstructsAll,vstructsRestart)  #may not yet include from structures.start
    pastStructsUpdate(vstructsFinished,atoms)  
    
    enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, niid, uncleOutput)
    enumerator.enumerate()
    ntot = enumerator.getNtot(os.getcwd()+'/enum') #number of all enumerated structures
    energiesLast = zeros((natoms,ntot),dtype=float) #energies of last iteration, sorted by structure name

    createEnumPastDir(atoms)
    
    converged = False
    iteration = 1

#======================================  Iteration loop =======================================
    while not converged:
        changed = False
        
        subprocess.call(['echo','\n========================================================'])
        subprocess.call(['echo','\t\tIteration ' + str(iteration)])
        subprocess.call(['echo','========================================================\n'])
        #since atoms may have changed, have to reinitialize this
        enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, niid, uncleOutput) 
        # Extract the pseudo-POSCARs from struct_enum.out
        extractor = Extractor.Extractor(atoms, uncleOutput, startFromExisting)           
        if iteration == 1:
            vstructsToStart = enumerator.chooseTrainingStructures(iteration,startFromExisting)
            vstructsToStart = extractor.checkPureInCurrent(iteration,vstructsToStart,vstructsFinished)           
        elif iteration > 1 and PriorOrIID == 'p':
            vstructsToStart = newStructsPrior  #from previous iteration 
        elif iteration > 1 and PriorOrIID == 'i':
            vstructsToStart = enumerator.chooseTrainingStructures(iteration,startFromExisting)   
        contcarToPoscar(vstructsRestart,atoms,iteration)
        vstructsToRun = joinLists(vstructsRestart,vstructsToStart)
        pastStructsUpdate(vstructsToStart,atoms)
        finalDir = extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsToStart,vstructsToRun)           
        # Create structures.in and structures.holdout files for each atom.
        os.chdir(maindir) #FIX this...shouldn't need it.
#        print "vstructsToStart IS SET FOR ITERATION 1"
#        if iteration == 1: vstructsToStart = [[],['1','3','435']]
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atoms, startFromExisting, iteration, finalDir, restartTimeout) 
        [newlyFinished, newlyFailed,newlyToRestart,vdata] = uncleFileMaker.makeUncleFiles(iteration, holdoutStructs,vstructsToRun,vdata) 
        #update the vstructs lists and past_structs files       
        print '[newlyFinished, newlyFailed, newlyToRestart]'
        print newlyFinished;print newlyFailed;print newlyToRestart
#        for iatom, atom in enumerate(atoms):
#            print atom, 'vstructsFinished in Main0', vstructsFinished[iatom]

        vstructsFinished = joinLists(vstructsFinished,newlyFinished)
        vstructsFailed = joinLists(vstructsFailed,newlyFailed)
        vstructsAll = joinLists(vstructsFinished,vstructsFailed)
        os.chdir(maindir) #FIX this...shouldn't need it.  


#        for iatom, atom in enumerate(atoms):
#            print atom, 'vstructsFinished in Main1', vstructsFinished[iatom]
        
        # Get all the structs that have been through VASP calculations for each atom. These
        # should be sorted by formation energy during the work done by makeUncleFiles()
        # TODO:  Check the precision of the energy per atom that I take from VASP and put
        # into structures.in.  See if it's coming from low-precision or normal-precision
        # rather than DOS.
                
        # Perform a fit to the VASP data in structures.in for each atom.
        fitter = Fitter.Fitter(atoms, mfitStructs, nfitSubsets, vstructsFinished,uncleOutput)
        fitter.makeFitDirectories()
        if iteration == 1 and startFromExisting.count(True) > 0: #at least one starts from existing
            fitter.holdoutFromIn(atoms,startFromExisting)
        fitter.fitVASPData(iteration)
    
        # Perform a ground state search on the fit for each atom.    
        gss = GSS.GSS(atoms, volRange, plotTitle, xlabel, ylabel, vstructsFinished,uncleOutput)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch(iteration)
        gss.makePlots(iteration)
        #I have to have everything in a fixed order,  to use energiesLast
        #get the priority of each structure in each atom
        priorities = gss.getGssInfo(iteration,vstructsFailed) #first structure listed is highest priority

        #test weighted energy difference since last iteration
        diffe = getDiffE(priorities,energiesLast,atoms)
        energiesLast = sort(priorities,order = ['struct'])['FE'] #to be used next iteration
        print 'energiesLast', energiesLast
        # Determine which atoms need to continue
#        diffMax = 0.01 # 10 meV convergence criterion
        diffMax = float(readfile(maindir + '/needed_files/diffMax')[0].strip()) #read from file called 'diffMax', so we can experiment with it during a long run
        print 'Weighted energy changes'
        atomnumbers = range(natoms)
        rmAtoms = []
        for iatom, atom in enumerate(atoms):
            print atom, diffe[iatom], 'eV'
            if diffe[iatom]<diffMax: #Converged
                print 'Atom has converged:', atom
                atoms.remove(atom)
                atomnumbers.remove(iatom)
                rmAtoms.append(iatom)
                print 'new atoms',atoms
#            print'BLOCKING remove atoms for no finished!'
            elif len(newlyFinished[iatom]) == 0 and len(vstructsToStart[iatom]) > 0:
                print 'Atom has only failed structures:', atom,'has been stopped.'
                atoms.remove(atom)
                atomnumbers.remove(iatom)
                rmAtoms.append(iatom)                
        natoms = len(atoms) #could be lower now

        if len(rmAtoms)>0:
            vstructsFinished = multiDelete(vstructsFinished,rmAtoms)
            vstructsFailed = multiDelete(vstructsFailed,rmAtoms)
            vstructsAll = multiDelete(vstructsAll,rmAtoms)
            vdata = delete(vdata,rmAtoms,axis=0)        
            energiesLast = delete(priorities,rmAtoms,axis=0)['FE'] #    priorities[ikeep,:]['FE']

#        print vdata.shape
#        print 'vdata0';print vdata[0,:]['struct'][:100]
#        vdata2 = deepcopy(vdata[atomnumbers[0],:]) #get this started
#        print 'vdata2';print vdata2[:]['struct'][:100]
#        print vdata2.shape
#        for i in atomnumbers[1:]:
#            vdata2 = vdata2.append(deepcopy(vdata[i,:]),axis=0)
#        print vdata2.shape
#        vdata = vdata2
#        print 'vdata';print vdata[0,:]['struct'][:100]


        newStructsPrior = getFromPriorities(priorities, priorNum, vstructsAll,atoms)
        if natoms == 0:
            subprocess.call(['echo','\n----------------- All atoms have converged! ---------------'])
            converged = True

        # Set holdoutStructs for the next iteration to the 100 (or less) lowest vasp structures
         
        for iatom in range(natoms):
            nfinished = count_nonzero(vdata[iatom,:100]['struct'])
            if  nfinished > 100:
                holdoutStructs[iatom] = vdata[iatom,:100]['struct'].tolist()
            else:
                holdoutStructs[iatom] = vdata[iatom,:nfinished]['struct'].tolist() 

        # Go to next iteration
        iteration += 1
        
#        if iteration == 4:
#            break
    
    uncleOutput.close()
#    lowestStructsFile.close()
#    lowestGssFile.close()
    # Should do some analysis after the loop has finished as well. """

    subprocess.call(['echo','\n---------- PROGRAM ENDED ----------\n'])
        

    
        
    
    
 
 
 
 
 
 
 
    