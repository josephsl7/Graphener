'''
Created on Aug 26, 2014

@author: eswens13
'''

import os, subprocess
from random import seed
from numpy import zeros
from copy import deepcopy

import Enumerator, Extractor, Structs2Poscar, JobManager, MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo


def readSettingsFile():
    currDir = os.getcwd()
    infile = open(currDir + '/needed_files/settings.in', 'r')
    inlines = []
    for line in infile:
        firstPart = line.strip().split()[0]
        firstChar = list(firstPart)[0]
        if firstChar != '#':
            inlines.append(line.strip())
    infile.close()
    
    atoms = []
    volRange = []
    clusterNums = []
    trainStructs = 0
    fitStructs = 0
    fitSubsets = 0
    plotTitle = "title"
    xlabel = "xlabel"
    ylabel = "ylabel"
    
    for line in inlines:
        if line.split()[0] == 'ATOMS:':
            i = 1
            adatoms = line.split()[1:]
            for adatom in adatoms:
                if adatom == '#':
                    break
                else:
                    atoms.append(adatom)
            
        elif line.split()[0] == 'VOL_RANGE:':
            low = int(line.split()[1])
            high = int(line.split()[2])
            volRange = [low, high]
            
        elif line.split()[0] == 'CLUSTER_NUMS:':
            parts = line.split()
            for i in xrange(1, 11):
                clusterNums.append(int(parts[i]))
        
        elif line.split()[0] == 'TRAINING_STRUCTS:':
            trainStructs = int(line.split()[1])
          
        elif line.split()[0] == 'FITTING_STRUCTS:':
            fitStructs = int(line.split()[1])
            
        elif line.split()[0] == 'STRUCT_SUBSETS:':
            fitSubsets = int(line.split()[1])
        
        elif line.split()[0] == 'PLOT_TITLE:':
            plotTitle = line.split('\'')[1]
        
        elif line.split()[0] == 'XLAB:':
            xlabel = line.split('\'')[1]
        
        elif line.split()[0] == 'YLAB:':
            ylabel = line.split('\'')[1]
    
    return [atoms, volRange, clusterNums, trainStructs, fitStructs, fitSubsets, plotTitle, xlabel, ylabel]

def contains(struct, structList):
    for i in xrange(len(structList)):
        if str(struct) == str(structList[i]):
            return True
    
    return False

def equals(alist, blist):
    clist = deepcopy(alist)
    dlist = deepcopy(blist)
    while len(clist) > 0:
        if len(dlist) == 0:
            return False
        if not contains(clist[0], dlist):
            return False
        else:
            clist.remove(clist[0])
            dlist.remove(clist[0])
    
    if len(dlist) > 0:
        return False
    else:
        return True
    
           
if __name__ == '__main__':
    seed()
    
    [atomList, volRange, clusterNums, trainingStructs, fitStructs, fitSubsets, plotTitle, xlabel, ylabel] = readSettingsFile()
    
    enumerator = Enumerator.Enumerator(atomList, volRange, clusterNums, trainingStructs)
    subprocess.call(['echo','\nEnumerating symmetrically unique structures. . .\n'])
    enumerator.enumerate()
    
    changed = True
    iter = 1
    prevStructs = []
    newStructs = []
    allStructs = []
    while changed:
        changed = False
        
        # Extract the structures from struct_enum.out
        extractor = Extractor.Extractor(atomList)
        if iter == 1:
            extractor.setTrainingStructs()
        elif iter > 1:
            extractor.setStructList(newStructs)
        #prevStructs = extractor.getStructList()
        extractor.extract()
    
        # Convert the extracted pseudo-POSCARs to VASP POSCAR files, make directories for them
        # and put the POSCARs in their corresponding directories.
        subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
        toPoscar = Structs2Poscar.Structs2Poscar(atomList, prevStructs)
        toPoscar.convertOutputsToPoscar()
     
        # Start VASP jobs and wait until they all complete or time out.
        manager = JobManager.JobManager(atomList)
        manager.runLowJobs(prevStructs)
        manager.runNormalJobs(prevStructs)
        manager.runDOSJobs(prevStructs) 
    
        # Create structures.in and structures.holdout files for each atom.
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atomList)
        uncleFileMaker.makeUncleFiles()
        
        # Get all the structs that have been through VASP calculations for each atom.
        structList = uncleFileMaker.getStructureList() 
        
        # Perform a fit to the VASP data in structures.in for each atom.
        fitter = Fitter.Fitter(atomList, fitStructs, fitSubsets)
        fitter.makeFitDirectories()
        fitter.fitVASPData()
    
        # Perform a ground state search on the fit for each atom.    
        gss = GSS.GSS(atomList, volRange, plotTitle, xlabel, ylabel)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch()
        if iter == 1:
            allStructs = gss.getAllGSSStructures()
        
        # Add the 100 structures with the lowest formation energy that nave not been through VASP
        # calculations to the newStructs list for each atom.
        newStructs = []
        added = zeros(len(atomList))
        for i in xrange(len(allStructs)):
            atomStructs = []
            for j in xrange(len(allStructs[i])):
                if added[i] >= 100:
                    break
                elif not contains(allStructs[i][j], structList):
                    atomStructs.append(allStructs[i][j])
                    added[i] += 1
            newStructs.append(atomStructs)
        
        # Remove the lowest 100 structs from allStructs for each atom.
        for i in xrange(len(allStructs)):
            for j in xrange(len(newStructs[i])):
                allStructs[i].remove(newStructs[i][j])
        
        # If one of the atoms has finished the loop, remove it from the atomList.
        toRemove = []
        for i in xrange(len(newStructs)):
            if equals(newStructs[i], prevStructs[i]):
                toRemove.append(i)
        
        if len(toRemove) == len(newStructs):
            changed = False
        else:
            for ind in toRemove:
                atomList.remove(atomList[ind])
                newStructs.remove(newStructs[ind])
                prevStructs.remove(prevStructs[ind])
        
        # Change the new structs at the end of the loop to the previous structs for the 
        # beginning of the next iteration.
        prevStructs = []
        for i in xrange(len(newStructs)):
            prevStructs.append(newStructs[i])
            
        # Keep track of which iteration we're on.
        iter += 1
        
    
    # After the loop has finished, make plots of the lowest formation energy structs
    # for each atom.
    gss.makePlots()
    
    # Should do some analysis after the loop has finished as well.
        

    
        
    
    
 
 
 
 
 
 
 
    