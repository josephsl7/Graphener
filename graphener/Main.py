'''
Created on Aug 26, 2014

@author: eswens13
'''

import os, subprocess
from random import seed
from numpy import zeros

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
            
if __name__ == '__main__':
    seed()
    
    [atomList, volRange, clusterNums, trainingStructs, fitStructs, fitSubsets, plotTitle, xlabel, ylabel] = readSettingsFile()
    
    """enumerator = Enumerator.Enumerator(atomList, volRange, clusterNums, trainingStructs)
    subprocess.call(['echo','\nEnumerating symmetrically unique structures. . .\n'])
    enumerator.enumerate()
    
    extractor = Extractor.Extractor(atomList)
    extractor.setTrainingStructs()"""
    
    changed = True
    iter = 1
    prevStructs = []
    newStructs = []
    allStructs = []
    while changed:
        changed = False
        """if iter > 1:
            extractor.setStructList(newStructs)
        prevStructs = extractor.getStructList()
        extractor.extract()
    
        subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
        toPoscar = Structs2Poscar.Structs2Poscar(atomList, prevStructs)
        toPoscar.convertOutputsToPoscar()
     
        manager = JobManager.JobManager(atomList)
        manager.runLowJobs(prevStructs)
        manager.runNormalJobs(prevStructs)
        manager.runDOSJobs(prevStructs)""" 
    
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atomList)
        uncleFileMaker.makeUncleFiles()
        structList = uncleFileMaker.getStructureList()
        
        fitter = Fitter.Fitter(atomList, fitStructs, fitSubsets)
        fitter.makeFitDirectories()
        fitter.fitVASPData()
    
        gss = GSS.GSS(atomList, volRange, plotTitle, xlabel, ylabel)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch()
        if iter == 1:
            allStructs = gss.getAllGSSStructures()
        
        added = zeros(len(atomList))
        for i in xrange(len(allStructs)):
            atomStructs = []
            for j in xrange(len(allStructs[i])):
                if added[i] >= 100:
                    break
                elif not contains(allStructs[i][j], structList):
                    atomStructs.append(allStructs[i][j])
                    allStructs[i].remove(allStructs[i][j])
                    added[i] += 1
            newStructs.append(atomStructs)
 
        #gss.makePlots()  Do this after the convergence loop has finished
        
        iter += 1
        
        # Need to call extractor.setStructList() before end of loop.
        # Check if the structList changed from what it was before. If so, reset 'changed' to True
        # and repeat the loop.
    
        
    
    
    