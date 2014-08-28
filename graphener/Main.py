'''
Created on Aug 26, 2014

@author: eswens13
'''

import os

import Enumerator, Extractor, Structs2Poscar, JobManager, Analyzer, DistanceInfo


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
    
    return [atoms, volRange, clusterNums, trainStructs, fitStructs, fitSubsets]
            

if __name__ == '__main__':
    [atomList, volRange, clusterNums, trainingStructs, fitStructs, fitSubsets] = readSettingsFile()
    
    enumerator = Enumerator.Enumerator(atomList, volRange)
    print "\nEnumerating symmetrically unique structures. . .\n"
    enumerator.enumerate()
    
    extractor = Extractor.Extractor(clusterNums, trainingStructs)
    extractor.extract()
    
    print "\nConverting outputs to VASP inputs. . .\n"
    toPoscar = Structs2Poscar.Structs2Poscar(atomList)
    toPoscar.convertOutputsToPoscar()
    
    manager = JobManager.JobManager(atomList)
    manager.runLowJobs()
    manager.runNormalJobs()
    manager.runDOSJobs()
    
    """print "\nAnalyzing movement during relaxation. . .\n"
    distanceInfo = DistanceInfo.DistanceInfo(atomList)
    distanceInfo.getDistanceInfo()
    
    print "\nAnalyzing convergence, CPU time, and energies of the structures. . .\n"
    analyzer = Analyzer.Analyzer(atomList)
    analyzer.analyze()"""
    
        
    
    
    