'''
Created on Aug 26, 2014

@author: eswens13
'''

#notes for Eric:
#removed "base".  replaced with startStructs for the boolean, and start.in for base.in

import os, subprocess,sys
from random import seed
from numpy import zeros,sqrt,std,amax,amin
from copy import deepcopy

import Enumerator, Extractor, Structs2Poscar, JobManager, MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo     

def contains(struct, alist):
    if len(alist) == 0:
        return False
    
    for i in xrange(len(alist)):
        if str(struct) == str(alist[i]):
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
            toRemove = clist[0]
            clist.remove(toRemove)
            dlist.remove(toRemove)
    
    if len(dlist) > 0:
        return False
    else:
        return True
    
def getDiffE(priority, eLast, atomList):
    '''Finds the L1 norm energy change between iterations, weighted by priority'''
    ediff = priority[:,:]['energy'] - eLast
    ediffL1 = zeros(len(atomList))
    for iatom in range(len(atomList)):
        priorsum = sum(priority[iatom,:]['prior'])
        ediff[iatom,:] = abs(ediff[iatom,:]) * priority[iatom,:]['prior']/priorsum
        ediffL1[iatom] = sqrt(sum(ediff[iatom,:]))
    return ediffL1

def parseStartStructures(atomList):
    """ Right now this method assumes that the structures.in file will have all of the pure
    structures and will have the following format for the structure identification line:
    graphene str #: (structure number) FE = 0.545, Concentration = .473
    or, for a pure structure:
    PURE M graphene str #: (structure number) FE = 0.0, Concentration = 1.0 """
    lastDir = os.getcwd()
    os.chdir(lastDir + '/needed_files')
    startFromExisting = []
    found = False
    for i in xrange(len(atomList)):
        if (os.path.exists('structures.start.' + atomList[i])):
            infile = open('structures.start.' + atomList[i],'r')
            lines = infile.readlines()
            infile.close()
               
            print lastDir + '/' + atomList[i] + '/past/past_structs.dat'
            outfile =  open(lastDir + '/' + atomList[i] + '/past/past_structs.dat', 'w')
            for j in xrange(len(lines)):
                if list(lines[j].strip().split()[0])[:2] == ['#','-']:
                    if j != len(lines) - 1:
                        structLine = lines[j + 1].strip().split()
                        if structLine[0].lower() == 'pure':
                            outfile.write(str(structLine[5]) + '\n')
                        else:
                            outfile.write(str(structLine[3]) + '\n')
            outfile.close()
            startFromExisting.append(True)
            found = True
        else:
            startFromExisting.append(False)
    
    os.chdir(lastDir)

    if not found:
        subprocess.call(['echo','\tNo structures.start files found!  Proceeding without them. . .'])
    return startFromExisting   

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
    startStructs = False
    PriorOrIID = 'g'
    trainStructs = 0
    fitStructs = 0
    fitSubsets = 0
    growNum = 0
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
            high = int(line.split()[1])
            volRange = [1, high]
            
        elif line.split()[0] == 'CLUSTER_NUMS:':
            parts = line.split()
            for i in xrange(1, 11):
                clusterNums.append(int(parts[i]))

        elif line.split()[0] == 'START_STRUCTS:': #read from a structures.in.start and fit these
            if str(line.split()[1]).lower() == 'true':
                startStructs = True
        
        elif line.split()[0] == 'PRIORITY/IID:':
            if str(line.split()[1]).lower() == 'i':
                PriorOrIID = 'i'

        elif line.split()[0] == 'TRAINING_STRUCTS:':
            trainStructs = int(line.split()[1])
          
        elif line.split()[0] == 'FITTING_STRUCTS:':
            fitStructs = int(line.split()[1])
            
        elif line.split()[0] == 'STRUCT_SUBSETS:':
            fitSubsets = int(line.split()[1])
        
        elif line.split()[0] == 'GROW_NUM:':
            growNum = int(line.split()[1])
        
        elif line.split()[0] == 'PLOT_TITLE:':
            plotTitle = line.split('\'')[1]
        
        elif line.split()[0] == 'XLAB:':
            xlabel = line.split('\'')[1]
        
        elif line.split()[0] == 'YLAB:':
            ylabel = line.split('\'')[1]
    
    return [atoms, volRange, clusterNums, startStructs, PriorOrIID, trainStructs, fitStructs, fitSubsets, growNum, plotTitle, xlabel, ylabel]

def writeFailedVasp(failedFile, failedStructs, iteration, atomList):
    try:
        failedFile.write('==============================================================\n')
        failedFile.write('\tIteration: ' + str(iteration) + '\n')
        failedFile.write('==============================================================\n')
        for i in xrange(len(failedStructs)):
            failedFile.write('\n******************** ' + atomList[i] + ' ********************\n')
            atomLength = len(failedStructs[i])
            if atomLength == 0:
                failedFile.write('\nNo structures have failed.\n')
            else:
                for j in xrange(len(failedStructs[i])):
                    if (j + 1) % 20 == 0 or j == len(failedStructs[i]) - 1:
                        failedFile.write(str(failedStructs[i][j]) + '\n')
                    else:
                        failedFile.write(str(failedStructs[i][j]) + ', ')

        
        failedFile.flush()
        os.fsync(failedFile.fileno())
    except IOError:
        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to failed_vasp file. ~~~~~~~~~~\n'])   

def writeLowestGSS(lowestGssFile, newStructs, iteration, atomList):
    try:
        lowestGssFile.write('\n==================================================================\n')
        lowestGssFile.write('\tIteration: ' + str(iteration) + '\n')
        lowestGssFile.write('==================================================================\n')
        for i in xrange(len(newStructs)):
            lowestGssFile.write('\n***************** ' + atomList[i] + ' *******************\n')
            atomLength = len(newStructs[i])
            if atomLength >= 100:
                for j in xrange(len(newStructs[i])):
                    if (j + 1) % 20 == 0 or j == 99:
                        lowestGssFile.write(str(newStructs[i][j]) + '\n')
                    else:
                        lowestGssFile.write(str(newStructs[i][j]) + ', ')
        
        lowestGssFile.flush()
        os.fsync(lowestGssFile.fileno())
    except IOError:
        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to lowest_gss file. ~~~~~~~~~~\n'])

def writeLowestVasp(lowestStructsFile, vaspStructs, iteration, atomList):
    try:
        lowestStructsFile.write('==============================================================\n')
        lowestStructsFile.write('\tIteration: ' + str(iteration) + '\n')
        lowestStructsFile.write('==============================================================\n')
        for i in xrange(len(vaspStructs)):
            lowestStructsFile.write('\n******************** ' + atomList[i] + ' ********************\n')
            atomLength = len(vaspStructs[i])
            if atomLength >= 100:
                for j in xrange(len(vaspStructs[i][:100])):
                    if (j + 1) % 20 == 0 or j == 99:
                        lowestStructsFile.write(str(vaspStructs[i][j]) + '\n')
                    else:
                        lowestStructsFile.write(str(vaspStructs[i][j]) + ', ')
            elif atomLength == 0:
                lowestStructsFile.write('\nNo structures submitted.\n')
            else:
                for j in xrange(len(vaspStructs[i][:atomLength])):
                    if (j + 1) % 20 == 0 or j == atomLength - 1:
                        lowestStructsFile.write(str(vaspStructs[i][j]) + '\n')
                    else:
                        lowestStructsFile.write(str(vaspStructs[i][j]) + ', ')                       
        lowestStructsFile.flush()
        os.fsync(lowestStructsFile.fileno())
        
    except IOError:
        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to lowest_vasp file. ~~~~~~~~~~\n'])

# -------------------------------------------- MAIN -----------------------------------------------
          
if __name__ == '__main__':
    dir = '/fslhome/bch/cluster_expansion/graphene/testtm2'
    print 'Starting in ', dir
    os.chdir(dir)

    
    
    seed()
    
    [atomList, volRange, clusterNums, startStructs, PriorOrIID, trainingStructs, fitStructs, fitSubsets, growNum, plotTitle, xlabel, ylabel] = readSettingsFile()
    uncleOutput = open('uncle_output.txt','w') # All output from UNCLE will be written to this file.
    
    startFromExisting = []
    if startStructs:
        startFromExisting = parseStartStructures(atomList)
    else:
        for i in xrange(len(atomList)):
            startFromExisting.append(False)
    
    enumerator = Enumerator.Enumerator(atomList, volRange, clusterNums, trainingStructs, uncleOutput)
    enumerator.enumerate()
    ntot = enumerator.getNtot(os.getcwd()+'/enum')
    eLast = zeros((ntot,len(atomList)),dtype=float) #energies of last iteration, sorted by structure name
    
    changed = True
    iteration = 1
    
    newStructs = []
    gssStructs = []
    pastStructs = []
    holdoutStructs = []
    lowestStructsFile = open('lowest_vasp.txt','w')
    lowestGssFile = open('lowest_gss.txt','w')
    failedFile = open('failed_vasp.txt','w')
    
    while changed:
        changed = False
        
        subprocess.call(['echo','\n========================================================'])
        subprocess.call(['echo','\t\tIteration ' + str(iteration)])
        subprocess.call(['echo','========================================================\n'])
        
        # Extract the pseudo-POSCARs from struct_enum.out
        extractor = Extractor.Extractor(atomList, uncleOutput, startFromExisting)
        if iteration == 1:
            enumerator.chooseTrainingStructures()
            pastStructs = extractor.getPastStructs()
            extractor.setStructsFromTraining(iteration, pastStructs)
        elif iteration > 1 and PriorOrIID == 'p':
            extractor.setStructsFromGSS(newStructs)
        elif iteration > 1 and PriorOrIID == 'i':
            enumerator.chooseTrainingStructures()
            pastStructs = extractor.getPastStructs()
            extractor.setStructsFromTraining(iteration, pastStructs)

        toCalculate = extractor.getStructList()

        extractor.extract()
    
        # Convert the extracted pseudo-POSCARs to VASP POSCAR files, make directories for them
        # and put the POSCARs in their corresponding directories.
        subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
        toPoscar = Structs2Poscar.Structs2Poscar(atomList, toCalculate)
        toPoscar.convertOutputsToPoscar()
     
        # Start VASP jobs and wait until they all complete or time out.
        manager = JobManager.JobManager(atomList)
        manager.runLowJobs(toCalculate)
        finalDir = '/'  #if want low precision runs only        
#        manager.runNormalJobs(toCalculate)
#        manager.runDOSJobs(toCalculate)
    
        # Create structures.in and structures.holdout files for each atom.
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atomList, startFromExisting, iteration)
        uncleFileMaker.finalDir = finalDir #bch
        uncleFileMaker.makeUncleFiles(iteration, holdoutStructs)
        
        # Get all the structs that have been through VASP calculations for each atom. These
        # should be sorted by formation energy during the work done by makeUncleFiles()
        # TODO:  Check the precision of the energy per atom that I take from VASP and put
        # into structures.in.  See if it's coming from low-precision or normal-precision
        # rather than DOS.
        [vaspStructs, failedStructs] = uncleFileMaker.getStructureList()
        structuresInLengths = uncleFileMaker.getStructuresInLengths() 
        
        # Perform a fit to the VASP data in structures.in for each atom.
        fitter = Fitter.Fitter(atomList, fitStructs, fitSubsets, structuresInLengths, uncleOutput)
        fitter.makeFitDirectories()
        fitter.fitVASPData(iteration)
    
        # Perform a ground state search on the fit for each atom.    
        gss = GSS.GSS(atomList, volRange, plotTitle, xlabel, ylabel, uncleOutput)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch(iteration)
        gss.makePlots(iteration)
        gssStructs = gss.getAllGSSStructures(iteration, failedStructs)
        priority = gss.getGssInfo(iteration)
        diffe = getDiffE(priority,eLast,atomList)
        print 'diffe',diffe
        eLast = priority[:,:]['energy']
        # Print the lowest energy structures that have been through VASP calculations to a file.
        writeLowestVasp(lowestStructsFile, vaspStructs, iteration, atomList)
        
        # Write the all the structures that have failed VASP calculations to a file.
        # TODO: Only write the failed structures that are unique to this iteration to the file.
        writeFailedVasp(failedFile, failedStructs, iteration, atomList)
        
        # Check the lowest 100 hundred structs from VASP against the lowest 100 structs from UNCLE
        # for each atom.  If they match, then that atom has converged and we remove it from the 
        # lists.
        removeAtoms = []
        removeGss = []
        removeVasp = []
        for i in xrange(len(vaspStructs)):
            atomLength = len(vaspStructs[i])
            if atomLength >= 100:
                if equals(vaspStructs[i][:100], gssStructs[i][:100]):
                    startFromExisting.remove(startFromExisting[i])
                    removeAtoms.append(atomList[i])
                    removeGss.append(gssStructs[i])
                    removeVasp.append(vaspStructs[i])
            elif atomLength != 0:
                # If there are not yet 100 structs that have converged in VASP, just check however
                # many have finished.  If there are no structures in the list, it is the first
                # iteration starting from an existing structures.in.start file so we need to go 
                # through another iteration before checking convergence.                
                if equals(vaspStructs[i][:atomLength], gssStructs[i][:atomLength]):
                    removeAtoms.append(atomList[i])
                    removeGss.append(gssStructs[i])
                    removeVasp.append(vaspStructs[i])
        
        for i in xrange(len(removeAtoms)):
            atomList.remove(removeAtoms[i])
        for i in xrange(len(removeGss)):
            gssStructs.remove(removeGss[i])
        for i in xrange(len(removeVasp)):
            vaspStructs.remove(removeVasp[i])
        
        # If all of the atoms have converged, exit the loop.  Else, keep going.
        if len(atomList) > 0:
            changed = True
                    
        if not changed:
            subprocess.call(['echo','\n----------------- The loop has converged! ---------------'])
            break

        # Add the the number of structures specified by growNum with the lowest formation energy 
        # that have not been through VASP calculations (converged or failed) to the newStructs 
        # list for each remaining atom.
        newStructs = []
        added = zeros(len(atomList))
        for i in xrange(len(gssStructs)):
            atomStructs = []
            for j in xrange(len(gssStructs[i])):
                if added[i] >= growNum:
                    break
                elif not contains(gssStructs[i][j], vaspStructs[i]):
                    atomStructs.append(str(gssStructs[i][j]))
                    added[i] += 1
            newStructs.append(atomStructs)
        
        # Print the new GSS structures to a file.
        writeLowestGSS(lowestGssFile, newStructs, iteration, atomList)
        
        # Set holdoutStructs to the vaspStructs that were used in the last iteration of the loop.
        if len(vaspStructs) > 100:
            holdoutStructs = deepcopy(vaspStructs[:100])
        else:
            holdoutStructs = deepcopy(vaspStructs[:len(vaspStructs)])

        # Keep track of which iteration we're on.
        iteration += 1
        
        if iteration == 4:
            break
    
    uncleOutput.close()
    lowestStructsFile.close()
    lowestGssFile.close()
    # Should do some analysis after the loop has finished as well. """

    subprocess.call(['echo','\n---------- PROGRAM ENDED ----------\n'])
        

    
        
    
    
 
 
 
 
 
 
 
    