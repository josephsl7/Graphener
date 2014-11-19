'''
Created on Aug 26, 2014

@author: eswens13
'''

#notes for Eric:

import os, subprocess,sys,re
from random import seed
from numpy import zeros,array,sqrt,std,amax,amin,int32
from copy import deepcopy

import Enumerator, Extractor, Structs2Poscar, JobManager, MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo     

def contains(struct, alist):
    if len(alist) == 0:
        return False
    
    for i in xrange(len(alist)):
        if str(struct) == str(alist[i]):
            return True
    
    return False

def createEnumPastDir(atoms):
    '''enumpast is a folder for storing past_structs.dat, and enumerating new iid structures 
    that need past_structs.dat. This creates or clears the folder and puts an empty past_structs.dat there'''
    for i,atom in enumerate(atoms):
        atomDir = os.getcwd() + '/' + atom
        vsDir = atomDir + '/enumpast'
        if os.path.isdir(vsDir): subprocess.call(['rm','-r' ,vsDir])                    
        subprocess.call(['mkdir', vsDir])                
        file = open(vsDir + '/past_structs.dat', 'w'); file.close()  #just create it. 

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

def getDiffE(priority, eLast, atoms):
    '''Finds the L1 norm energy change between iterations, weighted by priority'''
    ediff = priority[:,:]['energy'] - eLast
    ediffL1 = zeros(len(atoms))
    for iatom in range(len(atoms)):
        priorsum = sum(priority[iatom,:]['prior'])
        ediff[iatom,:] = abs(ediff[iatom,:]) * priority[iatom,:]['prior']/priorsum
        ediffL1[iatom] = sqrt(sum(ediff[iatom,:]))
    return ediffL1

def parseStartStructures(atoms,startFromExisting,vdata,vstructsFinished):
    """ Writes past_structures.dat that vasp needs, and returns the structures in structures.in, 
    
    This method assumes that the structures.in file will have all of the pure
    structures and will have the following format for the structure identification line:
    graphene str #: (structure number) FE = 0.545, Concentration = .473
    or, for a pure structure:
    PURE M graphene str #: (structure number) FE = 0.0, Concentration = 1.0 """
    istruct = 0
    Nend = vstructsFinished
    startList = []
    lastDir = os.getcwd()
    for i,atom in enumerate(atoms):
        if startFromExisting[i]: #chosen by user in settings file
            subList = []
            startFile = lastDir + '/needed_files/structures.start.' + atoms[i]       
            if os.path.exists(startFile):      
                infile = open(startFile,'r')
                lines = infile.readlines()
                infile.close() 
                outfile =  open(lastDir + '/' + atom + '/enumpast/past_structs.dat', 'w')
                for j in xrange(len(lines)):
                    if list(lines[j].strip().split()[0])[:2] == ['#','-']:
                        if j != len(lines) - 1:
                            structLine = lines[j + 1].strip().split()
                            istruct += 1
                            if structLine[0].lower() == 'pure':
                                outfile.write(str(structLine[5]) + '\n')
                                subList.append(structLine[5])
                                vdata[istruct]['struct'] = structLine[5]
                            else:
                                outfile.write(str(structLine[3]) + '\n')
                                subList.append(structLine[3])
                                vdata[istruct]['struct'] = structLine[3]
                            natomsList = lines[j+6].strip().split()
                            natoms = natomsList[0]+natomsList[1]
                            vdata[istruct]['conc'] = float(natomsList[1]/natoms)
                            vdata[istruct]['energy'] = lines[j+8+natoms].strip().split()[0]                    
                outfile.close()
            else: 
                subprocess.call(['echo',"ERROR: structures.start."+ atoms[i] + ' is missing from needed_files/'])       
        startList.append(subList)
    os.chdir(lastDir)
    return startList, vdata  
 
def readfile(filepath):
        file1 = open(filepath,'r')
        lines = file1.readlines()
        file1.close()
        return lines           

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
#    startFromExisting = []
    runTypes = ['low']
    PriorOrIID = 'p'
    ntrainStructs = 0
    mfitStructs = 0
    nfitSubsets = 0
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

        elif line.split()[0] == 'START_FROM_EXISTING:': #read from a structures.in.start and fit these.  Needs to come after the ATOMS list in the settings lin
            if line.split()[1].lower() == 'n': #"No" or "None" means no atoms will start from existing calcs.
               startFromExisting = [False]*len(atoms) 
            elif line.split()[1].lower() == 'a': #"All" 
               startFromExisting = [True]*len(atoms)
            else:
                try:  
                    parselist = re.search(':(.+?)#', line).group(1).lower().split()[:len(atoms)] #take only enough to match number of atoms, as it might be an old list
                    startFromExisting = [True if i == 't' else False for i in parselist]
                except AttributeError:
                    print "Can't parse START_STRUCTS T's and F's.  Defaulting to no start structs." 
                    startFromExisting = [False]*len(atoms)               
        elif line.split()[0] == 'PRIORITY/IID:':
            if str(line.split()[1]).lower() == 'i':
                PriorOrIID = 'i'

        elif line.split()[0] == 'RUN_TYPES:':
            try:
                runTypes = re.search(':(.+?)#', line).group(1).lower().split() #extracts only text between START_STRUCTS and #
            except AttributeError:
                print "Can't parse RUN_TYPES. Defaulting to LOW precision runs only."             

        elif line.split()[0] == 'N_TRAINING_STRUCTS:':
            ntrainStructs = int(line.split()[1])
          
        elif line.split()[0] == 'M_FITTING_STRUCTS:':
            mfitStructs = int(line.split()[1])
            
        elif line.split()[0] == 'N_STRUCT_SUBSETS:':
            nfitSubsets = int(line.split()[1])
        
        elif line.split()[0] == 'GROW_NUM:':
            growNum = int(line.split()[1])
        
        elif line.split()[0] == 'PLOT_TITLE:':
            plotTitle = line.split('\'')[1]
        
        elif line.split()[0] == 'XLAB:':
            xlabel = line.split('\'')[1]
        
        elif line.split()[0] == 'YLAB:':
            ylabel = line.split('\'')[1]
    
    return [atoms, volRange, clusterNums, startFromExisting, runTypes, PriorOrIID, ntrainStructs, mfitStructs, nfitSubsets, growNum, plotTitle, xlabel, ylabel]

def writeFailedVasp(failedFile, newlyFailed, iteration, atoms):
    try:
        failedFile.write('==============================================================\n')
        failedFile.write('\tIteration: ' + str(iteration) + '\n')
        failedFile.write('==============================================================\n')
        for i in xrange(len(newlyFailed)):
            failedFile.write('\n******************** ' + atoms[i] + ' ********************\n')
            atomLength = len(newlyFailed[i])
            if atomLength == 0:
                failedFile.write('\nNo structures have failed.\n')
            else:
                for j in xrange(len(newlyFailed[i])):
                    if (j + 1) % 20 == 0 or j == len(newlyFailed[i]) - 1:
                        failedFile.write(str(newlyFailed[i][j]) + '\n')
                    else:
                        failedFile.write(str(newlyFailed[i][j]) + ', ')

        
        failedFile.flush()
        os.fsync(failedFile.fileno())
    except IOError:
        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to failed_vasp file. ~~~~~~~~~~\n'])   

def writeLowestGSS(lowestGssFile, newStructs, iteration, atoms):
    try:
        lowestGssFile.write('\n==================================================================\n')
        lowestGssFile.write('\tIteration: ' + str(iteration) + '\n')
        lowestGssFile.write('==================================================================\n')
        for i in xrange(len(newStructs)):
            lowestGssFile.write('\n***************** ' + atoms[i] + ' *******************\n')
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

def writeLowestVasp(lowestStructsFile, vStructsFinished, iteration, atoms):
    try:
        lowestStructsFile.write('==============================================================\n')
        lowestStructsFile.write('\tIteration: ' + str(iteration) + '\n')
        lowestStructsFile.write('==============================================================\n')
        for i in xrange(len(vStructsFinished)):
            lowestStructsFile.write('\n******************** ' + atoms[i] + ' ********************\n')
            atomLength = len(vStructsFinished[i])
            if atomLength >= 100:
                for j in xrange(len(vStructsFinished[i][:100])):
                    if (j + 1) % 20 == 0 or j == 99:
                        lowestStructsFile.write(str(vStructsFinished[i][j]) + '\n')
                    else:
                        lowestStructsFile.write(str(vStructsFinished[i][j]) + ', ')
            elif atomLength == 0:
                lowestStructsFile.write('\nNo structures submitted.\n')
            else:
                for j in xrange(len(vStructsFinished[i][:atomLength])):
                    if (j + 1) % 20 == 0 or j == atomLength - 1:
                        lowestStructsFile.write(str(vStructsFinished[i][j]) + '\n')
                    else:
                        lowestStructsFile.write(str(vStructsFinished[i][j]) + ', ')                       
        lowestStructsFile.flush()
        os.fsync(lowestStructsFile.fileno())
        
    except IOError:
        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to lowest_vasp file. ~~~~~~~~~~\n'])
                
def extractToVasp(iteration,runTypes,atoms):
    ''' Convert the extracted pseudo-POSCARs to VASP POSCAR files, make directories for them
     and put the POSCARs in their corresponding directories. Run VASP'''
    extractor.setStructsFromTraining(iteration, pastStructs)
    vstructsCurrent = extractor.getStructList()
#    initialHoldout(vstructsCurrent,atoms)
    extractor.extract()
    subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
    toPoscar = Structs2Poscar.Structs2Poscar(atoms, vstructsCurrent)
    toPoscar.convertOutputsToPoscar()
    # Start VASP jobs and wait until they all complete or time out.
    manager2 = JobManager.JobManager(atoms)
    if runTypes ==['low']:
        manager2.runLowJobs(vstructsCurrent)
        finalDir = '/'
    elif runTypes ==['low','normal']:
        manager2.runLowJobs(vstructsCurrent)
        manager2.runNormalJobs(vstructsCurrent) 
        finalDir = '/normal'          
    elif runTypes ==['low','normal','dos']:
        manager2.runLowJobs(vstructsCurrent)
        manager2.runNormalJobs(vstructsCurrent)
        manager2.runDOSJobs(vstructsCurrent) 
        finalDir = '/DOS' 
    else:
        print 'Your RUN_TYPES is ', runTypes
        sys.exit('The only supported RUN_TYPES are "low", "low normal" and "low normal DOS"')   
    return finalDir

# -------------------------------------------- MAIN -----------------------------------------------
          
if __name__ == '__main__':
    dir = '/fslhome/bch/cluster_expansion/graphene/testtm2'
    print 'Starting in ', dir
    os.chdir(dir)

    
    seed()
    maxvstructs = 20000 #maximum number of structures for each atom
    vdata = zeros(maxvstructs,dtype = [('struct', int32),('conc', float), ('energy', float), ('FE', float),('BE', float),('HFE', float)]) #data from vasp
    vstructsAll = [] #every structure vasp has attempted, a list for each atom
    vstructsFinished = [] #every structure vasp has finished, a list for each atom
    vstructFailed = [] #every structure vasp has failed, a list for each atom
    vstructsCurrent = [] #the structures vasp has attempted or will attempt this iteration, a list for each atom
#    vstructs = [vstructsAll,vstructsFinished,vstructFailed,vstructsCurrent]
    
    [atoms, volRange, clusterNums, startFromExisting, runTypes, PriorOrIID, ntrainStructs, mfitStructs, nfitSubsets, growNum, plotTitle, xlabel, ylabel] = readSettingsFile()
    uncleOutput = open('uncle_output.txt','w') # All output from UNCLE will be written to this file.
   
    if not os.path.isdir('single_atoms'):
        manager1 = JobManager.JobManager(atoms)
        manager1.runSingleAtoms()
    if not os.path.isdir('hex_monolayer_refs'):
        manager1 = JobManager.JobManager(atoms)
        manager1.runHexMono()
    
    createEnumPastDir(atoms)
    
    enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, ntrainStructs, uncleOutput)
    enumerator.enumerate()
    ntot = enumerator.getNtot(os.getcwd()+'/enum')
    eLast = zeros((ntot,len(atoms)),dtype=float) #energies of last iteration, sorted by structure name
    
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
        extractor = Extractor.Extractor(atoms, uncleOutput, startFromExisting)
        if iteration == 1:
            if startFromExisting.count(False) == 0: #all start from existing
                vstructsFinished,vdata = parseStartStructures(atoms,startFromExisting,vdata,vstructsFinished)
                vstructsAll = vstructsFinished
#                pastStructs = extractor.getPastStructs() 
                finalDir = []
            else: #at least one atom needs calculations
                enumerator.chooseTrainingStructures(iteration,startFromExisting)
#                pastStructs = extractor.getPastStructs()
                extractor.setStructsFromTraining(iteration, vtructsAll)
                finalDir = extractToVasp(iteration,runTypes,atoms)
        elif iteration > 1 and PriorOrIID == 'p':
            extractor.setStructsFromGSS(newStructs)
            finalDir = extractToVasp(iteration,runTypes,atoms)
        elif iteration > 1 and PriorOrIID == 'i':
            enumerator.chooseTrainingStructures(iteration,startFromExisting)
#            pastStructs = extractor.getPastStructs()
            finalDir = extractToVasp(iteration,runTypes,atoms)
        else:
            finalDir = '' #starting exclusively from structures.start
   
        # Create structures.in and structures.holdout files for each atom.
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atoms, startFromExisting, iteration, finalDir)
        print holdoutStructs
        [newlyFinished, newlyFailed, vdata] = uncleFileMaker.makeUncleFiles(iteration, holdoutStructs,vstructsCurrent) 
        # Get all the structs that have been through VASP calculations for each atom. These
        # should be sorted by formation energy during the work done by makeUncleFiles()
        # TODO:  Check the precision of the energy per atom that I take from VASP and put
        # into structures.in.  See if it's coming from low-precision or normal-precision
        # rather than DOS.
#        [newlyFinished, newlyFailed, vdata] = uncleFileMaker.getNewVaspResults(vstructsCurrent,vdata) 
        #BCH replace the above.. Do it all in analyzeNewVasp in makeUncleFiles
        #This will be done above:
#            1. vstructs_current gets analyzed and sorted.  Append to finished and failed
#            3. Put energies and conc from vasp into vdata
        2. Append all of vstructs_current to past_structs.dat file
                
        # Perform a fit to the VASP data in structures.in for each atom.
        fitter = Fitter.Fitter(atoms, mfitStructs, nfitSubsets, vstructsFinished,uncleOutput)
        fitter.makeFitDirectories()
        fitter.fitVASPData(iteration)
    
        # Perform a ground state search on the fit for each atom.    
        gss = GSS.GSS(atoms, volRange, plotTitle, xlabel, ylabel, uncleOutput)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch(iteration)
        gss.makePlots(iteration)
        gssStructs = gss.getAllGSSStructures(iteration, vstructFailed)
        priority = gss.getGssInfo(iteration)
        diffe = getDiffE(priority,eLast,atoms)
        print 'diffe',diffe
        eLast = priority[:,:]['energy']
        # Print the lowest energy structures that have been through VASP calculations to a file.
        writeLowestVasp(lowestStructsFile, vStructsFinished, iteration, atoms)
        
        # Write the all the structures that have failed VASP calculations to a file.
        # TODO: Only write the failed structures that are unique to this iteration to the file.
        writeFailedVasp(failedFile, vstructFailed, iteration, atoms)
        
        # Check the lowest 100 hundred structs from VASP against the lowest 100 structs from UNCLE
        # for each atom.  If they match, then that atom has converged and we remove it from the 
        # lists.
        removeAtoms = []
        removeGss = []
        removeVasp = []
        for i in xrange(len(vStructsFinished)):
            atomLength = len(vStructsFinished[i])
            if atomLength >= 100:
                if equals(vStructsFinished[i][:100], gssStructs[i][:100]):
                    startFromExisting.remove(startFromExisting[i])#??? This should only be used during first interation anyway? 
                    removeAtoms.append(atoms[i])
                    removeGss.append(gssStructs[i])
                    removeVasp.append(vStructsFinished[i])
            elif atomLength != 0:
                # If there are not yet 100 structs that have converged in VASP, just check however
                # many have finished.  If there are no structures in the list, it is the first
                # iteration starting from an existing structures.in.start file so we need to go 
                # through another iteration before checking convergence.                
                if equals(vStructsFinished[i][:atomLength], gssStructs[i][:atomLength]):
                    removeAtoms.append(atoms[i])
                    removeGss.append(gssStructs[i])
                    removeVasp.append(vStructsFinished[i])
        
        for i in xrange(len(removeAtoms)):
            atoms.remove(removeAtoms[i])
        for i in xrange(len(removeGss)):
            gssStructs.remove(removeGss[i])
        for i in xrange(len(removeVasp)):
            vStructsFinished.remove(removeVasp[i])
        
        # If all of the atoms have converged, exit the loop.  Else, keep going.
        if len(atoms) > 0:
            changed = True
                    
        if not changed:
            subprocess.call(['echo','\n----------------- The loop has converged! ---------------'])
            break

        # Add the the number of structures specified by growNum with the lowest formation energy 
        # that have not been through VASP calculations (converged or failed) to the newStructs 
        # list for each remaining atom.
        newStructs = []
        added = zeros(len(atoms))
        for i in xrange(len(gssStructs)):
            atomStructs = []
            for j in xrange(len(gssStructs[i])):
                if added[i] >= growNum:
                    break
                elif not contains(gssStructs[i][j], vStructsFinished[i]):
                    atomStructs.append(str(gssStructs[i][j]))
                    added[i] += 1
            newStructs.append(atomStructs)
        
        # Print the new GSS structures to a file.
        writeLowestGSS(lowestGssFile, newStructs, iteration, atoms)
        
        # Set holdoutStructs to the vStructsFinished that were used in the last iteration of the loop.
        if len(vStructsFinished) > 100:
            holdoutStructs = deepcopy(vStructsFinished[:100])
        else:
            holdoutStructs = deepcopy(vStructsFinished[:len(vStructsFinished)])

        # Keep track of which iteration we're on.
        iteration += 1
        
        if iteration == 4:
            break
    
    uncleOutput.close()
    lowestStructsFile.close()
    lowestGssFile.close()
    # Should do some analysis after the loop has finished as well. """

    subprocess.call(['echo','\n---------- PROGRAM ENDED ----------\n'])
        

    
        
    
    
 
 
 
 
 
 
 
    