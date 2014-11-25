'''
Created on Aug 26, 2014

@author: eswens13
'''

import os, subprocess,sys,re
from random import seed
from numpy import zeros,array,sqrt,std,amax,amin,int32,sort,count_nonzero,delete
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

def extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsCurrent):
    ''' Convert the extracted pseudo-POSCARs to VASP POSCAR files, make directories for them
     and put the POSCARs in their corresponding directories. Run VASP'''
#    vstructsCurrent = extractor.getStructList()

    extractor.extract(vstructsCurrent)
    if len(vstructsCurrent)>0:
        subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
        toPoscar = Structs2Poscar.structs2Poscar(atoms, vstructsCurrent)
        toPoscar.convertOutputsToPoscar()
    # Start VASP jobs and wait until they all complete or time out.
            # Start VASP jobs and wait until they all complete or time out.
        manager2 = JobManager.JobManager(atoms)
        if runTypes ==['low']:
    #        print'block manager2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11'
            manager2.runLowJobs(vstructsCurrent)
        elif runTypes ==['low','normal']:
            manager2.runLowJobs(vstructsCurrent)
            manager2.runNormalJobs(vstructsCurrent)          
        elif runTypes ==['low','normal','dos']:
            manager2.runLowJobs(vstructsCurrent)
            manager2.runNormalJobs(vstructsCurrent)
            manager2.runDOSJobs(vstructsCurrent) 
        else:
            print 'Your RUN_TYPES is ', runTypes
            sys.exit('The only supported RUN_TYPES are "low", "low normal" and "low normal DOS"') 
    if runTypes ==['low']:
        finalDir = '/'
    elif runTypes ==['low','normal']:
        finalDir = '/normal'          
    elif runTypes ==['low','normal','dos']:
        finalDir = '/DOS' 
    else:
        print 'Your RUN_TYPES is ', runTypes
        sys.exit('The only supported RUN_TYPES are "low", "low normal" and "low normal DOS"')   
    return finalDir

def getFromPriorities(priorities, N, vstructsAll, natoms):
    '''take N structures new to vasp of highest priority'''
    newstructs = [[]]*natoms
    nNew = 0 
    istruct = 0 
    for iatom in range(natoms):
        sublist = []
        while nNew < N:
            struct = str(priorities[iatom,istruct]['struct'])
            if not struct in vstructsAll[iatom]:
                sublist.append(struct)
                nNew += 1
            istruct += 1
        newstructs[iatom] = sublist
    return newstructs
            
def getDiffE(priority, energiesLast, atoms):
    '''Finds the L1 norm energy change between iterations, weighted by priority'''
    ediff = priority[:,:]['FE'] - energiesLast
    ediffL1 = zeros(len(atoms))
    for iatom in range(len(atoms)):
        priorsum = sum(priority[iatom,:]['prior'])
        ediffL1[iatom] = sqrt(sum(abs(ediff[iatom,:]) * priority[iatom,:]['prior']/priorsum))
    return ediffL1

def joinLists(list1,list2):
    '''Joins lists of the [[sublist1],[sublist2],[sublist3]].  List1 and 2 must have the
    same length, but can different length sublists'''
    list3=[]
    for i in range(len(list1)):
        list3.append(list1[i]+list2[i])
    return list3

def multiDelete(list_, args):
    indexes = sorted(args, reverse=True)
    for index in indexes:
        del list_[index]
    return list_

def parseStartStructures(atoms,startFromExisting,vstructsFinished,vdata):
    """ Writes past_structures.dat that vasp needs, and returns the structures in structures.in, 
    
    This method assumes that the structures.in file will have all of the pure
    structures and will have the following format for the structure identification line:
    graphene str #: (structure number) FE = 0.545, Concentration = .473
    or, for a pure structure:
    PURE M graphene str #: (structure number) FE = 0.0, Concentration = 1.0 """
    
    istruct = 0
    startList = []
    lastDir = os.getcwd()
    for iatom,atom in enumerate(atoms):
        if startFromExisting[iatom]:  
            subList = []
            startFile = lastDir + '/needed_files/structures.start.' + atoms[iatom]
            subprocess.call(['cp',startFile, lastDir + '/' + atom + '/structures.start'])       
            if os.path.exists(startFile):      
                infile = open(startFile,'r')
                lines = infile.readlines()
                infile.close() 
                outfile =  open(lastDir + '/' + atom + '/enumpast/past_structs.dat', 'w')
                for j in xrange(len(lines)):
                    if list(lines[j].strip().split()[0])[:2] == ['#','-']:
                        if j != len(lines) - 1:
                            structLine = lines[j + 1].strip().split()
                            if structLine[0].lower() == 'pure':
                                outfile.write(str(structLine[5]) + '\n')
                                subList.append(structLine[5])
                                vdata[iatom,istruct]['struct'] = structLine[5]
                            else:
                                outfile.write(str(structLine[3]) + '\n')
                                subList.append(structLine[3])
                                vdata[iatom,istruct]['struct'] = structLine[3]
                            natomsList = [int(i) for i in lines[j+6].strip().split()]
                            natoms = natomsList[0]+natomsList[1]
                            vdata[iatom,istruct]['natoms'] = natoms
                            vdata[iatom,istruct]['conc'] = float(natomsList[1]/float(natoms))
                            vdata[iatom,istruct]['energy'] = float(lines[j+9+natoms].strip().split()[0])  
                            istruct += 1                  
                outfile.close()
            else: 
                subprocess.call(['echo',"ERROR: structures.start."+ atoms[iatom] + ' is missing from needed_files/'])       
        startList.append(subList)
        
        
    os.chdir(lastDir)
    return startList, vdata
  
def pastStructsUpdate(vstructsCurrent,atoms):
    for iatom,atom in enumerate(atoms):
        file = open( atom +'/enumpast/past_structs.dat','a')
        if len(vstructsCurrent[iatom])>0:
            for struct in vstructsCurrent[iatom]:
                file.write(str(struct)+'\n')
        file.close
     
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
    priorNum = 0
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

        elif line.split()[0] == 'START_FROM_EXISTING:': #read from a structures.start and fit these.  Needs to come after the ATOMS list in the settings lin
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
        
        elif line.split()[0] == 'PRIOR_NUM:':
            priorNum = int(line.split()[1])
        
        elif line.split()[0] == 'PLOT_TITLE:':
            plotTitle = line.split('\'')[1]
        
        elif line.split()[0] == 'XLAB:':
            xlabel = line.split('\'')[1]
        
        elif line.split()[0] == 'YLAB:':
            ylabel = line.split('\'')[1]
    
    return [atoms, volRange, clusterNums, startFromExisting, runTypes, PriorOrIID, ntrainStructs, mfitStructs, nfitSubsets, priorNum, plotTitle, xlabel, ylabel]

def writeFailedVasp(failedFile, newlyFailed, iteration, atoms):
    try:
        failedFile.write('==============================================================\n')
        failedFile.write('\tIteration: ' + str(iteration) + '\n')
        failedFile.write('==============================================================\n')
        for iatom in xrange(len(atoms)):
            failedFile.write('\n******************** ' + atoms[iatom] + ' ********************\n')
            atomLength = len(newlyFailed[iatom])
            if atomLength == 0:
                failedFile.write('\nNo structures have failed.\n')
            else:
                for j in xrange(len(newlyFailed[iatom])):
                    if (j + 1) % 20 == 0 or j == len(newlyFailed[iatom]) - 1:
                        failedFile.write(str(newlyFailed[iatom][j]) + '\n')
                    else:
                        failedFile.write(str(newlyFailed[iatom][j]) + ', ')

        
        failedFile.flush()
        os.fsync(failedFile.fileno())
    except IOError:
        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to failed_vasp file. ~~~~~~~~~~\n'])   

#def writeLowestGSS(lowestGssFile, newStructs, iteration, atoms):
#    try:
#        lowestGssFile.write('\n==================================================================\n')
#        lowestGssFile.write('\tIteration: ' + str(iteration) + '\n')
#        lowestGssFile.write('==================================================================\n')
#        for i in xrange(len(newStructs)):
#            lowestGssFile.write('\n***************** ' + atoms[i] + ' *******************\n')
#            atomLength = len(newStructs[i])
#            if atomLength >= 100:
#                for j in xrange(len(newStructs[i])):
#                    if (j + 1) % 20 == 0 or j == 99:
#                        lowestGssFile.write(str(newStructs[i][j]) + '\n')
#                    else:
#                        lowestGssFile.write(str(newStructs[i][j]) + ', ')
#        
#        lowestGssFile.flush()
#        os.fsync(lowestGssFile.fileno())
#    except IOError:
#        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to lowest_gss file. ~~~~~~~~~~\n'])

#def writeLowestVasp(lowestStructsFile, vStructsFinished, iteration, atoms):
#    try:
#        lowestStructsFile.write('==============================================================\n')
#        lowestStructsFile.write('\tIteration: ' + str(iteration) + '\n')
#        lowestStructsFile.write('==============================================================\n')
#        for i in xrange(len(vStructsFinished)):
#            lowestStructsFile.write('\n******************** ' + atoms[i] + ' ********************\n')
#            atomLength = len(vStructsFinished[i])
#            if atomLength >= 100:
#                for j in xrange(len(vStructsFinished[i][:100])):
#                    if (j + 1) % 20 == 0 or j == 99:
#                        lowestStructsFile.write(str(vStructsFinished[i][j]) + '\n')
#                    else:
#                        lowestStructsFile.write(str(vStructsFinished[i][j]) + ', ')
#            elif atomLength == 0:
#                lowestStructsFile.write('\nNo structures submitted.\n')
#            else:
#                for j in xrange(len(vStructsFinished[i][:atomLength])):
#                    if (j + 1) % 20 == 0 or j == atomLength - 1:
#                        lowestStructsFile.write(str(vStructsFinished[i][j]) + '\n')
#                    else:
#                        lowestStructsFile.write(str(vStructsFinished[i][j]) + ', ')                       
#        lowestStructsFile.flush()
#        os.fsync(lowestStructsFile.fileno())
#        
#    except IOError:
#        subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to lowest_vasp file. ~~~~~~~~~~\n'])
                


# -------------------------------------------- MAIN -----------------------------------------------
          
if __name__ == '__main__':
    maindir = '/fslhome/bch/cluster_expansion/graphene/testtm2'
#    maindir = os.getcwd()
    print 'Starting in ', maindir
    os.chdir(maindir)

    
    seed()
    maxvstructs = 20000 #maximum number of structures for each atom
    [atoms, volRange, clusterNums, startFromExisting, runTypes, PriorOrIID, ntrainStructs, mfitStructs, nfitSubsets, priorNum, plotTitle, xlabel, ylabel] = readSettingsFile()
    uncleOutput = open('uncle_output.txt','w') # All output from UNCLE will be written to this file.
    natoms = len(atoms)
    vstructsAll = [[]]*natoms #every structure vasp has attempted before this iteration, a list for each atom
    vstructsFinished = [[]]*natoms #every structure vasp has finished before this iteration, a list for each atom
    vstructsFailed = [[]]*natoms #every structure vasp has failed before this iteration, a list for each atom
    vstructsCurrent = [[]]*natoms #the structures vasp has attempted or will attempt this iteration, a list for each atom   
    newlyFinished = [[]]*natoms
    newlyFailed = [[]]*natoms
    newStructsPrior = [[]]*natoms  #structures from previous iteration with highest priority for vasp calculation
    holdoutStructs = [[]]*natoms
    #in vdata,unlike the other arrays, the struct field needs to be an integer, so I can count how many finished structures there are by count_nonzero
    vdata = zeros((natoms,maxvstructs),dtype = [('struct', int32),('conc', float), ('energy', float), ('natoms', int),('FE', float),('BE', float),('HFE', float)]) #data from vasp


    if not os.path.isdir('single_atoms'):
        manager1 = JobManager.JobManager(atoms)
        manager1.runSingleAtoms()
    if not os.path.isdir('hex_monolayer_refs'):
        manager1 = JobManager.JobManager(atoms)
        manager1.runHexMono()
    
    
    enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, ntrainStructs, uncleOutput)
    enumerator.enumerate()
    ntot = enumerator.getNtot(os.getcwd()+'/enum') #number of all enumerated structures
    energiesLast = zeros((natoms,ntot),dtype=float) #energies of last iteration, sorted by structure name

    createEnumPastDir(atoms)
    
    converged = False
    iteration = 1

    while not converged:
        changed = False
        
        subprocess.call(['echo','\n========================================================'])
        subprocess.call(['echo','\t\tIteration ' + str(iteration)])
        subprocess.call(['echo','========================================================\n'])
        #since atoms may have changed, have to reinitialize this
        enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, ntrainStructs, uncleOutput) 
        # Extract the pseudo-POSCARs from struct_enum.out
        extractor = Extractor.Extractor(atoms, uncleOutput, startFromExisting)
        if iteration == 1:
            vstructsFinished,vdata = parseStartStructures(atoms,startFromExisting,vstructsFinished,vdata,)
            vstructsCurrent = enumerator.chooseTrainingStructures(iteration,startFromExisting)
            vstructsCurrent = extractor.checkPureInCurrent(iteration,vstructsCurrent)
            finalDir = extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsCurrent)
        elif iteration > 1 and PriorOrIID == 'p':
            vstructsCurrent = newStructsPrior  #from previous iteration 
#            extractor.setStructsFromGSS(vstructsCurrent)
            finalDir = extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsCurrent)
        elif iteration > 1 and PriorOrIID == 'i':
            vstructsCurrent = enumerator.chooseTrainingStructures(iteration,startFromExisting)
            finalDir = extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsCurrent)
#        else:
#            finalDir = '' #starting exclusivelmy from structures.start
        

        # Create structures.in and structures.holdout files for each atom.
        os.chdir(maindir) #FIX this...shouldn't need it.
 
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atoms, startFromExisting, iteration, finalDir) 
        [newlyFinished, newlyFailed, vdata] = uncleFileMaker.makeUncleFiles(iteration, holdoutStructs,vstructsCurrent,vdata) 
        #update the vstructs lists and past_structs files       
        vstructsFinished = joinLists(vstructsFinished,newlyFinished)
        vstructsFailed = joinLists(vstructsFailed,newlyFailed)
        vstructsAll = joinLists(vstructsFinished,vstructsFailed)
        os.chdir(maindir) #FIX this...shouldn't need it.  
        pastStructsUpdate(vstructsCurrent,atoms)
        
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
        gss = GSS.GSS(atoms, volRange, plotTitle, xlabel, ylabel, uncleOutput)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch(iteration)
        gss.makePlots(iteration)
        #get the priority of each structure in each atom
        priorities = gss.getGssInfo(iteration,vstructsFailed) #first structure listed is highest priority
        #test weighted energy difference since last iteration
        diffe = getDiffE(priorities,energiesLast,atoms)
        energiesLast = priorities[:,:]['FE'] #to be used next iteration
        # Determine which atoms need to continue
        diffMax = 0.001 # 100 meV convergence criterion
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
            elif len(newlyFinished[i]) == 0 and len(vstructsCurrent[i] > 0):
                print 'Atom has only failed structures:', atom,'has been stopped.'
                atoms.remove(atom)
                atomnumbers.remove(iatom)
                rmAtoms.append(iatom)                
        natoms = len(atoms) #could be lower now
        #Keep only info on atoms that are continuing, so atom indices are right
#        del vstructsFinished[rmAtoms]
#        vstructsFinished2 = [[]]*natoms
#        vstructsFailed2 = [[]]*natoms
#        vstructsAll2 = [[]]*natoms
#        energiesLast2 = zeros((natoms,ntot),dtype=float) 
#        for i,ikeep in enumerate(atomnumbers):  
#            vstructsFinished2[i] = vstructsFinished[ikeep]
#            vstructsFailed2[i] = vstructsFailed[ikeep]
#            vstructsAll2[i] = vstructsAll2[ikeep]
#             #.tolist()           
#        vstructsFinished = vstructsFinished2
#        vstructsFailed = vstructsFailed2
#        energiesLast = energiesLast2

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


        newStructsPrior = getFromPriorities(priorities, priorNum, vstructsAll,natoms)

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
        

    
        
    
    
 
 
 
 
 
 
 
    