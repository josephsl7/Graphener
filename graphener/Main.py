'''
Created on Aug 26, 2014

@author: eswens13
'''

import os, subprocess,sys,re, time
from random import seed
from numpy import zeros,array,sqrt,std,amax,amin,int32,sort,count_nonzero,delete,mod
from copy import deepcopy

from comMethods import joinLists,structuresInWrite,writeLatticeVectors,readfile,writefile,\
                    convergeCheck,finishCheck,getNSW,getSteps

import Enumerator, Extractor, StructsToPoscar, JobManager, MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo     

def initializeStructs(atoms,restartTimeout):
    '''If structures.in exists in the atom folder, then parse that for starting structures.  If not, then leave lists empty for that atom'''
    lastDir = os.getcwd()
    natoms = len(atoms)
    subprocess.call(['echo','Checking structure folders for existing vasp data. . . \n'])
 
    nexistsStructsIn = 0
    nfoldersOK = 0
    #find starting methhod
    for iatom, atom in enumerate(atoms):
        nstruct = 0
        atomDir = lastDir + '/' + atom
        pureHdir =  atomDir + '/1'
        pureMdir =  atomDir + '/3'
#        os.chdir(atomDir); os.system('rm structures.in');os.chdir(lastDir) ######## remove me...used for testing
        if os.path.exists(atomDir + '/structures.in') and os.stat(atomDir + '/structures.in').st_size > 0: 
            nexistsStructsIn += 1
        if os.path.exists(pureHdir) and os.path.exists(pureMdir):
            for item in os.listdir(atomDir):
                itempath = atomDir + '/' + item
                if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
                    nstruct += 1
                    if nstruct == 3:
                        break #need at least 3 structures to make a fit 
        if nstruct == 3: nfoldersOK += 1
#        print 'atom nexistsStructsIn,nstruct ',nexistsStructsIn nstruct 
    if nexistsStructsIn == natoms: 
        startMethod = 'structures.in' 
        subprocess.call(['echo','Starting all from structures.in within atomic folders'])          
    elif nfoldersOK == natoms: 
        startMethod = 'struct folders'
        subprocess.call(['echo','Starting all from structure folders'])           
    elif 0 < nexistsStructsIn < natoms:
        subprocess.call(['echo','Some atoms folders have structures.in, and some do not.  They must all be alike'])
        sys.exit('Stop')
    elif 0 < nfoldersOK < natoms:
        subprocess.call(['echo','Some atoms folders have existing structure folders, and some do not.  They must all be alike to start with them.'])
        sys.exit('Stop')
    elif nexistsStructsIn == 0 and nfoldersOK == 0:
        startMethod = 'empty folders'
        subprocess.call(['echo','Starting all from iid structures'])
    else:
        subprocess.call(['echo','Error in choosing start method'])
        sys.exit('Stop')    
    #read data
    if startMethod == 'structures.in':
        [vstructsFinished, vstructsRestart, vstructsFailed, vdata] = parseStructsIn(atoms,restartTimeout)          
    elif startMethod == 'struct folders': 
           [vstructsFinished, vstructsRestart, vstructsFailed, vdata] = readInitialFolders(atoms,restartTimeout)           
    else:
        natoms = len(atoms)
        vstructsFinished = [[]]*natoms #every structure vasp has finished before this iteration, a list for each atom
        vstructsFailed = [[]]*natoms #every structure vasp has failed before this iteration, a list for each atom
        vstructsRestart = [[]]*natoms #those that timed out before starting or in previous iteration.  Must have folders in the atom directory.
        maxvstructs = 20000 #maximum number of structures for each atom
        #in vdata,unlike the other arrays, the struct field needs to be an integer, so I can count how many finished structures there are by count_nonzero
        vdata = zeros((natoms,maxvstructs),dtype = [('struct', int32),('conc', float), ('energy', float), ('natoms', int),('FE', float),('BE', float),('HFE', float)]) #data from vasp
   
    return  vstructsFinished, vstructsRestart, vstructsFailed, startMethod, vdata               
    os.chdir(lastDir)

        
def readInitialFolders(atoms,restartTimeout,):
    '''assigns all struct folders to either finished, restart, or failed.  Restart is optional.
    Initializes and writes structures.in and vdata.  Uncle formation energy is calculated (the other energies
    will be calculated in vdataToPlotFiles)
    
    Reasons not to restart (restart all others).  See slurmProblem():
    1. latest slurm file shows segfault or memory problem
    2. finished but not converged (too many ionic steps)  
    ''' 
    lastDir = os.getcwd()
    subprocess.call(['echo','Reading vasp data from structure folders . . . \n'])
    natoms = len(atoms)
    vstructsFinished = [[]]*natoms #every structure vasp has finished before this iteration, a list for each atom
    vstructsFailed = [[]]*natoms #every structure vasp has failed before this iteration, a list for each atom
    vstructsRestart = [[]]*natoms #those that timed out before starting or in previous iteration.  Must have folders in the atom directory.
    maxvstructs = 20000 #maximum number of structures for each atom
    #in vdata,unlike the other arrays, the struct field needs to be an integer, so I can count how many finished structures there are by count_nonzero
    vdata = zeros((natoms,maxvstructs),dtype = [('struct', int32),('conc', float), ('energy', float), ('natoms', int),('FE', float),('BE', float),('HFE', float)]) #data from vasp
    
    for iatom, atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom     
        #pure energies
        pureHdir =  atomDir + '/1'
        pureMdir =  atomDir + '/3'
        pureHenergy = float(readfile(pureHdir + '/OSZICAR')[-1].split()[2])/2.0 #2.0: per atom
        pureMenergy = float(readfile(pureMdir + '/OSZICAR')[-1].split()[2])/2.0 
        subprocess.call(['echo','\n\n{}'.format(atom)])
        subprocess.call(['echo','\tpure H system energy {}'.format(pureHenergy)]) 
        subprocess.call(['echo','\tpure M system energy {}'.format(pureMenergy)])     
        # all structures (including pure
        finishedStructs = []
        restartStructs = []
        failedStructs = []
        structlist = []
        ifinished = -1
        for item in os.listdir(atomDir):
            itempath = atomDir + '/' + item
            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
                if hasVaspFiles(itempath):
                    structlist.append(item)
                else:
                    subprocess.call(['rm','-r',itempath])  #delete folder     
        for istruct,struct in enumerate(structlist):
            failed = False
            vaspDir = atomDir + '/'+ str(struct)
            if mod(istruct+1,100) == 0 or istruct+1 ==len(structlist): 
                subprocess.call(['echo','\tChecked {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
            if os.stat(vaspDir + '/POSCAR').st_size > 0: 
                try:
                    atomCounts = setAtomCounts(vaspDir) #also changes any "new"-format POSCAR/CONTCAR back to old (removes the 6th line if it's text
                except:
                    failed = True
            if not failed and finishCheck(vaspDir) and not convergeCheck(vaspDir, getNSW(vaspDir)):
                failed = True 
            elif not failed and finishCheck(vaspDir) and convergeCheck(vaspDir, getNSW(vaspDir)): 
                finishedStructs.append(struct)
                ifinished += 1
                vdata[iatom,ifinished]['struct'] = struct
                natoms =  atomCounts[0] + atomCounts[1] 
                vdata[iatom,ifinished]['natoms'] = natoms
                #metal concentration
                conc = 0.0               
                if atomCounts[0] == 0:
                    conc = 1.0
                else:
                    conc = float(float(atomCounts[1])/natoms)
                vdata[iatom,ifinished]['conc'] = conc
                #energy and formation energy
                structEnergy = float(readfile(vaspDir + '/OSZICAR')[-1].split()[2])
                structEnergy = structEnergy/float(sum(atomCounts)) #per atom                     
                vdata[iatom,ifinished]['energy'] = structEnergy 
                formationEnergy = structEnergy - (conc * pureMenergy + (1.0 - conc) * pureHenergy)
                vdata[iatom,ifinished]['FE'] = formationEnergy
            elif restartTimeout and not failed and not slurmProblem(vaspDir):
                restartStructs.append(struct)  
            else:
                failedStructs.append(struct)
        vstructsFinished[iatom] =  finishedStructs 
        vstructsRestart[iatom] = restartStructs
        vstructsFailed[iatom] = failedStructs 
        subprocess.call(['echo','\tAtom {}: found {} finished, {} to restart, {} failed\n'.format(atom,len(finishedStructs),len(restartStructs),len(failedStructs))])
        #sort by FE
        nfinished = count_nonzero(vdata[iatom,:]['struct'])
        vdata[iatom,:nfinished] = sort(vdata[iatom,:nfinished],order = ['FE']) #must sort only the filled ones 
        
        #write structures.in, (or just the header if there is no existing data)
        structuresInWrite(atomDir,vdata[iatom,:nfinished]['struct'], vdata[iatom,:nfinished]['FE'],vdata[iatom,:nfinished]['conc'],vdata[iatom,:nfinished]['energy'],'w')
    os.chdir(lastDir)    
    return  vstructsFinished, vstructsRestart, vstructsFailed, vdata               

def parseStructsIn(atoms,vstructsFinished):
    """ Returns the structures, energies, concentrations from structures.start, 
    
    This method assumes that the structures.in file will have all of the pure
    structures and will have the following format for the structure identification line:
    graphene str #: (structure number) FE = 0.545, Concentration = .473
    or, for a pure structure:
    PURE M graphene str #: (structure number) FE = 0.0, Concentration = 1.0 """

    lastDir = os.getcwd()    
    natoms = len(atoms)
    vstructsFinished = [[]]*natoms #every structure vasp has finished before this iteration, a list for each atom
    vstructsFailed = [[]]*natoms #every structure vasp has failed before this iteration, a list for each atom
    vstructsRestart = [[]]*natoms #those that timed out before starting or in previous iteration.  Must have folders in the atom directory.
    maxvstructs = 20000 #maximum number of structures for each atom
    vdata = zeros((natoms,maxvstructs),dtype = [('struct', int32),('conc', float), ('energy', float), ('natoms', int),('FE', float),('BE', float),('HFE', float)]) #data from vasp
    for iatom,atom in enumerate(atoms):
        istruct = 0
        subList = []
        atomDir = lastDir + '/' + atom
        #parse structures.in
        startFile = atomDir + '/structures.in'
        subprocess.call(['cp',startFile,startFile + '_start']) #just for reference later
        lines = readfile(startFile)          
        for j in xrange(len(lines)):
            if list(lines[j].strip().split()[0])[:2] == ['#','-']:
                if j != len(lines) - 1:
                    structLine = lines[j + 1].strip().split()
                    if structLine[0].lower() == 'pure':
                        subList.append(structLine[5])
                        vdata[iatom,istruct]['struct'] = structLine[5]
                    else:
                        subList.append(structLine[3])
                        vdata[iatom,istruct]['struct'] = structLine[3]
                    natomsList = [int(i) for i in lines[j+6].strip().split()]
                    natoms = natomsList[0]+natomsList[1]
                    vdata[iatom,istruct]['natoms'] = natoms
                    vdata[iatom,istruct]['conc'] = float(natomsList[1]/float(natoms))
                    vdata[iatom,istruct]['energy'] = float(lines[j+9+natoms].strip().split()[0])  
                    istruct += 1                        
        vstructsFinished[iatom] = subList
        #get unfinished and failed
        structlist = []
        restartStructs = []
        failedStructs = []
        for item in os.listdir(atomDir):
            itempath = atomDir + '/' + item
            if item not in vstructsFinished[iatom] and os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
                if hasVaspFiles(itempath):
                    structlist.append(item)
                else:
                    subprocess.call(['rm','-r',itempath]) #delete folder     
        for istruct,struct in enumerate(structlist): #only the ones not in structures.in
            failed = False
            vaspDir = atomDir + '/'+ str(struct)
            if mod(istruct+1,100) == 0 or istruct+1 ==len(structlist): 
                subprocess.call(['printf','\tChecked {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
            
            if os.stat(vaspDir + '/POSCAR').st_size > 0: 
                try:
                    atomCounts = setAtomCounts(vaspDir) #also changes any "new"-format POSCAR/CONTCAR back to old (removes the 6th line if it's text
                except:
                    failed = True
            if not failed and finishCheck(vaspDir) and not convergeCheck(vaspDir, getNSW(vaspDir)):
                failed = True 
            elif not failed and finishCheck(vaspDir) and convergeCheck(vaspDir, getNSW(vaspDir)): 
                subprocess.call(['echo','\tAtom {}: found finished structure {} not in structures.in. Appending to restart list.'.format(atom,struct)])                    
                restartStructs.append(struct) 
            elif restartTimeout and not failed and not slurmProblem(vaspDir):
                restartStructs.append(struct)  
            else:
                failedStructs.append(struct)
        vstructsRestart[iatom] = restartStructs
        vstructsFailed[iatom] = failedStructs 
        subprocess.call(['echo','\n\tAtom {}: found {} finished, {} to restart, {} failed after reading structures.in'.format(atom,len(vstructsFinished[iatom]),len(restartStructs),len(failedStructs))])
        nfinished = count_nonzero(vdata[iatom,:]['struct'])
        vdata[iatom,:nfinished] = sort(vdata[iatom,:nfinished],order = ['FE']) #must sort only the filled ones  
    os.chdir(lastDir)
    return  vstructsFinished, vstructsRestart, vstructsFailed, vdata               

def hasVaspFiles(dir):
    return os.path.exists(dir+'/KPOINTS') and os.path.exists(dir+'/INCAR') and os.path.exists(dir+'/POSCAR') and os.path.exists(dir+'/POTCAR')

    
def checkStructsIn(list1,list2):
    '''makes sure that items in list1 are in list 2''' 
    for iatom, atom in enumerate(atoms):
        for item in list1[iatom]:
            if item not in list2[iatom]: 
                list2[iatom].append(item)
    return list2

def contains(struct, alist):
    if len(alist) == 0:
        return False
    for i in xrange(len(alist)):
        if str(struct) == str(alist[i]):
            return True   
    return False

def contcarToPoscar(structs,atoms,iteration):
    '''Prepares structures for restart, copying CONTCAR to POSCAR'''
    lastDir = os.getcwd()
    for iatom, atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom
        for struct in structs[iatom]:
            structDir = atomDir + '/' + str(struct)
            os.chdir(structDir)
            subprocess.call(['cp','POSCAR','POSCAR_{}'.format(iteration)] )
            if os.path.exists(structDir + '/CONTCAR'):
                if os.stat(structDir + '/CONTCAR').st_size > 0:
                    subprocess.call(['cp',structDir + '/CONTCAR',structDir + '/POSCAR'])    
    os.chdir(lastDir) 

def createEnumPastDir(atoms):
    '''enumpast is a folder for storing past_structs.dat, and enumerating new iid structures 
    that need past_structs.dat. This creates/clears the folder and puts an empty past_structs.dat there'''
    for i,atom in enumerate(atoms):
        atomDir = os.getcwd() + '/' + atom
        epDir = atomDir + '/enumpast'
        if os.path.isdir(epDir): subprocess.call(['rm','-r' ,epDir])                    
        subprocess.call(['mkdir', epDir])                
        file = open(epDir + '/past_structs.dat', 'w'); file.close()  #just create it. 

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

def extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsToStart,vstructsRestart):
    ''' Convert the extracted pseudo-POSCARs to VASP POSCAR files, make directories for them
     and put the POSCARs in their corresponding directories. Run VASP'''
#    vstructsToStart = extractor.getStructList()

    extractor.extract(vstructsToStart)
    if len(vstructsToStart)>0:
        subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
#        print 'BLOCKING toPoscar, .convert'
        toPoscar = StructsToPoscar.structsToPoscar(atoms, vstructsToStart)
        toPoscar.convertOutputsToPoscar()
    # Start VASP jobs and wait until they all complete or time out.
            # Start VASP jobs and wait until they all complete or time out.
        manager2 = JobManager.JobManager(atoms,ediffg)
        if runTypes ==['low']:
#            subprocess.call(['echo','Warning: BLOCKING RUN MANAGER for testing' ])
            manager2.runLowJobs(vstructsToStart,vstructsRestart)
        elif runTypes ==['low','normal']:
            manager2.runLowJobs(vstructsToStart,vstructsRestart)
            manager2.runNormalJobs(vstructsToStart,vstructsRestart)          
        elif runTypes ==['low','normal','dos']:
            manager2.runLowJobs(vstructsToStart,vstructsRestart)
            manager2.runNormalJobs(vstructsToStart,vstructsRestart)
            manager2.runDOSJobs(vstructsToStart,vstructsRestart) 
        else:
            subprocess.call(['echo','Your RUN_TYPES is {}'.format(runTypes)]) 
            sys.exit('The only supported RUN_TYPES are "low", "low normal" and "low normal DOS"') 
    if runTypes ==['low']:
        finalDir = '/'
    elif runTypes ==['low','normal']:
        finalDir = '/normal'          
    elif runTypes ==['low','normal','dos']:
        finalDir = '/DOS' 
    else:
        subprocess.call(['echo','Your RUN_TYPES is {}'.format(runTypes)]) 
        sys.exit('The only supported RUN_TYPES are "low", "low normal" and "low normal DOS"')   
    return finalDir

def getFromPriorities(priorities, N, vstructsAll, atoms):
    '''take N structures new to vasp of highest priority'''
    newstructs = [[]]*len(atoms)
    for iatom,atom in enumerate(atoms):
        nNew = 0 
        istruct = 0 
        sublist = []
        while nNew < N:
            struct = str(priorities[iatom,istruct]['struct'])
            if not struct in vstructsAll[iatom]:
                sublist.append(struct)
                nNew += 1
            istruct += 1
        newstructs[iatom] = sublist
    return newstructs
            
def getDiffE(priorities, energiesLast, atoms):
    '''Finds the L1 norm energy change between iterations, weighted by priority
       Must have priority sorted by structure name (done before called), and 
       energiesLast the same'''

    ediffL1 = zeros(len(atoms))
    for iatom in range(len(atoms)):
        ediffw = 0.0
        for istruct in range(len(priorities[iatom,:]['struct'])):#can't just enumerate because structures start at 1, not 0.
            ediff  = abs(priorities[iatom,istruct]['FE'] - energiesLast[iatom,istruct])
#            if ediff>1e-8:print'ediff', istruct, priorities[iatom,istruct]['struct'], ediff
            ediffw += ediff * priorities[iatom,istruct]['prior'] #weighted
        priorsum = sum(priorities[iatom,:]['prior'])
        ediffL1[iatom] = sqrt(ediff/priorsum)
#        ediffL1[iatom] = sqrt(sum(abs(ediff) * priorities[iatom,:]['prior']/priorsum))
    return ediffL1

def getNSW(dir): 
    proc = subprocess.Popen(['grep','-i','NSW',dir+'/INCAR'],stdout=subprocess.PIPE) 
    return int(proc.communicate()[0].split('=')[-1])


def multiDelete(list_, args):
    indexes = sorted(args, reverse=True)
    for index in indexes:
        del list_[index]
    return list_
 
def pastStructsUpdate(structs,atoms):
    lastDir = os.getcwd()
    for iatom,atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom
        pastStructs  = [line.strip() for line in readfile(atom +'/enumpast/past_structs.dat')]
        for struct in structs[iatom]:
            if struct not in pastStructs: pastStructs.append(struct)
        file = open( atomDir +'/enumpast/past_structs.dat','w')
        for struct in pastStructs:
            file.write(str(struct)+'\n')
        file.close       

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
    runTypes = ['low']
    PriorOrIID = 'p'
    niid = 0
    mfitStructs = 0
    nfitSubsets = 0
    priorNum = 0
    plotTitle = "title"
    xlabel = "xlabel"
    ylabel = "ylabel"
    restartTimeout = False 
    
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

#        elif line.split()[0] == 'START_FROM_EXISTING:': #read from a structures.start and fit these.  Needs to come after the ATOMS list in the settings lin
#            if line.split()[1][0].lower() == 'n': #"No" or "None" means no atoms will start from existing calcs.
#               startMethod = [False]*len(atoms) 
#            elif line.split()[1][0].lower() == 'a': #"All" 
#               startMethod = [True]*len(atoms)
#            else:
#                try:  
#                    parselist = re.search(':(.+?)#', line).group(1).lower().split()[:len(atoms)] #take only enough to match number of atoms, as it might be an old list
#                    startMethod = [True if i == 't' else False for i in parselist]
#                    if len(startMethod) < len(atoms): 
#                        subprocess.call(['echo','START_FROM_EXISTING length is too short for the number of atoms! Stopping'])
#                        sys.exit('Stop')
#                except AttributeError:
#                    print "Can't parse START_STRUCTS T's and F's.  Defaulting to no start structs." 
#                    startMethod = [False]*len(atoms)               
        elif line.split()[0] == 'PRIORITY/IID:':
            if str(line.split()[1]).lower() == 'i':
                PriorOrIID = 'i'

        elif line.split()[0] == 'RUN_TYPES:':
            try:
                runTypes = re.search(':(.+?)#', line).group(1).lower().split() #extracts only text between START_STRUCTS and #
            except AttributeError:
                subprocess.call(['echo',"Can't parse RUN_TYPES. Defaulting to LOW precision runs only."])            

        elif line.split()[0] == 'N_IID:':
            niid = int(line.split()[1])
          
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
            
        elif line.split()[0] == 'RESTART_TIMEOUT:':
            if line.split()[1][0].lower() == 'y': 
               restartTimeout = True
            
        elif line.split()[0] == 'EDIFFG:': #ionic relaxation energy max diff, for INCAR, per atom.  May also be used or scaled for ediff, electronic relax energy max diff
            ediffg = float(line.split()[1])
    
    return [atoms, volRange, clusterNums, runTypes, PriorOrIID, niid, mfitStructs, nfitSubsets, priorNum, plotTitle, xlabel, ylabel,restartTimeout,ediffg]

def removeStructs(list1,list2):
    '''Remove items from list1 that might be in list2, and return list2'''
    for iatom, atom in enumerate(atoms):
        for item in list1[iatom]:
            if item in list2[iatom]: list2[iatom].remove(item)
    return list2
#
#def timeoutCheck(dir):
#    '''Checks the latest slurm file for the words TIME LIMIT. 
#       Or if it doesn't have a slurm, it means it neverstarted, so start those'''
#    slurmlist = []
#    structfiles = os.listdir(dir)
#    failed = False
#    for file in structfiles:
#        filePath = dir + '/' + file
#        if 'slurm' in file and os.stat(filePath).st_size > 0:
#            slurmlist.append(file)
#    if len(slurmlist)==0: 
#        return True
#    else: #only use the last slurm
#        slurmlist.sort()
#        lastslurm = slurmlist[-1]
#        for line in readfile(dir + '/' + lastslurm):
#            if 'TIME LIMIT' in line: 
#                return True
#        return False

def slurmProblem(dir):
    '''Checks all slurm files for key strings that indicate we shouldn't restart. 
       Or if it doesn't have a slurm, it means it neverstarted, so start those'''
    slurmlist = []
    problems = ['memory','sigsev','segmentation']
    structfiles = os.listdir(dir)
    exceed = False
    for file in structfiles:
        filePath = dir + '/' + file
        if 'slurm' in file and os.stat(filePath).st_size > 0:
            slurmlist.append(file)
    for slurm in slurmlist:
        for line in readfile(dir + '/' + slurm):
            for problem in problems:
                if problem in line: 
                    return True        
    return False                         

def setAtomCounts(poscarDir):
    """ Retrieves the number of H atoms and the number of M atoms from the POSCAR file and sets 
        the corresponding members. 
        Also fixes "new" POSCAR/CONTCAR format (comes from CONTCAR) back to old for UNCLE use (removes the 6th line if it's text """
    atomCounts = []
    fixPOSCAR = False
    poscarLines = readfile(poscarDir + '/POSCAR')
    counts = poscarLines[5].strip().split() 
    if not counts[0][0].isdigit(): #have the "new" POSCAR format that gives the atom types in text on line 6 (5-python)
        fixPOSCAR = True
        counts = poscarLines[6].strip().split()  
    if len(counts) == 3:
        atomCounts.append(int(counts[1]))
        atomCounts.append(int(counts[2]))
    elif len(counts) == 2:
        if poscarLines[0].split()[1] == 'H':
            atomCounts.append(int(counts[1]))
            atomCounts.append(0)
        elif poscarLines[0].split()[1] == 'M':
            atomCounts.append(0)
            atomCounts.append(int(counts[1]))
    natoms = sum(atomCounts)
    if fixPOSCAR:
        del(poscarLines[5])
        writefile(poscarLines[:7+natoms*2],poscarDir + '/POSCAR') #:7+natoms*2 is because CONTCAR includes velocity lines that uncle doesn't want. The factor of 2 is because carbon atoms are not included in natoms      
    return atomCounts

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


# -------------------------------------------- MAIN -----------------------------------------------
'''All atom folders must start with the same method:  
1) structures.in from previous run must exist in each atom folder
or
2) structures.in must be deleted from each atomfolder. Structure folders must exist from previous runs (including pure cases) in atom folders, which can be a mix of finished and unfinished structures.  
s from previous run must exist
or
3) if all the above are missing from all folders, each atom will start from scratch with iid structures'''          
if __name__ == '__main__':
#    maindir = '/fslhome/bch/cluster_expansion/graphene/testtm3'  
#    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1'
    maindir = os.getcwd()

    subprocess.call(['echo','Starting in ' + maindir])
    
    os.chdir(maindir)
    pathMax = maindir + '/needed_files/diffMax'
    if os.path.exists(pathMax):
        diffMax = float(readfile(pathMax)[0].strip()) 
        subprocess.call(['echo','\ndiffMax read from file (convergence criterium): '+str(diffMax) + '\n'])
    else:
        subprocess.call(['echo','\nMissing diffMax file with the convergence criterium diffMax in eV. Stopping'])
        sys.exit('Stop')
    seed()

    [atoms, volRange, clusterNums, runTypes, PriorOrIID, niid, mfitStructs, nfitSubsets, priorNum, plotTitle, xlabel, ylabel,restartTimeout,ediffg] = readSettingsFile()
    uncleOutput = open('uncle_output.txt','w') # All output from UNCLE will be written to this file.
    natoms = len(atoms)
    vstructsAll = [[]]*natoms #every structure vasp has attempted before this iteration, a list for each atom
    vstructsToStart = [[]]*natoms #the structures vasp must create folders for and start this iteration, a list for each atom   
    vstructsRestart = [[]]*natoms
    newlyFinished = [[]]*natoms
    newlyFailed = [[]]*natoms
    newlyToRestart = [[]]*natoms
    newStructsPrior = [[]]*natoms  #structures from previous iteration with highest priority for vasp calculation
    holdoutStructs = [[]]*natoms

    if not os.path.isdir('single_atoms'):
        manager1 = JobManager.JobManager(atoms,ediffg)
        manager1.runSingleAtoms()
    if not os.path.isdir('hex_monolayer_refs'):
        manager1 = JobManager.JobManager(atoms,ediffg)
        manager1.runHexMono()
    #assign all existing struct folders to either finished, restart, or failed
    [vstructsFinished,vstructsRestart0,vstructsFailed,startMethod,vdata] = initializeStructs(atoms,restartTimeout)    
    vstructsAll = joinLists(vstructsFinished,vstructsFailed)
    vstructsAll = joinLists(vstructsAll,vstructsRestart0)  
  
    enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, niid, uncleOutput)
    subprocess.call(['echo','Warning: BLOCKING ENUMERATOR to save time' ])
#    enumerator.enumerate() #comment this out for testing
    ntot = enumerator.getNtot(os.getcwd()+'/enum') #number of all enumerated structures
    energiesLast = zeros((natoms,ntot),dtype=float) #energies of last iteration, sorted by structure name

    createEnumPastDir(atoms)
    pastStructsUpdate(vstructsAll,atoms)  
   
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
        extractor = Extractor.Extractor(atoms, uncleOutput, startMethod)           
        if iteration == 2: #bring in restarts from initial folders
            vstructsRestart = joinLists(vstructsRestart0,vstructsRestart)
        if iteration == 1: 
            if startMethod == 'empty folders': 
                vstructsToStart = enumerator.chooseTrainingStructures(iteration,startMethod)
                vstructsToStart = extractor.checkPureInCurrent(iteration,vstructsToStart,vstructsFinished)
            vstructsToRun = vstructsToStart #no restarts in first iteration
        elif iteration > 1 and PriorOrIID == 'p':
            vstructsToStart = newStructsPrior  #from previous iteration 
            contcarToPoscar(vstructsRestart,atoms,iteration) 
            vstructsToRun = joinLists(vstructsRestart,vstructsToStart)
        elif iteration > 1 and PriorOrIID == 'i':
            vstructsToStart = enumerator.chooseTrainingStructures(iteration,startMethod)   
            contcarToPoscar(vstructsRestart,atoms,iteration)
            vstructsToRun = joinLists(vstructsRestart,vstructsToStart)
 
        pastStructsUpdate(vstructsToStart,atoms)
        finalDir = extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsToStart,vstructsRestart)           
        os.chdir(maindir) #FIX this...shouldn't need it.
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atoms, startMethod, iteration, finalDir, restartTimeout) 
        [newlyFinished, newlyToRestart, newlyFailed,vdata] = uncleFileMaker.makeUncleFiles(iteration, holdoutStructs,vstructsToRun,vdata) 
#        print 'newlyFinished',newlyFinished
        vstructsFinished = joinLists(vstructsFinished,newlyFinished)
        vstructsFailed = joinLists(vstructsFailed,newlyFailed)
        vstructsAll = joinLists(vstructsFinished,vstructsFailed)
        os.chdir(maindir) #FIX this...shouldn't need it.  


#        for iatom, atom in enumerate(atoms):
#            print atom, 'vstructsFinished in Main1', vstructsFinished[iatom]
        
        # Get all the structs that have been through VASP calculations for each atom. These
        # should be sorted by formation energy during the work done by makeUncleFiles()
                
        # Perform a fit to the VASP data in structures.in for each atom.
        fitter = Fitter.Fitter(atoms, mfitStructs, nfitSubsets, vstructsFinished,uncleOutput)
        if iteration == 1: 
            fitter.makeFitDirectories()
            if startMethod != 'empty folders': fitter.holdoutFromIn(atoms)
        fitter.fitVASPData(iteration)
    
        # Perform a ground state search on the fit for each atom.    
        gss = GSS.GSS(atoms, volRange, plotTitle, xlabel, ylabel, vstructsFinished,uncleOutput)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch(iteration)
        gss.makePlots(iteration)
        #get the priority of each structure in each atom
        priorities = gss.getGssInfo(iteration,vstructsFailed) #first structure listed is highest priority
        #choose new structures while still sorted by priority
        newStructsPrior = getFromPriorities(priorities, priorNum, vstructsAll,atoms)

        #test weighted energy difference since last iteration
        # First sort priorities by structure name so we have an unchanging order
        for iatom, atom in enumerate(atoms): 
            priorities[iatom,:] = sort(priorities[iatom,:],order = ['struct'])

        for iatom, atom in enumerate(atoms):       
            print 'priorities for atom, iteration',atom, iteration
            print priorities[iatom,:10]['struct']
            print priorities[iatom,:10]['FE']
            print 'energiesLast'
            print energiesLast[iatom,:10]
        diffe = getDiffE(priorities,energiesLast,atoms)
        energiesLast = priorities['FE'] #to be used next iteration, sorted by structure
#        print 'new energiesLast', energiesLast
        
        # Determine which atoms need to continue
#        diffMax = 0.01 # 10 meV convergence criterion
        diffMax = float(readfile(maindir + '/needed_files/diffMax')[0].strip()) #read from file called 'diffMax', so we can experiment with it during a long run
        subprocess.call(['echo','Weighted energy changes']) 
        atomnumbers = range(natoms)
        rmAtoms = []
        for iatom, atom in enumerate(atoms):
            subprocess.call(['echo','\t{} {:8.6f} eV'.format(atom,diffe[iatom])])
            if diffe[iatom]<diffMax: #Converged
                subprocess.call(['echo','Atom {} has converged'.format(atom)])
                atoms.remove(atom)
                atomnumbers.remove(iatom)
                rmAtoms.append(iatom)
            elif len(newlyFinished[iatom]) + len(newlyToRestart[iatom]) == 0 and len(vstructsToStart[iatom]) > 0:
                subprocess.call(['echo','Atom {} did not finish any new structures, and has been stopped.:'.format(atom)])
                atoms.remove(atom)
                atomnumbers.remove(iatom)
                rmAtoms.append(iatom)                
        natoms = len(atoms) #could be lower now
        if len(rmAtoms)>0:
            vstructsFinished = multiDelete(vstructsFinished,rmAtoms)
            vstructsFailed = multiDelete(vstructsFailed,rmAtoms)
            newlyToRestart = multiDelete(newlyToRestart,rmAtoms)
            vstructsAll = multiDelete(vstructsAll,rmAtoms)
            vdata = delete(vdata,rmAtoms,axis=0)        
            energiesLast = delete(priorities,rmAtoms,axis=0)['FE'] #    priorities[ikeep,:]['FE']
            newStructsPrior = multiDelete(newStructsPrior,rmAtoms)
        if natoms == 0:
            subprocess.call(['echo','\n----------------- All atoms have finished ---------------'])
            converged = True
        else:       
            if natoms >1:
                subprocess.call(['echo','At end of iteration {} are continuing:'.format(' '.join(atoms))])
            else:
                subprocess.call(['echo','At end of iteration {} is continuing:'.format(' '.join(atoms))])
        #---------- prep for next iteration --------------
        vstructsRestart = newlyToRestart

        # Set holdoutStructs for the next iteration to the 100 (or less) lowest vasp structures
         
        for iatom in range(natoms):
            nfinished = count_nonzero(vdata[iatom,:100]['struct'])
            if  nfinished > 100:
                holdoutStructs[iatom] = vdata[iatom,:100]['struct'].tolist()
            else:
                holdoutStructs[iatom] = vdata[iatom,:nfinished]['struct'].tolist() 

        # Go to next iteration
        iteration += 1
        
    uncleOutput.close()


    subprocess.call(['echo','\n---------- PROGRAM ENDED ----------\n'])