'''
'''

import os, subprocess,sys,re,time
from random import seed
from numpy import zeros,array,sqrt,std,amax,amin,int32,sort,count_nonzero,delete,mod
from numpy import dot,rint,fromfile
from numpy.linalg import inv,det
from copy import deepcopy

from comMethods import *

import Enumerator, Extractor, StructsToPoscar, JobManager, MakeUncleFiles, Fitter, GSS, \
        Analyzer, MovementInfo, PlotStructures 

from PlotStructures import plotStructsByPrior,collateStructsConc, collateStructsHFE     

def initializeStructs(atoms,restartTimeout,rmStructIn,pureMetal):
    ''' '''
    lastDir = os.getcwd()
    natoms = len(atoms)
    subprocess.call(['echo','Checking structure folders for existing vasp data. . . \n'])
    nexistsStructsIn = 0
    nfoldersOK = 0
    #find starting method
    for iatom, atom in enumerate(atoms):
        nstruct = 0
        atomDir = lastDir + '/' + atom
        pureHdir =  atomDir + '/1'
        pureMdir =  atomDir + '/' + pureMetal       
        if rmStructIn and os.path.exists(atomDir + '/structures.in'): 
            os.chdir(atomDir)
            os.system('mv structures.in structures.in.old')
            subprocess.call(['echo','Not using structures.in: moved to structures.in.old'])          
            os.chdir(lastDir) 
        if os.path.exists(atomDir + '/structures.in') and os.stat(atomDir + '/structures.in').st_size > 0: 
            nexistsStructsIn += 1
        if os.path.exists(pureHdir) and os.path.exists(pureMdir) and finishCheck(pureHdir) and finishCheck(pureMdir):
            for item in os.listdir(atomDir):
                itempath = atomDir + '/' + item
                if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
                    nstruct += 1
                    if nstruct == 3:
                        break #need at least 3 structures to make a fit 
        if nstruct == 3: nfoldersOK += 1

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
        subprocess.call(['echo','Some atoms folders have at least three existing structure folders, and some do not.  They must all be alike to start with them.'])
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
        vdata = zeros((natoms,maxvstructs),dtype = [('struct', int32),('conc', float), \
            ('energy', float), ('nadatoms', int), ('nCarbon', int),('FE', float),('BE', float),('HFE', float)]) #data from vasp
   
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
    vdata = zeros((natoms,maxvstructs),dtype = [('struct', int32),('conc', float), \
        ('energy', float), ('nadatoms', int), ('nCarbon', int),('FE', float),('BE', float),('HFE', float)]) #data from vasp
    for iatom, atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom     
        #pure energies
        pureHdir =  atomDir + '/1'
        pureMdir =  atomDir + '/' + pureMetal
        pureAtomCounts = setAtomCounts(pureHdir) #in this case gives the number of adatom sites
        pureHenergy = float(readfile(pureHdir + '/OSZICAR')[-1].split()[2])/float(pureAtomCounts[1])  
        pureMenergy = float(readfile(pureMdir + '/OSZICAR')[-1].split()[2])/float(pureAtomCounts[1])
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
            if mod(istruct+1,100) == 0 or istruct+1 ==len(structlist): subprocess.call(['echo','\tChecked {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
            if os.stat(vaspDir + '/POSCAR').st_size > 0: 
                try:
                    atomCounts = setAtomCounts(vaspDir) #also changes any "new"-format POSCAR/CONTCAR back to old (removes the 6th line if it's text
                except:
                    failed = True
            if not failed and finishCheck(vaspDir) and not convergeCheck(vaspDir, getNSW(vaspDir)):
                failed = True 
            if (not os.path.exists(vaspDir + '/OUTCAR')) or outcarWarn(vaspDir): 
                failedStructs.append(struct)
                subprocess.call(['echo','\tOUTCAR warning for struct {}: failed'.format(struct)])
            elif not failed and not energyDropCheck(vaspDir):
                subprocess.call(['echo','\tEnergy rose unphysically for struct {}: failed'.format(struct)])
                failed = True               
            elif not failed and finishCheck(vaspDir) and convergeCheck(vaspDir, getNSW(vaspDir)): 
                finishedStructs.append(struct)
                ifinished += 1
                vdata[iatom,ifinished]['struct'] = struct
                nCarbon =  atomCounts[0]
                nadatoms =  float(atomCounts[1] + atomCounts[2]) 
                vdata[iatom,ifinished]['nadatoms'] = nadatoms
                vdata[iatom,ifinished]['nCarbon'] = nCarbon
                #metal concentration
                conc = 0.0               
                if atomCounts[1] == 0:
                    conc = 1.0
                else:
                    conc = float(float(atomCounts[2])/nadatoms)#metal conc
                vdata[iatom,ifinished]['conc'] = conc
                #energy and formation energy
                structEnergy = float(readfile(vaspDir + '/OSZICAR')[-1].split()[2])
                structEnergy = structEnergy/float(nadatoms) #per adatom                     
                vdata[iatom,ifinished]['energy'] = structEnergy 
                formationEnergy = structEnergy - (conc * pureMenergy + (1.0 - conc) * pureHenergy)
                vdata[iatom,ifinished]['FE'] = formationEnergy
            elif not failed and restartTimeout and not slurmProblem(vaspDir):
                restartStructs.append(struct)  
                subprocess.call(['echo','\tStruct {}: restart'.format(struct)])
            else:
                failedStructs.append(struct)
                subprocess.call(['echo','\tStruct {}: failed'.format(struct)])

        vstructsFinished[iatom] =  finishedStructs 
        vstructsRestart[iatom] = restartStructs
        vstructsFailed[iatom] = failedStructs 
        subprocess.call(['echo','\tAtom {}: found {} finished, {} to restart, {} failed\n'.format(atom,len(finishedStructs),len(restartStructs),len(failedStructs))])
        #sort by FE
        nfinished = count_nonzero(vdata[iatom,:]['struct'])
        vdata[iatom,:nfinished] = sort(vdata[iatom,:nfinished],order = ['FE']) #must sort only the filled ones 
        
        #write structures.in, (or just the header if there is no existing data)
        structuresWrite('all',atomDir,vdata[iatom,:nfinished]['struct'], vdata[iatom,:nfinished]['FE'],\
                        vdata[iatom,:nfinished]['conc'],vdata[iatom,:nfinished]['energy'],'.in','w')
    os.chdir(lastDir)    
    return  vstructsFinished, vstructsRestart, vstructsFailed, vdata               

def parseStructsIn(atoms,vstructsFinished):
    """ Returns the structures, energies, concentrations from structures.in, 
    
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
    vdata = zeros((natoms,maxvstructs),dtype = [('struct', int32),('conc', float), \
        ('energy', float), ('nadatoms', int), ('nCarbon', int),('FE', float),('BE', float),('HFE', float)]) #data from vasp
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
                    vectorLines = [line.strip().split() for line in lines[j+3:j+6]]                      
                    lattVec1 = [float(vectorLines[0][0]), float(vectorLines[0][1]), float(vectorLines[0][2])]
                    lattVec2 = [float(vectorLines[1][0]), float(vectorLines[1][1]), float(vectorLines[1][2])]
                    lattVec3 = [float(vectorLines[2][0]), float(vectorLines[2][1]), float(vectorLines[2][2])]
                    LV = zeros((3,3),dtype =float)
                    LV[:,0] = lattVec1
                    LV[:,1] = lattVec2
                    LV[:,2] = lattVec3 
                    cellVol = det(LV)
                    PLV = array(  [[0.00000000,    2.13128850,   -1.23050000], #primitive lattice vectors
                                   [0.00000000,    2.13128850,    1.23050000], 
                                   [1000.00000,    0.00000000,    0.00000000]])
                    primCellVol = det(PLV)
                    dxC = PLV[0,0]*2*0.333333333333333 #distance between 2 C atoms
                    nCarbon = 2 * int(rint(cellVol/primCellVol))       
                    structLine = lines[j + 1].strip().split()
                    if structLine[0].lower() == 'pure':
                        subList.append(structLine[5])
                        vdata[iatom,istruct]['struct'] = structLine[5]
                    else:
                        subList.append(structLine[3])
                        vdata[iatom,istruct]['struct'] = structLine[3]
                    nadatomsList = [int(i) for i in lines[j+6].strip().split()]
                    nadatoms = nadatomsList[0]+nadatomsList[1]
                    vdata[iatom,istruct]['nadatoms'] = nadatoms 
                    vdata[iatom,istruct]['nCarbon'] = nCarbon
                    vdata[iatom,istruct]['conc'] = float(nadatomsList[1]/float(nadatoms))
                    vdata[iatom,istruct]['energy'] = float(lines[j+9+nadatoms].strip().split()[0])  
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
            if (not os.path.exists(vaspDir + '/OUTCAR')) or outcarWarn(vaspDir):  
                failedStructs.append(struct)
                subprocess.call(['echo','\tOUTCAR warning for struct {}: failed'.format(struct)])
#                restartStructs.append(struct)
#                subprocess.call(['echo','\tOUTCAR warning for struct {}: restart'.format(struct)])
            elif not failed and not energyDropCheck(vaspDir):
                subprocess.call(['echo','\tEnergy rose unphysically for struct {}: failed'.format(struct)])
                failed = True                
            elif not failed and finishCheck(vaspDir) and not convergeCheck(vaspDir, getNSW(vaspDir)):# Don't restart any that have finished but converged.
                failed = True 
            elif not failed and finishCheck(vaspDir) and convergeCheck(vaspDir, getNSW(vaspDir)): 
                subprocess.call(['echo','\tAtom {}: found finished structure {} not in structures.in. Appending to restart list.'.format(atom,struct)])                    
                restartStructs.append(struct) 
            elif not failed and restartTimeout and not slurmProblem(vaspDir):
                #begin cluge
                os.chdir(vaspDir)
#                subprocess.call(['echo', 'submitting restarts in parse...']) 
#                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
#                jobid = proc.communicate()[0].split()[3]
#                subprocess.call(['echo', 'Submitted job ' + jobid])
                #end cluge
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
        
def enumerationDone(dir,nTotClusters,nTotStructs):
    '''Checks for the output files of enumeration'''
    print nTotClusters,nTotStructs,nTotClusters*nTotStructs
    if os.path.exists(dir + '/enum_PI_matrix.out'):
        subprocess.call(['echo','\nChecking enum_PI_matrix.out for correct size. . .'])
        lines = readfile(dir + '/enum_PI_matrix.out')
        n = len(lines)
        del lines #huge list!
        if n == nTotStructs: # and == nTotClusters: #The number of columns (clusters) doesn't make sent to me yet.
            subprocess.call(['echo','\tOK\n'])
            if os.path.exists(dir + '/clusters.out') and os.path.exists(dir + '/struct_enum.out'): 
                return True
        else:
            subprocess.call(['echo','\tNeeds to be calculated\n'])
    else:
        return False

def existAllJ1out(atoms):
    '''Checks if all atoms have the file fits/J.1.out'''
    existAll = True
    for atom in atoms:
        atomDir = os.getcwd() + '/' + atom
        if not os.path.exists('{}/fits/J.1.out'.format(atomDir)):
            existAll = False
            break
    return existAll    
    

def extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsToStart,vstructsRestart):
    ''' Convert the extracted pseudo-POSCARs to VASP POSCAR files, make directories for them
     and put the POSCARs in their corresponding directories. Run VASP'''
#    vstructsToStart = extractor.getStructList()
   
    if len(flat(vstructsToStart))>0:
        extractor.extract(vstructsToStart) # Extract the pseudo-POSCARs from struct_enum.out
        subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
#        print 'BLOCKING toPoscar, .convert'
        toPoscar = StructsToPoscar.structsToPoscar(atoms, vstructsToStart)
        toPoscar.convertOutputsToPoscar()
    # Start VASP jobs and wait until they all complete or time out.       
    vstructsToRun = joinLists([vstructsToStart,vstructsRestart])
    if len(flat(vstructsToRun))>0:
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

def getFromPriorities(priorities, vstructsAll, atoms,nNew):
    '''take N structures new to vasp of highest priority'''
    newstructs = [[]]*len(atoms)
    for iatom,atom in enumerate(atoms):
        nNewStruct = 0 
        istruct = 0 
        sublist = []
        subprocess.call(['echo','\tAtom {}, nNew {}'.format(atom,nNew[iatom])])
        while nNewStruct < nNew[iatom]:
            struct = str(priorities[iatom,istruct]['struct'])
            if not struct in vstructsAll[iatom]:
                sublist.append(struct)
                nNewStruct += 1
            istruct += 1
        subprocess.call(['echo','\tAtom {}, \tlowest priority new struct: {}, {}'.format(atom,str(priorities[iatom,istruct-1]['struct']),str(priorities[iatom,istruct-1]['prior']))])
            
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

def getnNew(atoms,vstructsFinished,vstructsRestart,niid,nPrior,nTotStructs,PriorOrIID):
    '''For each atom, chooses how many new structures to use for the next iteration'''
    natoms = len(atoms)
    nNew = [0]*natoms
    for iatom, atom in enumerate(atoms):
        nFinished = len(vstructsFinished[iatom])
        nRestart = len(vstructsRestart[iatom])
        if PriorOrIID == 'i':
            nNew[iatom] = min(niid,nTotStructs-(nFinished+nRestart)) 
        else:
            nNew[iatom] = min(nPrior,nTotStructs-(nFinished+nRestart))
    return nNew

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
        if len(line.strip().split()) > 0: #in case a blank line is in there
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
    maxiid = 2e6 
    mfitStructs = 0
    nfitSubsets = 0
    nPrior = 0
    plotTitle = "title"
    xlabel = "xlabel"
    ylabel = "ylabel"
    restartTimeout = False 
    rmStructIn = False
    maxE = 100.0 #eV, max from vasp to include in fit
    graphsOnly = False
    maxIter = 100
    distribute = True 
    extendpath = '' # builds run off of previous (usually smaller) run with symbolic links
    
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

        elif line.split()[0] == 'MAX_IID:':
            maxiid = int(line.split()[1])
         
        elif line.split()[0] == 'M_FITTING_STRUCTS:':
            mfitStructs = int(line.split()[1])
            
        elif line.split()[0] == 'N_STRUCT_SUBSETS:':
            nfitSubsets = int(line.split()[1])
        
        elif line.split()[0] == 'PRIOR_NUM:':
            nPrior = int(line.split()[1])
        
        elif line.split()[0] == 'PLOT_TITLE:':
            plotTitle = line.split('\'')[1]
        
        elif line.split()[0] == 'XLAB:':
            xlabel = line.split('\'')[1]
        
        
        elif line.split()[0] == 'YLAB:':
            ylabel = line.split('\'')[1]
            
        elif line.split()[0] == 'RESTART_TIMEOUT:':
            if line.split()[1][0].lower() == 'y': 
               restartTimeout = True
               
        elif line.split()[0] == 'REMOVE_STR.IN:':
            if line.split()[1][0].lower() == 'y': 
               rmStructIn = True               
            
        elif line.split()[0] == 'EDIFFG:': #ionic relaxation energy max diff, for INCAR, per atom.  May also be used or scaled for ediff, electronic relax energy max diff
            ediffg = float(line.split()[1])
            
        elif line.split()[0] == 'MAX_E:': #ionic relaxation energy max diff, for INCAR, per atom.  May also be used or scaled for ediff, electronic relax energy max diff
            maxE = float(line.split()[1])

        elif line.split()[0] == 'GRAPHS_ONLY:':
            if line.split()[1][0].lower() == 'y': 
               graphsOnly = True
               
        elif line.split()[0] == 'PURE_METAL:':
            pureMetal = line.split()[1] # structure number (string) of the pure metal case

        elif line.split()[0] == 'MAX_ITER:':
            maxIter = line.split()[1] # maximum number of iterations

        elif line.split()[0] == 'DISTRIBUTE:': # Whether to distribute over atom jobs the tasks such as iid choosing, fitting, and ground state search
            if line.split()[1][0].lower() == 'y': 
               distribute = True  
            elif line.split()[1][0].lower() == 'n': 
               distribute = False 
        elif line.split()[0] == 'EXTEND_PATH:':
            extendpath = line.split()[1]      
    
    return [atoms, volRange, clusterNums, runTypes, PriorOrIID, niid, maxiid, mfitStructs, nfitSubsets, nPrior, plotTitle, xlabel, \
            ylabel,restartTimeout, rmStructIn, ediffg, maxE, graphsOnly, pureMetal, maxIter, distribute, extendpath]

def removeStructs(list1,list2):
    '''Remove items from list1 that might be in list2, and return list2'''
    for iatom, atom in enumerate(atoms):
        for item in list1[iatom]:
            if item in list2[iatom]: list2[iatom].remove(item)
    return list2

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
    """ Retrieves the number of C, H and M atoms from the POSCAR file and sets 
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
        atomCounts.append(int(counts[0]))
        atomCounts.append(int(counts[1]))
        atomCounts.append(int(counts[2]))
    elif len(counts) == 2:
        atomCounts.append(int(counts[0]))
        if poscarLines[0].split()[1] == 'H':
            atomCounts.append(int(counts[1]))
            atomCounts.append(0)
        elif poscarLines[0].split()[1] == 'M':
            atomCounts.append(0)
            atomCounts.append(int(counts[1]))
    natoms = sum(atomCounts)
    if fixPOSCAR:
        del(poscarLines[5])
        writefile(poscarLines[:7+natoms],poscarDir + '/POSCAR') #:7+natoms is because CONTCAR includes velocity lines that uncle doesn't want. The factor of 2 is because carbon atoms are not included in natoms      
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
# -------------------------------------------- MAIN -----------------------------------------------
# -------------------------------------------- MAIN -----------------------------------------------
# -------------------------------------------- MAIN -----------------------------------------------


'''All atom folders must start with the same method:  
1) structures.in from previous run must exist in each atom folder
or
2) structures.in must be absent from each atomfolder. Structure folders must exist from previous runs (including pure cases) in atom folders, which can be a mix of finished and unfinished structures.  
s from previous run must exist
or
3) if all the above are missing from all folders, each atom will start from scratch with iid structures'''          
if __name__ == '__main__':
    maindir = os.getcwd() # leave this as default!
    
# override default maindir   

#    maindir = '/fslhome/bch/cluster_expansion/graphene/hollowTivac.v8'  
#    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1'
#    maindir = '/fslhome/bch/cluster_expansion/graphene/top.tm_row1.v8'
#    maindir = '/fslhome/bch/cluster_expansion/graphene/top.tm_row1.v15' 
#    maindir = '/fslhome/bch/cluster_expansion/graphene/hollowTiH.v8'
#    maindir = '/fslhome/bch/cluster_expansion/graphene/hollowTiH.v15'
#    maindir = '/fslhome/bch/cluster_expansion/graphene/hollowTopTiH.v8'

    subprocess.call(['echo','Starting in ' + maindir])
    #make sure the latest version of uncle is used
    subprocess.call(['rm',maindir+'/needed_files/uncle.x'])
    subprocess.call(['ln','-s','/fslhome/bch/bin/uncle',maindir+'/needed_files/uncle.x'])
    os.chdir(maindir)
    [atoms, volRange, clusterNums, runTypes, PriorOrIID, niid, maxiid, mfitStructs, nfitSubsets, nPrior, plotTitle, xlabel,\
             ylabel,restartTimeout, rmStructIn, ediffg, maxE, graphsOnly, pureMetal, maxIter, distribute, extendpath] = readSettingsFile()
    #build new larger run from previous smaller run, if maindir is empty and a path is given
#    if extendpath != '' and len(os.listdir(maindir))==0:
    nTotClusters = sum(clusterNums)    
    pathMax = maindir + '/needed_files/diffMax'
    if os.path.exists(pathMax):
        diffMax = float(readfile(pathMax)[0].strip()) 
        subprocess.call(['echo','\ndiffMax read from file (convergence criterium): '+str(diffMax) + '\n'])
    else:
        subprocess.call(['echo','\nMissing diffMax file with the convergence criterium diffMax in eV. Stopping'])
        sys.exit('Stop')
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
    [vstructsFinished,vstructsRestart0,vstructsFailed,startMethod,vdata] = initializeStructs(atoms,restartTimeout,rmStructIn,pureMetal)    
    vstructsAll = joinLists([vstructsFinished,vstructsFailed,vstructsRestart0])
       
    enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, uncleOutput, distribute)
    nTotStructs = enumerator.getNtot(os.getcwd()+'/enum') #number of all enumerated structures
#    subprocess.call(['echo','Warning: BLOCKING ENUMERATOR to save time' ])
    if not enumerationDone(maindir + '/enum',nTotClusters,nTotStructs):
        enumerator.enumerate()
    energiesLast = zeros((natoms,nTotStructs),dtype=float) #energies of last iteration, sorted by structure name

    createEnumPastDir(atoms)
    pastStructsUpdate(vstructsAll,atoms)  
    nNew = getnNew(atoms,vstructsFinished,vstructsRestart,niid,niid,nTotStructs,PriorOrIID) #how many new structures to get for each atom.  If have empty folders, niid is needed for both methods
   
    converged = False
    iteration = 1

#======================================  Iteration loop =======================================
    while not converged:
        
        subprocess.call(['echo','\n========================================================'])
        subprocess.call(['echo','\t\tIteration ' + str(iteration)])
        subprocess.call(['echo','========================================================\n'])
        #since atoms may have changed, have to reinitialize these
        enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, uncleOutput, distribute) 
        extractor = Extractor.Extractor(atoms, uncleOutput, startMethod,pureMetal)   
        # If enough iid structs have been submitted, go to priority method        
        nAllmin = min([len(vstructsAll[iatom]) for iatom in range(len(atoms))])
        if PriorOrIID == 'i' and nAllmin >= maxiid: 
            PriorOrIID = 'p'
        if iteration == 1: 
            if startMethod == 'empty folders': 
                vstructsToStart = enumerator.chooseTrainingStructures(iteration,startMethod,nNew,nTotStructs)
                vstructsToStart = extractor.checkPureInCurrent(iteration,vstructsToStart,vstructsFinished)
            vstructsToRun = vstructsToStart #no restarts in first iteration
        elif iteration > 1 and PriorOrIID == 'p':
            vstructsToStart = newStructsPrior  #from previous iteration 
            contcarToPoscar(vstructsRestart,atoms,iteration) 
            vstructsToRun = joinLists([vstructsRestart,vstructsToStart])
        elif iteration > 1 and PriorOrIID == 'i':
#            subprocess.call(['echo','Warning: BLOCKING chooseTrainingStructures (iter>1) to save time' ])
            vstructsToStart = enumerator.chooseTrainingStructures(iteration,startMethod,nNew,nTotStructs)   
            contcarToPoscar(vstructsRestart,atoms,iteration)
            vstructsToRun = joinLists([vstructsRestart,vstructsToStart])
        pastStructsUpdate(vstructsToStart,atoms)
        
#        #remove these lines
#        subprocess.call(['echo','Warning: Remove these two lines that convert Restart into Start !!!!!' ])
#        vstructsToStart = vstructsRestart
#        vstructsRestart = [[]]*natoms
#        #remove above

        finalDir = extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsToStart,vstructsRestart)           
        os.chdir(maindir) #test this...shouldn't need it.
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atoms, startMethod, iteration, finalDir, restartTimeout,pureMetal) 
        [newlyFinished, newlyToRestart, newlyFailed,vdata] = uncleFileMaker.makeUncleFiles(iteration, holdoutStructs,vstructsToRun,vdata) 
        vstructsFinished = joinLists([vstructsFinished,newlyFinished])
        vstructsFailed = joinLists([vstructsFailed,newlyFailed])
        vstructsAll = joinLists([vstructsFinished,vstructsFailed,newlyToRestart])
        os.chdir(maindir) #test this...shouldn't need it.  

        # Get all the structs that have been through VASP calculations for each atom. These
        # should be sorted by formation energy during the work done by makeUncleFiles()
                
        # Perform a fit to the VASP data in structures.in for each atom.
        if iteration > 1 or not existAllJ1out(atoms):
            fitter = Fitter.Fitter(atoms, mfitStructs, nfitSubsets, vstructsFinished,uncleOutput,distribute)
            if iteration == 1: 
                fitter.makeFitDirectories()
            fitter.writeHoldout(50,vstructsFinished,vdata)
#            subprocess.call(['echo','Warning: BLOCKING FITS to save time' ])
            fitter.fitVASPData(iteration,maxE)
    
        # Perform a ground state search on the fit for each atom.    
        gss = GSS.GSS(atoms, volRange, plotTitle, xlabel, ylabel, vstructsFinished,uncleOutput,pureMetal,finalDir,distribute)
        gss.makeGSSDirectories()
#        subprocess.call(['echo','Warning: BLOCKING GSS to save time' ])   
        gss.performGroundStateSearch(iteration)
        gss.makePlots(iteration)
        #get the priority of each structure in each atom
        priorities = gss.getGssInfo(iteration,vstructsFailed) #first structure listed is highest priority
        minPrior = 0.01
        plotStructsByPrior(atoms,minPrior,iteration)
        collateStructsConc(atoms,minPrior,iteration)
        NInPlot = 400
        collateStructsHFE(atoms,minPrior,NInPlot,iteration)              
        move = MovementInfo.MovementInfo(atoms,pureMetal,iteration,False) #instance
        move.getMovementInfo()        
        if graphsOnly: sys.exit('Done with graphs. Stopping')

                #---------- prep for next iteration --------------
        vstructsRestart = newlyToRestart
        #determine how many new structures to get for each atom, for any method
        if iteration == 1: #bring in restarts from initial folders
            vstructsRestart = joinLists([vstructsRestart0,vstructsRestart])
        nNew = getnNew(atoms,vstructsFinished,vstructsRestart,niid,nPrior,nTotStructs,PriorOrIID) #how many new structures to get for each atom
        #choose new structures while still sorted by priority
        newStructsPrior = getFromPriorities(priorities,vstructsAll,atoms,nNew)       
        # First sort priorities by structure name so we have an unchanging order
        for iatom, atom in enumerate(atoms): 
            priorities[iatom,:] = sort(priorities[iatom,:],order = ['struct'])
#        for iatom, atom in enumerate(atoms):       
#            print 'priorities for atom, iteration',atom, iteration
#            print priorities[iatom,:10]['struct']
#            print priorities[iatom,:10]['FE']
#            print 'energiesLast'
#            print energiesLast[iatom,:10]
        diffe = getDiffE(priorities,energiesLast,atoms)
        energiesLast = priorities['FE'] #to be used next iteration, sorted by structure
#        print 'new energiesLast', energiesLast
        
        # Determine which atoms need to continue
#        diffMax = 0.01 # 10 meV convergence criterion
        diffMax = float(readfile(maindir + '/needed_files/diffMax')[0].strip()) #read from file called 'diffMax', so we can experiment with it during a long run
        subprocess.call(['echo','Weighted energy changes']) 
        rmAtoms = []
        atomsCopy = deepcopy(atoms)
        for iatom, atom in enumerate(atoms):
            subprocess.call(['echo','\t{} {:8.6f} eV'.format(atom,diffe[iatom])])
            if diffe[iatom]<diffMax: #Converged
                subprocess.call(['echo','Atom {} has converged'.format(atom)])
                atomsCopy.remove(atom)
                rmAtoms.append(iatom)
            elif len(newlyFinished[iatom]) + len(vstructsRestart[iatom]) == 0 and len(vstructsToStart[iatom]) > 0:
                subprocess.call(['echo','\tAtom {} did not finish any new structures, and has been stopped.'.format(atom)])
                print 'nFinished,nrestart,nstart',len(newlyFinished[iatom]) , len(vstructsRestart[iatom]), len(vstructsToStart[iatom])
                atomsCopy.remove(atom)
                rmAtoms.append(iatom)
            elif nNew[iatom] == 0 and len(vstructsRestart[iatom]) == 0:
                subprocess.call(['echo','\tAtom {} has no structures to run, and has been stopped.'.format(atom)])
                atomsCopy.remove(atom)
                rmAtoms.append(iatom)
        atoms = atomsCopy
#        print 'rmAtoms',rmAtoms
#        print 'atoms',atoms
#        print 'vstructsRestart',vstructsRestart
                              
        natoms = len(atoms) #could be lower now
        if len(rmAtoms)>0:
            vstructsFinished = multiDelete(vstructsFinished,rmAtoms)
            vstructsFailed = multiDelete(vstructsFailed,rmAtoms)
            vstructsRestart = multiDelete(vstructsRestart,rmAtoms)
#            vstructsRestart0 = multiDelete(vstructsRestart0,rmAtoms)
            vstructsAll = multiDelete(vstructsAll,rmAtoms)
            vdata = delete(vdata,rmAtoms,axis=0)        
            energiesLast = delete(priorities,rmAtoms,axis=0)['FE'] #    priorities[ikeep,:]['FE']
            nNew = multiDelete(nNew,rmAtoms)
            newStructsPrior = multiDelete(newStructsPrior,rmAtoms)
#        print 'vstructsRestart',vstructsRestart
        if natoms == 0:
            subprocess.call(['echo','\n----------------- All atoms have finished ---------------'])
            converged = True
        elif iteration == maxIter:
            subprocess.call(['echo','\n----------------- Finished maximum number of iterations ---------------'])
            converged = True            
        else:       
            if natoms >1:
                subprocess.call(['echo','At end of iteration, {} are continuing.'.format(' '.join(atoms))])
            else:
                subprocess.call(['echo','At end of iteration, {} is continuing.'.format(' '.join(atoms))])


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