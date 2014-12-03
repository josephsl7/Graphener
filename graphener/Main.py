'''
Created on Aug 26, 2014

@author: eswens13
'''

import os, subprocess,sys,re
from random import seed
from numpy import zeros,array,sqrt,std,amax,amin,int32,sort,count_nonzero,delete,mod
from copy import deepcopy

import Enumerator, Extractor, StructsToPoscar, JobManager, MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo     

def checkInitialFolders(atoms,restartTimeout):
    '''assigns all struct folders to either finished, restart, or failed.  Restart is optional.
    Initializes and writes structures.in and vdata.  Uncle formation energy is calculated (the other energies
    will be calculated in vdataToPlotFiles)''' 

    lastDir = os.getcwd()
    subprocess.call(['echo','Checking structure folders for existing vasp data. . . \n'])
    natoms = len(atoms)
    vstructsFinished = [[]]*natoms #every structure vasp has finished before this iteration, a list for each atom
    vstructsFailed = [[]]*natoms #every structure vasp has failed before this iteration, a list for each atom
    vstructsRestart = [[]]*natoms #those that timed out before starting or in previous iteration.  Must have folders in the atom directory.
    startFromExisting = [False]*natoms 
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
        print '\n\n',atom
        print '\tpure H system energy',pureHenergy 
        print '\tpure M system energy',pureMenergy 
        # all structures (including pure
        finishedStructs = []
        restartStructs = []
        failedStructs = []
        structlist = []
        ifinished = -1
        for item in os.listdir(atomDir):
            itempath = atomDir + '/' + item
            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
                structlist.append(item)
        for istruct,struct in enumerate(structlist):
            failed = False
            if mod(istruct+1,10) == 0: print '\tChecking',istruct+1,'of',len(structlist), 'structures in', atom 
            vaspDir = atomDir + '/'+ struct
            try:
                atomCounts = setAtomCounts(vaspDir) #also changes andy "new"-format POSCAR/CONTCAR back to old (removes the 6th line if it's text
            except:
                failedStructs.append(struct)
                failed = True
            if not failed and finishCheck(vaspDir) and convergeCheck(vaspDir, getNSW(vaspDir)): 
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
            elif not failed and restartTimeout and timeoutCheck(vaspDir):
                restartStructs.append(struct)  
            else:
                failedStructs.append(struct)
        if len(finishedStructs) > 0 : startFromExisting[iatom] = True
        vstructsFinished[iatom] =  finishedStructs 
        vstructsRestart[iatom] = restartStructs
        vstructsFailed[iatom] = failedStructs 
        subprocess.call(['echo','\n\tAtom {}: found {} finished, {} to restart, {} failed'.format(atom,len(finishedStructs),len(restartStructs),len(failedStructs))])
        #sort by FE
        nfinished = count_nonzero(vdata[iatom,:]['struct'])
        vdata[iatom,:nfinished] = sort(vdata[iatom,:nfinished],order = ['FE']) #must sort only the filled ones 
        
        #write structures.in, (or just the header if there is no existing data)
        structuresInWrite(atomDir,vdata[iatom,:nfinished]['struct'], vdata[iatom,:nfinished]['FE'],vdata[iatom,:nfinished]['conc'],vdata[iatom,:nfinished]['energy'])
    
    return  vstructsFinished, vstructsRestart, vstructsFailed, startFromExisting, vdata               
    os.chdir(lastDir)

def structuresInWrite(atomDir, structlist, FElist, conclist,energylist):
    '''Goes back to makestr.x in case POSCAR has been changed (by overwriting with CONTCAR for example)
       Also writes this info as POSCAR_orig in the structure folder'''
    lastDir = os.getcwd()
    structuresInFile = open(atomDir + '/'+ 'structures.in', 'w')                                   
    structuresInFile.write("peratom\nnoweights\nposcar\n"); structuresInFile.flush()  #header
    os.chdir(atomDir)
#    print 'os.listdir(os.getcwd())',    os.listdir(os.getcwd())
    subprocess.call(['ln','-s','../enum/struct_enum.out'])
    subprocess.call(['rm vasp.0*'],shell=True) #start clean
    for istruct,struct in enumerate(structlist):
        vaspDir = atomDir + '/'+ str(struct)
        structuresInFile.write("#------------------------------------------------\n")
#        cmdstr = '../needed_files/makestr.x struct_enum.out ' + str(struct)
#        print cmdstr
#        os.system(cmdstr)
        subprocess.call(['../needed_files/makestr.x','struct_enum.out',str(struct)])
        subprocess.call(['mv vasp.0* {}'.format(vaspDir+'/POSCAR_orig')],shell=True)
        poscar = readfile(vaspDir+'/POSCAR_orig')           
        idString = 'graphene str #: ' + str(struct)
        structuresInFile.write(idString + " FE = " + str(FElist[istruct]) + ", Concentration = " + str(conclist[istruct]) + "\n")
        structuresInFile.write("1.0\n")
        writeLatticeVectors(poscar[2:5],structuresInFile)
        structuresInFile.writelines(poscar[5:])
        structuresInFile.write("#Energy:\n")
        structuresInFile.write(str(energylist[istruct]) + "\n")         
    structuresInFile.close() 
    os.chdir(lastDir)

def writeLatticeVectors(vecLines,outfile):
    """ Gets the lattice vectors from the first structure in the newlyFinished and sets the 
        corresponding member components. """  
    vec1 = vecLines[0].strip().split()
    vec2 = vecLines[1].strip().split()
    vec3 = vecLines[2].strip().split()
       
    vec1comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec1] #to handle ******
    vec2comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec2]
    vec3comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec3]
    
    vec1z = vec1comps[0]
    vec1x = vec1comps[1]
    vec1y = vec1comps[2]
    
    vec2z = vec2comps[0]
    vec2x = vec2comps[1]
    vec2y = vec2comps[2]
    
    vec3z = vec3comps[0]
    vec3x = vec3comps[1]
    vec3y = vec3comps[2]
        
#    vec1x = vec1comps[0]
#    vec1y = vec1comps[1]
#    if vec1comps[2] == 15.0:
#        vec1z = 1000.0
#    else:
#        vec1z = vec1comps[2]
#    
#    vec2x = vec2comps[0]
#    vec2y = vec2comps[1]
#    if vec2comps[2] == 15.0:
#        vec2z = 1000.0
#    else:
#        vec2z = vec2comps[2]
#    
#    vec3x = vec3comps[0]
#    vec3y = vec3comps[1]
#    if vec3comps[2] == 15.0:
#        vec3z = 1000.0
#    else:
#        vec3z = vec3comps[2]
    outfile
    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec1z, vec1x, vec1y))
    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec2z, vec2x, vec2y))
    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec3z, vec3x, vec3y))
    
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
    '''Prepares structures for restart'''
    lastDir = os.getcwd()
    restartStructs = [[]]*len(atoms)
    for iatom, atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom
        for struct in structs[iatom]:
            structDir = atomDir + '/' + str(struct)
            os.chdir(structDir)
            subprocess.call(['cp','POSCAR','POSCAR_{}'.format(iteration)] )
            if os.stat('CONTCAR').st_size > 0:
                subprocess.call(['cp','CONTCAR','POSCAR'] )   
    os.chdir(lastDir)

def convergeCheck(folder, NSW):
    """ Tests whether force convergence is done by whether the last line of OSZICAR (the last
        ionic relaxation step) is less than NSW."""
    try:
        value = getSteps(folder)
        return value < NSW  
    except:
        return False  

def createEnumPastDir(atoms):
    '''enumpast is a folder for storing past_structs.dat, and enumerating new iid structures 
    that need past_structs.dat. This creates/clears the folder and puts an empty past_structs.dat there'''
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

def extractToVasp(iteration,runTypes,atoms,vstructsAll,vstructsToStart,vstructsToRun):
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
        manager2 = JobManager.JobManager(atoms)
        if runTypes ==['low']:
#            print 'BLOCKING RUN MANAGER'
            manager2.runLowJobs(vstructsToStart,vstructsToRun)
        elif runTypes ==['low','normal']:
            manager2.runLowJobs(vstructsToStart,vstructsToRun)
            manager2.runNormalJobs(vstructsToStart,vstructsToRun)          
        elif runTypes ==['low','normal','dos']:
            manager2.runLowJobs(vstructsToStart,vstructsToRun)
            manager2.runNormalJobs(vstructsToStart,vstructsToRun)
            manager2.runDOSJobs(vstructsToStart,vstructsToRun) 
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
    
def finishCheck(folder):
    """ Tests whether VASP is done by finding "Voluntary" in last line of OUTCAR.  The input
        parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """   
    lastfolder = os.getcwd()
    os.chdir(folder)        
    proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'], stdout=subprocess.PIPE)
    newstring = proc.communicate()
    os.chdir(lastfolder)            
    return newstring[0].find('Voluntary') > -1

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
        ediff = priorities[iatom,:]['FE'] - energiesLast[iatom]
        priorsum = sum(priorities[iatom,:]['prior'])
        ediffL1[iatom] = sqrt(sum(abs(ediff) * priorities[iatom,:]['prior']/priorsum))
    return ediffL1

def getNSW(dir): 
    proc = subprocess.Popen(['grep','-i','NSW',dir+'/INCAR'],stdout=subprocess.PIPE) 
    return int(proc.communicate()[0].split('=')[-1])

def getSteps(folder):
    """ Returns the number of ionic relaxation steps that VASP performed, as an integer. """
    lastfolder = os.getcwd()
    os.chdir(folder)
    if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
        os.chdir(lastfolder) 
        return -9999
    oszicar = open('OSZICAR', 'r')
    laststep = oszicar.readlines()[-1].split()[0]
    oszicar.close()
    os.chdir(lastfolder)  
    try:
        value = int(laststep)
        return value
    except:
        return 9999 

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

#def parseStartStructures(atoms,startFromExisting,vstructsFinished,vdata):
#    """ Returns the structures, energies, concentrations from structures.start, 
#    
#    This method assumes that the structures.in file will have all of the pure
#    structures and will have the following format for the structure identification line:
#    graphene str #: (structure number) FE = 0.545, Concentration = .473
#    or, for a pure structure:
#    PURE M graphene str #: (structure number) FE = 0.0, Concentration = 1.0 """
#
#    lastDir = os.getcwd()    
#    istruct = 0
#    startList = [[]]*len(atoms)
#
#    for iatom,atom in enumerate(atoms):
#        subList = []
#        if startFromExisting[iatom]:  
#            startFile = lastDir + '/needed_files/structures.start.' + atoms[iatom]
#            subprocess.call(['cp',startFile, lastDir + '/' + atom + '/structures.start'])       
#            if os.path.exists(startFile):      
#                infile = open(startFile,'r')
#                lines = infile.readlines()
#                infile.close()            
#                for j in xrange(len(lines)):
#                    if list(lines[j].strip().split()[0])[:2] == ['#','-']:
#                        if j != len(lines) - 1:
#                            structLine = lines[j + 1].strip().split()
#                            if structLine[0].lower() == 'pure':
#                                subList.append(structLine[5])
#                                vdata[iatom,istruct]['struct'] = structLine[5]
#                            else:
#                                subList.append(structLine[3])
#                                vdata[iatom,istruct]['struct'] = structLine[3]
#                            natomsList = [int(i) for i in lines[j+6].strip().split()]
#                            natoms = natomsList[0]+natomsList[1]
#                            vdata[iatom,istruct]['natoms'] = natoms
#                            vdata[iatom,istruct]['conc'] = float(natomsList[1]/float(natoms))
#                            vdata[iatom,istruct]['energy'] = float(lines[j+9+natoms].strip().split()[0])  
#                            istruct += 1                  
#            else: 
#                subprocess.call(['echo',"ERROR: structures.start."+ atom + ' is missing from needed_files/'])       
#        startList[iatom] = subList    
#    os.chdir(lastDir)
#    return startList, vdata
  
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
#               startFromExisting = [False]*len(atoms) 
#            elif line.split()[1][0].lower() == 'a': #"All" 
#               startFromExisting = [True]*len(atoms)
#            else:
#                try:  
#                    parselist = re.search(':(.+?)#', line).group(1).lower().split()[:len(atoms)] #take only enough to match number of atoms, as it might be an old list
#                    startFromExisting = [True if i == 't' else False for i in parselist]
#                    if len(startFromExisting) < len(atoms): 
#                        subprocess.call(['echo','START_FROM_EXISTING length is too short for the number of atoms! Stopping'])
#                        sys.exit('Stop')
#                except AttributeError:
#                    print "Can't parse START_STRUCTS T's and F's.  Defaulting to no start structs." 
#                    startFromExisting = [False]*len(atoms)               
        elif line.split()[0] == 'PRIORITY/IID:':
            if str(line.split()[1]).lower() == 'i':
                PriorOrIID = 'i'

        elif line.split()[0] == 'RUN_TYPES:':
            try:
                runTypes = re.search(':(.+?)#', line).group(1).lower().split() #extracts only text between START_STRUCTS and #
            except AttributeError:
                print "Can't parse RUN_TYPES. Defaulting to LOW precision runs only."             

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
    
    return [atoms, volRange, clusterNums, runTypes, PriorOrIID, niid, mfitStructs, nfitSubsets, priorNum, plotTitle, xlabel, ylabel,restartTimeout]

def removeStructs(list1,list2):
    '''Remove items from list1 that might be in list2, and return list2'''
    for iatom, atom in enumerate(atoms):
        for item in list1[iatom]:
            if item in list2[iatom]: list2[iatom].remove(item)
    return list2

def searchTimeout(atoms,type,iteration):
    '''Finds structures with slurms that show time out on a previous run.
      Option 'all' is only implemented for the first iteration'''
    lastDir = os.getcwd()
    subprocess.call(['echo','Searching for timed out structures. . . '])
    restartStructs = [[]]*len(atoms)
    for iatom, atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom 
        sublist = []
        dirlist = []
        for item in os.listdir(atomDir):
            itempath = atomDir + '/' + item
            if os.path.isdir(itempath) and item[0].isdigit(): #look at dirs whose names are numbers
                dirlist.append(item)
        for struct in dirlist:
            structDir = atomDir + '/' + struct
            if not finishCheck(structDir): 
                slurmlist = []               
                structfiles = os.listdir(structDir)
                for file in structfiles:
                    filePath = structDir + '/' + file
                    if 'slurm' in file and os.stat(filePath).st_size > 0:
                        slurmlist.append(file)
                if iteration == 1 and type.lower() == 'all' and len(slurmlist)>0: #search all slurm files, in case you corrupted the last one by cancelling the run
                    for slurm in slurmlist:
                        for line in readfile(structDir + '/' + slurm):
                            if 'TIME LIMIT' in line: 
                                sublist.append(struct)
                                break                    
                elif len(slurmlist)>0: #only use the last slurm
                    slurmlist.sort()
                    lastslurm = slurmlist[-1]
                    for line in readfile(structDir + '/' + lastslurm):
                        if 'TIME LIMIT' in line: 
                            sublist.append(struct)
                            break
        restartStructs[iatom] = sublist
        print '\tFor {} found {} structures to restart'.format(atom,len(restartStructs[iatom]))                      
    os.chdir(lastDir)
    return restartStructs

def timeoutCheck(dir):
    '''checks the latest slurm file for the words TIME LIMIT'''
    slurmlist = []               
    structfiles = os.listdir(dir)
    for file in structfiles:
        filePath = dir + '/' + file
        if 'slurm' in file and os.stat(filePath).st_size > 0:
            slurmlist.append(file)
    if len(slurmlist)>0: #only use the last slurm
        slurmlist.sort()
        lastslurm = slurmlist[-1]
        for line in readfile(dir + '/' + lastslurm):
            if 'TIME LIMIT' in line: 
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
        writefile(poscarLines[:7+natoms],poscarDir + '/POSCAR') #:7+natoms is because CONTCAR includes velocity lines that uncle doesn't want          
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

def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

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

    [vstructsFinished,vstructsRestart,vstructsFailed,startFromExisting,vdata] = checkInitialFolders(atoms,restartTimeout) #assigns all struct folders to either finished, restart, or failed''' 
    print 'vstructsFinished';    print vstructsFinished
    print 'vstructsRestart'; print vstructsRestart
    print 'vstructsFailed';print vstructsFailed
    vstructsAll = joinLists(vstructsFinished,vstructsFailed)
    vstructsAll = joinLists(vstructsAll,vstructsRestart)  #may not yet include from structures.start
    pastStructsUpdate(vstructsFinished,atoms)  
    print'BLOCKING ENUMERATOR'   
    enumerator = Enumerator.Enumerator(atoms, volRange, clusterNums, niid, uncleOutput)
#    enumerator.enumerate()
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
        [newlyFinished, newlyToRestart, newlyFailed,vdata] = uncleFileMaker.makeUncleFiles(iteration, holdoutStructs,vstructsToRun,vdata) 
        #update the vstructs lists and past_structs files       
        print '[newlyFinished,newlyToRestart, newlyFailed ]'
        print newlyFinished;print newlyToRestart; print newlyFailed
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
        if iteration == 1: 
            fitter.makeFitDirectories()
            if startFromExisting.count(True) > 0: #at least one starts from existing
                fitter.holdoutFromIn(atoms,startFromExisting)
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
        for iatom in range(natoms): priorities[iatom,:] = sort(priorities[iatom,:],order = ['struct'])
        print 'priorities'
        print priorities[iatom,:10]['struct']
        print priorities[iatom,:10]['FE']
        diffe = getDiffE(priorities,energiesLast,atoms)
        energiesLast = priorities['FE'] #to be used next iteration, sorted by structure
        print 'new energiesLast', energiesLast
        
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

            elif len(newlyFinished[iatom]) + len(newlyToRestart[iatom])   == 0 and len(vstructsToStart[iatom]) > 0:
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
            newStructsPrior = multiDelete(newStructsPrior,rmAtoms)
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


        if natoms == 0:
            subprocess.call(['echo','\n----------------- All atoms have finished ---------------'])
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
        

    
        
    
    
 
 
 
 
 
 
 
    