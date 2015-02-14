import os, subprocess,sys,re,time
from numpy import mod, floor,zeros,array,sqrt,std,amax,amin,int32,sort,count_nonzero,delete
from sched import scheduler

def collatePlotsGSS(self,plotType,iteration):  
    '''Creates an HTML page with the plots and labels, from plots in /gss
    plotType: first part of plot file name'''
    lastDir = os.getcwd()
    plotsDir = lastDir + '/plots'
    if not os.path.exists(plotsDir):
        subprocess.call(['mkdir',plotsDir]) 
    subDir = plotsDir + '/plots{}'.format(plotType)
    if not os.path.exists(subDir):
        subprocess.call(['mkdir',subDir])        
    nRow = 3  # number of plots in row
    width = 1350
    height  = 900
    collatefile  = open(subDir +'/plots{}_{}.htm'.format(plotType,iteration),'w')
    collatefile.write(' <html>\n <HEAD>\n<TITLE> {} </TITLE>\n</HEAD>\n'.format(plotType))
    collatefile.write(' <BODY>\n <p style="font-size:20px"> <table border="1">\n <tr>\n') #start row

    for iImage,atom in enumerate(self.atoms):
        path =  '../../{}/gss/{}_{}.png'.format(atom,plotType,iteration)  
        atomtext = atom.split('_')[0]       
        collatefile.write('<td><p><img src="{}" width "{}" height "{}" ></p><p>{}</p></td>\n'.format(path,width,height,''))#Image and element under it
        if mod(iImage+1,nRow) == 0: 
            collatefile.write('</tr>\n<tr>\n') #end of row, begin new
    collatefile.write(' </tr></table> \n') #end of row and table                
    collatefile.write(' </BODY> </html>') #end of file 
    collatefile.close()  

def convergeCheck(folder, NSW):
    """ Tests whether force convergence is done by whether the last line of OSZICAR (the last
        ionic relaxation step) is less than NSW."""
    try:
        value = getSteps(folder)
        return value < NSW  
    except:
        return False  

def energyDropCheck(dir):
    '''tests whether the energies have dropped in OSZICAR...rising energies are unphysical 
    and show a numerical convergence problem. The factor like 0.99 allows for very small energy rises only'''
    lines = readfile(dir+'/OSZICAR') 
    energies = []
    for line in lines:
        try:
            if 'F=' in line and len(line.split())>2:
                energies.append(float(line.split()[2]))
        except:
            'do nothing' # this allows to continue if file is corrupt
    if len(energies) > 2:
        return energies[-1] <= 0.99*energies[0] 
    return False
    
def finishCheck(folder):
    """ Tests whether VASP is done by finding "Voluntary" in last line of OUTCAR.  The input
        parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """   
    lastfolder = os.getcwd()
    os.chdir(folder)        
    proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'], stdout=subprocess.PIPE)
    newstring = proc.communicate()
    os.chdir(lastfolder)            
    return newstring[0].find('Voluntary') > -1

def flat(list):
    '''flattens a list of lists (no further depth allowed'''
    return sum(list,[])
    
def getNSW(dir): 
    if os.path.exists(dir+'/INCAR'):
        proc = subprocess.Popen(['grep','-i','NSW',dir+'/INCAR'],stdout=subprocess.PIPE) 
        return int(proc.communicate()[0].split('=')[-1])
    else:
        return 1

def getPureEs(self, iatom):
    lastDir = os.getcwd()
    dir = lastDir + '/' + self.atoms[iatom]
    os.chdir(dir)

    pureEnergies = []
    pureStructs = getPureStructs(lastDir + '/enum')

    for i, iDir in enumerate(pureStructs):
        pureDir = dir + '/' + str(iDir)
    
        if os.path.exists(pureDir):
            self.setAtomCounts(pureDir)
            self.setEnergy(pureDir)
            pureEnergies.append(float(self.energy))
            subprocess.call(['echo','Pure ' + self.atoms[iatom].split(',')[i] +' energy: {}'.format(self.energy) ]) 
        else:
            subprocess.call(['echo','Missing pure ' + self.atoms[iatom].split(',')[i] + ' energy folder'])

#                    if etest != 999999:
#                        self.pureMenergy = etest
#                    else:
    self.pureEnergies = pureEnergies
    os.chdir(lastDir)

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

def hexMonolayerEnergies(self,dir1,iteration): 
    file = open(dir1 +'/hex_monolayer_refs/hex_energies','w')
    self.hexE = {}
    self.hexE['Vc'] = 0.0
    if iteration == 1: subprocess.call(['echo', '\nReading hexagonal monolayer energies\n'])
    for iatom,atom in enumerate(self.atoms):
        for element in atom.split(','):
            if element.find('Vc') == -1:
                dir2 = dir1 + '/hex_monolayer_refs'+'/'+element
                if finishCheck(dir2) and convergeCheck(dir2, getNSW(dir2)) and energyDropCheck(dir2): #finalDir
                    if iteration == 1: subprocess.call(['echo','{} monolayer (per atom): {:8.4f} '.format(element,self.getEnergy(dir2))])
                    file.write('{} monolayer (per atom): {:8.4f} \n'.format(element,self.getEnergy(dir2)))
                    self.hexE[element] = self.getEnergy(dir2) 
                else:
                    file.write('{} monolayer not converged \n'.format(element))
    os.chdir(dir1)  

def isequal(x,y):
    eps = 5.0e-5
    return abs(x-y)<eps

def joinLists(toJoin):
    '''Joins lists in toJoin of the form [[sublist1],[sublist2],[sublist3]...].  Lists must have the
    same length, but can different length sublists'''    
    joinedList = [[]]*len(toJoin[0])
    for isublist in range(len(toJoin[0])):
        sublist=[]
        for ilist in range(len(toJoin)):
            sublist += toJoin[ilist][isublist]
        joinedList[isublist] = sublist
    return joinedList

def outcarWarn(dir):
    lines = readfile(dir+'/OUTCAR')

    for i,line in enumerate(lines):
        if 'W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!' in line\
          and 'The number of bands has been changed' not in lines[i+2]: #this is OK 
            return True
    return False
  
def readfile(filepath):
    file1 = open(filepath,'r')
    lines = file1.readlines()
    file1.close()
    return lines
    
def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return
            
def setAtomCounts(self, poscarDir):
    """ Retrieves the number of C, H and M atoms from the POSCAR file and sets 
        the corresponding members. 
        Also fixes "new" POSCAR/CONTCAR format (comes from CONTCAR) back to old for UNCLE use (removes the 6th line if it's text """
    self.atomCounts = []

    vacancyNum = -1
    vacancies = 0

    fixPOSCAR = False
    poscarLines = readfile(poscarDir + '/POSCAR')
    counts = poscarLines[5].strip().split() 

    if not counts[0][0].isdigit(): # then we have the "new" POSCAR format that gives the atom types in text on line 6 (5-python)
        fixPOSCAR = True
        counts = poscarLines[6].strip().split()  

    case = int(poscarLines[0][poscarLines[0].find('Case:') + 6])

    present = poscarLines[0][poscarLines[0].find('atoms:') + 6:poscarLines[0].find('Case')].strip()
    countnum = 1

    self.atomCounts.append(int(counts[0]))

    for i in range(case):
        if present.find(str(i+1)) != -1:
            self.atomCounts.append(int(counts[countnum]))
            countnum = countnum + 1
        else:
            self.atomCounts.append(0)

    if poscarLines[0].find('Vacancies') != -1:
        vacancyNum = int(poscarLines[0].strip().split()[-1])
        vacancies = int(poscarLines[0].strip().split()[-4])
        self.atomCounts[vacancyNum] = vacancies

    natoms = sum(self.atomCounts)
    if fixPOSCAR:
        del(poscarLines[5])
        writefile(poscarLines[:7+natoms],poscarDir + '/POSCAR') #:7+natoms is because CONTCAR includes velocity lines that uncle doesn't want         

def setEnergy(self, directory):  
    """ Retrieves the energy of the structure from the OSZICAR file and sets the corresponding 
        member. """   
    try:
        oszicar = open(directory + self.finalDir + '/OSZICAR','r')
        energy = oszicar.readlines()[-1].split()[2]
        oszicar.close()
    except:
        energy = 0
    
    energy = float(energy)
    peratom = energy / float(sum(self.atomCounts[1:]))      
    self.energy = str(peratom)

def singleAtomsEnergies(self,dir1,iteration): 
    self.singleE = {}
    self.singleE['Vc'] = 0.0
    if iteration == 1: subprocess.call(['echo', '\nReading single atom energies\n'])
    file = open(dir1 +'/single_atoms/single_atom_energies','w')
    for iatom,atom in enumerate(self.atoms):
        for element in atom.split(','):
            if element.find('Vc') == -1:
                dir2 = dir1 + '/single_atoms'+'/'+element
                if self.electronicConvergeFinish(dir2): 
                    if iteration == 1: subprocess.call(['echo', 'Energy of {} atom: {:8.4f} \n'.format(element,self.getEnergy(dir2))])
                    file.write('{} atom: {:12.8f} \n'.format(element,self.getEnergy(dir2)))
                    self.singleE[element] = self.getEnergy(dir2)
    file.close()  
    os.chdir(dir1)      


def structuresWrite(howmany,atomDir,structlist, FElist,conclist,energylist,outType,writeType):
    '''Goes back to makestr.x in case POSCAR has been changed (by overwriting with CONTCAR for example)
       Also writes this info as POSCAR_orig in the structure folder.  outType is '.in' or '.holdout'.  writeType is either "w" or "a"
       depending on whether you are starting the file or appending to it.'''
    if howmany == 'all':
        N = len(structlist)
    else:
        N = min(howmany,len(structlist)) 
    lastDir = os.getcwd()
    existed = os.path.exists(atomDir + '/' + 'structures{}'.format(outType))
    outFile = open(atomDir + '/'+ 'structures{}'.format(outType), writeType)                                   
    if not existed:
        outFile.write("peratom\nnoweights\nposcar\n"); outFile.flush()  #header
    os.chdir(atomDir)
    subprocess.call(['ln','-s','../enum/struct_enum.out'])
#    subprocess.call(['rm vasp.0*'],shell=True) #start clean
    for istruct,struct in enumerate(structlist[:N]): #just take first N for now.  Can change to a slice later
        vaspDir = atomDir + '/'+ str(struct)
        outFile.write("#------------------------------------------------\n")

        subprocess.call(['../needed_files/makestr.x','struct_enum.out',str(struct)])
        vfile = 'vasp.' + '0'*(6-len(str(struct))) + str(struct)
        poscar = readfile(atomDir+'/'+vfile)
        if os.path.exists(vaspDir): subprocess.call(['cp', vfile, vaspDir+'/puedoPOSCAR'])
        subprocess.call(['rm', vfile])           
        idString = 'graphene str #: ' + str(struct)
        outFile.write(idString + " FE = " + str(FElist[istruct]) + ", Concentration = " + str(conclist[istruct]) + "\n")
        outFile.write("1.0\n")
        writeLatticeVectors(poscar[2:5],outFile)
        outFile.writelines(poscar[5:])
        outFile.write("#Energy:\n")
        outFile.write(str(energylist[istruct]) + "\n")         
    outFile.close() 
    os.chdir(lastDir)

def getPureStructs(enumDir):
    enumfile = open(enumDir + '/struct_enum.out', 'r')
    lines = [line for line in enumfile]
    enumfile.close()

    pureStructs = []

    for i, line in enumerate(lines):
        if line.find('#tot') != -1:
            start = i
            break
    end = 200
    if len(lines) < end:
        end = len(lines)

    for line in lines[start:end]:
        label = line.strip().split()[-1]
        if len([num for num in range(10) if str(num) in label]) == 1:
            pureStructs.append(int(line.strip().split()[0]))
    return pureStructs

def writeLatticeVectors(vecLines,outfile):
    """ Gets the lattice vectors and writes them"""  
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

    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec1z, vec1x, vec1y))
    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec2z, vec2x, vec2y))
    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec3z, vec3x, vec3y))   
    
def parallelJobFiles(atoms,subdir,walltime,mem,execString,atomStrings):
    """ Create job files for parallelizing over atoms.
    Walltime in decimal hours. Mem in Gb."""
    lastDir = os.getcwd()
    hrs = int(floor(walltime)); mints = int(mod(walltime,1)*60)
    for iatom, atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom
#        try:               
        os.chdir(atomDir + '/' + subdir)
        jobFile = open('job','w')
        jobFile.write('#!/bin/bash\n')
        jobFile.write('#SBATCH --time={}:{}:00\n'.format(hrs,mints))
        jobFile.write('#SBATCH --ntasks=1\n')
        jobFile.write('#SBATCH --mem-per-cpu={}G\n'.format(mem))
        jobFile.write('#SBATCH --mail-user=hess.byu@gmail.com\n')
        jobFile.write('#SBATCH --mail-type=FAIL\n\n')
        jobFile.write('#SBATCH --mail-type=end\n\n')
        jobFile.write("#SBATCH --job-name={}.{}\n\n".format(atom,subdir)) 
        jobFile.write(execString + ' ' + atomStrings[iatom] + ' > {}job.out'.format(subdir))
        jobFile.close()
#        except:
#            subprocess.call(['echo','\n~~~~~~~~~~ Failed while writing atom job files for ' + atom + '! ~~~~~~~~~~\n'])
    os.chdir(lastDir)

def parallelAtomsSubmit(atoms,subdir):
    lastDir = os.getcwd()
    jobIds = []    
    for iatom, atom in enumerate(atoms):
        atomDir = lastDir + '/' + atom
        try:               
            os.chdir(atomDir + '/' + subdir)
            proc = subprocess.Popen(['sbatch', 'job'], stdout=subprocess.PIPE)
            newid = proc.communicate()[0].split()[3]
            jobIds.append(newid)
        except:
            subprocess.call(['echo','\n~~~~~~~~~~ Failed while submitting atom job files for ' + atom + '! ~~~~~~~~~~\n'])
        os.chdir(lastDir)
    return jobIds

def parallelAtomsWait(jobIds):
    """ Waits for the training structures to be generated for all of the atoms before moving on.
        Checks every thirty seconds to see if the jobs have finished. """
    s = scheduler(time.time, time.sleep)    
    finished = False
    start_time = time.time()
    event_time = start_time
    subprocess.call(['echo','\n\tWaiting for above atom jobs to finish...\n'])
    while not finished:
        event_time += 10 #check every x seconds 
        s.enterabs(event_time, 1, doNothing, ([1]))
        s.run()
        finished = reportFinished(jobIds)

def doNothing(param1):
    """ This is just to satisfy the scheduler's requirements to call a method in the enterabs()
        command. """
    pass

def reportFinished(jobIds):
    """ Reports whether the job ids have dropped from the queue. """
    devnull = open(os.devnull, 'w')
    for jobid in jobIds:
        proc = subprocess.Popen(['squeue', '--job', jobid], stdout=subprocess.PIPE, stderr=devnull)
        output = proc.communicate()[0].split()
        if len(output) != 8 and len(output) != 0:   # It will list either all the headers or
            return False                            # sometimes an error and outputs nothing.
                                                    # The error in this case is an "invalid
                                                    # job id" error because the job is no
    return True                                     # longer on the supercomputer.  

def trimSmall(list_mat):
    low_values_indices = abs(list_mat) < 1.0e-3
    list_mat[low_values_indices] = 0.0
    return list_mat

def contains(self, struct, alist):
    """ Returns True if the list 'alist' contains the structure 'struct', False otherwise. """
    if len(alist) <= 0:
        return False
    for elem in alist:
        if str(struct) == str(elem):
            return True
    return False

def elConvergeCheck(self,folder,NELM):  
    """Tests electronic convergence is done by whether the electronic step is less than NELM."""
    try:
        value = self.getElSteps(folder)
        return value < NELM #True/False
    except:
        return False #True/False
                    
def electronicConvergeFinish(self,dir): 
    '''Test requires electronic convergence AND vasp finishing'''
    #get NELM, the max electronic steps allowed in the run.  Using first directory in dirslist
    proc = subprocess.Popen(['grep','-i','NELM',dir+'/INCAR'],stdout=subprocess.PIPE)
    result =  proc.communicate()[0]
    NELM = int(result.split('=')[1].split()[0])
    return self.elConvergeCheck(dir,NELM) and finishCheck(dir)   

def getElSteps(self,folder): 
    '''number of electronic steps for isolated runs'''
    lastfolder = os.getcwd()
    os.chdir(folder)
    try:
        oszicar = open('OSZICAR','r') 
        laststep = oszicar.readlines()[-2].split()[1] # Next to last line, 2nd entry
        oszicar.close()
        os.chdir(lastfolder) 
        value = int(laststep)
        return value         
    except:
        os.chdir(lastfolder)         
        return 9999    
    
def getEnergy(self,dir): 
    lines = readfile(dir+'/OSZICAR')
    if len(lines[-1].split())>1:
        energy = float(lines[-1].split()[2])  #Free energy
    else: 
        energy = 0.0
    return energy
