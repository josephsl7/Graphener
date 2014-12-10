import os, subprocess,sys,re

def convergeCheck(folder, NSW):
    """ Tests whether force convergence is done by whether the last line of OSZICAR (the last
        ionic relaxation step) is less than NSW."""
    try:
        value = getSteps(folder)
        return value < NSW  
    except:
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
            
def setAtomCounts(self,poscarDir):
    """ Retrieves the number of H atoms and the number of M atoms from the POSCAR file and sets 
        the corresponding members. 
        Also fixes "new" POSCAR/CONTCAR format (comes from CONTCAR) back to old for UNCLE use (removes the 6th line if it's text """
    self.atomCounts = []
    fixPOSCAR = False
    poscarLines = readfile(poscarDir + '/POSCAR')
    counts = poscarLines[5].strip().split() 
    if not counts[0][0].isdigit(): # then we have the "new" POSCAR format that gives the atom types in text on line 6 (5-python)
        fixPOSCAR = True
        counts = poscarLines[6].strip().split()  
    if len(counts) == 3:
        self.atomCounts.append(int(counts[1]))
        self.atomCounts.append(int(counts[2]))
    elif len(counts) == 2:
        if poscarLines[0].split()[1] == 'H':
            self.atomCounts.append(int(counts[1]))
            self.atomCounts.append(0)
        elif poscarLines[0].split()[1] == 'M':
            self.atomCounts.append(0)
            self.atomCounts.append(int(counts[1]))
    natoms = sum(self.atomCounts)
    if fixPOSCAR:
        del(poscarLines[5])
        writefile(poscarLines[:7+natoms*2],poscarDir + '/POSCAR') #:7+natoms is because CONTCAR includes velocity lines that uncle doesn't want          

def structuresInWrite(atomDir, structlist, FElist, conclist,energylist,writeType):
    '''Goes back to makestr.x in case POSCAR has been changed (by overwriting with CONTCAR for example)
       Also writes this info as POSCAR_orig in the structure folder'''
    lastDir = os.getcwd()
    structuresInFile = open(atomDir + '/'+ 'structures.in', writeType)                                   
    structuresInFile.write("peratom\nnoweights\nposcar\n"); structuresInFile.flush()  #header
    os.chdir(atomDir)
    subprocess.call(['ln','-s','../enum/struct_enum.out'])
    subprocess.call(['rm vasp.0*'],shell=True) #start clean
    for istruct,struct in enumerate(structlist):
        vaspDir = atomDir + '/'+ str(struct)
        structuresInFile.write("#------------------------------------------------\n")

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

    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec1z, vec1x, vec1y))
    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec2z, vec2x, vec2y))
    outfile.write("  %12.8f  %12.8f  %12.8f\n" % (vec3z, vec3x, vec3y))   
    