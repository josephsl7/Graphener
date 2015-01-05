'''Run DOS for all structure folders in maindir'''
import os, subprocess,sys,re
from comMethods import convergeCheck,finishCheck,getNSW,getSteps,readfile,writefile
from numpy import rint

def makeDOSDirectories(structlist,atom,NkDiv,walltime):
    """ After the low-precision relaxation, creates a directory for the Density of States
        run and populates it with the files from the low-precision run. 
        Copies the CONTCAR to the DOS POSCAR. """  
    topDir = os.getcwd()
    atomDir = topDir
    toStart = []
    for struct in structlist:
        structDir = atomDir + '/' + struct
        dosDir = atomDir + '/' + struct + '/DOS'
        if not os.path.exists(structDir):
            print structDir
            subprocess.call(['echo','No folder for {}, struct {}'.format(atom,struct)])
        elif finishCheck(structDir) and  convergeCheck(structDir, getNSW(structDir)):
            if not os.path.exists(dosDir):
                subprocess.call(['echo','Preparing {} struct {} for DOS run'.format(atom,struct)]) 
                os.chdir(structDir)
                subprocess.call(['mkdir', dosDir])
                subprocess.call(['ln','-s','/fslhome/bch/bin/vasp533',dosDir+'/vasp533'])
                subprocess.call(['cp',structDir+'/CONTCAR',dosDir+ '/POSCAR'])
                subprocess.call(['cp',structDir+'/POTCAR',dosDir+ '/POTCAR'])
                makeDOS_INCAR(dosDir)
                makeDOSJobFile(dosDir,atom+struct,walltime)
                makeKPOINTS(dosDir, NkDiv, NkDiv) 
                toStart.append(struct)
        else:
            subprocess.call(['echo','Structure ' + struct + ' has not finished or not converged'])
        os.chdir(topDir)
    return toStart

def makeDOS_INCAR(dosDir):
    """ Creates an INCAR file for the Density of States run in VASP.  The notable changes are:
            IBRION=-1 -- This tells VASP not to move the ions.
            NSW=0     -- This tells VASP that there will be no ionic relaxation steps.
            LORBIT=10 -- This creates the PROCAR file which can be used to project onto the 
                         C, H, and M atoms. """
            
    incar = open(dosDir + '/INCAR','w')
    
    incar.write("IBRION=-1\n")
    incar.write("NSW=0\n")
    incar.write("Algo=Normal\n")
    incar.write("PREC=Normal\n")
    incar.write("EDIFF=5E-4\n")
    incar.write("ISMEAR=0\n")
    incar.write("ISPIN=2\n")
    incar.write("NPAR=4\n")
    incar.write("LREAL=Auto\n")
    incar.write("LORBIT=11\n")
    incar.write("SIGMA=0.1\n")
    incar.write("LWAVE=.TRUE.\n")
    incar.write("LCHARG=.TRUE.\n")
    incar.close()  
    
def makeDOSJobFile(dosDir,name,walltime):
    
    jobFile = open(dosDir +'/job','w')
    
    jobFile.write("#!/bin/bash\n\n")
    jobFile.write("#SBATCH --time={}:00:00\n".format(walltime))
    jobFile.write("#SBATCH --ntasks=16\n")
    jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
    jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
    jobFile.write("#SBATCH --mail-type=END\n")  
    jobFile.write("#SBATCH --job-name=dos%s\n" % name)
    jobFile.write("#SBATCH --mail-type=FAIL\n")
    jobFile.write("\nmpiexec vasp533 > vasp.out\n")

    jobFile.close()  

def makeKPOINTS(dir, num1, num2):
    """ Creates a KPOINTS file based on the input parameters num1 and num2. It specifies that 
        the job will have num1 x num2 kpoints. For example, if we wanted to specify an 8x8 
        kpoints mesh, we would call makeKPOINTS(8, 8). """

    kpoints = open(dir + '/KPOINTS','w')

    kpoints.write("Automatic mesh\n")
    kpoints.write("0\n")
    kpoints.write("Gamma\n")
    kpoints.write(str(num1) + ' ' + str(num2) + ' 1\n')
    kpoints.write('0 0 0')

    kpoints.close()

def lastPriorFile(dir):
#    '''Return the latest priorities file name. 
    priorlist = []
    structfiles = os.listdir(dir)
    for file in structfiles:
        filePath = dir + '/' + file
        if 'priorities' in file and os.stat(filePath).st_size > 0:
            priorlist.append(file)
    if len(priorlist)==0: 
        return 'none'
    else: #only use the last prior
        priorlist.sort()
        lastprior = priorlist[-1]
        return lastprior

def startJobs(toStart,atomDir):
    """ Submits all the VASP jobs for structures in 'vstructsToRun' to the supercomputer for 
        low-precision relaxation and records their job IDs. """
    lastDir = os.getcwd()
    for structure in toStart:
        os.chdir(atomDir + '/' + structure + '/DOS')
        proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
        jobid = proc.communicate()[0].split()[3]
        subprocess.call(['echo', 'Submitted job ' + jobid])
    os.chdir(lastDir)
        
    
#======================================= Script =====================================
#======================================= Script =====================================
#======================================= Script =====================================

maindir = '/fslhome/bch/cluster_expansion/graphene/analysis/top.tm_row1/Ti/'
 
##    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1.continue'
#    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1'
#maindir = os.getcwd()

os.chdir(maindir)

walltime = 2 #hrs (will use nearest integer hour)
NkDiv = 30

subprocess.call(['echo','Starting in ' + maindir])
#os.chdir(maindir)
structlist = []
for item in os.listdir(maindir):
    itempath = maindir + '/' + item
    if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
        structlist.append(item)
sstring = ''; 
for struct in structlist: 
    sstring += struct + ' '
subprocess.call(['echo','found structures '+ sstring])
dirname = maindir.split('/')[-2]
toStart = makeDOSDirectories(structlist,dirname,NkDiv,int(rint(walltime)))
startJobs(toStart,maindir)
print "Done submitting jobs"

