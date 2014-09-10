import os, subprocess

class RunVasp:
    """ This class is responsible for preparing the directories, retrieving the needed files, and
        submitting VASP jobs to the supercomputer.  It keeps track of the SLURM job ids of all the
        jobs that are currently running from a particular instance of the class. """
    
    def __init__(self, atomList):
        """ CONSTRUCTOR """
        
        self.atomList = atomList
        self.neededFilesFolder = os.getcwd() + '/needed_files/'
        
        self.currJobIds = []
        
    def makeNormalDirectories(self):
        topDir = os.getcwd()
        for element in self.atomList:
            elementDir = topDir + '/' + element
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in os.listdir(os.getcwd()):
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        os.chdir(structDir)
                        subprocess.call(['mkdir', 'normal'])
                        subprocess.call(['cp','CONTCAR','DOSCAR','EIGENVAL',
                                         'IBZKPT','KPOINTS','vasp533',
                                         'OSZICAR','OUTCAR','PCDAT',
                                         'POSCAR','POTCAR','REPORT',
                                         'vasprun.xml','job','XDATCAR','normal'])
                        self.makeNormalINCAR()
                        subprocess.call(['cp','normal/CONTCAR','normal/POSCAR'])
                        os.chdir(elementDir)
            else:
                print "The directory " + elementDir + " does not exist."
            
            os.chdir(topDir)       
    
    def makeDOSDirectories(self):
        """ Creates the directories with the needed files for a Density of States run in VASP.
            The directories are created as sub-directories of the original structure directories. """
            
        topDir = os.getcwd()
        for element in self.atomList:
            elementDir = topDir + '/' + element
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in os.listdir(os.getcwd()):
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        if self.hasFinished(structDir + '/normal'):
                            os.chdir(structDir)
                            subprocess.call(['mkdir', 'DOS'])
                            subprocess.call(['cp','normal/CONTCAR','normal/DOSCAR','normal/EIGENVAL',
                                             'normal/IBZKPT','normal/KPOINTS','normal/vasp533',
                                             'normal/OSZICAR','normal/OUTCAR','normal/PCDAT',
                                             'normal/POSCAR','normal/POTCAR','normal/REPORT',
                                             'normal/vasprun.xml','normal/XDATCAR','DOS'])
                            self.makeDOS_INCAR()
                            self.makeDOSJobFile()
                            subprocess.call(['cp','DOS/CONTCAR','DOS/POSCAR'])
                            os.chdir(elementDir)
                        else:
                            print "Structure " + structure + " did not converge. Skipping. . ."
            else:
                print "The directory " + elementDir + " does not exist."
            
            os.chdir(topDir)
    
    def hasFinished(self, direc):
        """ Tests whether Vasp is done by finding "Voluntary" in last line of OUTCAR.  The input
        parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """
        lastfolder = os.getcwd()
        os.chdir(os.path.abspath(direc))
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)    
        return newstring[0].find('Voluntary') > -1 #True/False
             
    def makeLowINCARs(self):
        """ Creates a standard INCAR file and puts it in each different structure's top directory. """
        
        dirList = self.atomList
        
        for direc in dirList:
            incar = open(direc + '/INCAR','w')
    
            incar.write("IBRION=2\n")
            incar.write("ISIF=4\n")
            incar.write("NSW=400\n")
            incar.write("Algo=VeryFast\n")
            incar.write("PREC=Low\n")
            incar.write("EDIFF=5E-4\n")
            incar.write("EDIFFG=5E-4\n")
            incar.write("ISMEAR=0\n")
            incar.write("ISPIN=2\n")
            incar.write("LREAL=Auto\n")
            incar.write("SIGMA=0.1\n")
            incar.write("LWAVE=.TRUE.\n")
            incar.write("LCHARG=.TRUE.\n")
    
            incar.close()

    def makeNormalINCAR(self):
        """ Creates a standard INCAR file and puts it in each different structure's top directory. """
        
        incar = open('normal/INCAR','w')
    
        incar.write("IBRION=2\n")
        incar.write("ISIF=4\n")
        incar.write("NSW=400\n")
        incar.write("Algo=VeryFast\n")
        incar.write("PREC=Normal\n")
        incar.write("EDIFF=5E-4\n")
        incar.write("EDIFFG=5E-4\n")
        incar.write("ISMEAR=0\n")
        incar.write("ISPIN=2\n")
        incar.write("LREAL=Auto\n")
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")
    
        incar.close()

    def makeDOS_INCAR(self):
        """ Creates an INCAR file for the Density of States run in VASP.  The notable changes are:
                IBRION=-1 -- This tells VASP not to move the ions.
                NSW=0     -- This tells VASP that there will be no ionic relaxation steps.
                LORBIT=10 -- This creates the PROCAR file which can be used to project onto the C, H, and M atoms. """
                
        incar = open('DOS/INCAR','w')
        
        incar.write("IBRION=-1\n")
        incar.write("ISIF=4\n")
        incar.write("NSW=0\n")
        incar.write("Algo=Normal\n")
        incar.write("PREC=Normal\n")
        incar.write("EDIFF=5E-4\n")
        incar.write("EDIFFG=5E-4\n")
        incar.write("ISMEAR=0\n")
        incar.write("ISPIN=2\n")
        incar.write("LREAL=Auto\n")
        incar.write("LORBIT=10\n")
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")
    
        incar.close()
        
    def makeKPOINTS(self, num1, num2):
        """ Creates a KPOINTS file based on the input parameters num1 and num2. It specifies that the job will have
        num1 x num2 kpoints. For example, if we wanted to specify an 8x8 kpoints scenario, we would call
        makeKPOINTS(8, 8). """
        
        dirList = self.atomList
        
        for direc in dirList:
            kpoints = open(direc + '/KPOINTS','w')
    
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write(str(num1) + ' ' + str(num2) + ' 1\n')
            kpoints.write('0 0 0')
    
            kpoints.close()
           
    def makePOTCARs(self):
        """ Creates a POTCAR file based on the atoms in the given input list. Concatenates the individual POTCAR
        files to make a single POTCAR file for the multi-atom structure. It is assumed that the atoms in the
        list are given in the correct order. """
    
        for atom in self.atomList:
            CPotcarDir = "/fslhome/eswens13/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/C/POTCAR"
            HPotcarDir = "/fslhome/eswens13/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/H/POTCAR"
            atomPotcarDir = "/fslhome/eswens13/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + "/POTCAR"
            if os.path.exists(atomPotcarDir):
                CPotcar = open(CPotcarDir, 'r')
                CLines = CPotcar.readlines()
                CPotcar.close()
                    
                HPotcar = open(HPotcarDir, 'r')
                HLines = HPotcar.readlines()
                HPotcar.close()
                    
                atomPotcar = open(atomPotcarDir,'r')
                atomLines = atomPotcar.readlines()
                atomPotcar.close()
            
                potcar = open(atom + '/POTCAR', 'w')
                for line in CLines:
                    potcar.write(line)
                        
                for line in HLines:
                    potcar.write(line)
                    
                for line in atomLines:
                    potcar.write(line)
                        
                potcar.close()
                    
            else:
                print "ERROR: Could not find a POTCAR file for \'" + atom + "\'"
                print "Removing POTCAR . . ."
                potcar.close()
                subprocess.call(['rm','POTCAR'])
                return
    
    def makePurePOTCARs(self):
        """ Some of the structures that need to be submitted to VASP for relaxation are what we call
            "pure" structures.  This means that (other than carbon atoms) the structure only contains 
            one other type of atom.  In the binary representation of the structure, this means that it
            is either all '1's or all '0's.  VASP gives a segmentation fault even if we tell it that 
            there are zero of one of the kinds of atoms.  It needs a POTCAR file that doesn't even mention
            the element that is not a part of the "pure" structure.  This method creates these POTCAR 
            files. """
        
        for atom in self.atomList:
            atomPotcarDir = "/fslhome/eswens13/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom
            
            if os.path.isdir(atomPotcarDir):
                purePotcar = open(atom + "/C" + atom + "_POTCAR", "w")
                
                CPotcar = open("/fslhome/eswens13/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/C/POTCAR", "r")
                CLines = CPotcar.readlines()
                CPotcar.close()
                
                atomPotcar = open(atomPotcarDir + "/POTCAR", "r")
                atomLines = atomPotcar.readlines()
                atomPotcar.close()
                
                for line in CLines:
                    purePotcar.write(line)
                
                for line in atomLines:
                    purePotcar.write(line)
                
                purePotcar.close()
                
                CHPotcar = open(atom + "/CH_POTCAR", "w")
                
                HPotcar = open("/fslhome/eswens13/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/H/POTCAR", "r")
                HLines = HPotcar.readlines()
                HPotcar.close()
                
                for line in CLines:
                    CHPotcar.write(line)
                
                for line in HLines:
                    CHPotcar.write(line)
                
                CHPotcar.close()
                                            
    def makeJobFiles(self):
        """ Creates a standard job file for submitting a VASP job on the supercomputer. """
    
        dirList = self.atomList
        
        for direc in dirList:
            jobFile = open(direc + '/job','w')
    
            jobFile.write("#!/bin/bash\n\n")
            jobFile.write("#SBATCH --time=06:00:00\n")
            jobFile.write("#SBATCH --ntasks=16\n")
            jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
            jobFile.write("#SBATCH --mail-user=erandemswens@gmail.com\n")
            jobFile.write("#SBATCH --mail-type=FAIL\n")
            jobFile.write("\nmpiexec vasp533 > vasp.out\n")
    
            jobFile.close()
    
    def makeDOSJobFile(self):
        
        jobFile = open('DOS/job','w')
        
        jobFile.write("#!/bin/bash\n\n")
        jobFile.write("#SBATCH --time=06:00:00\n")
        jobFile.write("#SBATCH --ntasks=16\n")
        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
        jobFile.write("#SBATCH --mail-user=erandemswens@gmail.com\n")
        jobFile.write("#SBATCH --mail-type=FAIL\n")
        jobFile.write("\nmpiexec vasp533 > vasp.out\n")
    
        jobFile.close()

    def copyVaspExec(self):
        """ Copies the vasp executable file to the current working directory. """
    
        dirList = self.atomList
        
        for direc in dirList:
            subprocess.call(['cp','/fslhome/eswens13/bin/vasp',direc + '/vasp533'])

    def fillDirectories(self):
        """ Fills all the directories with the needed files for VASP to run, namely POSCAR, POTCAR, KPOINTS, 
            INCAR, job, and the VASP executable file. """
            
        for atom in self.atomList:
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + atom
            
            os.chdir(atomDir)
            contents = os.listdir(atomDir)
            structures = []
            for item in contents:
                if os.path.isdir(item):
                    structures.append(item)
            
            for structure in structures:
                structureDir = os.path.abspath(structure)
                subprocess.call(['cp', 'KPOINTS', 'INCAR', 'job', 'vasp533', structureDir])
                
                poscar = open(structureDir + '/POSCAR','r')
                poscarLines = [line.strip() for line in poscar]
                poscar.close()
                
                if poscarLines[0].split()[1] == 'H':
                    subprocess.call(['cp','CH_POTCAR',structureDir + '/POTCAR'])
                elif poscarLines[0].split()[1] == 'M':
                    subprocess.call(['cp','C' + atom + '_POTCAR',structureDir + '/POTCAR'])
                else:
                    subprocess.call(['cp','POTCAR',structureDir])
            
            os.chdir(lastDir)

    def startJobs(self):
        """ Submits all the VASP jobs to the supercomputer for initial ionic relaxation. """
        
        dirList = self.atomList
        self.clearCurrentJobIds()
        
        for direc in dirList:
            lastDir = os.getcwd()
            newDir = lastDir + '/' + direc
            
            contents = os.listdir(newDir)
            os.chdir(newDir)
            
            structures = []
            for item in contents:
                if os.path.isdir(item):
                    structures.append(item)
            
            for structure in structures:
                os.chdir(structure)
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
                print "Submitted job " + jobid
                self.currJobIds.append(jobid)
                os.chdir(newDir)
            
            os.chdir(lastDir)

    def startNormalJobs(self):
        dirList = self.atomList
        self.clearCurrentJobIds()
        
        for direc in dirList:
            lastDir = os.getcwd()
            newDir = lastDir + '/' + direc
            
            contents = os.listdir(newDir)
            os.chdir(newDir)
            
            structures = []
            for item in contents:
                if os.path.isdir(item + '/normal'):
                    structures.append(item)
            
            for structure in structures:
                os.chdir(structure + '/normal')
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
                print "Submitted job " + jobid
                self.currJobIds.append(jobid)
                os.chdir(newDir)
            
            os.chdir(lastDir)
 
    def startDOSJobs(self):
        """ Submits all the VASP jobs to the supercomputer for the Density of States run. """
        
        topDir = os.getcwd()
        self.clearCurrentJobIds()
        
        for element in self.atomList:
            elementDir = topDir + '/' + element
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                
                for structure in os.listdir(os.getcwd()):
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        os.chdir(structDir)
                        
                        dosDir = structDir + '/DOS'
                        if os.path.isdir(dosDir):
                            os.chdir(dosDir)
                            proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                            jobid = proc.communicate()[0].split()[3]
                            print "Submitted job " + jobid
                            self.currJobIds.append(jobid)
                        
                        os.chdir(structDir)
                    
                    os.chdir(elementDir)
            else:
                print "The directory " + elementDir + " does not exist."
            
            os.chdir(topDir)
    
    def clearCurrentJobIds(self):
        self.currJobIds = []
    
    def getCurrentJobIds(self):
        return self.currJobIds
    
    def prepareForVasp(self):
        self.makeLowINCARs()
        self.makePurePOTCARs()
        self.makePOTCARs()
        self.makeKPOINTS(8, 8)
        self.makeJobFiles()
        self.copyVaspExec()
        
        self.fillDirectories()
    
    def run(self, runNum):       
        if runNum == 1:
            self.startJobs()
    
        elif runNum == 2:
            self.makeNormalDirectories()
            self.startNormalJobs()
           
        elif runNum == 3:
            self.makeDOSDirectories()
            self.startDOSJobs()








        
        