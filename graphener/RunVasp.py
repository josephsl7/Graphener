import os, subprocess


class RunVasp:
    """ This class is responsible for preparing the directories, retrieving the needed files, and
        submitting VASP jobs to the supercomputer.  It keeps track of the SLURM job ids of all the
        jobs that are currently running from a particular instance of the class. """
    
    from comMethods import setAtomCounts
    def __init__(self, atoms,ediffg):
        """ CONSTRUCTOR """
        
        
        self.atoms = atoms
        self.neededFilesFolder = os.getcwd() + '/needed_files/'
        
        self.currJobIds = []
        self.ediffg = ediffg 
        

#    def addStructName(self,nameadd):
#        jobfile = open('job','r')
#        lines = jobfile.readlines()
#        jobfile.close()
#        jobfile = open('job','w')
#        for line in lines:
#            if 'job-name' in line: 
#                jobfile.write(line.strip('\n') + '_'+ nameadd + '\n')
#            else:
#                jobfile.write(line)
#        jobfile.close()

    def clearCurrentJobIds(self):
        """ Clears the list of current job IDs. """
        self.currJobIds = []
      
    def linkVaspExec(self):
        """ Copies the vasp executable file to the current working directory. """
        dirList = self.atoms
        for direc in dirList:
            subprocess.call(['ln','-s','/fslhome/bch/bin/vasp533',direc + '/vasp533']) 

    def fillDirectories(self, vstructsToStart):
        """ Fills all the directories in 'vstructsToStart' with the needed files for VASP to run, namely
            POSCAR, POTCAR, KPOINTS, INCAR, a SLURM job file, and the VASP executable file. """
        for iatom,atom in enumerate(self.atoms):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + atom
            
            os.chdir(atomDir)
            structs = []
            for item in vstructsToStart[iatom]:
                if os.path.isdir(item):
                    structs.append(item)
            
            for struct in structs:
                structDir = os.path.abspath(struct)
                self.makeJobFiles(structDir,atom+struct)  
                self.makeLowINCARs(structDir)  
                subprocess.call(['cp','-P','KPOINTS','vasp533', structDir])
                poscar = open(structDir + '/POSCAR','r')
                poscarLines = [line.strip() for line in poscar]
                poscar.close()
                
                if poscarLines[0].split()[1] == 'H':
                    subprocess.call(['cp','CH_POTCAR',structDir + '/POTCAR'])
                elif poscarLines[0].split()[1] == 'M':
                    subprocess.call(['cp','C' + atom + '_POTCAR',structDir + '/POTCAR'])
                else:
                    subprocess.call(['cp','POTCAR',structDir])            
            os.chdir(lastDir)

    def getCurrentJobIds(self):
        """ Returns the list of current job IDs. """
        return self.currJobIds  

    def makeDOSDirectories(self, vstructsToStart):
        """ After the normal-precision relaxation, creates a directory for the Density of States
            run and populates it with the files from the normal-precision run. Copies the normal
            CONTCAR to the DOS POSCAR. """  
        topDir = os.getcwd()
        for iatom,atom in enumerate(self.atoms):
            elementDir = topDir + '/' + atom
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in vstructsToStart[iatom]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        normalDir = structDir + '/normal'
                        if os.path.isdir(normalDir) and finishCheck(normalDir) and convergeCheck(normalDir, getNSW(normalDir)) and energyDropCheck(normalDir):
                            os.chdir(structDir)
                            subprocess.call(['mkdir', 'DOS'])
                            subprocess.call(['cp','normal/CONTCAR','normal/DOSCAR','normal/EIGENVAL',
                                             'normal/IBZKPT','normal/KPOINTS','normal/vasp533',
                                             'normal/OSZICAR','normal/OUTCAR','normal/PCDAT',
                                             'normal/POSCAR','normal/POTCAR','normal/REPORT',
                                             'normal/vasprun.xml','normal/XDATCAR','DOS'])
                            self.makeDOS_INCAR()
                            self.makeDOSJobFile(atom+structure) #bch 
                            subprocess.call(['cp','DOS/CONTCAR','DOS/POSCAR'])
                            os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir)

    def makeDOS_INCAR(self):
        """ Creates an INCAR file for the Density of States run in VASP.  The notable changes are:
                IBRION=-1 -- This tells VASP not to move the ions.
                NSW=0     -- This tells VASP that there will be no ionic relaxation steps.
                LORBIT=10 -- This creates the PROCAR file which can be used to project onto the 
                             C, H, and M atoms. """
                
        incar = open('DOS/INCAR','w')
        
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

    def makeDOSJobFile(self,name):
        
        jobFile = open('DOS/job','w')
        
        jobFile.write("#!/bin/bash\n\n")
        jobFile.write("#SBATCH --time=06:00:00\n")
        jobFile.write("#SBATCH --ntasks=16\n")
        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
        jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
        jobFile.write("#SBATCH --mail-type=END\n")  
        jobFile.write("#SBATCH --job-name=%s\n" % name)
        jobFile.write("#SBATCH --mail-type=FAIL\n")
        jobFile.write("\nmpiexec vasp533 > vasp.out\n")
    
        jobFile.close()          

    def makeJobFiles(self,dir,name):
        """ Creates a standard job file for submitting a VASP job on the supercomputer. 
        Done one structure at a time so we can have the structure in the name. """
        jobFile = open(dir + '/job','w')   
        jobFile.write("#!/bin/bash\n\n")
        jobFile.write("#SBATCH --time=24:00:00\n")
#            jobFile.write("#SBATCH --time=00:00:30\n")
        jobFile.write("#SBATCH --ntasks=16\n")
        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
        jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")              
        jobFile.write("#SBATCH --mail-type=FAIL\n")
        jobFile.write("#SBATCH --mail-type=END\n") 
        jobFile.write("#SBATCH --job-name=%s\n\n" % name) 
        jobFile.write("mpiexec vasp533 > vasp.out\n")    
        jobFile.close()

    def makeKPOINTS(self, num1, num2):
        """ Creates a KPOINTS file based on the input parameters num1 and num2. It specifies that 
            the job will have num1 x num2 kpoints. For example, if we wanted to specify an 8x8 
            kpoints mesh, we would call makeKPOINTS(8, 8). """
        dirList = self.atoms
        
        for direc in dirList:
            kpoints = open(direc + '/KPOINTS','w')
    
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write(str(num1) + ' ' + str(num2) + ' 1\n')
            kpoints.write('0 0 0')
    
            kpoints.close()


    def makeLowINCARs(self,dir):
        """ Creates a standard INCAR file and puts it in each different structure's top 
            directory. Ediffs increase with atom count, since it's a measure
            for total energy change, not per atom"""
        self.setAtomCounts(dir)
        natoms = sum(self.atomCounts)
        incar = open(dir + '/INCAR','w')    
        incar.write("IBRION=2\n")
        incar.write("ISIF=4\n")
        incar.write("NSW=400\n")
        incar.write("Algo=VeryFast\n")
        incar.write("PREC=Low\n") 
        incar.write("EDIFF={}\n".format(float(self.ediffg)*natoms))
        incar.write("EDIFFG={}\n".format(float(self.ediffg)*natoms))
        incar.write("ISMEAR=0\n")
        incar.write("ISPIN=2\n")
        incar.write("NPAR=4\n")
        incar.write("LREAL=Auto\n")
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")    
        incar.close()

    def makeNormalDirectories(self, vstructsToStart):
        topDir = os.getcwd()
        for iatom,atom in enumerate(self.atoms):
            elementDir = topDir + '/' + atom
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in vstructsToStart[iatom]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir) and finishCheck(structDir) and convergeCheck(structDir, getNSW(structDir)) and energyDropCheck(normalDir):
                        subprocess.call(['mkdir', 'normal'])
                        subprocess.call(['cp','-P','CONTCAR','KPOINTS','vasp533',
                                         'POTCAR','job','normal'])
                        subprocess.call(['mv','CHG','CHGCAR','WAVECAR','normal'])                 
                        self.makeNormalINCAR()
                        subprocess.call(['cp','normal/CONTCAR','normal/POSCAR'])
                        os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir) 

    def makeNormalINCAR(self):
        """ Creates a standard INCAR file for normal-precision relaxation. """ 
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
        incar.write("NPAR=4\n")
        incar.write("LREAL=Auto\n")
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")   
        incar.close()


    def makePOTCARs(self):
        """ Creates a POTCAR file for each atom in the member 'atoms'. Concatenates the 
            individual POTCAR files to make a single POTCAR file for the multi-atom structure. """
        for atom in self.atoms:
            CPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/C/POTCAR"
            HPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/H/POTCAR"
            atomPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + "/POTCAR"
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
                subprocess.call(['echo','ERROR: Could not find a POTCAR file for \'' + atom + '\''])
                subprocess.call(['echo','Removing POTCAR . . .'])
                potcar.close()
                subprocess.call(['rm','POTCAR'])
                return

    def makePurePOTCARs(self):
        """ Some of the structures that need to be submitted to VASP for relaxation are what we 
            call "pure" structures.  This means that (other than carbon atoms) the structure only 
            contains one other type of atom.  In the binary representation of the structure, this 
            means that it is either all '1's or all '0's.  VASP gives a segmentation fault when we 
            try to run a pure structure with a POTCAR that contains atoms that are not in the pure 
            structure, even if we tell it that there are zero of one of the kinds of atoms.  It 
            needs a POTCAR file that doesn't even mention the element that is not a part of the 
            "pure" structure. This method creates these POTCAR files. """       
        for atom in self.atoms:
            atomPotcarDir = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom
            purePotcar = open(atom + "/C" + atom + "_POTCAR", "w")                
            CPotcar = open("/fslhome/bch/hessgroup/vaspfiles/src/potpaw_PBE/C/POTCAR", "r")
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
            HPotcar = open("/fslhome/bch/hessgroup/vaspfiles/src/potpaw_PBE/H/POTCAR", "r")
            HLines = HPotcar.readlines()
            HPotcar.close()           
            for line in CLines:
                CHPotcar.write(line)           
            for line in HLines:
                CHPotcar.write(line)            
            CHPotcar.close()  
        
    def makeRunHexMono(self): #bch all
        topDir = os.getcwd()
        if not os.path.isdir('hex_monolayer_refs'): os.mkdir('hex_monolayer_refs')
        os.chdir('hex_monolayer_refs')
#        os.system('rm -r -f *')
        for atom in self.atoms:
            os.mkdir(atom)
            os.chdir(atom)
            atomPotcar = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + '/POTCAR'
            if os.path.exists(atomPotcar):
                subprocess.call(['cp', atomPotcar, '.'])
            else:
                system.exit('Failed to read POTCAR in makeSingleDirectories')                
            jobFile = open('job','w')
            jobFile.write("#!/bin/bash\n\n")
            jobFile.write("#SBATCH --time=03:00:00\n")
            jobFile.write("#SBATCH --ntasks=1\n")
            jobFile.write("#SBATCH --mem-per-cpu=4G\n")
            jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")              
            jobFile.write("#SBATCH --mail-type=FAIL\n")
            jobFile.write("#SBATCH --mail-type=END\n") 
            jobFile.write("#SBATCH --job-name=hexm%s\n" % atom)           
            jobFile.write("\nmpiexec vasp533 > vasp.out\n") 
            jobFile.close()
            incar = open('INCAR','w')
            incar.write("IBRION=2\n")
            incar.write("ISIF=4\n")
            incar.write("NSW=400\n")
            incar.write("PREC=High\n")
            incar.write("EDIFF=1E-6\n")
            incar.write("ISPIN=2\n")
            incar.write("LWAVE=.FALSE.\n")
            incar.write("LCHARG=.FALSE.\n")         
            incar.close()            
            kpoints = open('KPOINTS','w')
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write('8 8 8\n')
            kpoints.write('0 0 0')
            kpoints.close()
            poscar = open('POSCAR','w') 
            poscar.write('Hexagonal metal monolayer\n')
            poscar.write('1.0\n') #scale same as graphene for now
            poscar.write('2.13128850  -1.23050000  0\n') 
            poscar.write('2.13128850   1.23050000 0\n')
            poscar.write('0.00000000   0.00000000  15.00000000\n')
            poscar.write('1\n')
            poscar.write('Cartesian\n')
            poscar.write('0 0 0\n')
            poscar.close()
            subprocess.call(['sbatch','job'])
            os.chdir('../')      
        os.chdir(topDir)   
    
    def makeRunSingleDirectories(self): #bch all
        topDir = os.getcwd()
        if not os.path.isdir('single_atoms'): os.mkdir('single_atoms')
        os.chdir('single_atoms')
#        os.system('rm -r -f *')
        for atom in self.atoms:
            os.mkdir(atom)
            os.chdir(atom)
            atomPotcar = "/fslhome/bch/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + '/POTCAR'
            if os.path.exists(atomPotcar):
                subprocess.call(['cp', atomPotcar, '.'])
            else:
                system.exit('Failed to read POTCAR in makeSingleDirectories')                
            jobFile = open('job','w')
            jobFile.write("#!/bin/bash\n\n")
            jobFile.write("#SBATCH --time=03:00:00\n")
            jobFile.write("#SBATCH --ntasks=1\n")
            jobFile.write("#SBATCH --mem-per-cpu=4G\n")
            jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")              
            jobFile.write("#SBATCH --mail-type=FAIL\n")
            jobFile.write("#SBATCH --mail-type=END\n") 
            jobFile.write("#SBATCH --job-name=isol_%s\n" % atom)           
            jobFile.write("\nmpiexec vasp533 > vasp.out\n") 
            jobFile.close()
            incar = open('INCAR','w')
            incar.write("IBRION=-1\n")
            incar.write("NELM=400\n")
            incar.write("PREC=High\n")
            incar.write("EDIFF=1E-6\n")
            incar.write("ISPIN=2\n")
            incar.write("LWAVE=.FALSE.\n")
            incar.write("LCHARG=.FALSE.\n")         
            incar.close()            
            kpoints = open('KPOINTS','w')
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write('1 1 1\n')
            kpoints.write('0 0 0')
            kpoints.close()
            poscar = open('POSCAR','w') 
            poscar.write('Isolated atom\n')
            poscar.write('1\n')
            poscar.write('20 0 0\n') 
            poscar.write('0 20 0\n')
            poscar.write('0.0 0.0 20.0\n')
            poscar.write('1\n')
            poscar.write('Cartesian\n')
            poscar.write('0 0 0\n')
            poscar.close()
            subprocess.call(['sbatch','job'])
            os.chdir('../')      
        os.chdir(topDir)      
    
    def prepareForVasp(self, vstructsToStart):
        """ Makes all of the files that could be copied to a first, low-precision VASP run for any 
            given structure.  This includes concatenating the POTCARS for the pure and non-pure
            cases. """
        self.makePurePOTCARs()
        self.makePOTCARs()
        self.makeKPOINTS(6, 6)
        self.linkVaspExec()
        self.fillDirectories(vstructsToStart)

    def prepareRestarts(self, vstructsRestart):
        """ INCAR settings and job walltime updates might be needed"""
        lastDir = os.getcwd()
        for iatom,atom in enumerate(self.atoms):
            atomDir = lastDir + '/' + atom
            os.chdir(atomDir)
            structures = []
            for item in vstructsRestart[iatom]:
                if os.path.isdir(item):
                    structures.append(item)            
            for structure in structures:
                structureDir = os.path.abspath(structure)
                self.makeLowINCARs(structureDir)             
                self.makeJobFiles(structureDir,atom+structure)    
        os.chdir(lastDir)
        
    def run(self, runNum, vstructsToStart,vstructsToRun):
        """ Starts the VASP runs (specified by 'runNum') for each of the structures in
            'vstructsToRun'. For runNum = 1, starts a low-precision run, runNum = 2, starts a 
            normal-precision run, runNum = 3 starts a DOS run. """
        if runNum == 1:
            self.startJobs(vstructsToRun)
    
        elif runNum == 2:
            self.makeNormalDirectories(vstructsToStart)
            self.startNormalJobs(vstructsToRun)
           
        elif runNum == 3:
            self.makeDOSDirectories(vstructsToStart)
            self.startDOSJobs(vstructsToRun)

    def startDOSJobs(self, vstructsToStart):
        """ Submits all the VASP jobs for structures in 'vstructsToStart' to the supercomputer for 
            Density of States calculations. Records their SLURM job IDs. """
        topDir = os.getcwd()
        self.clearCurrentJobIds()
        
        for iatom,atom in enumerate(self.atoms):
            elementDir = topDir + '/' + atom
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                
                for structure in vstructsToStart[iatom]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        os.chdir(structDir)
                        
                        dosDir = structDir + '/DOS'
                        if os.path.isdir(dosDir):
                            os.chdir(dosDir)
                            proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                            jobid = proc.communicate()[0].split()[3]
                            subprocess.call(['echo','Submitted job ' + jobid])
                            self.currJobIds.append(jobid)
                        
                        os.chdir(structDir)
                    
                    os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir)

    def startJobs(self, vstructsToRun):
        """ Submits all the VASP jobs for structures in 'vstructsToRun' to the supercomputer for 
            low-precision relaxation and records their job IDs. """
        self.clearCurrentJobIds()
        for iatom,atom in enumerate(self.atoms):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + atom
            
            os.chdir(atomDir)
            
            structures = []
            for item in vstructsToRun[iatom]:
                if os.path.isdir(item):
                    structures.append(item)
            
            for structure in structures:
                os.chdir(structure)
                #self.addStructName(structure) #bch adds structure to name
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
#                subprocess.call(['echo', 'Submitted job ' + jobid])
                self.currJobIds.append(jobid)
                os.chdir(atomDir)
            os.chdir(lastDir)
        if len(self.currJobIds)>0:
            subprocess.call(['echo', 'Submitted {}jobs, ranging from ID {} to ID {}.'.format(len(self.currJobIds),self.currJobIds[0],self.currJobIds[-1])])            

    def startNormalJobs(self, vstructsToRun):
        """ Submits all the VASP jobs for structures in 'vstructsToRun' to the supercomputer for 
            normal-precision relaxation and records their job IDs. """
        self.clearCurrentJobIds()
        
        for iatom,atom in enumerate(self.atoms):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + atom
            
            os.chdir(atomDir)
            
            structures = []
            for item in vstructsToRun[iatom]:
                if os.path.isdir(item + '/normal'):
                    structures.append(item)
            
            for structure in structures:
                os.chdir(structure + '/normal')
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
                subprocess.call(['echo','Submitted job ' + jobid])
                self.currJobIds.append(jobid)
                os.chdir(atomDir)
            
            os.chdir(lastDir)




        









