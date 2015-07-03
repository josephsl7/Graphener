'''
'''

from numpy import sqrt,log10,array, mod,dot, transpose, std, size, zeros, where, int32, mean,square, sign,argmax,argmin
from numpy.linalg.linalg import inv, norm, det
import os, subprocess,sys
from copy import deepcopy
from comMethods import readfile,writefile,trimSmall,rms,finishCheck,convergeCheck,\
    energyDropCheck,getNSW
            
import StructsToPoscar
import matplotlib.pyplot as plt

#class ChargeDiff:
#    '''Runs jobs with all the C atoms removed and also with all the adatoms removed.  
#    Calculates the difference between the sum of those CHG files and the normal run'''
#    from comMethods import setAtomCounts,startJobs,getEnergy, elConvergeCheck, getElSteps,collatePlotsGSS
#    def __init__(self, atoms):
#        self.atoms = atoms
#        self.nC = 0 
#        self.nH = 0 
#        self.nM = 0
#        self.jobIds = []
# 
#    def createFolders(self,atom,atomDir,walltime,type):
#        '''Copies from previous run folders, but modifies the POSCARs with deltaz for adatoms, and changes INCAR so atoms don't relax'''
#        topDir = os.getcwd()
#        structlist = []
#        toStart = []
#        for item in os.listdir(atomDir):
#            itempath = atomDir + '/' + item
#            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
#                structlist.append(item)
#        sstring = ''; 
#        for struct in structlist: 
#            sstring += struct + ' '
#        subprocess.call(['echo','Found structures '+ sstring])
#        analysisDir = '{}/analysis'.format(atomDir)       
#        for istruct,struct in enumerate(structlist):
#            if mod(istruct+1,100) == 0 or istruct+1 ==len(structlist): subprocess.call(['echo','\tchgDiff prep for {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
#            structDir = atomDir + '/' + struct
#            planarDir = atomDir + '/' + struct + '/chgDiff_{}'.format(type)
#            if not os.path.exists(structDir):
#                subprocess.call(['echo','No vasp folder for {}, struct {}'.format(atom,struct)])
#            elif struct != '1' and finishCheck(structDir) and  convergeCheck(structDir, getNSW(structDir)) and energyDropCheck(structDir):
#                #remove this if statement when you want to start over on chgDiff calcs. 
##                if not os.path.exists(planarDir) or not finishCheck(planarDir): #not doing electronic convergence check for this precision...good within 1-5%% of BE
#                    if os.path.exists(planarDir): subprocess.call(['rm','-r',planarDir]) #start over for now
#                    os.mkdir(planarDir)
#                    subprocess.call(['ln','-s','/fslhome/bch/bin/vasp533',planarDir+'/vasp533'])
#                    subprocess.call(['cp',structDir+'/CONTCAR',planarDir+ '/POSCAR'])
#                    subprocess.call(['cp',structDir+'/INCAR',planarDir])
##                    subprocess.call(['ln','-s',structDir+'/CHG',planarDir]) #Use only for perturbation (force)
##                    subprocess.call(['ln','-s',structDir+'/CHGCAR',planarDir]) #Use only for perturbation (force)
##                    subprocess.call(['ln','-s',structDir+'/WAVECAR',planarDir]) #Use only for perturbation (force)
#                    subprocess.call(['ln','-s',structDir+'/POTCAR',planarDir])
#                    subprocess.call(['ln','-s',structDir+'/KPOINTS',planarDir])
#                    self.makePlanarJobFile(planarDir,atom+struct,walltime)
#                    self.modifyINCAR(planarDir)
#                    self.modifyPOSCAR1Atom(planarDir,self.getdMCatomIndex(struct,type,analysisDir))
#                    toStart.append(struct)
#            else:
#                print 
#                subprocess.call(['echo','Structure ' + struct + ' has not finished, converged or has an energy rise.'])
#            os.chdir(topDir)
#        return toStart            
#
#    def electronicConvergeFinish(self,dir): 
#        '''Test requires electronic convergence AND vasp finishing'''
##        get NELM, the max electronic steps allowed in the run. 
#        proc = subprocess.Popen(['grep','-i','NELM',dir+'/INCAR'],stdout=subprocess.PIPE)
#        result =  proc.communicate()[0]
#        NELM = int(result.split('=')[1].split()[0])
#        return self.elConvergeCheck(dir,NELM) and finishCheck(dir)
#
#    def getChgDiff(self, atom, atomDir, type):
#        ''' Calculates the difference between the sum of those CHG files and the normal run.
#        Saves this difference file for each struct as 'chgDiff'.   Also calculates chgDiffavg, the sum(abs(chgDiff))/Ngrid a
#        global measure of the charge changes due to the adatoms. '''    
#        topDir = os.getcwd()
#        structlist = []
#        toStart = []
#        for item in os.listdir(atomDir):
#            itempath = atomDir + '/' + item
#            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
#                structlist.append(item)
#        sstring = ''; 
#        for struct in structlist: 
#            sstring += struct + ' '
#        subprocess.call(['echo','Found structures '+ sstring]) 
#        file = open(atomDir + '/analysis/chgDiff_dz_{}_{}'.format(self.deltaz,type),'w' )      
#        for istruct,struct in enumerate(sorted(structlist)):
#            if mod(istruct+1,100) == 0 or istruct+1 == len(structlist): subprocess.call(['echo','\tChgDiff analyzed for {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
#            structDir = atomDir + '/' + struct
#            planarDir = atomDir + '/' + struct + '/chgDiff_{}'.format(type)
#            if not os.path.exists(structDir):
#                subprocess.call(['echo','No vasp folder for {}, struct {}'.format(atom,struct)])
#            elif not os.path.exists(planarDir):
#                subprocess.call(['echo','No planar BE folder for {}, struct {}'.format(atom,struct)]) 
#            elif finishCheck(planarDir) and  os.path.getsize(structDir + '/chgDiff/POSCAR') > 0:                
#                self.setAtomCounts(planarDir)
#                BE = (self.getEnergy(structDir) - self.getEnergy(planarDir))/sum(self.atomCounts[1:])#/self.deltaz
#                file.write('{} {:8.4f}\n'.format(struct,BE))
#            else:
#                subprocess.call(['echo','Structure ' + struct + ' has not finished the planar BE vasp calculation.'])
#            os.chdir(topDir)
#        file.close()
#
#    def makeChgDiffJobFile(self,planarDir,name,walltime):
#        jobFile = open(planarDir +'/job','w')
#        jobFile.write("#!/bin/bash\n\n")
#        jobFile.write("#SBATCH --time={}:00:00\n".format(walltime))
#        jobFile.write("#SBATCH --ntasks=16\n")
#        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
#        jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
#        jobFile.write("#SBATCH --mail-type=END\n")
#        jobFile.write("#SBATCH --mail-type=FAIL\n")  
#        jobFile.write("#SBATCH --job-name=PBE%s\n" % name)
#        jobFile.write("\nmpiexec vasp533 > vasp.out\n")   
#        jobFile.close() 
#
#    def modifyINCAR(self,planarDir):
#        '''Changes it so it has no ionic relaxation.  Ignore CHG and WAVECAR files'''
#        lines = readfile(planarDir + '/INCAR')
#        lines2 = deepcopy(lines)
#        for i,line in enumerate(lines):
#            if 'IBRION' in line:
#                lines2[i] = 'IBRION=2\n'  #testing with a Cr top structure 2, find we must have IBRION 2 and NSW 1.  NOT -1 and 0 as expected for static run!!!
#            elif 'NSW' in line:
#                lines2[i] = 'NSW=1\n'  #testing with a Cr top structure 2, find we must have IBRION 2 and NSW 1.  NOT -1 and 0 as expected for static run!!!
#            elif 'LCHARG' in line:
#                lines2[i] = 'LCHARG=.FALSE.\n' 
#            elif 'LWAV' in line:
#                lines2[i] = 'LWAV=.FALSE.\n'
#        writefile(lines2,planarDir + '/INCAR') 
#                          
##    def modifyPOSCAR(self,planarDir):
##        '''all adatoms are displaced deltaz away from the plane.  Since this POSCAR
##        came from a CONTAR, we are working in direct coordinates'''
##        lines = readfile(planarDir + '/POSCAR')
##        lines2 = deepcopy(lines)
##        lattvecs = lines[2:5]
##        lattvecs = [line.strip().split() for line in lattvecs]
##        zLV = max([abs(float(lattvecs[i][2])) for i in range(3)]) #z repeat distance in ang 
##        izLV = argmax([abs(float(lattvecs[i][2])) for i in range(3)])                
##        self.setAtomCounts(planarDir)
##        self.nC = self.atomCounts[0] 
##        self.nH = self.atomCounts[1] 
##        self.nM = self.atomCounts[2]
##        self.nTot = sum(self.atomCounts)
###        zLV2 = 25
##        zLV2 = zLV #keep it unchanged for now
##        dz = self.deltaz/zLV2
###        lines2[2+izLV] = '{}   {}   {}\n'.format(lattvecs[izLV][0],lattvecs[izLV][1],zLV2)
##        for i in range(self.nTot):
##            pos = [float(part) for part in lines[8 + i].strip().split()]
##            if i < self.nC:
##                if pos[2] > 0.5: #don't want Vasp's convention to shift atoms to all positive coordinates
##                    pos[2] -= 1
##                lines2[8 + i] = '{}     {}     {}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV)
##            else:
##                if pos[2] > 0.5:  
##                    pos[2] -= 1
##                lines2[8 + i] = '{}     {}     {}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV + sign(pos[2])*dz)
##        writefile(lines2,planarDir + '/POSCAR')
#        
#        
#    def modifyPOSCAR1Atom(self,planarDir,index):
#            '''Moves one metal atom with either maximum or minimum (given by type)metal-carbon distance, by deltaz.  Since this POSCAR
#            came from a CONTAR, we are working in direct coordinates'''
#            lines = readfile(planarDir + '/POSCAR')
#            lines2 = deepcopy(lines)
#            lattvecs = lines[2:5]
#            lattvecs = [line.strip().split() for line in lattvecs]
#            zLV = max([abs(float(lattvecs[i][2])) for i in range(3)]) #z repeat distance in ang 
#            izLV = argmax([abs(float(lattvecs[i][2])) for i in range(3)])                
#            self.setAtomCounts(planarDir)
#            self.nC = self.atomCounts[0] 
#            self.nH = self.atomCounts[1] 
#            self.nM = self.atomCounts[2]
#            self.nTot = sum(self.atomCounts)
#    #        zLV2 = 25
#            zLV2 = zLV #keep it unchanged for now
#            dz = self.deltaz/zLV2
#    #        lines2[2+izLV] = '{}   {}   {}\n'.format(lattvecs[izLV][0],lattvecs[izLV][1],zLV2)
#            for i in range(self.nTot):
#                pos = [float(part) for part in lines[8 + i].strip().split()]
#                if i < self.nC: #just rewrite, no shift for these atoms
#                    if pos[2] > 0.5: #don't want Vasp's convention to shift atoms to all positive coordinates
#                        pos[2] -= 1
#                    lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV)
#                elif i-self.nH-self.nC == index: #this is the metal atom we want to move
#                    if pos[2] > 0.5:  
#                        pos[2] -= 1
#                    lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV + sign(pos[2])*dz)
#                else: #just rewrite, no shift for these atoms
#                    if pos[2] > 0.5:  
#                        pos[2] -= 1 
#                    lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV)
#            writefile(lines2,planarDir + '/POSCAR')
#        
#    def chgDiffruns(self, walltime):  
#        '''Script for this process'''
#        topDir = os.getcwd()
#        for iatom, atom in enumerate(self.atoms):
#            atomDir = topDir + '/' + atom       
#            toStart = self.createFolders(atom,atomDir,walltime,'max')
#            self.startJobs(toStart,atomDir,'chgDiff_max') 
#            toStart = self.createFolders(atom,atomDir,walltime,'min')
#            self.startJobs(toStart,atomDir,'chgDiff_min') 
#            
#    def chgDiffanalysis(self):
#        topDir = os.getcwd()
#        for iatom, atom in enumerate(self.atoms):
#            atomDir = topDir + '/' + atom       
#            self.getChgDiff(atom,atomDir,'max')
#            self.getChgDiff(atom,atomDir,'min')
#
#    def chgDiffplot(self):
#        topDir = os.getcwd()
#        for iatom, atom in enumerate(self.atoms):
#            atomDir = topDir + '/' + atom       
#            self.plotChgDiff(atom,atomDir,'max')
#            self.plotChgDiff(atom,atomDir,'min')                
#            self.collatePlotsGSS('HFE_chgDiff_max','xx')
#            self.collatePlotsGSS('HFE_chgDiff_min','xx')
#            
#    def plotChgDiff(self,atom,atomDir,type):
#        '''Plot the vaspHFE vs energy, and make the marker color reflect the planar BE'''
#        lines = readfile(atomDir + '/analysis/analysis.csv')
#        planarLines = readfile(atomDir + '/analysis/chgDiff_dz_{}_{}'.format(self.deltaz,type))
#        structs = []
#        conc = []
#        HFE = []
#        chgDiff = []
##        concPos = []; concNeg = []
##        HFEPos = []; HFENeg = []
##        chgDiffPos = []; chgDiffNeg = [
#        ikeep = -1
#        for ip,pline in enumerate(planarLines):
#            if mod(ip+1,100) == 0 or ip+1 ==len(planarLines): subprocess.call(['echo','\tChgDiff plot data read for {} of {} structures in {}'.format(ip+1,len(planarLines),atom)])            
#            pstruct = int(pline.split()[0])
#            BE = float(pline.split()[1])
#            
#            for line in lines[1:]:
#                struct = int(line.split(',')[0])
#                if struct == pstruct:
#                    if BE<0:
#                        conc.append(float(line.split(',')[2]))
#                        HFE.append(float(line.split(',')[3]))
#                        chgDiff.append(BE)
#                        break
#                    else:
#                        subprocess.call(['echo','Struct {} for atom {} has positive planar BE {:7.4f}...not included in plot\n'.format(struct,atom,BE)])  
##        chgDiff = log10(abs(array(chgDiff)))
#
#        #plot
#        fig = plt.figure()
#        plt.xlabel('{} concentration x'.format(atomDir.split('/')[-1]))
#        plt.ylabel('Energy (eV)')
##        plt.title('FE colored by planar force (eV/$\AA$) at {:2.1f} $\AA$'.format(self.deltaz/2.0))        
#        plt.title('FE colored by planar BE (eV), displaced {:2.1f} $\AA$'.format(self.deltaz))
##        plt.scatter(conc, HFE, c=info , cmap='autumn')
#        plt.scatter(conc, HFE, c=chgDiff , cmap='jet')
#        plt.colorbar()
#        plt.xlim(0, 1.0)
##        plt.ylim(plt.ylim()[0], min(max(HFE),2.0))
#        plt.ylim(min(HFE)-0.1, min(max(HFE),2.0))
#        plt.show()
#        fig.savefig(atomDir + '/gss/HFE_chgDiff_{}'.format(type)) #xx because iteration number is not used
#        plt.close(fig)   

class PlanarBE:
    '''Finds the binding energy of the adatoms to the carbon plane by 1) creating run folders 
    with POSCAR altered so all adatoms are displaced deltaz away from the plane, 
    2) running those jobs without ionic relaxation  3) calculating the difference (per adatom) in this energy
    and the fully relaxed energy.'''
    from comMethods import setAtomCounts,startJobs,getEnergy, elConvergeCheck, getElSteps,collatePlotsGSS
    def __init__(self, atoms):
        self.atoms = atoms
#        self.deltaz = 0.2  # perturbation   
        self.deltaz = 5.0
        self.newLVz = 25.0   
        self.nC = 0 
        self.nH = 0 
        self.nM = 0
        self.jobIds = []
 
    def createFolders(self,atom,atomDir,walltime):
        '''Copies from previous run folders, but modifies the POSCARs with deltaz for adatoms, and changes INCAR so atoms don't relax'''
        topDir = os.getcwd()
        structlist = []
        toStart = []
        for item in os.listdir(atomDir):
            itempath = atomDir + '/' + item
            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
                structlist.append(item)
        sstring = ''; 
        for struct in structlist: 
            sstring += struct + ' '
        subprocess.call(['echo','Found structures '+ sstring])
        analysisDir = '{}/analysis'.format(atomDir)       
        for istruct,struct in enumerate(structlist):
            if mod(istruct+1,100) == 0 or istruct+1 ==len(structlist): subprocess.call(['echo','\tPlanarBE prep for {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
            structDir = atomDir + '/' + struct
            planarDir = atomDir + '/' + struct + '/planarBE'
            if not os.path.exists(structDir):
                subprocess.call(['echo','No vasp folder for {}, struct {}'.format(atom,struct)])
            elif struct != '1' and finishCheck(structDir) and  convergeCheck(structDir, getNSW(structDir)) and energyDropCheck(structDir):
                #remove this if statement when you want to start over on planarBE calcs. 
#                if not os.path.exists(planarDir) or not finishCheck(planarDir): #not doing electronic convergence check for this precision...good within 1-5%% of BE
                    if os.path.exists(planarDir): subprocess.call(['rm','-r',planarDir]) #start over for now
                    os.mkdir(planarDir)
                    subprocess.call(['ln','-s','/fslhome/bch/bin/vasp533',planarDir+'/vasp533'])
                    subprocess.call(['cp',structDir+'/CONTCAR',planarDir+ '/POSCAR'])
                    subprocess.call(['cp',structDir+'/INCAR',planarDir])
#                    subprocess.call(['ln','-s',structDir+'/CHG',planarDir]) #Use only for perturbation (force)
#                    subprocess.call(['ln','-s',structDir+'/CHGCAR',planarDir]) #Use only for perturbation (force)
#                    subprocess.call(['ln','-s',structDir+'/WAVECAR',planarDir]) #Use only for perturbation (force)
                    subprocess.call(['ln','-s',structDir+'/POTCAR',planarDir])
                    subprocess.call(['ln','-s',structDir+'/KPOINTS',planarDir])
                    self.makePlanarJobFile(planarDir,atom+struct,walltime)
                    self.modifyINCAR(planarDir)
                    self.modifyPOSCAR(planarDir)
                    toStart.append(struct)
            elif struct == '1':
                'do nothing'
            else:
                print 
                subprocess.call(['echo','Structure ' + struct + ' has not finished, converged or has an energy rise.'])
            os.chdir(topDir)
        return toStart            

    def electronicConvergeFinish(self,dir): 
        '''Test requires electronic convergence AND vasp finishing'''
#        get NELM, the max electronic steps allowed in the run. 
        proc = subprocess.Popen(['grep','-i','NELM',dir+'/INCAR'],stdout=subprocess.PIPE)
        result =  proc.communicate()[0]
        NELM = int(result.split('=')[1].split()[0])
        return self.elConvergeCheck(dir,NELM) and finishCheck(dir)

    def getPlanarBE(self, atom, atomDir):
        '''Calculates the difference (per adatom) between the old energy and the new energy
    and the fully relaxed energy. Writes to file planarBE of the form "struct,  BE"in the analysis folder '''    
        topDir = os.getcwd()
        structlist = []
        toStart = []
        for item in os.listdir(atomDir):
            itempath = atomDir + '/' + item
            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
                structlist.append(item)
        sstring = ''; 
        for struct in structlist: 
            sstring += struct + ' '
        subprocess.call(['echo','Found structures '+ sstring])
        if not os.path.exists(atomDir + '/analysis'): os.mkdir(atomDir + '/analysis') 
        file = open(atomDir + '/analysis/planarBE_dz_{}'.format(self.deltaz),'w' )      
        for istruct,struct in enumerate(sorted(structlist)):
            if mod(istruct+1,100) == 0 or istruct+1 == len(structlist): subprocess.call(['echo','\tPlanarBE analyzed for {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
            structDir = atomDir + '/' + struct
            planarDir = atomDir + '/' + struct + '/planarBE'
            if not os.path.exists(structDir):
                subprocess.call(['echo','No vasp folder for {}, struct {}'.format(atom,struct)])
            elif not os.path.exists(planarDir):
                subprocess.call(['echo','No planar BE folder for {}, struct {}'.format(atom,struct)]) 
            elif finishCheck(planarDir) and  os.path.getsize(structDir + '/planarBE/POSCAR') > 0:                
                self.setAtomCounts(planarDir)
                BE = (self.getEnergy(structDir) - self.getEnergy(planarDir))/sum(self.atomCounts[1:])#/self.deltaz
                file.write('{} {:8.4f}\n'.format(struct,BE))
            else:
                subprocess.call(['echo','Structure ' + struct + ' has not finished the planar BE vasp calculation.'])
            os.chdir(topDir)
        file.close()

    def makePlanarJobFile(self,planarDir,name,walltime):
        jobFile = open(planarDir +'/job','w')
        jobFile.write("#!/bin/bash\n\n")
        jobFile.write("#SBATCH --time={}:00:00\n".format(walltime))
        jobFile.write("#SBATCH --ntasks=16\n")
        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
        jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
        jobFile.write("#SBATCH --mail-type=END\n")
        jobFile.write("#SBATCH --mail-type=FAIL\n")  
        jobFile.write("#SBATCH --job-name=PBE%s\n" % name)
        jobFile.write("\nmpiexec vasp533 > vasp.out\n")   
        jobFile.close() 

    def modifyINCAR(self,planarDir):
        '''Changes it so it has no ionic relaxation.  Ignore CHG and WAVECAR files'''
        lines = readfile(planarDir + '/INCAR')
        lines2 = deepcopy(lines)
        for i,line in enumerate(lines):
            if 'IBRION' in line:
                lines2[i] = 'IBRION=2\n'  #testing with a Cr top structure 2, find we must have IBRION 2 and NSW 1.  NOT -1 and 0 as expected for static run!!!
            elif 'NSW' in line:
                lines2[i] = 'NSW=1\n'  #testing with a Cr top structure 2, find we must have IBRION 2 and NSW 1.  NOT -1 and 0 as expected for static run!!!
#            elif 'LCHARG' in line:
#                lines2[i] = 'LCHARG=.FALSE.\n' 
#            elif 'LWAV' in line:
#                lines2[i] = 'LWAV=.FALSE.\n'
        writefile(lines2,planarDir + '/INCAR') 
                          
    def modifyPOSCAR(self,planarDir):
        '''all adatoms are displaced deltaz away from the plane.  Since this POSCAR
        came from a CONTAR, we are working in direct coordinates'''
        lines = readfile(planarDir + '/POSCAR')
        lines2 = deepcopy(lines)
        lattvecs = lines[2:5]
        lattvecs = [line.strip().split() for line in lattvecs]
        zLV = max([abs(float(lattvecs[i][2])) for i in range(3)]) #z repeat distance in ang 
        izLV = argmax([abs(float(lattvecs[i][2])) for i in range(3)])                
        self.setAtomCounts(planarDir)
        self.nC = self.atomCounts[0] 
        self.nH = self.atomCounts[1] 
        self.nM = self.atomCounts[2]
        self.nTot = sum(self.atomCounts)
        zLV2 = self.newLVz
#        zLV2 = zLV #keep it unchanged for now
        dz = self.deltaz/zLV2
        lines2[2+izLV] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(float(lattvecs[izLV][0]),float(lattvecs[izLV][1]),zLV2)
        for i in range(self.nTot):
            pos = [float(part) for part in lines[8 + i].strip().split()]
            if i < self.nC:
                if pos[2] > 0.5: #don't want Vasp's convention to shift atoms to all positive coordinates
                    pos[2] -= 1
                lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV)
            else:
                if pos[2] > 0.5:  
                    pos[2] -= 1
                lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV + sign(pos[2])*dz)
        writefile(lines2,planarDir + '/POSCAR')
                
    def planarBEruns(self, walltime):  
        '''Script for this process'''
        topDir = os.getcwd()
        for iatom, atom in enumerate(self.atoms):
            atomDir = topDir + '/' + atom       
            toStart = self.createFolders(atom,atomDir,walltime)
            self.startJobs(toStart,atomDir,'planarBE') 
            
    def planarBEanalysis(self):
        topDir = os.getcwd()
        for iatom, atom in enumerate(self.atoms):
            atomDir = topDir + '/' + atom       
            self.getPlanarBE(atom,atomDir)
            self.getPlanarBE(atom,atomDir)

    def planarBEplot(self):
        topDir = os.getcwd()
        for iatom, atom in enumerate(self.atoms):
            atomDir = topDir + '/' + atom       
            self.plotPlanarBE(atom,atomDir)
            self.plotPlanarBE(atom,atomDir)                
            self.collatePlotsGSS('HFE_planarBE','xx')
            
    def plotPlanarBE(self,atom,atomDir):
        '''Plot the vaspHFE vs energy, and make the marker color reflect the planar BE'''
        lines = readfile(atomDir + '/analysis/analysis.csv')
        planarLines = readfile(atomDir + '/analysis/planarBE_dz_{}'.format(self.deltaz))
        structs = []
        conc = []
        HFE = []
        planarBE = []
#        concPos = []; concNeg = []
#        HFEPos = []; HFENeg = []
#        planarBEPos = []; planarBENeg = [
        ikeep = -1
        for ip,pline in enumerate(planarLines):
            if mod(ip+1,100) == 0 or ip+1 ==len(planarLines): subprocess.call(['echo','\tPlanarBE plot data read for {} of {} structures in {}'.format(ip+1,len(planarLines),atom)])            
            pstruct = int(pline.split()[0])
            BE = float(pline.split()[1])
            
            for line in lines[1:]:
                struct = int(line.split(',')[0])
                if struct == pstruct:
                    if BE<0:
                        conc.append(float(line.split(',')[2]))
                        HFE.append(float(line.split(',')[3]))
                        planarBE.append(BE)
                        break
                    else:
                        subprocess.call(['echo','Struct {} for atom {} has positive planar BE {:7.4f}...not included in plot\n'.format(struct,atom,BE)])  
#        planarBE = log10(abs(array(planarBE)))

        #plot
        fig = plt.figure()
        plt.xlabel('{} concentration x'.format(atomDir.split('/')[-1]))
        plt.ylabel('Energy (eV)')
#        plt.title('FE colored by planar force (eV/$\AA$) at {:2.1f} $\AA$'.format(self.deltaz/2.0))        
        plt.title('FE colored by planar BE (eV), displaced {:2.1f} $\AA$'.format(self.deltaz))
#        plt.scatter(conc, HFE, c=info , cmap='autumn')
        plt.scatter(conc, HFE, c=planarBE , cmap='jet')
        plt.colorbar()
        plt.xlim(0, 1.0)
#        plt.ylim(plt.ylim()[0], min(max(HFE),2.0))
        plt.ylim(min(HFE)-0.1, min(max(HFE),2.0))
        plt.show()
        fig.savefig(atomDir + '/gss/HFE_planarBE') #xx because iteration number is not used
        plt.close(fig)   


#class PlanarBE:
#    '''Finds the binding energy of the adatoms to the carbon plane by 1) creating run folders 
#    with POSCAR altered so all adatoms are displaced deltaz away from the plane, 
#    2) running those jobs without ionic relaxation  3) calculating the difference (per adatom) in this energy
#    and the fully relaxed energy.'''
#    from comMethods import setAtomCounts,startJobs,getEnergy, elConvergeCheck, getElSteps,collatePlotsGSS
#    def __init__(self, atoms):
#        self.atoms = atoms
##        self.deltaz = 0.2  # perturbation   
#        self.deltaz = 5   
#        self.nC = 0 
#        self.nH = 0 
#        self.nM = 0
#        self.jobIds = []
# 
#    def createFolders(self,atom,atomDir,walltime,type):
#        '''Copies from previous run folders, but modifies the POSCARs with deltaz for adatoms, and changes INCAR so atoms don't relax'''
#        topDir = os.getcwd()
#        structlist = []
#        toStart = []
#        for item in os.listdir(atomDir):
#            itempath = atomDir + '/' + item
#            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
#                structlist.append(item)
#        sstring = ''; 
#        for struct in structlist: 
#            sstring += struct + ' '
#        subprocess.call(['echo','Found structures '+ sstring])
#        analysisDir = '{}/analysis'.format(atomDir)       
#        for istruct,struct in enumerate(structlist):
#            if mod(istruct+1,100) == 0 or istruct+1 ==len(structlist): subprocess.call(['echo','\tPlanarBE prep for {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
#            structDir = atomDir + '/' + struct
#            planarDir = atomDir + '/' + struct + '/planarBE_{}'.format(type)
#            if not os.path.exists(structDir):
#                subprocess.call(['echo','No vasp folder for {}, struct {}'.format(atom,struct)])
#            elif struct != '1' and finishCheck(structDir) and  convergeCheck(structDir, getNSW(structDir)) and energyDropCheck(structDir):
#                #remove this if statement when you want to start over on planarBE calcs. 
##                if not os.path.exists(planarDir) or not finishCheck(planarDir): #not doing electronic convergence check for this precision...good within 1-5%% of BE
#                    if os.path.exists(planarDir): subprocess.call(['rm','-r',planarDir]) #start over for now
#                    os.mkdir(planarDir)
#                    subprocess.call(['ln','-s','/fslhome/bch/bin/vasp533',planarDir+'/vasp533'])
#                    subprocess.call(['cp',structDir+'/CONTCAR',planarDir+ '/POSCAR'])
#                    subprocess.call(['cp',structDir+'/INCAR',planarDir])
##                    subprocess.call(['ln','-s',structDir+'/CHG',planarDir]) #Use only for perturbation (force)
##                    subprocess.call(['ln','-s',structDir+'/CHGCAR',planarDir]) #Use only for perturbation (force)
##                    subprocess.call(['ln','-s',structDir+'/WAVECAR',planarDir]) #Use only for perturbation (force)
#                    subprocess.call(['ln','-s',structDir+'/POTCAR',planarDir])
#                    subprocess.call(['ln','-s',structDir+'/KPOINTS',planarDir])
#                    self.makePlanarJobFile(planarDir,atom+struct,walltime)
#                    self.modifyINCAR(planarDir)
#                    self.modifyPOSCAR1Atom(planarDir,self.getdMCatomIndex(struct,type,analysisDir))
#                    toStart.append(struct)
#            elif struct == '1':
#                'do nothing'
#            else:
#                print 
#                subprocess.call(['echo','Structure ' + struct + ' has not finished, converged or has an energy rise.'])
#            os.chdir(topDir)
#        return toStart            
#
#    def electronicConvergeFinish(self,dir): 
#        '''Test requires electronic convergence AND vasp finishing'''
##        get NELM, the max electronic steps allowed in the run. 
#        proc = subprocess.Popen(['grep','-i','NELM',dir+'/INCAR'],stdout=subprocess.PIPE)
#        result =  proc.communicate()[0]
#        NELM = int(result.split('=')[1].split()[0])
#        return self.elConvergeCheck(dir,NELM) and finishCheck(dir)
#
#    def getdMCatomIndex(self,struct,type,analysisDir):
#        '''returns the metal atom index (among all the atoms in a structure) that has 
#        the max or min dMC.'''
#        lines = readfile('{}/dMCmaxMin'.format(analysisDir))
#        for i,line in enumerate(lines):
#            data = line.split()
#            if data[0] == struct:
#                if type == 'min':
#                    return int(data[1])
#                elif type == 'max':
#                    return int(data[3])
#                break
#
#    def getPlanarBE(self, atom, atomDir, type):
#        '''Calculates the difference (per adatom) between the old energy and the new energy
#    and the fully relaxed energy. Writes to file planarBE of the form "struct,  BE"in the analysis folder '''    
#        topDir = os.getcwd()
#        structlist = []
#        toStart = []
#        for item in os.listdir(atomDir):
#            itempath = atomDir + '/' + item
#            if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
#                structlist.append(item)
#        sstring = ''; 
#        for struct in structlist: 
#            sstring += struct + ' '
#        subprocess.call(['echo','Found structures '+ sstring]) 
#        file = open(atomDir + '/analysis/planarBE_dz_{}_{}'.format(self.deltaz,type),'w' )      
#        for istruct,struct in enumerate(sorted(structlist)):
#            if mod(istruct+1,100) == 0 or istruct+1 == len(structlist): subprocess.call(['echo','\tPlanarBE analyzed for {} of {} structures in {}'.format(istruct+1,len(structlist),atom)])
#            structDir = atomDir + '/' + struct
#            planarDir = atomDir + '/' + struct + '/planarBE_{}'.format(type)
#            if not os.path.exists(structDir):
#                subprocess.call(['echo','No vasp folder for {}, struct {}'.format(atom,struct)])
#            elif not os.path.exists(planarDir):
#                subprocess.call(['echo','No planar BE folder for {}, struct {}'.format(atom,struct)]) 
#            elif finishCheck(planarDir) and  os.path.getsize(structDir + '/planarBE/POSCAR') > 0:                
#                self.setAtomCounts(planarDir)
#                BE = (self.getEnergy(structDir) - self.getEnergy(planarDir))/sum(self.atomCounts[1:])#/self.deltaz
#                file.write('{} {:8.4f}\n'.format(struct,BE))
#            else:
#                subprocess.call(['echo','Structure ' + struct + ' has not finished the planar BE vasp calculation.'])
#            os.chdir(topDir)
#        file.close()
#
#    def makePlanarJobFile(self,planarDir,name,walltime):
#        jobFile = open(planarDir +'/job','w')
#        jobFile.write("#!/bin/bash\n\n")
#        jobFile.write("#SBATCH --time={}:00:00\n".format(walltime))
#        jobFile.write("#SBATCH --ntasks=16\n")
#        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")
#        jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
#        jobFile.write("#SBATCH --mail-type=END\n")
#        jobFile.write("#SBATCH --mail-type=FAIL\n")  
#        jobFile.write("#SBATCH --job-name=PBE%s\n" % name)
#        jobFile.write("\nmpiexec vasp533 > vasp.out\n")   
#        jobFile.close() 
#
#    def modifyINCAR(self,planarDir):
#        '''Changes it so it has no ionic relaxation.  Ignore CHG and WAVECAR files'''
#        lines = readfile(planarDir + '/INCAR')
#        lines2 = deepcopy(lines)
#        for i,line in enumerate(lines):
#            if 'IBRION' in line:
#                lines2[i] = 'IBRION=2\n'  #testing with a Cr top structure 2, find we must have IBRION 2 and NSW 1.  NOT -1 and 0 as expected for static run!!!
#            elif 'NSW' in line:
#                lines2[i] = 'NSW=1\n'  #testing with a Cr top structure 2, find we must have IBRION 2 and NSW 1.  NOT -1 and 0 as expected for static run!!!
#            elif 'LCHARG' in line:
#                lines2[i] = 'LCHARG=.FALSE.\n' 
#            elif 'LWAV' in line:
#                lines2[i] = 'LWAV=.FALSE.\n'
#        writefile(lines2,planarDir + '/INCAR') 
#                          
##    def modifyPOSCAR(self,planarDir):
##        '''all adatoms are displaced deltaz away from the plane.  Since this POSCAR
##        came from a CONTAR, we are working in direct coordinates'''
##        lines = readfile(planarDir + '/POSCAR')
##        lines2 = deepcopy(lines)
##        lattvecs = lines[2:5]
##        lattvecs = [line.strip().split() for line in lattvecs]
##        zLV = max([abs(float(lattvecs[i][2])) for i in range(3)]) #z repeat distance in ang 
##        izLV = argmax([abs(float(lattvecs[i][2])) for i in range(3)])                
##        self.setAtomCounts(planarDir)
##        self.nC = self.atomCounts[0] 
##        self.nH = self.atomCounts[1] 
##        self.nM = self.atomCounts[2]
##        self.nTot = sum(self.atomCounts)
###        zLV2 = 25
##        zLV2 = zLV #keep it unchanged for now
##        dz = self.deltaz/zLV2
###        lines2[2+izLV] = '{}   {}   {}\n'.format(lattvecs[izLV][0],lattvecs[izLV][1],zLV2)
##        for i in range(self.nTot):
##            pos = [float(part) for part in lines[8 + i].strip().split()]
##            if i < self.nC:
##                if pos[2] > 0.5: #don't want Vasp's convention to shift atoms to all positive coordinates
##                    pos[2] -= 1
##                lines2[8 + i] = '{}     {}     {}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV)
##            else:
##                if pos[2] > 0.5:  
##                    pos[2] -= 1
##                lines2[8 + i] = '{}     {}     {}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV + sign(pos[2])*dz)
##        writefile(lines2,planarDir + '/POSCAR')
#        
#        
#    def modifyPOSCAR1Atom(self,planarDir,index):
#            '''Moves one metal atom with either maximum or minimum (given by type)metal-carbon distance, by deltaz.  Since this POSCAR
#            came from a CONTAR, we are working in direct coordinates'''
#            lines = readfile(planarDir + '/POSCAR')
#            lines2 = deepcopy(lines)
#            lattvecs = lines[2:5]
#            lattvecs = [line.strip().split() for line in lattvecs]
#            zLV = max([abs(float(lattvecs[i][2])) for i in range(3)]) #z repeat distance in ang 
#            izLV = argmax([abs(float(lattvecs[i][2])) for i in range(3)])                
#            self.setAtomCounts(planarDir)
#            self.nC = self.atomCounts[0] 
#            self.nH = self.atomCounts[1] 
#            self.nM = self.atomCounts[2]
#            self.nTot = sum(self.atomCounts)
#    #        zLV2 = 25
#            zLV2 = zLV #keep it unchanged for now
#            dz = self.deltaz/zLV2
#    #        lines2[2+izLV] = '{}   {}   {}\n'.format(lattvecs[izLV][0],lattvecs[izLV][1],zLV2)
#            for i in range(self.nTot):
#                pos = [float(part) for part in lines[8 + i].strip().split()]
#                if i < self.nC: #just rewrite, no shift for these atoms
#                    if pos[2] > 0.5: #don't want Vasp's convention to shift atoms to all positive coordinates
#                        pos[2] -= 1
#                    lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV)
#                elif i-self.nH-self.nC == index: #this is the metal atom we want to move
#                    if pos[2] > 0.5:  
#                        pos[2] -= 1
#                    lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV + sign(pos[2])*dz)
#                else: #just rewrite, no shift for these atoms
#                    if pos[2] > 0.5:  
#                        pos[2] -= 1 
#                    lines2[8 + i] = '{:16.12f} {:16.12f} {:16.12f}\n'.format(pos[0],pos[1],pos[2]*zLV2/zLV)
#            writefile(lines2,planarDir + '/POSCAR')
#        
#    def planarBEruns(self, walltime):  
#        '''Script for this process'''
#        topDir = os.getcwd()
#        for iatom, atom in enumerate(self.atoms):
#            atomDir = topDir + '/' + atom       
#            toStart = self.createFolders(atom,atomDir,walltime,'max')
#            self.startJobs(toStart,atomDir,'planarBE_max') 
#            toStart = self.createFolders(atom,atomDir,walltime,'min')
#            self.startJobs(toStart,atomDir,'planarBE_min') 
#            
#    def planarBEanalysis(self):
#        topDir = os.getcwd()
#        for iatom, atom in enumerate(self.atoms):
#            atomDir = topDir + '/' + atom       
#            self.getPlanarBE(atom,atomDir,'max')
#            self.getPlanarBE(atom,atomDir,'min')
#
#    def planarBEplot(self):
#        topDir = os.getcwd()
#        for iatom, atom in enumerate(self.atoms):
#            atomDir = topDir + '/' + atom       
#            self.plotPlanarBE(atom,atomDir,'max')
#            self.plotPlanarBE(atom,atomDir,'min')                
#            self.collatePlotsGSS('HFE_planarBE_max','xx')
#            self.collatePlotsGSS('HFE_planarBE_min','xx')
#            
#    def plotPlanarBE(self,atom,atomDir,type):
#        '''Plot the vaspHFE vs energy, and make the marker color reflect the planar BE'''
#        lines = readfile(atomDir + '/analysis/analysis.csv')
#        planarLines = readfile(atomDir + '/analysis/planarBE_dz_{}_{}'.format(self.deltaz,type))
#        structs = []
#        conc = []
#        HFE = []
#        planarBE = []
##        concPos = []; concNeg = []
##        HFEPos = []; HFENeg = []
##        planarBEPos = []; planarBENeg = [
#        ikeep = -1
#        for ip,pline in enumerate(planarLines):
#            if mod(ip+1,100) == 0 or ip+1 ==len(planarLines): subprocess.call(['echo','\tPlanarBE plot data read for {} of {} structures in {}'.format(ip+1,len(planarLines),atom)])            
#            pstruct = int(pline.split()[0])
#            BE = float(pline.split()[1])
#            
#            for line in lines[1:]:
#                struct = int(line.split(',')[0])
#                if struct == pstruct:
#                    if BE<0:
#                        conc.append(float(line.split(',')[2]))
#                        HFE.append(float(line.split(',')[3]))
#                        planarBE.append(BE)
#                        break
#                    else:
#                        subprocess.call(['echo','Struct {} for atom {} has positive planar BE {:7.4f}...not included in plot\n'.format(struct,atom,BE)])  
##        planarBE = log10(abs(array(planarBE)))
#
#        #plot
#        fig = plt.figure()
#        plt.xlabel('{} concentration x'.format(atomDir.split('/')[-1]))
#        plt.ylabel('Energy (eV)')
##        plt.title('FE colored by planar force (eV/$\AA$) at {:2.1f} $\AA$'.format(self.deltaz/2.0))        
#        plt.title('FE colored by planar BE (eV), displaced {:2.1f} $\AA$'.format(self.deltaz))
##        plt.scatter(conc, HFE, c=info , cmap='autumn')
#        plt.scatter(conc, HFE, c=planarBE , cmap='jet')
#        plt.colorbar()
#        plt.xlim(0, 1.0)
##        plt.ylim(plt.ylim()[0], min(max(HFE),2.0))
#        plt.ylim(min(HFE)-0.1, min(max(HFE),2.0))
#        plt.show()
#        fig.savefig(atomDir + '/gss/HFE_planarBE_{}'.format(type)) #xx because iteration number is not used
#        plt.close(fig)   

class Analysis:
    '''Includes movement of atoms and magnetic moments'''
    from comMethods import setAtomCounts,collatePlotsGSS
    def __init__(self, atoms,pureMetal,iteration,plotOnly):
        """ CONSTRUCTOR """
    
        self.atoms = atoms
        self.plotOnly = plotOnly  
        self.disStructList = []        
        self.inPlaneMove = []
        self.normalMove = []
        self.totalMove = []
        self.displace = None
        self.atomCounts = [] 
        self.nC = 0
        self.nH = 0
        self.nM = 0
        self.ntot = 0
        self.origCPos = []
        self.origMPos = []
        self.origHPos = []
        self.nomapOrigPos = [] #need to easily find NNs
        self.origPos = [] #all atoms' positions of current structure
        self.relaxPos = [] #all atoms' positions of current structure
        self.relaxCPos = []
        self.relaxMPos = []
        self.relaxHPos = []
        self.CCexpansion = 0.0
        self.origLattVecs = None 
        self.relaxLattVecs = None
        self.indexLVz = None
        self.hfe = []
        self.output = []
        self.pureMetal = pureMetal
        self.iteration = iteration
        self.struct = ''
        self.structDir = ''
        self.buckleAvg = 0.0
        self.conc = 0.0
        self.moveTotal = 0.0
        self.volFactor = 0
        self.moveInPlane = 0.0
        self.moveMInPlane = 0.0
        self.moveHInPlane = 0.0
        self.moveNormal = 0.0       
        self.closestMC = []
        self.magneticMom = 0.0

    def atomType(self,i):
        if 0 <= i < self.nC:
            atomtype = 'C'
        elif self.nC <= i < self.nC + self.nH:
            atomtype = 'H'
        elif self.nC + self.nH <= i < self.nTot:
            atomtype = 'M'
        return atomtype
        
    def correctRelaxPos(self):
        '''Contcar moves any position outside of cell into cell.  
         We want to change that, because we can't easily find true small movements from initial positions
         that are close to the cell edge'''
        iplanar = [0,1,2]
        del iplanar[self.indexLVz]
        for i,relaxpos in enumerate(self.relaxPos):
            origpos =  self.origPos[i,:]           
            dmin = 100
            for j in range(-1,2):
                for k in range(-1,2): 
                    shift = j*self.relaxLattVecs[iplanar[0],:] + k*self.relaxLattVecs[iplanar[1],:]              
                    dtest = norm(relaxpos + shift - origpos)
                    if dtest<dmin:
                        dmin = dtest
                        bestshift = shift
            relaxpos = relaxpos + bestshift
            self.relaxPos[i,:] = relaxpos
        self.relaxCPos = self.relaxPos[:self.nC,:]
        self.relaxHPos = self.relaxPos[self.nC:self.nC+self.nH,:]
        self.relaxMPos = self.relaxPos[-self.nM:,:]
               
    def correctz(self,directvec):
        '''Translates direct vectors with z positions farther than 0.5 to the cell edge to closer than that,
        because we like to keep z centered around zero'''
        eps = 0.01
        if directvec[self.indexLVz] > 0.5 :
            directvec[self.indexLVz] -= 1
        elif directvec[self.indexLVz] < -0.5:
            directvec[self.indexLVz] += 1    
        return directvec
    
    def getOrigPos(self):
        '''These are really the original positions mapped into the relaxed unit cell.  Will be in Cartesian coordinates. Go back to psuedoPOSCAR because POSCAR
        might have been overwritten by CONTCAR.  '''
        
        subprocess.call(['../needed_files/makestr.x','../enum/struct_enum.out',str(self.struct)])
        vfile = 'vasp.' + '0'*(6-len(str(self.struct))) + str(self.struct)
        toPoscar = StructsToPoscar.structsToPoscar([],[])
        toPoscar.convertOne(vfile) 
        subprocess.call(['rm', vfile])         
        poscarLines = readfile('POSCAR')
        subprocess.call(['rm', 'POSCAR'])       
        lattvecs = poscarLines[2:5]
        lattvecs = [line.strip().split() for line in lattvecs]
        self.origLattVecs = zeros((3,3),dtype = float)
        for i,vec in enumerate(lattvecs):
            newvec = [float(comp) for comp in vec]
            self.origLattVecs[i,:] = newvec
        positionLines = poscarLines[7:7 + self.nTot] # from vfile, these are already in direct lattice coordinates
        # Convert to Direct coordinates in terms of the !! NEW !! lattice vectors and then back to
        # Cartesian coordinates.
        positions = []
        for line in positionLines:
            position = line.strip().split()
            positions.append([float(comp) for comp in position])
        self.nomapOrigPos = array(positions) #need these to find NNpairs
        origPos = []
        for position in positions:
            #convert to direct
            directPos = dot(inv(transpose(self.origLattVecs)), transpose(position)) # Change to direct coordinates in the old cell
            origPos.append(dot(transpose(self.relaxLattVecs), transpose(directPos)))
        self.origPos = array(origPos)
        self.origCPos = self.origPos[:self.nC]
        self.origHPos = self.origPos[self.nC:self.nC+self.nH]
        self.origMPos = self.origPos[:-self.nM]
        #find C-C expansion vs starting lattice (%)
        self.CCexpansion = 100*(sqrt(abs(det(self.relaxLattVecs)/det(self.origLattVecs)\
                        *self.origLattVecs[self.indexLVz,2]/self.relaxLattVecs[self.indexLVz,2]))\
                        -1.0)
#        if self.struct == '3077':
#            print '3077',det(self.relaxLattVecs)/det(self.origLattVecs),self.origLattVecs[self.indexLVz,2],self.relaxLattVecs[self.indexLVz,2]
#            print 'origLattVecs', self.origLattVecs
#            print det(self.origLattVecs)   
    def getRelaxPos(self):
        '''Will be in Cartesian coordinates. Also gets atomCounts'''
        self.setAtomCounts(self.structDir)
        self.nC = self.atomCounts[0] 
        self.nH = self.atomCounts[1] 
        self.nM = self.atomCounts[2]
        self.nTot = sum(self.atomCounts) 
        self.volFactor = self.nC/2
        self.relaxLattVecs = zeros((3,3),dtype = float)
        contcarLines = readfile(self.structDir + '/CONTCAR')
        lattvecs = contcarLines[2:5]
        lattvecs = [line.strip().split() for line in lattvecs]    
        for i,vec in enumerate(lattvecs):
            newvec = [float(comp) for comp in vec]
            self.relaxLattVecs[i,:] = newvec
            cart = dot(transpose(self.relaxLattVecs), transpose(newvec))
#            if abs(newvec[2]) > max(abs(cart[1]),abs(cart[1])):
#                self.indexLVz = i #which LVector is in the z dir
        self.indexLVz = argmax([abs(float(lattvecs[i][2])) for i in range(3)])   
        positionLines = contcarLines[8:8 + self.nTot]     
        relaxPos = []
        for line in positionLines:
            newPosition = line.strip().split()
            newPosition = [float(pos) for pos in newPosition]
            newPosition = self.correctz(newPosition)
            # Already direct from Contcar.  Convert to Cartesian coordinates
            relaxPos.append(dot(transpose(self.relaxLattVecs), transpose(newPosition)))
        self.relaxPos = array(relaxPos)

    def getMove(self):
        '''Movements of atoms'''         
        self.displace = self.relaxPos - self.origPos
        self.moveInPlane = zeros(self.nTot,dtype = float)
        self.moveNormal = zeros(self.nTot,dtype = float)
        self.moveTotal = zeros(self.nTot,dtype = float)
        self.moveMInPlane = zeros(self.nTot,dtype = float)
        self.moveHInPlane = zeros(self.nTot,dtype = float)
        
        for i in range(len(self.displace)):
            self.moveInPlane[i] = norm(self.displace[i,:2]) #only x,y parts
            self.moveNormal[i] = abs(self.displace[i,2])
            self.moveTotal[i] = norm(self.displace[i,:])
            if self.atomType(i) =='M': 
                self.moveMInPlane[i] = self.moveInPlane[i]
            elif self.atomType(i) =='H':
                self.moveHInPlane[i] = self.moveInPlane[i] 
            
    def getMCBondingInfo(self):
        '''Find distances to nearest C atom''' 
        self.closestMC = []     
        for Mvec in self.relaxMPos:
            trialLengths = []
            for Cvec in self.relaxCPos:
                trialLengths.append(norm(Mvec-Cvec))
            self.closestMC.append(min(trialLengths)) #closest for this M atom 
        self.closestMC = array(self.closestMC)
           
    def getBucklingInfo(self): 
        NNpairs = self.getNNPairs(self.nomapOrigPos[:self.nC])   
        buckleDistances = []
        for pair in NNpairs:  
            buckleDistances.append(abs(self.relaxCPos[pair[0],2]-self.relaxCPos[pair[1],2])) #z components 
        self.buckleAvg =  sum(buckleDistances) / len(buckleDistances)
    
    def getNNPairs(self, Cpos):
        eps = 1e-4
        pairs = []
        for i in range(len(Cpos)):
            for j in range(len(Cpos)):
                if i != j:
                    d = norm(Cpos[i,:2] - Cpos[j,:2]) #take off the z components...looking for inplane distances
                    if abs(d - 1.42085899)<= eps and not self.NNPairExists([i,j], pairs):
                        pairs.append([i,j])       
        return pairs         
    
    def NNPairExists(self, pair, pairlist):
        if len(pairlist) == 0:
            return False
        else:
            for x in pairlist:
                if x[0] == pair[0] and x[1] == pair[1]:
                    return True
                elif x[0] == pair[1] and x[1] == pair[0]:
                    return True
            return False
    
    def getStructList(self,atomDir):
        structList = []
        hlines = readfile(atomDir + '/gss/vaspHFE_{}.out'.format(self.iteration))
        for line in hlines:
            structList.append(line.split()[0])
        self.hfe = zeros((len(structList)),dtype = [('struct', int32),('HFE', float)])
        for i,line in enumerate(hlines):
            self.hfe[i]['struct'] = line.split()[0]
            self.hfe[i]['HFE'] = line.split()[2]
        return structList
    
    def getAnalysis(self):
        topDir = os.getcwd()
        header =  "Structure, vol, conc, HFE, magnetic moment, CC-exp, moved_rms, max moved_para, rms M moved_para, max M moved_para, \
                rms H moved_para, max H moved_para, rms moved_perp, max moved_perp, min d_MC, rms d_MC, max d_MC, ave buckle"       
        self.output = [item.strip() for item in header.split(',')]
        if not self.plotOnly:
            for atom in self.atoms:
                atomDir = topDir + '/' + atom
                analysisDir = atomDir + '/analysis'
                if not os.path.exists(analysisDir): os.mkdir(analysisDir)
                csvfile = open(analysisDir +'/analysis.csv','w')
                csvfile.write(header + "\n")
                MCmaxMinfile = open(analysisDir +'/dMCmaxMin','w')
                structList = self.getStructList(atomDir)
                subprocess.call(['echo','************************************'])
                subprocess.call(['echo','    Analysis, ' + atom])
                subprocess.call(['echo','************************************'])
                os.chdir(atomDir)
                for istruct, structure in enumerate(structList):
                        if mod(istruct+1,100) == 0 or istruct+1 ==len(structList): subprocess.call(['echo','\tAnalyzed {} of {} structures for {}'.format(istruct+1,len(structList),atom)])
        #                print structure
                        self.struct = structure
                        self.structDir = atomDir + '/' + structure
                        self.getRelaxPos()
                        self.getOrigPos()
                        self.correctRelaxPos()   
                        self.getMove()
                        self.getBucklingInfo()
                        self.getMagnetization()
                        self.writeInfoToFile()
                        self.writeToCSV(csvfile)
                        self.writeMCminMax(MCmaxMinfile)
                csvfile.close()
                MCmaxMinfile.close()
        plots = ['max M moved_para','max H moved_para','max d_MC','rms d_MC','CC-exp',\
                 'rms M moved_para','rms H moved_para','magnetic moment'] 
        for atom in self.atoms:
            subprocess.call(['echo','****************************'])
            subprocess.call(['echo','    Analysis plots, ' + atom])
            subprocess.call(['echo','****************************'])
            atomDir = topDir + '/' + atom 
            for plot in plots:      
                self.plotHFEcoloredInfo(atomDir,plot)
        os.chdir(topDir)
        for plot in plots:      
            self.collatePlotsGSS('HFE_'+ plot.replace(' ','_',),self.iteration)
            
    def getMagnetization(self):       
        '''Magnetic moments of each structure from OUTCAR'''
        lastdir = os.getcwd()
        os.chdir(lastdir + '/' + self.struct)
        goodlines = []
        lines = readfile(lastdir + '/' + self.struct + '/OUTCAR')
        for line in lines:
            if 'number of electron' in line:
                goodlines.append(line)   
        self.magneticMom = float(goodlines[-1].split(' ')[-1]) 
        if self.nM > 0: 
            self.magneticMom = abs(self.magneticMom/float(self.nM))
        os.chdir(lastdir)  

    def plotHFEcoloredInfo(self,atomDir,infoType):
        '''Plot the vaspHFE vs energy, and make the marker color reflect the parameter infoType, 
        which is a string that is part of the analysis.csv header, the substrings of which are
        listed in self.output'''
        lines = readfile(atomDir + '/analysis/analysis.csv')
        structs = []
        conc = []
        HFE = []
        info = []
        iinfo = self.output.index(infoType)
#        print iinfo
        for line in lines[1:]:
#            print infoType,line
            structure = int(line.split(',')[0])
            structs.append(structure)
            conc.append(float(line.split(',')[2]))
            HFE.append(float(line.split(',')[3]))
            info.append(float(line.split(',')[iinfo]))
        if 1 in structs and infoType in ['min d_MC','max d_MC','rms d_MC']: #to keep pure H structure with zero M from setting the scale on the colorbar    
            infotemp = deepcopy(info)
            infotemp.pop(structs.index(1))
            info[structs.index(1)] = min(infotemp)
#        if infoType == 'CC-exp':
#                    print 'ccexp',min(info),argmin(info),structs[argmin(info)]
        #plot
        fig = plt.figure()
        plt.xlabel('{} concentration x'.format(atomDir.split('/')[-1]))
        plt.ylabel('Energy (eV)')
        plt.title('Formation energy colored by {}'.format(infoType))
#        plt.scatter(conc, HFE, c=info , cmap='autumn')
        plt.scatter(conc, HFE, c=info , cmap='jet')
#        print 'plot',infoType
#        print info
        plt.colorbar()
        plt.xlim(0, 1.0)
#        plt.ylim(plt.ylim()[0], min(max(HFE),2.0))
        plt.ylim(min(HFE)-0.1, min(max(HFE),2.0))
        plt.show()
        fig.savefig(atomDir + '/gss/HFE_{}_{}'.format(infoType.replace(' ','_',),self.iteration))
        plt.close(fig)
        self.meanMoveMInPlane = 0.0
        self.meanMoveHInPlane = 0.0
        self.meanMoveNormal = 0.0  

    def writeMCminMax(self,MCmaxMinfile):
        '''For each structure writes the indices of the metal atoms with minimum and maximum bond lengths to C atom.  
        Used by planarBE routines.  Format of a line: structure, index of min dMC, min dMC, index of max dMC, max dMC'''
        MCmaxMinfile.write('{:10d} {:6d} {:6.3f} {:6d} {:6.3f}\n'.format(int(self.struct),argmin(self.closestMC),min(self.closestMC)\
                                                                       ,argmax(self.closestMC),max(self.closestMC)))

    def writeToCSV(self,csvfile):
        # "Structure, vol, conc, moved_rms, Max moved_para, Max moved_perp, Min d_MC, Max d_MC, Ave Buckle dz\n")
#        HFEindex = where(self.hfe['struct']==self.struct)
#        print 'testing'
#        print self.closestMC
#        print 'square'
#        print square(self.closestMC)
#        print 'mean'
#        print mean(self.closestMC)
        HFEindex = self.hfe['struct'].tolist().index(int(self.struct))
        csvfile.write("%s, %d, %7.3f, %7.3f, %7.3f,%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f\n" %
                      (self.struct, self.volFactor , self.conc , self.hfe[HFEindex]['HFE'], self.magneticMom, self.CCexpansion, rms(self.moveTotal), max(self.moveInPlane),  \
                       rms(self.moveMInPlane), max(self.moveMInPlane), rms(self.moveHInPlane), max(self.moveHInPlane), \
                       rms(self.moveNormal), max(self.moveNormal), min(self.closestMC), rms(self.closestMC),max(self.closestMC), self.buckleAvg))

    def writeInfoToFile(self):
        outfile = open(self.structDir + '/analysis.dat','w')       
        outfile.write("*************************************************\n")
        outfile.write("      MOVEMENT and MAGNETIC INFO FOR STRUCTURE " + self.struct + "\n")
        outfile.write("*************************************************\n\n")   
        outfile.write("Volume Factor: " + str(self.volFactor) + "\t")       
        self.conc = float(float(self.nM) / float(self.nH + self.nM))        
        outfile.write("Metal concentration: " + str(self.conc) + "\n\n")        
        outfile.write("RMS distance moved in relaxation: " + str(rms(self.moveTotal)) + "\n\n")        
        outfile.write("Maximum distance moved in the plane: " + str(max(self.moveInPlane)) + "\n\n")        
        outfile.write("Maximum distance moved perpendicular to the plane: " + str(max(self.moveNormal)) + "\n\n")        
        if self.struct == '1':
            self.closestMC = [0.0]
            outfile.write("Minimum M-C bond length: NONE \n\n")
            outfile.write("Maximum M-C bond length: NONE \n\n")
        else:
            self.getMCBondingInfo()
            outfile.write("Minimum M-C bond length: " + str(min(self.closestMC)) + "\n\n")
            outfile.write("Maximum M-C bond length: " + str(max(self.closestMC)) + "\n\n")
        outfile.write("Average buckling distance: " + str(self.buckleAvg) + "\n")
        outfile.write("\nOriginal Pos : \t\t\t Relaxed Position: \t\t Movement\t\t delta-x \t delta-y \t delta-z\n")
        outfile.write("(mapped to new cell)\n")
        for i in range(len(self.origPos)):
            o = self.origPos[i,:]
            r = self.relaxPos[i,:]
            outfile.write("%s  %7.3f %7.3f %7.3f \t %7.3f %7.3f %7.3f \t %7.3f \t\t %7.3f \t %7.3f \t %7.3f\n" \
                          % (self.atomType(i), o[0], o[1], o[2], r[0], r[1], r[2],\
                              self.moveTotal[i], self.displace[i,0], self.displace[i,1], self.displace[i,2]))        
        outfile.write("Magnetic moment: {} \n".format(self.magneticMom))
        outfile.close()
        
        
    
