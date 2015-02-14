'''A hack to create and run jobs for structure folders that replace H atoms with a vacancy.  
Assumes that the POSCAR is the original starting one (hasn't been overwritten by CONTCAR) '''
import os, subprocess,sys,re
from comMethods import *
from numpy import rint

def changePOSCAR(struct,poscarDir,):
    """ Removes all reference to H atoms, including their positions 
        the corresponding members. 
        Also fixes "new" POSCAR/CONTCAR format (comes from CONTCAR) back to old for UNCLE use (removes the 6th line if it's text """
    if struct not in ['1','2']:
        poscarLines = readfile(poscarDir + '/POSCAR')
        counts = poscarLines[5].strip().split()   
        if len(counts) == 3:
            nC = int(counts[0])
            nH = int(counts[1])
            nM = int(counts[2])
    #    elif len(counts) == 2: #just leave these alone...these are the pure cases: will fail, fix manually
    #        atomCounts.append(int(counts[0]))
    #        if poscarLines[0].split()[1] == 'H':
    #            atomCounts.append(int(counts[1]))
    #            atomCounts.append(0)
    #        elif poscarLines[0].split()[1] == 'M':
    #            atomCounts.append(0)
    #            atomCounts.append(int(counts[1]))
#        natoms1 = nC + nH + nM
#        natoms2 = nC + nM
        poscarLines[5] = '{} {}\n'.format(nC,nM)
        del poscarLines[7+nC:7+nC+nH] #remove H atoms
        writefile(poscarLines,poscarDir + '/POSCAR') #:7+natoms is because CONTCAR includes velocity lines that uncle doesn't want. The factor of 2 is because carbon atoms are not included in natoms      
    elif struct == '1': # it's pure graphene'
        poscarLines =['PURE graphene str #: 1\n', '1.0\n','  2.13128850  -1.23050000   0.00000000\n' ,\
                '  2.13128850   1.23050000   0.00000000\n','0.00000000   0.00000000  15.00000000\n',\
                '2\n','Cartesian\n','  1.42085899   0.00000000   0.22856000\n',\
                '  2.84171799   0.00000000  -0.22856000\n']
        writefile(poscarLines,poscarDir + '/POSCAR') 
    #structure 2 (metal pure is OK)

#======================================= Script =====================================
#======================================= Script =====================================
#======================================= Script =====================================
#maindir = '/fslhome/bch/cluster_expansion/graphene/'
maindir = '/fslhome/bch/cluster_expansion/graphene/hollowTiH.v8/'
finaldir = '/fslhome/bch/cluster_expansion/graphene/hollowTivac.v8/'
#finaldir = '/fslhome/bch/cluster_expansion/graphene/vac.top.tm_row1.v15/'
 
##    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1.continue'
#    maindir = '/fslhome/bch/cluster_expansion/graphene/tm_row1'
#maindir = os.getcwd()

#os.chdir(maindir)
#if not os.path.exists(finaldir):
#    subprocess.call(['mkdir',finaldir])


atomlist = []
for item in os.listdir(maindir):
    itempath = maindir + '/' + item
    if os.path.isdir(itempath) and item[0].isupper(): #look only at dirs whose names are capitalized. These are the atoms
        atomlist.append(item)
for atom in atomlist:
    print
    print atom
    atomdir = maindir + '/' + atom
    finalatomdir = finaldir + '/' + atom
    subprocess.call(['mkdir',finalatomdir])    
    structlist = []
    for item in os.listdir(atomdir):
        itempath = atomdir + '/' + item
        if os.path.isdir(itempath) and item[0].isdigit(): #look only at dirs whose names are numbers
            structlist.append(item)        
    for struct in structlist:
        if int(struct) <= 191:
            print struct
            structdir = atomdir + '/' + struct
            finalstructdir = finalatomdir + '/' + struct
            subprocess.call(['mkdir',finalstructdir])
            subprocess.call(['cp',structdir+'/POSCAR', finalstructdir])
            changePOSCAR(struct,finalstructdir) 
            subprocess.call(['cp',atomdir + '/' + 'C' + atom + '_POTCAR',finalstructdir + '/POTCAR'])
            subprocess.call(['cp',structdir+'/KPOINTS', finalstructdir])
            subprocess.call(['cp',structdir+'/INCAR', finalstructdir])
            subprocess.call(['cp',structdir+'/job', finalstructdir])
            os.chdir(finalstructdir)
            proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
            jobid = proc.communicate()[0].split()[3]
            subprocess.call(['echo','Submitted job ' + jobid])

print "Done submitting jobs"

