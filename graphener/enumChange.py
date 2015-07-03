import subprocess, time
from numpy import *

def changeEnumFile():
        enumFile = 'struct_enum.out_MAKE'
        subprocess.call(['mv',enumFile, enumFile + '_OLD'])

        oldfile = open(enumFile + '_OLD','r')
        oldlines = [line for line in oldfile]
        oldfile.close()

        enumFile = enumFile + '_NEW'
        
        newfile = open(enumFile, 'w')
        startline = 0
        for i in xrange(len(oldlines)):
            if i == 1:
                newfile.write('bulk\n')
            else:
                newfile.write(oldlines[i])
            if oldlines[i].strip().split()[0] == 'start':
                startline = i + 1
                break

        vec1 = oldlines[2].strip().split()[:3]
        vec2 = oldlines[3].strip().split()[:3]
        vec3 = oldlines[4].strip().split()[:3]
       
        vec1comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec1] #to handle ******
        vec2comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec2]
        vec3comps = [float(comp) if comp[0]!='*' else 1000 for comp in vec3]
    
        #Parent lattice vectors
        pLV = array([[vec1comps[0],vec2comps[0],vec3comps[0]], [vec1comps[1],vec2comps[1],vec3comps[1]], [vec1comps[2],vec2comps[2],vec3comps[2]]])
        nD = int(oldlines[5].strip().split()[0]) #Number of sites per cell
        pBas = array([[float(comp) for comp in line.strip().split()[:3]] for line in oldlines[6:6 + nD]]) #Base vectors for sites in the primitive cell

        structNum = 1
        enumVcNum = 1
        newHNF = True
        for i in xrange(startline,len(oldlines)):
            if i != startline:
                if oldlines[i-1].strip().split()[8:26] != oldlines[i].strip().split()[8:26]:
                    newHNF = True
            if newHNF == True:
                S = array([int(x) for x in oldlines[i].strip().split()[8:11]]) #SNF vector
                [a,b,c,d,e,f] = [int(x) for x in oldlines[i].strip().split()[11:17]] #From HNF matrix
                L = array([[int(x) for x in oldlines[i].strip().split()[17:20]], [int(x) for x in oldlines[i].strip().split()[20:23]], [int(x) for x in oldlines[i].strip().split()[23:26]]]) #Left transform
                aBas = array(zeros((3, nD*a*c*f))) #The base vectors for sites in the structure
                gIndx = [0]*(nD*a*c*f) #Mapping function showing the site that goes with each structure in the labeling

                ic = -1  #counter
                for iD in range(nD): # Loop over the number at sites/parent cell (the d set)
                    for z1 in range(a):
                        for z2 in range((b*z1)/a, c+(b*z1)/a):
                            for z3 in range(z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c):
                                ic = ic + 1
                                aBas[:,ic]=dot(pLV,array([z1,z2,z3]))+pBas[iD,:] # Atomic basis vector in Cartesian coordinates
                                greal = dot(L,array([z1,z2,z3])) # Map position into the group.
                                g = greal.astype(int) # Convert the g-vector from real to integer
                                g = fmod(g,S) # Bring the g-vector back into the first tile.
                                gIndx[ic] = iD*S[0]*S[1]*S[2]+g[0]*S[1]*S[2]+g[1]*S[2]+g[2]

                print aBas
                forbid = []
                HNF = array([[a,0,0],[b,c,0],[d,e,f]])
                LV = dot(pLV,HNF)


                for vecNum1 in range(nD*a*c*f - 1):
                    for vecNum2 in range(vecNum1 + 1, nD*a*c*f):
                        for i1 in range(-1,2):
                            for i2 in range(-1,2):
                                if linalg.norm(aBas[:,vecNum2]-aBas[:,vecNum1] + i1*LV[:,1] + i2*LV[:,2]) < 1.9:
                                    forbid.append([int(gIndx[vecNum1]), int(gIndx[vecNum2])])
                newHNF = False

            oldline = oldlines[i]
            label = oldline.strip().split()[-1]
            skip = False
            if enumVcNum != -1 and len(set([char for char in label])) != 1:
                for pair in forbid:
                    if label[pair[0]] != str(enumVcNum) and label[pair[1]] != str(enumVcNum):
                        skip = True
                        break
            if skip == True:
                continue
            else:
                newline = '%11i%10i%8i%9i%9i%12i%4i%6i%4i%3i%3i%5i%3i%3i%3i%3i%3i%7i%5i%5i%5i%5i%5i%5i%5i%5i' % tuple([structNum] + [int(x) for x in oldline.strip().split()[1:-1]]) + '   ' + oldline.strip().split()[-1] + '\n'
                newfile.write(newline)
                structNum += 1
        
        newfile.close()

changeEnumFile()
