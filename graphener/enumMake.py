import subprocess, time
from numpy import *
from math import sqrt, ceil
from copy import deepcopy

enumFile = 'struct_enum.out'

newfile = open(enumFile, 'w')

newfile.write('''graphene  
surf
 1000.00000       0.0000000       0.0000000        # a1 parent lattice vector
  0.0000000       2.1312885      -1.2305000        # a2 parent lattice vector
  0.0000000       2.1312885       1.2305000        # a3 parent lattice vector
    6 # Number of points in the multilattice
 0.10000000E+00   0.0000000       0.0000000        # d01 d-vector, labels: 0/1
 0.10000000E+00   1.4208590       0.0000000        # d02 d-vector, labels: 0/1
 0.10000000E+00   2.8417180       0.0000000        # d03 d-vector, labels: 0/1
 0.10000000E+00   2.1312885       0.0000000        # d04 d-vector, labels: 0/1
 0.10000000E+00   1.0656443     -0.61525000        # d05 d-vector, labels: 0/1
 0.10000000E+00   1.0656443      0.61525000        # d06 d-vector, labels: 0/1
 2-nary case
   1   4 # Starting and ending cell sizes for search
 0.10000000E-02 # Epsilon (finite precision parameter)
Concentration check:
    F
full list of labelings (including incomplete labelings) is used
Equivalency list: 1  2  3  4  5  6
start   #tot      HNF     Hdegn   labdegn   Totdegn   #size idx    pg    SNF             HNF                 Left transform                          labeling
          1         1       1        1        1          16   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   000000
          2         1       1        1        1          16   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   011111
          3         1       1        1        1          16   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   101111
          4         1       1        2        2           8   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   110111
          5         1       1        1        1          12   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   111011
          6         1       1        3        3          14   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   111101
          7         1       1        3        3          15   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   111110
          8         1       1        1        1          16   1    24   1  1  1    1  0  1  0  0  1      1    0    0    0    1    0    0    0    1   111111
''')

structs = [[1,1,'010111'], [1,1,'011011'], [1,1,'011111'], [1,1,'101111'], [1,1,'110111'], [1,1,'111011'], [1,1,'111101'], [1,1,'111110'], [1,1,'111111']]

nD = 6
pLV = array([[1000, 0, 0],[0, 2.1312885, 2.1312885],[0,-1.2305000,1.2305000]])
pBas = array([[0.10000000,0.0000000,0.0000000], [0.10000000,1.4208590,0.0000000], [0.10000000, 2.8417180, 0.0000000], [0.10000000, 2.1312885, 0.0000000], [0.10000000,1.0656443,-0.61525000], [0.10000000,3.1969328,-0.61525000]])

structNum = 9
enumVcNum = 1
HNFNum = 1

def getHNFs(volume):
    HNFs = []
    Rs = [array([[1., 0., 0.], [0., -0.5, 0.8660254], [0., 0.86602541, 0.5]]), array([[1., 0., 0.], [ 0., -1.,  0.], [ 0.,  0., -1.]]), array([[1., 0., 0.], [0., -0.5, -0.8660254 ], [0., 0.86602541, -0.5]]), array([[1., 0., 0.], [0., 0.5, 0.8660254], [0., 0.86602541, -0.5]]), array([[1., 0., 0.], [0., -1., 0.], [0., 0., 1.]]), array([[1., 0., 0.], [0., -0.5, 0.8660254], [0., -0.86602541, -0.5]]), array([[1., 0., 0.], [0., 0.5, -0.8660254], [0., 0.86602541, 0.5]]), array([[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]]), array([[1., 0., 0.], [0., -0.5, -0.8660254], [0., -0.86602541, 0.5]]), array([[1., 0., 0.], [0., 0.5, 0.8660254], [0., -0.86602541, 0.5]]), array([[1., 0., 0.], [0., 0.5, -0.8660254], [0., -0.86602541, -0.5]])]
    pLV = array([[1000, 0, 0],[0, 2.1312885, 2.1312885],[0,-1.2305000,1.2305000]])

    a = 1
    for c in range(1, volume + 1):
        if float(c) > float(volume)**0.5 or float(volume) / float(c) != int(float(volume) / float(c)):
            continue
        f = int(volume / c)
        for b in range(0, 1):
            for d in range(0, 1):
                for e in range(0, f):
                    HNFs.append(array([[a,0,0],[b,c,0],[d,e,f]]))
    #print HNFs
    validHNFs = [HNFs[0]]
    for HNF in HNFs[1:]:
        skip = False
        for testHNF in deepcopy(validHNFs):
            #print 'HNF', HNF
            #print 'test HNF', testHNF
            for R in Rs:
                L1 = array([[1,0,0],[0,1,0],[0,-(testHNF[2,1]),1]])
                L2 = array([[1,0,0],[0,1,0],[0,-(HNF[2,1]),1]])
                #H = linalg.inv(pLV*testHNF)*linalg.inv(R)*(pLV*HNF)
                #H = linalg.inv((pLV*L1)*testHNF)*linalg.inv(R)*(pLV*L2*HNF)
                H = dot(dot(linalg.inv(dot(pLV,testHNF)),linalg.inv(R)),dot(pLV,HNF))
                #print 'R', R
                #print 'H', H
                valid = False
                for n in range(3):
                    for m in range(3):
                        #if not equal(mod(H[n,m], 1), 0):
                        if mod(H[n,m], 1) >= 0.0001 and mod(H[n,m], 1) <= 0.9999:
                            valid = True
                if valid == False:
                    skip = True
                    #break
            if skip == True:
                break
        if skip == False:
            validHNFs.append(HNF)
    return validHNFs
  

#for n in range(2,9): 
#    print n
#    HNFs = getHNFs(n) 
#    for HNF in HNFs:
#        print HNF              

for volume in range(2, 9):
    subprocess.call(['echo','\tGenerating structures for volume {}.\n'.format(volume)])
    HNFs = getHNFs(volume)
    for HNF in HNFs:
        a = HNF[0,0]
        b = HNF[1,0]
        c = HNF[1,1]
        d = HNF[2,0]
        e = HNF[2,1]
        f = HNF[2,2]
        S = array([a, c, f])
        L = array([[1,0,0],[0,1,0],[0,-e,1]])

        HNFNum += 1

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

        forbid = []
        LV = dot(pLV,HNF)

        for vecNum1 in range(nD*a*c*f - 1):
            for vecNum2 in range(vecNum1 + 1, nD*a*c*f):
                for i1 in range(-1,2):
                    for i2 in range(-1,2):
                        if linalg.norm(aBas[:,vecNum2]-aBas[:,vecNum1] + i1*LV[:,1] + i2*LV[:,2]) < 1.9:
                            forbid.append([int(gIndx[vecNum1]), int(gIndx[vecNum2])])

        chooseStructs = [0]*volume
        finished = False
        while not finished:
            chooseStructs[volume - 1] += 1
            structsReady = False
            while not structsReady:
                if 9 in chooseStructs:
                    i = chooseStructs.index(9)
                    chooseStructs[i] = 0
                    if i > 0:
                        chooseStructs[i-1] += 1
                    else:
                        finished = True
                        break
                else:
                    structsReady = True

            skip = False
            transStructs = []
            smallFactor = c
            largeFactor = f
            for j in range(0, volume - smallFactor + 1, smallFactor):
                if j != 0:
                    transStructs.append([chooseStructs[(i + j) % volume] for i in range(0, volume)])
            for j in range(0, volume - largeFactor + 1, largeFactor):
                if j != 0:
                    transStructs.append([chooseStructs[(i + j) % volume] for i in range(0, volume)])
            for testStruct in transStructs:
                if int(''.join([str(x) for x in testStruct])) < int(''.join([str(x) for x in chooseStructs])):
                    skip = True
                    break

            if skip == True:
                continue

            if finished == True:
                break
            instructs = [structs[i][2] for i in chooseStructs]

            label = ''
            for i in range(6):
                label += ''.join([instruct[i] for instruct in instructs])
            if enumVcNum != -1 and len(set([char for char in label])) != 1:
                for pair in forbid:
                    if label[pair[0]] != str(enumVcNum) and label[pair[1]] != str(enumVcNum):
                        skip = True
                        break
                if skip == True:
                    continue
                else:
                    newline = '%11i%10i%8i%9i%9i%12i%4i%6i%4i%3i%3i%5i%3i%3i%3i%3i%3i%7i%5i%5i%5i%5i%5i%5i%5i%5i' % tuple([structNum, HNFNum, 0,0,0,0, volume, 0, S[0],S[1],S[2], a,b,c,d,e,f, 1,0,0,0,1,0,0,L[2,1],1]) + '   ' +  label + '\n'
                    newfile.write(newline)
                    structNum += 1

newfile.close()

subprocess.call(['./remove_duplicates_enum_file.x'])

#Re-number the structures
infile = open('struct_enum_unique.out', 'r')
inlines = infile.readlines()
infile.close()

outfile = open('struct_enum_unique2.out', 'w')
i = 0
for line in inlines:
    outfile.write(line)
    i += 1
    if line.strip().split()[0] == 'start':
        break

for n, line in enumerate(inlines[i:]):
    newline = '%11i' % (n+1) + line[11:]
    outfile.write(newline)

outfile.close()
    












