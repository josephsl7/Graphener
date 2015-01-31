'''Run DOS for all structure folders in maindir'''
import os, subprocess,sys,re
from comMethods import convergeCheck,finishCheck,getNSW,getSteps,readfile,writefile
from numpy import rint
import PlotStructures
from PlotStructures import collateStructsHFE
   
#======================================= Script =====================================
#======================================= Script =====================================
#======================================= Script =====================================

maindir = '/fslhome/bch/cluster_expansion/graphene/top.tm_row1.v15/'
os.chdir(maindir)
atoms = []
for item in os.listdir(maindir):
    itempath = maindir + '/' + item
    if os.path.isdir(itempath) and item[0].isupper(): #look only at dirs whose names are capitalized. These are the atoms
        atoms.append(item)
print atoms
NInPlot = 400
iteration = 4
minPrior = 0.01
collateStructsHFE(atoms,minPrior,NInPlot,iteration)  

print "Done"