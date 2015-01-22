'''Run DOS for all structure folders in maindir'''
import os, subprocess,sys,re
from comMethods import convergeCheck,finishCheck,getNSW,getSteps,readfile,writefile
from numpy import rint
import PlotStructures
from PlotStructures import plotByPrior
   
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
minPrior = 10.0
iteration = 4
plotByPrior(atoms,minPrior,iteration)

print "Done"

