'''
Created on Aug 29, 2014

@author: eswens13
'''
import os, subprocess, sys
from numpy import amax, amin, zeros, sort, array, floor, exp, ceil, median, int32
from copy import deepcopy
from comMethods import *
class GSS:
    """ This class performs the UNCLE ground state search for the lowest formation energy
        structures. From this, we get a list of every structure that was enumerated and its
        corresponding formation energy as predicted by UNCLE. We use this list to decide which new
        structures to add into the model. We also keep track of the list and the plots of the list
        from iteration to iteration. """

    def __init__(self, atoms, volRange, plotTitle, xlabel, ylabel, vstructsFinished, uncleOutput):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.volRange = volRange
        self.plotTitle = plotTitle
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.vstructsFinished = vstructsFinished
        
        self.enumFolder = os.getcwd() + '/enum/'
        self.fitsDir = None
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
        self.Ncs = []
        self.priorities = []

    def contains(self, struct, alist):
        """ Returns true if 'struct' is found in 'alist', false otherwise. """
        if len(alist) == 0:
            return False
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
        return False

#    def getAllGSSStructures(self, iteration, newlyFailed):
#        """ Returns a list of all the structures sorted by their predicted formation energies.
#            It does this for each metal atom that has been specified by the user so this will 
#            actually return a list of lists--a list for each atom. """
#        allStructs = []
#        for n in xrange(len(self.atoms)):
#            atomStructs = []
#            structsByEnergy = []
#            gssFile = os.getcwd() + '/' + self.atoms[n] + '/gss/gss_' + str(iteration) + '.out'
#            infile = open(gssFile, 'r')
#            
#            for i,line in enumerate(infile):
#                if i >= 2:
#                    formEnergy = float(line.strip().split()[7])
#                    struct = int(line.strip().split()[0])                    
#                    # Do not include structures that failed VASP calculations.
#                    if not self.contains(struct, newlyFailed[n]):
#                        structsByEnergy.append([formEnergy, struct])
#            infile.close()
#            
#            structsByEnergy.sort()
#            
#            for struct in structsByEnergy:
#                atomStructs.append(str(struct[1]))
#            
#            allStructs.append(atomStructs)
#            
#        return allStructs
    
    def getGssInfo(self,iteration,vstructsFailed):   
        '''Get structure,Uncle formation energy, concentration into gssInfo.
        Write out gssFailedVasp.out, which lists structures and their model formation energy so they 
        can be identified on the gss plot as a different color'''
        lastDir = os.getcwd()
        dir = lastDir
        self.getNcs(iteration) #number at each concentration
        numberCs = len(self.Ncs)
        Ntot = sum(self.Ncs) #total number of structures              
        if iteration==1: subprocess.call(['echo',  'Number of concentrations: '+ str(numberCs)])     
        self.priorities = zeros((len(self.atoms),Ntot),dtype = [('struct', int32),('FE', float), ('prior', float)])
#        e_cutoff = zeros(numberCs,dtype = float)
        for iatom, atom in enumerate(self.atoms):
            if len(self.vstructsFinished[iatom]) > 1:
                atomDir = dir + '/' + self.atoms[iatom] + '/gss'
                gfvfile = open(atomDir + '/gssFailedVasp.out','w')
                gssInfo = zeros((Ntot),dtype = [('struct', 'S10'), ('conc', float), ('FE', float)]) #columns:  struct, concentration, energy
                lines = readfile(atomDir + '/gss_' + str(iteration) + '.out')
                for i,line in enumerate(lines[2:]): #first 2 lines are headers
                    struct = line.strip().split()[0]
                    conc = float(line.strip().split()[2]) #metal concentration
                    formEnergy = float(line.strip().split()[7])                
                    gssInfo[i-2]['struct'] = struct
                    gssInfo[i-2]['conc'] = conc
                    gssInfo[i-2]['FE'] = formEnergy
                    if struct in vstructsFailed[iatom]: 
                        gfvfile.write('{:10s}{:12.8f}{:12.8f}\n'.format(struct,conc,formEnergy))
                gfvfile.close()            
                gssInfo = sort(gssInfo, order=['conc','FE']) #sorts low to high
                emin = amin(gssInfo[:]['FE'])
                emed = median(gssInfo[:]['FE'])
                width_percent = 0.001   #weight the bottom % strongly
                iplace = 0
                for ic,Nc in enumerate(self.Ncs):
                    imin = iplace #place of lowest energy for this concentration
                    eminc = gssInfo[imin]['FE'] #lowest energy for this concentration
                    width = ceil(width_percent*Nc) + 2*width_percent*Nc*max(0,(emed-eminc)/(emed-emin)) #last term favors the lowest energy structures globally.  Is zero when eminc is at or above the median energy for this atom
                    for n in range(Nc):
                        istr = iplace+n
                        self.priorities[iatom,istr]['struct'] = gssInfo[istr]['struct']
                        self.priorities[iatom,istr]['FE'] = gssInfo[istr]['FE']
                        self.priorities[iatom,istr]['prior'] = 100 * exp(-(istr-imin)/width) 
                    iplace += Nc
    #            self.priorities = sort(self.priorities, order=['prior']) # sorted low to high 
                #Note:  if I do the sort above, the sort works in the sense that when I print the elements,
                # they are sorted.  BUT, when I write them to a file they are all zero!  Some memory/pointer problem? 
                # so I will use the work-around with linux sort the written file, then read them back in 
                priorfile = open(atomDir+'/temp','w')
                priorfile.write('structure,priority,concentration,FEnergy\n')
                for i in range(Ntot):
                    priorfile.write('{:10d}{:12.6f}{:8.4f}{:10.6f}\n'.format( \
                        self.priorities[iatom,i]['struct'],self.priorities[iatom,i]['prior'], \
                           gssInfo[i]['conc'],gssInfo[i]['FE']))               
                priorfile.close()
                os.chdir(atomDir)
                os.system('sort -g -r -k 2 temp > '+ 'priorities_{}.out'.format(iteration)) #puts highest priorites at top
                lines = readfile(atomDir+'/priorities_{}.out'.format(iteration))
                for i, line in enumerate(lines[:-1]): #skip last line which is header (now footer) as sorted
                    data = line.strip().split()
                    self.priorities[iatom,i]['struct'] = data[0]
                    self.priorities[iatom,i]['prior'] = data[1]
                    self.priorities[iatom,i]['FE'] = data[3]
                os.chdir(lastDir)                  
        os.chdir(lastDir)
        return self.priorities                
            
    def getNcs(self,iteration): #bch
        '''Find the number of structures at each concentration''' 
        lastDir = os.getcwd()
        dir = os.getcwd() + '/' + self.atoms[0] + '/gss/'
        os.chdir(dir)
        os.system('sort -g -k 3 gss_1.out > tempout') #sorts by column 3, metal concentration
        lines = readfile('tempout');os.system('rm tempout')         
        conc_old = 1.0 #starts with highest concentration first 
        Nc = 0  
        if iteration==1: concfile = open(lastDir + '/concentrations.info' , 'w')
        for line in lines[2:]:
            conc = float(line.strip().split()[1])
            if conc == conc_old:
                Nc += 1
            else: 
                if iteration==1: concfile.write('{:8.3f}{:10d}\n'.format(conc_old,Nc))
                self.Ncs.append(Nc) #for the previous concentration           
                Nc = 1
                conc_old = conc
        self.Ncs.append(1) #for the pure H concentration
        if iteration==1: concfile.close()
        os.chdir(lastDir)
                        
    def makeGSSDirectories(self):
        """ Creates the 'gss' directories for each different metal atom and populates them with 
            the files that UNCLE needs to perform a ground state search.  These files are 
            struct_enum.out, lat.in, J.1.out, groundstatesearch.in, and gss_plot.gp. """
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                self.fitsDir = atomDir + '/fits'
                gssDir = atomDir + '/gss'
                if not os.path.exists(gssDir):
                    subprocess.call(['mkdir',gssDir])
                subprocess.call(['cp',self.enumFolder + '/struct_enum.out', gssDir])
                subprocess.call(['cp',self.enumFolder + '/lat.in', gssDir])
                subprocess.call(['cp',self.fitsDir + '/J.1.out', gssDir])
           
                inlines = readfile(self.neededFilesDir + '/groundstatesearch.in')              
                outfile = open(gssDir + '/groundstatesearch.in','w')
                for i in xrange(len(inlines)):
                    if i == 4:
                        outfile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + "\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
                
                inlines = readfile(self.neededFilesDir + '/gss_plot.gp')
                outfile = open(gssDir + '/gss_plot.gp','w')
                for i in xrange(len(inlines)):
                    if i == 3:
                        outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atom) + "\"\n")#bch
                    elif i == 4:
                        outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                    elif i == 5:
                        outfile.write("set title \"" + self.plotTitle + " (" + atom + ")\"\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
                

    def makePlots(self, iteration):
        """ Creates the plots of the predicted energies of all the structures that have been 
            enumerated. Adds the iteration number onto the end of the filenames for the plots and
            the lists. """
        lastDir = os.getcwd() 
        for iatom, atom in enumerate(self.atoms):
            if len(self.vstructsFinished[iatom]) > 1:
                subprocess.call(['echo', '\nMaking plots for ' + atom + '. . .\n'])
                gssDir = os.getcwd() + '/' + atom + '/gss'
                os.chdir(gssDir)                                   
                subprocess.call(['mv', '../vaspFE.out', '.'])
                subprocess.call(['mv', '../vaspBE.out', '.'])
                subprocess.call(['mv', '../vaspHFE.out', '.']) 
                if os.path.isdir(gssDir):           
                    #Binding energy
                    inlines = readfile(self.neededFilesDir + '/BE_plot.gp')
                    outfile = open(gssDir + '/BE_plot.gp','w')
                    for i in xrange(len(inlines)):
                        if i == 3:
                            outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atom) + "\"\n")#bch
                        elif i == 4:
                            outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                        elif i == 5:
                            outfile.write("set title \"" + 'Binding energy vs graphene'+ " (" + atom + ")\"\n")
                        else:
                            outfile.write(inlines[i])
                    outfile.close()
            
                    #Hexagonal formation energy    
                    data = readfile(gssDir+'/vaspHFE.out')
                    ydata = [float(line.strip().split()[1]) for line in data]
                    ymax = min(amax(ydata),1.0) #don't let yrange get over 1 eV
                    inlines = readfile(self.neededFilesDir + '/HFE_plot.gp')
                    outfile = open(gssDir + '/HFE_plot.gp','w')
                    for i in xrange(len(inlines)):
                        if i == 3:
                            outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atom) + "\"\n")#bch
                        elif i == 4:
                            outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                        elif i == 5:
                            outfile.write("set title \"" + 'Formation energy vs H2, metal hex monolayer'+ " (" + atom + ")\"\n")
                        elif 'plot "' in inlines[i]:         
                            outfile.write('set yrange [:{}]\n'.format(ymax))
                            outfile.write(inlines[i])
                        else:
                            outfile.write(inlines[i])
                    outfile.close()
                
                    subprocess.call(['gnuplot', 'gss_plot.gp'])
                    subprocess.call(['gnuplot', 'BE_plot.gp'])
                    subprocess.call(['gnuplot', 'HFE_plot.gp'])
                    subprocess.call(['mv','gss.out','gss_' + str(iteration) + '.out'])
                    subprocess.call(['mv','gss.pdf','gss_' + str(iteration) + '.pdf'])
                    os.chdir(lastDir)
    
    def performGroundStateSearch(self, iteration):
        """ Performs the ground state search with the current fit from UNCLE. """
        lastDir = os.getcwd()
        natoms = len(self.atoms)
        subdir = 'gss'
        for iatom, atom in enumerate(self.atoms):
            if len(self.vstructsFinished[iatom]) > 1:
                gssDir = os.getcwd() + '/' + atom + '/gss/'
                if os.path.isdir(gssDir):
                    subprocess.call(['echo','\nPerforming ground state search for ' + atom + '. . .\n'])
                    os.chdir(gssDir)
                    if os.path.exists('gss.out'): subprocess.call(['rm','gss.out'])
                    os.chdir(lastDir)
        if natoms ==1:
            os.chdir(lastDir + '/' + self.atoms[0]  + '/' + subdir)
            subprocess.call([self.uncleExec, '21'], stdout=self.uncleOut)             
            os.chdir(lastDir)   
        else:#parallelize the atom jobs
            #make job files
            mem = '16' #Gb
            walltime = 1.0 #hrs
            execString = self.uncleExec + ' 21'
            parallelJobFiles(self.atoms,subdir,walltime,mem,execString)
            #submit jobs
            jobIds = parallelAtomsSubmit(self.atoms[1:],subdir)
            #use this job to calculate the first atom:
            os.chdir(lastDir + '/' + self.atoms[0]  + '/' + subdir)
            subprocess.call(['echo','\tThis job calculating the first atom: {}. Submitted jobs for the others.\n'.format(self.atoms[0])])
            subprocess.call([self.uncleExec, '21'], stdout=self.uncleOut)             
            os.chdir(lastDir)
            #wait for others
            parallelAtomsWait(jobIds)
            os.chdir(lastDir)             
                          







        