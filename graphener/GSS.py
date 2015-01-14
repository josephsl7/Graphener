'''
Created on Aug 29, 2014


'''
import os, subprocess, sys
from numpy import amax, amin, zeros, sort, array, floor, exp, ceil, median, int32, mod
from copy import deepcopy
from comMethods import *
class GSS:
    """ This class performs the UNCLE ground state search for the lowest formation energy
        structures. From this, we get a list of every structure that was enumerated and its
        corresponding formation energy as predicted by UNCLE. We use this list to decide which new
        structures to add into the model. We also keep track of the list and the plots of the list
        from iteration to iteration. """
    from comMethods import setAtomCounts,hexMonolayerEnergies,singleAtomsEnergies,getPureEs,\
                    contains,elConvergeCheck,electronicConvergeFinish,getElSteps,getEnergy,setEnergy

    def __init__(self, atoms, volRange, plotTitle, xlabel, ylabel, vstructsFinished, uncleOutput,pureMetal,finalDir):
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
        
        self.pureHenergy = 0.0
        self.pureMenergy = 0.0
        self.pureMetal = pureMetal
    
        self.energy = 0.0   
        self.singleE = [] 
        self.hexE = [] 
        self.finalDir = finalDir 

    def collate_plots(self,plotType,iteration):  
        '''Creates an HTML page with the plots and labels. plotType: gss,BE,HFE'''
        lastDir = os.getcwd()
        plots1Dir = lastDir + '/plots'
        if not os.path.exists(plots1Dir):
            subprocess.call(['mkdir',plots1Dir]) 
        plots2Dir = plots1Dir + '/plots{}'.format(plotType)
        if not os.path.exists(plots2Dir):
            subprocess.call(['mkdir',plots2Dir])        
        nRow = 5  # number of plots in row
        width = 1350
        height  = 900
        collatefile  = open(plots2Dir +'/plots{}_{}.htm'.format(plotType,iteration),'w')
        collatefile.write(' <html>\n <HEAD>\n<TITLE> {} </TITLE>\n</HEAD>\n'.format(plotType))
        collatefile.write(' <BODY>\n <p style="font-size:20px"> <table border="1">\n <tr>\n') #start row
    #    images = []
        iImage = 0
        for atom in self.atoms:
            path = lastDir + '/' + atom + '/gss/{}_{}.png'.format(plotType,iteration)
            iImage += 1  
            atomtext = atom.split('_')[0] 
            name = '{}{}.png'.format(atomtext,plotType)       
            subprocess.call(['cp',path,plots2Dir + '/{}'.format(name)])
    #        images.append()
            collatefile.write('<td><p><img src="{}" width "{}" height "{}" ></p><p>{}</p></td>\n'.format(name,width,height,''))#Image and element under it
            if mod(iImage,nRow) == 0: 
                collatefile.write('</tr>\n<tr>\n') #end of row, begin new
        collatefile.write(' </tr></table> \n') #end of row and table                
        collatefile.write(' </BODY> </html>') #end of file 
        collatefile.close()  
        
    def contains(self, struct, alist):
        """ Returns true if 'struct' is found in 'alist', false otherwise. """
        if len(alist) == 0:
            return False
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
        return False
    
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
        self.priorities = zeros((len(self.atoms),Ntot),dtype = [('struct', int32),('FE', float), ('conc', float), ('prior', float)])
#        e_cutoff = zeros(numberCs,dtype = float)
        for iatom, atom in enumerate(self.atoms):
            atomDir = lastDir + '/' + atom
            if len(self.vstructsFinished[iatom]) > 1:
                gssDir = dir + '/' + self.atoms[iatom] + '/gss'
                gfvfile = open(gssDir + '/gssFailedVasp.out','w')
                gssInfo = zeros((Ntot),dtype = [('struct', 'S10'), ('conc', float), ('FE', float)]) #columns:  struct, concentration, energy
                lines = readfile(gssDir + '/gss_' + str(iteration) + '.out')
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
                        self.priorities[iatom,istr]['conc'] = gssInfo[istr]['conc']
                        self.priorities[iatom,istr]['prior'] = 100 * exp(-(istr-imin)/width)  
                    iplace += Nc
                #sort highest to lowest
                self.priorities[iatom,:] = sort(self.priorities[iatom,:], order=['prior']) # sorted low to high 
                self.priorities[iatom,:]= self.priorities[iatom,::-1] # reversed: now highest to lowest
                os.chdir(lastDir)                  
                priorfile = open(atomDir+'/priorities_{}.out'.format(iteration),'w')  
                priorfile.write('structure,priority,concentration,FEnergy\n')
                for i in range(Ntot):
                    print i, self.priorities[iatom,i]['struct'],self.priorities[iatom,i]['prior'],gssInfo[i]['conc'],gssInfo[i]['FE']
                    priorfile.write('{:7d}   {:10.6f}{:8.4f}{:10.6f}\n'.format( \
                    self.priorities[iatom,i]['struct'],self.priorities[iatom,i]['prior'], \
                    self.priorities[iatom,i]['conc'],self.priorities[iatom,i]['FE']))
                priorfile.close()
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
#                        outfile.write("set title \"" + self.plotTitle + " (" + atom + ")\"\n")
                        outfile.write("set title '' \n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
                

    def makePlots(self, iteration):
        """ Creates the plots of the predicted energies of all the structures that have been 
            enumerated. Adds the iteration number onto the end of the filenames for the plots and
            the lists. """
        lastDir = os.getcwd() 
        self.writeUncleBE_HFE(iteration)
        for iatom, atom in enumerate(self.atoms):
            atomtext = atom.split('_')[0] #remove _suffixes
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
                            outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atomtext) + "\"\n")#bch
                        elif i == 4:
                            outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                        elif i == 5:
#                            outfile.write("set title \"" + 'Binding energy vs graphene'+ " (" + atom + ")\"\n")
                            outfile.write("set title \"\" \n")
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
                            outfile.write("set xlabel \"" + self.xlabel.replace('Metal',atomtext) + "\"\n")#bch
                        elif i == 4:
                            outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                        elif i == 5:
#                            outfile.write("set title \"" + 'Formation energy vs H2, metal hex monolayer'+ " (" + atomtext + ")\"\n")
                         outfile.write("set title \"\" \n")

                        elif 'plot "' in inlines[i]:         
#                            outfile.write('set yrange [:{}]\n'.format(ymax))
                            outfile.write(inlines[i])
                        else:
                            outfile.write(inlines[i])
                    outfile.close()
                
                    subprocess.call(['gnuplot', 'gss_plot.gp'])
                    subprocess.call(['gnuplot', 'BE_plot.gp'])
                    subprocess.call(['gnuplot', 'HFE_plot.gp'])
                    
                    subprocess.call(['convert -density 300 gss.pdf -resize 1800x2700 gss_' + str(iteration) + '.png'],shell = True)
                    subprocess.call(['convert -density 300 BE.pdf -resize 1800x2700 BE_' + str(iteration) + '.png'],shell = True)
                    subprocess.call(['convert -density 300 HFE.pdf -resize 1800x2700 HFE_' + str(iteration) + '.png'],shell = True)
#                    subprocess.call(['convert','-density','300','gss.pdf','-resize','1800x2700', 'gss_' + str(iteration) + '.png'],shell = True)
#                    subprocess.call(['convert','-density','300','BE.pdf','resize','1800x2700','BE_' + str(iteration) + '.png'])
#                    subprocess.call(['convert','-density','300','HFE.pdf','resize','1800x2700','HFE_' + str(iteration) + '.png'])                 

                    subprocess.call(['cp','gss.out','gss_' + str(iteration) + '.out'])
                    subprocess.call(['cp','vaspBE.out','vaspBE_' + str(iteration) + '.out'])
                    subprocess.call(['cp','vaspFE.out','vaspFE_' + str(iteration) + '.out'])
                    subprocess.call(['cp','vaspHFE.out','vaspHFE_' + str(iteration) + '.out'])
                    subprocess.call(['cp','uncleBE.out','uncleBE_' + str(iteration) + '.out'])
                    subprocess.call(['cp','uncleHFE.out','uncleHFE_' + str(iteration) + '.out'])

                    os.chdir(lastDir)
        self.collate_plots('gss',iteration)
        self.collate_plots('BE',iteration)
        self.collate_plots('HFE',iteration)
    
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
            walltime = 2.0 #hrs
            execString = self.uncleExec + ' 21'
            atomStrings = ['']*natoms
            parallelJobFiles(self.atoms,subdir,walltime,mem,execString,atomStrings) 
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

    def writeUncleBE_HFE(self,iteration):
        '''Using vasp atomic,moleculat and vasp hex monolayer reference energies, find the fitted
        values for binding energy (vs atomic) and HFE (vs hexagonal monolayer and H2)'''
        lastDir = os.getcwd()
        self.singleAtomsEnergies(os.getcwd(),iteration)
        self.hexMonolayerEnergies(os.getcwd(),iteration)   
        eIsolatedH = -1.1073   
        eIsolatedC = -1.3179 
        eH2 = -6.7591696/2.0 
        energyGraphene = -18.456521 #for 2 C atoms         
 
        for iatom, atom in enumerate(self.atoms):
            self.getPureEs(iatom)
            gssDir = lastDir + '/' + atom + '/gss'
            lines = readfile(gssDir + '/gss.out')
#            lines = readfile(gssDir + '/gss_' + str(iteration) + '.out')
            uncleBEfile = open(gssDir + '/uncleBE.out','w')  
            uncleHFEfile = open(gssDir + '/uncleHFE.out','w')
#            #get the volume factor of structure 3
#            vol3 = lines[5].strip().split()[3]
#            nsites = 3-vol3 #number of sites in smallest structure, either 1 or 2
            for i,line in enumerate(lines[2:]): #first 2 lines are headers
                struct = int(line.strip().split()[0])
                x = float(line.strip().split()[2]) #metal concentration
                vol  = int(line.strip().split()[3])
                nH  = int(line.strip().split()[4])
                nmetal  = int(line.strip().split()[5])
                nadatoms = nH + nmetal
                
                ncarbon = vol * 2  #always 2 C atoms in smallest cell
                
                FE = float(line.strip().split()[7]) #formation energy
                structEnergy = FE + x*self.pureMenergy + (1-x)*self.pureHenergy #uncle fit energy per adatom
                structEnergy = structEnergy * nadatoms #total
                bindEnergy = (structEnergy - ncarbon*energyGraphene/2.0 - nH*eIsolatedH - nmetal*self.singleE[iatom])/ float(nadatoms) #2 atoms in graphene 
                uncleBEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,x,bindEnergy))
                hexFormationEnergy = (structEnergy - energyGraphene*ncarbon/2.0  - nmetal *self.hexE[iatom] - nH*eH2)/float(nadatoms)
                uncleHFEfile.write('{:10d} {:12.8f} {:12.8f}\n'.format(struct,x,hexFormationEnergy))    
            uncleBEfile.close()
            uncleHFEfile.close() 




        