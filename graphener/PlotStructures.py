import sys, os, subprocess
import matplotlib.pyplot as plt
import numpy as np
from comMethods import *
import StructsToPoscar
from matplotlib.lines import Line2D 
from copy import deepcopy

def collateStructsConc(atoms,minPrior,iteration):  
    '''Creates an HTML page with structure plots for each concentration, with the highest priority
    structs at the top.  Labels include atom, conc, FE, priority. Store in /struct_plots within each atom folder.
     '''
    lastDir = os.getcwd()         
    nRow = 4  # number of plots in row
    width = 500
    height  = 500
    iImage = 0
    for iatom, atom in enumerate(atoms):
        print atom
        atomDir = lastDir + '/' + atom 
        plotsDir = atomDir + '/struct_plots'       
        if not os.path.exists(plotsDir):
            subprocess.call(['mkdir',plotsDir])
        plines = readfile(atomDir + '/priorities_{}.out'.format(iteration))
        if iatom == 0: #get list of concentrations 
            conclist = []
            for line in plines[1:]:
                conc = line.split()[2]
                if conc not in conclist:
                    conclist.append(conc)
            subprocess.call(['echo','Found {} concentrations'.format(len(conclist))])
        for conc in conclist:
            print '   ',conc
            collatefile  = open(plotsDir +'/{}_{}.htm'.format(conc,iteration),'w')
            collatefile.write(' <html>\n <HEAD>\n<TITLE> {} Conc {}, Iteration {} </TITLE>\n</HEAD>\n'.format(atom,conc,iteration))
            collatefile.write(' <BODY>\n <p style="font-size:20px"> <table border="1">\n <tr>\n') #start row
            iImage = 0 
            for line in plines[1:]:
                struct = line.split()[0] 
                conc2 = line.split()[2]
                prior = line.split()[1]
                path = '/structs/{}.png'.format(struct)
                if conc2 == conc and float(prior) >= minPrior and os.path.exists(plotsDir+'/'+path):
                    FE = line.split()[3]
                    label = '  <b>{}</b>    {}: <b>FE</b> {}, prior {}, <b>conc</b> {}'.format(struct,atom,round(FE,2),prior,conc)           
                    collatefile.write('<td><p><img src="../../{}" width "{}" height "{}" ></p><p>{}</p><p></p></td>\n'.format(path,width,height,label,))#Image and element under it
                    iImage += 1 
                    if mod(iImage,nRow) == 0: 
                        collatefile.write('</tr>\n<tr>\n') #end of row, begin new
            collatefile.write(' </tr></table> \n') #end of row and table                
            collatefile.write(' </BODY> </html>') #end of file 
            collatefile.close()

def collateStructsHFE(atoms,minPrior,NInPlot,iteration):  
    '''Creates an HTML page with structure plots ordered by hexagonal layer formation energy HFE, with the lowest (most negative) HFE first
    structs at the top. If minPrior >0, then it will give only the higher priority ones with lower HFE.  Labels include atom, conc, vHFE, vFE, priority. Store in /struct_plots within each atom folder/
    Calculates an error between vasp and uncle HFEs for each plot, and globally for these selected structures'''
    
    lastDir = os.getcwd()         
    nRow = 4  # number of plots in row
    width = 500
    height  = 500
    iImage = 0
    for iatom, atom in enumerate(atoms):
        print atom
        atomDir = lastDir + '/' + atom 
        plotsDir = atomDir + '/struct_plots'
        if not os.path.exists(plotsDir):     
            subprocess.call(['mkdir',plotsDir])
        plines = readfile(atomDir + '/priorities_{}.out'.format(iteration))
        hlines = readfile(atomDir + '/gss/vaspHFE_{}.out'.format(iteration))
        flines = readfile(atomDir + '/gss/vaspFE_{}.out'.format(iteration))
        nVaspCalc = len(hlines)
        cdata = zeros(len(hlines),dtype = [('struct', int32),('conc', float), \
            ('uFE', float),('vFE', float),('vHFE', float), ('fiterr', float),('prior', float)]) #vFE is formation energy from vasp, uFE is from uncle fit 
        for istruct,line in enumerate(hlines):
            cdata[istruct]['struct'] = int(line.split()[0])
            cdata[istruct]['vHFE'] = line.split()[2]
        structlist = cdata['struct'].tolist()
        for line in flines:
            struct = int(line.split()[0])
            cdata[structlist.index(struct)]['vFE'] = line.split()[2]
            cdata[structlist.index(struct)]['fiterr'] = cdata[structlist.index(struct)]['uFE']-cdata[structlist.index(struct)]['vFE']
        for i,line in enumerate(plines[1:]):
            struct = int(line.split()[0]) 
            try: #is istruct in structlist?
                istruct = structlist.index(struct)
                cdata[istruct]['prior'] = line.split()[1]
                cdata[istruct]['conc'] = line.split()[2]
                cdata[istruct]['uFE'] = line.split()[3]
            except:
                'do nothing'
        vcdata = sort(cdata,order = ['vHFE']) #lowest first
#        ucdata = sort(cdata,order = ['uFE']) 
        collatefile  = open(plotsDir +'/HFEsort_{}.htm'.format(iteration),'w')
        collatefile.write(' <html>\n <HEAD>\n<TITLE> {} HFE sort, Iteration {} </TITLE>\n</HEAD>\n'.format(atom,iteration))
        collatefile.write(' <BODY>\n <p style="font-size:20px"> <table border="1">\n <tr>\n') #start row
        iImage = 0
        istruct = 0
        err = 0.0
        err2 = 0.0
        errOr = 0.0
        N = min(NInPlot,nVaspCalc) 
        while iImage <=N and istruct <= nVaspCalc-1:
            struct = vcdata[istruct]['struct']
            path = 'structs/{}.png'.format(struct)
            if vcdata[istruct]['prior'] >= minPrior:
                if os.path.exists(lastDir+'/' + path):
#                    orderErr = errOrder(vcdata,ucdata,istruct)
                    orderErr = errOrderConc(vcdata,istruct)
    #                subprocess.call(['echo','Found plot for {}'.format(struct)])
                    label1 = '  <b>{}</b>  {}: <b>vHFE</b> {}, <b>vFE</b> {}, <b>OrderErr</b> {}'\
                        .format(struct,atom,round(vcdata[istruct]['vHFE'],3),round(vcdata[istruct]['vFE'],3),\
                                round(orderErr,4))           
                    label2 = '  <b>uPrior</b> {}, <b>Conc</b> {}'\
                        .format(round(vcdata[istruct]['prior'],3),vcdata[istruct]['conc']) 
                    collatefile.write('<td><p><img src="../../{}" width "{}" height "{}" ></p><p>{}</p><p>{}</p><p></p></td>\n'.format(path,width,height,label1,label2))#Image and element under it
                    iImage += 1 
                    if mod(iImage,nRow) == 0: 
                        collatefile.write('</tr>\n<tr>\n') #end of row, begin new
                    err += vcdata[istruct]['fiterr']**2
                    err2 += vcdata[istruct]['fiterr']
                    errOr += orderErr
                else:
                    subprocess.call(['echo','No plot found for structure {}, prior {}'.format(struct,vcdata[istruct]['prior'])])
            istruct += 1
        err = sqrt(err/iImage)
        err2 = err2/iImage
        errOr = errOr/iImage
        subprocess.call(['echo','RMS error for these structures: {}'.format(round(err,4))])
        subprocess.call(['echo','Average order error of fit vs calc for these structures: {}'.format(round(errOr,4))])
        subprocess.call(['echo','Average shift of fit vs calc for these structures: {}\n'.format(round(err2,4))])
        
        collatefile.write(' </tr></table> \n') #end of row and table                
        collatefile.write(' <p><b>RMS error for these structures:</b> {},</p>'.format(round(err,4))) #end of file 
        collatefile.write(' </BODY> </html>') #end of file 
        collatefile.close()

def errOrder(vcdata,ucdata,vindex):
    '''errOrder = |uncleFE(vaspindex) - uncleFE(uncleindex)|, where uncleFE is an ordered list containing only the vasp structures we have calculated.
    uncleindex is the index of the structure in uncleFE. vaspindex is the index of the structure in cdata, which must be ordered with the lowest energy at the top'''
    struct = vcdata[vindex]['struct']
    ustructlist = ucdata['struct'].tolist()
    uindex = ustructlist.index(struct)   
    return abs(ucdata[vindex]['uFE'] - ucdata[uindex]['uFE'])

def errOrderConc(vcdata,vindex):
    '''errOrder = |uncleFE(vaspindex) - uncleFE(uncleindex)|, where uncleFE is an ordered list containing only the vasp structures we have calculated.
    uncleindex is the index of the structure in uncleFE. vaspindex is the index of the structure in cdata, which must be ordered with the lowest energy at the top, 
    Need to reduce both lists to only those of the same concentration '''
   
    struct = vcdata[vindex]['struct']
    conc = vcdata[vindex]['conc']
    indicesKeep = []
    for istruct in range(len(vcdata)):
        conc2 = vcdata[istruct]['conc'] 
        if conc2 == conc:
            indicesKeep.append(istruct)
    vcdataRed =  zeros(len(indicesKeep),dtype = [('struct', int32),('vFE', float),('uFE', float)])  
    for istruct,ikeep in enumerate(indicesKeep):
        vcdataRed[istruct]['struct'] = vcdata[ikeep]['struct']
        vcdataRed[istruct]['uFE'] = vcdata[ikeep]['uFE']
        vcdataRed[istruct]['vFE'] = vcdata[ikeep]['vFE']
    vcdataRed = sort(vcdataRed,order = ['vFE'])
    ucdataRed = sort(vcdataRed,order = ['uFE'])
    ustructlist = ucdataRed['struct'].tolist()
    vstructlist = vcdataRed['struct'].tolist()
    uindex = ustructlist.index(struct)  
    vindex = vstructlist.index(struct)  
    return abs(ucdataRed[vindex]['uFE'] - ucdataRed[uindex]['uFE'])

def plotByPrior(atoms,minPrior,iteration):  
    ''''''
    lastDir = os.getcwd()
    structsDir = lastDir + '/structs'
    if not os.path.exists(structsDir):
        subprocess.call(['mkdir',structsDir])
    done = []
    #find the existing plots
    os.chdir(structsDir)
    for item in os.listdir(structsDir):
        if '.png' in item and item[0].isdigit(): 
            done.append(int(item.split('.')[0]))
    toPlot = []
    skipped = []

#    print done
    os.chdir(lastDir)
    for iatom, atom in enumerate(atoms):
        priorPath = lastDir + '/' + atom + '/priorities_{}.out'.format(iteration)
        lines = readfile(priorPath)
        for line in lines[1:]:
            prior = float(line.split()[1])
            struct = int(line.split()[0])
#            if prior >0: print struct,prior
            if os.path.isdir(lastDir + '/' + atom + '/' + str(struct)): # want to plot only calculated ones for now 
                if prior > minPrior and struct not in done:
                    toPlot.append(struct)
                elif prior > minPrior:
                    skipped.append(struct)
    toPlot = set(toPlot) #remove duplicates
    skipped = set(skipped)
    print 'Number of structures in priority range (priority > {}):  {}'.format(minPrior,len(toPlot)+len(skipped))
    print 'Number of structures already in structs folder:   {}'.format(len(done))
    print 'Number of structures left to plot:  {}'.format(len(toPlot))
#    print 'toPlot',toPlot
    plot_structs(toPlot)
    os.chdir(lastDir)

    
def plot_structs(structs):  
    '''If plot of structure is not in structsdir, then create it'''
    lastDir = os.getcwd()
    structsDir = lastDir + '/structs'
    if not os.path.exists(structsDir):
        subprocess.call(['mkdir',structsDir])
    done = []
    #find the existing plots
    os.chdir(structsDir)
    for item in os.listdir(structsDir):
        if '.png' in item and item[0].isdigit(): 
            done.append(int(item.split('.')[0]))
    #make new plots
    toPoscar = StructsToPoscar.structsToPoscar([],[])
    for i,struct in enumerate(structs):
        if struct not in done:
            #need the full POSCAR if we want to show empty C sites
            subprocess.call(['../needed_files/makestr.x','../enum/struct_enum.out',str(struct)])
            vfile = 'vasp.' + '0'*(6-len(str(struct))) + str(struct)
            toPoscar.convertOne(vfile)
            print i+1,struct
            plotSize = 23.5 #determines number of atoms in plot... Roughly N/4 hexagons across           
            plotter = PlotGraphene('.','POSCAR','-p','',plotSize)
            plotter.fillByVecs(int(plotSize*2.0))            
            plotter.addLines(int(plotSize*2.0))
            plotter.saveFigure()
            subprocess.call(['mv', 'POSCAR_plot.png', str(struct)+ '.png'])
            subprocess.call(['rm',vfile])             
    os.chdir(lastDir) 
    
class PlotGraphene:
    """ This plots either a POSCAR file or an output file (vasp.000xxx) from UNCLE.  When calling this class from
        the command line, call it as follows:
        
            python plotGraphene.py  /path/to/folder/  filename  [ -p | -u ]  [-z]
        
        The "-p" or "-u" tags tell us whether we are reading a POSCAR from VASP or an output file 
        from UNCLE.  The default tag (if none is supplied) is "-p".  The /path/to/folder/ is the 
        path to the folder containing the file we want to read and the filename is the name of the 
        actual file we want to read.  The -z tag enables differentiation between points above and 
        below the plane. """
    from comMethods import setAtomCounts    
    def __init__(self, poscarDir, poscarFile, kindTag, zTag, plotSize):

        if kindTag == "-u":
            self.uncle = True
        else:
            self.uncle = False
        
        self.poscarDir = poscarDir
        os.chdir(self.poscarDir)
        
        self.poscarFile = self.poscarDir + '/' + poscarFile
        
        if zTag == "-z" or zTag == "-Z":
            self.zFunc = True
        else:
            self.zFunc = False
        
        self.direct = False
        
        self.lattVec1 = []
        self.lattVec2 = []
        self.lattVec3 = []
        
        self.atomCounts = []
        self.xshift = 4.26257704
        self.yshift = 2.461
       
        self.origXlist = []
        self.origYlist = []
        self.origZlist = []
        
        self.xlist = []
        self.ylist = []
        self.zlist = []
         
        self.readPOSCAR()
         
        self.Ccirclelist = []
        self.Hcirclelist = []
        self.Mcirclelist = []
        self.plotSize = plotSize #determines number of atoms in plot... Roughly N/4 hexagons across   
        self.figure = None
#        self.ax = None
        self.initializeFigure()
    
        self.Hnum = 0
        self.Mnum = 0


        
    def readPOSCAR(self):
        
        poscarFile = open(self.poscarFile, 'r')
        poscarLines = poscarFile.readlines()
        
        # Extract the lattice vectors.
        if self.uncle:
            # We first need to find where the "z-vector" is.  It will contain a '********' in the 
            # first column.
            vecInds = []
            for i in [2,3,4]:
                zvec = False
                comps = poscarLines[i].strip().split()
                for j in xrange(3):
                    stringComp = comps[j]
                    if list(stringComp)[0] == '*':
                        zvec = True
                
                if not zvec:
                    vecInds.append(i)
                
            scrambledVec = []
            scrambledVec = poscarLines[vecInds[0]].strip().split()
            self.lattVec1 = [float(scrambledVec[1]), float(scrambledVec[2]), float(scrambledVec[0])]
            
            scrambledVec = poscarLines[vecInds[1]].strip().split()
            self.lattVec2 = [float(scrambledVec[1]), float(scrambledVec[2]), float(scrambledVec[0])]
        else:
            # We first need to find where the "z-vector" is.  We will know it is the z-vector if
            # its components are [0, 0, 15].  The Converter class makes this possible.
            vecInds = []
            for i in [2,3,4]:
                zvec = False
                comps = poscarLines[i].strip().split()
                if float(comps[0]) == 0.0 and float(comps[1]) == 0.0 and float(comps[2]) == 15.0:
                    zvec = True
                
                if not zvec:
                    vecInds.append(i)

            self.lattVec1 = poscarLines[vecInds[0]].strip().split()
            self.lattVec1 = [float(comp) for comp in self.lattVec1]
            
            self.lattVec2 = poscarLines[vecInds[1]].strip().split()
            self.lattVec2 = [float(comp) for comp in self.lattVec2]
        
        # Get the number of each type of atom.
        if self.uncle:
            nonCatomCounts = poscarLines[5].strip().split()
            nonCatomCounts = [int(count) for count in nonCatomCounts]
            
            self.atomCounts = [sum(nonCatomCounts), nonCatomCounts[0], nonCatomCounts[1]]
        else:
            self.setAtomCounts('.') # this is from comMethods...handles the pure cases right 
#            self.atomCounts = poscarLines[5].strip().split()
#            self.atomCounts = [int(count) for count in self.atomCounts]
        
        # Determine whether to compute atom positions in direct or cartesian coordinates.
        if poscarLines[6].strip().lower() == "d" or poscarLines[6].strip().lower() == "direct":
            self.direct = True
                
        positionLines = poscarLines[7:]
        
        if self.direct:
            if self.uncle:
                # If it is an UNCLE file, they don't explicitly state the C positions, so we run
                # this twice.  Once for the C positions and once for the H and M positions.
                self.getPositionsFromDirectCoordinates(positionLines, vecInds)
                self.getPositionsFromDirectCoordinates(positionLines, vecInds)
            else:
                self.getPositionsFromDirectCoordinates(positionLines)
        else:
            self.getPositionsFromCartesianCoordinates(positionLines)        
            
    def getPositionsFromDirectCoordinates(self, positionLines, vecInds):
        # Convert the vecInds to be zero-based.
        vecInds = [ind - 2 for ind in vecInds]
        
        zcoord = 1.0
        for line in positionLines:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[vecInds[0]] * self.lattVec1[0], position[vecInds[0]] * self.lattVec1[1], position[vecInds[0]] * self.lattVec1[2]]
            comp2 = [position[vecInds[1]] * self.lattVec2[0], position[vecInds[1]] * self.lattVec2[1], position[vecInds[1]] * self.lattVec2[2]]
            
            self.origXlist.append(float(comp1[0] + comp2[0]))
            self.xlist.append(float(comp1[0] + comp2[0]))
            self.origYlist.append(float(comp1[1] + comp2[1]))
            self.ylist.append(float(comp1[1] + comp2[1]))
            self.origZlist.append(zcoord)
            self.zlist.append(zcoord)
            
            # TODO:  Need to correct this buckling.
            if zcoord == 1.0:
                zcoord = -1.0
            else:
                zcoord = 1.0
        
    def getPositionsFromCartesianCoordinates(self, positionLines):
        for line in positionLines:
            positions = line.strip().split()
            self.origXlist.append(float(positions[0]))
            self.xlist.append(float(positions[0]))
            self.origYlist.append(float(positions[1]))
            self.ylist.append(float(positions[1]))
            self.origZlist.append(float(positions[2]))
            self.zlist.append(float(positions[2]))
    
    def initializeFigure(self):
        plt.clf() #don't carry over anything from previous figure!
        self.figure = plt.gcf()
        self.figure.gca().set_aspect('equal')
        plt.axis([0, self.plotSize, 0, self.plotSize]) 

    def periodicByVecs(self, vec1num, vec2num):
        self.Ccirclelist = []
        self.Hcirclelist = []
        self.Mcirclelist = []
        
        for i in range(0,sum(self.atomCounts)):
            self.xlist[i] = self.origXlist[i] + (vec1num * self.lattVec1[0]) + (vec2num * self.lattVec2[0])
            self.ylist[i] = self.origYlist[i] + (vec1num * self.lattVec1[1]) + (vec2num * self.lattVec2[1])
            
        numOfC = self.atomCounts[0]
        numOfH = self.atomCounts[1]
        numOfM = self.atomCounts[2]
        
        Cxlist = self.xlist[0:numOfC]
        Cylist = self.ylist[0:numOfC]
        Czlist = self.zlist[0:numOfC]

        Hxlist = self.xlist[numOfC:numOfC + numOfH]
        Hylist = self.ylist[numOfC:numOfC + numOfH]
        Hzlist = self.zlist[numOfC:numOfC + numOfH]

        Mxlist = self.xlist[numOfC + numOfH:sum(self.atomCounts)]
        Mylist = self.ylist[numOfC + numOfH:sum(self.atomCounts)]
        Mzlist = self.zlist[numOfC + numOfH:sum(self.atomCounts)]

        eps = 0.001

        for i in range(0,numOfC):
            if 0-eps <= Cxlist[i] <= self.plotSize + eps  and 0-eps <= Cylist[i] <= self.plotSize + eps: 
                self.Ccirclelist = self.Ccirclelist + [plt.Circle((Cxlist[i], Cylist[i]), .7, color='#2B65EC')]    
        for i in range(len(self.Ccirclelist)):
            self.figure.gca().add_artist(self.Ccirclelist[i])
        
        for i in range(0,numOfH):
            if 0-eps <= Hxlist[i] <= self.plotSize + eps  and 0-eps <= Hylist[i] <= self.plotSize + eps:
                if Hzlist[i] < 0 and self.zFunc:
                    self.Hcirclelist = self.Hcirclelist + [plt.Circle((Hxlist[i], Hylist[i]), .5, facecolor='#FFD801', edgecolor='black', lw = 3)]
                else:
                    self.Hcirclelist = self.Hcirclelist + [plt.Circle((Hxlist[i], Hylist[i]), .5, color='#FFD801')]      
        for i in range(len(self.Hcirclelist)):
            self.figure.gca().add_artist(self.Hcirclelist[i])
 
        for i in range(0,numOfM):
            if 0-eps <= Mxlist[i] <= self.plotSize + eps  and 0-eps <= Mylist[i] <= self.plotSize + eps:
                if Mzlist[i] < 0 and self.zFunc:
                    self.Mcirclelist = self.Mcirclelist + [plt.Circle((Mxlist[i], Mylist[i]), .5, facecolor='r', edgecolor='black', lw = 3)]
                else:
                    self.Mcirclelist = self.Mcirclelist + [plt.Circle((Mxlist[i], Mylist[i]), .5, color='r')]                   
        for i in range(len(self.Mcirclelist)):
            self.figure.gca().add_artist(self.Mcirclelist[i])
        
    def fillByVecs(self, num):
        if num == 0:
            self.periodicByVecs(0, 0)
        else:
            for i in xrange(-1, num + 1):
                for j in xrange(-1, num + 1):
                    self.periodicByVecs(i, j)
                    self.periodicByVecs(-i, j)
                    self.periodicByVecs(i, -j)
                    self.periodicByVecs(-i, -j)
    
    def addLines(self,reps):
        if self.lattVec1[0] == 0:
            slope1 = 'vertical'
        else:
            slope1 = self.lattVec1[1] / self.lattVec1[0]
        
        if self.lattVec2[0] == 0:
            slope2 = 'vertical'
        else:
            slope2 = self.lattVec2[1] / self.lattVec2[0]
#        print 'slope1, slope2',slope1, slope2
            
#        xpoints = array([-self.plotSize*2,-self.plotSize*2]) #make big enough that lines will cover larger structs

        for i in range(reps):
            self.shiftLineByVec1(i, slope2)
            self.shiftLineByVec1(-i, slope2)
            self.shiftLineByVec2(i, slope1)
            self.shiftLineByVec2(-i, slope1)

    def getIntercepts(self,slope,b):
        '''Find four intercepts of a line with the lines that define the square.  Only two of these
        will lie on the edges of the square'''
        a = self.plotSize
        points = []
        #x constant intercepts
#        print slope
#        print b,slope*a+b,-b/slope,(a-b)/slope
        if slope == 1.0 and b == 0: #diagonal line:
            points = [[0,0],[a,a]]
        elif  slope == -1.0 and b == a: #diagonal line:
            points = [[a,0],[0,a]]
        else:
            if 0 < b <= a:
                points.append([0,b])
            if 0 < slope*a+b <= a:
                points.append([a,slope*a+b])
            #y=constant intercepts
            if 0 <= -b/slope < a:
                points.append([-b/slope,0])
            if 0 <= (a-b)/slope < a:
                points.append([(a-b)/slope,a])         
        return points
            
    def shiftLineByVec1(self, n, slope):
        if slope == 0:
            dy = -n * self.lattVec1[1]
            if 0 <= dy <= self.plotSize:
                plt.axhline(y=dy,color='k')
        elif slope == 'vertical':
            dx = n * self.lattVec1[0]
            if 0 <= dx <= self.plotSize:
                plt.axvline(x=dx,color='k')
        else:
            b = n*(self.lattVec1[1] - slope*self.lattVec1[0])
            pts = self.getIntercepts(slope,b)
            if len(pts) == 2:
                xPoints = np.arange(pts[0][0],pts[1][0], 0.1*np.sign(pts[1][0]-pts[0][0])) 
                yPoints = slope * xPoints  + b   
                self.figure.gca().add_artist(Line2D([pts[0][0],pts[1][0]],[pts[0][1],pts[1][1]], color='k'))


    def shiftLineByVec2(self, n, slope):
        if slope == 0:
            dy = -n * self.lattVec2[1]
            if 0 <= dy <= self.plotSize:
                plt.axhline(y=dy,color='k')
        elif slope == 'vertical':
            dx = n * self.lattVec2[0]
            if 0 <= dx <= self.plotSize:
                plt.axvline(x=dx,color='k')
        else:
            b = n*(self.lattVec2[1] - slope*self.lattVec2[0])
            pts = self.getIntercepts(slope,b)
            if len(pts) == 2:
                xPoints = np.arange(pts[0][0],pts[1][0], 0.1*np.sign(pts[1][0]-pts[0][0])) 
                yPoints = slope * xPoints  + b
                self.figure.gca().add_artist(Line2D([pts[0][0],pts[1][0]],[pts[0][1],pts[1][1]], color='k'))
                
    def saveFigure(self):
        plt.plot()           
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.savefig(self.poscarFile + "_plot.png", bbox_inches='tight')
  
    