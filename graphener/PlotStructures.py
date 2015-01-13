import sys, os, subprocess
import matplotlib.pyplot as plt
import numpy as np

def collate_structs(self,plotType,iteration):  
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
    
def plot_structs(structs):  
    '''If plot of structure is not in structsdir, then create it'''
    lastDir = os.getcwd()
    structsDir = lastDir + '/structs'
    if not os.path.exists(structsDir):
        subprocess.call(['mkdir',structsDir])
    done = []
#    PlotGraphene()
    #find the existing plots
    os.chdir(structsDir)
    for item in os.listdir(structsDir):
        if '.png' in item and item[0].isdigit(): 
            done.append(item)
    #make new plots
    for struct in structs:
        if struct not in done:
            #make psuedoPoscar
            subprocess.call(['../needed_files/makestr.x','../enum/struct_enum.out',str(struct)])
            vfile = 'vasp.' + '0'*(6-len(str(struct))) + str(struct)
            print vfile                      
            plotter = PlotGraphene('.',vfile,'-u','')
            plotter.fillByVecs(10)
            plotter.addLines()
            plotter.saveFigure()             
#            PlotGraphene('.',vfile,'-u','')
#    subprocess.call(['rm vasp.0*'],shell=True)
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
        
    def __init__(self, poscarDir, poscarFile, kindTag, zTag):
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
        
        self.figure = None
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
            self.atomCounts = poscarLines[5].strip().split()
            self.atomCounts = [int(count) for count in self.atomCounts]
        
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
        self.figure = plt.gcf()
        self.figure.gca().set_aspect('equal')
        plt.axis([0, 15, 0, 15])
       
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
    
        for i in range(0,numOfC):
            self.Ccirclelist = self.Ccirclelist + [plt.Circle((Cxlist[i], Cylist[i]), .7, color='#2B65EC')]
    
        for i in range(0,numOfC):
            self.figure.gca().add_artist(self.Ccirclelist[i])
            
        for i in range(0,numOfH):
            if Hzlist[i] < 0 and self.zFunc:
                self.Hcirclelist = self.Hcirclelist + [plt.Circle((Hxlist[i], Hylist[i]), .5, facecolor='#FFD801', edgecolor='black', lw = 3)]
                self.Hnum += 1
            else:
                self.Hcirclelist = self.Hcirclelist + [plt.Circle((Hxlist[i], Hylist[i]), .5, color='#FFD801')]
                self.Hnum += 1
        
        for i in range(0,numOfH):
            self.figure.gca().add_artist(self.Hcirclelist[i])

        for i in range(0,numOfM):
            if Mzlist[i] < 0 and self.zFunc:
                self.Mcirclelist = self.Mcirclelist + [plt.Circle((Mxlist[i], Mylist[i]), .5, facecolor='r', edgecolor='black', lw = 3)]
                self.Mnum += 1
            else:
                self.Mcirclelist = self.Mcirclelist + [plt.Circle((Mxlist[i], Mylist[i]), .5, color='r')]
                self.Mnum += 1
                
        for i in range(0,numOfM):
            self.figure.gca().add_artist(self.Mcirclelist[i])
    
    def fillByVecs(self, num):
        if num == 0:
            self.periodicByVecs(0, 0)
        else:
            for i in xrange(0, num + 1):
                for j in xrange(0, num + 1):
                    self.periodicByVecs(i, j)
                    self.periodicByVecs(-i, j)
                    self.periodicByVecs(i, -j)
                    self.periodicByVecs(-i, -j)
    
    def addLines(self):
        if self.lattVec1[0] == 0:
            slope1 = 'vertical'
        else:
            slope1 = self.lattVec1[1] / self.lattVec1[0]
        
        if self.lattVec2[0] == 0:
            slope2 = 'vertical'
        else:
            slope2 = self.lattVec2[1] / self.lattVec2[0]
        
        xPoints = np.arange(-20, 20, .1)
        
        for i in range(10):
            self.shiftLineByVec1(i, xPoints, slope2)
            self.shiftLineByVec1(-i, xPoints, slope2)
            self.shiftLineByVec2(i, xPoints, slope1)
            self.shiftLineByVec2(-i, xPoints, slope1)
    
    def shiftLineByVec1(self, ntimes, xPoints, slope):
        if slope == 0:
            plt.plot(xPoints, (slope * xPoints) + (-ntimes * self.lattVec1[1]), color='k')
        elif slope == 'vertical':
            plt.plot(((ntimes * self.lattVec1[0]), (ntimes * self.lattVec1[0])), (0, 15), color='k')
        else:
            plt.plot(xPoints, 
                    (slope * xPoints) + (slope * (-ntimes * self.lattVec1[0])) + (ntimes * self.lattVec1[1]), color='k')
    
    def shiftLineByVec2(self, ntimes, xPoints, slope):
        if slope == 0:
            plt.plot(xPoints, (slope * xPoints) + (-ntimes * self.lattVec2[1]), color='k')
        elif slope == 'vertical':
            plt.plot(((ntimes * self.lattVec2[0]), (ntimes * self.lattVec2[0])), (0, 15), color='k')
        else:
            plt.plot(xPoints, 
                     (slope * xPoints) + (slope * (-ntimes * self.lattVec2[0])) + (ntimes * self.lattVec2[1]), color='k')
        
    def saveFigure(self):
        plt.plot()
        #plt.gca().axes.get_xaxis().set_visible(False)
        #plt.gca().axes.get_yaxis().set_visible(False)
        plt.savefig(self.poscarFile + "_plot.png", bbox_inches='tight')
  
    