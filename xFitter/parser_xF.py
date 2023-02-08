import numpy as np
import os
from distutils.dir_util import copy_tree

PDF="PDF4LHC21"
suffix = ".txt"
expIDs = ["FASERv_14",
          "FASERv_-14",
          "FASERv2_14",
          "FASERv2_-14"]
grids = ['nu_A_1-XSFPFCC.pineappl.lz4',   #neutrino
         'nub_A_1-XSFPFCC.pineappl.lz4',  #antineutrino
         'nu_A_1-XSFPFCC.pineappl.lz4',
         'nub_A_1-XSFPFCC.pineappl.lz4']

subdirs=['datafiles/','lhc/','fpf/','neutrinoDIS/','pseudodata/']
mdir=''  #To contain master data directory path, subdir "grids" & "prel"

#Produce subdirectories for tables
for subdir in subdirs:
    mdir+=subdir
    if not os.path.exists(mdir):
        os.mkdir(mdir)
if not os.path.exists(mdir+'/'+PDF):
    os.mkdir(mdir+'/'+PDF)
if not os.path.exists(mdir+'/'+PDF+'/prel/'):
    os.mkdir(mdir+'/'+PDF+'/prel/')
if not os.path.exists(mdir+'grids/'):
    copy_tree('../theory/', mdir+'grids/')
#TODO to be most useful, should also untar the file.
#However, python tarfile is unsafe so you better do it manually

#Param  sDir     Subdirectory under mdir where to write the tables
#       expID    The filename part with experiment name and particle  
#       index    Dataset must be assigned a unique index
#       x,Q2     Vec. of x and Q^2 values for each data point
#       xsec     Vec. of cross-sections
#       stat     Vec. of statistical uncertainties for each data point
#       uncor    ^-||- uncorrelated uncertainties
#       sStr     Names of sources of systematic uncertainty
#       syst     Matrix of syst. unc. values
def xFtableWriter(sdir,expID,index,x,Q2,xsec,stat,uncor,sStr,syst):
    #LaTeX and ROOT style strings for experiment identification
    exptag = expID.replace("v",    "$\\nu$"            )\
                  .replace("_14", " $\\nu_\\mu$"       )\
                  .replace("_-14"," $\\bar{\\nu}_\\mu$")
    expstr = exptag.replace("$",""  )\
                   .replace("\\","#" )#\
                   #.replace(" ", ", ")

    #Lightweight test mode: write predictions as a constant value table
    test=False

    expname="FASERv2" if "v2" in expID else "FASERv"
    f = open(mdir+sdir+expID+"-thexp.dat", "w")
    f.write("* "+exptag+" CC DIS pseudodata\n")
    f.write("&Data\n")
    f.write("   Name = '"+exptag+"'\n")
    f.write("   IndexDataset = "+str(index)+"\n") #Give set a unique ID
    f.write("   Reaction = 'CC nup'\n\n")
    f.write("   NData = "+str(len(x))+"\n")
    NCol=1+2+1+2+len(sStr)
    f.write("   NColumn = "+str(NCol)+"\n")
    f.write("   ColumnType = 'Flag',2*'Bin','Sigma',"+str(2+len(sStr))\
                                                     +"*'Error'\n")
    f.write("   ColumnName = 'binFlag','x','Q2','Sigma','stat','uncor'")
    for s in sStr:
		f.write(",'"+s+"'")
    f.write(    "\n")
    f.write("   TheoryType = 'expression\'\n")     
    f.write("   TermType   = 'reaction'\n")
    if test:
        f.write("   TermName   = 'K'\n")
        f.write("   TermSource = 'KFactor'\n")
        f.write("   TermInfo   = 'FileName="+mdir+"th/"+expID\
                                            +"-th.dat:FileColumn=4'\n")
        f.write("   TheorExpr  = 'K'\n\n")
    else:
        f.write("   TermName   = 'P'\n")
        f.write("   TermSource = 'PineAPPL'\n")
        f.write("   TermInfo   = 'GridName="+mdir+"grids/"\
                                            +expID.replace('-','m')\
                                            +'/grids-xsecs_A1/grids/'\
                                            +gridname+"'\n")
        f.write("   TheorExpr  = 'P'\n\n")
    f.write("   Percent = 2*True\n")
    f.write("&End\n")
    f.write("&PlotDesc\n")
    f.write("   PlotN = 4\n")
    f.write("   PlotDefColumn = 'Q2'\n")
    f.write("   PlotVarColumn = 'x'\n")
    f.write("   PlotDefValue = 3., 5., 11., 110., 1100.\n")
    f.write("   PlotOptions(1)  = 'Experiment:"\
                                   +expstr+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 4 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("   PlotOptions(2)  = 'Experiment:"\
                                   +expstr+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 10 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("   PlotOptions(3)  = 'Experiment:"\
                                   +expstr+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 100 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("   PlotOptions(4)  = 'Experiment:"\
                                   +expstr+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 1000 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("&End\n")
    f.write("*binFlag          x              Q2")
    f.write("           Sigma            stat           uncor")
    for su in sStr:        
		f.write("%16s" % su)
    f.write("\n")
    for ix in range(len(x)):
        f.write("1  ")
        f.write("{:16.6f}".format(x[ix]))
        f.write("{:16.6f}".format(Q2[ix]))
        f.write("{:16.6e}".format(xsec[ix]))
        f.write("{:16.6f}".format(100.*stat[ix]))
        f.write("{:16.6f}".format(100.*uncor[ix]))
        for isu,su in enumerate(sStr):  #Systematic uncertainties
            f.write("{:16.6f}".format(100.*syst[isu][ix]))			
        #TODO rm below if obsolete
        ##Uncertainties, given in %
        ##However, would require correlation matrix, use uncor instead
        #f.write("{:16.6f}".format(100.*ElUnc[ix]))
        #f.write("{:16.6f}".format(100.*thetaUnc[ix]))
        #f.write("{:16.6f}".format(100.*EhUnc[ix]))
        f.write("\n")        
        
    f.close()
# END xFtableWriter

# BEGIN "main"
for iexp,expID in enumerate(expIDs):

	#True:  write tables for preliminary xFitter run
	#False: write final tables to be used as pseudodata in fits.
    writePrel = False

    gridname = grids[iexp]
    infile = "../results/binned_events_" + expID + suffix
    #Read the predictions / pseudodata template
    f = open(infile, "r")
    lines=f.readlines()
    lines.pop(0)    #Remove 2 lines of header
    lines.pop(0)
    xlo,xhi,xav=[],[],[]
    Q2lo,Q2hi,Q2av=[],[],[]
    sigma=[]
    stat,EnuUnc,ElUnc,EhUnc,thetaUnc=[],[],[],[],[]
    for l in lines:
        xlo.append(       float(l.split()[0]))
        xhi.append(       float(l.split()[1]))
        xav.append(       float(l.split()[2]))
        Q2lo.append(      float(l.split()[3]))
        Q2hi.append(      float(l.split()[4]))
        Q2av.append(      float(l.split()[5]))
        sigma.append(10.**float(l.split()[9]))  #Cnvrt from log10(sigma)
        stat.append(      float(l.split()[11])/float(l.split()[10]))
        #Assume perfect detector, uncertainties only in energies and lepton angle
        #Also assume perfect flavor and sign recognition
        #Estimates based on J. Rojo's and T. Ariga's slides in FPF5 meeting
        #EnuUnc.append(0.3)   #Neutrino E unc ~ 30%
        ElUnc.append(0.3)    #Lepton E unc ~ 30%
        EhUnc.append(0.3)    #Hadronic FS E unc ~30-50% => let's be optimistic
        #theta unc. = 1 mrad, tan(theta) < 0.5, so O(theta unc)>1%. 
        #Should be bin-dependent, but take universal 10% for simplicity
        #Optimistically relative unc. > 0.001/0.46 and 
        thetaUnc.append(0.1)
    f.close()

    #Vary pseudodata
    sigmavar=[]
    fcorr=1.  #Data w/ correlated systematics would be more constraining
    fred=1.   #Reduction factor: improvements in detector accuracy etc
    uncor=[]
    for isig,sig in enumerate(sigma):
        #Gather unc estimates into an uncorrelated unc for fits,
        #estimating correlations would be out of scope for now
        expUnc  = (fcorr*fred*ElUnc[   isig])**2
        expUnc += (fcorr*fred*EhUnc[   isig])**2
        expUnc += (fcorr*fred*thetaUnc[isig])**2
        #expUnc += (fcorr*fred*EnuUnc[  isig])**2 #FIXME check if this should be there
        uncor.append(np.sqrt(expUnc))

        #Form pseudodata: vary M.Fieg's results by the estimated exp. unc.
        expUnc += stat[isig]**2  #Also stat. unc. factors into the variations
        expUnc = np.sqrt(expUnc)
#        sigmavar.append(sig*(1. + expUnc*np.random.normal()))
        sigmavar.append(sig) #FIXME DEBUG

    #Write xFitter tables for prel run
    if writePrel:
        xFtableWriter(PDF+"/prel/",expID,137+iexp, \
                      xav,Q2av,sigmavar,stat,uncor,[],[])
    else:
        #Read xFitter output from preliminary run
        f = open("PDF_profiling/"+PDF+"/prel/"+expID.replace('-','m')\
                                 +"/output/fittedresults.txt","r")
        lines=f.readlines()
        #Skip header, data starts at line 9
        for i in range(8):
            lines.pop(0)
        xsec=[]
        for l in lines:
            xsec.append(float(l.split()[6]))
        f.close()
        xsecvar = xsec  #FIXME vary to form pseudodata
        xFtableWriter(PDF+"/",expID,137+iexp, \
                      xav,Q2av,xsecvar,stat,uncor,[],[]) #TODO syst unc matrix
# END "main"
