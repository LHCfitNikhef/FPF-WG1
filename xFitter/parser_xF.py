import numpy as np
import os
from distutils.dir_util import copy_tree

suffix = ".txt"
expIDs = ["FASERv_14",
          "FASERv_-14",
          "FASERv2_14",
          "FASERv2_-14"]
grids = ['nu_A_1-XSFPFCC.pineappl.lz4',   #neutrino
         'nub_A_1-XSFPFCC.pineappl.lz4',  #antineutrino
         'nu_A_1-XSFPFCC.pineappl.lz4',
         'nub_A_1-XSFPFCC.pineappl.lz4']
exptags = ["FASER$\\nu$ $\\nu_\\mu$", "FASER$\\nu$ $\\bar{\\nu}_\\mu$", 
           "FASER$\\nu$2 $\\nu_\\mu$","FASER$\\nu$2 $\\bar{\\nu}_\\mu$"]
expstrs = ["FASER#nu, #nu_{#mu}","FASER#nu, #bar{#nu}_{#mu}",
           "FASER#nu2, #nu_{#mu}","FASER#nu2, #bar{#nu}_{#mu}"]

subdirs=['datafiles/','lhc/','fpf/','neutrinoDIS/','pseudodata/']
datadir=''

#Produce subdirectories for tables
for subdir in subdirs:
    datadir+=subdir
    if not os.path.exists(datadir):
        os.mkdir(datadir)
if not os.path.exists(datadir+'grids/'):
    copy_tree('../theory/', datadir+'grids/')
#TODO to be most useful, should also untar the file.
#However, python tarfile is unsafe so you better do it manually

#Read results, write tables
for iexp,expID in enumerate(expIDs):
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
    for x in lines:
        xlo.append(float(x.split()[0]))
        xhi.append(float(x.split()[1]))
        xav.append(float(x.split()[2]))
        Q2lo.append(float(x.split()[3]))
        Q2hi.append(float(x.split()[4]))
        Q2av.append(float(x.split()[5]))
        sigma.append(10.**float(x.split()[9]))  #Convert from log10(sigma)
        stat.append(float(x.split()[11])/float(x.split()[10]))
        #Assume perfect detector, uncertainties only in energies and lepton angle
        #Also assume perfect flavor and sign recognition
        #Estimates based on J. Rojo's and T. Ariga's slides in FPF5 meeting
        #EnuUnc.append(0.3)   #Neutrino E unc ~ 30%
        ElUnc.append(0.3)    #Lepton E unc ~ 30%
        EhUnc.append(0.3)    #Hadronic FS E unc ~30-50% => let's be optimistic
        #theta unc. = 1 mrad, tan(theta) < 0.5, so O(theta unc)>1%. 
        #Should be bin-dependent, but let's take universal 10% for simplicity
        #Optimistically relative unc. > 0.001/0.46 and 
        thetaUnc.append(0.1)
    f.close()

    #Vary pseudodata
    sigmavar=[]
    fcorr=1.  #Data w/ correlated systematics would be more constraining
    fred=1.   #Reduction factor, expect improvements in detector accuracy etc
    uncor=[]
    for isig,sig in enumerate(sigma):
        #Gather unc. estimates into an uncorrelated unc. for PDF fits,
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
        sigmavar.append(sig) #DEBUG

#    #Test functionality to write predictions as a constant value table
#    f = open(datadir+"/th/"+expID+"-th.dat", "w")
#    for ix, x in enumerate(xlo):
#        f.write("{:16.6f}".format(Q2lo[ix]))
#        f.write("{:16.6f}".format(Q2hi[ix]))
#        f.write("{:16.6f}".format(xlo[ix]))
#        f.write("{:16.6e}".format(sigma[ix]))  
#        f.write("\n")
#    f.close()

    #Write xFitter table
    test=False
    expname="FASERv2" if "v2" in exptags[iexp] else "FASERv"
    f = open(datadir+expID+"-thexp.dat", "w")
    f.write("* "+exptags[iexp]+"CC DIS pseudodata\n")
    f.write("&Data\n")
    f.write("   Name = '"+exptags[iexp]+"'\n")
    f.write("   IndexDataset = "+str(137+iexp)+"\n") #Assign unique IDs to datasets
    f.write("   Reaction = 'CC nup'\n\n")
    f.write("   NData = "+str(len(xlo))+"\n")
    f.write("   NColumn = 6\n")
    f.write("   ColumnType = 'Flag',2*'Bin','Sigma',2*'Error'\n")
    f.write("   ColumnName = 'binFlag','x','Q2','Sigma','stat','uncor'\n")
    f.write("   TheoryType = 'expression\'\n")
    f.write("   TermType   = 'reaction'\n")
    if test:
        f.write("   TermName   = 'K'\n")
        f.write("   TermSource = 'KFactor'\n")
        f.write("   TermInfo   = 'FileName="+datadir+"th/"+expID+"-th.dat:FileColumn=4'\n")
        f.write("   TheorExpr  = 'K'\n\n")
    else:
        f.write("   TermName   = 'P'\n")
        f.write("   TermSource = 'PineAPPL'\n")
        f.write("   TermInfo   = 'GridName="+datadir+"grids/"+expID.replace('-','m')+'/grids-xsecs_A1/grids/'+gridname+"'\n")
        f.write("   TheorExpr  = 'P'\n\n")
    f.write("   Percent = 2*True\n")
    f.write("&End\n")
    f.write("&PlotDesc\n")
    f.write("   PlotN = 4\n")
    f.write("   PlotDefColumn = 'Q2'\n")
    f.write("   PlotVarColumn = 'x'\n")
    f.write("   PlotDefValue = 3., 5., 11., 110., 1100.\n")
    f.write("   PlotOptions(1)  = 'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 4 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("   PlotOptions(2)  = 'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 10 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("   PlotOptions(3)  = 'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 100 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("   PlotOptions(4)  = 'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2} > 1000 @Ymin:0.01@Xlog@Ylog'\n")
    f.write("&End\n")
    f.write("*binFlag          x              Q2")
    f.write("           Sigma            stat           uncor\n")
    for ix, x in enumerate(xlo):
        f.write("1  ")
        f.write("{:16.6f}".format(xav[ix]))
        f.write("{:16.6f}".format(Q2av[ix]))
        f.write("{:16.6e}".format(sigmavar[ix]))
        f.write("{:16.6f}".format(100.*stat[ix]))
        f.write("{:16.6f}".format(100.*uncor[ix]))
        ##Uncertainties, given in %
        ##However, would require correlation matrix, use uncor instead
        #f.write("{:16.6f}".format(100.*ElUnc[ix]))
        #f.write("{:16.6f}".format(100.*thetaUnc[ix]))
        #f.write("{:16.6f}".format(100.*EhUnc[ix]))
        f.write("\n")        
        
    f.close()
