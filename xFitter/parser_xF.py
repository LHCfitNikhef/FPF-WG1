import numpy as np
import os

suffix = ".txt"
expIDs = ["FASERv_14","FASERv_-14", "FASERv2_14","FASERv2_-14"]
exptags = ["FASER$\\nu$ $\\nu_\\mu$", "FASER$\\nu$ $\\bar{\\nu}_\\mu$", 
           "FASER$\\nu$2 $\\nu_\\mu$","FASER$\\nu$2 $\\bar{\\nu}_\\mu$"]
expstrs = ["FASER#nu, #nu_{#mu}","FASER#nu, #bar{#nu}_{#mu}",
           "FASER#nu2, #nu_{#mu}","FASER#nu2, #bar{#nu}_{#mu}"]

#FIXME Const. toy estimates of unc. scale, to get chi2/N ~ 1
#TODO  make these bin-by-bin & experimentally motivated!
expUncs = [0.07, 0.15, 0.015, 0.025]

subdirs=["datafiles/","lhc/","fpf/","neutrinoDIS/","pseudodata/"]
datadir=""
for subdir in subdirs:
    datadir+=subdir
    if not os.path.exists(datadir):
        os.mkdir(datadir)
#FIXME for now, theory grids replaced by constant table directly from Max's results
if not os.path.exists(datadir+"th/"):
    os.mkdir(datadir+"th/")

for iexp,expID in enumerate(expIDs):
    infile = "../results/binned_events_" + expID + suffix
    #Read the predictions / pseudodata template
    f = open(infile, "r")
    lines=f.readlines()
    lines.pop(0)    #Remove 2 lines of header
    lines.pop(0)
    xlo=[]
    Q2lo=[]
    Q2hi=[]
    sigma=[]
    stat=[]
    for x in lines:
        xlo.append(float(x.split()[0]))  #Split at whitespace by default
        Q2lo.append(float(x.split()[3]))
        Q2hi.append(float(x.split()[4]))
        sigma.append(10**float(x.split()[9]))  #Convert from log10(sigma)?
        stat.append(100.*float(x.split()[11])/float(x.split()[10]))
    f.close()

    #Vary pseudodata
    sigmavar=[]
    for isig,sig in enumerate(sigma):
        sigmavar.append(abs(sig*(1. + expUncs[iexp]*np.random.normal())))
    
    #Write predictions as a constant value table
    f = open(datadir+"/th/"+expID+"-th.dat", "w")
    for ix, x in enumerate(xlo):
        f.write("{:16.6f}".format(Q2lo[ix]))
        f.write("{:16.6f}".format(Q2hi[ix]))
        f.write("{:16.6f}".format(xlo[ix]))
        f.write("{:16.6e}".format(sigma[ix]))  
        f.write("\n")
    f.close()

    #Write xFitter table
    f = open(datadir+expID+"-thexp.dat", "w")
    f.write("* "+exptags[iexp]+"CC DIS pseudodata\n")
    f.write("&Data\n")
    f.write("   Name = \'"+exptags[iexp]+"\'\n")
    f.write("   IndexDataset = "+str(137+iexp)+"\n") #Assign unique IDs to datasets
    f.write("   Reaction = \'CC nup\'\n\n")
    f.write("   NData = "+str(len(xlo))+"\n")
    f.write("   NColumn = 6\n")
    f.write("   ColumnType = 'Flag',3*\'Bin\',\'Sigma\',\'Error\'\n")
    f.write("   ColumnName = 'binFlag',\'x\',\'Q2min\',\'Q2max\',\'Sigma\',\'stat\'\n\n")
    f.write("   TheoryType = \'expression\'\n")
    f.write("   TermName   = \'K\'\n")
    f.write("   TermType   = \'reaction\'\n")
    f.write("   TermSource = \'KFactor\'\n")
    #FIXME below must eventually point to actual prediction grids
    f.write("   TermInfo   = \'FileName="+datadir+"th/"+expID+"-th.dat:FileColumn=4\'\n")
    f.write("   TheorExpr  = \'K\'\n\n")
    f.write("   Percent = 1*True\n")
    f.write("&End\n")
    f.write("&PlotDesc\n")
    f.write("   PlotN = 4\n")
    f.write("   PlotDefColumn = \'Q2min\'\n")
    f.write("   PlotVarColumn = \'x\'\n")
    f.write("   PlotDefValue = 4., 10., 100., 1000.\n")
    f.write("   PlotOptions(1)  = \'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2}_{min} = 4 @Ymin:0.01@Xlog@Ylog\'\n")
    f.write("   PlotOptions(2)  = \'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2}_{min} = 10 @Ymin:0.01@Xlog@Ylog\'\n")
    f.write("   PlotOptions(3)  = \'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2}_{min} = 100 @Ymin:0.01@Xlog@Ylog\'\n")
    f.write("   PlotOptions(4)  = \'Experiment:"+expstrs[iexp]+"@ExtraLabel:#nu p (CC)")
    f.write(  " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}")
    f.write(  " @Title:Q^{2}_{min} = 1000 @Ymin:0.01@Xlog@Ylog\'\n")
    f.write("&End\n")
    f.write("*binFlag        xlo           Q2min           Q2max")
    f.write("           Sigma            stat\n")
    for ix, x in enumerate(xlo):
        f.write("1  ")
        f.write("{:16.6f}".format(xlo[ix]))
        f.write("{:16.6f}".format(Q2lo[ix]))
        f.write("{:16.6f}".format(Q2hi[ix]))
        f.write("{:16.6e}".format(sigmavar[ix]))
        f.write("{:16.6f}".format(stat[ix]))
        f.write("\n")
    f.close()
