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
#       binF     Flags (int 1 or 0) if bin used in fit or not
#       x,Q2     Vec. of x and Q^2 values for each data point
#       xsec     Vec. of cross-sections
#       stat     Vec. of statistical uncertainties for each data point
#       uncor    Vec. of uncorrelated uncertainties
#       sStr     Names of correlated systematic uncertainty sources
#       syst     Matrix of corr. syst. unc. values
def xFtableWriter(sdir,expID,index,binF,x,Q2,xsec,stat,uncor,sStr,syst):
    #LaTeX and ROOT style strings for experiment identification
    exptag = expID.replace("v",    "$\\nu$"              )\
                  .replace("_14", " $\\nu_{\\mu}$"       )\
                  .replace("_-14"," $\\bar{\\nu}_{\\mu}$")
    expstr = exptag.replace("$",""  )\
                   .replace("\\","#" )

    #Check if uncor and stat uncs given
    useUncor = len(uncor) > 0
    useStat  = len(stat)  > 0

    #Lightweight test mode: write predictions as a constant value table
    test=False

    if not os.path.exists(mdir+sdir):
        os.mkdir(mdir+sdir)

    expname="FASERv2" if "v2" in expID else "FASERv"
    f = open(mdir+sdir+expID+"-thexp.dat", "w")
    f.write("* "+exptag+" CC DIS pseudodata\n")
    f.write("&Data\n")
    f.write("   Name = '"+exptag+"'\n")
    f.write("   IndexDataset = "+str(index)+"\n") #Give set a unique ID
    f.write("   Reaction = 'CC nup'\n\n")
    f.write("   NData = "+str(len(x))+"\n")
    Nerr = len(sStr)
    if useStat:
        Nerr = Nerr+1
    if useUncor:
        Nerr = Nerr+1
    Ncol=1+2+1+Nerr
    f.write("   NColumn = "+str(Ncol)+"\n")
    f.write("   ColumnType = 'Flag',2*'Bin','Sigma'")
    if Nerr>0:
        f.write(","+str(Nerr)+"*'Error'")
    f.write("\n")
    f.write("   ColumnName = 'binFlag','x','Q2','Sigma'")
    if useStat:
        f.write(",'stat'")
    if useUncor:
        f.write(",'uncor'")
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
                                            +expID.replace('-','m')+'/'\
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
    f.write("*binFlag          x              Q2           Sigma")
    if useStat:
        f.write("            stat")
    if useUncor:
        f.write("           uncor")
    for su in sStr:        
        f.write("%16s" % su)
    f.write("\n")
    ixsec=0
    for ix in range(len(x)):
        f.write("  "+str(binF[ix]))
        f.write("{:16.6f}".format(x[ix]))
        f.write("{:16.6f}".format(Q2[ix]))
        if binF[ix]==1:
            f.write("{:16.6e}".format(xsec[ixsec]))
            ixsec += 1
        else:
            f.write("        0.000000")
        if useStat:
            f.write("{:16.6f}".format(100.*stat[ix])) #Turn into %
        if useUncor:
            f.write("{:16.6f}".format(100.*uncor[ix])) #Turn into %
        for isu,su in enumerate(sStr):  #Systematic uncertainties
            f.write("{:16.6f}".format(100.*syst[isu][ix]))          
        f.write("\n")                
    f.close()
# END xFtableWriter

# BEGIN "main"
for iexp,expID in enumerate(expIDs):

    #Flags for including uncertainties
    useStat  = True
    useUncor = True
    #Syst. unc are included if a non-empty array is provided

    #True:  write tables for preliminary xFitter run
    #False: write final tables to be used as pseudodata in fits.
    writePrel = False

    #Init
    binF=[]
    xlo,xhi,xav=[],[],[]
    Q2lo,Q2hi,Q2av=[],[],[]
    sigma,stat,uncor=[],[],[]
 
    gridname = grids[iexp]
    infile = "../results/binned_sysevents_" + expID + suffix
    #Read the predictions / pseudodata template
    f = open(infile, "r")
    lines=f.readlines()
    lines.pop(0)    #Remove 2 lines of header
    lines.pop(0)
    for l in lines:
        lsplit = l.split()
        binF.append(1)
        xlo.append(float(lsplit[0]))
        xhi.append(float(lsplit[1]))
        xavstr=lsplit[2]
        if xavstr=="nan":
            xav.append(0.5*(xlo[-1]+xhi[-1]))
        else:
            xav.append(float(xavstr))
        Q2lo.append(float(lsplit[3]))
        Q2hi.append(float(lsplit[4]))
        Q2avstr = lsplit[5]
        if Q2avstr=="nan":
            Q2av.append(0.5*(Q2lo[-1]+Q2hi[-1]))
        else:
            Q2av.append(float(Q2avstr))
        sigstr = lsplit[9]
        if sigstr=="nan":
            sigma.append(0.)
            binF[-1]=0
        else:
            sigma.append(10.**float(sigstr))  #Convert from log10(sigma)
        Nevtstr = lsplit[10]
        if useStat:
            if Nevtstr=="0.00000":
                stat.append(0.)             
                binF[-1]=0
            else:
                stat.append(float(lsplit[11])/float(Nevtstr))
        if useUncor:
            if Nevtstr=="0.00000":
                uncor.append(0.)             
                binF[-1]=0
            else:
                uncor.append(float(lsplit[12])/float(Nevtstr))
    f.close()

    uncorALT=np.multiply(0.5,uncor)  #Optimistic scen., smaller unc.
    xsecvar=[]
    xsecvarALT=[]
    fcorr=1.  #Data w/ correlated systematics would be more constraining
    fred=1.   #Reduction factor: improvements in detector accuracy etc

    #Write xFitter tables for prel run
    if writePrel:
        xFtableWriter(PDF+"/prel/",expID,137+iexp, \
                      binF,xav,Q2av,sigma,stat,uncor,[],[])

    #Write xFitter tables for final run
    else:
 
        #Read xFitter output from preliminary run
        f = open("PDF_profiling/"+PDF+"/prel/"+expID.replace('-','m')\
                                 +"/output/fittedresults.txt","r")
        lines=f.readlines()
        for i in range(8):  #Skip header, data starts at line 9
            lines.pop(0)
        xsec=[]
        for l in lines:
            xsec.append(float(l.split()[6]))
        f.close()
 
        #Form pseudodata: vary results by the estimated exp. unc.
        for ixs,xs in enumerate(xsec):
            expUnc    = 0
            expUncALT = 0
            if useStat:
                expUnc    += stat[ixs]**2
                expUncALT += stat[ixs]**2
            if useUncor:
                expUnc    +=    uncor[ixs]**2
                expUncALT += uncorALT[ixs]**2
            expUnc    = np.sqrt(expUnc)
            expUncALT = np.sqrt(expUncALT)
            xsecvar.append(   xs*(1. +    expUnc*np.random.normal()))
            xsecvarALT.append(xs*(1. + expUncALT*np.random.normal()))
        xFtableWriter(PDF+"/",expID,137+iexp, \
                      binF,xav,Q2av,xsecvar,stat,uncor,[],[])
        xFtableWriter(PDF+"/uncorHalf/",expID,137+iexp, \
                      binF,xav,Q2av,xsecvarALT,stat,uncorALT,[],[])
# END "main"
