#!/usr/bin/env python
# coding: utf-8

import os,sys
import lhapdf
import numpy as np
import matplotlib.pyplot as py
from matplotlib import gridspec
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
from pylab import *
import scipy
from scipy.integrate import dblquad

print("\n *********************************************************")
print("      Impact of FPF data                               ")
print(" ***********************************************************\n")

q = 100 # GeV
stringQ=r'$Q=100~{\rm GeV}$'
filelabel='q100gev-ratios'
nx = 1000
xmin = 1e-3
xmax= 0.95
# Reduce verbosity of LHAPDF
lhapdf.setVerbosity(0)
# max number of replicas
nrepmax=1000
# number of flavours to be plotted
nfl=6
# Set x grid
X = np.logspace(log(xmin),log(xmax),nx)
# number of pdf sets
nset=2

nrep=np.zeros(nset, dtype='int')
nrep_max = 1000

pdfset=["230309-ern-002","230718-jcmwfaser-001"]

for iset in range(nset):

    p=lhapdf.getPDFSet(pdfset[iset])
    nrep[iset]=int(p.get_entry("NumMembers"))-1
    
    if(iset==0):
        fit1 = np.zeros((nrep[iset],nfl,nx))
    if(iset==1):
        fit2 = np.zeros((nrep[iset],nfl,nx))
    if(iset==2):
        fit3 = np.zeros((nrep[iset],nfl,nx))

    print(pdfset[iset])
    print("nrep = ",nrep[iset])

    # Run over replicas
    for i in range(1,nrep[iset]+1):
        p=lhapdf.mkPDF(pdfset[iset],i)
        lhapdf.setVerbosity(0)

        #print("irep = ",i)
        
        # Run over x arrat
        for k in range(nx):
            
            x = X[k]
            q2 = pow(q,2.0)

            # run over flavours
            for ifl in range(nfl):

                if(ifl==0):
                    # uV
                    if(iset==0):
                        fit1[i-1][ifl][k] = p.xfxQ(2,x,q) - p.xfxQ(-2,x,q)
                    if(iset==1):
                        fit2[i-1][ifl][k] = p.xfxQ(2,x,q) - p.xfxQ(-2,x,q)

                elif(ifl==1):
                    # dV
                    if(iset==0):
                        fit1[i-1][ifl][k] = p.xfxQ(1,x,q) - p.xfxQ(-1,x,q)
                    if(iset==1):
                        fit2[i-1][ifl][k] = p.xfxQ(1,x,q) - p.xfxQ(-1,x,q)

                elif(ifl==2):
                    # gluon
                    if(iset==0):
                        fit1[i-1][ifl][k] = p.xfxQ(0,x,q)
                    if(iset==1):
                        fit2[i-1][ifl][k] = p.xfxQ(0,x,q) 

                elif(ifl==3):
                    # singlet
                    if(iset==0):
                        fit1[i-1][ifl][k] = p.xfxQ(2,x,q) + p.xfxQ(-2,x,q) +\
                            p.xfxQ(1,x,q) + p.xfxQ(-1,x,q) +\
                            p.xfxQ(3,x,q) + p.xfxQ(-3,x,q)
                    if(iset==1):
                        fit2[i-1][ifl][k] = p.xfxQ(2,x,q) + p.xfxQ(-2,x,q) +\
                            p.xfxQ(1,x,q) + p.xfxQ(-1,x,q) +\
                            p.xfxQ(3,x,q) + p.xfxQ(-3,x,q)

                elif(ifl==4):
                    # s^+
                    if(iset==0):
                        fit1[i-1][ifl][k] = p.xfxQ(3,x,q) + p.xfxQ(-3,x,q) 
                    if(iset==1):
                        fit2[i-1][ifl][k] = p.xfxQ(3,x,q) + p.xfxQ(-3,x,q) 

                elif(ifl==5):
                    # c^+
                    if(iset==0):
                        fit1[i-1][ifl][k] = p.xfxQ(4,x,q) + p.xfxQ(-4,x,q) 
                    if(iset==1):
                        fit2[i-1][ifl][k] = p.xfxQ(4,x,q) + p.xfxQ(-4,x,q) 
                                 
                # end run over sets 
print("PDF arrays succesfully filled")

#---------------------------------------------------------------------
# Compute central values and uncertainties
#---------------------------------------------------------------------

p1_high = np.nanpercentile(fit1,84,axis=0)
p1_low = np.nanpercentile(fit1,16,axis=0)
p1_mid = ( p1_high + p1_low )/2.
p1_error = ( p1_high - p1_low )/2.

p2_high = np.nanpercentile(fit2,84,axis=0)
p2_low = np.nanpercentile(fit2,16,axis=0)
p2_mid = ( p2_high + p2_low )/2.
p2_error = ( p2_high - p2_low )/2.

#p3_high = np.nanpercentile(fit3,84,axis=0)
#p3_low = np.nanpercentile(fit3,16,axis=0)
#p3_mid = ( p3_high + p3_low )/2.
#p3_error = ( p3_high - p3_low )/2.

#---------------------------------------------------------------------
# Plot absolute SFs
#---------------------------------------------------------------------

py.clf()
ncols,nrows=3,2
py.figure(figsize=(ncols*5,nrows*3.5))
gs = gridspec.GridSpec(nrows,ncols)
rescolors = py.rcParams['axes.prop_cycle'].by_key()['color']

# pdflabels
labelpdf=[r"$xu_V(x,Q)/xu_V^{(\rm ref)}(x,Q)$",\
          r"$xd_V(x,Q)/xd_V^{(\rm ref)}(x,Q)$",\
          r"$xg(x,Q)/xg^{(\rm ref)}(x,Q)$",\
          r"$x\Sigma(x,Q)/x\Sigma^{(\rm ref)}(x,Q)$",\
          r"$xs^+(x,Q)/xs^+_{(\rm ref)}(x,Q)$",\
          r"$xc^+(x,Q)/xc^+_{(\rm ref)}(x,Q)$"]

icount=0
for ifl in range(nfl):

    ax = py.subplot(gs[icount])

    norm = p1_mid[ifl]
    
    p1=ax.plot(X,p1_mid[ifl]/norm,ls="dotted",color=rescolors[0])
    ax.fill_between(X,p1_high[ifl]/norm,p1_low[ifl]/norm,color=rescolors[0],alpha=0.2)
    p2=ax.fill(np.NaN,np.NaN,color=rescolors[0],alpha=0.2)

    norm = p2_mid[ifl]

    p3=ax.plot(X,p2_mid[ifl]/norm,ls="dashed",color=rescolors[1])
    ax.fill_between(X,p2_high[ifl]/norm,p2_low[ifl]/norm,color=rescolors[1],alpha=0.2)
    p4=ax.fill(np.NaN,np.NaN,color=rescolors[1],alpha=0.2)

#    p5=ax.plot(X,p3_mid[ifl],ls="solid")
#    ax.fill_between(X,p3_high[ifl],p3_low[ifl],color=rescolors[2],alpha=0.2)
#    p6=ax.fill(np.NaN,np.NaN,color=rescolors[2],alpha=0.2)


    #ax.set_xscale('linear')
    ax.set_xscale('log')    
    ax.set_xlim(0.004,0.7)
    ax.set_yscale('linear')
    ax.set_ylim(0.8,1.25)
    if(ifl==0):
        ax.set_ylim(0.93,1.07)
    if(ifl==1):
        ax.set_ylim(0.90,1.10)
    if(ifl==2):
        ax.set_ylim(0.95,1.08)
    if(ifl==3):
        ax.set_ylim(0.97,1.03)
    if(ifl==4):
        ax.set_ylim(0.91,1.09)
    if(ifl==5):
        ax.set_ylim(0.82,1.15)
     
    ax.tick_params(which='both',direction='in',labelsize=14,right=True)
    ax.tick_params(which='major',length=8)
    ax.tick_params(which='minor',length=6)
    ax.set_ylabel(labelpdf[ifl],fontsize=16)
    if(ifl==1):
        ax.set_xlabel(r'$x$',fontsize=18)
    
    if(ifl==1):
        ax.text(0.44,0.12,stringQ,fontsize=19,transform=ax.transAxes)

    # Add the legend
    if(ifl==2):
        ax.legend([(p1[0],p2[0]),(p3[0],p4[0])],\
                  [r"${\rm Current~Data~(NNPDF4.0)}$",\
                   r"${\rm Current~Data+FASER}\nu 2$"], \
                  frameon=True,loc=2,prop={'size':15})
        
    icount = icount + 1

py.tight_layout(pad=1.4, w_pad=1.0, h_pad=1.0)
py.savefig('FASERnu2-'+filelabel+'.pdf')
print('output plot: FASERnu2-'+filelabel+'.pdf')

exit()




