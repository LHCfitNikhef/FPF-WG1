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
import scipy.integrate as integrate
import scipy.special as special

#---------------------------------------------------------
#---------------------------------------------------------
# General plot settings
nmx = 30
mxmin = 9 # GeV
mxmax = 7000 # GeV
sqrts = 14000 # GeV, center of mass energy
s = math.pow(sqrts,2.0)
# Integration settings
epsrel=5e-3

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

nset=3
pdfset=["230309-ern-002","230824-tr-fpf-002","230824-jcm-fpf-001_iterated"]

pdfsetlab=[r"${\rm NNPDF4.0}$",\
                   r"${\rm NNPDF4.0}+{\rm FPF~(stat\,only)}$",\
                   r"${\rm NNPDF4.0}+{\rm FPF~(stat+sys)}$"]
error_option=["mc_1sigma","mc_1sigma","mc_1sigma"]
filelabel="-FPFall"

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

def lumi_wrap(y):

    # number of active flavours
    nfl=4 # up to charm (inclusive)
    
    # Define kinematics 
    x1 = y
    x2 = tau/x1
    q = mx[imx]

    # Check x1 and x2 are physical
    if( x1 < 0 or x1 > 1 or x2 < 0 or x2 > 1 ):
        print("invalid x1, x2 = ",x1,x2)
        exit()
    if( q < 1 or q > 1e6 ):
        print("invalid q = ",q)
        exit()
    
    # Evaluate lumis
    if(lumitype=="gg"): # gluon-gluon luminosity
        xpdf1 = p.xfxQ(0,x1,q)
        xpdf2 = p.xfxQ(0,x2,q) 
        integrand = (1.0/x1)*(xpdf1/x1)*(xpdf2/x2)
    elif(lumitype=="qg"): # quark-gluon luminosity
        integrand = 0
        xsinglet = 0
        for ifl in range(1,nfl+1):
            xsinglet = xsinglet + ( p.xfxQ(ifl,x1,q) + p.xfxQ(-ifl,x1,q) )
        xpdf2 = p.xfxQ(0,x2,q) 
        integrand = (1.0/x1)*(xsinglet/x1)*(xpdf2/x2)
    elif(lumitype=="qq"): # quark-antiquark
        integrand = 0
        for ifl in range(1,nfl+1):
            xpdf1 = p.xfxQ(ifl,x1,q)
            xpdf2 = p.xfxQ(-ifl,x2,q) 
            integrand = integrand + (1.0/x1)*(xpdf1/x1)*(xpdf2/x2)
    elif(lumitype=="q2"): # quark-quark
        integrand = 0
        for ifl in range(1,nfl+1):
            for jfl in range(1,nfl+1):
                xpdf1 = p.xfxQ(ifl,x1,q)
                xpdf2 = p.xfxQ(jfl,x2,q) 
            integrand = integrand + (1.0/x1)*(xpdf1/x1)*(xpdf2/x2)
    elif(lumitype=="ss"): # strange-antistrange
        xpdf1 = p.xfxQ(3,x1,q)
        xpdf2 = p.xfxQ(-3,x2,q) 
        integrand = (1.0/x1)*(xpdf1/x1)*(xpdf2/x2)
    elif(lumitype=="sg"): # strange-gluon
        xpdf1 = p.xfxQ(3,x1,q)+p.xfxQ(-3,x1,q)
        xpdf2 = p.xfxQ(0,x2,q) 
        integrand = (1.0/x1)*(xpdf1/x1)*(xpdf2/x2)
    else:
        print("invalid type of lumi = ",lumitype)
        exit()

    # Return integrand
    return integrand

     
#---------------------------------------------------------
#---------------------------------------------------------
# Common plotting settings

# Reduce verbosity of LHAPDF
lhapdf.setVerbosity(0)
# number of lumi combinations
nlumi=3
lumichannel=["gg","qq","q2"]
# Set mx grid (log spaced)
mx = np.logspace(np.log10(mxmin),np.log10(mxmax),nmx)

print("mx = ",mx)

# Threshold for ratios
eps=1e-60

print("\n Filling lumi arrays from LHAPDF \n")

nrep = np.zeros(nset, dtype='int')
# run over PDF sets
for iset in range(nset):
    # Initialise PDF set
    p=lhapdf.getPDFSet(pdfset[iset])
    print(p.description)
    nrep[iset] = int(p.get_entry("NumMembers"))-1
    print("nrep[iset] =", nrep[iset])
    # Arrays to store LHAPDF results
    if(iset==0):
        lumi1 = np.zeros((nlumi,nmx,nrep[iset]))
        lumi1_cv = np.zeros((nlumi,nmx))
    elif(iset==1):
        lumi2 = np.zeros((nlumi,nmx,nrep[iset]))
        lumi2_cv = np.zeros((nlumi,nmx))
    elif(iset==2):
        lumi3 = np.zeros((nlumi,nmx,nrep[iset]))
        lumi3_cv = np.zeros((nlumi,nmx))
        
    # Run over replicas
    for irep in range(1,nrep[iset]+1):
        p=lhapdf.mkPDF(pdfset[iset],irep)
        print(pdfset[iset], "  ",irep)

        # run over flavours
        for ilumi in range(nlumi):
                        
            # Run over x arrat
            for imx in range(nmx):
                
                tau = mx[imx]**2/s
                if(tau < 0 or tau > 1):
                    print("invalid value of tau = ",tau)
                    exit()
                # Define type of luminosity    
                lumitype=lumichannel[ilumi]
                lumi = integrate.quad(lambda y: lumi_wrap(y), tau, 1,epsrel=epsrel)
                # Divide by center of mass energy
                if(iset==0):
                    lumi1[ilumi,imx,irep-1] = lumi[0]/s
                elif(iset==1):
                    lumi2[ilumi,imx,irep-1] = lumi[0]/s
                elif(iset==2):
                    lumi3[ilumi,imx,irep-1] = lumi[0]/s
               

    # Central values
    p=lhapdf.mkPDF(pdfset[iset],0)
    print(pdfset[iset], "  ",irep)
    # run over flavours
    for ilumi in range(nlumi):
        
        # Run over x arrat
        for imx in range(nmx):
            
            tau = mx[imx]**2/s
            if(tau < 0 or tau > 1):
                print("invalid value of tau = ",tau)
                exit()
            # Define type of luminosity    
            lumitype=lumichannel[ilumi]
            lumi = integrate.quad(lambda y: lumi_wrap(y), tau, 1,epsrel=epsrel)
            # Divide by center of mass energy
            if(iset==0):
                lumi1_cv[ilumi,imx] = lumi[0]/s
            elif(iset==1):
                lumi2_cv[ilumi,imx] = lumi[0]/s
            elif(iset==2):
                lumi3_cv[ilumi,imx] = lumi[0]/s
             
# end run over sets 
print("lumi arrays succesfully filled")

#---------------------------------------------------------------------
# Compute central values and uncertainties
#---------------------------------------------------------------------

for iset in range(nset):

    # MC PDF sets, 68% CL intervals
    if(error_option[iset]=="mc_68cl"):

        if(iset==0):
            lumi1_high = np.nanpercentile(lumi1,84,axis=2)
            lumi1_low = np.nanpercentile(lumi1,16,axis=2)
            lumi1_mid = ( lumi1_high + lumi1_low )/2.
            lumi1_error = ( lumi1_high - lumi1_low )/2.

        elif(iset==1):
            lumi2_high = np.nanpercentile(lumi2,84,axis=2)
            lumi2_low = np.nanpercentile(lumi2,16,axis=2)
            lumi2_mid = ( lumi2_high + lumi2_low )/2.
            lumi2_error = ( lumi2_high - lumi2_low )/2.
            
        elif(iset==2):
            lumi3_high = np.nanpercentile(lumi3,84,axis=2)
            lumi3_low = np.nanpercentile(lumi3,16,axis=2)
            lumi3_mid = ( lumi3_high + lumi3_low )/2.
            lumi3_error = ( lumi3_high - lumi3_low )/2.
            
        else:
            print("invalid plotting option")
            exit()

    # MC PDF sets, 90% CL intervals
    elif(error_option[iset]=="mc_90cl"):

        if(iset==0):
            lumi1_high = np.nanpercentile(lumi1,95,axis=2)
            lumi1_low = np.nanpercentile(lumi1,5,axis=2)
            lumi1_mid = ( lumi1_high + lumi1_low )/2.
            lumi1_error = ( lumi1_high - lumi1_low )/2.

        elif(iset==1):
            lumi2_high = np.nanpercentile(lumi2,95,axis=2)
            lumi2_low = np.nanpercentile(lumi2,5,axis=2)
            lumi2_mid = ( lumi2_high + lumi2_low )/2.
            lumi2_error = ( lumi2_high - lumi2_low )/2.
            
        elif(iset==2):
            lumi3_high = np.nanpercentile(lumi3,95,axis=2)
            lumi3_low = np.nanpercentile(lumi3,5,axis=2)
            lumi3_mid = ( lumi3_high + lumi3_low )/2.
            lumi3_error = ( lumi3_high - lumi3_low )/2.

        else:
            print("invalid plotting option")
            exit()
     
    # MC PDF sets, one-sigma and mean
    elif(error_option[iset]=="mc_1sigma"):

        if(iset==0):
            lumi1_high = np.mean(lumi1,axis=2) + np.std(lumi1,axis=2)
            lumi1_low = np.mean(lumi1,axis=2) - np.std(lumi1,axis=2)
            lumi1_mid = np.mean(lumi1,axis=2)
            lumi1_error= np.std(lumi1,axis=2)

        elif(iset==1):
            lumi2_high = np.mean(lumi2,axis=2) + np.std(lumi2,axis=2)
            lumi2_low = np.mean(lumi2,axis=2) - np.std(lumi2,axis=2)
            lumi2_mid = np.mean(lumi2,axis=2)
            lumi2_error= np.std(lumi2,axis=2)

        elif(iset==2):
            lumi3_high = np.mean(lumi3,axis=2) + np.std(lumi3,axis=2)
            lumi3_low = np.mean(lumi3,axis=2) - np.std(lumi3,axis=2)
            lumi3_mid = np.mean(lumi3,axis=2)
            lumi3_error= np.std(lumi3,axis=2)

        else:
            print("invalid plotting option")
            exit()
            

    # Asymmetric Hessian with then normalisation to one-sigma
    elif(error_option[iset]=="ct" or error_option[iset]=="mmht" ):

        if(iset==0):
            lumi1_mid = np.mean(lumi1,axis=2)
            lumi1_error = np.std(lumi1,axis=2)
            neig = int(nrep[iset]/2) # Number of eigenvectors
            for imx in range(nmx):
                for ilumi in range(nlumi):
                    lumi1_mid[ilumi][imx]=lumi1_cv[ilumi][imx]
                    lumi1_error[ilumi][imx]=0 # initialisation
                    for ieg in range(neig):
                        lumi1_error[ilumi][imx] = lumi1_error[ilumi][imx] \
                            + (lumi1[ilumi,imx,2*ieg+1] - lumi1[ilumi,imx,2*ieg] )**2.0
                    lumi1_error[ilumi][imx]=math.sqrt(lumi1_error[ilumi][imx])/2
                    if(error_option[iset]=="ct"):
                        lumi1_error[ilumi][imx]=lumi1_error[ilumi][imx] / 1.642 # from 90% to 68% CL
            lumi1_high = lumi1_mid + lumi1_error 
            lumi1_low = lumi1_mid - lumi1_error

        elif(iset==1):
            lumi2_mid = np.mean(lumi2,axis=2)
            lumi2_error = np.std(lumi2,axis=2)
            neig = int(nrep[iset]/2) # Number of eigenvectors
            for imx in range(nmx):
                for ilumi in range(nlumi):
                    lumi2_mid[ilumi][imx]=lumi2_cv[ilumi][imx]
                    lumi2_error[ilumi][imx]=0 # initialisation
                    for ieg in range(neig):
                        lumi2_error[ilumi][imx] = lumi2_error[ilumi][imx] \
                            + (lumi2[ilumi,imx,2*ieg+1] - lumi2[ilumi,imx,2*ieg] )**2.0
                    lumi2_error[ilumi][imx]=math.sqrt(lumi2_error[ilumi][imx])/2
                    if(error_option[iset]=="ct"):
                        lumi2_error[ilumi][imx]=lumi2_error[ilumi][imx] / 1.642 # from 90% to 68% CL
            lumi2_high = lumi2_mid + lumi2_error 
            lumi2_low = lumi2_mid - lumi2_error
            
        elif(iset==2):
            lumi3_mid = np.mean(lumi3,axis=2)
            lumi3_error = np.std(lumi3,axis=2)
            neig = int(nrep[iset]/2) # Number of eigenvectors
            for imx in range(nmx):
                for ilumi in range(nlumi):
                    lumi3_mid[ilumi][imx]=lumi3_cv[ilumi][imx]
                    lumi3_error[ilumi][imx]=0 # initialisation
                    for ieg in range(neig):
                        lumi3_error[ilumi][imx] = lumi3_error[ilumi][imx] \
                            + (lumi3[ilumi,imx,2*ieg+1] - lumi3[ilumi,imx,2*ieg] )**2.0
                    lumi3_error[ilumi][imx]=math.sqrt(lumi3_error[ilumi][imx])/2
                    if(error_option[iset]=="ct"):
                        lumi3_error[ilumi][imx]=lumi3_error[ilumi][imx] / 1.642 # from 90% to 68% CL
            lumi3_high = lumi3_mid + lumi3_error 
            lumi3_low = lumi3_mid - lumi3_error
        else:
            print("invalid option")
            exit()

    # symmetric Hessian 
    elif(error_option[iset]=="symmhessian" ):

        if(iset==1):
            lumi2_mid = np.mean(lumi2,axis=2)
            lumi2_error = np.std(lumi2,axis=2)
            neig = nrep[iset] # Number of eigenvectors
            for imx in range(nmx):
                for ilumi in range(nlumi):
                    lumi2_mid[ilumi][imx]=lumi2_cv[ilumi][imx]
                    lumi2_error[ilumi][imx]=0 # initialisation
                    for ieg in range(neig):
                        lumi2_error[ilumi][imx] = lumi2_error[ilumi][imx] \
                            + (lumi2[ilumi,imx,ieg] - lumi2_cv[ilumi][imx])**2.0
                    lumi2_error[ilumi][imx]=math.sqrt(lumi2_error[ilumi][imx])
            lumi2_high = lumi2_mid + lumi2_error 
            lumi2_low = lumi2_mid - lumi2_error
        elif(iset==2):
            lumi3_mid = np.mean(lumi3,axis=2)
            lumi3_error = np.std(lumi3,axis=2)
            neig = nrep[iset] # Number of eigenvectors
            for imx in range(nmx):
                for ilumi in range(nlumi):
                    lumi3_mid[ilumi][imx]=lumi3_cv[ilumi][imx]
                    lumi3_error[ilumi][imx]=0 # initialisation
                    for ieg in range(neig):
                        lumi3_error[ilumi][imx] = lumi3_error[ilumi][imx] \
                            + (lumi3[ilumi,imx,ieg] - lumi3_cv[ilumi][imx])**2.0
                    lumi3_error[ilumi][imx]=math.sqrt(lumi3_error[ilumi][imx])
            lumi3_high = lumi3_mid + lumi3_error 
            lumi3_low = lumi3_mid - lumi3_error
        else:
            print("invalid option")
            exit()
               
    else:
        print("Invalid error option = ",error_option[iset])
        exit()
    
#----------------------------------------------------------------------


#*****************************************************************************
#*****************************************************************************

#---------------------------------------------------------------------
# Plot lumis
#---------------------------------------------------------------------
print("\n ****** Plotting relative luminosities ******* \n")

ncols,nrows=3,1
py.figure(figsize=(ncols*5,nrows*3.5))
gs = gridspec.GridSpec(nrows,ncols)
rescolors = py.rcParams['axes.prop_cycle'].by_key()['color']

# lumilabels
ylabel=[r"$\mathcal{L}_{gg}/\mathcal{L}_{gg}^{(\rm ref)}$",\
        r"$\mathcal{L}_{q\bar{q}}/\mathcal{L}_{q\bar{q}}^{(\rm ref)}$",\
        r"$\mathcal{L}_{qq}/\mathcal{L}_{qq}^{(\rm ref)}$"]

# y axis ranges
yranges=[[0.90,1.15],[0.95,1.06],[0.90,1.15]]

for ilumi in range(nlumi):

    norm = abs(lumi1_mid[ilumi]+eps) # avoid crossing zero
    # print(norm)

    # Only interested in relative reduction of uncertainties

    ax = py.subplot(gs[ilumi])
    p1=ax.plot(mx,norm/norm,lumi1_mid[ilumi],ls="solid",color=rescolors[0],lw=2)
    ax.fill_between(mx,lumi1_high[ilumi]/norm,lumi1_low[ilumi]/norm,color=rescolors[0],alpha=0.2)
    p2=ax.fill(np.NaN,np.NaN,color=rescolors[0],alpha=0.2)

    norm = abs(lumi3_mid[ilumi]+eps) # avoid crossing zero
    p3=ax.plot(mx,lumi3_high[ilumi]/norm,ls="dashed",color=rescolors[1])
    p4=ax.plot(mx,lumi3_low[ilumi]/norm,ls="dashed",color=rescolors[1])
    

    norm = abs(lumi2_mid[ilumi]+eps) # avoid crossing zero
    p5=ax.plot(mx,lumi2_high[ilumi]/norm,ls="solid",color=rescolors[2])
    p6=ax.plot(mx,lumi2_low[ilumi]/norm,ls="solid",color=rescolors[2])

    ax.set_xscale('log')
    ax.set_xlim(mxmin,mxmax)
    ax.set_ylim(yranges[ilumi][0],yranges[ilumi][1])
    ax.tick_params(which='both',direction='in',labelsize=12,right=True)
    ax.tick_params(which='major',length=7)
    ax.tick_params(which='minor',length=4)
    if (ilumi==0):
        ax.legend([(p1[0],p2[0]),(p5[0]),(p3[0])],\
                  [pdfsetlab[0],pdfsetlab[1],pdfsetlab[2]], frameon=True,loc=2,prop={'size':12})
    ax.set_xlabel(r'$m_X~({\rm GeV})$',fontsize=17)
    ax.set_ylabel(ylabel[ilumi],fontsize=19)
    if(ilumi==0):
        string=r'$\sqrt{s}=14'+'~{\\rm TeV}$'
        ax.text(0.10,0.08,string,fontsize=18,transform=ax.transAxes)

py.tight_layout(pad=1, w_pad=1, h_pad=1.0)
py.savefig('lumi'+filelabel+'.pdf')
print('output plot: lumi'+filelabel+'.pdf')



