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
nmx =15
mxmin = 2 # GeV
mxmax = 301 # GeV
sqrts = 100000 # GeV, center of mass energy
s = math.pow(sqrts,2.0)
# Integration settings
epsrel=1e-3

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

nset=2
pdfset=["NNPDF31_nnlo_as_0118","NNPDF31_nnlo_as_0118"]
pdfsetlab=[r"${\rm NNPDF3.1}$",\
           r"${\rm NNPDF3.1+FPF(charm)}$"]
filelabel="-fpf-nnpdf31-electron-FCC100"

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
    elif(lumitype=="pp"): # photon-photon
        xpdf1 = p.xfxQ(22,x1,q)
        xpdf2 = p.xfxQ(22,x2,q) 
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
nlumi=2
lumichannel=["gg","qq"]
# Set mx grid (log spaced)
mx = np.logspace(np.log10(mxmin),np.log10(mxmax),nmx)
print("mx = ",mx)

# Read the weights, 2% uncertainties
weights_31_electron = np.array([4.20959881738763e-11, 4.464387872762778e-15, 3.1710968727146294, 5.001055598481901e-40, 1.9233942360058736, 0.00522789677087585, 1.9381615486918123e-22, 6.948282303461209e-24, 0.06077148947256093, 0.02717385269034323, 2.7237343721995346e-07, 6.863857770437248e-08, 0.43768685340824637, 4.1745073051576e-56, 2.1674571660035514e-15, 0.0019915940908779625, 1.7514222525832988e-19, 1.56326409800828e-07, 2.99190186192526e-56, 6.345200573203433, 2.431594158442075e-72, 1.4876699734956801e-170, 1.0398008089287094e-28, 0.0795695277747303, 6.2966481209422165e-59, 8.266158483455245, 7.729730861058763e-11, 2.5349671453327466e-07, 7.3220973698555865, 0.05385970028056135, 3.9746962420487293, 5.4610775493535926e-08, 4.799952740474192e-09, 9.005135679581303e-17, 0.004779185280047421, 1.2188576570720227e-89, 0.08015797401432266, 2.673495381955513e-49, 1.4748764295804374e-54, 0.0413821314569009, 5.423206610433918e-26, 2.6072926891915743e-74, 0.3941385012788396, 1.0340939078530515e-69, 2.1442787185174665e-24, 4.74158005792525e-26, 5.609637823505153, 2.1392949939154007e-28, 2.013309843721551e-12, 0.022065932719304307, 2.9751983232020177e-17, 4.3592296923194025e-07, 0.0003895547202838083, 1.2813193941409293e-180, 1.8898708875105342e-08, 1.4826239836861598e-33, 1.3871613281381383e-116, 0.0, 8.72768736314932e-34, 1.7490925716947977e-123, 6.815676196382348, 4.673600067043014e-21, 0.0, 5.8300490016459126, 7.111216628408524, 2.1175426478333522e-44, 1.4045438863414828, 1.2308862814967925, 0.4863808327603043, 4.844882729083789, 9.753122227482136e-39, 1.0919961285872696e-10, 0.026631523859397675, 0.10380861410599444, 0.014768734579598864, 7.39147817289468e-35, 3.260880607408281e-05, 1.3250863804661024e-144, 3.7915201995068145e-37, 6.782391495693524e-08, 0.0002041913500827248, 0.0019915940908779625, 9.539232006441736, 0.0021576736648114066, 4.332549527854453e-05, 5.530851694716551e-12, 3.1590279147376088e-27, 5.734348239464141e-26, 4.6137199966092865e-07, 9.768635107801382, 3.548206954550099e-69, 0.06833888681115605, 8.511266347478142e-37, 9.594516110871194, 5.322736088990339, 3.942276229475236e-123, 8.770779990212123e-26, 0.011790387766815947, 2.3744309148936213e-40, 1.648933560894321e-34])

# Check normalisation
sum  =0
for irep in range(100):
    sum = sum + weights_31_electron[irep]
if(abs(sum - 100) > 1e-5):
    print("weights not normalised")
    exit()

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
            
# end run over sets 
print("lumi arrays succesfully filled")

#---------------------------------------------------------------------
# Compute central values and uncertainties
#---------------------------------------------------------------------

lumi1_high = np.mean(lumi1,axis=2) + np.std(lumi1,axis=2)
lumi1_low = np.mean(lumi1,axis=2) - np.std(lumi1,axis=2)
lumi1_mid = np.mean(lumi1,axis=2)
lumi1_error= np.std(lumi1,axis=2)

lumi2_high = np.mean(lumi2,axis=2) + np.std(lumi2,axis=2)
lumi2_low = np.mean(lumi2,axis=2) - np.std(lumi2,axis=2)
lumi2_mid = np.mean(lumi2,axis=2)
lumi2_error= np.std(lumi2,axis=2)

cv = np.zeros((nlumi,nmx))
err = np.zeros((nlumi,nmx))
for imx in range(nmx):
    for ilumi in range(nlumi):
        cv[ilumi][imx] = 0.0
        err[ilumi][imx] = 0.0
        for irep in range(100):
            print(irep, " ",imx," ",ilumi)
            cv[ilumi][imx] = cv[ilumi][imx] + weights_31_electron[irep] * lumi2[ilumi][imx][irep]/100
            err[ilumi][imx] = err[ilumi][imx] + weights_31_electron[irep] * pow(lumi2[ilumi][imx][irep],2.0)/100
        print("err 1 = ",err[ilumi][imx])
        err[ilumi][imx] = pow( abs( err[ilumi][imx] - pow(cv[ilumi][imx],2.0)), 0.5)
        print("err final = ",err[ilumi][imx])
        #print(ilumi," ",imx," ",cv[ilumi][imx])
        #print(cv)
        #if(imx>1 and ilumi > 0):
        #    exit()

lumi2_mid = cv
lumi2_error = err
lumi2_high = cv + err
lumi2_low = cv -err


#*****************************************************************************
#*****************************************************************************

#---------------------------------------------------------------------
# Plot lumis
#---------------------------------------------------------------------
print("\n ****** Plotting relative luminosities ******* \n")

ncols,nrows=2,1
py.figure(figsize=(ncols*5,nrows*3.5))
gs = gridspec.GridSpec(nrows,ncols)
rescolors = py.rcParams['axes.prop_cycle'].by_key()['color']

#lumichannel=["gg","qg","qq"]
# lumilabels
ylabel=[r"$\mathcal{L}_{gg}/\mathcal{L}_{gg}^{(\rm ref)}$",\
        r"$\mathcal{L}_{q\bar{q}}/\mathcal{L}_{q\bar{q}}^{(\rm ref)}$"]

# y axis ranges
yranges=[[0.5,1.5],[0.5,1.5],[0.5,1.5]]

for ilumi in range(nlumi):

    norm = abs(lumi1_mid[ilumi]) # avoid crossing zero
    # print(norm)

    ax = py.subplot(gs[ilumi])
    p1=ax.plot(mx,lumi1_mid[ilumi]/norm,ls="solid",color=rescolors[0])
    ax.fill_between(mx,lumi1_high[ilumi]/norm,lumi1_low[ilumi]/norm,color=rescolors[0],alpha=0.2)
    p2=ax.fill(np.NaN,np.NaN,color=rescolors[0],alpha=0.2)

    norm = abs(lumi2_mid[ilumi]) # avoid crossing zero
    p3=ax.plot(mx,lumi2_mid[ilumi]/norm,ls="dashed",color=rescolors[1])
    ax.fill_between(mx,lumi2_high[ilumi]/norm,lumi2_low[ilumi]/norm,color=rescolors[1],alpha=0.2)
    p4=ax.fill(np.NaN,np.NaN,color=rescolors[1],alpha=0.2)

    ax.set_xscale('log')
    ax.set_xlim(mxmin,mxmax)
    ax.set_ylim(yranges[ilumi][0],yranges[ilumi][1])
    ax.tick_params(which='both',direction='in',labelsize=12,right=True)
    ax.tick_params(which='major',length=7)
    ax.tick_params(which='minor',length=4)
    if( ilumi ==1 ):
        ax.legend([(p1[0],p2[0]),(p3[0],p4[0])],\
                  [pdfsetlab[0],pdfsetlab[1]], frameon=True,loc="best",prop={'size':13})

    ax.set_xlabel(r'$m_X~({\rm GeV})$',fontsize=17)
    ax.set_ylabel(ylabel[ilumi],fontsize=19)
    if(ilumi==0):
        string=r'$\sqrt{s}=100'+'~{\\rm TeV}$'
        ax.text(0.50,0.10,string,fontsize=20,transform=ax.transAxes)


py.tight_layout(pad=1, w_pad=1, h_pad=1.0)
py.savefig('lumi'+filelabel+'.pdf')
print('output plot: lumi'+filelabel+'.pdf')



