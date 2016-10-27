# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:16:35 2016

@author: virginiaprice

All the custom functions I've made for everything

List as follows:

aper_reduce -- averages over multiple galaxies in GALEX dataset
devline -- Can calculate standard deviation from linear trend
flux_fit -- Defines a fit function to fit GALEX aperture data
fuv_petr -- finds fuv flux in petrosian radius
ha_sb -- Calculates average galactic surface brightness inside petrosian radius
hub_dist -- Converts from arcsec to kpc

"""

import numpy as np
from astropy.table import Table,join
from scipy.optimize import curve_fit


    

def aper_reduce(aper,dataorerror):
    objids = np.unique(aper['objID'])
    zeros = np.zeros((len(objids),len(aper.colnames)))
    zeros[:,0] = objids

    aper_red = Table(zeros,names=aper.colnames)

    
    if dataorerror=='uncorr':
        for k in xrange(len(aper_red)):
            curr=aper[aper['objID']==aper_red['objID'][k]]
            avg_flx = np.mean(curr['fuv_flux'])
            aper_red['fuv_flux'][k]=avg_flx
    else:
        aper_red.remove_columns(['plate','fiberid','MJD'])
        if dataorerror == 'data':
            basestr = 'FUV_FLUX_APER_'
        elif dataorerror=='error':
            basestr = 'FUV_FLUXERR_APER_'
        
    
                
        
        for i in xrange(len(aper_red)):      #loop through all object IDs
    
            curr = aper[aper['objID']==aper_red['objID'][i]]    #picks out all objects with the same objID
            curr = curr[curr[basestr + '1'] >= 0.] # gets rid of any -99 flux error values because we don't need that sort of negativity
            avg_aperj = []      # start an array for apertures 1-6
            for j in xrange(6):
                stselect = basestr + str(j+1)
                
                avg_aperj=np.mean(curr[stselect])  
                aper_red[stselect][i]=avg_aperj

            
    return aper_red
    
    
def devline(SB,ha2fuv):
    logSB = SB
    logha2fuv = np.log10(ha2fuv)
    # Fit a line to the log of the data
    a = np.polyfit(logSB,logha2fuv,1)
    line = np.poly1d(a)
    slope = a[0]
    
    x = np.linspace(min(logSB),max(logSB),len(logSB))   # line arrays to plot
    y = line(x)

    
    subtr = logha2fuv-line(logSB)   # Subtract away linear dependence
    dev=np.std(subtr)               # Caluculate standard deviation

    pltdata = Table([logSB,logha2fuv,subtr,x,y],names=('logSB','logratio','subtracted','fitx','fity'))
    
    return pltdata,dev,slope
    

''' Returns fuv flux in the same aperture radius as the Ha petrosian radius'''  
def flux_fit(x,a,b,c):
    return (c-np.exp(a-b*x))
    
    
def fuv_petr(aper_red,aper_red_err,PETR):
    FUV_PETR = join(aper_red,PETR,keys='objID')
    Tbl_err = join(FUV_PETR,aper_red_err,keys='objID')
    FUV_PETR['FUV_petr'] = np.zeros(len(FUV_PETR))
    FUV_PETR['FUV_fib'] = np.zeros(len(FUV_PETR))

    # aperture radii are constant so define them here for fast code   
    ap_x = np.array([3.,4.5,7.5,12.,18.,25.5])/2.        # Divide by 2 since GALEX uses diameter for aperture, not radius
    for i in xrange(len(FUV_PETR)):      
        ap_curr = np.array([FUV_PETR['FUV_FLUX_APER_1'][i],FUV_PETR['FUV_FLUX_APER_2'][i],FUV_PETR['FUV_FLUX_APER_3'][i],FUV_PETR['FUV_FLUX_APER_4'][i],FUV_PETR['FUV_FLUX_APER_5'][i],FUV_PETR['FUV_FLUX_APER_6'][i]])
        ap_cerr = np.array([Tbl_err['FUV_FLUXERR_APER_1'][i],Tbl_err['FUV_FLUXERR_APER_2'][i],Tbl_err['FUV_FLUXERR_APER_3'][i],Tbl_err['FUV_FLUXERR_APER_4'][i],Tbl_err['FUV_FLUXERR_APER_5'][i],Tbl_err['FUV_FLUXERR_APER_6'][i]])
        try:

            apfit,apfit_cov = curve_fit(flux_fit,ap_x,ap_curr)      # fit to function defined by flux_fit above

            a = flux_fit(FUV_PETR['petroRad_r'][i],apfit[0],apfit[1],apfit[2])  # Calculate FUV flux in petrosian radius of galaxy
            FUV_PETR['FUV_petr'][i] = a*1.4e-15  # convert from cps to erg/cm^2/s/angstrom  
            FUV_PETR['FUV_fib'][i] = FUV_PETR['FUV_FLUX_APER_1'][i]*1.4e-15 #Gives FUV fiber flux in 3" aperture (same units as above)
            
            xrng = np.linspace(0,max(ap_x),200)
            curv = flux_fit(xrng,apfit[0],apfit[1],apfit[2])
            
            
#            if i==0:
#                plt.figure()
#                plt.plot(xrng,curv,c='r',ls='dashdot')
#                plt.scatter(ap_x,ap_curr)
#                plt.errorbar(ap_x,ap_curr,yerr=ap_cerr,fmt='none')
#                plt.xlabel(r'Fiber radius (arcsec)')
#                plt.ylabel('FUV flux (counts per second)')
#                plt.plot((FUV_PETR['petroRad_r'][i], FUV_PETR['petroRad_r'][i]), (0, .6), 'k-')
#                plt.scatter(FUV_PETR['petroRad_r'][i],a,marker='x',c='g')
#                plt.show()
            
        except:
            FUV_PETR['FUV_petr'][i]=-99.
    # Clean up output table a little bit    
    FUV_PETR = FUV_PETR[FUV_PETR['FUV_petr']>0.]
    FUV_PETR.keep_columns(('objID','plate','fiberid','MJD','FUV_petr','FUV_fib'))
    return FUV_PETR
    
    


def sb(Ha,PETR,r_fibflux,FUV_petr,correct):
# Inputs:
# Ha - Table with SDSS H-alpha fluxes (erg/cm^2/s)
# PETR - Table with petrosian fluxes & radii in r-band
# r_fibflux - table with fiber fluxes in r-band
# correct - string input, y/n -- correct units or no?
    Ha_PETR = join(Ha,PETR,keys=('plate','fiberid'))
    Ha_SB = join(Ha_PETR,r_fibflux,keys=('objID','plate','fiberid','MJD'))
    SB = join(Ha_SB,FUV_petr,keys=('objID','plate','fiberid','MJD'))
    
    SB['Ha_petr'] = np.zeros(len(SB)) 
    
    SB['Ha_SB'] = np.zeros(len(SB)) 
    SB['r_SB'] = np.zeros(len(SB))
    SB['FUV_SB'] = np.zeros(len(SB)) 

            
    for i in xrange(len(SB)):
        SB['Ha_petr'][i] = SB['petroMag_r'][i]/SB['fiberMag_r'][i]*SB['Ha_flux'][i]                                                                            # note conversion from arsec --> radians
        
        if correct == 'y':
          
            SB['Ha_SB'][i] = 2.*SB['Ha_petr'][i]/SB['petroRad_r'][i]**2*8.1e46
            
            SB['r_SB'][i] = 2.*SB['petroMag_r'][i]/SB['petroRad_r'][i]**2*8.1e46
            SB['FUV_SB'][i] = 2.*SB['FUV_petr'][i]/SB['petroRad_r'][i]**2*8.1e46
        else:
            SB['Ha_SB'][i] = 2.*SB['Ha_flux'][i]/SB['petroRad_r'][i]**2*8.1e46
    return SB

def datrange(data):
    print(min(data),max(data))