# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:08:21 2016

@author: Virginia Price
"""

# import necessary packages #
# ------------------------- #
import numpy as np
from astropy.table import Table,join
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy.optimize import curve_fit
from myfunctions import *


# Have some conversion factors #
# ---------------------------- #

# GALEX counts per second to erg/cm^2/s/ang
cps2inst = 1.4e-15

# Ha flux order of magntitude (erg/cm^2/s)
Haord = 1.e-17






# read data in from file #
# ---------------------- #
if 'xover_raw' in locals():
    # Do nothing
    print('Data already loaded')
else:
#    ''' Here we see all the raw data files before compilation into xover, in case xover must be regenerated'''


    ''' RAW DATASETS imported directly from the interwebs'''
    
    NO= ascii.read('/users/virginiaprice/Documents/Research/OQ-project/Data/kias1033_5_P-MJD-F_MPAJHU_Zall_stellarMass_BPT_SFR_NSA.txt',
    include_names=('plate','fiberID','MJD','vflag','redshift','absmag','Z12logOH','Z12logOppH','N12logNH','logNO','BPTclass','t3'))
    NO.rename_column('fiberID','fiberid')
        
    NO = NO[~np.isnan(NO['logNO'])]     #Delete any galaxies w/o metallicity measurements
    NO = NO[NO['Z12logOH'] > 6.5]
    NO = NO[NO['BPTclass']==1.]
    NO = NO[NO['t3']<3.]   
    
    
    # Flux in units of maggys
    PETR = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/PETR_r.csv',delimiter=',')
    PETR=PETR[PETR['petroMag_r']>0.]
    
    '''not actually used any more'''
    # Calibrated fuv flux in uJy #
    FUV_raw = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/gal_fluxdata.csv',delimiter=',')
    FUV = FUV_raw[FUV_raw['fuv_flux']>0.]
    FUV['fuv_flux']=FUV['fuv_flux']*1e-17 #convert from uJy to erg cm^-2 s^-1 A^-1
    '''end'''
    
    Ha = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/Ha_data_comp.txt')
    Ha = Ha[Ha['Ha_flux']>0.]        # Currently in 1e-17 erg/cm^2/s
    Ha['Ha_flux'] = Ha['Ha_flux']*1.e-17    # Convert to erg/cm^2/s


    
    # Flux in units of maggys
    r_fibflux = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/red_fibmags.csv',delimiter=',')
    r_fibflux = r_fibflux[r_fibflux['fiberMag_r']>0]
    
    
    # GALEX aperture data    
    aper = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/gal_aperdata.csv',delimiter=',')
    aper_err = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/gal_aperdata_err.csv',delimiter=',')
    
    
    
    '''Altered datasets  -- defined by "xover" filenames'''
    # Metallicities and SFR data
    MT = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/compiled/xover_data_v0.txt')    
    
    # Metallicities, SFR, and flux data (corrected for aperture)
    aper_cor = ascii.read('/Users/virginiaprice/Documents/Research/OQ-project/Data/compiled/xover_data_v1.txt')



# Put data into tables #
# -------------------- #
# these are the old metallicity tables
'''if 'NO_dat' in locals():
    # Do nothing
    print('Data already sorted')
else:
     Define column labels
    NO_cols = ('index','MPA_index','plate','fiberid','MJD','ra','dec','vflag','absmag','Z12logOH','Z12logOppOH','N12logNH','logNO')
    FR_cols = ('plate','fiberid','MJD','ra','dec','FUV_SFR','Halpha_SFR')

    NO_dat = Table(NO[:,[0,33,30,32,31,1,2,29,5,35,37,39,41]],names = NO_cols)
    FR_dat = Table(FR[:,[0,1,2,5,6,10,11]], names = FR_cols)
    NO_dat = NO_dat[~np.isnan(NO_dat['logNO'])]     #Delete any galaxies w/o metallicity measurements
    NO_dat = NO_dat[NO_dat['Z12logOH'] > 6.5]'''
    
  
print('Done loading data')

# Merge uncorrected & corrected data
# First reduce the GALEX data set by averaging over fluxes
ap_red = aper_reduce(aper,'data')
ap_red_err = aper_reduce(aper_err,'error')

# Find FUV flux within petrosian galaxy radius
FUV_petr = fuv_petr(ap_red,ap_red_err,PETR)

# Perform aperture correction in H-alpha
SB = sb(Ha,PETR,r_fibflux,FUV_petr,'y')

flux_data = join(NO,SB,keys=('plate','fiberid','MJD'))
flux_data['flux_ratio'] = flux_data['Ha_petr']/flux_data['FUV_petr']

# Add a few more cuts








# Use the magic function 'join' to join the datasets #
# ------------------------------------------------------ #
#xover = join(NO_dat,FR_dat,keys=('plate','fiberid','MJD'))
#fluxmatch = join(r_fibflux,aper,keys=('objID'))
#fxmt = join(fluxmatch,PETR,keys=('objID'))
#print(str(len(xover['index'])) + ' galaxies in sample')

#Make void, wall arrays for more readable calling
#walls = xover[xover['vflag'] == 0.]
#voids = xover[xover['vflag'] == 1.]
#other = xover[xover['vflag'] == 2.]

#Ratios of FUV/Halpha rate
#wallsfrratio = np.divide(walls['FUV_sfr'],walls['Halpha_sfr'])
#voidsfrratio = np.divide(voids['FUV_sfr'],voids['Halpha_sfr'])



