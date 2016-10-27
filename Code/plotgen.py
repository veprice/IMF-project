# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:00:27 2016

@author: virginiaprice
"""
''' LIST OF GENERATED PLOTS SO FAR
Uncomment to generate plots from script'''
#
#''' ##################################### '''
'''# Basic metallicity ratios vs SFR'''
#
'''Halpha'''
'''------------------------'''
#'''First find correlation coefficients between matrices'''
#cr_wH=np.corrcoef(walls['logNO'],walls['Halpha_SFR'])[0,1] 
#cr_vH = np.corrcoef(voids['logNO'],voids['Halpha_SFR'])[0,1]
#
#'''Now plot the thing'''
#plt.figure()
#plt.scatter(walls['logNO'],walls['Halpha_SFR'],c='b')
#plt.scatter(voids['logNO'],voids['Halpha_SFR'],c='r')
#plt.xlabel(r'Metallicity Ratio ($\log(N/O)$)')
#plt.ylabel(r'$SFR_{H\alpha}$')
#plt.legend(['Wall','Void'])
#plt.text(0.05,9.5,r'$\rho_{wall} =$ '+str(cr_wH))
#plt.text(0.05,8.5,r'$\rho_{void} =$ '+str(cr_vH))
#
#plt.show()
#
'''FUV'''
'''------------------------'''
#
#'''First find correlation coefficients between matrices'''
#cr_wuv=np.corrcoef(walls['logNO'],walls['FUV_SFR'])[0,1]
#cr_vuv=np.corrcoef(voids['logNO'],voids['FUV_SFR'])[0,1]
#
#''' Now plot the thing '''
#plt.figure()
#plt.scatter(walls['logNO'],walls['FUV_SFR'],c='b')
#plt.scatter(voids['logNO'],voids['FUV_SFR'],c='r')
#plt.xlabel(r'Metallicity Ratio ($\log(N/O)$)')
#plt.ylabel(r'$SFR_{FUV}$')
#plt.legend(['Wall','Void'])
#plt.text(0.05,9.5,r'$\rho_{wall} =$ '+str(cr_wuv))
#plt.text(0.05,8.5,r'$\rho_{void} =$ '+str(cr_vuv))
#
#plt.show()

'''Ratio'''
'''------------------------'''
'''First find correlation coefficients between matrices'''
cr_w=np.corrcoef(walls['logNO'],np.log10(wallsfrratio))[0,1]
cr_v=np.corrcoef(voids['logNO'],np.log10(voidsfrratio))[0,1]

# Now plot the thing
plt.figure()
plt.scatter(walls['logNO'],wallsfrratio,c='b')
plt.scatter(voids['logNO'],voidsfrratio,c='r')
plt.xlabel(r'Metallicity Ratio ($\log(N/O)$)')
plt.ylabel(r'$SFR_{FUV}/SFR_{H\alpha}$ (log)')
plt.yscale('log')
plt.ylim(1e-3,1e2)
plt.legend(['Wall','Void'])
#plt.text(-0.05,50,r'$\rho_{wall} =$ '+str(cr_w))
#plt.text(-0.05,30,r'$\rho_{void} =$ '+str(cr_v))


#''' ##################################### '''
''' Histograms of abs magnitudes '''
#plt.figure()
#plt.hist(xover['absmag'],80,color='g')
#plt.xlabel('Absolute Magnitude')
#plt.ylabel('Number of Galaxies')
#plt.show()
#
#plt.figure()
#plt.hist(walls['absmag'],80,color='b')
#plt.xlabel('Absolute Magnitude')
#plt.ylabel('Number of Wall Galaxies')
#plt.show()
#
#plt.figure()
#plt.hist(voids['absmag'],50,color='r')
#plt.xlabel('Absolute Magnitude')
#plt.ylabel('Number of Void Galaxies')
#plt.show()


''' Reproduction of Muerer et al '''
plt.figure()
plt.scatter(np.log10(fd_c['Ha_SB']),np.log10(fd_c['flux_ratio']),c='red')
plt.xlabel(r'$log(\Sigma_{H\alpha}$ [$\mathrm{W\;kpc}^{-2})$]')
plt.ylabel(r'$log(F_{H\alpha}/f_{FUV}$ $[\mathrm{\AA}]$)')