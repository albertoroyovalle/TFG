#!/usr/bin/env python
# coding: utf-8

# """
# 11-04-2021:
# Para adaptar el codigo getting_sdss_profiles es necesario:
# 
# 1)eliminar el bucle original (for ii,gal in enumerate(gals): [1]).
# 
# 2)De ese bucle se leen datos de un catalogo. Los que necesitamos en este programa son:
# 
# gals=np.array(["79071"])
# gal="79071" -----------------> La razon de escribir asi gals y gal es para que no de muchos problemas el programa si el dia de mañana necesito implementar el bucle [1]
# field="hudf" 
# zz=0.6218
# pa=155
# ar=0.8
# 
# 3) Buscar y reemplazar [ii] (eliminando todos los terminos)
# 
# 4)A la hora de leer los datos de las bandas H J I B... rescordar que la salida de este programa nos daba documentos del tipo data_H.cat. 
# El programa originalmente lee:  tt = Table.read("./sb_profiles/"+gal_sb_cat,names=["dist_arcsec","dist_kpc","sb","sb_err"],format="ascii.commented_header")
# tener en cuenta que gal_sb_cat = str(gal)+"_"+field[ii]+"_"+band+".cat" (lo que difiere de nuestros archivos data_H.cat). Esto se arregla sustituendo gal_sb_cat por
# gal_sb_cat="data_"+band+".cat"
# 
# 
# """

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import pdb #for debugging purposes only
from astropy.table import Table
from scipy.interpolate import interp1d
from os import listdir
from os.path import isfile, join
import sys
###########################################
from my_python_library_v2 import correcting_by_extinction, correcting_by_inclination
from functs_for_truncations import obtaining_best_redshifts, obtaining_best_ar_and_pa


# In[2]:


lambda_sdss= np.array([4718.9,6185.2,7499.7,8961.5]) #http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php taking the lambda_mean
bands_sdss = np.array(["g","r","i","z"])
path_output="./sdss_profiles/"


gals=np.array(["79071"])
gal="79071"
field="hudf" ###CAMBIAR

zz=0.6218
pa=155
ar=0.8

if not os.path.exists(path_output):
    os.mkdir(path_output)


#subset of objects_to_analyse


# In[3]:


if field == "goodss":
    lambda_obs = np.array([4359.4,6035.8,8128.7,10651.0,12576.2,15436.3]) #http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php taking the lambda_mean
    bands_obs  = np.array(["B","V","I","Y","J","H"])
else:
    lambda_obs = np.array([6035.8,8128.7,12576.2,15436.3]) #http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php taking the lambda_mean
    bands_obs  = np.array(["V","I","J","H"])
flag_exists = np.full_like(bands_obs,True,dtype=bool)
mags     = { band : np.array([]) for band in bands_obs }
mags_err = { band : np.array([]) for band in bands_obs }
print("bands_obs: ", bands_obs)
print(flag_exists, "\n \n",mags,"\n \n",mags_err)
print(lambda_obs)


# In[5]:


#I load in memory all the surface brightness profiles
for jj,band in enumerate(bands_obs):
    
    #gal_sb_cat = str(gal)+"_"+field+"_"+band+".cat" Cambié esto ya que la entrada de mis funciones era distinta
    gal_sb_cat="data_"+band+".cat"
    if os.path.exists("./sb_profiles/"+gal_sb_cat):
        tt = Table.read("./sb_profiles/"+gal_sb_cat,names=["dist_arcsec","dist_kpc","sb","sb_err"],format="ascii.commented_header")
        vector_dist_arcsec = np.array(tt["dist_arcsec"])
        vector_dist_kpc    = np.array(tt["dist_kpc"])
        mags[band]         = np.array(tt["sb"])
        mags_err[band]     = np.array(tt["sb_err"])
        if np.isfinite(mags[band]).any() == False: flag_exists[jj] = False #if the image is empty, flag it
    else:
        flag_exists[jj] = False

    if flag_exists[jj] == True:
        
        mags[band] = correcting_by_extinction(mags[band],field,band)
    
        mags[band] = correcting_by_inclination(mags[band], ar)
print(flag_exists)
print(mags)


# 

# In[6]:



bands_obs  = bands_obs[flag_exists] #I only take the bands that exist
lambda_obs = lambda_obs[flag_exists]
    
lambda_restframe = lambda_obs/(1.+zz) #de-redshifting wavelengths
    
mags_corr = { band: np.full(vector_dist_arcsec.size,float("NaN")) for band in bands_obs }
mags_err2 = { band: np.full(vector_dist_arcsec.size,float("NaN")) for band in bands_obs }
    
for jj,band in enumerate(bands_obs):
    mags_corr[band] = mags[band] - 10.*np.log10(1.+zz) #correcting by cosmological dimming a la Giavalisco+04 and Ribeiro+16 --> and as you already use lambda restframe, you need to use (1+z)^4
    #I don't work in flux as before because the galaxy SED is in a logarithmic space (magnitudes) and there it is the place where I should approximate by linear relations
    #before flux_in_Jy_restframe = 10.**( (mags_corr-8.9)/(-2.5) ) #from mag=-2.5log(flux_in_Jy)-2.5log(10^-23)-48.6 Oke & Gunn 1983 ;I will do the extrapolation in flux and not in magnitude just in case there was any problem for interpolating logarithms
    mags_err2[band] = mags_err[band]

print(bands_sdss)
print(flag_exists)


# In[7]:


for jj,band in enumerate(bands_sdss):
    
    #you take the two filters that are in between lambda_sdss
    comp_wav = lambda_restframe - lambda_sdss[jj]
    filter_sign_wav = comp_wav <= 0.
    if filter_sign_wav.any() == False: filter_sign_wav[1]  = True  #to take the first two filters
    if filter_sign_wav.all() == True:  filter_sign_wav[-1] = False #to take the last two filters
    for kk,sign in enumerate(filter_sign_wav):
        if kk == 0:
            last_sign = sign
        else:
            if sign != last_sign:
                ind1 = kk-1
                ind2 = kk
                break
    
    """if you take the smallest difference in wavelength
    dist_wav = np.abs(lambda_restframe - lambda_sdss[jj])
    ind1 = np.where( dist_wav == np.min(dist_wav) )
    if ind1 == []: continue
    ind1 = int(ind1[0]) #[0] because it returns a tuple, int because it returns an array
    dist_wav[ind1] = 9e99
    ind2 = np.where( dist_wav == np.min(dist_wav) )
    if ind2 == []: continue
    ind2 = int(ind2[0]) #[0] because it returns a tuple, int because it returns an array
    """
    
    band1 = bands_obs[ind1]
    band2 = bands_obs[ind2]
    
    mag1 = mags_corr[band1]
    err_mag1 = mags_err2[band1]
    mag2 = mags_corr[band2]
    err_mag2 = mags_err2[band2]
        
    mags_new       = np.zeros_like(vector_dist_arcsec)
    mags_new_error = np.zeros_like(vector_dist_arcsec)
    
    mags_to_analyse=np.array([mag1,mag2])
    wavelengths_to_analyze=np.array([np.full_like(vector_dist_arcsec,lambda_restframe[ind1]),np.full_like(vector_dist_arcsec,lambda_restframe[ind2])])
    
    for kk in range(len(vector_dist_arcsec)):    
        mags_to_analyse=np.array([mag1[kk],mag2[kk]])
        wavelengths_to_analyze=np.array([lambda_restframe[ind1],lambda_restframe[ind2]])

        # Interpolate the filter response to data wavelength
        interp   = interp1d(wavelengths_to_analyze,mags_to_analyse,fill_value="extrapolate") #with two datapoints I can only interpolate linearly
        mag_new = interp(lambda_sdss[jj])

        #again for the magnitudes without the psf contribution
        mags_new[kk] = mag_new #before mags_new[kk] = -2.5*np.log10(flux_new)+8.9
        mags_new_error[kk] = np.sqrt( (err_mag1[kk]**2.) + (err_mag2[kk]**2.) )
        # ~ mags_new_error[kk] = err_mag1 #I take the error of the closest point in wavelength to the new obtained mag

    #to remove previous files that may have a different name
    old_data = [ff for ff in listdir(path_output) if isfile(join(path_output, ff)) == True and ff.startswith(str(gal)+"_"+field+"_sdss_"+band) == True]
    if old_data != []: 
        for old_file in old_data: 
            os.remove(path_output+old_file)
    
    data = Table( [vector_dist_arcsec, vector_dist_kpc, mags_new, mags_new_error], names = ["vector_dist_arcsec","vector_dist_kpc","sb_sdss","sb_sdss_error"] )
    data.write(path_output+str(gal)+"_"+field+"_sdss_"+band+"_profile_"+band1+band2+".cat", format = "ascii.commented_header",overwrite=True)


# In[ ]:




