#!/usr/bin/env python
# coding: utf-8

# PROGRAMA PARA LAS CORRECCIONES COSMOLÓGICAS PARA LAS GALAXIAS EN LOS FILTROS SDSS:
# 
# Nota: Para otras galaxias en filtros Hudf (hubble) la correccion esta implementada en el codigo getting_sdss_profiles.py

# In[22]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib import rc,rcParams
import matplotlib as mpl
from os import listdir
from os.path import isfile, join
import os
import pdb
import sys
from functs_for_truncations import obtaining_best_redshifts
from my_python_library_v2 import calculate_mass_profile
from my_python_library_v2 import correcting_by_extinction
from my_python_library_v2 import correcting_by_inclination


# In[23]:


path_output="./sdss_profiles/"
path_sdss_profiles="./sdss_profiles/"
bands_obs = np.array(["g","r","i","z"])
full_catalog = True #True if you want to display the mass-to-light ratios in the output catalog

gals=np.array(["ngc1277"])
gal="ngc1277"
field="kids" ###CAMBIAR

zz=0.01690
pa=57.73
ar=1.9047619047619047


# In[24]:


files_in_sdss_profiles_folder = [ff for ff in listdir(path_sdss_profiles) if isfile(join(path_sdss_profiles, ff))]

if not os.path.exists(path_output):
    os.mkdir(path_output)


# In[25]:


mags     = { band : np.array([]) for band in bands_obs }
mags_err = { band : np.array([]) for band in bands_obs }


# In[26]:


#I load in memory all the surface brightness profiles
for jj,band in enumerate(bands_obs):
    
    #gal_sb_cat = str(gal)+"_"+field+"_"+band+".cat" Cambié esto ya que la entrada de mis funciones era distinta
    gal_sb_cat="data_"+band+".cat"
    if os.path.exists("./sdss_profiles/"+"ngc1277_hudf_sdss_"+band+".cat"):
        tt = Table.read("./sdss_profiles/"+"ngc1277_hudf_sdss_"+band+".cat",names=["dist_arcsec","dist_kpc","sb","sb_err"],format="ascii.commented_header")
        vector_dist_arcsec = np.array(tt["dist_arcsec"])
        vector_dist_kpc    = np.array(tt["dist_kpc"])
        mags[band]         = np.array(tt["sb"])
        mags_err[band]     = np.array(tt["sb_err"])


print(mags)


# In[27]:


for jj,band in enumerate(bands_obs):
    
    mags[band] = correcting_by_extinction(mags[band],field,band,False)
    """
    #https://ned.ipac.caltech.edu/extinction_calculator
    
    A esta página la gusta que metas las coordenadas del tipo (03h19m51.49s, +41d34m24.7s)
    RA:  h m s  ====>  Horas,          Minutos          Segundos 
    DEC: d m s  ====>  grados (° o d), minutes (′ o m), and seconds (″ o s)
    
    
    Ascension recta  AR->   15º=1hora----------> 360º=24h (ojo, NO 12h)
    Declinación      DEC->
    
    Hacer esto para las galaxias SDSS:
    
    79010
    
    ngc1277
    

    """
    
    
    
    mags[band] = correcting_by_inclination(mags[band], ar)
    
    
    print(band,mags[band])


# In[28]:


mags_corr = { band: np.full(vector_dist_arcsec.size,float("NaN")) for band in bands_obs }
mags_err2 = { band: np.full(vector_dist_arcsec.size,float("NaN")) for band in bands_obs }
    


# In[29]:


for jj,band in enumerate(bands_obs):
    mags_corr[band] = mags[band] - 7.5*np.log10(1.+zz) #correcting by cosmological dimming a la Giavalisco+04 and Ribeiro+16 --> and as you already use lambda restframe, you need to use (1+z)^4
    #I don't work in flux as before because the galaxy SED is in a logarithmic space (magnitudes) and there it is the place where I should approximate by linear relations
    #before flux_in_Jy_restframe = 10.**( (mags_corr-8.9)/(-2.5) ) #from mag=-2.5log(flux_in_Jy)-2.5log(10^-23)-48.6 Oke & Gunn 1983 ;I will do the extrapolation in flux and not in magnitude just in case there was any problem for interpolating logarithms
    mags_err2[band] = mags_err[band]


# In[30]:


mags_corr["g"]


# In[33]:


bands_sdss = np.array(["g","r","i","z"])
for jj,band in enumerate(bands_sdss):
    data = Table( [vector_dist_arcsec, vector_dist_kpc, mags_corr[band], mags_err2[band]], names = ["vector_dist_arcsec","vector_dist_kpc","sb","sb_error"] )
    data.write(path_output+"ngc1277_hudf_sdss_correct_"+band+".cat", format = "ascii.commented_header",overwrite=True)


# In[ ]:




