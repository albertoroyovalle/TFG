#!/usr/bin/env python
# coding: utf-8

# In[27]:


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


# """
# 13-04-2021
# Para adaptar el programa es necesario:
# 
# 1) Eliminar el bucle original (for ii,gal in enumerate(gals): [1]).
# 
# 2) De ese bucle se leen datos de un catalogo. Los que necesitamos en este programa son:
# 
#     gals=np.array(["79071"])
#     gal="79071" -----------------> La razon de escribir asi gals y gal es para que no de muchos problemas el programa si el dia de mañana necesito implementar el bucle [1]
#     field="hudf" *
#     zz=0.6218
#     pa=155
#     ar=0.8
# 
# *Estrictamente solo necesitamos el dato field.(Aun asi yo introduzco el resto por si en algun caso los necesitaramos)
# 
# 3) Buscar y reemplazar [ii] (eliminando todos los terminos)---------> todos los [ii] pertenecen a field
# 
# """
# 
# """
# to run this program, be sure that the sb_profiles_B.cat (and the like) do not have nan values
# """
# 
# RECORDAR: 1Msun/kpc²=1Msun/(kpc²*1000000 pc²/Kpc²) aL REPRESENTAR

# In[28]:


path_output="./mass_profiles/"
path_sdss_profiles="./sdss_profiles/"
bands_sdss = np.array(["g","r","i","z"])
full_catalog = True #True if you want to display the mass-to-light ratios in the output catalog

gals=np.array(["ngc1277"])
gal="ngc1277"
field="hudf" ###CAMBIAR

zz=0.6218
pa=155
ar=0.8


# In[29]:


files_in_sdss_profiles_folder = [ff for ff in listdir(path_sdss_profiles) if isfile(join(path_sdss_profiles, ff))]

if not os.path.exists(path_output):
    os.mkdir(path_output)


# In[30]:


sb_profiles = { band : np.array([]) for band in bands_sdss }
err_sb_profiles = { band : np.array([]) for band in bands_sdss }


# In[31]:


for jj,band in enumerate(bands_sdss):
              
    try:
        filename = [ss for ss in files_in_sdss_profiles_folder if ss.startswith("ngc1277_hudf_sdss_correct_"+band) == True][0] #[0] for not being a list
    except:
        dist_arcsec           = np.array([float("NaN")])
        dist_kpc              = np.array([float("NaN")])
        sb_profiles[band]     = np.array([float("NaN")])
        err_sb_profiles[band] = np.array([float("NaN")])
        continue
    
    tt = Table.read(path_sdss_profiles+filename,names=["dist_arcsec","dist_kpc","sb","err_sb"],format='ascii.commented_header')
    dist_arcsec = tt["dist_arcsec"]
    dist_kpc    = tt["dist_kpc"]
    sb = tt["sb"]
    err_sb = tt["err_sb"]

    sb_profiles[band] = sb
    err_sb_profiles[band] = err_sb
    print("filenames: ",filename)
    print(sb)


# In[ ]:





# In[32]:


#obtaining the mass profiles
error_in_color_gr = np.sqrt( (err_sb_profiles["g"]**2.)+(err_sb_profiles["r"]**2.) )
error_in_color_gi = np.sqrt( (err_sb_profiles["g"]**2.)+(err_sb_profiles["i"]**2.) )
error_in_color_gz = np.sqrt( (err_sb_profiles["g"]**2.)+(err_sb_profiles["z"]**2.) )
error_in_color_ri = np.sqrt( (err_sb_profiles["r"]**2.)+(err_sb_profiles["i"]**2.) )
error_in_color_rz = np.sqrt( (err_sb_profiles["r"]**2.)+(err_sb_profiles["z"]**2.) )
error_in_color_iz = np.sqrt( (err_sb_profiles["i"]**2.)+(err_sb_profiles["z"]**2.) )

mass_profile_g_gr, error_mass_profile_g_gr, MtoL_g_gr, dMtoL_g_gr = calculate_mass_profile(sb_profiles["g"], sb_profiles["g"]-sb_profiles["r"], error_in_color_gr, "g", "g-r")
mass_profile_g_gi, error_mass_profile_g_gi, MtoL_g_gi, dMtoL_g_gi = calculate_mass_profile(sb_profiles["g"], sb_profiles["g"]-sb_profiles["i"], error_in_color_gi, "g", "g-i")
mass_profile_g_gz, error_mass_profile_g_gz, MtoL_g_gz, dMtoL_g_gz = calculate_mass_profile(sb_profiles["g"], sb_profiles["g"]-sb_profiles["z"], error_in_color_gz, "g", "g-z")
mass_profile_g_ri, error_mass_profile_g_ri, MtoL_g_ri, dMtoL_g_ri = calculate_mass_profile(sb_profiles["g"], sb_profiles["r"]-sb_profiles["i"], error_in_color_ri, "g", "r-i")
mass_profile_g_rz, error_mass_profile_g_rz, MtoL_g_rz, dMtoL_g_rz = calculate_mass_profile(sb_profiles["g"], sb_profiles["r"]-sb_profiles["z"], error_in_color_rz, "g", "r-z")
mass_profile_g_iz, error_mass_profile_g_iz, MtoL_g_iz, dMtoL_g_iz = calculate_mass_profile(sb_profiles["g"], sb_profiles["i"]-sb_profiles["z"], error_in_color_iz, "g", "i-z")

mass_profile_r_gr, error_mass_profile_r_gr, MtoL_r_gr, dMtoL_r_gr = calculate_mass_profile(sb_profiles["r"], sb_profiles["g"]-sb_profiles["r"], error_in_color_gr, "r", "g-r")
mass_profile_r_gi, error_mass_profile_r_gi, MtoL_r_gi, dMtoL_r_gi = calculate_mass_profile(sb_profiles["r"], sb_profiles["g"]-sb_profiles["i"], error_in_color_gi, "r", "g-i")
mass_profile_r_gz, error_mass_profile_r_gz, MtoL_r_gz, dMtoL_r_gz = calculate_mass_profile(sb_profiles["r"], sb_profiles["g"]-sb_profiles["z"], error_in_color_gz, "r", "g-z")
mass_profile_r_ri, error_mass_profile_r_ri, MtoL_r_ri, dMtoL_r_ri = calculate_mass_profile(sb_profiles["r"], sb_profiles["r"]-sb_profiles["i"], error_in_color_ri, "r", "r-i")
mass_profile_r_rz, error_mass_profile_r_rz, MtoL_r_rz, dMtoL_r_rz = calculate_mass_profile(sb_profiles["r"], sb_profiles["r"]-sb_profiles["z"], error_in_color_rz, "r", "r-z")
mass_profile_r_iz, error_mass_profile_r_iz, MtoL_r_iz, dMtoL_r_iz = calculate_mass_profile(sb_profiles["r"], sb_profiles["i"]-sb_profiles["z"], error_in_color_iz, "r", "i-z")

mass_profile_i_gr, error_mass_profile_i_gr, MtoL_i_gr, dMtoL_i_gr = calculate_mass_profile(sb_profiles["i"], sb_profiles["g"]-sb_profiles["r"], error_in_color_gr, "i", "g-r")
mass_profile_i_gi, error_mass_profile_i_gi, MtoL_i_gi, dMtoL_i_gi = calculate_mass_profile(sb_profiles["i"], sb_profiles["g"]-sb_profiles["i"], error_in_color_gi, "i", "g-i")
mass_profile_i_gz, error_mass_profile_i_gz, MtoL_i_gz, dMtoL_i_gz = calculate_mass_profile(sb_profiles["i"], sb_profiles["g"]-sb_profiles["z"], error_in_color_gz, "i", "g-z")
mass_profile_i_ri, error_mass_profile_i_ri, MtoL_i_ri, dMtoL_i_ri = calculate_mass_profile(sb_profiles["i"], sb_profiles["r"]-sb_profiles["i"], error_in_color_ri, "i", "r-i")
mass_profile_i_rz, error_mass_profile_i_rz, MtoL_i_rz, dMtoL_i_rz = calculate_mass_profile(sb_profiles["i"], sb_profiles["r"]-sb_profiles["z"], error_in_color_rz, "i", "r-z")
mass_profile_i_iz, error_mass_profile_i_iz, MtoL_i_iz, dMtoL_i_iz = calculate_mass_profile(sb_profiles["i"], sb_profiles["i"]-sb_profiles["z"], error_in_color_iz, "i", "i-z")

mass_profile_z_gr, error_mass_profile_z_gr, MtoL_z_gr, dMtoL_z_gr = calculate_mass_profile(sb_profiles["z"], sb_profiles["g"]-sb_profiles["r"], error_in_color_gr, "z", "g-r")
mass_profile_z_gi, error_mass_profile_z_gi, MtoL_z_gi, dMtoL_z_gi = calculate_mass_profile(sb_profiles["z"], sb_profiles["g"]-sb_profiles["i"], error_in_color_gi, "z", "g-i")
mass_profile_z_gz, error_mass_profile_z_gz, MtoL_z_gz, dMtoL_z_gz = calculate_mass_profile(sb_profiles["z"], sb_profiles["g"]-sb_profiles["z"], error_in_color_gz, "z", "g-z")
mass_profile_z_ri, error_mass_profile_z_ri, MtoL_z_ri, dMtoL_z_ri = calculate_mass_profile(sb_profiles["z"], sb_profiles["r"]-sb_profiles["i"], error_in_color_ri, "z", "r-i")
mass_profile_z_rz, error_mass_profile_z_rz, MtoL_z_rz, dMtoL_z_rz = calculate_mass_profile(sb_profiles["z"], sb_profiles["r"]-sb_profiles["z"], error_in_color_rz, "z", "r-z")
mass_profile_z_iz, error_mass_profile_z_iz, MtoL_z_iz, dMtoL_z_iz = calculate_mass_profile(sb_profiles["z"], sb_profiles["i"]-sb_profiles["z"], error_in_color_iz, "z", "i-z")

#to make everything of the same size
if np.isnan(sb_profiles["g"][0]) == True:
    mass_profile_g_gr = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_g_gi = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_g_gz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_g_ri = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_g_rz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_g_iz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_g_gr = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_g_gi = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_g_gz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_g_ri = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_g_rz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_g_iz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_g_gr = np.full_like(dist_arcsec,float("NaN"))
    MtoL_g_gi = np.full_like(dist_arcsec,float("NaN"))
    MtoL_g_gz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_g_ri = np.full_like(dist_arcsec,float("NaN"))
    MtoL_g_rz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_g_iz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_g_gr = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_g_gi = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_g_gz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_g_ri = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_g_rz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_g_iz = np.full_like(dist_arcsec,float("NaN"))
if np.isnan(sb_profiles["r"][0]) == True:
    mass_profile_r_gr = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_r_gi = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_r_gz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_r_ri = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_r_rz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_r_iz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_r_gr = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_r_gi = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_r_gz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_r_ri = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_r_rz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_r_iz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_r_gr = np.full_like(dist_arcsec,float("NaN"))
    MtoL_r_gi = np.full_like(dist_arcsec,float("NaN"))
    MtoL_r_gz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_r_ri = np.full_like(dist_arcsec,float("NaN"))
    MtoL_r_rz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_r_iz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_r_gr = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_r_gi = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_r_gz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_r_ri = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_r_rz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_r_iz = np.full_like(dist_arcsec,float("NaN"))
if np.isnan(sb_profiles["i"][0]) == True:
    mass_profile_i_gr = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_i_gi = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_i_gz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_i_ri = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_i_rz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_i_iz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_i_gr = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_i_gi = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_i_gz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_i_ri = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_i_rz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_i_iz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_i_gr = np.full_like(dist_arcsec,float("NaN"))
    MtoL_i_gi = np.full_like(dist_arcsec,float("NaN"))
    MtoL_i_gz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_i_ri = np.full_like(dist_arcsec,float("NaN"))
    MtoL_i_rz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_i_iz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_i_gr = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_i_gi = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_i_gz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_i_ri = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_i_rz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_i_iz = np.full_like(dist_arcsec,float("NaN"))
if np.isnan(sb_profiles["z"][0]) == True:
    mass_profile_z_gr = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_z_gi = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_z_gz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_z_ri = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_z_rz = np.full_like(dist_arcsec,float("NaN"))
    mass_profile_z_iz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_z_gr = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_z_gi = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_z_gz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_z_ri = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_z_rz = np.full_like(dist_arcsec,float("NaN"))
    error_mass_profile_z_iz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_z_gr = np.full_like(dist_arcsec,float("NaN"))
    MtoL_z_gi = np.full_like(dist_arcsec,float("NaN"))
    MtoL_z_gz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_z_ri = np.full_like(dist_arcsec,float("NaN"))
    MtoL_z_rz = np.full_like(dist_arcsec,float("NaN"))
    MtoL_z_iz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_z_gr = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_z_gi = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_z_gz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_z_ri = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_z_rz = np.full_like(dist_arcsec,float("NaN"))
    dMtoL_z_iz = np.full_like(dist_arcsec,float("NaN"))
    
#I write the mass profiles in a text catalog
if full_catalog == True:
    data = Table( [dist_arcsec, dist_kpc, 
                   mass_profile_g_gr, error_mass_profile_g_gr, mass_profile_g_gi, error_mass_profile_g_gi, mass_profile_g_gz, error_mass_profile_g_gz, 
                   mass_profile_g_ri, error_mass_profile_g_ri, mass_profile_g_rz, error_mass_profile_g_rz, mass_profile_g_iz, error_mass_profile_g_iz, 
                   mass_profile_r_gr, error_mass_profile_r_gr, mass_profile_r_gi, error_mass_profile_r_gi, mass_profile_r_gz, error_mass_profile_r_gz, 
                   mass_profile_r_ri, error_mass_profile_r_ri, mass_profile_r_rz, error_mass_profile_r_rz, mass_profile_r_iz, error_mass_profile_r_iz, 
                   mass_profile_i_gr, error_mass_profile_i_gr, mass_profile_i_gi, error_mass_profile_i_gi, mass_profile_i_gz, error_mass_profile_i_gz, 
                   mass_profile_i_ri, error_mass_profile_i_ri, mass_profile_i_rz, error_mass_profile_i_rz, mass_profile_i_iz, error_mass_profile_i_iz, 
                   mass_profile_z_gr, error_mass_profile_z_gr, mass_profile_z_gi, error_mass_profile_z_gi, mass_profile_z_gz, error_mass_profile_z_gz, 
                   mass_profile_z_ri, error_mass_profile_z_ri, mass_profile_z_rz, error_mass_profile_z_rz, mass_profile_z_iz, error_mass_profile_z_iz,
                   MtoL_g_gr, dMtoL_g_gr, MtoL_g_gi, dMtoL_g_gi, MtoL_g_gz, dMtoL_g_gz, 
                   MtoL_g_ri, dMtoL_g_ri, MtoL_g_rz, dMtoL_g_rz, MtoL_g_iz, dMtoL_g_iz, 
                   MtoL_r_gr, dMtoL_r_gr, MtoL_r_gi, dMtoL_r_gi, MtoL_r_gz, dMtoL_r_gz, 
                   MtoL_r_ri, dMtoL_r_ri, MtoL_r_rz, dMtoL_r_rz, MtoL_r_iz, dMtoL_r_iz, 
                   MtoL_i_gr, dMtoL_i_gr, MtoL_i_gi, dMtoL_i_gi, MtoL_i_gz, dMtoL_i_gz, 
                   MtoL_i_ri, dMtoL_i_ri, MtoL_i_rz, dMtoL_i_rz, MtoL_i_iz, dMtoL_i_iz, 
                   MtoL_z_gr, dMtoL_z_gr, MtoL_z_gi, dMtoL_z_gi, MtoL_z_gz, dMtoL_z_gz, 
                   MtoL_z_ri, dMtoL_z_ri, MtoL_z_rz, dMtoL_z_rz, MtoL_z_iz, dMtoL_z_iz],                                                          
                   names = ["dist_arcsec","dist_kpc",
                            "mass_profile_g_gr","err_mass_profile_g_gr","mass_profile_g_gi","err_mass_profile_g_gi","mass_profile_g_gz","err_mass_profile_g_gz",
                            "mass_profile_g_ri","err_mass_profile_g_ri","mass_profile_g_rz","err_mass_profile_g_rz","mass_profile_g_iz","err_mass_profile_g_iz",
                            "mass_profile_r_gr","err_mass_profile_r_gr","mass_profile_r_gi","err_mass_profile_r_gi","mass_profile_r_gz","err_mass_profile_r_gz",
                            "mass_profile_r_ri","err_mass_profile_r_ri","mass_profile_r_rz","err_mass_profile_r_rz","mass_profile_r_iz","err_mass_profile_r_iz",
                            "mass_profile_i_gr","err_mass_profile_i_gr","mass_profile_i_gi","err_mass_profile_i_gi","mass_profile_i_gz","err_mass_profile_i_gz",
                            "mass_profile_i_ri","err_mass_profile_i_ri","mass_profile_i_rz","err_mass_profile_i_rz","mass_profile_i_iz","err_mass_profile_i_iz",
                            "mass_profile_z_gr","err_mass_profile_z_gr","mass_profile_z_gi","err_mass_profile_z_gi","mass_profile_z_gz","err_mass_profile_z_gz",
                            "mass_profile_z_ri","err_mass_profile_z_ri","mass_profile_z_rz","err_mass_profile_z_rz","mass_profile_z_iz","err_mass_profile_z_iz",
                            "MtoL_g_gr","dMtoL_g_gr","MtoL_g_gi","dMtoL_g_gi","MtoL_g_gz","dMtoL_g_gz",
                            "MtoL_g_ri","dMtoL_g_ri","MtoL_g_rz","dMtoL_g_rz","MtoL_g_iz","dMtoL_g_iz",
                            "MtoL_r_gr","dMtoL_r_gr","MtoL_r_gi","dMtoL_r_gi","MtoL_r_gz","dMtoL_r_gz",
                            "MtoL_r_ri","dMtoL_r_ri","MtoL_r_rz","dMtoL_r_rz","MtoL_r_iz","dMtoL_r_iz",
                            "MtoL_i_gr","dMtoL_i_gr","MtoL_i_gi","dMtoL_i_gi","MtoL_i_gz","dMtoL_i_gz",
                            "MtoL_i_ri","dMtoL_i_ri","MtoL_i_rz","dMtoL_i_rz","MtoL_i_iz","dMtoL_i_iz",
                            "MtoL_z_gr","dMtoL_z_gr","MtoL_z_gi","dMtoL_z_gi","MtoL_z_gz","dMtoL_z_gz",
                            "MtoL_z_ri","dMtoL_z_ri","MtoL_z_rz","dMtoL_z_rz","MtoL_z_iz","dMtoL_z_iz"] )
    data.write(path_output+"mass_profiles_full_"+gal+"_"+field+".cat", format = "ascii.commented_header",overwrite=True)
else:
    data = Table( [dist_arcsec, dist_kpc, 
                   mass_profile_g_gr, error_mass_profile_g_gr, mass_profile_g_gi, error_mass_profile_g_gi, mass_profile_g_gz, error_mass_profile_g_gz, 
                   mass_profile_g_ri, error_mass_profile_g_ri, mass_profile_g_rz, error_mass_profile_g_rz, mass_profile_g_iz, error_mass_profile_g_iz, 
                   mass_profile_r_gr, error_mass_profile_r_gr, mass_profile_r_gi, error_mass_profile_r_gi, mass_profile_r_gz, error_mass_profile_r_gz, 
                   mass_profile_r_ri, error_mass_profile_r_ri, mass_profile_r_rz, error_mass_profile_r_rz, mass_profile_r_iz, error_mass_profile_r_iz, 
                   mass_profile_i_gr, error_mass_profile_i_gr, mass_profile_i_gi, error_mass_profile_i_gi, mass_profile_i_gz, error_mass_profile_i_gz, 
                   mass_profile_i_ri, error_mass_profile_i_ri, mass_profile_i_rz, error_mass_profile_i_rz, mass_profile_i_iz, error_mass_profile_i_iz, 
                   mass_profile_z_gr, error_mass_profile_z_gr, mass_profile_z_gi, error_mass_profile_z_gi, mass_profile_z_gz, error_mass_profile_z_gz, 
                   mass_profile_z_ri, error_mass_profile_z_ri, mass_profile_z_rz, error_mass_profile_z_rz, mass_profile_z_iz, error_mass_profile_z_iz],                                                          
                   names = ["dist_arcsec","dist_kpc",
                            "mass_profile_g_gr","err_mass_profile_g_gr","mass_profile_g_gi","err_mass_profile_g_gi","mass_profile_g_gz","err_mass_profile_g_gz",
                            "mass_profile_g_ri","err_mass_profile_g_ri","mass_profile_g_rz","err_mass_profile_g_rz","mass_profile_g_iz","err_mass_profile_g_iz",
                            "mass_profile_r_gr","err_mass_profile_r_gr","mass_profile_r_gi","err_mass_profile_r_gi","mass_profile_r_gz","err_mass_profile_r_gz",
                            "mass_profile_r_ri","err_mass_profile_r_ri","mass_profile_r_rz","err_mass_profile_r_rz","mass_profile_r_iz","err_mass_profile_r_iz",
                            "mass_profile_i_gr","err_mass_profile_i_gr","mass_profile_i_gi","err_mass_profile_i_gi","mass_profile_i_gz","err_mass_profile_i_gz",
                            "mass_profile_i_ri","err_mass_profile_i_ri","mass_profile_i_rz","err_mass_profile_i_rz","mass_profile_i_iz","err_mass_profile_i_iz",
                            "mass_profile_z_gr","err_mass_profile_z_gr","mass_profile_z_gi","err_mass_profile_z_gi","mass_profile_z_gz","err_mass_profile_z_gz",
                            "mass_profile_z_ri","err_mass_profile_z_ri","mass_profile_z_rz","err_mass_profile_z_rz","mass_profile_z_iz","err_mass_profile_z_iz"] )
    data.write(path_output+"mass_profiles_"+gal+"_"+field+".cat", format = "ascii.commented_header",overwrite=True)


# In[ ]:




