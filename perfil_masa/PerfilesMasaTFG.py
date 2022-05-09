#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import pdb
from astropy.cosmology import FlatLambdaCDM
from matplotlib import rc, rcParams
import matplotlib as mpl


# In[2]:


def plot_options():
    #starting the plotting
    rc('axes', linewidth=5)
    
    rc('font', weight='bold')
    
    rc('font', size=17)
    
    
    # set tick width
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['xtick.major.width'] = 4
    mpl.rcParams['xtick.minor.size'] = 5
    mpl.rcParams['xtick.minor.width'] = 2
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['ytick.major.width'] = 4
    mpl.rcParams['ytick.minor.size'] = 5
    mpl.rcParams['ytick.minor.width'] = 2
    mpl.rcParams['xtick.direction'] = "in"
    mpl.rcParams['ytick.direction'] = "in"
    mpl.rcParams['xtick.major.pad'] = 7 #padding between the x-axis and the numbers corresponding to the major ticks
    #for more params mpl.rcParams.keys()
    


# 

# In[3]:


tt = Table.read("mass_profiles_full_ngc1277_hudf.cat",names=["dist_arcsec","dist_kpc",
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
                            "MtoL_z_ri","dMtoL_z_ri","MtoL_z_rz","dMtoL_z_rz","MtoL_z_iz","dMtoL_z_iz"],format="ascii.commented_header")


# In[4]:


names=np.array(["dist_arcsec","dist_kpc",
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
                            "MtoL_z_ri","dMtoL_z_ri","MtoL_z_rz","dMtoL_z_rz","MtoL_z_iz","dMtoL_z_iz"])


# In[5]:




zz=0.01690
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(zz)
print("Kpc/arcsec: ", kpc_per_arcsec)
print("Arcsec/kpc: ", 1/kpc_per_arcsec)


# In[6]:


c=0
for i,a in enumerate(names):
    globals()[names[i]]=np.array(tt[names[i]])
   
    if a.startswith("mass_profile_"):
        print(names[i])
        c=c+1
        
print(c)


# In[7]:


band=np.array(["g","r","i","z"])

legends_1=np.array(["g_gr","g_gi","g_gz","g_ri","g_rz","g_iz"])
colors=np.array(["pink","blue","green","yellow","brown","red"])


# In[8]:


print("Mediana=massprofile/número de masses profiles")
print("c: ",c)
mediana=(mass_profile_g_gr+mass_profile_g_gi+
         mass_profile_g_gz+mass_profile_g_ri+
         mass_profile_g_rz+mass_profile_g_iz+
         mass_profile_r_gr+mass_profile_r_gi+
         mass_profile_r_gz+mass_profile_r_ri+
         mass_profile_r_rz+mass_profile_r_iz+
         mass_profile_i_gr+mass_profile_i_gi+
         mass_profile_i_gz+mass_profile_i_ri+
         mass_profile_i_rz+mass_profile_i_iz+
         mass_profile_z_gr+mass_profile_z_gi+
         mass_profile_z_gz+mass_profile_z_ri+
         mass_profile_z_rz+mass_profile_z_iz)/c
print("Mediana: ", mediana)


# In[12]:


fig = plt.figure(figsize=(22,22))
plot_options()
ax = fig.add_subplot(2,2,1)
ax.scatter(dist_arcsec,(MtoL_g_gr),color=colors[0],label="gr")
ax.scatter(dist_arcsec,(MtoL_g_gi),color=colors[1],label="gi")
ax.scatter(dist_arcsec,(MtoL_g_gz),color=colors[2],label="gz")
ax.scatter(dist_arcsec,(MtoL_g_ri),color=colors[3],label="ri")
ax.scatter(dist_arcsec,(MtoL_g_rz),color=colors[4],label="rz")
ax.scatter(dist_arcsec,(MtoL_g_iz),color=colors[5],label="iz")

ax.errorbar(dist_arcsec,(MtoL_g_gr), color=colors[0])
ax.errorbar(dist_arcsec,(MtoL_g_gi), color=colors[1])
ax.errorbar(dist_arcsec,(MtoL_g_gz), color=colors[2])
ax.errorbar(dist_arcsec,(MtoL_g_ri), color=colors[3])
ax.errorbar(dist_arcsec,(MtoL_g_rz), color=colors[4])
ax.errorbar(dist_arcsec,(MtoL_g_iz), color=colors[5])



ax.set_xlim(0,50)
ax.set_ylim(0,14)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_g_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("dist_kpc")
plt.grid()

ax.legend()  
ax.set_ylabel("M/L ratio")
ax.set_xlabel("dist_arcsec")
#plt.yscale('log')
plt.title("Perfil M/L g",fontweight="bold")
#plt.savefig("Perfil_Masa_g.pdf")
ax = fig.add_subplot(2,2,2)
ax.scatter(dist_arcsec,(MtoL_r_gr),color=colors[0],label="gr")
ax.scatter(dist_arcsec,(MtoL_r_gi),color=colors[1],label="gi")
ax.scatter(dist_arcsec,(MtoL_r_gz),color=colors[2],label="gz")
ax.scatter(dist_arcsec,(MtoL_r_ri),color=colors[3],label="ri")
ax.scatter(dist_arcsec,(MtoL_r_rz),color=colors[4],label="rz")
ax.scatter(dist_arcsec,(MtoL_r_iz),color=colors[5],label="iz")

ax.errorbar(dist_arcsec,(MtoL_r_gr), color=colors[0])
ax.errorbar(dist_arcsec,(MtoL_r_gi), color=colors[1])
ax.errorbar(dist_arcsec,(MtoL_r_gz), color=colors[2])
ax.errorbar(dist_arcsec,(MtoL_r_ri), color=colors[3])
ax.errorbar(dist_arcsec,(MtoL_r_rz), color=colors[4])
ax.errorbar(dist_arcsec,(MtoL_r_iz), color=colors[5])



ax.set_xlim(0,50)
ax.set_ylim(0,14)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_g_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("dist_kpc")
plt.grid()

ax.legend()  
ax.set_ylabel("M/L ratio")
ax.set_xlabel("dist_arcsec")
#plt.yscale('log')
plt.title("Perfil M/L r",fontweight="bold")


#Representaciones principales
ax = fig.add_subplot(2,2,3)
ax.scatter(dist_arcsec,(MtoL_i_gr),color=colors[0],label="gr")
ax.scatter(dist_arcsec,(MtoL_i_gi),color=colors[1],label="gi")
ax.scatter(dist_arcsec,(MtoL_i_gz),color=colors[2],label="gz")
ax.scatter(dist_arcsec,(MtoL_i_ri),color=colors[3],label="ri")
ax.scatter(dist_arcsec,(MtoL_i_rz),color=colors[4],label="rz")
ax.scatter(dist_arcsec,(MtoL_i_iz),color=colors[5],label="iz")

ax.errorbar(dist_arcsec,(MtoL_i_gr), color=colors[0])
ax.errorbar(dist_arcsec,(MtoL_i_gi), color=colors[1])
ax.errorbar(dist_arcsec,(MtoL_i_gz), color=colors[2])
ax.errorbar(dist_arcsec,(MtoL_i_ri), color=colors[3])
ax.errorbar(dist_arcsec,(MtoL_i_rz), color=colors[4])
ax.errorbar(dist_arcsec,(MtoL_i_iz), color=colors[5])



ax.set_xlim(0,50)
ax.set_ylim(0,14)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_g_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("dist_arcsec")
plt.grid()

ax.legend()  
ax.set_ylabel("M/L ratio")
ax.set_xlabel("dist_kpc")
#plt.yscale('log')
plt.title("Perfil M/L i",fontweight="bold")
ax = fig.add_subplot(2,2,4)

ax.scatter(dist_arcsec,(MtoL_z_gr),color=colors[0],label="gr")
ax.scatter(dist_arcsec,(MtoL_z_gi),color=colors[1],label="gi")
ax.scatter(dist_arcsec,(MtoL_z_gz),color=colors[2],label="gz")
ax.scatter(dist_arcsec,(MtoL_z_ri),color=colors[3],label="ri")
ax.scatter(dist_arcsec,(MtoL_z_rz),color=colors[4],label="rz")
ax.scatter(dist_arcsec,(MtoL_z_iz),color=colors[5],label="iz")

ax.errorbar(dist_arcsec,(MtoL_z_gr), color=colors[0])
ax.errorbar(dist_arcsec,(MtoL_z_gi), color=colors[1])
ax.errorbar(dist_arcsec,(MtoL_z_gz), color=colors[2])
ax.errorbar(dist_arcsec,(MtoL_z_ri), color=colors[3])
ax.errorbar(dist_arcsec,(MtoL_z_rz), color=colors[4])
ax.errorbar(dist_arcsec,(MtoL_z_iz), color=colors[5])



ax.set_xlim(0,50)
ax.set_ylim(0,14)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_g_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("dist_kpc")
plt.grid()

ax.legend()  
ax.set_ylabel("M/L ratio")
ax.set_xlabel("dist_arcsec")
#plt.yscale('log')
plt.title("Perfil M/L z",fontweight="bold")



plt.savefig("Perfiles_Mass_to_light_ratio_ngc1277.pdf")


# In[1]:


print("ax.get_xlim()[1]*kpc_per_arcsec.value",ax.get_xlim()[1]*kpc_per_arcsec.value)


# Hacer mediana.

# In[64]:


fig = plt.figure(figsize=(25,25))
plot_options()
ax = fig.add_subplot(2,2,1)
ax.scatter(dist_arcsec,(mass_profile_g_gr/1e6),color=colors[0],label="mass_profile_g_gr")
ax.scatter(dist_arcsec,(mass_profile_g_gi/1e6),color=colors[1],label="mass_profile_g_gi")
ax.scatter(dist_arcsec,(mass_profile_g_gz/1e6),color=colors[2],label="mass_profile_g_gz")
#ax.scatter(dist_arcsec,(mass_profile_g_ri/1e6),color=colors[3],label="mass_profile_g_ri")
ax.scatter(dist_arcsec,(mass_profile_g_rz/1e6),color=colors[4],label="mass_profile_g_rz")
#ax.scatter(dist_arcsec,(mass_profile_g_iz/1e6),color=colors[5],label="mass_profile_g_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)

ax.errorbar(dist_arcsec,(mass_profile_g_gr/1e6), yerr = (err_mass_profile_g_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_g_gi/1e6), yerr = (err_mass_profile_g_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_g_gz/1e6), yerr = (err_mass_profile_g_gz/1e6), color=colors[2])
#ax.errorbar(dist_arcsec,(mass_profile_g_ri/1e6), yerr = (err_mass_profile_g_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_g_rz/1e6), yerr = (err_mass_profile_g_rz/1e6), color=colors[4])
#ax.errorbar(dist_arcsec,(mass_profile_g_iz/1e6), yerr = (err_mass_profile_g_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")


ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_g_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa g",fontweight="bold")
#plt.savefig("Perfil_Masa_g.pdf")
ax = fig.add_subplot(2,2,2)
ax.scatter(dist_arcsec,(mass_profile_r_gr/1e6),color=colors[0],label="mass_profile_r_gr")
ax.scatter(dist_arcsec,(mass_profile_r_gi/1e6),color=colors[1],label="mass_profile_r_gi")
ax.scatter(dist_arcsec,(mass_profile_r_gz/1e6),color=colors[2],label="mass_profile_r_gz")
#ax.scatter(dist_arcsec,(mass_profile_r_ri/1e6),color=colors[3],label="mass_profile_r_ri")
ax.scatter(dist_arcsec,(mass_profile_r_rz/1e6),color=colors[4],label="mass_profile_r_rz")
#ax.scatter(dist_arcsec,(mass_profile_r_iz/1e6),color=colors[5],label="mass_profile_r_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)

ax.errorbar(dist_arcsec,(mass_profile_r_gr/1e6), yerr = (err_mass_profile_r_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_r_gi/1e6), yerr = (err_mass_profile_r_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_r_gz/1e6), yerr = (err_mass_profile_r_gz/1e6), color=colors[2])
#ax.errorbar(dist_arcsec,(mass_profile_r_ri/1e6), yerr = (err_mass_profile_r_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_r_rz/1e6), yerr = (err_mass_profile_r_rz/1e6), color=colors[4])
#ax.errorbar(dist_arcsec,(mass_profile_r_iz/1e6), yerr = (err_mass_profile_r_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")

ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_r_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa r",fontweight="bold")
#Representaciones principales
ax = fig.add_subplot(2,2,3)
ax.scatter(dist_arcsec,(mass_profile_i_gr/1e6),color=colors[0],label="mass_profile_i_gr")
ax.scatter(dist_arcsec,(mass_profile_i_gi/1e6),color=colors[1],label="mass_profile_i_gi")
ax.scatter(dist_arcsec,(mass_profile_i_gz/1e6),color=colors[2],label="mass_profile_i_gz")
#ax.scatter(dist_arcsec,(mass_profile_i_ri/1e6),color=colors[3],label="mass_profile_i_ri")
ax.scatter(dist_arcsec,(mass_profile_i_rz/1e6),color=colors[4],label="mass_profile_i_rz")
#ax.scatter(dist_arcsec,(mass_profile_i_iz/1e6),color=colors[5],label="mass_profile_i_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)
#Errores
ax.errorbar(dist_arcsec,(mass_profile_i_gr/1e6), yerr = (err_mass_profile_i_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_i_gi/1e6), yerr = (err_mass_profile_i_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_i_gz/1e6), yerr = (err_mass_profile_i_gz/1e6), color=colors[2])
#ax.errorbar(dist_arcsec,(mass_profile_i_ri/1e6), yerr = (err_mass_profile_i_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_i_rz/1e6), yerr = (err_mass_profile_i_rz/1e6), color=colors[4])
#ax.errorbar(dist_arcsec,(mass_profile_i_iz/1e6), yerr = (err_mass_profile_i_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")
#Establecer limites en los ejes
ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_i_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa i",fontweight="bold")
#
ax = fig.add_subplot(2,2,4)

ax.scatter(dist_arcsec,(mass_profile_z_gr/1e6),color=colors[0],label="mass_profile_z_gr")
ax.scatter(dist_arcsec,(mass_profile_z_gi/1e6),color=colors[1],label="mass_profile_z_gi")
ax.scatter(dist_arcsec,(mass_profile_z_gz/1e6),color=colors[2],label="mass_profile_z_gz")
#ax.scatter(dist_arcsec,(mass_profile_z_ri/1e6),color=colors[3],label="mass_profile_z_ri")
ax.scatter(dist_arcsec,(mass_profile_z_rz/1e6),color=colors[4],label="mass_profile_z_rz")
#ax.scatter(dist_arcsec,(mass_profile_z_iz/1e6),color=colors[5],label="mass_profile_z_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)

ax.errorbar(dist_arcsec,(mass_profile_z_gr/1e6), yerr = (err_mass_profile_z_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_z_gi/1e6), yerr = (err_mass_profile_z_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_z_gz/1e6), yerr = (err_mass_profile_z_gz/1e6), color=colors[2])
#ax.errorbar(dist_arcsec,(mass_profile_z_ri/1e6), yerr = (err_mass_profile_z_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_z_rz/1e6), yerr = (err_mass_profile_z_rz/1e6), color=colors[4])
#ax.errorbar(dist_arcsec,(mass_profile_z_iz/1e6), yerr = (err_mass_profile_z_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")

ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_z_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa z",fontweight="bold")



plt.savefig("Perfiles_Masa_ngc1277_Mejorado.pdf")


# In[65]:


fig = plt.figure(figsize=(25,25))
plot_options()
ax = fig.add_subplot(2,2,1)
ax.scatter(dist_arcsec,(mass_profile_g_gr/1e6),color=colors[0],label="mass_profile_g_gr")
ax.scatter(dist_arcsec,(mass_profile_g_gi/1e6),color=colors[1],label="mass_profile_g_gi")
ax.scatter(dist_arcsec,(mass_profile_g_gz/1e6),color=colors[2],label="mass_profile_g_gz")
ax.scatter(dist_arcsec,(mass_profile_g_ri/1e6),color=colors[3],label="mass_profile_g_ri")
ax.scatter(dist_arcsec,(mass_profile_g_rz/1e6),color=colors[4],label="mass_profile_g_rz")
ax.scatter(dist_arcsec,(mass_profile_g_iz/1e6),color=colors[5],label="mass_profile_g_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)

ax.errorbar(dist_arcsec,(mass_profile_g_gr/1e6), yerr = (err_mass_profile_g_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_g_gi/1e6), yerr = (err_mass_profile_g_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_g_gz/1e6), yerr = (err_mass_profile_g_gz/1e6), color=colors[2])
ax.errorbar(dist_arcsec,(mass_profile_g_ri/1e6), yerr = (err_mass_profile_g_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_g_rz/1e6), yerr = (err_mass_profile_g_rz/1e6), color=colors[4])
ax.errorbar(dist_arcsec,(mass_profile_g_iz/1e6), yerr = (err_mass_profile_g_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")


ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_g_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa g",fontweight="bold")
#plt.savefig("Perfil_Masa_g.pdf")
ax = fig.add_subplot(2,2,2)
ax.scatter(dist_arcsec,(mass_profile_r_gr/1e6),color=colors[0],label="mass_profile_r_gr")
ax.scatter(dist_arcsec,(mass_profile_r_gi/1e6),color=colors[1],label="mass_profile_r_gi")
ax.scatter(dist_arcsec,(mass_profile_r_gz/1e6),color=colors[2],label="mass_profile_r_gz")
ax.scatter(dist_arcsec,(mass_profile_r_ri/1e6),color=colors[3],label="mass_profile_r_ri")
ax.scatter(dist_arcsec,(mass_profile_r_rz/1e6),color=colors[4],label="mass_profile_r_rz")
ax.scatter(dist_arcsec,(mass_profile_r_iz/1e6),color=colors[5],label="mass_profile_r_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)

ax.errorbar(dist_arcsec,(mass_profile_r_gr/1e6), yerr = (err_mass_profile_r_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_r_gi/1e6), yerr = (err_mass_profile_r_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_r_gz/1e6), yerr = (err_mass_profile_r_gz/1e6), color=colors[2])
ax.errorbar(dist_arcsec,(mass_profile_r_ri/1e6), yerr = (err_mass_profile_r_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_r_rz/1e6), yerr = (err_mass_profile_r_rz/1e6), color=colors[4])
ax.errorbar(dist_arcsec,(mass_profile_r_iz/1e6), yerr = (err_mass_profile_r_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")

ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_r_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa r",fontweight="bold")
#Representaciones principales
ax = fig.add_subplot(2,2,3)
ax.scatter(dist_arcsec,(mass_profile_i_gr/1e6),color=colors[0],label="mass_profile_i_gr")
ax.scatter(dist_arcsec,(mass_profile_i_gi/1e6),color=colors[1],label="mass_profile_i_gi")
ax.scatter(dist_arcsec,(mass_profile_i_gz/1e6),color=colors[2],label="mass_profile_i_gz")
ax.scatter(dist_arcsec,(mass_profile_i_ri/1e6),color=colors[3],label="mass_profile_i_ri")
ax.scatter(dist_arcsec,(mass_profile_i_rz/1e6),color=colors[4],label="mass_profile_i_rz")
ax.scatter(dist_arcsec,(mass_profile_i_iz/1e6),color=colors[5],label="mass_profile_i_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)
#Errores
ax.errorbar(dist_arcsec,(mass_profile_i_gr/1e6), yerr = (err_mass_profile_i_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_i_gi/1e6), yerr = (err_mass_profile_i_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_i_gz/1e6), yerr = (err_mass_profile_i_gz/1e6), color=colors[2])
ax.errorbar(dist_arcsec,(mass_profile_i_ri/1e6), yerr = (err_mass_profile_i_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_i_rz/1e6), yerr = (err_mass_profile_i_rz/1e6), color=colors[4])
ax.errorbar(dist_arcsec,(mass_profile_i_iz/1e6), yerr = (err_mass_profile_i_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")
#Establecer limites en los ejes
ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_i_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa i",fontweight="bold")
#
ax = fig.add_subplot(2,2,4)

ax.scatter(dist_arcsec,(mass_profile_z_gr/1e6),color=colors[0],label="mass_profile_z_gr")
ax.scatter(dist_arcsec,(mass_profile_z_gi/1e6),color=colors[1],label="mass_profile_z_gi")
ax.scatter(dist_arcsec,(mass_profile_z_gz/1e6),color=colors[2],label="mass_profile_z_gz")
ax.scatter(dist_arcsec,(mass_profile_z_ri/1e6),color=colors[3],label="mass_profile_z_ri")
ax.scatter(dist_arcsec,(mass_profile_z_rz/1e6),color=colors[4],label="mass_profile_z_rz")
ax.scatter(dist_arcsec,(mass_profile_z_iz/1e6),color=colors[5],label="mass_profile_z_iz")
ax.scatter(dist_arcsec,mediana/1e6,color="black",label="Mediana del perfil de masa",linewidths=2)

ax.errorbar(dist_arcsec,(mass_profile_z_gr/1e6), yerr = (err_mass_profile_z_gr/1e6), color=colors[0])
ax.errorbar(dist_arcsec,(mass_profile_z_gi/1e6), yerr = (err_mass_profile_z_gi/1e6), color=colors[1])
ax.errorbar(dist_arcsec,(mass_profile_z_gz/1e6), yerr = (err_mass_profile_z_gz/1e6), color=colors[2])
ax.errorbar(dist_arcsec,(mass_profile_z_ri/1e6), yerr = (err_mass_profile_z_ri/1e6), color=colors[3])
ax.errorbar(dist_arcsec,(mass_profile_z_rz/1e6), yerr = (err_mass_profile_z_rz/1e6), color=colors[4])
ax.errorbar(dist_arcsec,(mass_profile_z_iz/1e6), yerr = (err_mass_profile_z_iz/1e6), color=colors[5])
ax.errorbar(dist_arcsec,mediana/1e6,color="black")

ax.set_xlim(0,50)
ax.set_ylim(1e-2,1e3)
plt.grid()

#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc,mass_profile_z_gr/1e6, alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

ax.legend()  
ax.set_ylabel("Densidad de Masa [M☉/pc]")
plt.yscale('log')
ax.set_xlabel("Distancia [arcsec]")
plt.title("Perfil de Masa z",fontweight="bold")



plt.savefig("Perfiles_Masa_ngc1277.pdf")


# In[ ]:





# In[ ]:





# In[ ]:




