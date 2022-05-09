#!/usr/bin/env python
# coding: utf-8

# 

# In[38]:


import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import pdb
from astropy.cosmology import FlatLambdaCDM
from matplotlib import rc, rcParams
import matplotlib as mpl


# In[39]:


def plot_options():
    #starting the plotting
    rc('axes', linewidth=5)
    
    rc('font', weight='bold')
    
    rc('font', size=20)
    
    
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
    


# In[40]:


def creating_galactocentric_distance_vector_in_kpc(step,end,transition=None,step2=None):
    """
    #step: distance step in kpc
    #end: end of distance vector in kpc
    #transition: [optional] transition distance in kpc for the two regimes in spatial resolution of the distance vector
    #step2: [optional] distance step in kpc for the second part of the distance vector
    """

    if transition == None:
        x0_1 = step
        x1_1 = end
        distance_in_kpc  = np.arange(x0_1,x1_1,step)
    else:
        x0_1 = step
        x1_1 = transition
        first_part_vect_dist  = np.arange(x0_1,x1_1,step)
        x0_2 = x1_1
        x1_2 = end
        second_part_vect_dist = np.arange(x0_2,x1_2,step2)
        distance_in_kpc = np.concatenate((first_part_vect_dist,second_part_vect_dist))

    return(distance_in_kpc)


# In[41]:


def creating_distance_vector_right_sampling(vector_dists):
    """
    the dist vector for the calculating the surface brightness does not start at 0. I should change the final distance vector for being representative of the step/2
    /*example*/ input: 1,2,3,4,5; output: 0.5,1.5,2.5,3.5,4.5
    """
    corrected_vector_dists = np.full_like(vector_dists, float("NaN") )
    
    for ii, element in enumerate(vector_dists):
        if ii == 0: previous_element = 0
        corrected_vector_dists[ii] = previous_element + ((element - previous_element)/2.)
        previous_element = element
    return(corrected_vector_dists)


# In[42]:


step_1 = 0.5 #kpc
end = 28. #kpc
transition_1_2 = 3. #kpc
step_2 = 0.5 #kpc
threshold = 5.
path_output = "./sb_profiles/"
sb_threshold = 27.5 #level apart from which we should be dominated by the sky
#


# In[43]:


bands=np.array(["u","g","r","i","z"])
colors=np.array(["purple","blue","green","red","black"])
lambdas=np.array(["358.0","475.4","620.4","769.8","966.5"])


# In[44]:


zz=0.01690
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
distance_in_kpc = creating_galactocentric_distance_vector_in_kpc(step_1,end,transition_1_2,step_2)
resolution = distance_in_kpc.size
distance_in_kpc_right_sampling = creating_distance_vector_right_sampling(distance_in_kpc)
kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(zz)
aa = distance_in_kpc/kpc_per_arcsec.value                                   #calculating the semi-major axis
distance_in_arcsec = np.array(aa)                                                  #copying and making them two independent arrays
distance_in_arcsec_right_sampling = distance_in_kpc_right_sampling / kpc_per_arcsec.value


# In[45]:


dist_arcsec={}
dist_kpc={}
sb={}
sb_err={}
for band in bands:
    tt = Table.read("data_"+band+".cat",
                    names=["vector_dist_arcsec","vector_dist_kpc","sb","sb_err"], 
                    format="ascii.commented_header")
    
    
    dist_arcsec[band]=np.array(tt["vector_dist_arcsec"])
    dist_kpc[band]=np.array(tt["vector_dist_kpc"])
    sb[band]=np.array(tt["sb"])
    sb_err[band]=np.array(tt["sb_err"])


# In[46]:


dist_ref_for_sky={}
sky_dist1={}
sky_dist2={}
flux_subtracted={}
num_pixs_sky={}
sigmasb={}


for band in bands:
    tt = Table.read("fluxes_subtracted_"+band+".cat",
                    names=["galaxy","dist_ref_for_sky","sky_dist1","sky_dist2","flux_subtracted","num_pixs_sky","3sigmasb"], 
                    format="ascii.commented_header")
    
    
    dist_ref_for_sky[band]=np.array(tt["dist_ref_for_sky"])
    sky_dist1[band]=np.array(tt["sky_dist1"])
    sky_dist2[band]=np.array(tt["sky_dist2"])
    flux_subtracted[band]=np.array(tt["flux_subtracted"])
    num_pixs_sky[band]=np.array(tt["num_pixs_sky"])
    sigmasb[band]=np.array(tt["3sigmasb"])


# In[47]:


fig = plt.figure(figsize=(10,10))
plot_options()
ax = fig.add_subplot(1,1,1)

for i,band in enumerate(bands):
    ax.scatter(dist_arcsec[band],sb[band],s=10,color=colors[i],label="Banda "+band+": "+lambdas[i]+" nm")
    ax.errorbar(dist_arcsec[band], sb[band], yerr = sb_err[band] , color=colors[i])
for i,band in enumerate(bands):   
    plt.axhline(sigmasb[band],color=colors[i],linestyle="--" )
    
    
    


ax.legend()    
#limites
min_dist = 0.   #X valor minimo
max_dist = 45.  #X valor maximo
mag_y0 = 31.
mag_y1 = 15.
#Límites
ax.set_xlim(min_dist,max_dist)
ax.set_ylim(mag_y0,mag_y1)
plt.grid()
#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc["u"],sb ["u"], alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()
#Titulos
ax.set_ylabel("Brillo [mag/arcesc²]")
ax.set_xlabel("Distancia [arcsec]")
#ax2.set_xlabel("dist_kpc")
plt.title("Perfiles Galaxia NGC1277 (sustraida)")

#plt.show()

plt.savefig("Perfil_ngc1277.pdf")


# In[48]:


for band in bands:
    print(sigmasb[band])
    print(band,"\n")


# In[49]:


for i,band in enumerate(bands):   
    print(band,sigmasb[band])


# In[50]:


for i,band in enumerate(bands):
    print("Bandas :", band)


# In[51]:


perfcolor={}
#Color u
perfcolor["u_g"]=sb["u"]-sb["g"]
perfcolor["u_r"]=sb["u"]-sb["r"]
perfcolor["u_i"]=sb["u"]-sb["i"]
perfcolor["u_z"]=sb["u"]-sb["z"]



#Color g
perfcolor["g_u"]=sb["g"]-sb["u"]
perfcolor["g_r"]=sb["g"]-sb["r"]
perfcolor["g_i"]=sb["g"]-sb["i"]
perfcolor["g_z"]=sb["g"]-sb["z"]




#Color r
perfcolor["r_u"]=sb["r"]-sb["u"]
perfcolor["r_g"]=sb["r"]-sb["g"]
perfcolor["r_i"]=sb["r"]-sb["i"]
perfcolor["r_z"]=sb["r"]-sb["z"]



#Color i
perfcolor["i_u"]=sb["i"]-sb["u"]
perfcolor["i_g"]=sb["i"]-sb["g"]
perfcolor["i_r"]=sb["i"]-sb["r"]
perfcolor["i_z"]=sb["i"]-sb["z"]



#Color z
perfcolor["z_u"]=sb["z"]-sb["u"]
perfcolor["z_g"]=sb["z"]-sb["g"]
perfcolor["z_r"]=sb["z"]-sb["r"]
perfcolor["z_i"]=sb["z"]-sb["i"]


# In[68]:


fig = plt.figure(figsize=(10,10))
plot_options()
ax = fig.add_subplot(1,1,1)

ax.scatter(dist_arcsec["g"],perfcolor["g_r"],color=colors[0],label="gr")
ax.errorbar(dist_arcsec["g"],perfcolor["g_r"],color=colors[0])


ax.scatter(dist_arcsec["g"],perfcolor["g_i"],color=colors[1],label="gi")
ax.errorbar(dist_arcsec["g"],perfcolor["g_i"],color=colors[1])


ax.scatter(dist_arcsec["g"],perfcolor["g_z"],color=colors[2],label="gz")
ax.errorbar(dist_arcsec["g"],perfcolor["g_z"],color=colors[2])


ax.scatter(dist_arcsec["g"],perfcolor["r_i"],color=colors[3],label="ri")
ax.errorbar(dist_arcsec["g"],perfcolor["r_i"],color=colors[3])


ax.scatter(dist_arcsec["g"],perfcolor["r_z"],color=colors[4],label="rz")
ax.errorbar(dist_arcsec["g"],perfcolor["r_z"],color=colors[4])


ax.scatter(dist_arcsec["g"],perfcolor["i_z"],color="yellow",label="iz")
ax.errorbar(dist_arcsec["g"],perfcolor["i_z"],color="yellow")


ax.legend()


#limites
min_dist = 0.   #X valor minimo
max_dist = 45.  #X valor maximo
mag_y0 = 0.
mag_y1 = 2.5

#Límites
ax.set_xlim(min_dist,max_dist)
ax.set_ylim(mag_y0,mag_y1)
plt.grid()




#Tercer eje
ax2 = ax.twiny()
ax2.plot(dist_kpc["u"],sb ["u"], alpha=0)
ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec.value)
ax2.set_xlabel("Distancia [kpc]")
plt.grid()

#Titulos
ax.set_ylabel("Magnitud del Color [mag/arcesc²]")
ax.set_xlabel("Distancia [arcsec]")
#ax2.set_xlabel("dist_kpc")
plt.title("Perfiles de Color Galaxia NGC1277")
plt.savefig("Perfil_Color_ngc1277.pdf")


# In[ ]:




