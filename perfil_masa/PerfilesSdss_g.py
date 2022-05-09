#!/usr/bin/env python
# coding: utf-8

#   Hsex--->MAG_ZEROPOINT	22.5
# 	PIXEL_SCALE	0.396
#  
#  La banda mas roja es la de z
#  X:251.4020    
#  Y:252.1227 
#  ellip:0.475  (interesa pasar a ar=1./(1.-ellip[target]) )
#  pa:57.73
#  
#  
# 

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pdb #for debugging purposes only
from astropy.cosmology import FlatLambdaCDM
import matplotlib as mpl
from astropy.stats import sigma_clipped_stats
import scipy.stats as stt
from astropy.table import Table
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib import rc,rcParams
import sys
from my_python_library_v2 import finding_coo_brightest_pix, resistant_mean, dist_ellipse, calculate_flux_profile, sb_from_flux
from functs_for_truncations import obtaining_best_ar_and_pa, find_right_image, obtaining_best_redshifts, find_right_mask
from detecting_central_obj import detecting_central_obj 
from my_python_library_v2 import calculate_limiting_mag
from my_python_library_v2 import calculate_mass_profile


# In[3]:


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


# In[4]:


def phot_options(field,band):
    if   band == 'u':
        pixel_scale = 0.396 #[arcsec/pix]
        #if   (field == "goodss") == True: zeropoint = 25.68
        #elif (field == "goodsn") == True: zeropoint = float("NaN")
        #elif (field == "cosmos") == True: zeropoint = float("NaN")
        #elif (field == "uds") == True:    zeropoint = float("NaN")
        #elif (field == "egs") == True:    zeropoint = float("NaN")
        #elif (field == "hudf") == True:   zeropoint = 25.673
        if (field == "sdss") == True:   zeropoint = 22.5
    elif band == 'g':
        pixel_scale = 0.396 #[arcsec/pix]
        #if   (field == "goodss") == True: zeropoint = 26.51 #26.493
        #elif (field == "goodsn") == True: zeropoint = 26.493
        #elif (field == "cosmos") == True: zeropoint = 26.493
        #elif (field == "uds") == True:    zeropoint = 26.49113
        #elif (field == "egs") == True:    zeropoint = 26.49113
        #elif (field == "hudf") == True:   zeropoint = 26.486
        if (field == "sdss") == True:   zeropoint = 22.5
    elif band == 'r':
        pixel_scale = 0.396 #[arcsec/pix]
        #if   (field == "goodss") == True: zeropoint = 25.94 #25.947
        #elif (field == "goodsn") == True: zeropoint = 25.947
        #elif (field == "cosmos") == True: zeropoint = 25.947
        #elif (field == "uds") == True:    zeropoint = 25.94333
        #elif (field == "egs") == True:    zeropoint = 25.94333
        #elif (field == "hudf") == True:   zeropoint = 25.654
        if (field == "sdss") == True:   zeropoint = 22.5
    elif band == 'i':
        pixel_scale = 0.396 #[arcsec/pix]
        #if   (field == "goodss") == True: zeropoint = 26.27
        #elif (field == "goodsn") == True: zeropoint = float("NaN")
        #elif (field == "cosmos") == True: zeropoint = float("NaN")
        #elif (field == "uds") == True:    zeropoint = float("NaN")
        #elif (field == "egs") == True:    zeropoint = float("NaN")
        #elif (field == "hudf") == True:   zeropoint = 26.269
        if (field == "sdss") == True:   zeropoint = 22.5
    elif band == 'z':
        pixel_scale = 0.396 #[arcsec/pix]
        #if   (field == "goodss") == True: zeropoint = 26.23 #26.2303
        #elif (field == "goodsn") == True: zeropoint = 26.2303
        #elif (field == "cosmos") == True: zeropoint = 26.2303
        #elif (field == "uds") == True:    zeropoint = 26.25
        #elif (field == "egs") == True:    zeropoint = 26.25
        #elif (field == "hudf") == True:   zeropoint = 26.230
        if (field == "sdss") == True:   zeropoint = 22.5#
        
    return(zeropoint, pixel_scale)


# In[5]:


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


# In[6]:


#opt_bands = ("V","B","I") #OJO: first V because I use it to calculate where to subtract the sky. it's a tuple because I want this to be an inmutable object
#nir_bands = ("Y","J","H")

#Sdss_bands = ("u","g","r","i","z")

bands = np.array(["g"])  ###CAMBIAR
b="g"                    ###CAMBIAR
field="sdss"             ###CAMBIAR

step_1 = 0.5 #kpc
end = 30. #kpc
transition_1_2 = 3. #kpc
step_2 = 0.5 #kpc
threshold = 5.
path_output = "./sb_profiles/"
#sb_threshold = 27.5 #level apart from which we should be dominated by the sky
sb_threshold = 24

gals=np.array(["79071"])
gal=np.array(["79071"]) #Seguir mirando esto...(Mirar si puedo dejar el bucle ||||| for ii,gal in enumerate(gals) )
 

   
zz=0.01690
img_filename="ngc1277_sdss_g.fits"     ####CAMBIAR
#rms_filename=
mask_filename="ngc1277_final_sdss_mask.fits" #### CAMBIAR
segmentation_cat_filename="first_z.cat"  ##### CAMBIAR
###Leo los datos de sextractor first.cat-----------


# In[7]:


#Leo los datos de sextractor first.cat--------------
#----------> El objetivo de esto es encontrar el objeto central(target) para despues determinar centroid=[(x_image,y_image)]        theta        ar(axis ratio=A/B)

    #Voy a necesitar:
#1 NUMBER                 Running object number                                     
#2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
#3 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
#4 MAG_APER               Fixed aperture magnitude vector                            [mag]
#5 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
#6 X_IMAGE                Object position along x                                    [pixel]
#7 Y_IMAGE                Object position along y                                    [pixel]
#8 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
#9 ELLIPTICITY            1 - B_IMAGE/A_IMAGE                                       
#10 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
#11 KRON_RADIUS            Kron apertures in units of A or B                         
#12 A_IMAGE                Profile RMS along major axis                               [pixel]


#OJO que en el caso de las bandas sdss Vamos a utilizar los centroides y datos de la banda mas roja: z


# In[8]:


try:
    objs,ra,dec,mag,magerr,x_image,y_image,flux_radius,ellipticity,theta,kron,a_image =np.loadtxt(segmentation_cat_filename, dtype = float, unpack = True)
except:
    file_log.write(' segmentation_empty\n')

sex_output = np.array([objs,ra,dec,mag,magerr,x_image,y_image,flux_radius,ellipticity,theta,kron,a_image],float)     ##Lo que nos da sextractor
           
if sex_output.ndim > 1:
    objs  = sex_output[0,:]
    xc    = sex_output[5,:]
    yc    = sex_output[6,:]
    ellip = sex_output[8,:]
    theta = sex_output[9,:]
    kron  = sex_output[10,:]
    a_image = sex_output[11,:]
else:
    objs  = np.array([sex_output[0]])
    xc    = np.array([sex_output[5]])
    yc    = np.array([sex_output[6]])
    ellip = np.array([sex_output[8]])
    theta = np.array([sex_output[9]])
    kron  = np.array([sex_output[10]])
    a_image = np.array([sex_output[11]])


# In[9]:


##trozo de codigo para:
##Detectar el centroide, pa, etc. 
#Se necesita el "first_"+banda+".cat" correspondiente

hdu = fits.open(img_filename)
hdr = hdu[0].header
x_size = hdr["NAXIS1"]
y_size = hdr["NAXIS2"]
target = detecting_central_obj(sex_output,x_size/2.,y_size/2.)   #Busco el objeto central(Target), función que se basa en buscar el elementro central de la imagen
xcentr=xc[target]
#xcentr=368#Coordenada x de nuestro target
ycentr=yc[target]   
#ycentr=366#Coordenada y de nuestro target
centroid=np.array([xcentr,ycentr])         ##Encuentro el centroide del target
pa=theta[target] 
#pa=155#Angulo de posicion del target
ar=1./(1.-ellip[target]) 
#ar=0.8#axis ratio-Razón de semiejes. Pag 33 Sextractor User's manual
#pdb.set_trace() Hasta aqui todo correcto
zeropoint,pixel_scale = phot_options(field,"z")      ##b=banda que utilizo
print("Zp, ",zeropoint,"Pixel scale: ",pixel_scale)
print("Centroide: ",xcentr,ycentr)
print("Position angle (pa): ",pa)
print("Axis ratio (ar) : ", ar)
print("ellip ", ellip[target])


# In[10]:


cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
distance_in_kpc = creating_galactocentric_distance_vector_in_kpc(step_1,end,transition_1_2,step_2)
resolution = distance_in_kpc.size
distance_in_kpc_right_sampling = creating_distance_vector_right_sampling(distance_in_kpc)


# In[11]:


if not os.path.exists(path_output):  #Creamos directorio donde volcaremos los resultados
    os.mkdir(path_output)


# In[12]:


for band in bands:
    flags                   = { band : np.zeros_like(gals, dtype=int)   for band in bands }                                                             
    vector_dist_ref_for_sky = { band : np.full_like (gals,float("NaN")) for band in bands }
    vector_sky_dist1        = { band : np.full_like (gals,float("NaN")) for band in bands }
    vector_sky_dist2        = { band : np.full_like (gals,float("NaN")) for band in bands }
    fluxes_subtracted       = { band : np.full_like (gals,float("NaN")) for band in bands }
    sky_errors_vector       = { band : np.full_like (gals,float("NaN"),dtype=float) for band in bands }
    num_pixs_aper_vector    = { band : np.full_like (gals,float("NaN"),dtype=float) for band in bands }
    #pdb.set_trace() Hasta aquí todo bien -----------> Diccionarios  ={"Banda": Nan}     np.full_like (gals,float("NaN")) crea un array igual que gals con solo elementos "NaN"


# In[13]:



kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(zz)
aa = distance_in_kpc/kpc_per_arcsec.value                                   #calculating the semi-major axis
distance_in_arcsec = np.array(aa)                                                  #copying and making them two independent arrays
distance_in_arcsec_right_sampling = distance_in_kpc_right_sampling / kpc_per_arcsec.value

#size=3.6*kpc_per_arcsec #10 arcsec * (kpc/arcsec)-----> distancia donde creemos que termina la galaxia
#print(size)
#type(size)
#No recuerdo para que hice esto de size, supongo que fuera para checkear el limite de la galaxia en kpc
#print(distance_in_kpc_right_sampling,distance_in_arcsec_right_sampling)
print(30*kpc_per_arcsec.value)

print("Kpc/arcsec: ", kpc_per_arcsec)
print("Arcsec/kpc: ", 1/kpc_per_arcsec)
print("\n")
print("Distance in arcsec: ",distance_in_arcsec, "\n")

print("Distance in kpc: ",distance_in_kpc)
print("\n")
print("Distance in kpc right sampling: ",distance_in_kpc_right_sampling)

print("\n","He utilizado la regla en sextractor el el radio a ojo es de unos 23 pixeles---|> 23 pix*0.396arcsec/pix=8.712 arcsec")


# OJO:
#     Photo options cambiar a banda utilizada siempre
#     if "band" == cambiar a la banda utilizada
#     
#     Despues de ver el primer perfil sin el flujo sustraido sleccionar vector ref for sky

# In[14]:


#### for jj,band in enumerate(bands):                                         ###Puedo dejar este bucle si en 143 bands="Z"  Mirar si las funciones incluyen los parametros de esta banda
#####Determinamos el zeropoint y el pixelscale (Se podria haber introducido directamente, yo he modificado la funcion para cuando field=None y band=("Z") nos de los valores zeropoint=30 pixel_scale=0.339)
    zeropoint, pixel_scale = phot_options(field,band)
    
##################Abrimos la imagen y leemos los datos, header...

    hdu = fits.open(img_filename)
    img = hdu[0].data
    hdr = hdu[0].header
    zero_eles = np.count_nonzero(np.logical_or(np.isnan(img),img==0)) #it actually counts the number of elements with zero values
    if zero_eles > img.size/2.: #if more than half the image is out of the border, I don't take the image
        continue
    """
    ###No disponemos de rms
    hdu = fits.open(rms_filename)                                    ###rms_filename=.....?????????
    rms = hdu[0].data
    """
    
###################Abrimos la mascara y leemos datos

    hdu = fits.open(mask_filename)
    mask = hdu[0].data
    """
    Primer perfil de luminosidas(y flujo)!!!!!!!!!!!!!!!
    Ya podemos determinar el flujo y demas cosas ¡necesarias para determinar el brillo superficial sb-----------------> Perfil de luminosidad en dos pasos:
    Determinar perfil de flujo con(calculate_flux_profile(distance_in_arcsec,centroid,ar,pa,img,mask,pixel_scale,zeropoint))
    Determinar el perfil de luminosidad a partir del de flujo:(sb_from_flux(flux,zeropoint,pixel_scale,sigma_flux,num_pixs_aper))
         
    Necesitaremos:
    distance_in_arcsec= distance_in_kpc/kpc_per_arcsec.value (ojo que array)
    centroid             = Centroide del target
    ar                       = Axis ratio--- Razon de semiejes
    pa                      = Position angle-- Angulo de posicion (Respecto x (Revisar))
    img                    = hdu[0].data (Datos imagen)
    mask                 = hdu[0].data (Datos mask)
    pixel_scale        =0.339
    zeropoint           =30
    """
    print("centroid: ",centroid,"\n", "ar: ",ar,"\n", "pa: ",pa,"\n", "img: ",img_filename,"\n", "mask ",mask_filename,"\n", "pixel_scale ",pixel_scale ,"\n","zeropoint: ",zeropoint)
    flux1,sigma_flux1,num_pixs_aper1  = calculate_flux_profile(distance_in_arcsec,centroid,ar,pa,img,mask,pixel_scale,zeropoint)   ##      ------> Recordar modificar la funcion en my_python_library_v2: def calculate_flux_profile(distance_vector,centroid,axis_ratio,position_angle,img,mask,pix_scale,zeropoint,gal=None,check=False):
    sb1, err_sb1 = sb_from_flux(flux1,zeropoint,pixel_scale,sigma_flux1,num_pixs_aper1)
    #pdb.set_trace() Hasta aqui todo bien
####################################
####################################
    """
    Dudas con esto, repasarlo bien
    """
#determining where to subtract the pedestal in the sky
    if band == "g":                                              #Cambiar band V por Z?
        for kk,sb_level in enumerate(sb1):
            if sb_level > sb_threshold:                    #sb_threshold definido al principio del programa
                break
        dist_ref_for_sky1 = distance_in_arcsec[kk]*2.
        dist_ref_for_sky2 = distance_in_arcsec[kk]*3.
        zp_none, pix_scale = phot_options(field,band)
        print("zpnone pixscale ",zp_none, pix_scale)
        max_dist_ref_for_sky = (hdr["NAXIS1"]*pix_scale)/2.
        if dist_ref_for_sky2 > max_dist_ref_for_sky: 
            dist_ref_for_sky2 = max_dist_ref_for_sky
            dist_ref_for_sky1 = max_dist_ref_for_sky-2.
        #pdb.set_trace() 
    

    
    #vector_ref_for_sky = np.array([dist_ref_for_sky1,dist_ref_for_sky2]) #Distances between where I measure the sky
    vector_ref_for_sky=np.array([40,50])
    print("vector ref",vector_ref_for_sky)
    #print(type(vector_ref_for_sky),vector_ref_for_sky)

    
    if ar < 0.3: #in order to measure a big chunk of the sky, because otherwise the large ellipticity prevents us to measure enough pixels in the annulus
        tmp_ar = 0.3 
    else:
        tmp_ar = ar
            
    """
    2º Perfil de flujo
    Preguntar de nuevo y apuntar en condiciones que es esto...
    He entendido que es algo asi como intentar "quitar el ruido de fondo" ----------------Once I subtract the sky level
    ruido: Ruido base()+pequeñas variaciones de ruido
    """
    flux2,sigma_flux2,num_pixs_aper2 = calculate_flux_profile(vector_ref_for_sky,centroid,tmp_ar,pa,img,mask,pixel_scale,zeropoint) #######PERFIL
    #pdb.set_trace()
    pixels_with_no_info = img == 0 #to preserve their zero value
    img = img - flux2[1] #as in Pohlen & Trujillo (2006) and Trujillo+19; [1] because it is the flux calculated at the right distance
    img[pixels_with_no_info] = 0
    vector_dist_ref_for_sky[band]= max_dist_ref_for_sky                                                                         ######AQUI LIO
    vector_sky_dist1       [band] = vector_ref_for_sky[0]
    vector_sky_dist2       [band] = vector_ref_for_sky[1]
    fluxes_subtracted      [band]= flux2[1] #[1] because it is the flux calculated at the right distance
    sky_errors_vector      [band]= sigma_flux2[1] #[1] because it is the flux calculated at the right distance
    num_pixs_aper_vector   [band] = num_pixs_aper2[1] #[1] because it is the flux calculated at the right distance
    
    print("fluxes subtracted",flux2[1])    
#=====================================================
    """
    Nuevo perfil de flujo con la imagen=img-flux[1]      Nueva imagen con el fondo sustraido
    """
    flux3,sigma_flux3,num_pixs_aper3  = calculate_flux_profile(distance_in_arcsec,centroid,ar,pa,img,mask,pixel_scale,zeropoint,gal)
    #none1,sigma_rms,num_pixs_rms  = calculate_flux_profile(distance_in_arcsec,centroid,ar,pa,rms,mask,pixel_scale,zeropoint,gal) #
    sb3, err_sb3 = sb_from_flux(flux3,zeropoint,pixel_scale,sigma_flux3,num_pixs_aper3,rms_bin=sky_errors_vector[band]/np.sqrt(num_pixs_aper_vector[band]))
    data = Table( [distance_in_arcsec_right_sampling, distance_in_kpc_right_sampling, sb3, err_sb3], names = ["vector_dist_arcsec","vector_dist_kpc","sb","sb_error"] )
    data.write(path_output+"data_"+band+".cat", format = "ascii.commented_header",overwrite=True)
    #pdb.set_trace() 
    #if the file with the subtracted sky exists, I update the values for this very galaxy    
    
    


# In[ ]:


print("Dist arcsec",dist_arcsec)
print("Dist Kpc", dist_kpc )
print("Flujo 1" , flux1)
print("sb1: " , sb1)


# In[ ]:


dist_arcsec=np.array([distance_in_arcsec_right_sampling]) 
dist_kpc=np.array([distance_in_kpc_right_sampling]) 
sb_Z=np.array([sb1]) 
sb_err_Z=np.array([err_sb1])
#kpc_per_arcsec


min_dist = 0.
max_dist = 60.
mag_y0 = 25.
mag_y1 = 17.

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.scatter(distance_in_arcsec_right_sampling, sb1)
ax.errorbar(distance_in_arcsec_right_sampling, sb1, yerr = err_sb1)
ax.set_xlim(min_dist,max_dist)
ax.set_ylim(mag_y0,mag_y1)
ax.yaxis.set_ticks_position('both')

ax2 = ax.twiny()
print(ax.get_xlim()[0]*kpc_per_arcsec)
print(ax.get_xlim()[1]*kpc_per_arcsec)
ax2.plot(distance_in_kpc_right_sampling, sb1, alpha=0)
#ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec,ax.get_xlim()[1]*kpc_per_arcsec)
ax.set_ylabel("sb")
ax.set_xlabel("dist_arsec")
ax2.set_xlabel("dist_kpc")
plt.title("Perfil Sin Flujo sustraido")
#plt.savefig("./perfimg/"+"Perfil Sin Flujo sustraido.jpg")
plt.show()


# In[ ]:


dist_arcsec=np.array([distance_in_arcsec_right_sampling]) 
dist_kpc=np.array([distance_in_kpc_right_sampling]) 
sb_Z=np.array([sb3]) 
sb_err_Z=np.array([err_sb3])
#kpc_per_arcsec


min_dist = 0.   #X valor minimo
max_dist = 60.  #X valor maximo
mag_y0 = 31.
mag_y1 = 16.
#Límites

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.scatter(distance_in_arcsec_right_sampling, sb3)
ax.errorbar(distance_in_arcsec_right_sampling, sb3, yerr = err_sb3)
ax.set_xlim(min_dist,max_dist)
ax.set_ylim(mag_y0,mag_y1)
ax.yaxis.set_ticks_position('both')

#plt.axhline(y=sb_limit-2.5*np.log10(3))

ax2 = ax.twiny()
print(ax.get_xlim()[0]*kpc_per_arcsec)
print(ax.get_xlim()[1]*kpc_per_arcsec)
ax2.plot(distance_in_kpc_right_sampling, sb3, alpha=0)
#ax2.set_xlim(ax.get_xlim()[0]*kpc_per_arcsec.value,ax.get_xlim()[1]*kpc_per_arcsec)
ax.set_ylabel("sb")
ax.set_xlabel("dist_arsec")
ax2.set_xlabel("dist_kpc")
plt.title("Perfil Con Flujo sustraido(CriterioBucle)")
#plt.savefig("./perfimg/"+"Perfil Con Flujo sustraido(CriterioBucle).jpg")
plt.show()


# In[ ]:


#max_dist_pix=5*(1/pixel_scale) ##Cambiar esto eh
max_dist_pix=(22/kpc_per_arcsec.value)*(1/pixel_scale)
print(max_dist_pix)
sb_limit=calculate_limiting_mag(img_filename,mask_filename,"joining_individual_pixels_to_create_apers","box",10,pixel_scale,1e4,zeropoint,ar,pa,max_dist_pix,xcentr,ycentr)
print("sigmaSb",sb_limit)
print("3simgaSb",sb_limit-2.5*np.log10(3))
sb3sigma=sb_limit-2.5*np.log10(3)


# In[ ]:


data = Table( [distance_in_arcsec_right_sampling, distance_in_kpc_right_sampling, sb1, err_sb1], names = ["vector_dist_arcsec","vector_dist_kpc","sb","sb_error"] )
data.write(path_output+"pre_data_g"+".cat", format = "ascii.commented_header",overwrite=True)


# In[ ]:


##Escribimos Flujos sustraidos (Para Jesus)
for band in bands:
    if os.path.isfile(path_output+"fluxes_subtracted_"+band+".cat") == True:
        tt_fluxes = Table.read(path_output+"fluxes_subtracted_"+band+".cat",names = ["galaxy",
                                                            "dist_ref_for_sky",
                                                            "sky_dist1","sky_dist2",
                                                            "flux_subtracted",
                                                            "num_pixs_sky","3sigmasb"],format = "ascii.commented_header")

        tt_fluxes["dist_ref_for_sky"]= np.array([vector_dist_ref_for_sky[band]])
        tt_fluxes["sky_dist1"]    = np.array(vector_sky_dist1 [band])
        tt_fluxes["sky_dist2"]    = np.array(vector_sky_dist2  [band])
        tt_fluxes["flux_subtracted"]= np.array([fluxes_subtracted   [band]])
        tt_fluxes["3sigmasb"]=np.array([sb3sigma])
    
    if os.path.isfile(path_output+"fluxes_subtracted_"+band+".cat") == False:
        data = Table([ gals,np.array([vector_dist_ref_for_sky[band]]),
                    np.array([vector_sky_dist1[band]]),
                    np.array([vector_sky_dist2[band]]),
                    np.array([fluxes_subtracted[band]]),
                    np.array([num_pixs_aper_vector[band]]),
                    np.array([sb3sigma])  ],
                    names = ["galaxy","dist_ref_for_sky","sky_dist1","sky_dist2","flux_subtracted","num_pixs_sky","3sigmasb"])
        data.write(path_output+"fluxes_subtracted_"+band+".cat", format = "ascii.commented_header",overwrite=True)
    else:
        tt_fluxes.write(path_output+"fluxes_subtracted_"+band+".cat", format = "ascii.commented_header",overwrite=True)        


# In[ ]:




