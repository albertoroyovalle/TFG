"""
Optical masks based on the I-band image, actually I only create the I-band mask
Near-infrared masks based on the H-band image, actually I only create the H-band mask
Optical and near-infrared mask have different parameters (see get_masking_options function)
I unmask objects at distances < too_close [kpc]
I mask objects at distances > too_far_away [kpc]
It is very important to set the variable hudf_flag to the right value. Depending its value, I create the masks (or not) for HUDF objects
OUTPUTS:
- object_status.txt -> Galaxies for which the program does not find the image gal+"_"+field[cont_gal]+"_"+band
- if hudf_flag == False: gal+"_"+field[cont_gal]+"_"+band+     "_mask_all_neigh.fits
- if hudf_flag == True: gal+"_"+field[cont_gal]+"_"+band+"_hudf_mask_all_neigh.fits"
"""

import numpy as np
import matplotlib.pyplot as plt
import pdb
import os 
from astropy.io import fits
from astropy.table import Table                                                       #Table.read for Creating diccionaries from the information of sextract
from astropy.cosmology import FlatLambdaCDM                           #?

#####Import functions from other codes####

from my_python_library_v2 import dist_ellipse 
#from functs_for_truncations import find_right_image, obtaining_best_redshifts    #Ask for this code: We dont need this function.

#from create_mask_tuned_nir import create_mask_tuned_nir         #Create mask
#from create_mask_tuned_opt import create_mask_tuned_opt       #Create mask  Diferences?
from detecting_central_obj import detecting_central_obj               #It detects the central object of the image(AKA the galaxy we are going to study)

####################################
####################################
####################################
def get_masking_options(band):      
    if band == "I":
        pix_scale = 0.03
        enlargemask = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size = 40
        enlargemask_far = 7 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 60
        
    elif band == "H":
        pix_scale = 0.06
        enlargemask = 3 #enlarge ellipses factor for the masking with ellipses function
        max_size = 20
        enlargemask_far = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 40
        
    elif band =="Z":
        pix_scale = 0.339
        enlargemask = 10 #enlarge ellipses factor for the masking with ellipses function
        max_size = 20
        enlargemask_far = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 40
    if band == "I":
        pix_scale = 0.396
        enlargemask = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size = 40
        enlargemask_far = 7 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 60
        
    elif band == "z":
        pix_scale = 0.396
        enlargemask = 3 #enlarge ellipses factor for the masking with ellipses function
        max_size = 20
        enlargemask_far = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 40
        
    elif band =="u":
        pix_scale = 0.396
        enlargemask = 10 #enlarge ellipses factor for the masking with ellipses function
        max_size = 20
        enlargemask_far = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 40
    elif band=="i":
        pix_scale = 0.396
        enlargemask = 3 #enlarge ellipses factor for the masking with ellipses function
        max_size = 20
        enlargemask_far = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 40
    
    elif band=="r":
        pix_scale = 0.396
        enlargemask = 3 #enlarge ellipses factor for the masking with ellipses function
        max_size = 20
        enlargemask_far = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 40
    
    elif band=="g":
        pix_scale = 0.396
        enlargemask = 3 #enlarge ellipses factor for the masking with ellipses function
        max_size = 20
        enlargemask_far = 5 #enlarge ellipses factor for the masking with ellipses function
        max_size_far = 40
    return(pix_scale, enlargemask, enlargemask_far, max_size, max_size_far)
        
def mask_obj(primary_mask,x_size,y_size,xc,yc,axis_ratio,position_angle,a_image,kron_radius,enlargemask,max_size):
    new_img = dist_ellipse([x_size,y_size], xc, yc, axis_ratio, position_angle)

    first_size = a_image*kron_radius
    if first_size > 5.:
        final_size = first_size*enlargemask 
    else:
        final_size = 5.*enlargemask
        
    if final_size > max_size: final_size = max_size

    np.putmask(new_img,new_img < final_size,1) #pixels I will mask
    np.putmask(new_img,new_img > final_size,0) #pixels I will not mask

    np.logical_or(primary_mask,new_img,primary_mask) #v1,v2,result

    return(primary_mask)

def unmask_obj(primary_mask,x_size,y_size,xc,yc,axis_ratio,position_angle,a_image,kron_radius,enlargemask,max_size):
    new_img = dist_ellipse([x_size,y_size], xc, yc, axis_ratio, position_angle)

    first_size = a_image*kron_radius
    if first_size > 5.:
        final_size = first_size*enlargemask 
    else:
        final_size = 5.*enlargemask

    if final_size > max_size: final_size = max_size

    np.putmask(new_img,new_img < final_size,0) #pixels I will mask
    np.putmask(new_img,new_img > final_size,1) #pixels I will not mask

    np.logical_and(primary_mask,new_img,primary_mask) #v1,v2,result

    return(primary_mask)
    
def create_mask_tuned(img_filename,sex_output,target,enlargemask,filter_objs,max_size):                       

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

    img = fits.open(img_filename)
    header = img[0].header
    y_size = header['NAXIS1']
    x_size = header['NAXIS2']

    primary_mask = np.zeros( (x_size,y_size), dtype=bool)

    #masking all objects
    for ii in range(len(objs)):
        primary_mask = mask_obj(primary_mask,x_size, y_size, xc[ii], yc[ii], 1./(1.-ellip[ii]), theta[ii],a_image[ii],kron[ii],enlargemask,max_size)

    #unmasking the central region
    primary_mask = unmask_obj(primary_mask,x_size,y_size,xc[target],yc[target],1./(1.-ellip[target]), theta[target],a_image[target],kron[target],enlargemask,max_size)
   
    #unmasking selected neighbours. It only works when filter_objs != from None
    if filter_objs!=None:
        for ii in range(len(filter_objs)):
            if (ii != target) & (filter_objs[ii]==True):
                primary_mask = unmask_obj(primary_mask,x_size,y_size,xc[ii],yc[ii],1./(1.-ellip[ii]), theta[ii],a_image[ii],kron[ii],enlargemask,max_size)
    #in order to be accepted by Cappellari's algorithms
    ones  = primary_mask == 1
    zeros = primary_mask == 0
    primary_mask[ones] = 0
    primary_mask[zeros] = 1
    
    return(primary_mask.astype(int))


#########################################################
#########################################################
#########################################################
def main():
    threshold_rms = None                                 #threshold for removing pixels in the rms image
    #bands = np.array(["I","H"])                           #all optical bands come from the I band, all NIR ones from the H band
    #bands=np.array(["I","H","Z"])
    
    #bands=np.array(["Z"])
    #bands=np.array(["i"])
    #bands=np.array(["r"])
    #bands=np.array(["g"])
    
    type_of_analysis = "single"
    #catalog_path = "./"
    #img_path = "./galaxy_images/"
    #rms_path = "./sigma_images/"
    path_output = "./masking_all_but_target/"
    #catalog = "galaxies_main_props_with_struct_params_and_no_bad_objects.cat"
    #number_cpus = 3
    too_close = 10. #kpc
    too_far_away = 50. #kpc
    segmentation_cat_filename="first_z.cat"
    img_filename=input("Image filename: ")
    bands=np.array([input("Band: ")])
    #img_filename="79071-viking-z.fits"



    #########################
    
    file_log = open("./object_status.txt", "w")
    
    ###Creation of a directory for our output data
    if not os.path.exists(path_output):
        os.mkdir(path_output)

            ####################################
            #################################### FOR in bands
    for cont_band,band in enumerate(bands):
        #reading the SExtractor output catalog=================================================================================
        
        try:
            objs,ra,dec,mag,magerr,x_image,y_image,flux_radius,ellipticity,theta,kron,a_image =np.loadtxt(segmentation_cat_filename, dtype = float, unpack = True)
        except:
            file_log.write(' segmentation_empty\n')
            continue
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
        #pdb.set_trace() ---------------------------> Hasta aqui todo correcto


#====================================================================================================
                
        #detecting the galaxy to be analyzed============================    Necesitare incluir tambien el .fits de la galaxia  Para leer los Naxis1 y Naxis2
        hdu = fits.open(img_filename)
        hdr = hdu[0].header
        x_size = hdr["NAXIS1"]
        y_size = hdr["NAXIS2"]
    
        target = detecting_central_obj(sex_output,x_size/2.,y_size/2.)                                                                        ###We detect the central object
        pix_scale,enlargemask,enlargemask_far,max_size,max_size_far = get_masking_options(band)
        #pdb.set_trace()
        #!!!ERROR!!!: TypeError: cannot unpack non-iterable NoneType object
        #ERROR solucionado

        #print("Error unmasking the neighbours")
        filter_objs=None
        mask = create_mask_tuned(img_filename,sex_output,target,enlargemask,filter_objs,max_size) 
        fits.writeto  (path_output+"_"+"_mask_all_neigh.fits",mask,hdr,overwrite=True)

if __name__ == "__main__":
    main()
