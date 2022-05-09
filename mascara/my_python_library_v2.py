import numpy as np
import os
import pdb
from astropy.io import fits
from astropy import wcs
from astropy.io import ascii
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from matplotlib import rc,rcParams
import matplotlib as mpl
from astropy.stats import sigma_clipped_stats
import scipy.stats as stt

def pix2coo(filename,xx,yy):

    hdu = fits.open(filename)
    ww = wcs.WCS(hdu[0].header)
    
    centroid = []
    centroid.append([float("NaN"),float("NaN")]) #when you use a single element, you need to add an empty coordinate
    centroid.append([xx,yy])
    centroid = np.array(centroid, dtype=float)

    empty, coos = ww.wcs_pix2world(centroid,1)
    
    ra  = coos[0]
    dec = coos[1]
    return(ra,dec)  #[1] to remove the empty coordinate
    
    
def coo2pix(filename,ra,dec):

    hdu = fits.open(filename)
    ww = wcs.WCS(hdu[0].header)
    
    coo = []
    coo.append([float("NaN"),float("NaN")]) #when you use a single element, you need to add an empty coordinate
    coo.append([ra,dec])
    coo = np.array(coo, dtype = float)
    
    empty, pixs = ww.wcs_world2pix(coo,1)

    pix_x = pixs[0]
    pix_y = pixs[1]
    return(pix_x,pix_y) #[1] to remove the empty coordinate
    
    
def finding_coo_brightest_pix(img_filename, zoom_center=50):
    """
    zoom center is the parameter that controls how much portion of the inner part of the image we take in order to look for the brightest pixel (just in case there are nearby objects that are brighter than our target)
    OUTPUT: In reality, the first position in the output corresponds to the y-coordinate in the image, and the second to the x
    """
    if os.path.isfile(img_filename) == True:
        #if the brightest pixel was the central (in counts, not in mag arcsec^-2), this function will only be return np.unravel_index(np.nanargmax( image ),image.shape)
        hdu = fits.open(img_filename)
        img = hdu[0].data
        hdr = hdu[0].header
        x_size = hdr["NAXIS1"]
        y_size = hdr["NAXIS2"]
        
        matrix = np.array(img)
        #looking in the inner part of the stamp
        offset_x = x_size/zoom_center
        offset_y = y_size/zoom_center
        submatrix = matrix[round((x_size/2)-offset_x):round((x_size/2)+offset_x),round((y_size/2)-offset_y):round((y_size/2)+offset_y)]
        if np.isnan(submatrix).all() == True: return(float("NaN"),float("NaN")) #if nothing in the central submatrix go out returning NaNs
        indices_max_value = np.unravel_index(np.nanargmax(submatrix),submatrix.shape)
        filter_pos = matrix == submatrix[indices_max_value[0],indices_max_value[1]] #I assume that there is only a single pixel with the brightest value
        pos = np.where(filter_pos)
        return(int(pos[0][0]),int(pos[1][0]))
        #np.nanargmin->gets the minimum element of the image without taking into account the NaN values
        #np.unravel_index->gets the two dimensional position from a coordinate in a flatten vector
    else:
        return(float("NaN"),float("NaN"))
        
        
def finding_coo_brightest_pix_without_file(img, dim_x, dim_y, zoom_center=50):
    """
    zoom center is the parameter that controls how much portion of the inner part of the image we take in order to look for the brightest pixel (just in case there are nearby objects that are brighter than our target)
    OUTPUT: In reality, the first position in the output corresponds to the y-coordinate in the image, and the second to the x
    """
    if np.isnan(img).all() != True: #if not all the elements are NaN
        #if the brightest pixel was the central (in counts, not in mag arcsec^-2), this function will only be return np.unravel_index(np.nanargmax( image ),image.shape)
        x_size = dim_x
        y_size = dim_y
        
        matrix = np.array(img)
        #looking in the inner part of the stamp
        offset_x = x_size/zoom_center
        offset_y = y_size/zoom_center
        submatrix = matrix[round((x_size/2)-offset_x):round((x_size/2)+offset_x),round((y_size/2)-offset_y):round((y_size/2)+offset_y)]
        if np.isnan(submatrix).all() == True: return(float("NaN"),float("NaN")) #if nothing in the central submatrix go out returning NaNs
        indices_max_value = np.unravel_index(np.nanargmax(submatrix),submatrix.shape)
        filter_pos = matrix == submatrix[indices_max_value[0],indices_max_value[1]] #I assume that there is only a single pixel with the brightest value
        pos = np.where(filter_pos)
        return(int(pos[0][0]),int(pos[1][0]))
        #np.nanargmin->gets the minimum element of the image without taking into account the NaN values
        #np.unravel_index->gets the two dimensional position from a coordinate in a flatten vector
    else:
        return(float("NaN"),float("NaN"))
        
    
def finding_coo_brightest_SB(img_filename, zoom_center=50):
    """
    zoom center is the parameter that controls how much portion of the inner part of the image we take in order to look for the brightest pixel (just in case there are nearby objects that are brighter than our target)
    """
    if os.path.isfile(img_filename) == True:
        #the brightest surface brightness = the pixel with the lowest value
        #if the brightest pixel was the central, this function will only be return np.unravel_index(np.nanargmin( image ),image.shape)
        hdu = fits.open(img_filename)
        img = hdu[0].data
        hdr = hdu[0].header
        x_size = hdr["NAXIS1"]
        y_size = hdr["NAXIS2"]
        
        matrix = np.array(img)
        #looking in the inner part of the stamp
        offset_x = x_size/zoom_center
        offset_y = y_size/zoom_center
        submatrix = matrix[round((x_size/2)-offset_x):round((x_size/2)+offset_x),round((y_size/2)-offset_y):round((y_size/2)+offset_y)]
        """
        if x_size == y_size:
            submatrix = matrix[round((x_size/2)-offset_x):round((x_size/2)+offset_x),round((y_size/2)-offset_y):round((y_size/2)+offset_y)]
        else:
            submatrix = matrix #in the case the stamp is not a square, just look in the whole stamp
        """
        if np.isnan(submatrix).all() == True: return(float("NaN"),float("NaN")) #if nothing in the central submatrix go out returning NaNs
        indices_min_value = np.unravel_index(np.nanargmin(submatrix),submatrix.shape)
        filter_pos = matrix == submatrix[indices_min_value[0],indices_min_value[1]] #I assume that there is only a single pixel with the brightest value
        pos = np.where(filter_pos)
        return(int(pos[0][0]),int(pos[1][0]))
        #np.nanargmin->gets the minimum element of the image without taking into account the NaN values
        #np.unravel_index->gets the two dimensional position from a coordinate in a flatten vector
    else:
        return(float("NaN"),float("NaN"))
        
    
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
    
    
def resistant_mean(vector, threshold = 5.):
    """
    this function is (numerically) the same as the resistant_mean function in IDL
    """

    vector = vector[~np.isnan(vector)]
    if vector.size == 0:
        return(float("NaN"),float("NaN"),float("NaN"))
    else:
        clipping_apers                         = stt.sigmaclip(vector,low=threshold,high=threshold)
        objs_after_clipping                    = clipping_apers[0] #objects after trimming
        mean_objs_after_clipping               = np.nanmean(objs_after_clipping)
        stddev_objs_after_clipping             = np.nanstd(objs_after_clipping) #this is the scatter or the standard deviation
        stddev_of_the_mean_objs_after_clipping = stddev_objs_after_clipping/np.sqrt(objs_after_clipping.size-1) #this is the standard deviation of the mean
        return(mean_objs_after_clipping,stddev_objs_after_clipping,stddev_of_the_mean_objs_after_clipping)
        

def dist_ellipse(n, xc, yc, ratio, pa=0): 
    """
    original implementation (like DIST_ELLIPSE IDL function)
    
    N = either  a scalar specifying the size of the N x N square output
              array, or a 2 element vector specifying the size of the
               M x N rectangular output array.
       XC,YC - Scalars giving the position of the ellipse center.   This does
               not necessarily have to be within the image
       RATIO - Scalar giving the ratio of the major to minor axis.   This
               should be greater than 1 for position angle to have its
               standard meaning.
    OPTIONAL INPUTS:
      POS_ANG - Position angle of the major axis in degrees, measured counter-clockwise
               from the Y axis.  For an image in standard orientation
               (North up, East left) this is the astronomical position angle.
               Default is 0 degrees.
    OUTPUT:
       IM - REAL*4 elliptical mask array, of size M x N.  THe value of each
               pixel is equal to the semi-major axis of the ellipse of center
                XC,YC, axial ratio RATIO, and position angle POS_ANG, which
               passes through the pixel.
    """

    ang = np.radians(pa + 90.)
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    nx = n[1]
    ny = n[0]
    x = np.arange(-xc,nx-xc)
    y = np.arange(-yc,ny-yc)

    im = np.empty(n)
    xcosang = x*cosang
    xsinang = x*sinang

    for i in range(0, ny):

        xtemp = xcosang + y[i]*sinang
        ytemp = -xsinang + y[i]*cosang
        im[i,:] = np.sqrt((xtemp*ratio)**2 + ytemp**2)

    return im
    
    
#function that gets all the pixels within a given square aperture
#OUTPUT
#pixels_x and pixels_y: to get the all the pixels you need to do the permutations of these two variables
#flag_ok: tells you whether the squares were drawn within the limits of your image (0=no, 1=yes)
def get_pix_within_square_aper(center_x,center_y,dim1,dim2,aper_size_pix):
    #the aperture size must be divided to get the "radius", but remember it is a square
    offset = math.floor(aper_size_pix/2.)

    x0 = center_x - offset
    y0 = center_y - offset
    x1 = center_x + offset
    y1 = center_y + offset

    #I get the pixels
    pixels_x = np.arange(x0,x1)
    pixels_y = np.arange(y0,y1)

    #if the pixels do not clash with the image borders, flag_ok is equal to 1
    if np.all( [pixels_x >= 0,pixels_y >= 0,pixels_x < dim1,pixels_y < dim2] ):
        flag_ok = 1
    else:
        flag_ok = 0

    #as I will return pixels, it makes sense they are integer numbers
    pixels_x = pixels_x.astype(int)
    pixels_y = pixels_y.astype(int)

    return(flag_ok, pixels_x, pixels_y)
    
    
def get_pixs_in_circ_aper(center_x,center_y,dim1,dim2,radius):
    """
    as it contains a double for loop, it takes ages if dealing with very large apertures. do not do them
    """
    if radius > 0:
        x0 = np.round(center_x-radius).astype(int)
        x1 = np.round(center_x+radius).astype(int)
        y0 = np.round(center_y-radius).astype(int)
        y1 = np.round(center_y+radius).astype(int)
        
    squared_radius = radius**2. #to speed up things
    
    pixs_x = np.array([],dtype=int)
    pixs_y = np.array([],dtype=int)
    for ii in range(x0,x1):
        for jj in range(y0,y1):
            if ((ii-center_x)**2.)+((jj-center_y)**2.) < squared_radius and ii > 0 and ii < dim1 and jj > 0 and jj < dim2:
                pixs_x = np.append(pixs_x,ii)
                pixs_y = np.append(pixs_y,jj)
    
    return((pixs_x,pixs_y))
    

def masking_regions(img,file_with_regions,dim1,dim2):
    if os.path.isfile(file_with_regions) == True:
        with open(file_with_regions,"r") as ff:
            for line in ff:
                if line.find("box") != -1:
                    params = line.split(",")
                    xc     = float(params[0][4:])
                    yc     = float(params[1])
                    lenght = float(params[2])
                    height = float(params[3])
                    x0 = np.round(xc-lenght/2).astype(int)
                    x1 = np.round(xc+lenght/2).astype(int)
                    y0 = np.round(yc-height/2).astype(int)
                    y1 = np.round(yc+height/2).astype(int)
                    if x0 < 0: x0 = 0
                    if x1 > dim1: x1 = dim1-1
                    if y0 < 0: y0 = 0
                    if y1 > dim2: y1 = dim2-1
                    img[y0:y1,x0:x1] = 0
                if line.find("circle") != -1:
                    params = line.split(",")
                    xc     = float(params[0][7:])
                    yc     = float(params[1])
                    radius = float(params[2][:-2])
                    if radius > 0:
                        pixels_x,pixels_y = get_pixs_in_circ_aper(xc,yc,img.shape[0],img.shape[1],radius)
                        for ii in range(len(pixels_x)):
                            img[pixels_y[ii],pixels_x[ii]] = 0

    return(img)
    
                
def unmasking_regions(img,file_with_regions,dim1,dim2):
    if os.path.isfile(file_with_regions) == True:
        with open(file_with_regions,"r") as ff:
            for line in ff:
                if line.find("box") != -1:
                    params = line.split(",")
                    xc     = float(params[0][4:])
                    yc     = float(params[1])
                    lenght = float(params[2])
                    height = float(params[3])
                    x0 = np.round(xc-lenght/2).astype(int)
                    x1 = np.round(xc+lenght/2).astype(int)
                    y0 = np.round(yc-height/2).astype(int)
                    y1 = np.round(yc+height/2).astype(int)
                    if x0 < 0: x0 = 0
                    if x1 > dim1: x1 = dim1-1
                    if y0 < 0: y0 = 0
                    if y1 > dim2: y1 = dim2-1
                    img[y0:y1,x0:x1] = 1
                if line.find("circle") != -1:
                    params = line.split(",")
                    xc     = float(params[0][7:])
                    yc     = float(params[1])
                    radius = float(params[2][:-2])
                    if radius > 0:
                        pixels_x,pixels_y = get_pixs_in_circ_aper(xc,yc,img.shape[0],img.shape[1],radius)
                        for ii in range(len(pixels_x)):
                            img[pixels_y[ii],pixels_x[ii]] = 1
    return(img)
    
    
def correcting_by_extinction(mags,field,band):
    #extinction quantities from https://ned.ipac.caltech.edu/extinction_calculator
    
    if   band == 'B':
        if   (field == "goodss") == True: mags = mags - 0.0315
    elif band == 'V':
        if   (field == "goodss") == True: mags = mags - 0.0215
        elif (field == "goodsn") == True: mags = mags - 0.029
        elif (field == "cosmos") == True: mags = mags - 0.0455
        elif (field == "uds") == True:    mags = mags - 0.050
        elif (field == "egs") == True:    mags = mags - 0.0235
    elif band == 'I':
        if   (field == "goodss") == True: mags = mags - 0.013
        elif (field == "goodsn") == True: mags = mags - 0.018
        elif (field == "cosmos") == True: mags = mags - 0.028
        elif (field == "uds") == True:    mags = mags - 0.031
        elif (field == "egs") == True:    mags = mags - 0.015
    elif band == 'Y':
        if   (field == "goodss") == True: mags = mags - 0.0085
    elif band == 'J':
        if   (field == "goodss") == True: mags = mags - 0.0065
        elif (field == "goodsn") == True: mags = mags - 0.0085
        elif (field == "cosmos") == True: mags = mags - 0.0135
        elif (field == "uds") == True:    mags = mags - 0.015
        elif (field == "egs") == True:    mags = mags - 0.007
    elif band == 'H':
        if   (field == "goodss") == True: mags = mags - 0.0045
        elif (field == "goodsn") == True: mags = mags - 0.006
        elif (field == "cosmos") == True: mags = mags - 0.0095
        elif (field == "uds") == True:    mags = mags - 0.01
        elif (field == "egs") == True:    mags = mags - 0.005

    return(mags)
    
    
def correcting_by_inclination(sb_profile, axis_ratio):
    #according to the recipe in Trujillo+20's Section 5.2 
    #values for z_0/h = 0.12
    alphas = np.array([2.845,-7.833,10.792,-8.482,2.679])
    
    delta_mu = 0.
    for ii,alpha in enumerate(alphas):
        delta_mu = delta_mu + alpha*(axis_ratio**ii)
        
    sb_profile = sb_profile + delta_mu

    return(sb_profile)
    
    
def plot_options():
    #starting the plotting
    rc('axes', linewidth=3)
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
    mpl.rcParams['axes.xmargin'] = 0
    #for more params mpl.rcParams.keys()
    
    
def obtaining_random_center_for_aper(dim1,dim2):
    center_x = np.random.uniform()
    center_y = np.random.uniform()
    center_x = center_x * dim1
    center_y = center_y * dim2
    center_x = math.floor(center_x)
    center_y = math.floor(center_y)    
    return(center_x,center_y)
    
    
def calculate_limiting_mag(filename_image,filename_mask,method,aper_shape,aper_size,pix_scale,no_apers,zeropoint,axis_ratio,position_angle,max_dist_pix,centroid_x,centroid_y):
    """
    aper_size in arcsec. if it is a box, it must be the size. if it is a circle, it must be the the diameter
    """
    threshold = 3 #[sigma] for the resistant mean
    no_apers = int(no_apers)
    position_angle = position_angle + 90. #+90 to follow the share the same reference frame
    
    #reading images===================================================
    img = fits.open(filename_image)
    image  = img[0].data
    header = img[0].header

    msk = fits.open(filename_mask)
    mask     = msk[0].data
    hdr_mask = msk[0].header

    dim1 = image.shape[0]
    dim2 = image.shape[1]
    #==================================================================

    #overwriting the mask in the science image
    pixels_to_avoid = mask == 0
    image    [pixels_to_avoid] = float("NaN")
    
    #creating mask for the central galaxy
    ellip_dist = dist_ellipse([dim1,dim2],centroid_x,centroid_y,1./axis_ratio,position_angle)
    pix_central_gal = ellip_dist < max_dist_pix
    if pix_central_gal.all() == True: #if all pixels should be flag according the masking criterion, only mask the very central ones
        image[int((dim1/2)-100):int((dim1/2)+100),int((dim2/2)-100):int((dim2/2)+100)] = float("NaN")
    else:
        image[pix_central_gal] = float("NaN")
    
    sum_of_each_aper = np.array([])
    sum_of_each_aper_rms = np.array([])
    
    if method == "regular_apertures":

        if aper_shape == "box":
            #for each aperture
            for ii in range(no_apers):
                if ii%1000 == 0: print(ii) #to check the algorithm is working
                
                while True: #for entering once at least
                    
                    center_x, center_y = obtaining_random_center_for_aper(dim1,dim2)

                    #I get the aperture
                    aper_size_pix = aper_size/pix_scale #with decimals
                    flag_ok, pixels_x, pixels_y = get_pix_within_square_aper(center_x,center_y,dim1,dim2,aper_size_pix)
                    if flag_ok == 0: continue #checking whether I am not clashing with the image borders

                    #I get the values within the positions I obtained
                    pixel_values     = np.array([]) 
                    pixel_values_rms = np.array([])
                    for pixel_x in pixels_x:
                        for pixel_y in pixels_y:
                            pixel_values     = np.append(pixel_values    , image[pixel_y,pixel_x]) #OJO x and y!!! 

                    #if none of the values is NaN
                    if np.all(np.isfinite(pixel_values)):
                        #ff.write('box('+str(int(center_x)+1)+','+str(int(center_y)+1)+','+str(int(pixels_x.size))+','+str(int(pixels_y.size))+',0)\n') #OJO x and y!!! 
                        sum_of_each_aper     = np.append(sum_of_each_aper    ,np.sum(pixel_values))
                        #+1 in the coordinates because of moving from Python (first pixel is [0,0]) to ds9 (first pixel is [1,1]) reference frame
                        break
                    else:
                        #flag_exit=0
                        continue
        
        elif aper_shape == "circ":
            #for each aperture
            for ii in range(no_apers):
                if ii%1000 == 0: print(ii) #to check the algorithm is working
                
                while True: #for entering once at least
                    
                    center_x, center_y = obtaining_random_center_for_aper(dim1,dim2)

                    #I get the aperture
                    aper_size_pix = aper_size/pix_scale #with decimals
                    flag_ok, pixels_x, pixels_y = get_pix_within_circ_aper(center_x,center_y,dim1,dim2,aper_size_pix)
                    if flag_ok == 0: continue #checking whether I am not clashing with the image borders

                    #I get the values within the positions I obtained
                    pixel_values     = np.array([]) 
                    pixel_values_rms = np.array([])
                    for jj in range(pixels_x.size):
                            pixel_values     = np.append(pixel_values    , image    [pixels_y[jj],pixels_x[jj]]) #OJO x and y!!! 
                    
                    #if none of the values is NaN
                    if np.all(np.isfinite(pixel_values)):
                        #ff.write('circle('+str(int(center_x)+1)+','+str(int(center_y)+1)+','+str(math.floor(aper_size_pix/2.))+')\n') #OJO x and y!!! 
                        sum_of_each_aper     = np.append(sum_of_each_aper    ,np.sum(pixel_values))
                        #+1 in the coordinates because of moving from Python (first pixel is [0,0]) to ds9 (first pixel is [1,1]) reference frame                          
                        break
                    else:
                        #flag_exit=0
                        continue  
            
    elif method == "joining_individual_pixels_to_create_apers":
        valid_pixs     = image[~np.isnan(image)]      
        
        aper_size_pix = aper_size/pix_scale
        if aper_shape == "box" : no_pixs_to_take = np.floor( aper_size_pix**2. )
        if aper_shape == "circ": no_pixs_to_take = np.floor( np.pi*((aper_size_pix/2.)**2.) )
        no_pixs_to_take = int(no_pixs_to_take)
        
        pixel_values     = np.array([]) 
        for ii in range(no_apers):
            values_pixs_in_virtual_aperture     = np.random.choice(valid_pixs    ,no_pixs_to_take,replace=True)
            sum_of_each_aper    = np.append(sum_of_each_aper    ,np.sum(values_pixs_in_virtual_aperture))

    mean_value_apers,std_value_apers,std_value_of_the_mean_apers = resistant_mean(sum_of_each_aper,threshold)
    SB_limit         = -2.5*np.log10(std_value_apers)+zeropoint+5.*np.log10(aper_size)
    """
    how Nacho does it?
    if aper_size < 1 arcsec
    SB_limit_standard_nacho = -2.5*np.log10(3.*std_value_apers)+5.*np.log10(aper_size)+zeropoint[cont_band] #3sigma in 10x10 arsec boxes
    if aper_size > 1 arcsec
    SB_limit_standard_nacho = -2.5*np.log10(3.*std_value_apers)+2.5*np.log10(aper_size)+zeropoint[cont_band] #3sigma in 10x10 arsec boxes
    because the area increases as aper_size*aper_size but the signal_to_noise only increases as aper_size
    """
    
    return(SB_limit) #3sigma add -2.5*np.log10(3), 5sigma add -2.5*np.log10(5)
    
    
def calculate_sb_profile(distance_vector,centroid,axis_ratio,position_angle,img,mask,pix_scale,zeropoint,gal,rms=None,check=False,only_flux=False):
    """
    INPUT
    distance_vector: follows the semi-major axis, in arcsec
    centroid: in pixels
    axis_ratio: I assume it is in the GALFIT reference frame
    position_angle: I assume it is in the GALFIT reference frame
    img: the galaxy image
    mask: the galaxy mask
    pix_scale: the pixel scale in arcsec/pix
    rms: the rms image
    check: whether I like to create a fits image with the outcome of dist_ellipse to check the axis ratio and position angle
    only_flux: the return value corresponds to the fluxes in the apertures, not the surface brightnesses
    OUTPUT
    surf_bright: the surface brightness profile
    surf_bright_error: the eror in the surface brightness profile
    """
    threshold = 3.
    
    if np.isnan(img.shape[0]) or np.isnan(img.shape[1]) or np.isnan(centroid[0]) or np.isnan(centroid[1]) or np.isnan(axis_ratio) or np.isnan(position_angle):
        if only_flux == False:
            surf_bright       = np.full_like(distance_vector,float("NaN"))
            surf_bright_error = np.full_like(distance_vector,float("NaN"))
        else:
            flux       = np.full_like(distance_vector,float("NaN"))
            sigma_flux = np.full_like(distance_vector,float("NaN"))    
            
    else:

        resolution = distance_vector.size
        position_angle = position_angle + 90. #+90 to follow the share the same reference frame
        ellip_dist = dist_ellipse([img.shape[0],img.shape[1]],centroid[0],centroid[1],1./axis_ratio,position_angle)
        if check == True: fits.writeto("dist_ellipse.fits",ellip_dist,overwrite=True) #draw the ellipses to confirm their axis ratio and position angle

        #not taking into account the pixels in the mask 
        pix_to_remove1 = mask == 0
        pix_to_remove2 = img  == 0
        pix_to_remove = np.logical_or(pix_to_remove1,pix_to_remove2)
        mask = np.array(mask,dtype=float) #necessary for adding the NaNs later on
        mask[pix_to_remove] = float("NaN") #for the mask to have NaNs and not zeros (to avoid confusions)

        flux       = np.array([])
        sigma_flux = np.array([])
        for ii in range(resolution):
            if ii == 0:
                ellip_annulus = np.logical_and(ellip_dist > 0.,                              ellip_dist < distance_vector[ii]/pix_scale)
            else:
                ellip_annulus = np.logical_and(ellip_dist > distance_vector[ii-1]/pix_scale, ellip_dist < distance_vector[ii]/pix_scale)

            image_in_annulus = img[ellip_annulus]*mask[ellip_annulus]
            image_in_annulus = image_in_annulus[~np.isnan(image_in_annulus)] #remove all NaNs
                
            mean_aperture,median_aperture,stddev_aperture = sigma_clipped_stats(image_in_annulus,sigma=threshold)
            flux = np.append(flux,mean_aperture)
            if rms is None:
                sigma_flux = np.append(sigma_flux,stddev_aperture/np.sqrt(image_in_annulus.size-1))
            else:
                rms_in_annulus = rms[ellip_annulus]*mask[ellip_annulus]
                rms_in_annulus = rms_in_annulus[~np.isnan(rms_in_annulus)] #remove all NaNs
                none1,none2,stddev_aperture = resistant_mean(image_in_annulus,threshold)
                none1,none2,stddev_rms      = resistant_mean(rms_in_annulus,  threshold)
                sigma_total = np.sqrt( (stddev_aperture**2.)+(stddev_rms**2.) ) #sigma_total = sqrt(sigma_aperture^2 + sigma_rms^2); regarding my discussion with Alex Borlaff
                sigma_flux = np.append(sigma_flux,sigma_total)

        surf_bright       =             -2.5*np.log10(  flux             )+zeropoint+5.*np.log10(pix_scale)
        surf_bright_error = np.absolute(-2.5*np.log10( (flux+sigma_flux) )+zeropoint+5.*np.log10(pix_scale) - surf_bright)
        
    if only_flux == False:
        return(surf_bright,surf_bright_error) #normal usage
    else:
        return(flux,sigma_flux)
        
    
def calculate_mass_profile(sb_profile,color,error_in_color,filter_name,color_name):
    
    #from Roediger & Courteau et al. (2015)=======
    mm = {"g-r": {"g": 2.029, "r": 1.629, "i": 1.438, "z": 1.306}, 
          "g-i": {"g": 1.379, "r": 1.110, "i": 0.979, "z": 0.886},
          "g-z": {"g": 1.116, "r": 0.900, "i": 0.793, "z": 0.716},
          "r-i": {"g": 4.107, "r": 3.325, "i": 2.925, "z": 2.634},
          "r-z": {"g": 2.322, "r": 1.883, "i": 1.655, "z": 1.483},
          "i-z": {"g": 5.164, "r": 4.201, "i": 3.683, "z": 3.283}}
    bb = {"g-r": {"g": -0.984, "r": -0.792, "i": -0.771, "z": -0.796}, 
          "g-i": {"g": -1.067, "r": -0.861, "i": -0.831, "z": -0.848},
          "g-z": {"g": -1.132, "r": -0.916, "i": -0.878, "z": -0.888},
          "r-i": {"g": -1.170, "r": -0.952, "i": -0.908, "z": -0.912},
          "r-z": {"g": -1.211, "r": -0.987, "i": -0.937, "z": -0.935},
          "i-z": {"g": -1.212, "r": -0.991, "i": -0.939, "z": -0.931}} 
    #===============================
    
    #Sun's absolute mag according to http://mips.as.arizona.edu/~cnaw/sun.html
    M_Sun = {"g": 5.11,"r": 4.65,"i": 4.53,"z": 4.50}
        
    if np.isnan(sb_profile[0]) == False: 

        color_with_errors = color + error_in_color

        log_M_over_L       = bb[color_name][filter_name] + (mm[color_name][filter_name] * color            ) -0.15  #0.15 for obtaining it in Kroupa IMF
        error_log_M_over_L = bb[color_name][filter_name] + (mm[color_name][filter_name] * color_with_errors) - 0.15 #0.15 for obtaining it in Kroupa IMF
        
        log_mass_profile       =       log_M_over_L - 0.4*(sb_profile-M_Sun[filter_name]) + 8.629 + 6. #6. from pc^2 to Kpc^2
        error_log_mass_profile = error_log_M_over_L - 0.4*(sb_profile-M_Sun[filter_name]) + 8.629 + 6. #6. from pc^2 to Kpc^2
        log_mass_profile       = log_mass_profile - 0.04       #to move from Kroupa to Chabrier mass according to Bruce+12
        error_log_mass_profile = error_log_mass_profile - 0.04 #to move from Kroupa to Chabrier mass according to Bruce+12

        return(10**log_mass_profile,np.abs( (10**log_mass_profile)-(10**error_log_mass_profile) ),10**log_M_over_L,10**error_log_M_over_L)
    else:
        return(float("NaN"),        float("NaN"),                                                 float("NaN"),    float("NaN"))
