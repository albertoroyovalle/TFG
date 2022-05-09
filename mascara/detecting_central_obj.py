import numpy as np

def detecting_central_obj(sex_output,x_central,y_central):

    if sex_output.ndim > 1:
        objs         = sex_output[0, :]
        x_sextractor = sex_output[5, :]
        y_sextractor = sex_output[6, :]
    else:
        objs         = np.array([sex_output[0]])
        x_sextractor = np.array([sex_output[5]])
        y_sextractor = np.array([sex_output[6]])        

    #detecting central obj
    dists = np.sqrt( ((x_central-x_sextractor)**2.) + ((y_central-y_sextractor)**2.) )
    filter_dists = dists == dists.min()
    target = objs[filter_dists]

    return(int(target)-1) #-1 because SExtractor first element is 1 and Python first element is 0
