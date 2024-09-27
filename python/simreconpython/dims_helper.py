import numpy as np
from skimage.transform import resize

def reshape_zdpyx_dpzyx(raw, nphases=5, ndirections=3):
    """
    Reshape z-dir-phase-y-x (fastest changing index to slowest) data into 
    dir-phase-z-y-x (fastest changing index to slowest) data.

    Parameters:
    raw (np.ndarray): The raw input data array.
    nphases (int): Number of phases. Default is 5.
    ndirections (int): Number of directions. Default is 3.

    Returns:
    np.ndarray: Reshaped 5D array with dimensions (ndirections, nphases, nz, height, width).
    """
    nz = raw.shape[0] // (nphases * ndirections)
    images = np.zeros((ndirections, nphases, nz, raw.shape[1], raw.shape[2]))

    for d in range(ndirections):
        for p in range(nphases):
            n = d * nphases + p
            blob = raw[n::15, :, :]
            images[d, p, :, :, :] = blob

    return images

def reshape_tile_data_to_match_raw_data(tile_data, nz, ny, nx, ndirs, nphases):
    """
    Each entry in the tile data dictation is a 2d array representing the parameter values for each tile
    this function is used to make a 3d array with the same shape as the processed data, where each tile is 
    filled with the corresponding parameter values so we can visualize the parameter values overlayed with the raw data
    
    Parameters:
    tile_data (np.ndarray): The tile data array.
    nz (int): Number of z slices.
    ny (int): Number of pixels in y direction.
    nx (int): Number of pixels in x direction.
    ndirs (int): Number of directions.
    nphases (int): Number of phases.
    """

    # first reshape xy to match the size of processed data, use order 0, so we don't interpolate, each parameter value is
    #repeated for all pixels in the tile
    temp2 = resize(tile_data, (tile_data.shape[0], ny, nx),order=0,anti_aliasing = True)
    temp2.shape

    # for each direction replicate the 2d array for each z slice
    temp3 = [np.stack([temp2[dir,:,:]]*nz, axis=0) for dir in range(ndirs)]
    temp3 = np.array(temp3)

    # for each direction replicate the 3d array for each phase
    temp3.shape
    temp4 = [np.stack([temp3[dir,:,:,:]]*nphases, axis=0) for dir in range(3)]
    temp4 = np.array(temp4)
    temp4.shape

    return temp4