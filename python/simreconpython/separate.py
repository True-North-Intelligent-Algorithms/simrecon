import numpy as np

def get_volume(data, phase, dir, num_phases, num_directions):
    """
    
    """
    print(dir*num_phases+phase)
    volume = data[dir*num_phases+phase::num_phases*num_directions,:,:].copy()

    return volume

def makematrix(nphases, norders):
    """
    
    """
    print('make matrix')

    phi = 2*np.pi/nphases

    sepmatrix = np.zeros((2*norders-1, nphases))

    for j in range(nphases):
      sepmatrix[0, j] = 1.0
      for order in range(1, norders):
        sepmatrix[2*order-1, j] = np.cos(j*order*phi)
        sepmatrix[2*order, j] = np.sin(j*order*phi)

    return sepmatrix

def apodize_2d(image, napodize):
    image_out = image.copy()
    ny, nx = image.shape
    for k in range(nx):
        diff = (image[ny-1, k] - image[0, k]) / 2
        for l in range(napodize):
            fact = 1 - np.sin(((l + 0.5) / napodize) * np.pi * 0.5)
            image_out[l, k] += diff * fact
            image_out[ny-1-l, k] -= diff * fact
    for l in range(ny):
        diff = (image[l, nx-1] - image[l, 0]) / 2
        for k in range(napodize):
            fact = 1 - np.sin(((k + 0.5) / napodize) * np.pi * 0.5)
            image_out[l, k] += diff * fact
            image_out[l, nx-1-k] -= diff * fact
    return image_out  # Optionally return the modified image as a flat array

def apodize_z(image, napodize):
    image_out = image.copy()
    nz, ny, nx = image.shape
    
    for k in range(nx):
        for l in range(ny):
            diff = (image[nz-1, l, k] - image[0, l, k])/2
            for m in range(napodize):
                delta = 0.01
                fact = (1- (napodize - m) / napodize)+delta
                image_out[m, l, k] = fact * image[m, l, k]
                image_out[nz-1-m, l, k] = fact * image[nz-1-m, l, k]
    return image_out  # Optionally return the modified image as a flat array

def apodize_3d(image, napodize):
    temp = image.copy().astype('float32')
    for n in range(temp.shape[0]):
        temp[n,:,:] = apodize_2d(temp[n,:,:], napodize)

    return temp 

def normalize(image):
    temp = image.copy()
    sum1 = temp[0,:,:].sum()
    for n in range(temp.shape[0]):
        sumn = temp[n,:,:].sum()
        temp[n,:,:] = temp[n,:,:] * sum1 / sumn

    return temp

def separate(matrix, data):
    """
    
    """
    separated = np.einsum('ij,jkl->ikl', matrix, data)
    return separated
