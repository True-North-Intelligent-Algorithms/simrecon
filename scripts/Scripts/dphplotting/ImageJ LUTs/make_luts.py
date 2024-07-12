# !/usr/bin/env python

'''
File to convert matplotlib cmaps into ImageJ LUTs
'''

import numpy as np
from matplotlib.cm import get_cmap

cmaps_to_convert = ('viridis', 'inferno', 'plasma', 'magma')


def make_cmap(cname):
    '''
    take a valid cmap name and convert and save it to a LUT
    '''
    # get the cmap
    cmap = get_cmap(cname)
    # pull the data
    data = np.array(cmap.colors)
    # convert to 8-bit RGB
    scaled_data = (data*255).astype(np.uint8)
    # save the data
    scaled_data.T.tofile(cname.capitalize()+'.lut')

if __name__ == "__main__":
    # run from command line?
    for cname in cmaps_to_convert:
        make_cmap(cname)
