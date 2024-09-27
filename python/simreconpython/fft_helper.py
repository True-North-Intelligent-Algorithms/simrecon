import numpy as np
from tnia.plotting.projections import show_xyz_slice_center

def centered_fft_magnitude(im, log):
    """Compute the magnitude of the FFT of an image, centered around the
    origin.   The magnitude is returned as a 2D array.

    Parameters
    ----------
    im : 2D array
        The image to be transformed.
    log : bool
        If True, the magnitude is returned as log(1 + magnitude).

    Returns
    -------
    2D array
        The magnitude of the FFT of the image.
    """
    fft_im = np.fft.fftshift(np.fft.fftn(im))
    # Compute the magnitude of the FFT.
    magnitude = np.abs(fft_im)
    if log:
        magnitude = np.log1p(magnitude)
    return magnitude

def show_centered_fft_magnitude(im, log=True):
    """Compute and display the magnitude of the FFT of an image, centered
    around the origin. The magnitude is displayed as a 2D array.

    Parameters
    ----------
    im : 2D array
        The image to be transformed.
    log : bool
        If True, the magnitude is displayed as log(1 + magnitude).
    """
    magnitude = centered_fft_magnitude(im, log)
    fig = show_xyz_slice_center(magnitude)
    return fig
