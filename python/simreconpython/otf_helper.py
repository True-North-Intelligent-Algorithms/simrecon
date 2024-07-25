import numpy as np
from tnia.plotting.plt_helper import imshow_multi2d

def get_otf_figures(otf_data):
    print('get_otf_figures has recieved otf of shape',otf_data.shape) 
    otf = [np.transpose(otf_data[i,:,:] ) for i in range(otf_data.shape[0])]
    otf_abs = [np.abs(otf[i]) for i in range(otf_data.shape[0])]
    fig = imshow_multi2d(otf_abs, ['otf1', 'otf2', 'otf3'], 1 ,3)
    
    real_shape = ( otf[0].shape[0], (otf[0].shape[1]-1)*2)
    psfs = [np.fft.irfftn(o, s=real_shape) for o in otf]
    psfs = [np.transpose(psf) for psf in psfs]
    otfs = [np.fft.fftn(psf) for psf in psfs]
    otfs_abs = [np.abs(np.fft.fftshift(otfs[i])) for i in range(otf_data.shape[0])]
    delta = 0
    for o in otfs_abs: o[o<delta] = delta 
    otfs_log = [np.log(otfs_abs[i]) for i in range(otf_data.shape[0])]
    psfs_shift = [np.fft.fftshift(psf) for psf in psfs]
    
    fig_otfs = imshow_multi2d(otfs_abs, [f'$otf_0$', f'$otf_1$', f'$otf_2$'], 1 ,3, xlabels = [f'$k_z$', '$k_z$', '$k_z$'], ylabels = [f'$k_r$', '$k_r$', '$k_r$'], height = 5)
    fig_otfs.suptitle('Separated OTF order 0 to 2')
    #fig.savefig(r'D:\Janelia\Data 2024-06-12\band_figures\separated_otfs.png')

    fig_otfs_log = imshow_multi2d(otfs_log, [f'$otf_0$', f'$otf_1$', f'$otf_2$'], 1 ,3, xlabels = [f'$k_z$', '$k_z$', '$k_z$'], ylabels = [f'$k_r$', '$k_r$', '$k_r$'], height = 5)
    fig_otfs_log.suptitle('Separated OTF order 0 to 2')
    
    fig_psfs = imshow_multi2d(psfs_shift, ['$psf_0$', '$psf_1$', '$psf_2'], 1 ,3, xlabels = [f'$z$', '$z$', '$z$'], ylabels = [f'$r$', '$r$', '$r$'], height = 5)
    fig_psfs.suptitle('PSF order 0 to 2')
    #fig.savefig(r'D:\Janelia\Data 2024-06-12\band_figures\separated_psfs.png')

    return fig_otfs, fig_otfs_log, fig_psfs
