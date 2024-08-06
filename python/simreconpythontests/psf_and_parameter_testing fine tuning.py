import os
import sys
import mrc as Mrc
from tnia.plotting.projections import show_xyz_slice_center, show_xyz_slice
sys.path.insert(1, r'./python')
from simreconpython.fft_helper import show_centered_fft_magnitude, centered_fft_magnitude
from simreconpython.otf_helper import get_otf_figures
import numpy as np
from tnia.nd.ndutil import centercrop
import matplotlib.pyplot as plt

computer = 'bnort'

if computer == 'default':
    sys.path.insert(1, 'Y:\Cryo_data2\Data Processing Notebooks')
    sys.path.insert(1, 'Y:\Cryo_data2\Data Processing Notebooks\Scripts')
elif computer == 'bnort':
    sys.path.insert(1, r'C:\Users\bnort\work\Janelia\code\\simrecon\scripts\Scripts')
    sys.path.insert(1, r'C:\Users\bnort\work\Janelia\code\\simrecon\scripts')
else:
    pass

from simrecon_utils import simrecon, split_process_recombine

def psf_and_parameter_testing():
    print(os.getcwd())
    print("PSF and parameter testing")
    
    #input_name = r'D:\Janelia\Data 2024-06-06\Wiener, gammaApo and SupressR parameter testing\488nm comparison Brian\experiments_July_9\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    input_name = r'D:\Janelia\Data 2024-06-30\488cm cell 5 good signal_20240627_131236  Fail fine tuning\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #input_name = r'D:\Janelia\Data 2024-06-30\NA Test\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #input_name = r'D:\Janelia\Data 2024-07-13\makemodel\makemodel.mrc'
    #input_name = r'D:\Janelia\Data 2024-06-30\488nm RPE-1 cell 2 2_20240626_165005  Fail\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    otf_name = r'C:\Users\bnort\work\Janelia\ims\computed_OTF_folder\488nmLinOTF0_mask.mrc'
    #otf_name = r'C:\Users\bnort\work\Janelia\ims\computed_OTF_folder\488 OTF Bead 8_20190919_141256.mrc'
    wl = 488

    #input_name = r'D:\Janelia\Data 2024-06-12\parameter_testing_july15\560 nm 655 40 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #input_name = r'D:\Janelia\Data 2024-06-10\parameter_test_July16\560 nm 655 40 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #otf_name = r'D:\Janelia\Data 2024-06-12\561nm OTF used\560 201909_19-20_best.mrc'

    #input_name = r'D:\Janelia\test_data\raw.mrc'
    #otf_name = r'D:\Janelia\test_data\otf.mrc'
    z_display_slice = 10


    # set default parameters
    base_kwargs = dict(
                    nphases=5,
                    ndirs=3,
                    angle0= 1.29,
                    negDangle=True,              
                    na= 0.85,
                    nimm= 1.0,
                    zoomfact= 2.0, 
                    background= 100.0,           
                    wiener= 0.001,
                    fastSIM=True,
                    otfRA= True,
                    dampenOrder0=True,
                    k0searchall=True,
                    equalizez=True,
                    preciseapo=True,
                    gammaApo=0.7,
                    suppressR=15.0,
                    nthreads = 8
                )
    gammaApo = 0.0
    suppressR = 1.0
    wiener = 0.001
    driftfix = False 
    nosuppress = False
    zpadto = 0 
    noOrder0 = False
    makemodel = False
    tiled = False
    tile_size = 32
    tile_overlap = 16
    apodizeoutput = 1
    #na = 1.2

    for gammaApo in [0.0311]:#[20.,2.,0.2, 0.02]:#[2, 0.2, 0.02]:
        for suppressR in [3.54]:#[10, 1, 0.1]:
            for wiener in [0.00363]:#[0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]:
                user_text = 'gApo_'+str(gammaApo)+'_supR_'+str(suppressR)+'_w_'+str(wiener)+'_wl_'+str(wl)

                psf_base_name = os.path.basename(otf_name)
                psf_base_name = os.path.splitext(psf_base_name)[0]

                user_text += '_'+psf_base_name

                if zpadto > 0:
                    user_text += '_zpadto_'+str(zpadto)

                if driftfix:
                    user_text += '_driftfix'

                if nosuppress:
                    user_text += '_nosuppress'

                if noOrder0:
                    user_text += '_noOrder0'
                
                if makemodel:
                    user_text += '_model'

                if tiled:
                    user_text += '_t_'+str(tile_size)+'_'+str(tile_overlap)

                if apodizeoutput == 2:
                    user_text += '_triangle'

                print(user_text)
                
                base_kwargs.update(dict(apodizeoutput=apodizeoutput, gammaApo=gammaApo, suppressR=suppressR, wiener=wiener, nosuppress=nosuppress, driftfix = driftfix, noOrder0=noOrder0, zpadto=zpadto, makemodel=makemodel))   # default Full frame Recon. parameters    
                sim_kwargs = dict(                                                                                                            
                input_file= input_name,
                otf_file= otf_name,
                ls= (wl/1000)/2/0.81)
                #ls = 0.2035)
                    
                sim_kwargs.update(base_kwargs)
                            
                #create processed file output name
                output_file = sim_kwargs["input_file"].replace(".mrc", '_proc'+'_' + user_text + ".mrc")
                sim_kwargs["output_file"] = output_file

                otf_mrc = Mrc.Mrc(otf_name)
                otf_data = otf_mrc.data
                otf_fig, otf_log_fig, psf_fig = get_otf_figures(otf_data)

                # does output file already exist? 
                if False: #os.path.exists(output_file):
                    print("Output file already exists.  Figures will be generated on existing file.")
                else:
                    if tiled: 
                        output_file, sim_output = split_process_recombine(sim_kwargs["input_file"], tile_size, tile_overlap, sim_kwargs)
                    else:
                        sim_output = simrecon(**sim_kwargs)
                    print(sim_output)
                    for i in range(len(sim_output)):
                        sim_output[i] = sim_output[i].replace('\r', '')
                    with open(output_file.replace(".mrc", ".txt"), "w") as myfile:
                        myfile.write(str(sim_kwargs))
                        myfile.write("\n".join(sim_output))
                    

                input_path = os.path.dirname(input_name)
                figure_path = os.path.join(input_path, user_text)

                if not os.path.exists(figure_path):
                    os.makedirs(figure_path)

                otf_fig.savefig(os.path.join(figure_path, user_text + "_otf.png"))
                otf_log_fig.savefig(os.path.join(figure_path, user_text + "_otf_log.png"))
                psf_fig.savefig(os.path.join(figure_path, user_text + "_psf.png"))    

                reconstructed_mrc = Mrc.Mrc(output_file)
                reconstructed = reconstructed_mrc.data
                print(reconstructed.shape)


                # Assuming 'reconstructed' is your data array
                # Create the histogram
                fig, ax = plt.subplots()
                ax.hist(reconstructed.flatten(), bins=256, color='blue', alpha=0.7)
                ax.set_title('Histogram of Reconstructed Data')
                ax.set_xlabel('Pixel Intensity')
                ax.set_ylabel('Frequency')

                # Save the histogram figure
                histogram_path = os.path.join(figure_path, user_text + "_reconstructed_histogram.png")
                fig.savefig(histogram_path)

                #fig = show_xyz_slice_center(reconstructed, sxy=1, sz=20)
                fig = show_xyz_slice(reconstructed, reconstructed.shape[2]//2, reconstructed.shape[1]//2, z_display_slice, sxy=1, sz=5)
                print(reconstructed.min(), reconstructed.max())
                fig.savefig(os.path.join(figure_path, user_text + "_xyz_center.png"))

                view = True

                yc = 697
                xc = 843
                reconstructed_crop = reconstructed[:, max(yc-512,0):min(yc+512, reconstructed.shape[1]-1), max(xc-512,0):min(xc+512, reconstructed.shape[2]-1)]
                fig = show_xyz_slice_center(reconstructed_crop, sxy=1, sz=5)
                fig.savefig(os.path.join(figure_path, user_text + "_xyz_center_crop.png"))
                
                yc = 527
                xc = 773
                size = 256
                reconstructed_crop = reconstructed[:, max(yc-size,0):min(yc+size, reconstructed.shape[1]-1), max(xc-size,0):min(xc+size, reconstructed.shape[2]-1)]
                fig = show_xyz_slice_center(reconstructed_crop, sxy=1, sz=5)
                fig.savefig(os.path.join(figure_path, user_text + "_xyz_center_crop2.png"))

     
                yc = 527
                xc = 818
                size = 128
                reconstructed_crop = reconstructed[:, max(yc-size,0):min(yc+size, reconstructed.shape[1]-1), max(xc-size,0):min(xc+size, reconstructed.shape[2]-1)]
                fig = show_xyz_slice_center(reconstructed_crop, sxy=1, sz=5)
                fig.savefig(os.path.join(figure_path, user_text + "_xyz_center_crop3.png"))
                fig.savefig(os.path.join(input_path, 'wiener_gA1.0_triangle', user_text + "_xyz_center_crop3.png"))

                center_line = reconstructed[:, yc, xc]
                # Create a new figure and axis
                fig, ax = plt.subplots()

                # Plot the center line
                ax.plot(center_line)

                # Optionally, set labels and title
                ax.set_xlabel('Position along the line')
                ax.set_ylabel('Intensity')
                ax.set_title('Center Line Intensity Profile')

                # Save the figure
                fig.savefig(os.path.join(figure_path, user_text + "_center_line_profile.png"))
                
                fft = np.fft.fftn(reconstructed)
                fft = np.fft.fftshift(fft)
                fft_app_mag = np.abs(fft)
                
                fft_app_mag_log = np.log(fft_app_mag + 1)

                xc = fft_app_mag.shape[1]//2
                yc = fft_app_mag.shape[2]//2
                #fft_app_mag[:,xc-10:xc+10,yc-10:yc+10] = 0
                fig = show_xyz_slice_center(fft_app_mag, sxy=1, sz=20, gamma=0.2)
                fig.savefig(os.path.join(figure_path, user_text + "_fft_center_mag.png"))

                fig = show_xyz_slice_center(fft_app_mag_log, sxy=1, sz=20)
                fig.savefig(os.path.join(figure_path, user_text + "_fft_center.png"))

if __name__ == "__main__":
    psf_and_parameter_testing()