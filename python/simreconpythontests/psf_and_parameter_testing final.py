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

from simrecon_utils import simrecon, split_process_recombine, process_txt_output, plot_params

def psf_and_parameter_testing():
    print(os.getcwd())
    print("PSF and parameter testing")
    
    #input_name = r'D:\Janelia\Data 2024-06-30\488cm cell 5 good signal_20240627_131236  Fail vector tuning zcrop experiment\cropped\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #input_name = r'D:\Janelia\Data 2024-06-30\488cm cell 5 good signal_20240627_131236  Fail vector tuning Best\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    input_name = r'D:\Janelia\Data 2024-06-30\filteroverlaps5\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #input_name = r'D:\Janelia\Data 2024-06-30\488nm RPE-1 cell 2 2_20240626_165005  Fail zcrop experiment\cropped\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    z_display_slice = 9 
    otf_name = r'C:\Users\bnort\work\Janelia\ims\computed_OTF_folder\488nmLinOTF0_mask.mrc'
    wl = 488

    #input_name = r'D:\Janelia\Data 2024-08-22\561nm cell 5 good signal_20240627_151300\560 nm 615 45 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'

    #otf_name = r'C:\Users\bnort\work\Janelia\otfs\computed_OTF_folder\532 OTF Bead 3_20190920_154920.mrc'

    #input_name = r'D:\Janelia\Data 2024-09-03\cell4\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #otf_name = r'C:\Users\bnort\work\Janelia\ims\computed_OTF_folder\488nmLinOTF0_mask.mrc'
    #wl = 488

    #input_name = r'D:\Janelia\Data 2024-09-06\LID701_ROI4_20231120_162836\560 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #otf_name = r'D:\Janelia\Data 2024-08-22\BEST PSFs\560 201909_19-20_best.mrc'
    #wl = 560

    #input_name = r'D:\Janelia\Data 2024-06-30\cell 2 Sep 23\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    input_name = r'D:\Janelia\Data 2024-06-30\cell 5 Oct 8 crop\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'

    #initial_results_name = r'D:\Janelia\Data 2024-06-30\cell 5 Oct 1\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.2_supR_0.01_w_0.0005_wl_488_488nmLinOTF0_mask_t_16_8_forceotfamp_1_1_tile16_pad8.txt'
    #initial_results_name = r'D:\Janelia\Data 2024-06-30\cell 5 Oct 1\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.3_supR_0.1_w_0.0005_wl_488_488nmLinOTF0_mask_t_8_4_forceotfamp_1_1_tile8_pad4.txt'
    
    #initial_results_name = r'D:\Janelia\Data 2024-06-30\cell 5 Oct 1\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.2_supR_0.1_w_0.0005_wl_488_488nmLinOTF0_mask_t_32_16_forceotfamp_1_1_tile32_pad16.txt'
    #input_name = r'D:\Janelia\Data 2024-06-30\cell 5 Oct 1 cropped\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #initial_results_name = r'D:\Janelia\Data 2024-06-30\cell 5 Oct 1 cropped\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.2_supR_0.5_w_0.0001_wl_488_488nmLinOTF0_mask_t_16_8_forceotfamp_1_1_tile16_pad8.txt'

    #input_name = r'D:\Janelia\Data 2024-09-03\cell 4 oct 2\560 nm 615 45 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #initial_results_name = r'D:\Janelia\Data 2024-09-03\cell 4 oct 2\560 nm 615 45 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.3_supR_1.0_w_0.0025_wl_488_560 201909_19-20_best_t_16_8_forceotfamp_1_1_tile16_pad8.txt'
    #otf_name = r'D:\Janelia\Data 2024-08-22\BEST PSFs\560 201909_19-20_best.mrc'

    #input_name = r'D:\Janelia\Data 2024-10-02\560cm cell 4 _20240627_124604\560 nm 615 45 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    #otf_name = r'D:\Janelia\Data 2024-08-22\BEST PSFs\560 201909_19-20_best.mrc'
    #initial_results_name = r'D:\Janelia\Data 2024-06-30\cell 5 Oct 1\560 nm 615 45 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.3_supR_0.1_w_0.001_wl_488_560 201909_19-20_best_t_16_8_forceotfamp_1_1_tile16_pad8.txt'

    input_name = r'D:\Janelia\summary_experiment\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'

    input_name = r'D:\Janelia\Data 2024-11-21\position 3 nice cells 561nm 4.74um stack\560 nm 655 40 filter 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc'
    otf_name = r'D:\Janelia\Data 2024-08-22\BEST PSFs\560 201909_19-20_best.mrc'
    wl = 560

    input_name = r'D:\Janelia\Data 2024-11-21\view_overlaps\488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1.mrc' 
    otf_name = r'C:\Users\bnort\work\Janelia\ims\computed_OTF_folder\488nmLinOTF0_mask.mrc'
    wl = 488
    
    # set default parameters
    base_kwargs = dict(
                    nphases=5,
                    ndirs=3,
                    angle0= 1.36,
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
    gammaApo = 0.5
    suppressR = 2.0
    wiener = 0.00025
    driftfix = False 
    nosuppress = False
    zpadto = 0 
    noOrder0 = False
    makemodel = False
    
    tiled = False
    
    tile_size = 16 
    tile_overlap = 8
    filter_tiles = False
    keeporder2 = False
    
    apodizeoutput = 1
    nofilteroverlaps = False 
    forcemodamp = False
    o1 = 0.8
    o2 = 0.3
    nok0search = False

    for gammaApo in [0.5]:#[20.,2.,0.2, 0.02]:#[2, 0.2, 0.02]:
        for suppressR in [2.0]:#[10, 1, 0.1]:
            for wiener in [0.0005]:#[0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]:
                
                user_text = 'order1_gApo_'+str(gammaApo)+'_supR_'+str(suppressR)+'_w_'+str(wiener)+'_wl_'+str(wl)

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
                elif apodizeoutput == 0:
                    user_text += '_noapo'

                if nofilteroverlaps==True:
                    user_text += '_nofilteroverlaps'

                if forcemodamp:
                    user_text += '_forcemodamp_'+str(o1)+'_'+str(o2)
                    forcemodamp = [o1, o2]
               
                if keeporder2:
                    user_text += '_keeporder2'

                if nok0search:
                    user_text=user_text+'_nok0search'

                if filter_tiles:
                    user_text=user_text+'_filter_tile2'

                    tile_limits = {}
                    tile_limits['spacing_min'] = 0.355
                    tile_limits['spacing_max'] = 0.365
                    tile_limits['spacing_default']=0.361
                    tile_limits['angle_min'] = 1.3
                    tile_limits['angle_max'] = 1.4
                    tile_limits['angle_default'] = 1.36
                    tile_limits['amp1_min'] = 0.98
                    tile_limits['amp1_max'] = 1.02
                    tile_limits['amp1_default'] = 1.0
                    tile_limits['amp2_min'] = 0.98
                    tile_limits['amp2_max'] = 1.02
                    tile_limits['amp2_default'] = 1.0
                    
                    user_text=user_text+'_filter_tile'
                    
                    tile_limits['spacing_min'] = 0.34
                    tile_limits['spacing_max'] = 0.38
                    tile_limits['spacing_default']=0.361
                    tile_limits['angle_min'] = 1.0
                    tile_limits['angle_max'] = 2.0
                    tile_limits['angle_default'] = 1.36
                    tile_limits['amp1_min'] = 0.98
                    tile_limits['amp1_max'] = 1.02
                    tile_limits['amp1_default'] = 1.0
                    tile_limits['amp2_min'] = 0.98
                    tile_limits['amp2_max'] = 1.02
                    tile_limits['amp2_default'] = 1.0
                else:
                    tile_limits = None

                print(user_text)
                
                if forcemodamp:
                    #base_kwargs.update(dict(nofilteroverlaps=nofilteroverlaps, apodizeoutput=apodizeoutput, gammaApo=gammaApo, suppressR=suppressR, wiener=wiener, nosuppress=nosuppress, driftfix = driftfix, noOrder0=noOrder0, zpadto=zpadto, makemodel=makemodel, otfamp=[1,1], forcemodamp=forcemodamp, nok0search=nok0search))   # default Full frame Recon. parameters    
                    base_kwargs.update(dict(nofilteroverlaps=nofilteroverlaps, apodizeoutput=apodizeoutput, gammaApo=gammaApo, suppressR=suppressR, wiener=wiener, nosuppress=nosuppress, driftfix = driftfix, noOrder0=noOrder0, zpadto=zpadto, makemodel=makemodel, forcemodamp=forcemodamp, nok0search=nok0search))   # default Full frame Recon. parameters    
                #elif forceotfamp:
                #    base_kwargs.update(dict(nofilteroverlaps=nofilteroverlaps, apodizeoutput=apodizeoutput, gammaApo=gammaApo, suppressR=suppressR, wiener=wiener, nosuppress=nosuppress, driftfix = driftfix, noOrder0=noOrder0, zpadto=zpadto, makemodel=makemodel, otfamp=forceotfamp))
                else:
                    base_kwargs.update(dict(nofilteroverlaps=nofilteroverlaps, apodizeoutput=apodizeoutput, gammaApo=gammaApo, suppressR=suppressR, wiener=wiener, nosuppress=nosuppress, driftfix = driftfix, noOrder0=noOrder0, zpadto=zpadto, makemodel=makemodel, keeporder2=keeporder2))   # default Full frame Recon. parameters    
                sim_kwargs = dict(                                                                                                            
                input_file= input_name,
                otf_file= otf_name,
                #ls = 0.315)
                ls= (wl/1000)/2/0.775)
                #ls = 0.2035)
                    
                sim_kwargs.update(base_kwargs)
                            
                #create processed file output name
                output_file = sim_kwargs["input_file"].replace(".mrc", '_proc'+'_' + user_text + ".mrc")
                sim_kwargs["output_file"] = output_file

                otf_mrc = Mrc.Mrc(otf_name)
                otf_data = otf_mrc.data
                otf_fig, otf_log_fig, psf_fig = get_otf_figures(otf_data)

                input_path = os.path.dirname(input_name)
                figure_path = os.path.join(input_path, user_text)
                
                if not os.path.exists(figure_path):
                    os.makedirs(figure_path)

                if tiled: 
                    output_file, sim_output = split_process_recombine(sim_kwargs["input_file"], tile_size, tile_overlap, sim_kwargs, tile_limits = tile_limits)

                else:
                    sim_output = simrecon(**sim_kwargs)
                print(sim_output)
                for i in range(len(sim_output)):
                    sim_output[i] = sim_output[i].replace('\r', '')
                
                if tiled:
                    with open(output_file.replace(".mrc", ".txt"), "w") as myfile:
                        myfile.write(str(sim_kwargs))
                        temp = "\n".join(sim_output)
                        myfile.write(temp)

                        txt_output = process_txt_output(temp)
                        fig, axs = plot_params(txt_output)
                        tile_params_name = os.path.join(figure_path, user_text + "_tile_params.png")

                        # Note: The figure will show the un-filtered tile parameter values
                        # In the future this should be changed but currently the code to parse the text file
                        # only works if the parameter search is performed.  This is a limitation of the current code.
                        fig.savefig(tile_params_name, dpi=300, bbox_inches="tight")


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
                #fig.savefig(os.path.join(input_path, 'wiener_gA1.0_triangle', user_text + "_xyz_center_crop3.png"))

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