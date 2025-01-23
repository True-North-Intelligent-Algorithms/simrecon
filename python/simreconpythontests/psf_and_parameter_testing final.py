import os
import re
import sys
import mrc as Mrc
from tnia.plotting.projections import show_xyz_slice_center, show_xyz_slice
sys.path.insert(1, r'./python')
from simreconpython.fft_helper import show_centered_fft_magnitude, centered_fft_magnitude
from simreconpython.otf_helper import get_otf_figures
import numpy as np
from tnia.nd.ndutil import centercrop
import matplotlib.pyplot as plt
import json

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

def psf_and_parameter_testing(json_name):
    print(os.getcwd())
    print("PSF and parameter testing")

    # Read the JSON file
    with open(json_name, 'r') as file:
        data = json.load(file)
    
    gammaApo = data['gammaApo'] 
    suppressR = data['suppressR']
    wiener = data['wiener']
    driftfix = data['driftfix']
    nosuppress = data['nosuppress']
    zpadto = data['zpadto']
    noOrder0 = data['noOrder0']
    makemodel = data['makemodel']

    tiled = data['tiled']
    tile_size = data['tile_size']
    tile_overlap = data['tile_overlap']
    filter_tiles = data['filter_tiles']

    keeporder2 = data['keeporder2']

    apodizeoutput = data['apodizeoutput']

    nofilteroverlaps = data['nofilteroverlaps']
    forcemodamp = data['forcemodamp']
    o1 = data['o1']
    o2 = data['o2']
    nok0search = data['nok0search']
    
    # Loop through each dataset and process it
    for dataset in data["datasets"]:

        input_name = dataset["input_name"]
        otf_name = dataset["otf_name"]
        wl = dataset["wl"]
        na = dataset["na"]

        print(f"Processing dataset:")
        print(f"Input Name: {input_name}")
        print(f"OTF Name: {otf_name}")
        print(f"Wavelength (wl): {wl}")
        print(f"Numerical Aperture (na): {na}")
        # Add your processing logic here

        z_display_slice = 9 
        
        prefix = ''

        # set default parameters
        base_kwargs = dict(
                        nphases=5,
                        ndirs=3,
                        angle0= 1.29,#1.36,
                        negDangle=True,              
                        na= na,
                        nimm= 1.0,
                        zoomfact= 2.0, 
                        background= 100.0,           
                        fastSIM=True,
                        otfRA= True,
                        dampenOrder0=True,
                        k0searchall=True,
                        equalizez=True,
                        preciseapo=True,
                        nthreads = 8
                    )
        '''
        gammaApo = 0.3
        suppressR = 10.0
        wiener = 0.001
        driftfix = False 
        nosuppress = False
        zpadto = 0 
        noOrder0 = False
        makemodel = False
        
        tiled = False
        tile_size = 64 
        tile_overlap = 32
        filter_tiles = False
        
        keeporder2 = False
        
        apodizeoutput = 1
        
        nofilteroverlaps = False
        forcemodamp = False
        o1 = 1.0
        o2 = 0.1
        nok0search = False
        '''
                   
        user_text = prefix+'gApo_'+str(gammaApo)+'_supR_'+str(suppressR)+'_w_'+str(wiener)+'_wl_'+str(wl)

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
            user_text=user_text+'_filter_tile3'

            tile_limits = {}
            tile_limits['spacing_min'] = 0.41
            tile_limits['spacing_max'] = 0.42
            tile_limits['spacing_default']=0.414
            tile_limits['angle_min'] = 1.3
            tile_limits['angle_max'] = 1.4
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
            base_kwargs.update(dict(nofilteroverlaps=nofilteroverlaps, apodizeoutput=apodizeoutput, gammaApo=gammaApo, suppressR=suppressR, wiener=wiener, nosuppress=nosuppress, driftfix = driftfix, noOrder0=noOrder0, zpadto=zpadto, makemodel=makemodel, forcemodamp=forcemodamp, nok0search=nok0search))   # default Full frame Recon. parameters    
        else:
            base_kwargs.update(dict(nofilteroverlaps=nofilteroverlaps, apodizeoutput=apodizeoutput, gammaApo=gammaApo, suppressR=suppressR, wiener=wiener, nosuppress=nosuppress, driftfix = driftfix, noOrder0=noOrder0, zpadto=zpadto, makemodel=makemodel, keeporder2=keeporder2))   # default Full frame Recon. parameters    
        
        base_kwargs.update(dict(                                                                                                            
        input_file= input_name,
        otf_file= otf_name,
        #ls = 0.315)
        ls= (wl/1000)/2/0.81
        #ls = 0.2035)
        ))
            
        base_kwargs.update(base_kwargs)

        #create processed file output name
        output_file = base_kwargs["input_file"].replace(".mrc", '_proc'+'_' + user_text + ".mrc")
        base_kwargs["output_file"] = output_file

        otf_mrc = Mrc.Mrc(otf_name)
        otf_data = otf_mrc.data
        otf_fig, otf_log_fig, psf_fig = get_otf_figures(otf_data)

        input_path = os.path.dirname(input_name)
        figure_path = os.path.join(input_path, user_text)

        if not os.path.exists(figure_path):
            os.makedirs(figure_path)

        if tiled: 
            output_file, sim_output = split_process_recombine(base_kwargs["input_file"], tile_size, tile_overlap, base_kwargs, tile_limits = tile_limits)

        else:
            sim_output = simrecon(**base_kwargs)

        print(sim_output)
        for i in range(len(sim_output)):
            sim_output[i] = sim_output[i].replace('\r', '')

        with open(output_file.replace(".mrc", ".txt"), "w") as myfile:
            myfile.write(str(base_kwargs))
            temp = "\n".join(sim_output)
            myfile.write(temp)

        if tiled:
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

        #yc = 697
        #xc = 843
        yc1 = reconstructed.shape[1]//2
        xc1 = reconstructed.shape[2]//2
        reconstructed_crop = reconstructed[:, max(yc1-512,0):min(yc1+512, reconstructed.shape[1]-1), max(xc1-512,0):min(xc1+512, reconstructed.shape[2]-1)]
        fig = show_xyz_slice_center(reconstructed_crop, sxy=1, sz=5)
        fig.savefig(os.path.join(figure_path, user_text + "_xyz_center_crop.png"))

        #yc2 = 527
        #xc2 = 773
        yc2 = reconstructed.shape[1]//2
        xc2 = reconstructed.shape[2]//2
        size = 256
        reconstructed_crop = reconstructed[:, max(yc2-size,0):min(yc2+size, reconstructed.shape[1]-1), max(xc2-size,0):min(xc2+size, reconstructed.shape[2]-1)]
        fig = show_xyz_slice_center(reconstructed_crop, sxy=1, sz=5)
        fig.savefig(os.path.join(figure_path, user_text + "_xyz_center_crop2.png"))


        #yc = 527
        #xc = 818
        yc3 = reconstructed.shape[1]//2
        xc3 = reconstructed.shape[2]//2
        size = 128
        reconstructed_crop = reconstructed[:, max(yc3-size,0):min(yc3+size, reconstructed.shape[1]-1), max(xc3-size,0):min(xc3+size, reconstructed.shape[2]-1)]
        fig = show_xyz_slice_center(reconstructed_crop, sxy=1, sz=5)
        fig.savefig(os.path.join(figure_path, user_text + "_xyz_center_crop3.png"))
        #fig.savefig(os.path.join(input_path, 'wiener_gA1.0_triangle', user_text + "_xyz_center_crop3.png"))

        center_line = reconstructed[:, yc1, xc1]
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
    json_name = r'D:\Janelia\gold_standard_set_final\repeat_old.json'
    psf_and_parameter_testing(json_name)