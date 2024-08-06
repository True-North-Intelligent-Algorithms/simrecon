import os
import sys
import mrc as Mrc
from tnia.plotting.projections import show_xyz_slice_center, show_xyz_slice
sys.path.insert(1, r'./python')
from simreconpython.fft_helper import show_centered_fft_magnitude, centered_fft_magnitude
import numpy as np

def explore_reconstructed():
    print(os.getcwd())
    print("Exploring reconstructed data")

    input_path = r"D:\Janelia\Data 2024-06-06\Wiener, gammaApo and SupressR parameter testing\488nm comparison Brian\CELL 4 - 1.0W 300ms updated HWPQWP pos_20240503_115909"
    figure_path = os.path.join(input_path, "figures2")

    if not os.path.exists(figure_path):
        os.makedirs(figure_path)

    reconstructed_names = [r'488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.3_supR_10_w_0.001_wl_488_jsPSF.mrc',
        r'488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.3_supR_10_w_0.001_wl_488_dhPSF.mrc',
        r'488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.3_supR_10_w_0.001_wl_488_bnPSF.mrc']
    
    reconstructed_names =  [r'488 nm 5 phases 0.81 NA React_All Linear SIM_cam1_1_proc_gApo_0.3_supR_10_w_0.001_wl_488_bnPSF4.mrc']
    
    for reconstructed_name in reconstructed_names:
        full_name = os.path.join(input_path, reconstructed_name)
        reconstructed_mrc = Mrc.Mrc(full_name)
        reconstructed = reconstructed_mrc.data
        print(reconstructed.shape)
        fig = show_xyz_slice_center(reconstructed, sxy=1, sz=20)
        print(reconstructed.min(), reconstructed.max())
        fig.savefig(os.path.join(figure_path, reconstructed_name + "_xyz_center.png"))
        
        fft = np.fft.fftn(reconstructed)
        fft = np.fft.fftshift(fft)
        fft_app_mag = np.abs(fft)
        fft_app_mag_log = np.log(fft_app_mag + 1)
        fig = show_xyz_slice_center(fft_app_mag_log, sxy=1, sz=20)
        fig.savefig(os.path.join(figure_path, reconstructed_name + "_fft_center.png"))


if __name__ == "__main__":
    explore_reconstructed()