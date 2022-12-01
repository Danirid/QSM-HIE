import nibabel as nib
import os
import numpy as np

def prepare_data(data_path, basename_image):
    # Files names for the multiple echoes and outputs
    names_images = []
    for i in range(1,9):
        names_images.append(basename_image + '_e'+str(i))
    suffix_image = '.nii'
    fname_out_magnitude = os.path.join(data_path, basename_image + '_magnitude' + suffix_image)
    fname_out_phase = os.path.join(data_path, basename_image + '_phase' + suffix_image)
    fname_out_real = os.path.join(data_path, basename_image + '_real' + suffix_image)
    fname_out_imag = os.path.join(data_path, basename_image + '_imag' + suffix_image)

    magnitude_echoes, phase_echoes = [], []
    real_echoes, imag_echoes = [], []
    for i, fname_echo in enumerate(names_images):
        fname_image = fname_echo + suffix_image
        path_image = os.path.join(data_path, fname_image)

        # Loading the images and modifying the header
        image = nib.load(path_image)
        data = image.get_fdata()
        if i == 0:
            hdr = image.header.copy()
            hdr['dim'][4] = 3

        # Calculating the magnitude and phase based on complex data
        real, imag = data[:, :, :, 0], data[:, :, :, 1]
        magnitude_echo = np.sqrt(np.power(real, 2) + np.power(imag, 2))
        phase_echo = np.arctan2(imag, real)

        # Adding to the buffer lists
        magnitude_echoes.append(magnitude_echo)
        phase_echoes.append(phase_echo)

        real_echoes.append(real)
        imag_echoes.append(imag)

    # Concatenating echoes into a single matrix
    magnitude = np.stack(magnitude_echoes, axis=-1)
    phase = np.stack(phase_echoes, axis=-1)
    real = np.stack(real_echoes, axis=-1)
    imag = np.stack(imag_echoes, axis=-1)

    # Saving the magnitude and phase images
    image_magnitude = nib.nifti1.Nifti1Image(magnitude.copy(), None, hdr)
    nib.save(image_magnitude, fname_out_magnitude)
    image_phase = nib.nifti1.Nifti1Image(phase.copy(), None, hdr)
    nib.save(image_phase, fname_out_phase)
    image_real = nib.nifti1.Nifti1Image(real.copy(), None, hdr)
    nib.save(image_real, fname_out_real)
    image_imag = nib.nifti1.Nifti1Image(imag.copy(), None, hdr)
    nib.save(image_imag, fname_out_imag)

path="/Users/noeedc/Documents/Stage2021/QSM-MEDI_toolbox/MEDI_data/03_Invivo_GE/AXL_QSM/nifti"
basename_image = os.listdir(path)[0].split('_e')[0]    
prepare_data(path, basename_image)