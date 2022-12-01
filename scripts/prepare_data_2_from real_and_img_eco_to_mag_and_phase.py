import nibabel as nib
import os
import numpy as np

def prep_data(data_path, basename_image):
# Files names for the multiple echoes and outputs
    name_images = [basename_image + '_e1', basename_image + '_e2', basename_image + '_e3',basename_image + '_e4',basename_image + '_e5',basename_image + '_e6',basename_image + '_e7',basename_image + '_e8',basename_image + '_e9']
    suffix_parts = {'real': '_real', 'imag': '_imaginary'}
    suffix_image = '.nii'
    fname_out_real = os.path.join(data_path, basename_image + '_real' + suffix_image)
    fname_out_imag = os.path.join(data_path, basename_image + '_imag' + suffix_image)
    fname_out_magnitude = os.path.join(data_path, basename_image + '_magnitude' + suffix_image)
    fname_out_phase = os.path.join(data_path, basename_image + '_phase' + suffix_image)
    hdr = None

    magnitude_echoes, phase_echoes = [], []
    real_echoes, imag_echoes = [], []
    for i, fname_echo in enumerate(name_images):
        fname_image_real = fname_echo + suffix_parts['real'] + suffix_image
        fname_image_imag = fname_echo + suffix_parts['imag'] + suffix_image
        path_image_real = os.path.join(data_path, fname_image_real)
        path_image_imag = os.path.join(data_path, fname_image_imag)

        # Loading the images and modifying the header
        image_real = nib.load(path_image_real)
        real = image_real.get_fdata()
        image_imag = nib.load(path_image_imag)
        imag = image_imag.get_fdata()
        if i == 0:
            hdr = image_real.header.copy()
            hdr['dim'][4] = 3

        # Calculating the magnitude and phase based on complex data
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

# path="/Users/noeedc/Documents/Stage2021/data_new"
# files = [f for f in os.listdir(path) if not f.startswith('.')]
# for f in files:
#     data_path = os.path.join(path,f)
#     basename_image = os.listdir(f)[0].split('_e')[0]    
#     prep_data(data_path, basename_image)


prep_data('/home/magic-chusj-2/Desktop/testData/MEDI_data/03_Invivo_GE/dcm2niix', 'AXL_QSM_MRB2_Stroke-Tramua_6.08_20111030211857_5')