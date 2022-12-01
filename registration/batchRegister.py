#!/usr/bin/env python
# coding: utf-8

####################################
# Filename: batchRegister.py
# Author: Noée Ducros-Chabot
# co-Author: Daniel Ridani 
# Date: 06/06/2021
# Python version: 3.8.3
# 
# Description: Script to batch register QSM data and apply transform on QSM results
# 
# Inputs
# ----------
# data_path : str
#     Path to folder containing datasets (refer to README.md for formating)
# template_path : str
#     Path to nifti image used as template for registration. For HIE datasets, dchp mean 40 ga masked was used.
#     See https://gin.g-node.org/BioMedIA/dhcp-volumetric-atlas-groupwise.git
#
####################################

import os
import glob2
import nibabel as nib
import pickle
from pytictoc import TicToc

from registration.register import register

data_path = '/path/to/data'
template_path = "/path/to/dhcp-volumetric-atlas-groupwise/mean/ga_40/average_t1_masked.nii"

# Register datasets: 
static_nii = nib.load(template_path)

datasets = sorted([os.path.abspath(os.path.join(data_path, p)) for p in os.listdir(data_path)])
for dataset in datasets:
    print('\n')
    print(dataset.split('/')[-1])
    # Register magnitude image
    moving_nii = nib.load(os.path.join(dataset, 'QSM_results/magMask.nii'))
    t = TicToc()
    print('Registration:')
    t.tic()
    warped_moving, mapping = register(static_nii, moving_nii)
    t.toc()

    # Saving registered magnitude
    outpath = os.path.join(dataset, 'Registration_results')
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    img = nib.Nifti1Image(warped_moving, static_nii.affine, header =static_nii.header )
    nib.save(img, os.path.join(outpath,"registered_mag"))

    # Saving DiffeomorphicMap in file
    with open(os.path.join(outpath,'transformation.pkl'), 'wb') as f:
        pickle.dump(mapping, f)

    # Applying transform to QSM results
    print('Transforming QSM results')
    qsm_files = glob2.glob(os.path.join(dataset, 'QSM_results','**/chi.nii'))
    for qsm_file in qsm_files: 
        qsm_nii = nib.load(qsm_file)
        qsm = qsm_nii.get_fdata()
        warped_qsm = mapping.transform(qsm)

        # Save registered qsm images
        folder_origin = os.path.dirname(qsm_file)
        folder = folder_origin.replace('QSM_results', 'Registration_results')
        if not os.path.exists(folder):
            os.makedirs(folder)

        img = nib.Nifti1Image(warped_qsm, static_nii.affine, header = static_nii.header)
        nib.save(img, os.path.join(folder, "registered_qsm.nii"))

        # Applying transform to erodedMask and saving results (used for CSvO2 extraction)
        # Mask may vary depending on algo used.
        erodedMask_nii = nib.load(os.path.join(folder_origin, 'erodedMask.nii'))
        erodedMask = erodedMask_nii.get_fdata()

        warped_mask = mapping.transform(erodedMask)
        img = nib.Nifti1Image(warped_mask, static_nii.affine, header = static_nii.header)
        nib.save(img, os.path.join(folder, "registered_mask.nii"))
