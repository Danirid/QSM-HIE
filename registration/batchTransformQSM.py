#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply transform (result from registration) to QSM results and eroded mask for all datasets from pickle files (containing DiffeomorphicMap). 
Inputs : directory & static_nii
"""
import os, sys 
import pickle
import nibabel as nib
import glob

directory = '/path/to/Data'
static_nii = nib.load("/path/to/dhcp-volumetric-atlas-groupwise/mean/ga_40/average_t1_masked.nii")

folders = sorted(glob.glob(os.path.join(directory, "*/")))
for folder in folders: 
    print(folder.split('/')[-2])
    trans = pickle.load(open(os.path.join(folder, 'Registration_results', 'transformation.pkl'),'rb'))
    
    # Loop through all methodes (mix of algorithms ex.resharp_rts)
    dir = os.path.join(folder,'QSM_results')
    methodes = [name for name in os.listdir(os.path.join(folder,'QSM_results')) if os.path.isdir(os.path.join(folder, 'QSM_results', name))]
    for methode in methodes: 
        qsm_nii = nib.load(os.path.join(folder, 'QSM_results', methode, "chi.nii"))
        qsm = qsm_nii.get_fdata()
        warped_qsm = trans.transform(qsm)
        img = nib.Nifti1Image(warped_qsm, static_nii.affine, header =static_nii.header )

        out_folder = os.path.join(folder,'Registration_results', methode)
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        nib.save(img, os.path.join(out_folder,"registered_qsm.nii"))

        mask_nii = nib.load(os.path.join(folder, 'QSM_results', methode, "erodedMask.nii"))
        mask = mask_nii.get_fdata()
        warped_mask = trans.transform(mask)
        warped_mask[warped_mask>0]=1 # Bc some pixel took intermediate values when transformed
        warped_mask = warped_mask.astype(bool)
        img = nib.Nifti1Image(warped_mask, static_nii.affine, header =static_nii.header )
        nib.save(img, os.path.join(out_folder, "registered_mask.nii"))        
        