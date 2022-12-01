#!/usr/bin/env python
# coding: utf-8

####################################
# Filename: batchRegister.py
# Author: Noée Ducros-Chabot 
# co-Author: Daniel Ridani
# Date: 06/06/2021
# Python version: 3.8.3
# 
# Registers magnitude image to template space
# 
# Inputs
# ----------
# static_path : str
#     Path to static nifti image (used as template for registration)
# moving_path : str
#     Path to nifti image to be registered. 
# output_path : str
#     Path to output folder where the registered image and the DiffeomorphicMap will be saved.
# qsm_path : str
#     Path to suceptibility map nifti image (QSM result) to be transformed. 
#     Leave blank ([] or False) for no qsm tranformation.
####################################

import numpy as np
import nibabel as nib
from dipy.viz import regtools
from dipy.align.imaffine import (transform_centers_of_mass,
                                 AffineMap,
                                 MutualInformationMetric,
                                 AffineRegistration)
from dipy.align.transforms import (TranslationTransform3D,
                                   RigidTransform3D,
                                   AffineTransform3D)
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.metrics import *
import matplotlib as plt
import pickle
import os

def register(static_nii, moving_nii):
    """
    Symmetric Diffeomorphic Registration of moving image to template image world.

    Parameters
    ----------
    static_nii : NiBabel image
        Template image.
    moving_nii : str
        Moving image. This image will be registered to template coordinate.

    Returns
    -------
    registered_img np.array
        Registered moving image 
    
    mapping dipy.align.imwarp.DiffeomorphicMap
        Diffeomorphic map

    """
    static = static_nii.get_fdata()
    static_grid2world = static_nii.affine

    moving = moving_nii.get_fdata()
    moving_grid2world = moving_nii.affine

    print('Rigid Registration')
    identity = np.eye(4)
    affine_map = AffineMap(identity, static.shape, static_grid2world, moving.shape, moving_grid2world)
    resampled = affine_map.transform(moving)

    c_of_mass = transform_centers_of_mass(static, static_grid2world, moving, moving_grid2world)

    nbins = 32
    sampling_prop = None
    metric = MutualInformationMetric(nbins, sampling_prop)

    level_iters = [10000, 1000, 100]
    sigmas = [3.0, 1.0, 0.0]
    factors = [4, 2, 1]
    affreg = AffineRegistration(metric=metric,level_iters=level_iters,sigmas=sigmas,factors=factors)

    transform = TranslationTransform3D()
    params0 = None
    starting_affine = c_of_mass.affine
    translation = affreg.optimize(static, moving, transform, params0,static_grid2world, moving_grid2world, starting_affine=starting_affine)

    transform = RigidTransform3D()
    params0 = None
    starting_affine = translation.affine
    rigid = affreg.optimize(static, moving, transform, params0,
                            static_grid2world, moving_grid2world,
                            starting_affine=starting_affine)

    print('Affine Registration')
    transform = AffineTransform3D()
    params0 = None
    starting_affine = rigid.affine
    affine = affreg.optimize(static, moving, transform, params0,
                            static_grid2world, moving_grid2world,
                            starting_affine=starting_affine)

    transformed = affine.transform(moving)

    print('Symmetric Diffeomorphic Registration')
    metric = CCMetric(3, sigma_diff = 5)
    level_iters = [10,10,5] 
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

    mapping = sdr.optimize(static, moving, static_grid2world, moving_grid2world, affine.affine)

    return mapping.transform(moving), mapping

if __name__ == "__main__":
    static_path = "/home/magic-chusj-2/Desktop/Registration/dhcp-volumetric-atlas-groupwise/mean/ga_40/average_t1_masked.nii"
    moving_path = "/home/magic-chusj-2/Desktop/QSM/QSM.m/Results/HIE_baby_001/magMask.nii"
    output_path = '/home/magic-chusj-2/Desktop/Registration/results/dipy/HIE_baby_001'
    qsm_path = "/home/magic-chusj-2/Desktop/QSM/QSM.m/Results/HIE_baby_001/resharp_rts/chi.nii"

    static_nii = nib.load(static_path)
    moving_nii = nib.load(moving_path)

    # Registration
    warped_moving, mapping = register(static_nii, moving_nii)

    # Saving registered magnitude image
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    img = nib.Nifti1Image(warped_moving, static_nii.affine, header =static_nii.header )
    nib.save(img, os.path.join(output_path,"registered_mag"))

    # Saving DiffeomorphicMap
    with open(os.path.join(output_path,'transformation.pkl'), 'wb') as f:
        pickle.dump(mapping, f)

    # Apply transformation to QSM image if exist:
    if qsm_path and os.path.isfile(qsm_path):
        qsm_nii = nib.load(qsm_path)
        qsm = qsm_nii.get_fdata()

        warped_qsm = mapping.transform(qsm)

        img = nib.Nifti1Image(warped_qsm, static_nii.affine, header =static_nii.header )
        nib.save(img, os.path.join(output_path,"registered_qsm"))