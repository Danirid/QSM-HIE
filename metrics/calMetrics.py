#!/usr/bin/env python
# -*- coding: utf-8 -*-

####################################
# Filename: batchRegister.py
# Author: Noée Ducros-Chabot 
# co-Author: Daniel Ridani
# Date of creation: August 6th, 2021
# Python version: 3.8.3
# 
# Description : This script calcuates CsVO2 and cartography of susceptibility from registered 
# qsm and structures/tissues label map and csv (from dchp atlas) for all datasets
#  
# Inputs
# ----------
# directory : str
#     Path to folder containing datasets (refer to README.md for formating)
# tissuesVol_path/structVol_path : str
#     Path to nii tissues/structure label map (for this project dhcp-volumetric-atlas-groupwise/mean/ga_40/tissues[or structures].nii.gz was used) 
# tissues_csv/structures_csv : str
#     Path to csv files containing labels for label maps regions (see github config files)
# outpath : str
#      Path to where csv file containing the metrics will be created
#
####################################
import os
import numpy as np
import matplotlib as plt
import nibabel as nib
import pandas as pd
import seaborn as sns
import glob

def calCSvO2(qsm, mask, perc = 99.95):
    thresh = np.percentile(qsm, perc)

    # Mean susceptibility of cerebral veins 
    x_vein = np.mean(qsm[qsm>thresh])
    
    # Mean cerebral susceptibility
    mask[qsm>thresh]=0 
    x_avg = np.mean(qsm[mask])

    # Calculating CSvO2
    dx = x_vein - x_avg
    dx_do = 4*np.pi*0.21
    HCT = 0.5
    CSvO2 = 1 - dx/(dx_do*HCT)
    return thresh, CSvO2

def formatDict(dict):
    # format dictionary values from str to list of int
    for key in dict: 
        dict[key] = dict[key].strip('][').split(',')
        dict[key] = [int(i) for i in dict[key]]


def QSM_map(qsm, label_map, label_dict, mask):
    # qsm : qsm image (nifti) to extract mean values (for certain zone)
    # label_map : label map image (nifti)
    # label_dict: dictionary containing label name and their values
    # mask : erroded mask, bool image to apply to label map 

    QSM_values = label_dict.copy()
    label_map[~mask]=0 
    for key in label_dict: 
        values = label_dict[key]
        if isinstance(values, list):
            if len(values)==1:
                x = np.mean(qsm[label_map == values[0]])
            elif len(values)>1: 
                means = []
                for value in values:
                    x = np.mean(qsm[label_map == value])
                    means.append(x)
                x = np.mean(means)
            else: 
                print(f'Check for error with {key} values')
                x = None
        else: 
            # Assumes that if not list than is a single value
            x = np.mean(qsm[label_map == values])
        QSM_values[key] = x
    return QSM_values

if __name__ == "__main__":
    directory = '/path/to/Data'
    tissuesVol_path =  '/path/to/dhcp-volumetric-atlas-groupwise/mean/ga_40/tissues.nii.gz'
    tissues_csv = os.path.join(os.getcwd(), 'config/tissues.csv')
    structVol_path = '/path/to/dhcp-volumetric-atlas-groupwise/mean/ga_40/structures.nii.gz'
    struct_csv = os.path.join(os.getcwd(), 'config/structures.csv')
    outpath = '/path/to/outfolder'

    #Loading label maps and their dictionaries (from csv file)
    tissues_nii = nib.load(tissuesVol_path)
    tissues = tissues_nii.get_fdata()
    tissues_label = pd.read_csv(tissues_csv, header = None, index_col=0, squeeze=True).to_dict()
    struct_nii = nib.load(structVol_path)
    struct = struct_nii.get_fdata()
    struct_label = pd.read_csv(struct_csv, header = None, index_col=0, squeeze=True).to_dict()
    formatDict(struct_label)

    IDs = []
    bckgRemoval = []
    inversion = []
    CSvO2s = []
    threshs = []
    tissues_values = {i:[] for i in tissues_label.keys()}
    struct_values = {i:[] for i in struct_label.keys()}

    folders = sorted(glob.glob(os.path.join(directory, "*/")))
    # Loop through all datasets
    for folder in folders:  
        print(folder.split('/')[-2])
        dir = os.path.join(folder,'Registration_results')
        methodes = [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]
        # Loop through all methodes (mix of algorithms ex.resharp_rts)
        for methode in methodes: 
            ID = folder.split('/')[-2]
            IDs.append(ID)

            # Creating dataframe
            bckgRemoval.append(methode.split('_')[0])
            inversion.append(methode.split('_')[-1])

            # Loading unregistered QSM data and mask
            qsm_nii = nib.load(os.path.join(folder, 'QSM_results', methode, "chi.nii"))
            qsm = qsm_nii.get_fdata()
            mask_nii = nib.load(os.path.join(folder, 'QSM_results', methode, "erodedMask.nii"))
            mask = mask_nii.get_fdata().astype(bool)

            # Calculating CSvO2 value
            thresh, CSvO2 = calCSvO2(qsm,mask)
            CSvO2s.append(CSvO2)
            threshs.append(thresh)

            # Loading registered QSM data
            qsm_nii = nib.load(os.path.join(dir, methode, "registered_qsm.nii"))
            qsm = qsm_nii.get_fdata()

            # Loading registered erroded mask
            erodedMask_nii = nib.load(os.path.join(dir, methode, "registered_mask.nii"))
            erodedMask = erodedMask_nii.get_fdata()
            erodedMask[erodedMask>0]=1 # Bc some pixel took intermediate values when transformed
            erodedMask = erodedMask.astype(bool)

            # Calculating QSM values from tissues label map
            x_tissues = QSM_map(qsm, tissues, tissues_label, erodedMask)
            x_structs = QSM_map(qsm, struct, struct_label, erodedMask)
            for key in x_tissues:
                tissues_values[key].append(x_tissues[key])
            for key in x_structs: 
                struct_values[key].append(x_structs[key])

    # Creating dataframe
    data = {'ID': IDs, 'bckgRemoval' : bckgRemoval, 'dipoleInv' : inversion, 'Thresh': threshs, 'CSvO2': CSvO2s}
    all_data = {**data , **struct_values, **tissues_values}
    df = pd.DataFrame(all_data).set_index('ID')

    if not os.path.exists(outpath):
            os.makedirs(outpath)

    df.to_csv(os.path.join(outpath,'metrics.csv'))