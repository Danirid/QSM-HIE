# QSM-HIE
This project aims to develop and compare an image reconstruction pipeline to assess brain oxygenation in neonatal asphyxia using quantitative susceptibility mapping (QSM)
The basic steps of this project consist of: 
- Applying QSM.m on a datasets. 
- Registering the QSM magnitudes of each subject to template space of **dhcp's 40 ga mean as template**.
- Applying the transforms to the susceptibility maps.
- Extracting metrics from registered QSM results.

## Getting started

### Prerequisites

* [QSM.m toolbox](https://github.com/kamesy/QSM.m.git)
* [dhcp-volumetric-atlas-groupwise](https://gin.g-node.org/BioMedIA/dhcp-volumetric-atlas-groupwise.git) (for registration template and cartography)
* Python packages (in config/requirements.txt)
 
### Repos architecture
```bash
QSM_HIE
├── README.md
├── batchQSM.m 		# Batch QSM reconstruction on data
├── batchRegister.py		# Performs registration and apply on QSM results on all subjects.
├── registration
├── metrics
├── scripts
└── config 			# python dependencies and label map legend used for cartography. 
```

To run the batch scripts, the data should follow the form: 
```
Data
├── HIE_baby_0XX
│   ├── DICOM			# The folder containing the dicom files must be named DICOM.
│   ├── *magnitude.nii.gz		# Nifti header used to save QSM reconstruction output.
│   │
│   ├── QSM_results 			# Will be created after running batchQSM.m
│   │    ├── magn.nii			
│   │    ├── maskMag.nii
│   │    ├── erodedMask.nii
│   │    ├── mask.nii			# If masks needs manual modification, rename mask_.nii (will be use in further batch_dicom runs)
│   │    ├── method_1 			#ex. resharp_rts
│   │    │     ├── chi.nii
│   │    │     └── erodedMask.nii
│   │    ├── method_2 			#ex. rehsarp_SSTGV
│   │    │     ├── chi.nii
│   │    │     └── erodedMask.nii
│   │    └── ... 
│   │
│   └── Registration_results 		# Will be created after running batchRegister.py
│        ├── registered_mag.nii
│        ├── transformation.pkl
│        ├── method_1
│        │     ├── registered_qsm.nii
│        │     └── registered_mask.nii
│        ├── method_2
│        │     ├── registered_qsm.nii
│        │     └── registered_mask.nii
│        └── ... 
└── ...

```
Note : Any nifti files in right coordinate systeme (with proper header) can be used instead of the *magntiude.nii.gz (need to modify batchQSM.m in consequence). [dcm2niix](https://github.com/rordenlab/dcm2niix) could be usefull to convert DICOM into nifti files.

### Usage
:warning: Warning :warning: In order for QSM.m toolbox to work on GE dicoms (like HIE baby dataset), phase needs to be inverse. This line must be added to the code: 
```
phas = data(:,:,:,:,2);
phas = -phas;
```
#### Batch use:
Run the following scripts in order to apply QSM reconstruction and registration on a dataset:
- batchQSM.m
- batch_register.py

Running both scripts will result in susceptibility maps registered in a template's coordinate system for every subject of the dataset. This allows for extraction of clinical metric and suceptibility cartography.

- Once the registration was performed, registration/batchTransformQSM.py can be used to batch apply transform to QSM results using pickle file (containing DiffeomorphicMap) without having to go through the whole registration process once more.

#### Single use: 
Scripts can also be run individually. 
- To apply QSM on a single subject refer to s cripts/exampleQSM.m (or QSM.m/examples/example2). 
- To register a single subject either use registration/register.py or registration/register.ipynb. The jupyter notebook is recommended for single use, as it displays intermediate registration steps as well as a basic overlay of the final diffeomorphic registration. 

#### Metric assessment: 
- metrics/calMetrics.py: Extracts CSvO2 values and cartography (mean susceptibility values for certain regions of the brain) from registered susceptibility maps of population. The results are saved in csv files. 
- metrics/displayMetrics.ipynb: Plots the metrics from the registered QSM. It can also display a comparison between different methods. This script serves more of an example as what can be extracted and displayed from the csv metrics file.

#### Scripts:
The scripts folder contains an series of scripts that could be useful to the process; such as biased field correction, nifti2gif, contrast enhancement, ect.

#### Files organized by Noée Ducros-Chabot
