% This code batch applies QSM reconstruction on a population (see README.md
% for proper format). 
% Author: Noée Ducros-Chabot 
% co-Author: Daniel Ridani
%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
% Inputs : 
% directory - path to directory containing datasets
% bckremoval - methode to be used for background removal. options: irsharp,
%              ismv, lbv, pdf, resharp (default), sharp, vsharp
% dipoleInv - methode used for dipole inversion. options : ilsqr, medi,
%             ndi, rts (default), sstgv, tikh, tkd, tsvd, tv.
% B0 - Field strength (default: 3) 
% f - Imaging frequency (default: 127.7993)
% niftifile - path to nifti file to copy header info (for output files)
%
% mag - 3D or 4D magnitude image
% phas - wrapped phase (3d/4d array)
% vsz - voxel size. 
% TEs - echos times in seconds.
% bdir - unit vector of B field direction.
% mask0 - mask of region of interest (3d array).
% mask1 - eroded mask of region of interest (3d array).
% uphas - unwrapped local phase (3d array).
% fl - local field (3d aray).
% x - susceptibility map.

%%
% add QSM.m to path
run('/path/to/QSM.m/addpathqsm.m'); % Add local path to QSM.m's addpathqsm.m file.
addpath('scripts');
%% Inputs 
% Path to data (refer to README.md for formating)
directory = '/path/to/dataset/'; 

% QSM algorithms to use
bckremoval = 'resharp'; 
dipoleInv = 'rts';
%% 
d = dir(directory);
dfolders = d([d(:).isdir]);
patients = dfolders(~ismember({dfolders(:).name},{'.','..'}));

% Constants for normalizing phase
B0 = 3;
f = 127.7993; % See ImagingFrequency/(0018,0084) dicom attribute. 
GYRO = 2*pi*f/B0;
%%
for patient = patients'
    % DICOM
    folder = patient.name;
    disp(newline);
    disp(folder);
    filename = fullfile(directory,folder, 'DICOM');
    
    niftifile = dir(fullfile(directory,folder, '*magnitude.nii.gz'));
    niftifile = fullfile(directory,folder,niftifile(1).name);

    outpath = fullfile(directory,folder,'QSM_results', append(bckremoval, '_',dipoleInv)); 
    
    if ~exist(outpath, 'dir')
           mkdir(outpath)
    end
    %%
    disp('Loading the data...');
    % Load dataset
    [hdr, data] = dicom(filename);
    mag = data(:,:,:,:,1);
    phas = data(:,:,:,:,2);
    phas = -phas;

    vsz = hdr.vsz;
    TEs = hdr.TEs;
    bdir = hdr.bdir;
    %% brain extraction using FSL's bet. won't work on windows.
    fMask = fullfile(outpath,'..', 'mask_.nii');
    if isfile(fMask)
        disp('Loading modified brain mask...');
        mask_ = niftiread(fMask);
        mask0 = logical(flip(mask_, 2));
    else
        disp('Creating brain mask...');
        % brain extraction using FSL's bet. won't work on windows.
        mask0 = generateMask(mag(:,:,:,end), vsz, '-m -n -f 0.5');
    end

    % erode mask to deal with boundary inconsistencies during brain extraction
    mask1 = erodeMask(mask0, 5);
    
    %% unwrap phase + background field removing
    disp('Laplacien unwraping and background field removing...');
    uphas = unwrapLaplacian(phas, mask1, vsz);

    % convert units
    for t = 1:size(uphas, 4)
        uphas(:,:,:,t) = uphas(:,:,:,t) ./ (B0 * GYRO * TEs(t));
    end
    uphas = average(uphas,4); %Combining TE images
    
    %% remove non-harmonic background fields
    disp('Removing non-harmonic background fields...');
    switch bckremoval
        case 'sharp'
            [fl, mask1] = sharp(uphas, mask1, vsz);
        case 'vsharp'
            [fl, mask1] = vsharp(uphas, mask1, vsz);
        case 'resharp'
            [fl, mask1] = resharp(uphas, mask1, vsz);
        case 'irsharp' 
            [fl, mask1] = irsharp(uphas, phas, mask1, vsz);
        case 'lbv'
            [fl] = lbv(uphas, mask1, vsz);
        case 'pdf'
            [fl] = pdf(uphas, mask1, vsz, [], bdir);
        case 'ismv'
            [fl, mask1] = ismv(uphas, mask1, vsz);
        otherwise
            warning('Did not recognize background removal methode. Using resharp.');
            [fl, mask1] = resharp(uphas, mask1, vsz);
    end
    %% dipole inversion
    disp('Dipole inversion...');
    switch dipoleInv
        case 'ilsqr'
            [x_ilsqr, xsa, xfs, xlsqr] = ilsqr(fl, logical(mask1), vsz, bdir);
        case 'medi'
            x = medi(fl, mask1, vsz, mag, [], bdir);
        case 'ndi'
            x = ndi(fl, mask1, vsz, [], bdir);
        case 'rts' 
            x = rts(fl, mask1, vsz, bdir);
        case 'tikh'
            x = tikh(fl, mask1, vsz, bdir);
        case 'tkd'
            x  = tkd(fl, mask1, vsz, bdir);
        case 'tsvd'
            x = tsvd(fl, mask1, vsz, bdir);
        case 'tv'
            x = tv(fl, mask1, vsz, bdir);
        otherwise
            warning('Did not recognize dipole inversion methode. Using rts.');
            x = rts(fl, mask1, vsz, bdir);
    end
    
    %% save resultss
    disp('Saving results...');
    
    if ~isfile(fullfile(outpath,'..', 'mag.nii'))
        writeNifti(niftifile, fullfile(outpath, '..','mag.nii'), flip(mean(mag,4),2));
    end
    if ~isfile(fullfile(outpath,'..', 'magMask.nii'))
        writeNifti(niftifile, fullfile(outpath,'..', 'magMask.nii'), flip(mean(mag,4).*mask0,2));
    end
    if ~isfile(fullfile(outpath,'..', 'mask.nii'))
        writeNifti(niftifile, fullfile(outpath,'..', 'mask.nii'), flip(int16(mask0),2));
    end
    
    writeNifti(niftifile, fullfile(outpath,'erodedMask.nii'), flip(int16(mask1),2));
    writeNifti(niftifile, fullfile(outpath, 'chi.nii'), flip(x,2));

end
