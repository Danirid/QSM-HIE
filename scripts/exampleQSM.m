%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
% Inputs :
% dicomDir - path to folder containing dicom files
% niftifile - path to nifti file to copy header info (for output files)
% bckremoval - methode to be used for background removal. options: irsharp,
%              ismv, lbv, pdf, resharp (default), sharp, vsharp
% dipoleInv - methode used for dipole inversion. options : ilsqr, medi,
%             ndi, rts (default), sstgv, tikh, tkd, tsvd, tv.
% 
% B0 - Field strength (to be changed manually) 
% f - Imaging frequency (to be changed manually)
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
run('addpathqsm.m');
addpath('scripts');
%%
% INPUTs
dicomDir = 'path/to/dicom';
niftifile ='path/to/nifti';
outpath = 'path/to/output';
fMask = 'path/to/mask'; %optional

% QSM algorithms to use
bckremoval = 'resharp'; 
dipoleInv = 'rts';

if ~exist(outpath, 'dir')
       mkdir(outpath)
end

% Constants for normalizing phase
B0 = 3;
GYRO = 2*pi*127.7993/B0; % 2*pi*f=y*B0.
%%
disp('Loading the data...');
% Load dataset
[hdr, data] = dicom(dicomDir);

mag = data(:,:,:,:,1);
phas = data(:,:,:,:,2);
phas = -phas;

% if magnitude and phase stored in separate files
%  [hdr, data] = dicom(filenameMag);
%  mag = data(:,:,:,:,1);
%  [hdr, data] = dicom(filenamePhase);
%  phas = data(:,:,:,:,2);

vsz = hdr.vsz;
TEs = hdr.TEs;
bdir = hdr.bdir;

%%
% brain extraction using FSL's bet. won't work on windows.
if isfile(fMask)
    disp('Loading modified brain mask...');
    mask_ = niftiread(fMask);
    mask0 = double(flip(mask_,2));
else
    disp('Creating brain mask...');
    mask0 = generateMask(mag(:,:,:,end), vsz, '-m -n -f 0.5'); %fsl necessary
    % mask0 = BET(mean(mag,4),size(mag(:,:,:,end)),vsz); % Add function
    % from MEDI toolbox in this toolbox (to extract brain mask without fsl).
end

%%
% erode mask to deal with boundary inconsistencies during brain extraction
mask1 = erodeMask(mask0, 5);
%%
% unwrap phase + background field removing
disp('Laplacien unwraping and background field removing...');
uphas = unwrapLaplacian(phas, mask1, vsz);

%convert units
for t = 1:size(uphas, 4)
    uphas(:,:,:,t) = uphas(:,:,:,t) ./ (B0 * GYRO * TEs(t));
end

uphas = mean(uphas,4); %combining TE images
%%
% remove non-harmonic background fields
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
%%
% dipole inversion
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
%%
% save mat-file
% save qsm.mat mag phas uphas x mask1 bdir vsz TEs
%%
% or nifti
disp('Saving results...');
writeNifti(niftifile, fullfile(outpath, 'chi.nii'), flip(x,2)); 
writeNifti(niftifile, fullfile(outpath, 'mag.nii'), flip(mean(mag,4), 2));
writeNifti(niftifile, fullfile(outpath, 'erodedMask.nii'), flip(mask1, 2));
writeNifti(niftifile, fullfile(outpath, 'magMask.nii'), flip(mean(mag,4).*mask0, 2));
writeNifti(niftifile, fullfile(outpath, 'mask.nii'), flip(mask0, 2));
