function [] = writeNifti(niftipath, outputpath, image)
% Reformating hdr based on input image
info = niftiinfo(niftipath);
info.Datatype = class(image);
if isa(image,'single')
    info.BitsPerPixel = 32; 
elseif isa(image, 'double')
    info.BitsPerPixel = 64; 
elseif isa(image, 'int16')
    info.BitsPerPixel = 16; 
elseif isa(image, 'int32')
    info.BitsPerPixel = 32;
else
    warning('BitsPerPixel of image not modified in nifti header')  
end
    
info.ImageSize = size(image);
info.Description = '';
if length(info.PixelDimensions)>3 && length(size(image))<=3
    info.PixelDimensions = info.PixelDimensions(1:end-1);
end

info.AdditiveOffset = 0;
info.MultiplicativeScaling = 1;

% Writing nifti image
niftiwrite(image, outputpath, info)

