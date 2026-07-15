function mask = createMask(im, threshold)
%createMask Create a binary mask from an image using a threshold.
%
%   Inputs:
%   im: Input image (2D or 3D)
%   threshold: Threshold value for creating the mask
%
%   Outputs:
%   mask: Binary mask (same size as input image)
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    im {mustBeNumeric, mustBeReal}
    threshold (1,1) {mustBeNumeric} = graythresh(im)
end

% Apply threshold
mask = im > threshold;

% Clean up mask
for iSlice = 1:size(mask, 3)
    mask(:,:,iSlice) = imfill(mask(:,:,iSlice), 'holes');
    % mask(:,:,iSlice) = bwareaopen(mask(:,:,iSlice), 100);
end
if size(mask,3) > 1
    se = strel('sphere', 3);
    mask = imopen(mask, se);
    mask = imclose(mask, se);
end
se = strel('disk', 3);
mask = imopen(mask, se);
mask = imclose(mask, se);


end
