function imWarped = warpImage(im, imSpacing, displ, displSpacing, interpMethod, extrapVal)
%warpImage Warp an image using a displacement field.
%
%   Inputs:
%   im: Input image as a 4D array [x,y,z,t]
%   imSpacing: Voxel spacing in each spatial dimension of the input image
%   displ: Displacement field as a 5D array [x,y,z,u,t]
%   displSpacing: Voxel spacing in each spatial dimension of the
%       displacement field
%   interpMethod (optional): Interpolation method (e.g., 'linear',
%       'nearest', 'spline'), defaults to 'linear'
%   extrapVal: Value to use for points outside the image domain, defaults
%       to 0
%
%   Outputs:
%   imWarped: Warped image (same size as displacement field, which can be
%       different from input image size)
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    im (:,:,:,:) {mustBeNumeric}
    imSpacing (1,3) {mustBePositive}
    displ (:,:,:,3,:) {mustBeNumeric}
    displSpacing (1,3) {mustBePositive}
    interpMethod {mustBeTextScalar} = 'linear'
    extrapVal {mustBeNumeric} = 0
end

assert(size(im,4) == size(displ,5), ...
    'Input image and displacement field must have the same number of time frames')

% Coordinate grids for interpolation
xVec = (-floor(size(displ,1)/2) : ceil((size(displ,1)-1)/2)-1) * displSpacing(1);
xVec = cast(xVec, 'like', displ);
yVec = (-floor(size(displ,2)/2) : ceil((size(displ,2)-1)/2)-1) * displSpacing(2);
yVec = cast(yVec, 'like', displ);
zVec = (-floor(size(displ,3)/2) : ceil((size(displ,3)-1)/2)-1) * displSpacing(3);
zVec = cast(zVec, 'like', displ);
[xGrid, yGrid, zGrid] = ndgrid(xVec, yVec, zVec);
x0 = cat(4, xGrid, yGrid, zGrid);

xVec = (-floor(size(im,1)/2) : ceil((size(im,1)-1)/2)-1) * imSpacing(1);
yVec = (-floor(size(im,2)/2) : ceil((size(im,2)-1)/2)-1) * imSpacing(2);
zVec = (-floor(size(im,3)/2) : ceil((size(im,3)-1)/2)-1) * imSpacing(3);

% Interpolate image on warped grid
imWarped = zeros(size(displ, [1:3,5]), 'like', im);
for iTime = 1:size(im, 4)
    xt = x0 + displ(:,:,:,:,iTime);
    imWarped(:,:,:,iTime) = interpn(xVec, yVec, zVec, im(:,:,:,iTime), ...
        xt(:,:,:,1), xt(:,:,:,2), xt(:,:,:,3), interpMethod, extrapVal);
end

end
