function im = maskImage(im, mask, clims, cmap, bgColor)
%maskImage Masks out the background of an image
%
%   Inputs:
%   im: Input image
%   mask: Binary mask to imprint on the image
%   clims: Color limits for the image
%   cmap: Colormap for the image
%   bgColor: Background color for the masked-out regions
%
%   Outputs:
%   im: Truecolor RGB image with the mask imprinted
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    im {mustBeNumeric}
    mask {mustBeNumericOrLogical}
    clims (1,2) {mustBeNumeric} = [0, 1]
    cmap (:,3) {mustBeNumeric} = gray(256)
    bgColor (1,3) {mustBeNumeric} = [0, 0, 0]
end

% Create colormap bins
cmapBins = linspace(clims(1), clims(2), size(cmap,1)+1);
cmapBins(1) = -Inf;
cmapBins(end) = Inf;

% Create truecolor RGB image given the colormap
imSz = size(im);
im = discretize(im, cmapBins);
im = ind2rgb(im(:), cmap);

% Imprint the mask
im(~mask(:),:) = repmat(bgColor, nnz(~mask), 1);
im = reshape(im, [imSz, 3]);

end
