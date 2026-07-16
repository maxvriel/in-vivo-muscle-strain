function plotContours(labels, colors, lineWidth)
%plotContours Plots the contours of all unique labels (except 0) in the
%segmentation labels image in the current Axes with the given color
%table and line width.
%
%   Inputs:
%   labels: Segmentation labels as integer indices
%   colors: Color table for the labels, with the labels as zero-based
%       indices
%   lineWidth: Width of the contour lines
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    labels (:,:) {mustBeNumeric, mustBeInteger, mustBeNonnegative}
    colors (:,3) {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(colors, 1)}
    lineWidth (1,1) {mustBeNumeric, mustBePositive} = 1
end

% Find label indices, excluding the background
lab = unique(labels(:));
lab = setdiff(lab, 0);

for iLab = 1:length(lab)
    % Plot boundaries of the current label (there can be multiple disconnected regions)
    boundary = bwboundaries(labels == lab(iLab), 'CoordinateOrder', 'xy');
    for ib = 1:length(boundary)
        plot(boundary{ib}(:,1), boundary{ib}(:,2), 'LineWidth', lineWidth, ...
            'Color', colors(lab(iLab)+1,:))
    end
end

end
