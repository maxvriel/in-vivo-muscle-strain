function [segmIm, segmLabels, segmSpacing, labelsInfo, labelsColors] = readSegmentation(folderName)
%readSegmentation Load segmentation data.
%
%   Inputs:
%   folderName: Path to the folder containing segmentation files
%
%   Outputs:
%   segmIm: Out-of-phase image used for segmentation
%   segmLabels: Segmentation labels as integer indices
%   segmSpacing: Spacing of the segmentation image in mm
%   labelsInfo: Table with information about the segmentation labels
%   labelsColors: Colors for each segmentation label
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    folderName {mustBeTextScalar, mustBeFolder}
end

% Read out-of-phase image and segmentation labels
segmImInfo = niftiinfo(fullfile(folderName, 'outphase.nii.gz'));
segmSpacing = segmImInfo.PixelDimensions;
segmIm = niftiread(segmImInfo);
segmLabels = niftiread(fullfile(folderName, 'segmentation.nii.gz'));

% Flip images and remove black areas
segmIm = flip(flip(flip(segmIm, 1), 2), 3);
segmIm = segmIm(:,73:216,:);
segmLabels = flip(flip(flip(segmLabels, 1), 2), 3);
segmLabels = segmLabels(:,73:216,:);

% Normalize image
segmIm = single(segmIm);
segmIm = segmIm / max(segmIm(:));

% Parse labels info
labelsInfo = readtable(fullfile(folderName, 'segmentation.txt'), 'CommentStyle', '#');
labelsInfo = renamevars(labelsInfo, 1:width(labelsInfo), ...
    {'IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'});
labelsColors = zeros(max(labelsInfo.IDX)+1, 3);
labelsColors(labelsInfo.IDX+1,:) = labelsInfo{:, {'R', 'G', 'B'}} / 255;

end
