% Creates Figure 8
% Figure 8: From left to right: Out-of-phase image from the Dixon scan with
% the contours of the segmented muscles overlaid on top for volunteer 1,
% mean octahedral shear strain values at the maximum deformation for each
% segmented muscle in dynamic scans 1 and 2, and muscle-wise differences in
% mean octahedral shear strain value. Note that the biggest decrease
% appears in the posterior part of the leg, where the hamstring muscles are
% located, which are the main knee flexors.
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl
clearvars
close all

% Set paths
setup

% Load colormap
load('colormaps/lajolla.mat')
load('colormaps/vik.mat')

%% Load segmentations
[segmIm, segmLabels, segmSpacing, labelsInfo, labelsColors] = readSegmentation('data/volunteer1/');

% Create mask for segmentation
segmMask = createMask(segmIm, prctile(segmIm(:), 65));
segmMask = imdilate(segmMask, strel('disk', 3));

%% Load data
recFiles = {'recon/volunteer1/dynamic1_recon.h5', ...
    'recon/volunteer1/dynamic2_recon.h5'};

masks = cell(size(recFiles));
segmLabelsReg = cell(size(recFiles));
strainTables = cell(size(recFiles));
for iRec = 1:length(recFiles)
    % Load reconstruction
    [dynIm, dt, FOV, vel] = readRecon(recFiles{iRec});
    
    % Integrate velocity field to displacement field
    spacing = FOV ./ size(dynIm, 1:3)';
    displ = velocityToDisplacement(vel, spacing, dt, 1);

    % Find time frame of maximum deformation
    [~, tMax] = max(vecnorm(reshape(displ, [], size(displ,5))));

    % Calculate octahedral shear strain
    oss = calcOctahedralShearStrain(displ(:,:,:,:,tMax), spacing);

    % Create mask for first time frame
    imInit = abs(dynIm(:,:,:,1));
    imInit = imInit / max(imInit(:));
    mask = createMask(imInit, prctile(imInit(:), 65));

    % Elastix registration
    elastixDir = sprintf('./registration/tmp/segmentation/volunteer1/rec%d', iRec);
    if ~isfolder(elastixDir)
        % Create mask to use during registration
        maskReg = mask;
        maskReg = imdilate(maskReg, strel('disk', 3));
        maskReg(:,:,[1:7, 58:64]) = false;
        
        % Run registration
        elastix(segmIm, segmSpacing, imInit, spacing, elastixDir, ...
            {'./registration/elastixParamsSpline1.txt', './registration/elastixParamsSpline2.txt'}, ...
            maskReg, segmMask);
    end

    % Warp segmentation labels
    displReg = transformix(fullfile(elastixDir, 'TransformParameters.1.txt'));
    labelsReg = warpImage(segmLabels, segmSpacing, displReg, spacing, 'nearest', 0);

    % Remove outer slices and select right leg
    labelsReg = labelsReg((1:64)+4,:,8:57);
    mask = mask((1:64)+4,:,8:57);
    oss = oss((1:64)+4,:,8:57);

    % Calculate mean strain values for each segmented muscle
    [strainInfo, strainLabels] = groupsummary(oss(:), labelsReg(:), 'mean');
    [~, labelIdx] = ismember(strainLabels, labelsInfo.IDX);
    strainTable = table(strainLabels, labelsInfo.LABEL(labelIdx), ...
        strainInfo(:,1), 'VariableNames', {'Index', 'Label', 'MeanStrain'});
    
    masks{iRec} = mask;
    segmLabelsReg{iRec} = labelsReg;
    strainTables{iRec} = strainTable;
end

% Remove outer slices
segmIm = segmIm(:,:,11:40);
segmLabels = segmLabels(:,:,11:40);

%% Show figure
clims = [0, 0.3];
climsDiff = [-0.25, 0.25];

widths = [0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1];
heights = [0.1, 0.5, 0.2];
axPos = createAxesPositions(widths, heights);

hFig = figure('Name', 'Strain Difference 1', 'Color', 'white', 'DefaultAxesFontSize', 14);

hFig.Position(1:3) = [1, 1, 1200];
hFig.Position(4) = hFig.Position(3)*sum(heights)/sum(widths);

hAx = axes('Position', axPos{1,1});
imshow(abs(segmIm(1:end/2,:,10)).')
hold on
plotContours(segmLabels(1:end/2,:,10).', labelsColors, 1.5)
daspect(hAx, [1./segmSpacing([1,2]), 1])
title('Segmentation', 'FontSize', 16)
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 1)/2-20+[-10,10], 'Y', size(segmIm, 2)-20*[1,1])
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 1)/2-20*[1,1], 'Y', size(segmIm, 2)-20+[-10,10])
text(size(segmIm, 1)/2-20+[0,0,-15,15], size(segmIm, 2)-20+[-15,15,0,0], ...
    {'A','P','R','L'}, 'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

[~, labelIdx1] = ismember(segmLabelsReg{1}, strainTables{1}.Index);
[~, labelIdx2] = ismember(segmLabelsReg{1}, strainTables{2}.Index);
labelIdx1 = labelIdx1 .* (segmLabelsReg{1} ~= 0);
labelIdx2 = labelIdx2 .* (segmLabelsReg{1} ~= 0);
strainR = zeros(size(segmLabelsReg{1}));
strainR(labelIdx1 > 0) = strainTables{1}.MeanStrain(labelIdx1(labelIdx1 > 0));
strainC = zeros(size(segmLabelsReg{1}));
strainC(labelIdx2 > 0) = strainTables{2}.MeanStrain(labelIdx2(labelIdx2 > 0));
diffIdx = labelIdx1 > 0 & labelIdx2 > 0;
meanStrain = zeros(size(segmLabelsReg{1}));
meanStrain(diffIdx) = (strainTables{2}.MeanStrain(labelIdx2(diffIdx)) - ...
    strainTables{1}.MeanStrain(labelIdx1(diffIdx)));

axes('Position', axPos{1,2})
imshow(maskImage(strainR(:,:,25).', strainR(:,:,25).' ~= 0, clims, lajolla, [0.5,0.5,0.5]))
title('Passive muscles', 'FontSize', 16)

hAx = axes('Position', axPos{1,3});
imshow(maskImage(strainC(:,:,25).', strainC(:,:,25).' ~= 0, clims, lajolla, [0.5,0.5,0.5]))
title('Isometric contractions', 'FontSize', 16)
clim(hAx, clims)
colormap(hAx, lajolla)
cbPos = combineAxesPositions(axPos(1,2:3));
cbPos(2) = cbPos(2) - 0.4*heights(end)/sum(heights);
cbPos(4) = 0.05;
cb = colorbar(hAx, 'Location', 'southoutside', 'Position', cbPos, 'FontSize', 12);
yl = ylabel(cb, 'OSS (-)', 'FontSize', 12);
yl.Position(2) = yl.Position(2) * 1.4;

hAx = axes('Position', axPos{1,4});
imshow(maskImage(meanStrain(:,:,25).', meanStrain(:,:,25).' ~= 0, climsDiff, vik, [0.5,0.5,0.5]))
title('Difference', 'FontSize', 16)
clim(hAx, climsDiff)
colormap(hAx, vik)
cbPos = axPos{1,4};
cbPos(2) = cbPos(2) - 0.4*heights(end)/sum(heights);
cbPos(4) = 0.05;
cb = colorbar(hAx, 'Location', 'southoutside', 'Position', cbPos, 'FontSize', 12);
yl = ylabel(cb, '\DeltaOSS (-)', 'FontSize', 12);
yl.Position(2) = yl.Position(2) * 1.4;

annotation('textbox', ...
    [0.5*(sum(axPos{1,2}([1,3]))+axPos{1,3}(1))-0.5*0.05, sum(axPos{1,2}([2,4]).*[1,0.5])-0.5*0.05, 0.05, 0.05], ...
    'String', char(8211), 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
annotation('textbox', ...
    [0.5*(sum(axPos{1,3}([1,3]))+axPos{1,4}(1))-0.5*0.05, sum(axPos{1,3}([2,4]).*[1,0.5])-0.5*0.05, 0.05, 0.05], ...
    'String', '=', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')

exportgraphics(hFig, './figures/figure8.eps', 'ContentType', 'vector', 'Padding', 'figure')
exportgraphics(hFig, './figures/figure8.png', 'Resolution', 300, 'Padding', 'figure')
