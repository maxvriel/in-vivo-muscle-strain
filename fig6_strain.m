% Creates Figure 7
% Figure 7: Top: Out-of-phase image from the Dixon scan with the contours
% of the segmented muscles overlaid on top for volunteer 1. The red lines
% in each plane show the slices shown in the other two orthogonal planes.
% Bottom: Octahedral shear strain values for dynamic scans 1 and 2 at the
% time point of maximum deformation. The segmentations are also given as
% anatomical reference. Note the absence of strain in the undeformed leg
% and the decreased strain values during the isometric contraction,
% especially in the posterior part of the leg. An animated version of this
% figure over time is available as Supplementary Video SV2. A: Anterior; P:
% Posterior; R: Right; L: Left; H: Head; F: Feet.
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl
clearvars
close all

% Set paths
setup

% Load colormap
load('colormaps/lajolla.mat')

%% Load segmentations
[segmIm, segmLabels, segmSpacing, labelsInfo, labelsColors] = readSegmentation('data/volunteer1/');

% Create mask for segmentation
segmMask = createMask(segmIm, prctile(segmIm(:), 65));
segmMask = imdilate(segmMask, strel('disk', 3));

%% Load data
recFiles = {'recon/volunteer1/dynamic1_recon.h5', 'recon/volunteer1/dynamic2_recon.h5'};

masks = cell(size(recFiles));
segmLabelsReg = cell(size(recFiles));
strains = cell(size(recFiles));
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

    % Remove outer slices
    labelsReg = labelsReg(:,:,8:57);
    mask = mask(:,:,8:57);
    oss = oss(:,:,8:57);
    
    masks{iRec} = mask;
    segmLabelsReg{iRec} = labelsReg;
    strains{iRec} = oss;
end

% Remove outer slices
segmIm = segmIm(:,:,11:40);
segmLabels = segmLabels(:,:,11:40);

%% Show figure
lineColor = 'r';

climsStrain = [0, 0.8];

widths = [0.2, 1, 0.1, 1, 0.1, 0.5, 0.6];
heights = [0.15, 0.5, 0.15, 0.5, 0.1, 0.5, 0.1];
axPos = createAxesPositions(widths, heights);

hFig = figure('Name', 'Strain', 'Color', 'white', 'DefaultAxesFontSize', 14);
hFig.Position(1:3) = [1, 1, 1200];
hFig.Position(4) = hFig.Position(3)*sum(heights)/sum(widths);

axes('Position', combineAxesPositions(axPos(1,:)), 'Visible', 'off')
title('Segmentation', 'FontSize', 16, 'Visible', 'on')

% Transverse view
hAx = axes('Position', axPos{1,1});
imshow(segmIm(:,:,16).')
hold on
plotContours(segmLabels(:,:,16).', labelsColors, 1.5)
daspect(hAx, [1./segmSpacing([1,2]), 1])
xline(98, 'Color', lineColor)
yline(54, 'Color', lineColor)
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 1)-20+[-10,10], 'Y', size(segmIm, 2)-20*[1,1])
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 1)-20*[1,1], 'Y', size(segmIm, 2)-20+[-10,10])
text(size(segmIm, 1)-20+[0,0,-15,15], size(segmIm, 2)-20+[-15,15,0,0], ...
    {'A','P','R','L'}, 'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

% Coronal view
hAx = axes('Position', axPos{1,2});
imshow(squeeze(segmIm(:,54,:)).')
hold on
plotContours(squeeze(segmLabels(:,54,:)).', labelsColors, 1.5)
daspect(hAx, [1./segmSpacing([1,3]), 1])
xline(98, 'Color', lineColor)
yline(16, 'Color', lineColor)
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 1)-20+[-10,10], 'Y', size(segmIm, 3)-5*[1,1])
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 1)-20*[1,1], 'Y', size(segmIm, 3)-5+[-2.5,2.5])
text(size(segmIm, 1)-20+[0,0,-15,15], size(segmIm, 3)-5+[-3.75,3.75,0,0], ...
    {'H','F','R','L'}, 'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

% Sagittal view
hAx = axes('Position', axPos{1,3});
imshow(squeeze(segmIm(98,:,:)).')
hold on
plotContours(squeeze(segmLabels(98,:,:)).', labelsColors, 1.5)
daspect(hAx, [1./segmSpacing([2,3]), 1])
xline(54, 'Color', lineColor)
yline(16, 'Color', lineColor)
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 2)-20+[-10,10], 'Y', size(segmIm, 3)-5*[1,1])
hl = annotation(hFig, 'line', 'Color', 'w');
set(hl, 'Parent', hAx, 'X', size(segmIm, 2)-20*[1,1], 'Y', size(segmIm, 3)-5+[-2.5,2.5])
text(size(segmIm, 2)-20+[0,0,-15,15], size(segmIm, 3)-5+[-3.75,3.75,0,0], ...
    {'H','F','A','P'}, 'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

% Show muscle names
legendIdx = [32,34,36,38,40,44,48,56,58,60,62,64,104].';
hP = fill(NaN(3,length(legendIdx)), NaN(3,length(legendIdx)), 0, 'EdgeColor', 'none');
set(hP, {'FaceColor'}, num2cell(labelsColors(legendIdx+1,:), 2))
[~, labelIdx] = ismember(legendIdx, labelsInfo.IDX);
idxLabels = labelsInfo.LABEL(labelIdx);
idxNames = regexpi(idxLabels, '\d+ (.*) (?:Left|Right)?', 'tokens', 'once');
idxNames = cat(1, idxNames{:});
idxNames = replace(idxNames, 'Satorius', 'Sartorius');
hLeg = legend(hP, idxNames, 'FontSize', 10);
hLeg.Position(1) = sum(axPos{1,3}([1,3])) + 0.1*widths(end)/sum(widths);
hLeg.Position(2) = axPos{1,3}(2) + 0.5*axPos{1,3}(4) - 0.5*hLeg.Position(4);

axes('Position', combineAxesPositions(axPos(2,:)), 'Visible', 'off')
title('Strain', 'FontSize', 16, 'Visible', 'on')

for i = 1:2
    % Transverse view
    axes('Position', axPos{i+1,1})
    imshow(maskImage(strains{i}(:,:,25).', masks{i}(:,:,25).', ...
        climsStrain, lajolla, [0.5,0.5,0.5]));
    hold on
    plotContours(segmLabelsReg{i}(:,:,25).', labelsColors, 1)
    yticks([])
    if i == 1
        ylabel({'Passive', 'muscles'}, 'FontSize', 14)
    elseif i == 2
        ylabel({'Isometric', 'contractions'}, 'FontSize', 14)
    end

    % Coronal view
    axes('Position', axPos{i+1,2})
    imshow(maskImage(squeeze(strains{i}(:,26,:)).', squeeze(masks{i}(:,26,:)).', ...
        climsStrain, lajolla, [0.5,0.5,0.5]));
    hold on
    plotContours(squeeze(segmLabelsReg{i}(:,26,:)).', labelsColors, 1)
    
    % Sagittal view
    axes('Position', axPos{i+1,3})
    imshow(maskImage(squeeze(strains{i}(46,:,:)).', squeeze(masks{i}(46,:,:)).', ...
        climsStrain, lajolla, [0.5,0.5,0.5]));
    hold on
    plotContours(squeeze(segmLabelsReg{i}(46,:,:)).', labelsColors, 1)
end

hAx = axes('Position', axPos{3,3}, 'Visible', 'off');
clim(hAx, climsStrain)
colormap(hAx, lajolla)
hCb = colorbar();
hCb.Position(1) = sum(axPos{3,3}([1,3])) + 0.1*widths(end)/sum(widths);
hCb.Position(2) = axPos{3,3}(2);
hCb.Position(3) = 0.1*widths(end)/sum(widths);
hCb.Position(4) = axPos{3,3}(4)+axPos{2,3}(4)+heights(5)/sum(heights);
ylabel(hCb, 'OSS (-)', 'FontSize', 14, 'Rotation', -90)
