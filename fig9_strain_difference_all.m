% Creates Figure 9
% Figure 9: Mean octahedral shear strain difference in the central
% transverse slice at the maximum deformation over each segmented muscle
% between the scan with passive muscles and the scan during the isometric
% contraction. Note that in all cases, the biggest decrease appears in the
% posterior part of the leg, where the hamstring muscles are located, which
% are the main knee flexors.
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

%% Load data
subjects = [1,2,4:9];
segmLabelsReg = cell(4,length(subjects));
strainTables = cell(4,length(subjects));
for iSubj = 1:length(subjects)
    % Load segmentation
    segmFolder = sprintf('data/volunteer%d/', subjects(iSubj));
    [segmIm, segmLabels, segmSpacing, labelsInfo, labelsColors] = readSegmentation(segmFolder);

    % Create mask for segmentation
    segmMask = createMask(segmIm, prctile(segmIm(:), 65));
    segmMask = imdilate(segmMask, strel('disk', 3));

    for iDyn1 = 1:4
        % Load reconstruction
        recFile = sprintf('recon/volunteer%d/dynamic%d_recon.h5', subjects(iSubj), iDyn1);
        [dynIm, dt, FOV, vel] = readRecon(recFile);
        
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
        elastixDir = sprintf('./registration/tmp/segmentation/volunteer%d/dynamicscan%d', subjects(iSubj), iDyn1);
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
        
        segmLabelsReg{iDyn1,iSubj} = labelsReg;
        strainTables{iDyn1,iSubj} = strainTable;
    end
end

%% Show figure
clims = [-0.25, 0.25];
cmap = vik;

widths = [0.15, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.4];
heights = [0.2, 0.5, 0.1, 0.5, 0.2, 0.5, 0.1, 0.5, 0.1];
axPos = createAxesPositions(widths, heights);

hFig = figure('Name', 'Strain Difference', 'Color', 'white', 'DefaultAxesFontSize', 14);

hFig.Position(1:3) = [1, 1, 800];
hFig.Position(4) = hFig.Position(3)*sum(heights)/sum(widths);

for iSubj = 1:length(subjects)
    for iCuffPos = 1:2
        iDyn1 = (iCuffPos-1)*2+1;
        iDyn2 = iDyn1+1;

        [~, labelIdx1] = ismember(segmLabelsReg{iDyn1,iSubj}, strainTables{iDyn1,iSubj}.Index);
        [~, labelIdx2] = ismember(segmLabelsReg{iDyn1,iSubj}, strainTables{iDyn2,iSubj}.Index);
        labelIdx1 = labelIdx1 .* (segmLabelsReg{iDyn1,iSubj} ~= 0);
        labelIdx2 = labelIdx2 .* (segmLabelsReg{iDyn1,iSubj} ~= 0);
        diffIdx = labelIdx1 > 0 & labelIdx2 > 0;
        meanStrain = zeros(size(segmLabelsReg{iDyn1,iSubj}));
        meanStrain(diffIdx) = (strainTables{iDyn2,iSubj}.MeanStrain(labelIdx2(diffIdx)) - ...
            strainTables{iDyn1,iSubj}.MeanStrain(labelIdx1(diffIdx)));

        iRow = (iSubj > size(axPos,2))*2 + iCuffPos;
        iCol = mod(iSubj-1, size(axPos,2)) + 1;
        axes('Position', axPos{iRow,iCol})

        imshow(maskImage(meanStrain(:,:,25).', meanStrain(:,:,25).' ~= 0, clims, vik, [0.5,0.5,0.5]))
        yticks([])

        if iCuffPos == 1
            title(sprintf('Volunteer #%d', iSubj), 'FontSize', 14)
        end
        if iCol == 1
            ylabel(sprintf('Cuff position %d', iCuffPos), 'FontSize', 14)
        end

    end
end

hAx = axes('Position', axPos{3,3}, 'Visible', 'off');
clim(hAx, clims)
colormap(hAx, vik)
cb = colorbar('FontSize', 11);
cb.Position = combineAxesPositions(axPos(1:4,end));
cb.Position(1) = cb.Position(1) + cb.Position(3) + 0.2*widths(end)/sum(widths);
cb.Position(3) = 0.015;
ylabel(cb, '\DeltaOSS (-)', 'FontSize', 12, 'Rotation', -90)

hAx = axes('Position', combineAxesPositions(axPos(1,:)), 'Visible', 'off');
ht = title(hAx, 'Strain difference', 'FontSize', 16, 'Visible', 'on');
ht.Position(2) = ht.Position(2) + 0.16;

annotation('line', [0, 1] + ([widths(1), -widths(end)] + 0*[-0.5,0.5]*widths(1))/sum(widths), ...
    axPos{2,1}(2) - 0.4*heights(5)/sum(heights) * [1,1], ...
    'LineWidth', 1.5)
