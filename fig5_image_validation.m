% Creates Figure 5
% Figure 5: Validation of the time-resolved dynamic image series. During
% the dynamic scan, the pressure was increased and decreased continuously,
% while the pressure was fixed during the 22-second acquisition of each
% validation image. The SSIM values were determined between images from
% corresponding time points. The graph at the bottom shows the SSIM values
% for all nine validation time points, across all eight volunteers.
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl
clearvars
close all

% Set paths
setup

%% Calculate SSIM values
subjects = [1,2,4:9];
nVal = 9;
ssimVals = zeros(nVal, length(subjects));
ssimRef = zeros(nVal, length(subjects));
plotSubject = 1;
for iSubj = 1:length(subjects)
    recFile = sprintf('recon/volunteer%d/dynamic1_recon.h5', subjects(iSubj));
    valFile = sprintf('data/volunteer%d/validation.mrd', subjects(iSubj));
    dynFile = sprintf('data/volunteer%d/dynamic1.mrd', subjects(iSubj));

    % Load reconstruction
    [dynIm, dt, FOV, vel] = readRecon(recFile);
    climsIm = [0, max(abs(dynIm(:)))];

    % Integrate velocity field to displacement field
    spacing = FOV ./ size(dynIm, 1:3)';
    displ = velocityToDisplacement(vel, spacing, dt, 1);

    % Find time frame of maximum deformation
    [~, tMax] = max(vecnorm(reshape(displ, [], size(displ,5))));

    % Load validation and dynamic data
    [valData, valHeader, valCsm] = readMrd(valFile);
    [dynData, dynHeader, dynCsm] = readMrd(dynFile);

    % Reconstruct validation images
    valIm = imageRecon(valData, valHeader, valCsm);

    % Remove slice oversampling
    nKz = dynHeader.encoding(1).encodedSpace.matrixSize.z;
    nZ = dynHeader.encoding(1).reconSpace.matrixSize.z;
    zIdx = floor((nKz-nZ)/2) + (1:nZ);
    dynIm = dynIm(:,:,zIdx,:);

    % Select right leg
    dynIm = dynIm((1:64)+4,:,:,:);
    valIm = valIm((1:64)+4,:,:,:);

    % Take magnitude of the images
    dynIm = abs(dynIm);
    valIm = abs(valIm);

    % Coil compression
    nVirtCoils = 8;

    allData = cat(1, valData.data{2:end});
    [~, ~, ccMat] = svd(allData, 'econ', 'vector');

    valData.data = cellfun(@(x) x*ccMat(:,1:nVirtCoils), valData.data, 'UniformOutput', false);
    valData.head.active_channels(:) = nVirtCoils;
    valHeader.acquisitionSystemInformation.receiverChannels = nVirtCoils;
    valCsm = tensorprod(valCsm, ccMat(:,1:nVirtCoils), 4, 1);

    % Find k0 projections for the validation data
    ky0 = valHeader.encoding(1).encodingLimits.kspace_encoding_step_1.center;
    kz0 = valHeader.encoding(1).encodingLimits.kspace_encoding_step_2.center;
    isNotNoise = ~valData.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
    isKy0 = valData.head.idx.kspace_encode_step_1 == ky0;
    isKz0 = valData.head.idx.kspace_encode_step_2 == kz0;
    valProj = valData.data(isNotNoise & isKy0 & isKz0);
    valProj = cat(3, valProj{:});
    valProj = fftshift(ifft(ifftshift(valProj,1), [], 1), 1) * sqrt(size(valProj,1));
    valProj = reshape(valProj, size(valProj,1)*nVirtCoils, []);

    % Coil compression (use same compression matrix)
    dynData.data = cellfun(@(x) x*ccMat(:,1:nVirtCoils), dynData.data, 'UniformOutput', false);
    dynData.head.active_channels(:) = nVirtCoils;
    dynHeader.acquisitionSystemInformation.receiverChannels = nVirtCoils;
    dynCsm = tensorprod(dynCsm, ccMat(:,1:nVirtCoils), 4, 1);

    % Find k0 projections for the dynamic data (use the second repetition
    % and average across four consecutive shots)
    ky0 = dynHeader.encoding(1).encodingLimits.kspace_encoding_step_1.center;
    kz0 = dynHeader.encoding(1).encodingLimits.kspace_encoding_step_2.center;
    isNotNoise = ~dynData.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
    isKy0 = dynData.head.idx.kspace_encode_step_1 == ky0;
    isKz0 = dynData.head.idx.kspace_encode_step_2 == kz0;
    isRep2 = dynData.head.idx.user(1, :) == 2;
    dynProj = dynData.data(isNotNoise & isKy0 & isKz0 & isRep2);
    dynProj = cat(3, dynProj{:});
    dynProj = reshape(dynProj, size(dynProj,1), nVirtCoils, 4, []);
    dynProj = mean(dynProj, 3);
    dynProj = fftshift(ifft(ifftshift(dynProj,1), [], 1), 1) * sqrt(size(dynProj,1));
    dynProj = reshape(dynProj, size(dynProj,1)*nVirtCoils, []);

    % Calculate correlation matrix
    valProjCentered = valProj - mean(valProj);
    dynProjCentered = dynProj - mean(dynProj);
    corrMat = valProjCentered' * dynProjCentered;
    corrMat = corrMat ./ vecnorm(valProjCentered).';
    corrMat = corrMat ./ vecnorm(dynProjCentered);

    % Determine the time frame of maximum correlation for each validation image
    frameIdx = zeros(nVal, 1);
    [~, frameIdx(1:floor(nVal/2))] = max(corrMat(1:floor(nVal/2),1:tMax), [], 2);
    [~, frameIdx(floor(nVal/2)+1)] = max(corrMat(floor(nVal/2)+1,:), [], 2);
    [~, frameIdx(floor(nVal/2)+2:end)] = max(corrMat(floor(nVal/2)+2:end,tMax:end), [], 2);
    frameIdx(floor(nVal/2)+2:end) = frameIdx(floor(nVal/2)+2:end) + tMax-1;

    % Sample dynamic image frames
    dynImSampled = dynIm(:,:,:,frameIdx);

    % Elastix registration
    valImReg = valIm;
    for iVal = 1:nVal
        frameDir = sprintf('./registration/tmp/validation/volunteer%d/frame%d', subjects(iSubj), iVal);
        if ~isfolder(frameDir)
            mask = createMask(dynImSampled(:,:,:,iVal), prctile(dynIm(:), 65));
            elastix(valIm(:,:,:,iVal), spacing, dynImSampled(:,:,:,iVal), ...
                spacing, frameDir, './registration/elastixParamsRigid.txt', mask);
        end
        [~, valImReg(:,:,:,iVal)] = transformix(fullfile(frameDir, 'TransformParameters.0.txt'));
    end

    % Calculate SSIM values
    ssimVals(:,iSubj) = squeeze(ssim(dynImSampled, valImReg, ...
        'DynamicRange', climsIm(2), 'DataFormat', 'SSSB'));
    ssimRef(:,iSubj) = squeeze(ssim(repmat(dynIm(:,:,:,1), 1, 1, 1, nVal), ...
        valImReg, 'DynamicRange', climsIm(2), 'DataFormat', 'SSSB'));

    % Save images for plotting
    if subjects(iSubj) == plotSubject
        plotIm = dynImSampled(:,:,26,5).';
        plotImRef = valImReg(:,:,26,5).';
    end
end

%% Show figure
valColor = [252, 96, 57]/255;
dyColor = [17, 145, 250]/255;

cmap = [0.1216    0.4667    0.7059;
        1.0000    0.4980    0.0549;
        0.1725    0.6275    0.1725;
        0.8392    0.1529    0.1569;
        0.5804    0.4039    0.7412;
        0.5490    0.3373    0.2941;
        0.8902    0.4667    0.7608;
        0.4980    0.4980    0.4980;
        0.7373    0.7412    0.1333;
        0.0902    0.7451    0.8118];

widths = [0.4, 1.5, 0.2, 1.5, 0.2];
heights = [0.4, 1, 0.4, 1, 0.2, 1.5, 0.3];
axPos = createAxesPositions(widths, heights);

hFig = figure('Name', 'Image Validation', 'Color', 'white', 'DefaultAxesFontSize', 12);
hFig.Position(1:3) = [1, 1, 600];
hFig.Position(4) = round(hFig.Position(3)*sum(heights)/sum(widths));

x = linspace(-pi/2, pi/2, 101).';
x = x(1:end-1);
y = [40*sin(x)+40;
    -40*sin(x)+40;
    0];

hAx = axes('Position', axPos{1,1});
hold on
plot(0.01*(0:length(y)-1), y, 'k', 'LineWidth', 2)
hold on
xline(1, 'Color', dyColor, 'LineWidth', 2, 'Alpha', 1)
xlim tight
ylim([0, 90])
xlabel('Time', 'FontSize', 12)
ylabel('Pressure (mmHg)', 'FontSize', 12)
xticks([])
yticks(0:20:80)
grid on
set(hAx, 'XGrid', 'off')
box off
title('Dynamic data', 'FontSize', 12)

hAx = axes('Position', axPos{2,1});
imshow(plotIm);
clim([0, Inf]);
fig1Pos = tightPosition(hAx);

annotation('rectangle', fig1Pos, 'Color', dyColor, 'LineWidth', 4, 'PickableParts', 'none');

y = [zeros(100,1);
    10*sin(x)+10;
    20*ones(100,1);
    10*sin(x)+30;
    40*ones(100,1);
    10*sin(x)+50;
    60*ones(100,1);
    10*sin(x)+70;
    80*ones(100,1);
    -10*sin(x)+70;
    60*ones(100,1);
    -10*sin(x)+50;
    40*ones(100,1);
    -10*sin(x)+30;
    20*ones(100,1);
    -10*sin(x)+10;
    zeros(100,1);
    0];

hAx = axes('Position', axPos{1,2});
patch([8;8;9;9], [0;100;100;0], valColor, 'FaceAlpha', 1, 'EdgeColor', 'none');
hold on
plot(0.01*(0:length(y)-1), y, 'k', 'LineWidth', 2)
xlim tight
ylim([0, 90])
xl = xlabel('Validation image index', 'FontSize', 12);
xl.Position(2) = xl.Position(2) - 4;
xticks((0:8)*2+0.5)
xticklabels(string(1:9))
yticks(0:20:80)
grid on
set(hAx, 'XGrid', 'off')
box off
title('Validation data', 'FontSize', 12)

hAx = axes('Position', axPos{2,2});
imshow(plotImRef);
clim([0, Inf]);
annotation('rectangle', tightPosition(hAx), 'Color', valColor, ...
    'LineWidth', 4, 'PickableParts', 'none');
fig2Pos = tightPosition(hAx);

annotation('doublearrow', [sum(fig1Pos([1,3])), fig2Pos(1)]+0.015*[1,-1], ...
    fig1Pos(2)+fig1Pos(4)*0.5*[1,1], 'LineWidth', 1.5)
annotation('textbox', [sum([fig1Pos([1,3])-0.1/2, fig2Pos(1)])/2, fig1Pos(2)+fig1Pos(4)/2, 0.1, 0.05], ...
    'String', 'SSIM', 'FontSize', 12, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'EdgeColor', 'none')

axes('Position', combineAxesPositions(axPos(3,:)))
hl = plot(ssimVals, 'LineWidth', 2);
set(hl, {'Color'}, num2cell(cmap(1:length(subjects),:), 2))
hold on
hl2 = plot(ssimRef, '--', 'LineWidth', 2);
set(hl2, {'Color'}, num2cell(cmap(1:length(subjects),:), 2))
xlim tight
xticks(1:size(ssimVals,1))
ylim([0.5, 1])
grid on
xlabel('Validation image index', 'FontSize', 12)
ylabel('SSIM', 'FontSize', 12)

hlVal = plot(NaN, NaN, 'k', 'LineWidth', 2);
hlRef = plot(NaN, NaN, 'k--', 'LineWidth', 2);
hlBlank = plot(NaN, NaN, 'w', 'LineWidth', 2);
legend([hl; hlBlank; hlVal; hlRef], ["Volunteer " + (1:length(subjects)), "", "Validation", "Static reference"], ...
    'Location', 'eastoutside')

axes('Position', combineAxesPositions(axPos(1,:)), 'Visible', 'off');
ht = title('Dynamic image validation', 'FontSize', 14, 'Visible', 'on');
ht.Position(2) = 1.2;

exportgraphics(hFig, './figure/figure5.eps', 'ContentType', 'vector', 'Padding', 'figure')
exportgraphics(hFig, './figure/figure5.png', 'Resolution', 300, 'Padding', 'figure')
