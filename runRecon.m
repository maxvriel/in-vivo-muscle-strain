% Run the joint reconstruction of the time-resolved images and velocity
% fields, as described in the manuscript "Time-resolved 3D imaging and
% strain analysis for in vivo muscle dynamics"
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl
clearvars
close all

% Set paths
setup

dataPath = './data';
reconPath = './recon';

%% Recon settings
subjectIdx = 1;
dynScanIdx = 1;

nResLevels = 2;
nIterations = [15, 5];
nReadoutsPerFrame = 4*32;
lambda = 100;
mu = 10;

saveResult = true;
showResult = true;

%% Load data
dataFile = sprintf('%s/volunteer%d/dynamic%d.mrd', dataPath, subjectIdx, dynScanIdx);
[data, header, csm] = readMrd(dataFile);

fieldOfView = header.encoding(1).encodedSpace.fieldOfView_mm;
fieldOfView = [fieldOfView.x, fieldOfView.y, fieldOfView.z];
timeStep = nReadoutsPerFrame * header.sequenceParameters.TR * 1e-3;

nx = header.encoding(1).encodedSpace.matrixSize.x;
ny = header.encoding(1).encodedSpace.matrixSize.y;
nz = header.encoding(1).encodedSpace.matrixSize.z;

%% Coil compression
nVirtCoils = 8;
calibSize = [24, 24, 24];
kxCalib = floor(nx/2) - floor(calibSize(1)/2) + (1:calibSize(1));
kyCalib = floor(ny/2) - floor(calibSize(2)/2) + (1:calibSize(2));
kzCalib = floor(nz/2) - floor(calibSize(3)/2) + (1:calibSize(3));
ccIdx = ~data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
ccIdx = ccIdx & ismember(data.head.idx.kspace_encode_step_1, kxCalib);
ccData = cellfun(@(d) d(kxCalib,:), data.data(ccIdx), 'UniformOutput', false);
ccData = cat(1, ccData{:});
[~, ~, ccMat] = svd(ccData, 'econ', 'vector');

data.data = cellfun(@(x) x*ccMat(:,1:nVirtCoils), data.data, 'UniformOutput', false);
data.head.active_channels(:) = nVirtCoils;
header.acquisitionSystemInformation.receiverChannels = nVirtCoils;
csm = tensorprod(csm, ccMat(:,1:nVirtCoils), 4, 1);

%% Create spline basis
% Create spatial coordinate grid
nDims = 3;
spacing = fieldOfView ./ [nx, ny, nz];
xVec = (-floor(nx/2) : ceil(nx/2)-1) * spacing(1);
yVec = (-floor(ny/2) : ceil(ny/2)-1) * spacing(2);
zVec = (-floor(nz/2) : ceil(nz/2)-1) * spacing(3);

% Create B-spline basis functions
splOrder = [4, 4, 4];
nKnots = [24, 12, 12];
splineBasis = cell(1,nDims);
xKnots = [-ones(1, splOrder(1)-2), linspace(-1,1,nKnots(1)), ones(1, splOrder(1)-2)] * (fieldOfView(1)/2);
yKnots = [-ones(1, splOrder(2)-2), linspace(-1,1,nKnots(2)), ones(1, splOrder(2)-2)] * (fieldOfView(2)/2);
zKnots = [-ones(1, splOrder(3)-2), linspace(-1,1,nKnots(3)), ones(1, splOrder(3)-2)] * (fieldOfView(3)/2);
nBx = length(xKnots) - splOrder(1);
nBy = length(yKnots) - splOrder(2);
nBz = length(zKnots) - splOrder(3);
Bx = zeros(nx, nBx);
By = zeros(ny, nBy);
Bz = zeros(nz, nBz);
for iSpline = 1:nBx
    controlPoints  = zeros(1, nBx);
    controlPoints(iSpline) = 1;
    Bx(:,iSpline) = evalBSpline(xKnots, controlPoints, xVec);
end
for iSpline = 1:nBy
    controlPoints = zeros(1, nBy);
    controlPoints(iSpline) = 1;
    By(:,iSpline) = evalBSpline(yKnots, controlPoints, yVec);
end
for iSpline = 1:nBz
    controlPoints = zeros(1, nBz);
    controlPoints(iSpline) = 1;
    Bz(:,iSpline) = evalBSpline(zKnots, controlPoints, zVec);
end

% Find spline coefficients that are in the foreground region
[~, maxBx] = max(Bx);
[~, maxBy] = max(By);
[~, maxBz] = max(Bz);
maxSplIdx = reshape(maxBz-1, 1, 1, nBz)*nx*ny + (maxBy-1)*nx + maxBx.';
csmRms = sqrt(sum(abs(csm).^2, 4));
fgIdx = csmRms(maxSplIdx) > 0;

%% Divide data into time frames
% For Dynamic scan 1, use the second pressure cycle
if dynScanIdx == 1
    userIdx = 2;
else
    userIdx = 0;
end

% Use repetition index to indicate time frames
readoutMask = ~data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT') & ...
    data.head.idx.user(1,:) == userIdx;
nReadouts = nnz(readoutMask);
nTime = ceil(nReadouts/nReadoutsPerFrame);
repIdx = repelem(1:nTime, nReadoutsPerFrame);
repIdx = repIdx(1:nReadouts);

data.head.idx.repetition(:) = 0;
data.head.idx.repetition(readoutMask) = repIdx;

%% Sort data into a k-space array and apply density compensation
% The density compensation is performed by the averaging in sortData
[sortedData, samplCount] = sortData(data, header, true, 1:nTime);
samplMask = samplCount > 0;

% Add virtual "zero" measurements in k-space corners as additional regularization
[Y,Z] = ndgrid(1:ny, 1:nz);
r = [ny, nz] / 2;
k0 = floor(r)+1;
shutter = (((Y-k0(1))/r(1)).^2+((Z-k0(2))/r(2)).^2) <= 1;
samplMask = reshape(samplMask, 1, ny*nz, 1, nTime);
samplMask(1, ~shutter(:), 1, :) = 1;
samplMask = reshape(samplMask, 1, ny, nz, 1, nTime);

%% Run joint reconstruction
velSplEst = zeros(nBx*nBy*nBz, nDims, nTime-1, 'single');
for iResLevel = 1:nResLevels
    % Determine data size for current resolution level
    samplFactor = 2^(nResLevels-iResLevel);
    nxRes = ceil(nx / samplFactor);
    nyRes = ceil(ny / samplFactor);
    nzRes = ceil(nz / samplFactor);
    spacingRes = fieldOfView ./ [nxRes, nyRes, nzRes];
    imEst = zeros(nxRes, nyRes, nzRes, nTime, 'single');

    % Subsample data and coil sensitivities
    xVecRes = (-floor(nxRes/2) : ceil(nxRes/2)-1) * spacingRes(1);
    yVecRes = (-floor(nyRes/2) : ceil(nyRes/2)-1) * spacingRes(2);
    zVecRes = (-floor(nzRes/2) : ceil(nzRes/2)-1) * spacingRes(3);
    kxIdxRes = floor(nx/2) - floor(nxRes/2) + (1:nxRes);
    kyIdxRes = floor(ny/2) - floor(nyRes/2) + (1:nyRes);
    kzIdxRes = floor(nz/2) - floor(nzRes/2) + (1:nzRes);
    dataRes = sortedData(kxIdxRes, kyIdxRes, kzIdxRes, :, :);
    samplMaskRes = samplMask(1, kyIdxRes, kzIdxRes, 1, :);
    csmRes = interpn(xVec, yVec, zVec, csm, xVecRes.', yVecRes, zVecRes);
    BxRes = interp1(xVec, Bx, xVecRes);
    ByRes = interp1(yVec, By, yVecRes);
    BzRes = interp1(zVec, Bz, zVecRes);

    % Convert spline coefficients to velocity field
    velocity = reshape(velSplEst, 1, nBx, nBy, nBz, nDims, nTime-1);
    velocity = pagemtimes(velocity, 'none', BxRes, 'transpose');
    velocity = reshape(velocity, nxRes, nBy, nBz, nDims, nTime-1);
    velocity = pagemtimes(velocity, 'none', ByRes, 'transpose');
    velocity = reshape(velocity, nxRes*nyRes, nBz, nDims, nTime-1);
    velocity = pagemtimes(velocity, 'none', BzRes, 'transpose');
    velocity = reshape(velocity, nxRes, nyRes, nzRes, nDims, nTime-1);

    for iIteration = 1:nIterations(iResLevel)
        ticIter = tic;
        fprintf('Resolution level %d/%d, outer iteration %d/%d ...\n', iResLevel, ...
            nResLevels, iIteration, nIterations(iResLevel));
        
        [lhsIm, rhsIm, preconFcn] = createModelsImages(lambda, dataRes, ...
            samplMaskRes, csmRes, velocity, fieldOfView, timeStep);

        % Uncomment for debugging
        % testOperatorFcn(lhsIm, [nxRes, nyRes, nzRes, nTime], ...
        %     [nxRes*nyRes*nzRes*nTime*nVirtCoils + nxRes*nyRes*nzRes*(nTime-1), 1]);

        imEst = lsqrWrapper(lhsIm, rhsIm, 1e-6, 500, preconFcn, [], imEst(:));
        imEst = reshape(imEst, [nxRes, nyRes, nzRes, nTime]);

        % Fix splines in the background to zero
        vEst = velSplEst(fgIdx,:,:);
        
        [lhsVel, rhsVel] = createModelsVelocity(lambda, mu, imEst, ...
            fieldOfView, timeStep, {BxRes,ByRes,BzRes}, fgIdx);

        % Uncomment for debugging
        % testOperatorFcn(lhsVel, [nnz(fgIdx), nDims, nTime-1], ...
        %     [2*(nxRes*nyRes*nzRes*(nTime-1) + prod(nSplines)*nDims*(nTime-1)), 1], false);

        vEst = lsqrWrapper(lhsVel, rhsVel, 1e-6, 1000, [], [], vEst(:));
        vEst = reshape(vEst, nnz(fgIdx), nDims, nTime-1);
        
        velSplEst(fgIdx,:,:) = vEst;

        % Convert spline coefficients to velocity field
        velocity = reshape(velSplEst, 1, nBx, nBy, nBz, nDims, nTime-1);
        velocity = pagemtimes(velocity, 'none', BxRes, 'transpose');
        velocity = reshape(velocity, nxRes, nBy, nBz, nDims, nTime-1);
        velocity = pagemtimes(velocity, 'none', ByRes, 'transpose');
        velocity = reshape(velocity, nxRes*nyRes, nBz, nDims, nTime-1);
        velocity = pagemtimes(velocity, 'none', BzRes, 'transpose');
        velocity = reshape(velocity, nxRes, nyRes, nzRes, nDims, nTime-1);

        tIter = toc(ticIter);
        fprintf('Done in %.2f s\n', tIter)
    end

end

%% Save reconstruction result
if saveResult == true
    subjFolder = sprintf('%s/volunteer%d', reconPath, subjectIdx);
    if ~isfolder(subjFolder)
        mkdir(subjFolder)
    end
    recFile = fullfile(subjFolder, sprintf('dynamic%d_recon.h5', dynScanIdx));
    
    writeRecon(recFile, imEst, timeStep, fieldOfView, velSplEst, Bx, By, Bz);
end

%% Show result
if showResult == true
    % Integrate velocity field to displacement field
    velocity = cat(5, zeros(nx, ny, nz, nDims, 1, 'like', velocity), velocity);
    displ = velocityToDisplacement(velocity, spacing, timeStep, 1);
    
    % Find time frame of maximum deformation
    [~, tMax] = max(vecnorm(reshape(displ, [], size(displ,5))));
    
    % Create mask for displacement field
    maskInit = createMask(abs(imEst(:,:,:,1)), prctile(abs(imEst(:)), 65));
    mask = createMask(abs(imEst(:,:,:,tMax)), prctile(abs(imEst(:)), 65));
    
    % Location of the shown slices
    xIdx = 40;
    yIdx = 39;
    zIdx = 25;
    
    % Create slices
    imSlices = {abs(imEst(:,:,zIdx,tMax)) .* mask(:,:,zIdx), ...
        squeeze(abs(imEst(:,yIdx,:,tMax)) .* mask(:,yIdx,:)), ...
        squeeze(abs(imEst(xIdx,:,:,tMax)) .* mask(xIdx,:,:))};
    displSlices = {squeeze(displ(:,:,zIdx,1:2,tMax) .* maskInit(:,:,zIdx))./reshape(spacing(1:2), 1, 1, 2), ...
        squeeze(displ(:,yIdx,:,[1,3],tMax) .* maskInit(:,yIdx,:))./reshape(spacing([1,3]), 1, 1, 2), ...
        squeeze(displ(xIdx,:,:,2:3,tMax) .* maskInit(xIdx,:,:))./reshape(spacing(2:3), 1, 1, 2)};
    
    % Create figure
    widths = [0.1, 1, 0.1, 1, 0.1, 0.5, 0.1];
    heights = [0.1, 0.5, 0.1];
    axPos = createAxesPositions(widths, heights);
    
    hFig = figure('Name', 'Displacements', 'Color', 'white', 'DefaultAxesFontSize', 14);
    hFig.Position(1:3) = [1, 1, 1200];
    hFig.Position(4) = hFig.Position(3)*sum(heights)/sum(widths);
    
    arrowColor = [255, 204, 0]/255;

    titles = {'Transverse', 'Coronal', 'Sagittal'};
    
    for iView = 1:3
        hAx = axes('Position', axPos{1,iView});
        imshow(imSlices{iView}.')
        clim([0, max(abs(imEst(:)))])
        hold on
        xVec = 3:4:size(imSlices{iView},1);
        yVec = 3:4:size(imSlices{iView},2);
        ux = displSlices{iView}(xVec,yVec,1).';
        uy = displSlices{iView}(xVec,yVec,2).';
        nzIdx = (ux ~= 0) | (uy ~= 0);
        [xGrid,yGrid] = meshgrid(xVec,yVec);
        quiver(hAx, xGrid(nzIdx), yGrid(nzIdx), ux(nzIdx), uy(nzIdx), ...
            'LineWidth', 0.5, 'Color', arrowColor, 'AutoScale', 'off')
        title(titles{iView})
    end
end

%% Local functions
function varargout = lsqrWrapper(varargin)
%lsqrWrapper Wrapper around lsqr to provide useful error messages inside
%the linear operator functions.

    varargout = cell(1, nargout);
    try
        [varargout{:}] = lsqr(varargin{:});
    catch ME
        if strcmp(ME.identifier, 'MATLAB:iterapp:InvalidInput')
            ME = ME.cause{1};
        end
        rethrow(ME)
    end
end
