% Creates Figure 4
% Figure 4: Three-dimensional images and displacement fields at the time
% point of maximum deformation, for all four dynamic scans, in three
% orthogonal planes. The pressure cuff was positioned in the initial
% position for the top two rows, and was rotated for the bottom two rows.
% The muscles were passive for rows 1 and 3, while an isometric contraction
% was performed for rows 2 and 4. Note that the deformation fields differ
% considerably when performing the isometric contraction or when rotating
% the cuff. An animated version of this figure over time is available as
% Supplementary Video SV1. A: Anterior; P: Posterior; R: Right; L: Left; H:
% Head; F: Feet.
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl
clearvars
close all

addpath(genpath(pwd))

%% Load data
recFiles = {'recon/volunteer6/dynamic1_recon.h5', 'recon/volunteer6/dynamic2_recon.h5', ...
    'recon/volunteer6/dynamic3_recon.h5', 'recon/volunteer6/dynamic4_recon.h5'};

images = cell(size(recFiles));
displacements = cell(size(recFiles));
dt = zeros(size(recFiles));
climsIm = [0,0];
for iRec = 1:length(recFiles)
    % Load images
    im = h5read(recFiles{iRec}, '/recon/images');
    im = im.real + 1i*im.imag;

    climsIm(2) = max(climsIm(2), max(abs(im(:))));

    % Load velocity field
    Bx = h5read(recFiles{iRec}, '/recon/velocity/splines_x');
    By = h5read(recFiles{iRec}, '/recon/velocity/splines_y');
    Bz = h5read(recFiles{iRec}, '/recon/velocity/splines_z');
    coeffs = h5read(recFiles{iRec}, '/recon/velocity/coefficients');
    dt = h5readatt(recFiles{iRec}, '/recon', 'dt');
    FOV = h5readatt(recFiles{iRec}, '/recon', 'FOV');
    spacing = FOV ./ size(im, 1:3)';

    nDims = 3;
    [nx, nBx] = size(Bx);
    [ny, nBy] = size(By);
    [nz, nBz] = size(Bz);
    nt = size(coeffs, 2);
    
    vel = reshape(coeffs, 1, nBx, nBy, nBz, nDims, nt);
    vel = pagemtimes(vel, 'none', Bx, 'transpose');
    vel = reshape(vel, nx, nBy, nBz, nDims, nt);
    vel = pagemtimes(vel, 'none', By, 'transpose');
    vel = reshape(vel, nx*ny, nBz, nDims, nt);
    vel = pagemtimes(vel, 'none', Bz, 'transpose');
    vel = reshape(vel, nx, ny, nz, nDims, nt);

    % Integrate velocity field to displacement field
    displ = velocityToDisplacement(vel, spacing, dt, 1);

    % Find time frame of maximum deformation
    [~, tMax] = max(vecnorm(reshape(displ, [], size(displ,5))));
    
    % Create mask for displacement field
    mask = createMask(abs(im(:,:,:,1)), prctile(abs(im(:)), 65));

    % Remove slice oversampling
    images{iRec} = abs(im(:,:,8:57,tMax));
    displacements{iRec} = displ(:,:,8:57,:,tMax) .* mask(:,:,8:57);
end

% Location of the shown slices
xIdx = 40;
yIdx = 39;
zIdx = 25;

%% Show figure
widths = [0.2, 0.1, 0.05, 1, 0.1, 1, 0.1, 0.5, 0.1];
heights = [0.1, repmat([0.5, 0.1], 1, 4)];
axPos = createAxesPositions(widths, heights);

hFig = figure('Name', 'Displacements', 'Color', 'white', 'DefaultAxesFontSize', 14);
hFig.Position(1:3) = [1, 1, 1200];
hFig.Position(4) = hFig.Position(3)*sum(heights)/sum(widths);

arrowColor = [255, 204, 0]/255;

ylabs = {{'Passive muscles', '(Dynamic scan 1)'}, {'Isometric contractions', '(Dynamic scan 2)'}, ...
    {'Passive muscles', '(Dynamic scan 3)'}, {'Isometric contractions', '(Dynamic scan 4)'}};

for iRec = 1:length(recFiles)
    % Transverse view
    hAx = axes('Position', axPos{iRec,2});
    imshow(images{iRec}(:,:,zIdx).')
    clim(climsIm)
    hold on
    xVec = 3:4:size(images{iRec},1);
    yVec = 3:4:size(images{iRec},2);
    ux = displacements{iRec}(xVec,yVec,zIdx,1).' / spacing(1);
    uy = displacements{iRec}(xVec,yVec,zIdx,2).' / spacing(2);
    nzIdx = (ux ~= 0) | (uy ~= 0);
    [xGrid,yGrid] = meshgrid(xVec,yVec);
    quiver(hAx, xGrid(nzIdx), yGrid(nzIdx), ux(nzIdx), uy(nzIdx), ...
        'LineWidth', 0.5, 'Color', arrowColor, 'AutoScale', 'off')
    ylabel(ylabs{iRec})
    if iRec == 1
        hl = annotation(hFig, 'line', 'Color', 'w', 'X', size(images{iRec},1) - 10 + [-5,5], 'Y', size(images{iRec},2) - 10*[1,1]);
        hl.Parent = hAx;
        hl = annotation(hFig, 'line', 'Color', 'w', 'X', size(images{iRec},1) - 10*[1,1], 'Y', size(images{iRec},2) - 10 + [-5,5]);
        hl.Parent = hAx;
        text(size(images{iRec},1)-10+[0,0,-7.5,7.5], size(images{iRec},2)-10+[-7.5,7.5,0,0], {'A','P','R','L'}, ...
            'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

        axes('Position', axPos{iRec,2}, 'Visible', 'off');
        title('Transverse', 'Visible', 'on', 'FontSize', 14)
    end

    % Coronal view
    hAx = axes('Position', axPos{iRec,3});
    imshow(squeeze(images{iRec}(:,yIdx,:)).')
    clim(climsIm)
    hold on
    xVec = 3:4:size(images{iRec},1);
    zVec = 3:4:size(images{iRec},3);
    ux = squeeze(displacements{iRec}(xVec,yIdx,zVec,1)).' / spacing(1);
    uz = squeeze(displacements{iRec}(xVec,yIdx,zVec,3)).' / spacing(3);
    nzIdx = (ux ~= 0) | (uz ~= 0);
    [xGrid,zGrid] = meshgrid(xVec,zVec);
    quiver(hAx, xGrid(nzIdx), zGrid(nzIdx), ux(nzIdx), uz(nzIdx), ...
        'LineWidth', 0.5, 'Color', arrowColor, 'AutoScale', 'off')
    if iRec == 1
        hl = annotation(hFig, 'line', 'Color', 'w', 'X', size(images{iRec},1) - 10 + [-5,5], 'Y', size(images{iRec},3) - 10*[1,1]);
        hl.Parent = hAx;
        hl = annotation(hFig, 'line', 'Color', 'w', 'X', size(images{iRec},1) - 10*[1,1], 'Y', size(images{iRec},3) - 10 + [-5,5]);
        hl.Parent = hAx;
        text(size(images{iRec},1)-10+[0,0,-7.5,7.5], size(images{iRec},3)-10+[-7.5,7.5,0,0], {'H','F','R','L'}, ...
            'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

        axes('Position', axPos{iRec,3}, 'Visible', 'off');
        title('Coronal', 'Visible', 'on', 'FontSize', 14)
    end

    % Sagittal view
    hAx = axes('Position', axPos{iRec,4});
    imshow(squeeze(images{iRec}(xIdx,:,:)).')
    clim(climsIm)
    hold on
    yVec = 3:4:size(images{iRec},2);
    zVec = 3:4:size(images{iRec},3);
    uy = squeeze(displacements{iRec}(xIdx,yVec,zVec,2)).' / spacing(2);
    uz = squeeze(displacements{iRec}(xIdx,yVec,zVec,3)).' / spacing(3);
    nzIdx = (uy ~= 0) | (uz ~= 0);
    [yGrid,zGrid] = meshgrid(yVec,zVec);
    quiver(hAx, yGrid(nzIdx), zGrid(nzIdx), uy(nzIdx), uz(nzIdx), ...
        'LineWidth', 0.5, 'Color', arrowColor, 'AutoScale', 'off')
    if iRec == 1
        hl = annotation(hFig, 'line', 'Color', 'w', 'X', size(images{iRec},2) - 10 + [-5,5], 'Y', size(images{iRec},3) - 10*[1,1]);
        hl.Parent = hAx;
        hl = annotation(hFig, 'line', 'Color', 'w', 'X', size(images{iRec},2) - 10*[1,1], 'Y', size(images{iRec},3) - 10 + [-5,5]);
        hl.Parent = hAx;
        text(size(images{iRec},2)-10+[0,0,-7.5,7.5], size(images{iRec},3)-10+[-7.5,7.5,0,0], {'H','F','A','P'}, ...
            'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

        axes('Position', axPos{iRec,4}, 'Visible', 'off');
        title('Sagittal', 'Visible', 'on', 'FontSize', 14)
    end
end

axes('Position', combineAxesPositions(axPos(1:2,:)), 'Visible', 'off')
ylabel('Cuff position 1', 'Visible', 'on')
annotation('line', 0.04*[1,1], [sum(axPos{1,1}([2,4])), axPos{2,1}(2)])
annotation('line', [0.04,0.05], axPos{2,1}(2)*[1,1])
annotation('line', [0.04,0.05], sum(axPos{1,1}([2,4]))*[1,1])
axes('Position', combineAxesPositions(axPos(3:4,:)), 'Visible', 'off')
ylabel('Cuff position 2', 'Visible', 'on')
annotation('line', 0.04*[1,1], [sum(axPos{3,1}([2,4])), axPos{4,1}(2)])
annotation('line', [0.04,0.05], axPos{4,1}(2)*[1,1])
annotation('line', [0.04,0.05], sum(axPos{3,1}([2,4]))*[1,1])

exportgraphics(hFig, './figures/figure4.eps', 'ContentType', 'vector', 'Padding', 'figure')
exportgraphics(hFig, './figures/figure4.png', 'Resolution', 300, 'Padding', 'figure')
