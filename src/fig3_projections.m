% Creates Figure 3 in the manuscript
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl
clearvars
close all
clc

addpath('./visualization/')

%% Load data
recFile = '../data/volunteer2/dynamic1_recon.h5';

im = h5read(recFile, '/recon/images');
im = abs(im.real + 1i*im.imag);

dt = h5readatt(recFile, '/recon', 'dt');

climsIm = [0, max(abs(im(:)))];

% Remove slice oversampling
im = im(:, :, 8:57, :);

xIdx = 40;
yIdx = 39;
zIdx = 25;
tIdx = 23;

%% Show figure
widths = [0.1, 2, 0.15, 1, 0.3];
heights = [0.1, 1, 0.15, 1, 0.1];
axPos = createAxesPositions(widths, heights);

hFig = figure('Name', 'Projections', 'Color', 'white', 'DefaultAxesFontSize', 14);

hFig.Position(1:3) = [1, 1, 800];
hFig.Position(4) = hFig.Position(3)*sum(heights)/sum(widths);

% Transverse view
hAx = axes('Position', axPos{1,1});
imshow(squeeze(im(:,:,zIdx,tIdx)).')
hold on
yline(yIdx, 'r')
xline(xIdx, 'r')
clim(climsIm)
hl = annotation(hFig, 'line', 'Color', 'w');
hl.Parent = hAx;
hl.X = size(im,1) - 10 + [-5,5];
hl.Y = size(im,2) - 10*[1,1];
hl = annotation(hFig, 'line', 'Color', 'w');
hl.Parent = hAx;
hl.X = size(im,1) - 10*[1,1];
hl.Y = size(im,2) - 10 + [-5,5];
text(size(im,1)-10+[0,0,-7.5,7.5], size(im,2)-10+[-7.5,7.5,0,0], {'A','P','R','L'}, ...
    'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

% Sagittal view
hAx = axes('Position', axPos{1,2});
imshow(squeeze(im(xIdx,:,:,tIdx)).')
hold on
yline(zIdx, 'r')
clim(climsIm)
hl = annotation(hFig, 'line', 'Color', 'w');
hl.Parent = hAx;
hl.X = size(im,2) - 10 + [-5,5];
hl.Y = size(im,3) - 10*[1,1];
hl = annotation(hFig, 'line', 'Color', 'w');
hl.Parent = hAx;
hl.X = size(im,2) - 10*[1,1];
hl.Y = size(im,3) - 10 + [-5,5];
text(size(im,2)-10+[0,0,-7.5,7.5], size(im,3)-10+[-7.5,7.5,0,0], {'H','F','A','P'}, ...
    'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center')

% RL,t-view
hAx = axes('Position', axPos{2,1});
imshow(squeeze(im(:,yIdx,zIdx,:)).')
daspect auto
clim(climsIm)
hold on
yline(tIdx, 'r');
ya = get(hAx, 'YAxis');
set(ya, 'Visible', 'on', 'FontSize', 14)
yticks([1, size(im,4)])
yticklabels({})
set(hAx, 'YAxisLocation', 'right')

%AP,t-view
hAx = axes('Position', axPos{2,2});
imshow(squeeze(im(xIdx,:,zIdx,:)).')
daspect auto
clim(climsIm)
hold on
yline(tIdx, 'r');
ya = get(hAx, 'YAxis');
set(ya, 'Visible', 'on', 'FontSize', 14)
yticks([1, size(im,4)])
yticklabels(compose('%.1f', [0, size(im,4)-1]*dt))
set(hAx, 'YAxisLocation', 'right')
hl = ylabel('Time (s)', 'Rotation', -90);
hl.Position(1) = 77.5;

annotation('arrow', sum(axPos{2,2}([1,3])) + 0.8*widths(3)/sum(widths)/2*[1,1], ...
    axPos{2,2}(2) + [axPos{2,2}(4), 0] + 0.25*axPos{2,2}(4)*[-1,1], 'LineWidth', 1.5)
