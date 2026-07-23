% Creates Figure 6
% Figure 6: Validation of the motion fields. Each frame of the dynamic
% image series is warped back to the undeformed state using the estimated
% displacement fields u. This image is then compared to the image of the
% first time frame using the SSIM. The graph shows the SSIM values for all
% time points across all eight volunteers.
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl
clearvars
close all

% Set paths
setup

%% Calculate SSIM values
subjects = [1,2,4:9];
ssimVals = zeros(50, length(subjects));
ssimRef = zeros(50, length(subjects));
plotSubject = 1;
for iSubj = 1:length(subjects)
    recFile = sprintf('recon/volunteer%d/dynamic1_recon.h5', subjects(iSubj));

    % Load reconstruction
    [dynIm, dt, FOV, vel] = readRecon(recFile);
    climsIm = [0, max(abs(dynIm(:)))];

    % Integrate velocity field to displacement field
    spacing = FOV ./ size(dynIm, 1:3)';
    displ = velocityToDisplacement(vel, spacing, dt, 1);

    % Find time frame of maximum deformation
    [~, tMax] = max(vecnorm(reshape(displ, [], size(displ,5))));

    % Warp each frame back to the first frame using the displacement fields
    warpedIm = warpImage(dynIm, spacing, displ, spacing, 'spline', 0);

    % Remove outer slices, select right leg, and take magnitude
    dynIm = abs(dynIm((1:64)+4,:,8:57,:));
    displ = displ((1:64)+4,:,8:57,:,:);
    warpedIm = abs(warpedIm((1:64)+4,:,8:57,:));

    % Calculate SSIM values
    ssimVals(:,iSubj) = squeeze(ssim(warpedIm, repmat(dynIm(:,:,:,1), 1, 1, 1, size(warpedIm, 4)), ...
        'DynamicRange', climsIm(2), 'DataFormat', 'SSSB'));
    ssimRef(:,iSubj) = squeeze(ssim(dynIm, repmat(dynIm(:,:,:,1), 1, 1, 1, size(warpedIm, 4)), ...
        'DynamicRange', climsIm(2), 'DataFormat', 'SSSB'));

    % Save images for plotting
    if subjects(iSubj) == plotSubject
        plotImInit = dynIm(:,:,26,1).';  
        plotImDyn = dynIm(:,:,26,tMax).';
        plotImWarp = warpedIm(:,:,26,tMax).';
        mask = createMask(dynIm(:,:,:,1), prctile(abs(dynIm(:)), 65));
        plotDispl = displ(:,:,26,:,tMax) .* mask(:,:,26) ./ reshape(spacing, 1, 1, 1, 3);
    end
end

%% Show figure
valColor = [252, 96, 57]/255;
initColor = [44, 191, 44]/255;
dyColor = [17, 145, 250]/255;
arrowColor = [255, 204, 0]/255;

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

hFig = figure('Color', 'white', 'DefaultAxesFontSize', 12);
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
xline(0.05, 'Color', initColor, 'LineWidth', 2, 'Alpha', 1)
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

hAx = axes('Position', axPos{1,2});
imshow(plotImDyn, []);
hold on
xVec = 3:4:size(plotImDyn,1);
yVec = 3:4:size(plotImDyn,2);
ux = plotDispl(xVec,yVec,1).';
uy = plotDispl(xVec,yVec,2).';
nzIdx = (ux ~= 0) | (uy ~= 0);
[xGrid,yGrid] = meshgrid(xVec,yVec);
quiver(hAx, xGrid(nzIdx), yGrid(nzIdx), ux(nzIdx), uy(nzIdx), ...
    'LineWidth', 1, 'Color', arrowColor, 'AutoScale', 'off')
annotation('rectangle', tightPosition(hAx), 'Color', dyColor, 'LineWidth', 4);

hAx = axes('Position', axPos{2,1});
imshow(plotImInit, []);
fig1Pos = tightPosition(hAx);
annotation('rectangle', fig1Pos, 'Color', initColor, 'LineWidth', 4);

hAx = axes('Position', axPos{2,2});
imshow(plotImWarp, []);
fig2Pos = tightPosition(hAx);

annotation('arrow', (axPos{1,2}(1)+0.5*axPos{1,2}(3))*[1,1], ...
    [axPos{1,2}(2), axPos{2,2}(2)+axPos{2,2}(4)]+0.015*[-1,1], ...
    'LineWidth', 1.5)
annotation('textbox', [axPos{1,2}(1)+0.5*axPos{1,2}(3), axPos{1,2}(2)-0.5*heights(3)/sum(heights)-0.05/2, 0.2, 0.05], ...
    'String', 'Warping', 'FontSize', 12, 'VerticalAlignment', 'middle', 'EdgeColor', 'none')

annotation('doublearrow', [sum(fig1Pos([1,3])), fig2Pos(1)]+0.015*[1,-1], ...
    fig1Pos(2)+fig1Pos(4)*0.5*[1,1], 'LineWidth', 1.5)
annotation('textbox', [sum([fig1Pos([1,3])-0.1/2, fig2Pos(1)])/2, fig1Pos(2)+fig1Pos(4)/2, 0.1, 0.05], ...
    'String', 'SSIM', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'EdgeColor', 'none')

hAx = axes('Position', combineAxesPositions(axPos(3,:)));
hl = plot((0:size(ssimVals,1)-1)*dt, ssimVals, 'LineWidth', 2);
set(hl, {'Color'}, num2cell(cmap(1:length(subjects),:), 2))
hold on
hl2 = plot((0:size(ssimRef,1)-1)*dt, ssimRef, '--', 'LineWidth', 2);
set(hl2, {'Color'}, num2cell(cmap(1:length(subjects),:), 2))
xlim tight
ylim([0.5, 1])
grid on
xlabel('Time (s)', 'FontSize', 12)
ylabel('SSIM', 'FontSize', 12)

hlVal = plot(NaN, NaN, 'k', 'LineWidth', 2);
hlRef = plot(NaN, NaN, 'k--', 'LineWidth', 2);
hlBlank = plot(NaN, NaN, 'w', 'LineWidth', 2);
legend([hl; hlBlank; hlVal; hlRef], ["Volunteer " + (1:length(subjects)), "", "Validation", "Static reference"], ...
    'Location', 'eastoutside')

axes('Position', combineAxesPositions(axPos(1,:)), 'Visible', 'off');
ht = title('Motion field validation', 'FontSize', 14, 'Visible', 'on');
ht.Position(2) = 1.2;
