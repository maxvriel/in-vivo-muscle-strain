function uField = velocityToDisplacement(velField, spacing, timeStep, nSteps)
%velocityToDisplacement Convert a velocity field to a displacement field.
%
%   Inputs:
%   velField: Velocity field as a 5D array [x,y,z,v,t]
%   spacing: Spatial spacing of the velocity field in each spatial 
%       dimension
%   timeStep: Time step between the velocity field frames in seconds
%   nSteps (optional): Number of steps for Eulerian integration; more steps
%       is more accurate but slower
%
%   Outputs:
%   uField: Displacement field as a 5D array [x,y,z,u,t]
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    velField (:,:,:,:,:) {mustBeNumeric}
    spacing (1,3) {mustBePositive, mustBeNonempty}
    timeStep (1,1) {mustBePositive}
    nSteps (1,1) {mustBePositive, mustBeInteger} = 1
end

% Coordinate grids for interpolation
xVec = cast((0:size(velField,1)-1) * spacing(1), 'like', velField);
yVec = cast((0:size(velField,2)-1) * spacing(2), 'like', velField);
zVec = cast((0:size(velField,3)-1) * spacing(3), 'like', velField);
[xGrid, yGrid, zGrid] = ndgrid(xVec, yVec, zVec);
x0 = cat(4, xGrid, yGrid, zGrid);

% Allocate output displacement field
uField = zeros(size(velField), 'like', velField);

xt = x0;
for iTime = 1:size(velField, 5)
    vt = velField(:,:,:,:,iTime);
    vInterp = griddedInterpolant({xVec,yVec,zVec}, vt, 'linear', 'none');
    % Use forward Eulerian integration steps
    % Optionally, use smaller time steps
    for iStep = 1:nSteps
        du = vInterp(xt(:,:,:,1), xt(:,:,:,2), xt(:,:,:,3));
        du(isnan(du)) = 0;
        xt = xt + du*(timeStep/nSteps);
    end
    % Store the displacement for the current time step
    uField(:,:,:,:,iTime) = xt - x0;
end

end
