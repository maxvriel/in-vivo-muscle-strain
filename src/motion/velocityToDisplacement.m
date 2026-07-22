function displ = velocityToDisplacement(velocity, spacing, timeStep, nSteps)
%velocityToDisplacement Convert a velocity field to a displacement field.
%
%   Inputs:
%   velocity: Velocity field as a 5D array [x,y,z,v,t]
%   spacing: Spatial spacing of the velocity field in each spatial 
%       dimension
%   timeStep: Time step between the velocity field frames
%   nSteps (optional): Number of steps for Eulerian integration; more steps
%       is more accurate but slower
%
%   Outputs:
%   displ: Displacement field as a 5D array [x,y,z,u,t]
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    velocity (:,:,:,:,:) {mustBeNumeric}
    spacing (1,3) {mustBePositive, mustBeNonempty}
    timeStep (1,1) {mustBePositive}
    nSteps (1,1) {mustBePositive, mustBeInteger} = 1
end

% Coordinate grids for interpolation
xVec = cast((0:size(velocity,1)-1) * spacing(1), 'like', velocity);
yVec = cast((0:size(velocity,2)-1) * spacing(2), 'like', velocity);
zVec = cast((0:size(velocity,3)-1) * spacing(3), 'like', velocity);
[xGrid, yGrid, zGrid] = ndgrid(xVec, yVec, zVec);
x0 = cat(4, xGrid, yGrid, zGrid);

% Allocate output displacement field
displ = zeros(size(velocity), 'like', velocity);

xt = x0;
for iTime = 1:size(velocity, 5)
    vt = velocity(:,:,:,:,iTime);
    vInterp = griddedInterpolant({xVec,yVec,zVec}, vt, 'linear', 'none');
    % Use forward Eulerian integration steps
    % Optionally, use smaller time steps
    for iStep = 1:nSteps
        du = vInterp(xt(:,:,:,1), xt(:,:,:,2), xt(:,:,:,3));
        du(isnan(du)) = 0;
        xt = xt + du*(timeStep/nSteps);
    end
    % Store the displacement for the current time step
    displ(:,:,:,:,iTime) = xt - x0;
end

end
