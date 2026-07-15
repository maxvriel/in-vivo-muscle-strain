function oss = calcOctahedralShearStrain(displ, spacing)
%calcOss Calculate octahedral shear strain from a displacement field
%
%   Inputs:
%   displ: Displacement field as a 5D array [x,y,z,u,t]
%   spacing: Voxel spacing in each dimension in meters
%
%   Outputs:
%   oss: Octahedral shear strain as a 3D array [x,y,z,t]
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    displ (:,:,:,:,:) {mustBeNumeric}
    spacing (1,3) {mustBePositive, mustBeNonempty}
end

oss = zeros(size(displ, [1:3,5]), 'like', displ);
for iTime = 1:size(displ, 5)
    % Calculate deformation gradient tensor
    F = zeros([3, 3, size(displ, 1:3)], 'like', displ);
    for iDim = 1:3
        % Calculate the gradient of the displacement field in each dimension
        u = reshape(displ(:,:,:,iDim,iTime), [1, 1, size(displ, 1:3)]);
        [~, ~, dudx, dudy, dudz] = gradient(u, 1, 1, spacing(1), spacing(2), spacing(3));
        F(iDim,1,:,:,:) = dudx;
        F(iDim,2,:,:,:) = dudy;
        F(iDim,3,:,:,:) = dudz;
        F(iDim,iDim,:,:) = F(iDim,iDim,:,:) + 1; % Add identity to diagonal
    end

    % Calculate the right Cauchy-Green deformation tensor C = F^T * F
    C = pagemtimes(F, 'transpose', F, 'none');

    % Calculate eigenvectors and eigenvalues of C
    % Since C is real and symmetric, its eigenvalues are real
    mu = pageeig(C);
    mu = real(mu);

    % Calculate principal strains of logarithmic strain tensor E = 0.5*ln(C)
    Ep = 0.5*log(mu);

    % Calculate octahedral shear strain
    s = 1/3*sqrt((Ep(1,:)-Ep(2,:)).^2 + (Ep(2,:)-Ep(3,:)).^2 + (Ep(1,:)-Ep(3,:)).^2);
    oss(:,:,:,iTime) = reshape(s, size(displ, 1:3));
end

end
