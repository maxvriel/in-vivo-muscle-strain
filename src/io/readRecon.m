function [images, timeStep, fieldOfView, velocity] = readRecon(fileName)
%readRecon Load reconstructed images and velocity fields from HDF5 file.
%
%   Inputs:
%   fileName: Path to the HDF5 file containing the reconstruction data
%
%   Outputs:
%   images: Reconstructed complex-valued images as a 4D array [x, y, z, t]
%   timeStep: Time step between image frames in seconds
%   fieldOfView: Field of view in mm
%   velocity: Velocity fields as a 5D array [x, y, z, v, t]
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    fileName {mustBeTextScalar, mustBeFile}
end

% Load images
images = h5read(fileName, '/recon/images');
images = images.real + 1i*images.imag;

timeStep = h5readatt(fileName, '/recon', 'time_step');
fieldOfView = h5readatt(fileName, '/recon', 'field_of_view');

if nargout >= 4
    % Load velocity field
    Bx = h5read(fileName, '/recon/velocity/splines_x');
    By = h5read(fileName, '/recon/velocity/splines_y');
    Bz = h5read(fileName, '/recon/velocity/splines_z');
    coeffs = h5read(fileName, '/recon/velocity/coefficients');

    nDims = 3;
    [nx, nBx] = size(Bx);
    [ny, nBy] = size(By);
    [nz, nBz] = size(Bz);
    nt = size(coeffs, 2);

    velocity = reshape(coeffs, 1, nBx, nBy, nBz, nDims, nt);
    velocity = pagemtimes(velocity, 'none', Bx, 'transpose');
    velocity = reshape(velocity, nx, nBy, nBz, nDims, nt);
    velocity = pagemtimes(velocity, 'none', By, 'transpose');
    velocity = reshape(velocity, nx*ny, nBz, nDims, nt);
    velocity = pagemtimes(velocity, 'none', Bz, 'transpose');
    velocity = reshape(velocity, nx, ny, nz, nDims, nt);
end

end
