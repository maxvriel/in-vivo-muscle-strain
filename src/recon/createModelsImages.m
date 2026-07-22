function [lhsModel, rhsVec, preconFcn] = createModelsImages(lambda, data, samplMask, csm, velocity, fieldOfView, timeStep)
%createModelsImages Create linear operator for the image reconstruction
%step of the joint reconstruction problem.
%
%   Inputs:
%   lambda: Regularization parameter for the motion model
%   data: Acquired k-space data as a 5D array [x, y, z, coil, time]
%   samplMask: Binary sampling mask as a 5D array [1, y, z, 1, time]
%   csm: Coil sensitivity maps as a 4D array [x, y, z, coil]
%   velocity: Velocity field as a 5D array [x, y, z, dim, time]
%   fieldOfView: Field of view in each spatial dimension
%   timeStep: Time step between the time frames
%
%   Outputs:
%   lhsModel: Function handle for the left-hand-side operator of the linear
%       system, with the image as input
%   rhsVec: Right-hand-side vector of the linear system
%   preconFcn: Function handle for the preconditioner
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

% Dimensions
[nx, ny, nz, nCoils, nTime] = size(data, 1:5);
nDims = size(velocity, 4);

% Prevent casting during optimization
samplMask = cast(samplMask, 'like', data);

% Create spatial derivative factor in k-space
Dkx = 1i*2*pi * reshape(-floor(nx/2) : ceil(nx/2)-1, nx, 1) * (1/fieldOfView(1));
Dky = 1i*2*pi * reshape(-floor(ny/2) : ceil(ny/2)-1, 1, ny) * (1/fieldOfView(2));
Dkz = 1i*2*pi * reshape(-floor(nz/2) : ceil(nz/2)-1, 1, 1, nz) * (1/fieldOfView(3));
Dk = cat(4, Dkx.*ones(1,ny,nz), Dky.*ones(nx,1,nz), Dkz.*ones(nx,ny,1));

% Calculate diagonal preconditioner
preconDiag = sum(abs(csm).^2, 4);
tempKernel = zeros(1, 1, 1, 1, nTime, 'like', csm);
tempKernel(floor(nTime/2)+[0,1,2]) = [-1, 2, -1] / timeStep^2;
tempKernel = ifftshift(tempKernel);
% Kernel is symmetric in time, so fft should be real
tempKernel = real(fft(tempKernel, [], 5));
preconDiag = preconDiag + lambda*tempKernel;
nzIdx = preconDiag ~= 0;
preconDiag(nzIdx) = 1 ./ sqrt(preconDiag(nzIdx));

% Precompute scaling factors
csm = csm / sqrt(ny*nz);

% Precompute fftshifts
csm = ifftshift(ifftshift(csm, 2), 3);
data = ifftshift(ifftshift(data, 2), 3);
samplMask = ifftshift(ifftshift(samplMask, 2), 3);
Dk = ifftshift(Dk, 1);
Dk = ifftshift(Dk, 2);
Dk = ifftshift(Dk, 3);

% Precompute IFFT along fully sampled readout direction
data = ifftshift(data, 1);
data = ifft(data, [], 1);
data = data * sqrt(nx);
data = fftshift(data, 1);

% Check if the velocity field contains nonzeros
anyVel = any(velocity(:));

lhsModel = @lhsFcn;
rhsVec = [sqrt(0.5)*data(:); zeros(nx*ny*nz*1*(nTime-1), 1, 'like', data)];
preconFcn = @applyPrecon;

% Use nested functions to define the linear operators and preconditioner

function y = lhsFcn(x, opt)
    % Left-hand-side operator function
    
    if nargin < 2 || strcmp(opt, 'notransp')
        % Pass through model operators
        y = [sqrt(0.5)*dataModel(x, 'notransp'); sqrt(0.5*lambda)*motionModel(x, 'notransp')];
    elseif strcmp(opt, 'transp')
        % Add model operators
        y = sqrt(0.5)*dataModel(x(1:nx*ny*nz*nCoils*nTime), 'transp') + ...
            sqrt(0.5*lambda)*motionModel(x(nx*ny*nz*nCoils*nTime+1:end), 'transp');
    else
        error('Invalid option. Use ''notransp'' or ''transp''.');
    end

end

function y = dataModel(x, opt)
    % Data model function

    if strcmp(opt, 'notransp')
        % Apply the data model to x

        x = reshape(x, nx, ny, nz, 1, nTime);

        x = ifftshift(ifftshift(x, 2), 3);

        % Apply coil sensitivity maps
        y = x .* csm;

        % Forward FFT
        y = fft(y, [], 2);
        y = fft(y, [], 3);

        % Apply sampling mask
        y = y .* samplMask;

        y = y(:);

    elseif strcmp(opt, 'transp')
        % Apply the adjoint of the data model to x

        x = reshape(x, nx, ny, nz, nCoils, nTime);

        % Apply sampling mask
        x = x .* samplMask;

        % Adjoint FFT
        x = conj(x);
        x = fft(x, [], 2);
        x = fft(x, [], 3);
        x = conj(x);

        % Apply coil sensitivity maps
        x = x .* conj(csm);
        y = sum(x, 4);

        y = fftshift(fftshift(y, 2), 3);

        y = y(:);

    else
        error('Invalid option. Use ''notransp'' or ''transp''.');
    end

end

function y = motionModel(x, opt)
    % Motion model function

    if strcmp(opt, 'notransp')
        % Apply the motion model to x

        x = reshape(x, nx, ny, nz, 1, nTime);

        if anyVel
            % Temporal interpolation
            y = (x(:,:,:,:,1:end-1) + x(:,:,:,:,2:end)) * 0.5;

            % Multiply with the velocity field
            y = y .* velocity;

            % Spatial derivative in the frequency domain
            y = fft(y, [], 1);
            y = fft(y, [], 2);
            y = fft(y, [], 3);

            y = y .* Dk;

            % Sum over dimensions
            y = sum(y, 4);

            y = ifft(y, [], 1);
            y = ifft(y, [], 2);
            y = ifft(y, [], 3);
        else
            y = zeros(1, 'like', x);
        end

        % Temporal derivative
        y = y + diff(x, 1, 5) / timeStep;

        y = y(:);

    elseif strcmp(opt, 'transp')
        % Apply the adjoint of the motion model to x
    
        x = reshape(x, nx, ny, nz, 1, nTime-1);

        if anyVel
            y = fft(x, [], 1);
            y = fft(y, [], 2);
            y = fft(y, [], 3);

            % Sum over dimensions
            y = repmat(y, 1, 1, 1, nDims, 1);

            % Spatial derivative in the frequency domain
            y = y .* conj(Dk);

            y = ifft(y, [], 1);
            y = ifft(y, [], 2);
            y = ifft(y, [], 3);

            % Multiply with the velocity field
            % Velocity is real, so conj(velocity) = velocity
            y = y .* velocity;
            y = sum(y, 4);

            % Temporal interpolation
            z = zeros(nx, ny, nz, 1, nTime, 'like', y);
            z(:,:,:,:,1:end-1) = y;
            z(:,:,:,:,2:end) = z(:,:,:,:,2:end) + y;
            z = z * 0.5;
        else
            z = zeros(1, 'like', x);
        end

        % Temporal derivative
        y = zeros(nx, ny, nz, 1, nTime, 'like', x);
        y(:,:,:,:,1:end-1) = -x;
        y(:,:,:,:,2:end) = y(:,:,:,:,2:end) + x;
        y = y / timeStep;

        y = y(:) + z(:);

    else
        error('Invalid option. Use ''notransp'' or ''transp''.');
    end

end

function x = applyPrecon(x, opt)
    % Preconditioner function
    
    if strcmp(opt, 'notransp')

        x = reshape(x, nx, ny, nz, 1, nTime);

        % FFT in time
        x = fft(x, [], 5);

        % Apply diagonal preconditioner
        x = x .* preconDiag;

        % Inverse FFT in time
        x = ifft(x, [], 5);

        x = x(:);
        
    elseif strcmp(opt, 'transp')
        
        x = reshape(x, nx, ny, nz, 1, nTime);

        % FFT in time
        x = fft(x, [], 5);

        % Apply diagonal preconditioner 
        x = x .* conj(preconDiag);

        % Inverse FFT in time
        x = ifft(x, [], 5);

        x = x(:);

    else
        error('Invalid option. Use ''notransp'' or ''transp''.');
    end

end

end
