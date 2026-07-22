function [lhsModel, rhsVec] = createModelsVelocity(lambda, mu, imEst, fieldOfView, timeStep, splineBasis, fgIdx)
%createModelsVelocity Create linear operator for the velocity estimation
%step of the joint reconstruction problem.
%
%   Inputs:
%   lambda: Regularization parameter for the motion model
%   mu: Regularization parameter for the spline coefficients regularization
%   imEst: Estimated images as a 4D array [x, y, z, time]
%   fieldOfView: Field of view in each spatial dimension
%   timeStep: Time step between the time frames
%   splineBasis: Cell array containing the B-spline basis functions for
%       each spatial dimension
%   fgIdx: Boolean mask indicating the splines that are not in the
%       background
%
%   Outputs:
%   lhsModel: Function handle for the left-hand-side operator of the linear
%       system, with the velocity spline coefficients as input
%   rhsVec: Right-hand-side vector of the linear system
%
% Copyright (c) 2026, UMC Utrecht 
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

% Dimensions
[nx, ny, nz, nTime] = size(imEst);
nBx = size(splineBasis{1}, 2);
nBy = size(splineBasis{2}, 2);
nBz = size(splineBasis{3}, 2);
nDims = 3;

% Create spatial derivative factor in k-space
Dkx = 1i*2*pi * reshape(-floor(nx/2) : ceil(nx/2)-1, nx, 1) * (1/fieldOfView(1));
Dky = 1i*2*pi * reshape(-floor(ny/2) : ceil(ny/2)-1, 1, ny) * (1/fieldOfView(2));
Dkz = 1i*2*pi * reshape(-floor(nz/2) : ceil(nz/2)-1, 1, 1, nz) * (1/fieldOfView(3));
Dk = cat(4, Dkx.*ones(1,ny,nz), Dky.*ones(nx,1,nz), Dkz.*ones(nx,ny,1));

% Temporal derivative
dmdt = diff(imEst, 1, 4) / timeStep;

% Temporal interpolation
mTempInterp = (imEst(:,:,:,1:end-1) + imEst(:,:,:,2:end)) * 0.5;
mTempInterp = reshape(mTempInterp, nx, ny, nz, 1, nTime-1);
mTempInterp = repmat(mTempInterp, 1, 1, 1, nDims, 1);

Bx = splineBasis{1};
By = splineBasis{2};
Bz = splineBasis{3};

% Precompute fftshifts
Dk = ifftshift(Dk, 1);
Dk = ifftshift(Dk, 2);
Dk = ifftshift(Dk, 3);

% Precompute FFT
mTempInterp = mTempInterp / sqrt(nx*ny*nz);
dmdt = fft(dmdt, [], 1);
dmdt = fft(dmdt, [], 2);
dmdt = fft(dmdt, [], 3);
dmdt = dmdt / sqrt(nx*ny*nz);

% Check if the images contains nonzeros
anyIm = any(imEst(:));

lhsModel = @lhsFcn;
rhsVec = [sqrt(0.5*lambda)*-dmdt(:); zeros(nBx*nBy*nBz*nDims*(nTime-1), 1, 'like', dmdt)];

% Force the solution to be real by splitting the real and imaginary rows
rhsVec = [real(rhsVec); imag(rhsVec)];

% Use nested functions to define the linear operators

function y = lhsFcn(x, opt)
    % Left-hand-side operator function
    
    if nargin < 2 || strcmp(opt, 'notransp')
        % Fill in the background spline coefficients with zeros
        xFull = zeros(nBx*nBy*nBz, nDims, nTime-1, 'like', x);
        xFull(fgIdx,:,:) = reshape(x, [], nDims, nTime-1);
        % Pass through model operators
        y = [sqrt(0.5*lambda)*motionModel(xFull, 'notransp'); sqrt(0.5*mu)*xFull(:)];
        % Force the solution to be real by splitting the real and imaginary rows
        y = [real(y); imag(y)];
    elseif strcmp(opt, 'transp')
        % Combine real and imaginary rows
        x = x(1:end/2) + 1i*x(end/2+1:end);
        % Add model operators
        y = real(sqrt(0.5*lambda)*motionModel(x(1:nx*ny*nz*(nTime-1)), 'transp')) + ...
            real(sqrt(0.5*mu)*x(nx*ny*nz*(nTime-1)+1:end));
        % Remove the background spline coefficients
        y = reshape(y, nBx*nBy*nBz, nDims, nTime-1);
        y = y(fgIdx,:,:);
        y = y(:);
    else
        error('Invalid option. Use ''notransp'' or ''transp''.');
    end

end

function y = motionModel(x, opt)
    % Motion model function

    if strcmp(opt, 'notransp')
        % Apply the motion model to x

        if anyIm
            y = reshape(x, nBx, nBy, nBz, nDims, nTime-1);

            % Convert spline coefficients to velocity field
            y = reshape(y, 1, nBx, nBy, nBz, nDims, nTime-1);
            y = pagemtimes(y, 'none', Bx, 'transpose');
            y = reshape(y, nx, nBy, nBz, nDims, nTime-1);
            y = pagemtimes(y, 'none', By, 'transpose');
            y = reshape(y, nx*ny, nBz, nDims, nTime-1);
            y = pagemtimes(y, 'none', Bz, 'transpose');
            y = reshape(y, nx, ny, nz, nDims, nTime-1);

            % Multiply with the images
            y = y .* mTempInterp;

            % Spatial derivative in the frequency domain
            y = fft(y, [], 1);
            y = fft(y, [], 2);
            y = fft(y, [], 3);

            y = y .* Dk;

            % Sum over dimensions
            y = sum(y, 4);
        else
            y = zeros(nx, ny, nz, 1, nTime-1, 'like', x);
        end

        y = y(:);

    elseif strcmp(opt, 'transp')
        % Apply the adjoint of the motion model to x
    
        if anyIm
            x = reshape(x, nx, ny, nz, 1, nTime-1);

            % Sum over dimensions
            y = repmat(x, 1, 1, 1, nDims, 1);

            % Spatial derivative in the frequency domain
            y = y .* conj(Dk);

            y = conj(y);
            y = fft(y, [], 1);
            y = fft(y, [], 2);
            y = fft(y, [], 3);
            y = conj(y);

            % Multiply with the images
            y = y .* conj(mTempInterp);

            % Convert spline coefficients to velocity field
            y = reshape(y, nx*ny, nz, nDims, nTime-1);
            y = pagemtimes(y, 'none', conj(Bz), 'none');
            y = reshape(y, nx, ny, nBz, nDims, nTime-1);
            y = pagemtimes(y, 'none', conj(By), 'none');
            y = reshape(y, 1, nx, nBy, nBz, nDims, nTime-1);
            y = pagemtimes(y, 'none', conj(Bx), 'none');
            y = reshape(y, nBx, nBy, nBz, nDims, nTime-1);
        else
            y = zeros(nBx, nBy, nBz, nDims, nTime-1, 'like', x);
        end

        y = y(:);

    else
        error('Invalid option. Use ''notransp'' or ''transp''.');
    end

end

end
