function images = imageRecon(data, header, coilMaps, dynIdx, userIdx)
%imageRecon Reconstruct k-space data to images by performing a
%coil-weighted zero-filled reconstruction.
%
% Inputs:
%   data: ISMRMRD data object
%   header: ISMRMRD header structure
%   coilMaps: Coil sensitivity maps
%   dynIdx: Dynamic indices to select, defaults to all
%   userIdx: User-defined indices to select, defaults to all
%
% Outputs:
%   images: Reconstructed images [x, y, z, dynamics]
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    data (1,1) ismrmrd.Acquisition
    header (1,1) struct
    coilMaps (:,:,:,:) {mustBeNumeric}
    dynIdx (1,:) {mustBeInteger, mustBeNonnegative} = ...
        0:header.encoding(1).encodingLimits.repetition.maximum
    userIdx (1,:) {mustBeInteger, mustBeNonnegative} = ...
        0:max(data.head.idx.user(1,:))
end

% Extract header information
nKx = header.encoding(1).encodedSpace.matrixSize.x;
nKy = header.encoding(1).encodedSpace.matrixSize.y;
nKz = header.encoding(1).encodedSpace.matrixSize.z;
nX = header.encoding(1).reconSpace.matrixSize.x;
nY = header.encoding(1).reconSpace.matrixSize.y;
nZ = header.encoding(1).reconSpace.matrixSize.z;
nDyn = length(dynIdx);
nCh = header.acquisitionSystemInformation.receiverChannels;

% Reconstruct images
% Only select the first slice, contrast, etc. in case there are multiple
images = zeros(nX, nY, nZ, nDyn, 'like', data.data{1});
for iDyn = 1:nDyn
    
    kspace = sortData(data, header, true, dynIdx(iDyn), userIdx);

    % Reconstruct image
    kspace = fftshift(ifft(ifftshift(kspace, 1), [], 1), 1) * sqrt(nKx);
    kspace = fftshift(ifft(ifftshift(kspace, 2), [], 2), 2) * sqrt(nKy);
    kspace = fftshift(ifft(ifftshift(kspace, 3), [], 3), 3) * sqrt(nKz);

    % Apply coil-weighted combination
    im = zeros(nKx, nKy, nKz, 'like', kspace);
    for iCh = 1:nCh
        im = im + conj(coilMaps(:,:,:,iCh)) .* kspace(:,:,:,iCh);
    end
    csmSoS = sum(abs(coilMaps).^2, 4);
    nzIdx = csmSoS ~= 0;
    im(nzIdx) = im(nzIdx) ./ csmSoS(nzIdx);

    % Crop or zero-pad and place in images array
    kIdx = {':', ':', ':'};
    imIdx = {':', ':', ':'};
    if nKx > nX
        kIdx{1} = floor((nKx-nX)/2) + (1:nX);
    elseif nX > nKx
        imIdx{1} = floor((nX-nKx)/2) + (1:nKx);
    end
    if nKy > nY
        kIdx{2} = floor((nKy-nY)/2) + (1:nY);
    elseif nY > nKy
        imIdx{2} = floor((nY-nKy)/2) + (1:nKy);
    end
    if nKz > nZ
        kIdx{3} = floor((nKz-nZ)/2) + (1:nZ);
    elseif nZ > nKz
        imIdx{3} = floor((nZ-nKz)/2) + (1:nKz);
    end
    images(imIdx{:}, iDyn) = im(kIdx{:});

end

end
