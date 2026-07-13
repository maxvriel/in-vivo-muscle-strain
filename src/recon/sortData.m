function kspace = sortData(data, header, densityComp, dynIdx, userIdx)
%sortData Sort k-space data in a 4D array and averages the data that has
% been acquired multiple times
%
% Inputs:
%   data: ISMRMRD data object
%   header: ISMRMRD header structure
%   densityComp (optional): Boolean indicating whether to apply density
%       compensation, defaults to true
%   dynIdx (optional): Dynamic indices to select, defaults to all
%   userIdx (optional): User-defined indices to select, defaults to all
%
% Outputs:
%   kspace: k-space array [kx, ky, kz, coils, dynamics]
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    data (1,1) ismrmrd.Acquisition
    header (1,1) struct
    densityComp (1,1) logical = true
    dynIdx (1,:) {mustBeInteger, mustBeNonnegative} = ...
        0:header.encoding(1).encodingLimits.repetition.maximum
    userIdx (1,:) {mustBeInteger, mustBeNonnegative} = ...
        0:max(data.head.idx.user(1,:))
end

% Extract header information
nKx = header.encoding(1).encodedSpace.matrixSize.x;
nKy = header.encoding(1).encodedSpace.matrixSize.y;
nKz = header.encoding(1).encodedSpace.matrixSize.z;
kxLim = header.encoding(1).encodingLimits.kspace_encoding_step_0;
kyLim = header.encoding(1).encodingLimits.kspace_encoding_step_1;
kzLim = header.encoding(1).encodingLimits.kspace_encoding_step_2;
nDyn = length(dynIdx);
nCh = header.acquisitionSystemInformation.receiverChannels;

% Allocate k-space array
kspace = zeros(nKx, nKy, nKz, nCh, nDyn, 'like', data.data{1});

for iDyn = 1:nDyn
    % Keep track of the number of times the line has been sampled
    % For CASPR data, this is not a constant over k-space
    count = zeros(1, nKy, nKz);

    % Select k-space readouts
    % Only select the first slice, contrast, etc. in case there are multiple
    acqIdx = find(~data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT') & ...
        data.head.idx.repetition == dynIdx(iDyn) & ...
        data.head.idx.slice == 0 & ...
        data.head.idx.contrast == 0 & ...
        data.head.idx.phase == 0 & ...
        data.head.idx.set == 0 & ...
        data.head.idx.segment == 0 & ...
        ismember(data.head.idx.user(1,:), userIdx));
    for iAcq = 1:length(acqIdx)    
        kx = floor(nKx/2) - kxLim.center + (kxLim.minimum:kxLim.maximum) + 1;
        ky = floor(nKy/2) - kyLim.center + data.head.idx.kspace_encode_step_1(acqIdx(iAcq)) + 1;
        kz = floor(nKz/2) - kzLim.center + data.head.idx.kspace_encode_step_2(acqIdx(iAcq)) + 1;
        kspace(kx,ky,kz,:,iDyn) = kspace(kx,ky,kz,:,iDyn) + ...
            reshape(data.data{acqIdx(iAcq)}, data.head.number_of_samples(acqIdx(iAcq)), 1, 1, nCh);
        count(1,ky,kz) = count(1,ky,kz) + 1;
    end

    % Sampling density compensation, preventing division by zero
    if densityComp
        kspace(:,:,:,:,iDyn) = kspace(:,:,:,:,iDyn) ./ max(count, 1);
    end
end

end
