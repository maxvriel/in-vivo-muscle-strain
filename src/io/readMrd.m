function [data, header, coilMaps] = readMrd(dataFile, decorrelate)
%readMrd Load ISMRMRD k-space data and coil sensitivity maps from an
%ISMRMRD data file.
%
% Inputs:
%   dataFile: Path to the ISMRMRD data file
%   decorrelate (optional): Boolean flag for coil decorrelation, defaults
%       to true
%
% Outputs:
%   data: ISMRMRD data, as an ismrmrd.Acquisition object
%   header: ISMRMRD header structure
%   coilMaps: Coil sensitivity maps
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    dataFile {mustBeTextScalar, mustBeFile}
    decorrelate (1,1) logical = true
end

% Read ISMRMRD data and header
dataset = ismrmrd.Dataset(dataFile);
header = ismrmrd.xml.deserialize(dataset.readxml());
data = dataset.readAcquisition();
dataset.close();

% Read coil sensitivity maps
coilMaps = h5read(dataFile, '/csm');
coilMaps = coilMaps.real + 1i*coilMaps.imag;

% Apply coil decorrelation
if decorrelate
    % Extract noise data and calculate decorrelation matrix
    noiseIdx = data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
    if nnz(noiseIdx) == 1
        noiseData = data.data{noiseIdx};
        noiseData = noiseData - mean(noiseData);
        noiseCov = conj(noiseData'*noiseData) / (size(noiseData, 1)-1);
        decorrMat = conj(chol(noiseCov));
    elseif nnz(noiseIdx) == 0
        warning('No noise data available for noise decorrelation')
        return
    elseif nnz(noiseIdx) > 1
        error('More than one noise measurement')
    end

    % Decorrelate k-space data
    for iAcq = 1:length(data.data)
        if ~data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT', iAcq)
            data.data{iAcq} = data.data{iAcq} / decorrMat;
        end
    end

    if nargout >= 3
        % Decorrelate coil sensitivity maps
        cmSz = size(coilMaps);
        coilMaps = reshape(coilMaps, prod(cmSz(1:3)), cmSz(4));
        coilMaps = coilMaps / decorrMat;
        coilMaps = reshape(coilMaps, cmSz);
    end
end
