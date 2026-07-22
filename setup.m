% Sets all relevent paths and environment variables to run the scripts in
% this repository.
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

addpath(genpath('./src'))
addpath('./colormaps')
% addpath('path/to/ismrmrd/matlab')
if ~isenv('ELASTIXPATH')
    % setenv('ELASTIXPATH', 'path/to/elastix')
end

% Check if ISMRMRD and elastix are available
if ~exist('ismrmrd.Dataset', 'class')
    error('ISMRMRD library cannot be found; add its path to setup.m')
end
if ~exist(fullfile(getenv('ELASTIXPATH'), 'bin', 'elastix'), 'file')
    error('Elastix cannot be found; make sure ELASTIXPATH is set correctly in setup.m')
end
