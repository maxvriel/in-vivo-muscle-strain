addpath(genpath('./src'))
addpath('./colormaps')
addpath('~/Code/external/ismrmrd/matlab')
if ~isenv('ELASTIXPATH')
    setenv('ELASTIXPATH', '/localscratch/mriel3/Software/elastix-latest')
end

% Check if ISMRMRD and elastix are available
if ~exist('ismrmrd.Dataset', 'class')
    error('ISMRMRD library cannot be found')
end
if ~exist(fullfile(getenv('ELASTIXPATH'), 'bin', 'elastix'), 'file')
    error('Elastix cannot be found; make sure ELASTIXPATH is set correctly')
end
