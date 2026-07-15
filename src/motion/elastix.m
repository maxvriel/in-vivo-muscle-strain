function [imWarped, deformField] = elastix(imMoving, spacingMoving, imFixed, spacingFixed, outDir, paramFiles, maskFixed, maskMoving)
%elastix Perform image registration with elastix
%
%   Inputs:
%   imMoving: Moving image
%   spacingMoving: Voxel spacing of the moving image
%   imFixed: Fixed image
%   spacingFixed: Voxel spacing of the fixed image
%   outDir: Output directory for elastix results; must be empty
%   paramFiles: Path(s) to elastix parameter file; specify multiple files
%       to run multiple registrations in succession
%   maskFixed (optional): Fixed image mask
%   maskMoving (optional): Moving image mask
%   
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    imMoving {mustBeNumeric}
    spacingMoving (1,:) {mustBePositive}
    imFixed {mustBeNumeric}
    spacingFixed (1,:) {mustBePositive}
    outDir {mustBeTextScalar}
    paramFiles (1,:) {mustBeFile}
    maskFixed {mustBeNumericOrLogical} = []
    maskMoving {mustBeNumericOrLogical} = []
end

assert(isenv('ELASTIXPATH'), 'Environment variable ELASTIXPATH must be set')

% Make sure output folder exists and is empty
if ~isfolder(outDir)
    mkdir(outDir)
end
existingFiles = dir(fullfile(outDir, '*'));
existingFiles = existingFiles(~[existingFiles.isdir]);
if ~isempty(existingFiles)
    error('Output directory is not empty')
end

% Create command to run elastix
cmd = sprintf('LD_LIBRARY_PATH=$ELASTIXPATH/lib && $ELASTIXPATH/bin/elastix -out %s', outDir);
cmd = append(cmd, char(join(append(' -p ', paramFiles), '')));
fixedFileName = fullfile(outDir, 'fixed.mhd');
writeMhd(fixedFileName, imFixed, spacingFixed);
movingFileName = fullfile(outDir, 'moving.mhd');
writeMhd(movingFileName, imMoving, spacingMoving);
cmd = sprintf('%s -f %s -m %s', cmd, fixedFileName, movingFileName);
if ~isempty(maskFixed)
    maskFileName = fullfile(outDir, 'fMask.mhd');
    writeMhd(maskFileName, maskFixed, spacingFixed);
    cmd = sprintf('%s -fMask %s', cmd, maskFileName);
end
if ~isempty(maskMoving)
    maskFileName = fullfile(outDir, 'mMask.mhd');
    writeMhd(maskFileName, maskMoving, spacingMoving);
    cmd = sprintf('%s -mMask %s', cmd, maskFileName);
end
% Run elastix
status = system(cmd);
if status ~= 0
    error('Elastix failed with status %d', status)
end

% Run transformix to get output image and/or deformation field
transfParamFile = fullfile(outDir, sprintf('TransformParameters.%d.txt', length(string(paramFiles))-1));
if nargout == 1
    imWarped = transformix(transfParamFile);
elseif nargout >= 2
    [imWarped, deformField] = transformix(transfParamFile);
end

end
