function [deformField, imWarped] = transformix(transfParamFile)
%transformix Retrieve deformation field from elastix registration, and
%optionally warp the moving image
%
%   Inputs:
%   transfParamFile: Path to transformation parameters file
%
%   Outputs:
%   deformField: Deformation field
%   imWarped: Warped moving image
%   
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    transfParamFile {mustBeTextScalar, mustBeFile}
end

assert(isenv('ELASTIXPATH'), 'Environment variable ELASTIXPATH must be set')

outDir = fileparts(transfParamFile);
if nargout >= 2
    movingFileName = fullfile(outDir, 'moving.mhd');
    if ~isfile(movingFileName)
        error('Moving image file %s not found in transformix directory', movingFileName)
    end
end

% Create command to run transformix
cmd = 'LD_LIBRARY_PATH=$ELASTIXPATH/lib';
cmd = sprintf('%s && $ELASTIXPATH/bin/transformix -tp %s', cmd, transfParamFile);
cmd = append(cmd, ' -def all');
cmd = sprintf('%s -out %s', cmd, outDir);
if nargout >= 2
    cmd = sprintf('%s -in %s', cmd, movingFileName);
end

% Run transformix
status = system(cmd);
if status ~= 0
    error('Transformix failed with status %d', status)
end

% Read output
deformField = readMhd(fullfile(outDir, 'deformationField.mhd'));
if nargout >= 2
    imWarped = readMhd(fullfile(outDir, 'result.mhd'));
end

end
