function data = readMhd(fileName)
% readMhd Read a MetaImage (.mhd) file
%
%   Inputs:
%   fileName: Path to the MetaImage header file (.mhd)
%
%   Outputs:
%   data: Image or vector field data
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    fileName {mustBeTextScalar, mustBeFile}
end

assert(endsWith(fileName, '.mhd'), 'Input file must be a .mhd file')

% Read header file
lines = readlines(fileName);

% Put header info in structure
info = struct();
for iLine = 1:length(lines)
    propval = split(lines(iLine), ' = ');
    if length(propval) == 2
        prop = propval(1);
        val = propval(2);
        valnum = str2double(split(val).');
        if ~any(isnan(valnum))
            val = valnum;
        end
        if strcmpi(val, 'True')
            val = true;
        elseif strcmpi(val, 'False')
            val = false;
        end
        info.(prop) = val;
    end
end

if ~isfield(info, 'ElementNumberOfChannels')
    info.ElementNumberOfChannels = 1;
end

% Read raw data file
folder = fileparts(fileName);
rawFileName = fullfile(folder, info.ElementDataFile);
if ~isfile(rawFileName)
    error('Cannot locate raw data file')
end
if info.CompressedData
    error('Cannot read compressed data')
end
precision = lower(erase(info.ElementType, 'MET_'));
if strcmp(precision, 'char')
    precision = 'uchar';
end
fid = fopen(rawFileName);
if fid == -1
    error('Could not open raw data file for reading')
end
data = fread(fid, sprintf('*%s', precision));
fclose(fid);

% Reshape data and convert mm to m
if info.ElementNumberOfChannels == 1
    data = reshape(data, info.DimSize);
else
    data = reshape(data, [info.ElementNumberOfChannels, info.DimSize]);
    data = permute(data, [2:info.NDims+1, 1]);
    data = data * 1e-3;
end

end
