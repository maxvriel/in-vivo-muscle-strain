function writeMhd(fileName, data, spacing, isVectorField, transfMatrix)
%writeMhd Write data to a MetaImage (.mhd) file.
%
%   Inputs:
%   fileName: Path to the MetaImage header file (.mhd)
%   data: Image or vector field data
%   spacing: Voxel spacing in each spatial dimension in mm
%   isVectorField (optional): Boolean indicating if the data is a vector
%       field, defaults to false
%   transfMatrix (optional): Transformation matrix between the image and
%       world coordinates, defaults to identity
%   
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    fileName {mustBeTextScalar}
    data {mustBeNumericOrLogical, mustBeReal}
    spacing (1,:) {mustBePositive}
    isVectorField (1,1) logical = false
    transfMatrix (:,:) = eye(ndims(data)-isVectorField)
end

assert(endsWith(fileName, '.mhd'), 'Input file must be a .mhd file')

nDims = ndims(data);
if isVectorField
    nDims = nDims - 1;
end
assert(ismember(nDims, [2, 3]), 'Data must be 2D or 3D')

if isa(data, 'single')
    precision = 'float';
elseif isa(data, 'double')
    precision = 'double';
elseif isa(data, 'uint8') || islogical(data)
    precision = 'uchar';
elseif isa(data, 'int8')
    precision = 'char';
elseif isa(data, 'uint16')
    precision = 'ushort';
elseif isa(data, 'int16')
    precision = 'short';
elseif isa(data, 'uint32')
    precision = 'uint';
elseif isa(data, 'int32')
    precision = 'int';
else
    error('Unsupported data type: %s', class(data))
end

offset = (spacing(1:nDims) .* -floor(size(data,1:nDims)/2)) * transfMatrix';

[folder, baseFileName] = fileparts(fileName);
if ~isfolder(folder)
    mkdir(folder)
end

% Create header structure (convert m to mm)
info.ObjectType = "Image";
info.NDims = nDims;
info.BinaryData = true;
info.BinaryDataByteOrderMSB = false;
info.CompressedData = false;
info.TransformMatrix = transfMatrix(:).';
info.Offset = offset;
info.CenterOfRotation = zeros(1, nDims);
info.ElementSpacing = spacing(1:nDims);
info.DimSize = size(data, 1:nDims);
if isVectorField
    info.ElementNumberOfChannels = size(data, nDims+1);
else
    info.ElementNumberOfChannels = 1;
end
info.ElementType = sprintf("MET_%s", upper(precision));
info.ElementDataFile = sprintf("%s.raw", baseFileName);

% Convert structure to string array
fields = fieldnames(info);
lines = repmat("", length(fields)+1, 1);
for iField = 1:length(fields)
    prop = fields{iField};
    val = info.(prop);
    if isempty(val)
        % Do not write anything
        continue
    elseif ischar(val) || isstring(val)
        % Write string
        val = join(compose('%s', val));
    elseif isnumeric(val)
        if all(mod(val,1) == 0)
            % Write integers
            val = join(compose('%d', val));
        else
            % Write floating points (only up to the 6th decimal)
            if any(abs(val) < 1e-5)
                warning('Possible loss of precision when writing property %s', prop)
            end
            val = join(compose('%.6f', val));
        end
    elseif islogical(val)
        % Write logical strings by indexing into array
        logicalString = ["False", "True"];
        val = join(compose('%s', logicalString(val+1)));
    end
    lines(iField) = sprintf("%s = %s", prop, val{1});
end

% Write header file
writelines(lines, fileName);

% Permute data
if isVectorField
    data = permute(data, [nDims+1, 1:nDims]);
end

% Write raw data file
if strcmp(precision, 'char')
    precision = 'uchar';
end
fid = fopen(fullfile(folder, info.ElementDataFile), 'w');
if fid == -1
    error('Could not open raw data file for writing')
end
fwrite(fid, data, precision);
fclose(fid);

end
