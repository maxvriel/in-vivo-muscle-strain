function axPos = createAxesPositions(widths, heights)
%createAxesPositions Calculate the Position property for Axes objects in a
%grid given a vector of sizes of the Axes and the gaps in between
%   Inputs:
%   widths: Vector of relative widths: [left, plot, gap, plot, ..., right]
%   heights: Vector of relative heights: [top, plot, gap, plot, ..., bottom]
%
%   Outputs:
%   axPos: Cell array with one Position vector for each Axes
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

assert(mod(length(widths), 2) == 1, 'Input widths must have an odd number of elements')
assert(mod(length(heights), 2) == 1, 'Input heights must have an odd number of elements')

nw = (length(widths)-1) / 2;
nh = (length(heights)-1) / 2;
axPos = cell(nh, nw);
norm = [sum(widths), sum(heights), sum(widths), sum(heights)];

for iw = 1:nw
    for ih = 1:nh
        left = sum(widths(1:(iw-1)*2+1));
        height = sum(heights(end-(nh-ih)*2:end));
        axPos{ih,iw} = [left, height, widths(iw*2), heights(end-(nh-ih)*2-1)] ./ norm;
    end
end

end
