function y = evalBSpline(knots, controlPoints, x)
%evalBSpline Evaluates a 1D B-spline using de Boor's algorithm.
%
%   Inputs:
%   knots: Knot vector, must be monotonically non-decreasing
%   controlPoints: Control points of the B-spline
%   x: Locations at which to evaluate the B-spline
%
%   Outputs:
%   y: Values of the B-spline at the given locations x
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    knots (:,1) {mustBeNumeric, mustBeReal, mustBeNonDecreasing}
    controlPoints (:,1) {mustBeNumeric}
    x {mustBeNumeric, mustBeReal}
end

nKnots = length(knots);
nCp = length(controlPoints);
order = length(knots) - length(controlPoints);
p = order-1;
assert(nKnots >= 2, 'At least two knots required')
assert(nKnots > nCp, 'Number of knots must be larger than number of control points')

xSz = size(x);
x = x(:);

% Find the interval index for all locations x
idx = discretize(x, knots);
validIdx = ~isnan(idx);
idx = idx(validIdx);
x = x(validIdx);

% Perform de Boor's algorithm
d = zeros(length(idx), p+1, 'like', controlPoints);
cpIdx = idx + (-p:0);
validCp = cpIdx >= 1 & cpIdx <= nCp;
d(validCp) = controlPoints(cpIdx(validCp));
for r = 1:p
    for j = p:-1:r
        i1 = j+idx-p;
        i2 = j+idx-r+1;
        
        alpha = zeros(length(x), 1, 'like', x);
        beta = zeros(length(x), 1, 'like', x);
        for ix = 1:length(x)
            if d(ix,j+1) ~= 0
                if knots(i2(ix)) == knots(i1(ix))
                    alpha(ix) = 0;
                else
                    alpha(ix) = (x(ix) - knots(i1(ix))) ./ (knots(i2(ix)) - knots(i1(ix)));
                end
            end
            if d(ix,j) ~= 0
                if knots(i2(ix)) == knots(i1(ix))
                    beta(ix) = 1;
                else
                    beta(ix) = (knots(i2(ix)) - x(ix)) ./ (knots(i2(ix)) - knots(i1(ix)));
                end
            end
        end

        d(:,j+1) = beta.*d(:,j) + alpha.*d(:,j+1);
    end
end

% Elements outside the spline range are zero
y = zeros(prod(xSz), 1, 'like', d);
y(validIdx) = d(:,end);

y = reshape(y, xSz);

end

function mustBeNonDecreasing(knots)
    %mustBeNonDecreasing Validation function for knot vector

    if ~all(diff(knots) >= 0)
        error('Knot values must be monotonically non-decreasing')
    end

end
