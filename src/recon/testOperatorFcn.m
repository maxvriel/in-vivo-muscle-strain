function testOperatorFcn(modelFcn, inputSize, outputSize, complex)
%testOperatorFcn Test the adjoint property of an operator.
%   For any operator A, the adjoint property states that:
%   <Ax, y> = <x, A^H y> for all x and y.
%
%   Inputs:
%   modelFcn: Function handle for the operator to test
%   inputSize: Input size of the operator
%   outputSize: Output size of the operator
%   complex (optional): Boolean indicating whether to use complex inputs,
%       defaults to true
%
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

arguments
    modelFcn (1,1) function_handle
    inputSize (1,:) {mustBePositive, mustBeInteger}
    outputSize (1,:) {mustBePositive, mustBeInteger}
    complex (1,1) logical = true
end

% Create random inputs
x = randn(inputSize);
y = randn(outputSize);
if complex
    x = x + 1i*randn(inputSize);
    y = y + 1i*randn(outputSize);
end

% Pass inputs through the model function and its adjoint
Ax = modelFcn(x, 'notransp');
r1 = Ax(:)'*y(:);
AHy = modelFcn(y, 'transp');
r2 = x(:)'*AHy(:);

% Check adjoint property
assert(abs(r1 - r2) < 1e3*eps(max([r1,r2])), 'Adjoint test failed')

end
