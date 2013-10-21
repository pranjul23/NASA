% Z = TENSPROD(X, Y, J, K)
%
% Tensor product: reduce along dimension J of X and K of Y.  (J and K
% default to 1.)  Z's first NDIMS(X)-1 dimensions correspond to the
% remaining dimensions of X; Z's last NDIMS(Y)-1 dimensions correspond to
% the remaining NDIMS(Y)-1 dimensions of Y.

% Copyright 2012 Geoff Gordon
%
% This file is part of SpecDS.
%
% SpecDS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% SpecDS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with SpecDS.  If not, see <http://www.gnu.org/licenses/>.

function z = tensprod(x, y, j, k)

% default arguments
if (nargin < 3)
    j = 1;
    k = 1;
end

% permute X, Y so that the reduction dimension is first
if (j ~= 1)
    x = permute(x, [j 1:j-1 j+1:ndims(x)]);
end
if (k ~= 1)
    y = permute(y, [k 1:k-1 k+1:ndims(y)]);
end

xdims = size(x);
ydims = size(y);
x = reshape(x, [xdims ones(1, ndims(y)-1)]);
y = reshape(y, [ydims(1) ones(1, ndims(x)-1) ydims(2:end)]);
x = repmat(x, [ones(size(xdims)) ydims(2:end)]);
y = repmat(y, [1 xdims(2:end) ones(1,ndims(y)-1)]);

zdims = [xdims(2:end) ydims(2:end)];
z = sum(x .* y, 1);
z = reshape(z, zdims);

