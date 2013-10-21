% DS = DISTANCES(XS, YS)
%
% Calculate squared distances between groups of points.  Suppose XS and YS
% are matrices of dimension K * N and K * M, respectively.  Consider the
% columns of XS and YS as points in K-dimensional space.  DS will then be
% an N * M matrix, with D(I,J) equal to the squared distance between
% XS(:,I) and YS(:,J).
%
% DS = DISTANCES(XS)
%
% With a single argument, assumes YS = XS, making a square matrix of
% inter-point squared distances.

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

function ds = distances(xs, ys)

if (nargin == 1)
    ys = xs;
end

if (size(xs,1) ~= size(ys,1))
    error('distances:dimension', 'XS and YS must be same height (%d, %d)\n', size(xs,1), size(ys,1));
end

n = size(xs, 2);
m = size(ys, 2);

ds = zeros(n, m);

for i = 1:n
    ds(i,:) = sum((ys - repmat(xs(:,i), 1, m)).^2,1);
end
