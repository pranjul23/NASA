% Y = MAKESHIFTS(X, k)
%
% Make a Hankel matrix by shifting X and stacking the results k times.

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

function Y = makeshifts(X, k)

[m, n] = size(X);
Y = zeros(m*k, n-k+1);
for i = 1:k
    Y(1+(i-1)*m:i*m, :) = X(:, i:n-k+i);
end
