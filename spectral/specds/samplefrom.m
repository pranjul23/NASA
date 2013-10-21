% IDX = SAMPLEFROM(D,N)
%
% Sample from a multinomial distribution.  D is a vector of probabilities,
% N is the number of samples desired.  IDX is a vector of N samples, each
% with distribution given by D, returned in sorted order.

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

function idx = samplefrom(d, n)

if (nargin < 2)
    n = 1;
end

d = d(:);
p = sort(rand(n,1));
cd = [cumsum(d); 5];			% 5 is a sentinel >1
idx = zeros(n,1);

thisidx = 1;
for i = 1:n
    while (p(i) > cd(thisidx))
        thisidx = thisidx + 1;
    end
    idx(i) = thisidx;
end


