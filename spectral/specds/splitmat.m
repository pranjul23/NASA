% SUBX = SPLITMAT(X, COLS, PADSTART, PADEND, TPS)
%
% Split a matrix X into a list of submatrices, based on the unique values
% of one or more columns.
%
% Output SUBX is a cell array, with one element per unique value of
% X(:,COLS).  Each element of SUBX consists of a subset of the rows of X
% for which X(:,COLS) is identical; the identical entries (corresponding to
% the columns in COLS) will be stripped off.
%
% If optional arguments PADSTART and PADEND are zero (the default), the
% total number of rows in all elements of SUBX will equal the number of
% rows of X.  Each element of SUBX will have the same number of columns,
% equal to the number of columns in X minus the number of elements of COLS.
%
% If PADSTART > 0, each element of SUBX will be padded at the beginning
% with PADSTART copies of the row [0 0 0 ... 1].  If PADEND > 0, each
% element of SUBX will be padded at the end with PADEND copies of the row
% [0 0 0 ... 1].  In both cases, the single 1 will be in a new column,
% distinct from all the original columns.  (If both PADSTART > 0 and PADEND
% > 0, we will have two new columns, the first for padding at the start,
% and the second for padding at the end.)
%
% If optional last argument TPS is true, transpose each element of SUBX.

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


function subX = splitmat(X, cols, padstart, padend, tps)

if (nargin < 3)
    padstart = 0;
end

if (nargin < 4)
    padend = 0;
end

if (nargin < 5)
    tps = false;
end

[~,~,seqidx] = unique(X(:,cols),'rows');
othercols = setdiff(1:size(X,2), cols);
numseqs = max(seqidx);
subX = cell(numseqs,1);
for i = 1:numseqs
    thisx = X(seqidx==i, othercols);
    if (padstart > 0)
        [n,m] = size(thisx);
        thisx = [zeros(padstart, m) ones(padstart,1); thisx zeros(n,1)]; 
    end
    if (padend > 0)
        [n,m] = size(thisx);
        thisx = [thisx zeros(n,1); zeros(padend, m) ones(padend,1)];
    end
    if (tps)
        subX{i} = thisx';
    else
        subX{i} = thisx;
    end
end

%#ok<*AGROW> (the size-changing of thisx is a feature not a bug)
