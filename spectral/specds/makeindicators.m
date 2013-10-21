% FUNCTION [IDX, KEY] = MAKEINDICATORS(XS, KEY, SPARS)
%
% Given a vector XS (either a 1D array of numbers, or a 1D cell array of
% strings), produce a corresponding array IDX of indicator variables based
% on the unique elements of XS.  IDX will have one row per element of XS,
% and its rows will reference the unique elements of XS.
%
% Optional second argument KEY lists the possible values for elements of
% XS.  If missing or empty, KEY defaults to a sorted list of the unique
% elements of XS.  Whether provided or not, KEY is returned as the second
% output value, so that we can use it in a future call to get matching
% indicator variables.
%
% There are two formats for IDX, selected by the optional third argument
% SPARS: if SPARS is true (the default), IDX will have one column for each
% element of KEY.  There will be a single 1 in each row of IDX: for row I,
% if XS(I)==KEY(J), then IDX(I,J)==1.  On the other hand, if SPARS is
% false, IDX will have just one column, and its value will be an index into
% KEY: that is, KEY(IDX(I))==XS(I).
%
% Rows of XS which are not matched in KEYS will result in all-zero rows of
% IDX.  Elements of KEYS which do not appear in XS will result in all-zero
% columns of IDX (in the sparse format), or index values that don't appear
% in IDX (in the dense format).
%
% Examples:
%
% >> full(makeindicators({'red', 'green', 'green', 'blue'}, {'red','green'}))
% 
% ans =
% 
%      1     0
%      0     1
%      0     1
%      0     0
%
% >> [idx,key] = makeindicators({'red', 'green', 'green', 'blue'}, [], 0)
% 
% idx =
% 
%      3
%      2
%      2
%      1
% 
% 
% key = 
% 
%     'blue'    'green'    'red'

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

function [idx, key] = makeindicators(xs, key, spars)

% check whether a key was passed in
if (nargin < 2 || isempty(key))
    % if no key passed in, build a key and a map KEY2XS showing which key
    % is in each row of XS
    [key,~,key2xs] = unique(xs);
else
    % if so, we only need to build KEY2XS
    [kk,~,kk2xs] = unique(xs);
    [~, loc] = ismember(kk, key);
    key2xs = loc(kk2xs);
end

% default the format argument
if (nargin < 3)
    spars = true;
end

% row and column indices of 1s in IDX, ignoring possibility of unmatched
% rows of XS
is = 1:length(key2xs);
js = key2xs;

% Delete entries of IS and JS corresponding to 0s in KEY2XS, since these
% are entries of XS with no matching entry in KEY
mask = js > 0;
is = is(mask);
js = js(mask);

% construct and return the indicator matrix
if (spars)
    idx = sparse(is, js, ones(size(js)), length(xs), length(key));
else
    idx = accumarray([is(:) ones(size(is(:)))], js, [length(xs) 1]);
end
