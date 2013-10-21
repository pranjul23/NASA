% STATES = QUICKSTATES(MODL, X, SKIP)
%
% A quick way to compute approximate predictive states.  Given a model MODL
% as returned by SPECDS, estimate predictive states for a sequence of
% observations X by simple linear projection of history.  Rows of X
% correspond to observation dimensions; columns correspond to time steps. X
% should have MODL.nfeats rows, and at least MODL.npast columns.
%
% X may also be a cell array of observation sequences, in which case we
% compute state estimates for each element X{I} separately.
%
% The optional final argument SKIP (defaults to 0) tells us to skip this
% many observations from the end of each sequence -- useful if we want to
% compute only those states for which we can test corresponding predictions
% against reality.
%
% The output is a cell array of matrices of state estimates. STATES{I} will
% have MODL.k+1 rows (corresponding to latent state dimensions, the last of
% which is constant) and SIZE(X{I},2)-MODL.npast+1-SKIP columns
% (corresponding to states immediately after observation numbers MODL.npast
% through SIZE(X{I},2)-SKIP).

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

function states = quickstates(modl, seqs, skip)

% default argument
if (nargin < 3)
    skip = 0;
end

% a single observation matrix is equivalent to a cell array containing just
% that matrix
if (~iscell(seqs))
    seq = seqs;
    seqs = cell(1);
    seqs{1} = seq;
end

% must provide at least one observation matrix
if (length(seqs) < 1)
    error('specds:empty', 'must provide at least one observation matrix');
end

% extract some parameters from the model
rrr = modl.past2state;
npast = modl.npast;
nfeat = modl.nobs;

% compute the states
states = cell(length(seqs),1);
for ss = 1:length(seqs)
    if (npast > size(seqs{ss},2))
        continue;
    end
    st = 0;
    for i = 1:npast
        st = st + rrr(:,(i-1)*nfeat+1:i*nfeat) * seqs{ss}(:, i:end-npast+i-skip);
    end
    states{ss} = [st; ones(1, size(st,2))];
end
