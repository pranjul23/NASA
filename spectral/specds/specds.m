% MODEL = SPECDS(X, LABELS, NPAST, NFUT, K, LAMBDA)
%
% Spectral learning of a dynamical system from observations X.  X may be a
% matrix, in which case each column of X is an observation from a single
% time step in a trajectory of the dynamical system.  X may also be a cell
% array of matrices, corresponding to multiple trajectories sampled from
% the same start state.
%
% LABELS selects between unsupervised and semisupervised learning.  If
% LABELS is [], we do unsupervised learning.  Else, if X is a cell array,
% then LABELS should be a cell array the same size as X; LABELS{I} should
% have either 1 column (corresponding to supervision at the level of whole
% trajectories) or the same number of columns as X{I} (corresponding to
% supervision at the level of time steps).  If X is a matrix, then LABELS
% may also be a matrix (with the same number of columns as X).
%
% Estimate states by using a window of NPAST steps in the past to predict a
% window of NFUT steps in the future.  Labels (if present) are included
% only in the future (and so they need not be present at test time).  Learn
% a dynamical system with (K+1)-dimensional states (of which the last
% dimension is a constant).  Regularization parameter LAMBDA > 0.
%
% Return value MODEL is a struct with fields:
%
%   trso: the transition parameters
%   tro: the observation parameters
%   sbar: the stationary state
%   s1: the initial state
%   svs: the proportion of variance explained by each state dimension
%   past2state: transformation from history to the learned latent state
%   futbasis: the transformation from the learned latent state to a vector
%     of predictions about the future
%   obsbasis: a dimension reduction for the observation vector (if needed)
%   npast, nfut, k, lambda: copied from the input
%   nobs: dimension of observation vector
%   nlabs: dimension of label vector (if present)
%
% Note that error in estimating s1 only goes to zero with the number of
% sequences (the length of X if X is a cell array, or 1 if X is a matrix);
% error in all other parameters goes to zero with the total number of
% observations (the total number of columns in all elements of X, often
% much larger).  So, s1 may be unusable even when other parameters are
% reasonable estimates.  Typical RMSE is O(1/sqrt(n)), where n is number of
% sequences or number of observations.

% to do:
%   conditioning variables / actions
%   assumption now is large T *or* high-dim observations; what if both?
%   eliminate common subexpressions in computation of covariance blocks
%   visualize learned model and learned states better?
%   save RRR and/or states

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


function model = specds(seqs, labels, npast, nfut, k, lambda)

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

% a single label matrix is equivalent to a cell array containing just that
% matrix
if (~isempty(labels) && ~iscell(labels))
    lab = labels;
    labels = cell(1);
    labels{1} = lab;
end

% must provide same number of label matrices as observation matrices
if (~isempty(labels) && length(labels) ~= length(seqs))
    error('specds:labels', 'if labels are provided, should have one label matrix per observation matrix');
end

% some parameters
nfeat = size(seqs{1},1);
if (isempty(labels))
    nlabs = 0;
else
    nlabs = size(labels{1},1);
end
maxlag = nfut+npast-1;

% check requested dimension
if (nfeat * npast < k || nfeat * nfut < k)
    error('specds:dim', 'requested dimension K=%d is larger than window dimension', k);
end

% collect all observations, and ensure that the number of features is <=
% the number of steps by making an orthogonal projection
allX = [seqs{:}];
nsteps = size(allX, 2);
if (nfeat > nsteps)
    [obsbasis,~,~] = svd(allX, 'econ');
    nfeat = nsteps;
    for ss = 1:length(seqs)
        seqs{ss} = obsbasis'*seqs{ss};
    end
else
    obsbasis = [];
end

% check label matrix sizes, and determine whether we are using
% trajectory-level or step-level supervision
nflab = 0;
if (~isempty(labels))
    obslevel = true;
    trajlevel = true;
    for ss = 1:length(seqs)
        sz = size(labels{ss}, 2);
        if (sz ~= 1)
            trajlevel = false;
        end
        if (sz ~= size(seqs{ss},2))
            obslevel = false;
        end
    end
    if (obslevel && trajlevel)
        trajlevel = false;
    end
    if (~(obslevel || trajlevel))
        error('specds:labsize', 'if labels provided, must have either one label per sequence or one label per step');
    end
    if (trajlevel)
        nflab = 1;
    else
        nflab = nfut;
    end
end

% initialize the output struct
model.npast = npast;
model.nfut = nfut;
model.k = k;
model.lambda = lambda;
model.obsbasis = obsbasis;
model.nlabs = nlabs;
model.nobs = size(seqs{1},1);

% go through all sequences to collect covariances
sigfp = zeros(nfeat*nfut,nfeat*npast);
siglp = zeros(nlabs*nflab,nfeat*npast);
sigpp = zeros(nfeat*npast);
count = 0;
skipped = 0;
for ss = 1:length(seqs)
    
    % check window size
    thisX = seqs{ss};
    thisnsteps = size(thisX,2);
    if (nfut+npast+1 > thisnsteps)
        skipped = skipped + 1;
        continue;
    else
        count = count + thisnsteps - maxlag;
    end
    
    % covariance between future and past
    for i = 1:nfut
        for j = 1:npast
            v = thisX(:,i+npast:end-maxlag+i+npast-1) * thisX(:,j:end-maxlag+j-1)';
            rs = nfeat*(i-1)+1:nfeat*i;
            cs = nfeat*(j-1)+1:nfeat*j;
            sigfp(rs,cs) = sigfp(rs,cs) + v;
        end
    end
    
    % covariance between labels and past
    if (~isempty(labels))
        thisL = labels{ss};
        for i = 1:nflab
            for j = 1:npast
                if (trajlevel)
                    v = thisL * sum(thisX(:,j:end-maxlag+j-1),2)';
                else
                    v = thisL(:,i+npast:end-maxlag+i+npast-1) * thisX(:,j:end-maxlag+j-1)';
                end                
                rs = nlabs*(i-1)+1:nlabs*i;
                cs = nfeat*(j-1)+1:nfeat*j;
                siglp(rs,cs) = siglp(rs,cs) + v;
            end
        end
    end
    
    % covariance of past
    for i = 1:npast
        v = thisX(:,i:end-maxlag+i-1) * thisX(:,i:end-maxlag+i-1)';
        rs = nfeat*(i-1)+1:nfeat*i;
        sigpp(rs,rs) = sigpp(rs,rs) + v;
        for j = i+1:npast
            v = thisX(:,j:end-maxlag+j-1) * thisX(:,i:end-maxlag+i-1)';
            cs = nfeat*(j-1)+1:nfeat*j;
            sigpp(rs,cs) = sigpp(rs,cs) + v';
            sigpp(cs,rs) = sigpp(cs,rs) + v;
        end
    end
end
if (skipped > 0)
    warning('specds:nsteps', 'skipped %d/%d observation sequences: must have at least %d time steps to allow npast=%d, nfut=%d', skipped, length(seqs), nfut+npast+1, npast, nfut);
end

% scale by sample size
sigpp = sigpp / count;
sigfp = sigfp / count;
siglp = siglp / count;

% Reduced-rank regression: SVD the covariance between future and whitened
% past to get a basis for the maximally predictable part of the future.
R = chol(sigpp+(lambda/(nsteps-maxlag))*eye(size(sigpp)));
[Ubas,svs,~] = svds([sigfp; siglp] / R, k);
UbasF = Ubas(1:size(sigfp,1),:);
UbasL = Ubas(size(sigfp,1)+1:end,:);
model.futbasis = UbasF*svs;
model.labelbasis = UbasL*svs;
futbinv = svs \ UbasF';
labbinv = svs \ UbasL';
svs = diag(svs);
if (length(svs) < k)
    warning('specds:svs', 'unable to extract %d svs; learned states are %d-dimensional', k, length(svs)+1);
    k = length(svs);
    model.k = k;
end
model.svs = svs;

% Regress from past to projected future to get a weight matrix which
% multiplies the past and produces predictive states -- we would prefer no
% regularization at all here, but to prevent numerical problems we use a
% tiny amount.
ridge = 1e-10*mean(abs(sigpp(:)))*speye(size(sigpp));
rrr = (futbinv * sigfp) / (sigpp + ridge);
model.past2state = rrr;

% Compute predictive states; estimate mean and initial state
states = quickstates(model, seqs, nfut+1);
s1 = 0;
ns1 = 0;
for ss = 1:length(states)
    if (isempty(states{ss}))
        continue;
    end
    s1 = s1 + states{ss}(:,1);
    ns1 = ns1 + 1;
end
states = cat(2,states{:});
model.sbar = mean(states,2);
model.s1 = s1 / ns1;


% states = cell(length(seqs),1);
% s1 = 0;
% ns1 = 0;
% for ss = 1:length(seqs)
%     if (nfut+npast+1 > size(seqs{ss},2))
%         continue;
%     end
%     st = 0;
%     for i = 1:npast
%         st = st + rrr(:,(i-1)*nfeat+1:i*nfeat) * seqs{ss}(:, i:end-npast-nfut+i-1);
%     end
%     states{ss} = st;
%     s1 = s1 + st(:,1);
%     ns1 = ns1 + 1;
% end
% states = cat(2,states{:});
% states = [states; ones(1, size(states,2))];
% model.sbar = mean(states,2);
% model.s1 = [s1 / ns1; 1];

% plot(states(1,:), states(2,:), 'x', 'linewidth', 2)

% Compute regression targets for future states
splus = cell(length(seqs), 1);
for ss = 1:length(seqs)
    if (nfut+npast+1 > size(seqs{ss},2))
        continue;
    end
    st = 0;
    for i = 1:nfut
        st = st + futbinv(:,(i-1)*nfeat+1:i*nfeat) * seqs{ss}(:,npast+i+1:end-nfut+i);
    end
    for i = 1:nflab
        if (trajlevel)
            st = st + repmat(labbinv * labels{ss}, 1, size(st,2));
        else
            st = st + labbinv(:,(i-1)*nlabs+1:i*nlabs) * labels{ss}(:,npast+i-1:end-nfut+i);
        end
    end
    splus{ss} = st;
end
splus = cat(2, splus{:});
splus = [splus; ones(1, size(splus,2))];

% Collect observations directly between states and splus
obs = cell(length(seqs),1);
for ss = 1:length(seqs)
    if (nfut+npast+1 > size(seqs{ss},2))
        continue;
    end
    obs{ss} = seqs{ss}(:,npast+1:end-nfut);
end
obs = cat(2, obs{:});

% matrix that regresses onto predictive state -- note no penalty on
% constant term since we know its variance exactly.
ridge = diag(sparse([ones(k,1); 0]));
projpast = states' / (states*states'+lambda*ridge);

% regressions to get dynamical system parameters
nsdim = size(states,1);
nodim = size(obs,1);
tro = zeros(nodim,nsdim,nodim);
trso = zeros(nsdim,nsdim,nodim);
for i = 1:nodim
    spluso = splus .* repmat(obs(i,:),nsdim,1);
    oo = obs .* repmat(obs(i,:),nodim,1);
    trso(:,:,i) = spluso * projpast;
    tro(:,:,i) = oo * projpast;
end
model.trso = trso;
model.tro = tro;


