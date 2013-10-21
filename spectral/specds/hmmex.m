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

%% make a k-state HMM, and set some parameters

% latent dimension k, window for past and future, ridge parameter lambda
k = 10;
npast = 10;
nfut = 5;
lambda = 1;

fprintf('Building HMM\n')
nstates = 5;
nobs = 7;
ths = (1:nstates)*2*pi/nstates;
statepos = [cos(ths); sin(ths)];
obspos = [zeros(1,nobs); ((1:nobs)*2/(nobs+1))-1];
% obspos = 2*rand(2,nobs)-1;  % comment in to randomize observations too
T = exp(-2*distances(statepos));
T = T / diag(sparse(sum(T)));
O = exp(-3*distances(obspos, statepos));
O = O / diag(sparse(sum(O)));
st1 = zeros(nstates,1);
st1(1) = 1;

%% Visualize the HMM

figure(1);
clf;

subplot('position', [.05 .2 .6 .6]);
basewidth = 5;
for i = 1:nstates
    for j = 1:nstates
        lw = sqrt(T(i,j))*basewidth;
        if (lw > .75 && i ~= j)
            line([statepos(1,i) statepos(1,j)], [statepos(2,i) statepos(2,j)], 'linewidth', lw, 'color', max(0,1 - [1 1 1]*lw*.5));
        end
    end
end
line(statepos(1,:), statepos(2,:), 'marker', 'o', 'markersize', 20, 'linestyle', 'none', 'linewidth', 2);
line(obspos(1,:), obspos(2,:), 'marker', 'x', 'markersize', 15, 'linestyle', 'none', 'linewidth', 2);
axis equal square off

subplot('position', [.7 .07 .2 .4]);
imagesc(T);
axis equal off

subplot('position', [.7 .5 .2 .4]);
imagesc(O);
axis equal off

set(gcf, 'Color', [1 1 1])
colormap default
pause(.5);

%% sample from the HMM

fprintf('Sampling from the HMM\n');
nsteps = 50000;
obseq = samplehmm(T, O, st1, nsteps);
obsmat = accumarray([obseq (1:nsteps)'], 1);

%% display matrix of shifted observations (visualization -- not used later)

XY = makeshifts(obsmat(:,101:300), npast+nfut);
X = XY(1:npast*nobs,:);
Y = XY(npast*nobs+1:end,:);

figure(1)
clf;
imagesc(kron(XY,ones(3)));
axis equal tight off
horizline(3*npast*nobs+.5, 'g', 1);
set(gca, 'FontSize', 20);
set(gcf, 'Color', [1 1 1])
colormap default
pause(.5)

%% learn the HMM

fprintf('Learning the HMM\n')
modl = specds(obsmat, [], npast, nfut, k, lambda);

%% compute some predictions from true model and learned model, compare

% true model
stat = ones(nstates,1)/nstates;
trupreds = zeros(nobs,nobs*nobs);
for i = 1:nobs
    for j = 1:nobs
        pr = O * T * diag(sparse(O(j,:))) * T * diag(sparse(O(i,:))) * stat;
        trupreds(:,(i-1)*nobs+j) = pr;
    end
end

% learned model
preds = zeros(nobs,nobs*nobs);
for i = 1:nobs
    for j = 1:nobs
        st = modl.sbar;
        ei = zeros(nobs,1);
        ei(i) = 1;
        ej = zeros(nobs, 1);
        ej(j) = 1;
        ocov = tensprod(modl.tro, st, 2, 1);
        pr = ocov(i,i);
        st = tensprod(modl.trso, st, 2, 1) * (ocov \ ei);
        ocov = tensprod(modl.tro, st, 2, 1);
        pr = pr * ocov(j,j);
        st = tensprod(modl.trso, st, 2, 1) * (ocov \ ej);
        pr = pr * diag(tensprod(modl.tro, st, 2, 1));
        preds(:,(i-1)*nobs+j) = pr;
    end
end

% comparison plot
figure(1)
clf;
subplot(3,1,1);
imagesc(kron(preds,ones(3))); 
axis equal tight off;
set(gca, 'FontSize', 20);
set(gcf, 'Color', [1 1 1])
title Predictions
subplot(3,1,2);
imagesc(kron(trupreds,ones(3))); 
axis equal tight off;
set(gca, 'FontSize', 20);
set(gcf, 'Color', [1 1 1])
colormap default
title('Ground Truth')
pause(.5)


%% compute and plot convergence of prediction error

fprintf('Making convergence plot\n')
ns = round(exp(log(300):.5:log(1000000)));
runs = 3;
errs = zeros(runs,length(ns));
for run = 1:runs;
    fprintf('Run %d\n', run);
    thisobseq = samplehmm(T, O, st1, max(ns));
    thisobsmat = accumarray([thisobseq (1:max(ns))'], 1);
    for iter = 1:length(ns);
        nn = ns(iter);
        fprintf('learning from first %d observations\n',nn);
        thismod = specds(thisobsmat(:,1:nn), [], npast, nfut, k, lambda);
        preds = zeros(nobs,nobs*nobs);
        for i = 1:nobs
            for j = 1:nobs
                st = thismod.sbar;
                ei = zeros(nobs,1);
                ei(i) = 1;
                ej = zeros(nobs, 1);
                ej(j) = 1;
                ocov = tensprod(thismod.tro, st, 2, 1);
                pr = ocov(i,i);
                st = tensprod(thismod.trso, st, 2, 1) * (ocov \ ei);
                ocov = tensprod(thismod.tro, st, 2, 1);
                pr = pr * ocov(j,j);
                st = tensprod(thismod.trso, st, 2, 1) * (ocov \ ej);
                pr = pr * diag(tensprod(thismod.tro, st, 2, 1));
                preds(:,(i-1)*nobs+j) = pr;
            end
        end
        errs(run, iter) = mean(abs(preds(:)-trupreds(:)));
    end
end

% scaling factor: error if we were to predict uniform
errscale = mean(abs(trupreds(:)-1/nobs^3));

% plot it
figure(1)
clf
loglog(ns, errs/errscale, ns, 30*ns.^-.5, 'k--', 'linewidth', 2)
legend({'Run 1', 'Run 2', 'Run 3', 'c/sqrt(n)'});
axis([min(ns) max(ns) 5e-3 1]);
set(gca, 'FontSize', 20);
set(gcf, 'Color', [1 1 1])
colormap default
