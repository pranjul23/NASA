% OB = SAMPLEHMM(T, O, S1, NSTEPS)
%
% Sample a sequence of NSTEPS observations from an HMM with transition
% matrix T, observation matrix O, and initial state S1.

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

function obseq = samplehmm(T, O, s1, nsteps)

st = s1;
obseq = zeros(nsteps,1);
% stseq = zeros(nstates,nsteps);
for t = 1:nsteps
%     stseq(:,t) = st;
    o = samplefrom(O*T*st);
    st = O(o,:)' .* (T*st);
    st = st/sum(st);
    obseq(t) = o;
end
