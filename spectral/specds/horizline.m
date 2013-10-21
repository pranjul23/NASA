% H = HORIZLINE(Y)
%
% Draw a horizontal line across the current plot at the specified height.
% Return a handle to the line.
%
% H = HORIZLINE(Y, COLOR, WIDTH, ...)
%
% Optional arguments allow specifying line color and thickness.  Defaults
% are black and 2.  Additional arguments are property-value pairs that get
% passed to line() -- e.g., 'LineStyle', ':'.

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


function h = horizline(y, color, width, varargin)

if (nargin < 2)
    color = 'k';
end

if (nargin < 3)
    width = 2;
end

ax = axis;

h = line([ax(1) ax(2)], [y y], 'Color', color, 'LineWidth', width, varargin{:});
