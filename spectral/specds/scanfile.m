% [DAT, ...] = SCANFILE(FILENAME, HEADERLINE, ...)
%
% Read a tab-delimited file of data with optional header lines, and convert
% to a numerical array.  After the HEADERLINE argument, all subsequent
% arguments are column specifier strings, as described below; each column
% specifier produces a single column of DAT.
%
% FILENAME is the path to the file.  HEADERLINE is a nonnegative integer
% specifying which line in the file is the header.  Lines before the header
% line are skipped; data begins immediately after the header.  If
% HEADERLINE==0, there is no header, and columns can only be specified by
% integer indices.
%
% A specifier can be:
%
%   A string that matches a field in the header line.  The corresponding
%   column of the file is converted to real numbers (using sscanf('%f') on
%   each row).
%
%   A positive integer column index.  The corresponding column is converted
%   to real numbers as above.
%
%   A string ID:XXX where XXX is a field in the header line.  The
%   corresponding column is converted to small integer indicators, one
%   value for each unique element.  (E.g., a column with entries APPLE
%   BANANA BANANA CHERRY would be converted to [1; 2; 2; 3].)  The key for
%   interpreting the indicators is returned in the corresponding output
%   argument. (E.g., in the example above, the key would be {'APPLE',
%   'BANANA', 'CHERRY'}.) 
%
%   A negative integer N.  The column ABS(N) is converted to small integer
%   indicators as above, and a key is returned in the corresponding output
%   argument.
%
%   A struct containing two fields: SPEC and KEY.  SPEC is an integer
%   (interpreted as a column index) or a string (interpreted as a field
%   name).  KEY is a key as defined above.  Elements not mentioned in KEY
%   are defaulted to 0.  (UNIMPLEMENTED)
%
% Optional output arguments after the first are keys for interpreting the
% indicator variables, as described above.  Keys for non-indicator
% specifiers are [].

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


function [dat, varargout] = scanfile(filename, headerlines, varargin)

% read the header line (or first line of file, if header is absent) and
% split into fields
fid = fopen(filename);
for i = 1:headerlines-1
    fscanf(fid, '%[^\n\r] ', 1);
end
headerline = fscanf(fid, '%[^\n\r] ', 1);
fclose(fid);
headers = regexp(headerline, '\t', 'split');
fmt = repmat('%s', 1, length(headers));

% scan into a cell array of strings (skipping the header line)
fid = fopen(filename);
X = textscan(fid, fmt, 'Delimiter', '\t', 'Headerlines', headerlines);
fclose(fid);

% catch special case of empty file
if (length(X) < 1)
    dat = [];
    return;
end

% initialize output arguments
nout = max(nargout,1) - 1;
varargout = cell(nout,1);

% process each column specifier
nrows = length(X{1});
ncols = length(varargin);
dat = zeros(nrows, ncols);
for i = 1:ncols
    
    % printable representation of column specifier
    spec = varargin{i};
    if (isnumeric(spec))
        specstr = sprintf('%f', spec);
    elseif (isstruct(spec))
        specstr = sprintf('%s', spec.col);
    elseif (ischar(spec))
        specstr = spec;
    else
        specstr = 'not a number, string, or struct';
    end
    
    % parse column specifier
    mkind = 0;
    key = [];
    if (isstruct(spec))
        mkind = 1;
        if (ischar(spec.col))
            col = find(strcmp(spec.col, headers), 1);
        else
            col = spec.col;
        end
        key = spec.key;
    elseif (isnumeric(spec))
        if (spec < 0)
            mkind = 1;
            col = abs(spec);
        elseif (spec==0)
            error('scanfile:zero', 'column indices must be nonzero');
        else
            col = spec;
        end
    elseif (ischar(spec) && ~isempty(regexp(spec, '^ID:', 'once')))
        mkind = 1;
        col = find(strcmp(spec(4:end), headers), 1);
    elseif (ischar(spec))
        col = find(strcmp(spec, headers), 1);
    else
        error('scanfile:specifier', 'did not understand column specifier: %s', specstr);
    end
    
    % check for column name missing or index out of range
    if (isempty(col) || col > length(headers))
        error('scanfile:colname', 'Missing column: %s\n', specstr);
    end
    
    % process column
    if (mkind)
        [idx, key] = makeindicators(X{col},key,0);
        dat(:,i) = idx;
        if (i <= nout)
            varargout{i} = key;
        end
    else
        for j = 1:nrows
            entry = sscanf(X{col}{j}, '%f');
            if (isempty(entry))
                entry = NaN;
            end
            dat(j,i) = entry;
        end
    end
end
