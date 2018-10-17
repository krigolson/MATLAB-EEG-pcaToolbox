function cmp=ep_compareStrings(string1,string2)
% cmp=ep_compareStrings(string1,string2) -
% Compares two strings to determine if one is larger than the other.
% If both strings are numbers, then will be converted to numbers
% and compared numerically rather than alphabetically.
% Leading and trailing blanks are dropped.
%
%Input:
%  string1        : A string.
%  string2        : A string.
%
%Output:
%   cmp           : -1 if string1 is smaller than string2, 0 if they are
%                   the same, and 1 if string1 is larger.
%History
%  by Joseph Dien (8/14/15)
%  jdien07@mac.com
%  based in part on code posted at: http://www.mathworks.com/matlabcentral/answers/39374-compare-two-strings-based-on-ascii-dictionary-order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     Copyright (C) 1999-2018  Joseph Dien
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


if ~ischar(string1) || ~ischar(string2)
    cmp=[];
    return;
end;

string1=strtrim(string1);
string2=strtrim(string2);

if length(string1)==0 || length(string2)==0
    cmp=[];
    return;
end;

num1=str2num(string1);
num2=str2num(string2);

if ~isempty(num1) && ~isempty(num2)
    if num1==num2
        cmp=0;
    elseif num1<num2
        cmp=-1;
    else
        cmp=1;
    end;
    return
end;

a=string1;
b=string2;

%code from "Geoff" at matlabcentral starts

% Force the strings to equal length
x = char({a;b});
% Subtract one from the other
d = x(1,:) - x(2,:);
% Remove zero entries
d(~d) = [];
if isempty(d)
    cmp = 0;
else
    cmp = d(1);
end

cmp=sign(cmp);

%code from "Geoff" at matlabcentral ends

if cmp==0
    if length(string1) > length(string2)
        cmp=1;
    elseif length(string1) < length(string2)
        cmp=-1;
    end;
end;
