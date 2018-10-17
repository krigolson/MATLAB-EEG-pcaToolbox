function time=ep_ms2time(ms);
%  time=ep_ms2time(ms);
%       Converts ms number to string expressing time in digital form
%
%Inputs:
%  ms   : milliseconds
%Outputs:
%  time : string output in digital format (hh:mm:ss:ms)

%History:
%  by Joseph Dien (1/19/14)
%  jdien07@mac.com
%
%
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

time='';

if ~isnumeric(ms)
    return
end;

hours=floor(ms/3600000);

ms=rem(ms,3600000);

minutes=floor(ms/60000);

ms=rem(ms,60000);

seconds=floor(ms/1000);

ms=rem(ms,1000);

time=[sprintf('%02.0f',hours) ':' sprintf('%02.0f',minutes) ':' sprintf('%02.0f',seconds) ':' sprintf('%03.0f',ms)];