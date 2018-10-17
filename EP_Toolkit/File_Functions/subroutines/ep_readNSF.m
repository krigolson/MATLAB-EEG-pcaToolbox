function [theData,theEvents]=ep_readNSF(fileName);
% ep_readNSF - [theData,theEvents]=ep_readNSF(fileName); -
% Reads the contents of an EGI .nsf data file.
%
%Input:
%	fileName  : The name of the file, including the suffix.
%
%Output:
%  theData : The cell array with the data.
%  theEvents : The cell array with the events.

%History
%  by Joseph Dien (2/26/11)
%  jdien07@mac.com
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

theData=[];
theEvents=[];

try
    load('-mat', fileName);
catch
    msg{1}=['The attempt to load in the file ' fileName 'resulted in the error:' lasterr];
    [msg]=ep_errorMsg(msg);
    return
end;

theVars=whos;

for i=1:length(theVars)
    if strcmp('ECI_TCPIP_',theVars(i).name(1:min(10,length(theVars(i).name))))
        if isempty(theEvents)
            eval(['theEvents=' theVars(i).name ';']);
        else
            msg{1}='There were multiple event arrays in the file.';
            [msg]=ep_errorMsg(msg);
            return
        end;
    else
        if ~strcmp('Impedances_',theVars(i).name(1:min(11,length(theVars(i).name)))) && ~strcmp('fileName',theVars(i).name) && ~strcmp('theData',theVars(i).name) && ~strcmp('theEvents',theVars(i).name)
            if isempty(theData)
                eval(['theData=' theVars(i).name ';']);
            else
                msg{1}='There were multiple data arrays in the file.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
    end;
end;