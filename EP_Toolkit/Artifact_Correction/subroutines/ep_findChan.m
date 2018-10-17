function [theChan theOrder] = ep_findChan(eloc, badChans, targetLocs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [theChan theOrder] = ep_findChan(eloc, badChans)
%
%	Given the electrode locations, figures out which non-bad channel is closest to the desired channel coordinates.
%
%Inputs
%  badChans:   list of bad channels.
%  eloc          : structure containing the channel names and locations.
%                 Assumes only electrodes present, no fiducials.
%   .X           : x-coordinate  
%   .Y           : y-coordinate 
%   .Z           : z-coordinate 
%
%Outputs
%	theChan      : The channel closest to Cz without being a bad channel.  Zero if they are all bad.
%   theOrder     : This channel is which closest to Cz (e.g., 3 means 3rd, and 1st and 2nd were therefore bad channels and so skipped).
%
% History:
%
% by Joseph Dien (9/7/16)
% jdien07@mac.com
%
% modified 12/18/17 JD
% Generalized function to find any given channel by coordinates.  Changed name of function as well.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Cz [0 -9.7989 107.9359]

nChan=length(eloc);
x=zeros(1,nChan);
y=zeros(1,nChan);
z=zeros(1,nChan);

for i=1:nChan
    if ~isempty(eloc(i).X) && ~isempty(eloc(i).Y) && ~isempty(eloc(i).Z)
        x(i)=eloc(i).X;
        y(i)=eloc(i).Y;
        z(i)=eloc(i).Z;
    else
        x(i)=NaN;
        y(i)=NaN;
        z(i)=NaN;
    end;
end;

if max(abs(x)) <= 1
    MNItransform=[ 0 -15 0 0.08 0 -1.571 102 93 100 ]; %assume eloc coordinates are from a .elp file
else
    MNItransform=[ 0 -15 4 0.05 0 -1.571 10.2 12 12.2 ]; %assume eloc coordinates are from a .sfp file
end;

MNIcoord(:,1) = x';
MNIcoord(:,2) = y';
MNIcoord(:,3) = z';

transmat=traditionaldipfit(MNItransform);
MNIcoord=transmat*[MNIcoord ones(size(MNIcoord,1),1)]';
MNIcoord=MNIcoord(1:3,:)';

x=MNIcoord(:,1)';
y=MNIcoord(:,2)';
z=MNIcoord(:,3)';

elecDistances = zeros(nChan,1);
for chan = 1:nChan
    if ~isnan(x(chan)) && ~isnan(y(chan)) && ~isnan(z(chan))
        elecDistances(chan)=sqrt((x(chan)-targetLocs(1))^2+(y(chan)-targetLocs(2))^2+(z(chan)-targetLocs(3))^2);
    else
        elecDistances(chan)=inf;
    end;
end;
[B IX]=sort(elecDistances);

theOrder=1;
while ismember(IX(theOrder),badChans)
    theOrder=theOrder+1;
end;
if theOrder > nChan
    theChan=0;
else
    theChan=IX(theOrder);
end;