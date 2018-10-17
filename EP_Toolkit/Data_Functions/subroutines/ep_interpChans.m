function [sevenDdataOut]=ep_interpChans(sevenDdataIn, oldElocs, newElocs);
% [sevenDdataOut]=ep_interpChans(sevenDdataIn, oldElocs, newElocs);
% Given data and their electrode coordinates, will generate the interpolated EEG corresponding to a set of new electrode coordinates.
%
%Input:
%    sevenDdataIn   : input 7D data matrix [channels, time points, cells/trials, subjects, factors, freqs, relations]
%    oldElocs      : The electrode location information, one for each channel (see readlocs header)
%                    eloc is the same length as the channels, with REG channels having a blank entry.
%    newElocs      : The electrode location information for the new interpolated channels.
%
%Outputs:
%    sevenDdataOut      : output 7D data matrix [channels, time points, cells/trials, subjects, factors, freqs, relations]

%History
%  by Joseph Dien (8/18/15)
%  jdien07@mac.com
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
%

sevenDdataOut=[];

if size(sevenDdataIn,7) > 1
    disp('Error: Cannot interpolate relational data.');
    return;
end;

maxRad=0.5;
GRID_SCALE=67;

%old electrode coordinates
[y,x] = pol2cart(([oldElocs.theta]/360)*2*pi,[oldElocs.radius]);  % transform electrode locations from polar to cartesian coordinates
y=-y; %flip y-coordinate so that nose is upwards.
plotrad = min(1.0,max([oldElocs.radius])*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
x = x*(maxRad/plotrad);
y = y*(maxRad/plotrad);

xmin = min(-maxRad,min(x));
xmax = max(maxRad,max(x));
ymin = min(-maxRad,min(y));
ymax = max(maxRad,max(y));

x=round(((x/(xmax-xmin))*GRID_SCALE)+ceil(GRID_SCALE/2));
y=round(((y/(ymax-ymin))*GRID_SCALE)+ceil(GRID_SCALE/2));

%new electrode coordinates
[y1,x1] = pol2cart(([newElocs.theta]/360)*2*pi,[newElocs.radius]);  % transform electrode locations from polar to cartesian coordinates
y1=-y1; %flip y-coordinate so that nose is upwards.
plotrad = min(1.0,max([newElocs.radius])*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
x1 = x1*(maxRad/plotrad);
y1 = y1*(maxRad/plotrad);

xmin1 = min(-maxRad,min(x1));
xmax1 = max(maxRad,max(x1));
ymin1 = min(-maxRad,min(y1));
ymax1 = max(maxRad,max(y1));

x1=round(((x1/(xmax1-xmin1))*GRID_SCALE)+ceil(GRID_SCALE/2));
y1=round(((y1/(ymax1-ymin1))*GRID_SCALE)+ceil(GRID_SCALE/2));

fprintf('%60s\n',' ' );
for iChan=1:length(x1)
    for iFreq=1:size(sevenDdataIn,6)
        for iFac=1:size(sevenDdataIn,5)
            for iSub=1:size(sevenDdataIn,4)
                for iCell=1:size(sevenDdataIn,3)
                    for iPoint=1:size(sevenDdataIn,2)
                        fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('Remapping point %d in channel %4d of %4d',iPoint, iChan, length(x1)))
                        [Xi,Yi,Zi] = griddata(x,y,sevenDdataIn(:,iPoint,iCell,iSub,iFac,iFreq),[1:GRID_SCALE]',[1:GRID_SCALE],'v4');
                        %v4 interpolates to the edge of the box.  With other interpolation options, if a bad channel is at the edge of the montage then it
                        %would just be NaN.
                        sevenDdataOut(iChan,iPoint,iCell,iSub,iFac,iFreq)=Zi(y1(iChan),x1(iChan));
                    end;
                end;
            end;
        end;
    end;
end;
fprintf('%60s\n',' ' );
