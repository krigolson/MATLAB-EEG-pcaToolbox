function elecDistances=ep_closestChans(eloc);
%  elecDistances=closestChans(eloc);
%       computes distance from every electrode to every other one
%
%Inputs:
%  eloc          : structure containing the channel names and locations.
%                 Assumes only electrodes present, no fiducials.
%   .X           : x-coordinate  
%   .Y           : y-coordinate 
%   .Z           : z-coordinate 
%   
%Outputs:
%  elecDistances : distances to other electrodes (electrodes, electrodes)
%

%History:
%  by Joseph Dien (2/8/09)
%  jdien07@mac.com
%
%  bugfix 4/22/11 JD
%  Sets channels with missing coordinates to infinite.

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

elecDistances = zeros(nChan);

for chan1 = 1:nChan-1
    for chan2 = chan1+1:nChan
        if ~isnan(x(chan1)) && ~isnan(y(chan1)) && ~isnan(z(chan1)) && ~isnan(x(chan2)) && ~isnan(y(chan2)) && ~isnan(z(chan2))
            elecDistances(chan1,chan2)=sqrt((x(chan1)-x(chan2))^2+(y(chan1)-y(chan2))^2+(z(chan1)-z(chan2))^2);
            elecDistances(chan2,chan1)=elecDistances(chan1,chan2);
        else
            elecDistances(chan1,chan2)=inf;
        end;
    end;
end;