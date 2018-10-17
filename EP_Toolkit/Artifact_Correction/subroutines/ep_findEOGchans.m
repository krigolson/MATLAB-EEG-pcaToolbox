function eog=ep_findEOGchans(eloc);
%  eog=ep_findEOGchans(eloc);
%       Finds EOG channels from electrode coordinates from .ced file.
%
%Inputs:
%  eloc          : structure containing the channel names and locations.
%                 Assumes only electrodes present, no fiducials.
%   .X           : x-coordinate
%   .Y           : y-coordinate
%   .Z           : z-coordinate
%
%Outputs:
%   eog.LUVEOG: left upper vertical EOG channel
%   eog.RUVEOG: right upper vertical EOG channel
%   eog.LLVEOG: left lower vertical EOG channel
%   eog.RLVEOG: right lower vertical EOG channel
%   eog.LHEOG: left horizontal EOG channel
%   eog.RHEOG: right horizontal EOG channel
%

%History:
%  by Joseph Dien (2/11/09)
%  jdien07@mac.com
%
%  bugfix 4/22/11 JD
%  Sets channels with missing coordinates to infinite to avoid crashes.

%standard EOG locations based on EGI 128 channel GSN200 net
% standard.LUVEOG=[9.2995, 4.0822, -1.7159];
% standard.RUVEOG=[9.2995, -4.0822, -1.7159];
% standard.LLVEOG=[7.1654, 3.3413, -8.0735];
% standard.RLVEOG=[7.1654, -3.3413, -8.0735];
% standard.LHEOG=[5.5373, 6.0007, -5.7386];
% standard.RHEOG=[5.5373, -6.0007, -5.7386];

standard.LUVEOG=[-49.0057183853233,78.6801264317438,-21.6480883280076;];
standard.RUVEOG=[48.9670795825849,78.7000559279483,-21.6490856340447;];
standard.LLVEOG=[-40.1104850480940,53.0647992321382,-98.0260280902464;];
standard.RLVEOG=[40.0807132886293,53.0811116190050,-98.0268443899529;];
standard.LHEOG=[-72.0199020626803,37.8961374881195,-68.7455379891984;];
standard.RHEOG=[71.9968949502115,37.9254331933456,-68.7470039963360;];

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
        elecDistances(chan)=sqrt((x(chan)-standard.LUVEOG(1))^2+(y(chan)-standard.LUVEOG(2))^2+(z(chan)-standard.LUVEOG(3))^2);
    else
        elecDistances(chan)=inf;
    end;
end;
[B IX]=sort(elecDistances);
eog.LUVEOG=IX(1);

elecDistances = zeros(nChan,1);
for chan = 1:nChan
    if ~isnan(x(chan)) && ~isnan(y(chan)) && ~isnan(z(chan))
        elecDistances(chan)=sqrt((x(chan)-standard.RUVEOG(1))^2+(y(chan)-standard.RUVEOG(2))^2+(z(chan)-standard.RUVEOG(3))^2);
    else
        elecDistances(chan)=inf;
    end;
end;
[B IX]=sort(elecDistances);
eog.RUVEOG=IX(1);

elecDistances = zeros(nChan,1);
for chan = 1:nChan
    if ~isnan(x(chan)) && ~isnan(y(chan)) && ~isnan(z(chan))
        elecDistances(chan)=sqrt((x(chan)-standard.LLVEOG(1))^2+(y(chan)-standard.LLVEOG(2))^2+(z(chan)-standard.LLVEOG(3))^2);
    else
        elecDistances(chan)=inf;
    end;
end;
[B IX]=sort(elecDistances);
eog.LLVEOG=IX(1);

elecDistances = zeros(nChan,1);
for chan = 1:nChan
    if ~isnan(x(chan)) && ~isnan(y(chan)) && ~isnan(z(chan))
        elecDistances(chan)=sqrt((x(chan)-standard.RLVEOG(1))^2+(y(chan)-standard.RLVEOG(2))^2+(z(chan)-standard.RLVEOG(3))^2);
    else
        elecDistances(chan)=inf;
    end;
end;
[B IX]=sort(elecDistances);
eog.RLVEOG=IX(1);

elecDistances = zeros(nChan,1);
for chan = 1:nChan
    if ~isnan(x(chan)) && ~isnan(y(chan)) && ~isnan(z(chan))
        elecDistances(chan)=sqrt((x(chan)-standard.LHEOG(1))^2+(y(chan)-standard.LHEOG(2))^2+(z(chan)-standard.LHEOG(3))^2);
    else
        elecDistances(chan)=inf;
    end;
end;
[B IX]=sort(elecDistances);
eog.LHEOG=IX(1);

elecDistances = zeros(nChan,1);
for chan = 1:nChan
    if ~isnan(x(chan)) && ~isnan(y(chan)) && ~isnan(z(chan))
        elecDistances(chan)=sqrt((x(chan)-standard.RHEOG(1))^2+(y(chan)-standard.RHEOG(2))^2+(z(chan)-standard.RHEOG(3))^2);
    else
        elecDistances(chan)=inf;
    end;
end;
[B IX]=sort(elecDistances);
eog.RHEOG=IX(1);

