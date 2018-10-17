function [err out]=ep_addSLAY(outfile,montage)
%function [err out]=ep_addSLAY(SLAYfile,outfile,montage)
%adds SLAY resource to EGIS files
%
%Inputs
%  SLAYfile	: The file with the all the SLAY resources stored in it.
%  outfile  : The EGIS file to have the SLAY added to it.
%  montage  : The electrode montage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%History
%  by Joseph Dien (11/6/08)
%  jdien07@mac.com

resourceID = 0;

switch montage
    case 'GSN200-128-21'
        resourceID = 128;
    case 'GSN200-128-1'
        resourceID = 129;
    case 'GSN200-128-R1'
        resourceID = 130;
    case 'GSN200-64-A1'
        resourceID = 131;
    case 'GSN200-64-I1'
        resourceID = 132;
    case 'GSN200-64-T1'
        resourceID = 133;
    case 'GSN200-64-A2'
        resourceID = 134;
    case 'preGSN-64'
        resourceID = 135;
    case 'GSN200-256-2'
        resourceID = 136;
    case 'GSN200-256-21'
        resourceID = 137;
    case 'Hydrocel-64-1'
        resourceID = 138;
    case 'Hydrocel-128-1'
        resourceID = 139;
    case 'Hydrocel-256-1'
        resourceID = 140;
end;

if resourceID == 0
    err=1;
    out='Montage unknown.';
else
    SLAYfile = 'SLAY.rtf'; %The file with the all the SLAY resources stored in it.
    thisFile=which('ep_addSLAY');
    [pathstr, name, ext] = fileparts(thisFile); %find the folder that has all the SLAY related function files.
    
    [err out]=system(['osascript ''' which('addSLAY.scpt') ''' ''' [pathstr '/' SLAYfile] ''' ''' outfile ''' ' num2str(resourceID) ' ' montage]);
end;    