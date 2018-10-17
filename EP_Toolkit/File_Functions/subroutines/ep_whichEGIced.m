function [cedFile]=ep_whichEGIced(montage);
%function [cedFile]=ep_whichEGIced(montage, nChans);
%identifies which .ced electrode coordinate file matches an EGI montage
%
%Inputs
%  montage  : The electrode montage.
%
%Outputs
%  cedFile  : The name of the matching .ced file to use.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%History
%  by Joseph Dien (2/7/09)
%  jdien07@mac.com
%
%  modified (4/17/09) JD
%  Dropped channel input.  It is now assumed that if there is one more channel in the .ced file than there is in the
%  dataset the last one is an implicit reference.

cedFile='';

switch montage
    case 'GSN200-128-21'
        cedFile='GSN129.ced';
    case 'GSN200-128-1'
        cedFile='GSN129.ced';
    case 'GSN200-128-R1'
        cedFile='GSN129.ced';
    case 'GSN200-64-A1'
        cedFile='GSN65v1_0.ced';
    case 'GSN200-64-I1'
        cedFile='GSN65v1_0.ced';
    case 'GSN200-64-T1'
        cedFile='GSN65v1_0.ced';
    case 'GSN200-64-A2'
        cedFile='GSN65v2_0.ced';
    case 'preGSN-64'
        cedFile='preGSN65.ced';
    case 'GSN200-256-2'
        cedFile='GSN257.ced';
    case 'GSN200-256-21'
        cedFile='GSN257.ced';
    case {'Hydrocel-64-1','HydroCel GSN 64 1.0'}
        cedFile='GSN-Hydrocel-65 1.0.ced';
    case {'Hydrocel-128-1','HydroCel GSN 128 1.0'}
        cedFile='GSN-Hydrocel-129.ced';
    case {'Hydrocel-256-1','HydroCel GSN 256 1.0'}
        cedFile='GSN-Hydrocel-257.ced';
    otherwise
        cedFile='';
end;

