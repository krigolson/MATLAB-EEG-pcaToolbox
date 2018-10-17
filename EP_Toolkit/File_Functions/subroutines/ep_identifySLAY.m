function [errInFile, montage]=ep_identifySLAY(infile);
%function [errInFile, montage]=ep_identifySLAY(infile);
%identifies SLAY resource in EGIS files
%
%Inputs
%  infile  : The EGIS file with the SLAY resource.
%Outputs
%  errInFile: Trouble reading in SLAY montage.
%  montage  : The electrode montage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%History
%  by Joseph Dien (11/6/08)
%  jdien07@mac.com
%
% bugfix 1/28/13 JD
% Fixed identifying all EGIS files as Hydrocel-128P under OS X 10.8.

[errInFile theSLAY]=system(['osascript ''' which('readSLAY.scpt') ''' ''' infile ''' ' num2str(128)]);
montage=theSLAY;
SLAYoffset=findstr(theSLAY,'SLAY');

thisFile=which('ep_identifySLAY');
[pathstr, name, ext] = fileparts(thisFile);

theID=0;
if errInFile ~= 1
    for i=1:29
        [err out{i}]=system(['osascript ''' which('readSLAY.scpt') ''' ''' [pathstr '/SLAY.rtf'] ''' ' num2str(i+127)]);
        TAGoffset=findstr(out{i},'SLAY');
        tag=out{i}(1+TAGoffset:50+TAGoffset);
        if tag==theSLAY(1+SLAYoffset:50+SLAYoffset)
            if theID
                theID=0;
                disp('Error: multiple matches for SLAY resource.');
                break
            else
                theID=i;
            end;
        end;
    end;
    
    switch theID
        case 0
            disp('Could not identify SLAY resource.');
            montage=[];
 %up through NS 4.1
        case 1
            montage='GSN200-128-21';
        case 2
            montage='GSN200-128-1';
        case 3
            montage='GSN200-128-R1';
        case 4
            montage='GSN200-64-A1';
        case 5
            montage='GSN200-64-I1';
        case 6
            montage='GSN200-64-T1';
        case 7
            montage='GSN200-64-A2';
        case 8
            montage='preGSN-64';
        case 9
            montage='GSN200-256-2';
        case 10
            montage='GSN200-256-21';
        case 11
            montage='Hydrocel-64-1'; %4.3 version
        case 12
            montage='Hydrocel-128-1'; %4.3 version?
        case 13
            montage='Hydrocel-256-1'; %4.3 version
 %after NS 4.1
        case 14
            montage='GSN200-128-21';
        case 15
            montage='GSN200-128-1';
        case 16
            montage='GSN200-128-R1';
        case 17
            montage='GSN200-64-A1';
        case 18
            montage='GSN200-64-I1';
        case 19
            montage='GSN200-64-T1';
        case 20
            montage='GSN200-64-A2';
        case 21
            montage='preGSN-64';
        case 22
            montage='GSN200-256';
        case 23
            montage='Hydrocel-64-1';
        case 24
            montage='Hydrocel-128-1';
        case 25
            montage='Hydrocel-256-1';
        case 26
            montage='GSN200-128-21';
        case 27
            montage='Hydrocel-128-1'; %4.1 version
        case 28
            montage='Hydrocel-256-1'; %another 4.3 version?
        case 29
            montage='Hydrocel-128-P'; %Premie
    end;
else
    if strcmp(montage(end-6:end-1),'(-192)')
        [montage]=ep_askForMontage;
        errInFile=0;
    end;
end;