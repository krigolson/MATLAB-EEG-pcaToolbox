function [montage]=ep_askForMontage;    
%  [montage]=ep_askForMontage; 
%       ask for montage.
%
%Outputs:
%  montage      : The electrode montage.

%History:
%  by Joseph Dien (11/5/08)
%  jdien07@mac.com

theSelection = menu('Choose a montage',...
        'Adult GSN200 128-channel 2.1',...
        'Adult GSN 128-channel 1.0',...
        'Adult Recumbent GSN 128-channel 1.0',...
        'Adult GSN 64-channel 1.0',...
        'Infant GSN 64-channel 1.0',...
        'Toddler GSN 64-channel 1.0',...
        'Adult GSN 64-channel 2.0',...
        'Adult preGSN 64-channel 1.0',...
        'Adult GSN 256-channel 2.0',...
        'Adult GSN 256-channel 2.1',...
        'Adult Hydrocel 64-channel 1.0',...
        'Adult Hydrocel 128-channel 1.0',...
        'Adult Hydrocel 256-channel 1.0',...
        'None of the above.');

    drawnow %make sure the menu goes away immediately after the button is pressed.
    
    switch theSelection
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
            montage='Hydrocel-64-1';
        case 12
            montage='Hydrocel-128-1';
        case 13
            montage='Hydrocel-256-1';
        case 14
            montage='';
    end;