function EPdataset=ep_saveEPdataset(EPdataset,EPdata,entryNumber,fileSaved);
%  EPdataset=ep_saveEPdataset(EPdataset,EPdata,entryNumber,fileSaved);
%       Adds a new dataset to the work directory.
%
%Inputs:
%  EPdataset      : Structured array with the list of files in the work directory
%     .EPwork     : The path of the work directory.
%     .dataName   : The list of dataset names.
%  EPdata         : Structured array with the data and accompanying information.  See readData.
%  entryNumber    : Where to enter the dataset in the dataset list.  Defaults to the end of the list.
%  saved          : Has the latest version of the file been saved?  Defaults to yes.
%
%Outputs:
%  EPdataset        : Structured array with the list of files in the work directory

%History:
%  by Joseph Dien (7/29/09)
%  jdien07@mac.com
%
%  bugfix 8/3/09 JD
%  Was dropping suffixes even though readData had already dropped the file suffix, resulting in truncation of file names
%  with multiple periods in the name.  Also, file names with multiple periods would confuse Matlab about the file type,
%  resulting in it saving the file as a ASCII file rather than as a .mat file.
%
%  modified 10/31/09 JD
%  Added more information to EPdataset to speed up main pane and preprocess pane.
%
% modified 6/15/10 JD
% Marks when a file isn't saved yet.
%
% bugfix 3/24/14 JD
% Fixed when loading in a new file that had the same name as multiple existing files, appending dashed number to prior dashed
% number instead of replacing it (e.g., name-1-2)
%
% modified 5/1/18 JD
% If trying to save a file and preferences are set to v6 or v7 and it is over 2GB, then instead of not saving it will now save in v7.3 format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    entryNumber=length(EPdataset.dataset)+1;
end;

if nargin<4
    fileSaved='yes';
end;

if length(entryNumber) > 1
    msg{1}='Multiple datasets were specified for deletion.  Only one at a time is allowed.';
    [msg]=ep_errorMsg(msg);
    return
end;

if isempty(EPdata.dataName)
    EPdata.dataName='PCA';
end;
sameName=1;
suffix=0;
baseName=EPdata.dataName;
tempName=baseName;
while sameName
    sameName=0;
    for i=1:length(EPdataset.dataset)
        if strcmp(EPdataset.dataset(i).dataName,tempName)
            sameName=1;
        end;
    end;
    if sameName
        suffix=suffix+1;
        tempName=[baseName '-' num2str(suffix)];
    end;
end;
EPdata.dataName=tempName;

newDataset=ep_addToEPworkCache(EPdata);

if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
    msg{1}='The work directory cannot be found.';
    [msg]=ep_errorMsg(msg);
    return
end;

EPdataset=ep_checkEPworkCache(EPdataset);

EPdataset.dataset=[EPdataset.dataset(1:entryNumber-1) newDataset EPdataset.dataset(entryNumber:end)];

EPdataset.dataset(entryNumber).saved = fileSaved;

eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset'' EPdataset']);

warning('');
warning('off','MATLAB:save:sizeTooBigForMATFile')
eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep EPdata.dataName '.mat'' EPdata;']);
warning('on')
if ~isempty(lastwarn)
    [msg,msgID] = lastwarn;
    if strcmp(msgID,'MATLAB:save:sizeTooBigForMATFile')
        save('-mat', [EPdataset.EPwork filesep 'EPwork' filesep EPdata.dataName '.mat'], 'EPdata','-v7.3');
    end;
end;



