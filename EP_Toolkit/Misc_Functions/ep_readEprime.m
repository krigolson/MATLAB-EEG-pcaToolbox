function [experimentFieldNames, experimentData, fieldNames, trialData]=ep_readEprime(fname, textFormat);
%function [experimentFieldNames, experimentData, fieldNames, trialData]=ep_readEprime(fname);
%Reads in an E-Prime text file.
%
%Inputs
%  fname                  : The text file with the E-Prime output in it. 
%  textFormat             : The type of text file (defaults to 'UTF-16').
%Outputs
%  experimentFieldNames   : The experiment information field names
%  experimentData         : The experiment information
%  fieldNames             : The trial information field names
%  trialData              : The trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%History
%  by Joseph Dien (11/4/08)
%  jdien07@mac.com
%
% bugfix and modified 10/1/09 JD
% Crash when more than one level to E-Prime text file.
% Can now handle trials with different field names.
%
% modified 2/8/18 JD
% Modified to handle different text file types by adding textFormat specification.
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

if nargin < 1
    fname=[];
end;

if nargin < 2
    textFormat='UTF-16';
end;

if isempty(fname)
    %[infid fname] = get_fid('r');
    [outFiles, pathname] = uigetfile('*.*');
    fname=[pathname filesep outFiles];
% else
%     infid = fopen(fname);
end

% if infid == -1
%     error('File error')
% end;

% C = textscan(infid, '%s','delimiter','\n');
% C = textscanu(fname,'UTF-16LE',9,13);
% fclose(infid);

[f,msg]=fopen(fname,'r','n',textFormat);
txt=fscanf(f,'%c');
C=textscan(txt,'%[^\n]','delimiter','\n');

if ~isempty(C{1}{size(C{1},1),:}) %make sure the final line is a blank.
    C{1}{size(C{1},1)+1,:}='';
end;

%Get experiment information
lineCounter=1;
experimentFieldNames=[];
experimentData=[];
while isempty(strfind(C{1}{lineCounter,:},'Experiment'))
    lineCounter=lineCounter+1;
    if lineCounter > size(C{1},1)
        error('ran out of file');
    end;
end;

while isempty(strfind(C{1}{lineCounter,:},'*** Header End ***'))
    thisLine=strtrim(C{1}{lineCounter,:});
    [thisName theData]=strtok(thisLine,': ');
    experimentFieldNames{size(experimentFieldNames,2)+1}=thisName;
    experimentData{size(experimentData,2)+1}=theData(3:end);
    lineCounter=lineCounter+1;
    if lineCounter > size(C{1},1)
        error('ran out of file');
    end;
end;

%Get trial field names
lineCounter=1;
fieldNames=[];
while isempty(strfind(C{1}{lineCounter,:},'*** LogFrame Start ***'))
    lineCounter=lineCounter+1;
    if lineCounter > size(C{1},1)
        error('ran out of file');
    end;
end;
trialStart=C{1}{lineCounter-1,:};

%Get data
lineCounter=1;
trialNum=0;
trialData=[];

while lineCounter < size(C{1},1)
    lineCounter=lineCounter+1;
    if ~isempty(strfind(C{1}{lineCounter,:},trialStart)) %found a new trial
        lineCounter=lineCounter+2;
        if lineCounter > size(C{1},1)
            break;
        end;
        trialNum=trialNum+1;
        while isempty(strfind(C{1}{lineCounter,:},'*** LogFrame End ***'))
            if lineCounter > size(C{1},1)
                break;
            end;
            thisLine=strtrim(C{1}{lineCounter,:});
            [thisName theData]=strtok(thisLine,': ');
            if isempty(find(strcmp(thisName,fieldNames))) %is it a new field name?
                fieldNames{end+1}=thisName;
            end;
            trialData{trialNum,find(strcmp(thisName,fieldNames))}=theData(3:end);
            lineCounter=lineCounter+1;
        end;
    end;
end;

