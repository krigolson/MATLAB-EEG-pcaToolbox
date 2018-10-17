function EPdata = ep_readEventText(EPdata, specFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EPdata = ep_readEventText(EPdata, specFileName) - reads in events text file
% and replaces events structure in EPdata with the new events information.
%
%	Reads in events text file for file formats that do not support such
%	information.
%
%Inputs
%   EPdata         : Structured array with the data and accompanying information.  See readData.
%	specFileName: file name for events text file, including path and suffix.
%
%Outputs
%   EPdata         : Structured array with the data and accompanying information.  See readData.
%
% History:
%
% by Joseph Dien (12/14)
% jdien07@mac.com
%
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

numSubs=length(EPdata.subNames);
numWaves=length(EPdata.cellNames);

if exist(specFileName,'file')
    specFid=fopen(specFileName);
    if specFid==-1
        disp(['Error: Unable to open the spec file: ' specFileName]);
    else
        tempVar=fgetl(specFid);
        delim='\t';
        if ~isempty(strfind(tempVar,',')) && isempty(strfind(tempVar,'\t'))
            delim=','; %if there are commas and no tabs, assume it is a comma-delimited file.
        elseif ~isempty(strfind(tempVar,' ')) && isempty(strfind(tempVar,'\t'))
            delim=' '; %if there are spaces and no tabs, assume it is a space-delimited file.
        end;
        
        numcols=length(regexp(tempVar,[delim '+']))+1; %determine number of columns based on number of delimiters (treating repeats as a single delimiter)
        if regexp(tempVar,[delim '$'])
            numcols=numcols-1; %if there is an extra tab at the end of the line, drop it.
        end;
        
        frewind(specFid);
        theData=textscan(specFid, repmat('%s',1,numcols),'Delimiter',delim);
        fclose(specFid);
        if isempty(theData)
            disp(['Error: Addition of event information failed.']);
        else
            numRows=length(theData{1})-1; %the first row are the spec labels
            numSpecs=length(theData);
            specNames={'subject';'cell';'type';'sample';'value';'duration'};
            keyNames={'code';'data';'datatype';'description'};
            numKeys=(numSpecs-length(specNames))/length(keyNames);
            if mod(numSpecs-length(specNames),length(keyNames)) ~=0
                disp(['Error: Number of spec fields in the events text file is incorrect.']);
            else
                badSpecNames=0;
                for iSpec=1:length(specNames)
                    if ~strcmp(theData{iSpec}{1},specNames{iSpec})
                        badSpecNames=iSpec;
                    end;
                end;
                if badSpecNames
                    disp(['Error: The spec name ' specNames{badSpecNames} ' in the events text file does not conform to the expected format.']);
                else
                    badSpecNames=0;
                    for iSpec=length(specNames)+1:4:numSpecs
                        if ~strcmp(theData{iSpec}{1}(1:length(keyNames{1})),keyNames{1})
                            badSpecNames=iSpec;
                        end;
                        if ~strcmp(theData{iSpec+1}{1}(1:length(keyNames{2})),keyNames{2})
                            badSpecNames=iSpec;
                        end;
                        if ~strcmp(theData{iSpec+2}{1}(1:length(keyNames{3})),keyNames{3})
                            badSpecNames=iSpec;
                        end;
                        if ~strcmp(theData{iSpec+3}{1}(1:length(keyNames{4})),keyNames{4})
                            badSpecNames=iSpec;
                        end;
                    end;
                    if badSpecNames
                        disp(['Error: The spec name ' specNames{badSpecNames} ' in the events text file does not conform to the expected format.']);
                    else
                        badSub=0;
                        badCell=0;
                        for iEvent=2:numRows+1
                            theSub=str2num(theData{1}{iEvent});
                            if isempty(theSub) || (theSub < 1) || (theSub > numSubs)
                                badSub=1;
                            end;
                            theCell=str2num(theData{2}{iEvent});
                            if isempty(theCell) || (theCell < 1) || (theCell > numWaves)
                                badCell=1;
                            end;
                        end;
                        if badSub || badCell
                            disp(['Error: The subject and/or cell numbers in the events text file have a problem.']);
                        else
                            newEvents=cell(numSubs,numWaves);
                            for iEvent=2:numRows+1
                                theSub=str2num(theData{1}{iEvent});
                                theCell=str2num(theData{2}{iEvent});
                                oneEvent=[];
                                oneEvent.type=theData{3}{iEvent};
                                oneEvent.sample=str2num(theData{4}{iEvent});
                                oneEvent.value=theData{5}{iEvent};
                                oneEvent.duration=str2num(theData{6}{iEvent});
                                if numKeys
                                    for iKey=1:numKeys
                                        theSpec=((iKey-1)*length(keyNames))+length(specNames)+1;
                                        if ~isempty(theData{theSpec}{iEvent})
                                            oneEvent.keys(iKey).code=theData{theSpec}{iEvent};
                                            oneEvent.keys(iKey).data=theData{theSpec+1}{iEvent};
                                            oneEvent.keys(iKey).datatype=theData{theSpec+2}{iEvent};
                                            oneEvent.keys(iKey).description=theData{theSpec+3}{iEvent};
                                        else
                                            oneEvent.keys(iKey)=struct('code','','data','','datatype','','description','');
                                        end;
                                    end;
                                else
                                    oneEvent.keys=struct('code','','data','','datatype','','description','');
                                end;
                                newEvents{theSub,theCell}(end+1)=oneEvent;
                            end;
                            events=newEvents;
                            disp('Replacing current contents of events structure with contents of the events text file.');
                        end;
                    end;
                end;
            end;
        end;
        EPdata.events=events;
    end;
end;

