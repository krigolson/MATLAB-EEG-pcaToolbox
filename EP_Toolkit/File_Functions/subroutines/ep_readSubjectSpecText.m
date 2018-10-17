function EPdata = ep_readSubjectSpecText(EPdata, specFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EPdata = ep_readSubjectSpecText(EPdata, specFileName) - reads in subject specs text file
% and replaces subject specs structure in EPdata with the new information.  Also subject names.
%
%Inputs
%   EPdata         : Structured array with the data and accompanying information.  See readData.
%	specFileName: file name for subject specs text file, including path and suffix.
%
%Outputs
%   EPdata         : Structured array with the data and accompanying information.  See readData.
%
% History:
%
% by Joseph Dien (11/16)
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
            disp(['Error: Addition of subject spec information failed.']);
        else
            numSpecs=length(theData)-1; %number of specs, not including the subject numbers
            specNames=cell(numSpecs,1);
            for iSpec=1:numSpecs
                specNames{iSpec}=theData{iSpec+1}{1};
            end;
            if (length(theData{1})-1) ~= numSubs
                disp(['Error: Number of subjects in the subject specs text file is incorrect.']);
            else
                subNames=cell(numSubs,1);
                for iSubject=1:numSubs
                    subNames{iSubject}=theData{1}{iSubject+1};
                end;
                subjectSpecs=cell(numSubs,numSpecs);
                for iSpec=1:numSpecs
                    for iSubject=1:min(numSpecs,length(theData{iSpec+1})-1)
                        subjectSpecs{iSubject,iSpec}=theData{iSpec+1}{iSubject+1};
                    end;
                end;
                disp('Replacing current contents of subject specs structure with contents of the subject specs text file.');
                EPdata.subjectSpecs=subjectSpecs;
                EPdata.subjectSpecNames=specNames;
                EPdata.subNames=subNames;
            end;
        end;
    end;
end;

