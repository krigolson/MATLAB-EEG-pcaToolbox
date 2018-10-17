function [outputLog] = ep_markBadDataChunks(inFile, startChunk, endChunk, badChans, theSubject);
% ep_markBadDataChunks(inFile, startChunk, endChunk, badChans, theSubject);
%	Marks bad channels and trials with a flat line interrupted by a huge spike for another program like NetStation to fix.
%
%Inputs
%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.
%	startChunk: starting chunk (usually 1)
%   endChunk:   ending chunk
%   badChans:   list of globally bad channels.  Will be set to a flat line.
%   theSubject: which subject of the file is being processed.
%
%   The input chunks include: dataChunk in EP data format.
%
%Outputs
%   outputLog: output messages from bad data process.
%
%   Updated output chunks: dataChunk in EP data format.
%
% History:
%
% by Joseph Dien (2/09)
% jdien07@mac.com
%
% modified 3/14/09 JD
% Changed to use EP format data to provide more flexibility with I/O functions.
%
% modified 4/17/09 JD
% Dropped eloc as separate input parameter (now part of data).
%
% bugfix 12/8/09 JD
% Fixed crash when there are bad channels.  Thanks to Alex Lamey.
%
% modified 2/11/10 JD
% Will now work with subject average files with multiple subjects.
% Gets bad channel and bad trial info from the data chunk rather than from the function call.
%
% bufix 3/11/14 JD
% Handles decimal sampling rates gracefully.

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

msg='Marking bad channels and trials.';
disp(msg);
outputLog{1}=msg;

trialCount=1;
for chunk = startChunk:endChunk
    eval(['load ''' deblank(inFile) '''-' num2str(chunk) '.mat']);
    if strcmp(dataChunk.dataType,'continuous')
        numTrials=floor(size(dataChunk.data,2)/ceil(dataChunk.Fs)); %excess time points are tacked onto final epoch
        numSamples = min(ceil(dataChunk.Fs),size(dataChunk.data,2)); %one second epochs
    else
        numTrials = length(dataChunk.cellNames);
        numSamples=length(dataChunk.timeNames);
    end;
    numChans=length(dataChunk.chanNames);
    spike=[1000 zeros(1,numSamples-1)];

    badChanNum=dataChunk.analysis.badChans;
    badTrialNum=dataChunk.analysis.badTrials;

    for trial=1:numTrials
        if badTrialNum(trial+trialCount-1)
            if strcmp(dataChunk.dataType,'continuous')
                if trial == numTrials %excess time points are tacked onto final epoch
                    dataChunk.data(:,(trial-1)*numSamples+1:end,1,theSubject)=repmat([1000 zeros(1,size(dataChunk.data,2)-(numSamples*(numTrials-1))-1)],numChans,1);
                else
                    dataChunk.data(:,(trial-1)*numSamples+1:trial*numSamples,1,theSubject)=repmat(spike,numChans,1);
                end;
            else
                dataChunk.data(:,:,trial,theSubject)=repmat(spike,numChans,1);
            end;
        else
            if strcmp(dataChunk.dataType,'continuous')
                if trial == numTrials %excess time points are tacked onto final epoch
                    theData=dataChunk.data(:,(trial-1)*numSamples+1:end,1,theSubject);
                else
                    theData=dataChunk.data(:,(trial-1)*numSamples+1:trial*numSamples,1,theSubject);
                end;
            else
                theData=dataChunk.data(:,:,trial,theSubject);
            end;
            trialBadChans=find(badChanNum(1,trialCount+trial-1,:));
            trialGoodChans=setdiff([1:numChans],trialBadChans);
            theData(trialBadChans,:)=repmat(spike,length(trialBadChans),1);
            if strcmp(dataChunk.dataType,'continuous')
                if trial == numTrials %excess time points are tacked onto final epoch
                    dataChunk.data(:,(trial-1)*numSamples+1:end,1,theSubject)=theData;
                else
                    dataChunk.data(:,(trial-1)*numSamples+1:trial*numSamples,1,theSubject)=theData;
                end;
            else
                dataChunk.data(:,:,trial,theSubject)=theData;
            end;
        end;
    end;
    eval (['save ''' inFile '''-' num2str(chunk) '.mat dataChunk']);
    trialCount=trialCount+numTrials;
end;
