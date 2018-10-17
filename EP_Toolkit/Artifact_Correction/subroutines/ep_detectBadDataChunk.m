function [outputLog, globalBadChans, totBadChanNum, totBadTrialNum] = ep_detectBadDataChunk(inFile, startChunk, endChunk, badDataCriteria, globalBadChans, editMode, theSubject);
% [outputLog, globalBadChans, totBadChanNum, totBadTrialNum] = ep_detectBadDataChunk(inFile, startChunk, endChunk, badDataCriteria, globalBadChans, editMode, theSubject);
%	Detects bad channels and bad trials.  Bad channels have a
%	minimum/maximum difference over threshold.  Also, channels whose maximum difference from all other electrodes is
%   no smaller than a certain threshold is declared
%   bad.  Trials with too many bad channels (including globally bad channels) declared bad trials.  Also,
%   there is an option for trials with bad channels next to another bad channel (but not globally bad channels) declared bad trials
%   because they usually represent trials with generalized artifacts rather than isolated bad channels.
%   Not enabled by default as it generally turns out to be too stringent.
%   Channels that are bad on too many of the good trials are declared globally bad.
%   Criteria settings determine whether each of these considerations are applied and to what degree.
%   Bad trials (that is cells) not registered for average files.
%
%Inputs
%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.
%	startChunk: starting chunk (usually 1)
%   endChunk:   ending chunk
%   badDataCriteria:  Criteria for detecting bad data.
%       .window:    moving average window for smoothing
%       .minmax:    difference from minimum to maximum for bad channel
%       .trialminmax:  difference from minimum to maximum for bad trial
%       .badnum:    percent of bad channels exceeded to declare bad trial, rounding down
%       .hminmax:   difference from minimum to maximum for bad horizontal EOG
%       .neighbors: number of electrodes considered to be neighbors
%       .badchan:   maximum microvolt difference allowed from best matching neighbor
%       .maxneighbor:   maximum microvolt difference allowed from best matching neighbor
%       .blink:     threshold correlation with blink template, 0 to 1
%       .detrend:   1 to detrend
%       .badtrials: percentage of good trials chan is bad to declare a channel globally bad
%       .replace:   1 to interpolate bad channels from neighbors.
%       .noadjacent:1 to not allow adjacent bad channels (trial or subject declared bad)
%   globalBadChans:   list of globally bad channels.  Will be set to a flat line.
%   editMode  : How to identify artifacts (automatic: use automatic criteria and enter marks into file; manual: use
%               existing marks in file; both: use existing marks in file and add additional ones based on automatic criteria)
%               (default: both)
%   theSubject: which subject of the file is being processed.
%
%   The input chunks include: dataChunk in EP data format.
%
%Outputs
%   outputLog: output messages from bad data process.
%   globalBadChans:   updated list of globally bad channels.
%   totBadChanNum:  Total list of bad channels (subject,trial,channel).
%   totBadTrialNum: Total list of bad trials.
%
%   Updated output chunks: dataChunk in EP data format.
%
% History:
%
% by Joseph Dien (2/09)
% jdien07@mac.com
%
% modified 3/19/09 JD
% Changed to use EP format data to provide more flexibility with I/O functions.
%
% modified 4/17/09 JD
% Treats flat channels as bad data unless they are a reference channel which is flat over entire dataset.
% Dropped eloc as separate input parameter (now part of data).
%
% modified 5/28/09 JD
% no longer zeroes out bad channels and trials, just notes them for later correction.
%
% bugfix 6/11/09 JD
% Was always applying no adjacent constraint even when option was turned off.
%
% bugfix 6/22/09 JD
% if no baseline set then will not try to subtract it (avoiding divide by zero error warning)
% Fix to badChanNum not being in right format, causing file not to be written out.
%
% modified & bugfix 9/18/09 JD
% Added support for multiple refChans to deal with mean mastoid data where the presence of the two reference channels (correlated -1)
% was causing ICA problems.  Trial specs coded only for single_trial data.
% Commented out setting bad trial code for "edit" field as was causing crash for non-EGIS files and
% not being used in any case.
% Was identifying wrong channel as bad when there was a flat channel and there was a reference channel with a lower
% number.
%
% bugfix 11/21/09 JD
% Replaced "union" commands with "unique" commands because certain situations caused the "union" command to crash in
% Matlab 2007.
%
% bugfix 12/3/09 JD
% Fix for bug where channel 1 would be marked as bad in too many trials if there was only one globally bad channel.
%
% modified 2/7/10 JD
% Changed bad channel field to negative numbers for still bad channels.
%
% modified 2/11/10 JD
% Will now work with subject average files with multiple subjects.
%
% bugfix 10/12/10 JD
% Fixed bad chan and bad trial fields not being recorded for all but last chunk when data are being processed as multiple chunks.
%
%  modified 10/12/10 JD
%  For continuous files, data now divided into one second epochs and can be artifact rejected in an epochwise fashion
%  in same fashion as segmented data.
%
% bugfix 4/18/11 JD
% Number of bad channels per trial in artifact correction log calculated incorrectly.  Artifact correction function itself not affected. 
% 
% modified 1/27/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
%
% bugfix 5/17/12 JD
% Bad channel field equals NaN in average files.
%
% modified 1/28/13 JD
% Added options for using manual edit marks instead of or in addition to automatic marking for bad channel and trials.
% Removed baseline correction as now performed as part of detrending step.
% 
% modified 9/22/13 JD
% Restricted blink correction to EEG channels.
%
% bugfix 12/28/13 JD
% Fixed crash when there are multiple global and trialwise bad channels when running under a version of Matlab earlier than 2013a.
%
%  bufix 3/11/14 JD
%  Handles decimal sampling rates gracefully.
%
% bugfix 4/8/14 JD
% Fixed crash when there are multiple global and trialwise bad channels.  Apparently Mathworks changed something again.
%
% bugfix 9/1/16 JD
% Fixed channels with large offsets being erroneously tagged as being bad channels due to filter edge artifact.
%
% modified 9/2/16 JD
% Check for difference from inverse of all other channels too since channels in sparser arrays can end up being the only opposite polarity waveform.
% Baseline correct with first sample for each epoch when applying maxneighbor criterion so simple offset does not trigger it.
%
% modified 4/26/17 JD
% Adjusted "too many bad chans" criterion so not rounded off and also based on total EEG channels rather than all channels.
%
% bugfix 4/2/18 JD
% Fixed not outputing to log correct mean number of bad channels for average files.
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

msg=['Detecting bad data.'];
disp(msg);
outputLog=[];

totBadChanNum=[];
totBadTrialNum=[];
badTrialOrigCount=0;
badTrialTooManyBadChanCount=0;
badTrialBadNeighborsCount=0;
for iChunk = startChunk:endChunk
    eval(['load ''' deblank(inFile) '''-' num2str(iChunk) '.mat']);
    elecDistances=ep_closestChans(dataChunk.eloc);
    if ~strcmp(dataChunk.reference.type,'AVG')
        refChan= dataChunk.reference.current;
    else
        refChan=[];
    end;
    
    numSamples=length(dataChunk.timeNames);
    numChans=length(dataChunk.chanNames);
    EEGchans=find(strcmp('EEG',dataChunk.chanTypes));
    nonRefChans=setdiff(EEGchans,refChan);
    
    moveAverageWindow=floor(badDataCriteria.window/(1000/ceil(dataChunk.Fs))); %samples width of window
    if strcmp(dataChunk.dataType,'continuous')
        numTrials=floor(size(dataChunk.data,2)/ceil(dataChunk.Fs)); %excess time points are tacked onto final epoch
        trialSize = min(ceil(dataChunk.Fs),size(dataChunk.data,2)); %one second epochs
    else
        numTrials = length(dataChunk.cellNames);
    end;
    if strcmp(editMode,'both')
        badChanNum{iChunk}=dataChunk.analysis.badChans;
        badTrialNum{iChunk}=dataChunk.analysis.badTrials;
    elseif strcmp(editMode,'automatic')
        badChanNum{iChunk}=zeros(size(dataChunk.analysis.badChans));
        badTrialNum{iChunk}=zeros(size(dataChunk.analysis.badTrials));
    else
        error('mode not supported');
    end;
    badTrialOrigCount=badTrialOrigCount+sum(dataChunk.analysis.badTrials);
    
    badChanNum{iChunk}(theSubject,:,globalBadChans)=-1;
    numRegressors=min(badDataCriteria.neighbors,numChans-1);
    neighbors=zeros(numChans,numRegressors);
    for chan=1:numChans
        [E IX]=sort(elecDistances(chan,:));
        neighbors(chan,:)=IX(2:numRegressors+1);
    end;
    
    chanDiffs=zeros(numChans,1);
    for iTrial=1:numTrials
        if strcmp(dataChunk.dataType,'continuous')
            theData=dataChunk.data(:,(iTrial-1)*trialSize+1:iTrial*trialSize,1,theSubject);
            if iTrial == numTrials
                theData=dataChunk.data(:,(iTrial-1)*trialSize+1:end,1,theSubject); %excess time points are tacked onto final epoch
            end;
        else
            theData=squeeze(dataChunk.data(:,:,iTrial,theSubject));
        end;
        
        %moving average window smoothing
        if (badDataCriteria.window == 0) || (moveAverageWindow == 0)
            theFilteredData=theData;
        else
            theFilteredData=filter((1/moveAverageWindow)*ones(moveAverageWindow,1),1,[repmat(theData(:,1),1,moveAverageWindow) theData]')';
            theFilteredData=theFilteredData(:,moveAverageWindow+1:end);
        end;
        
        %channel bad if too much change over course of trial
        trialBadChans=intersect(find(abs(max(theFilteredData')-min(theFilteredData'))>=badDataCriteria.minmax)',EEGchans);
        %channel bad if flat and not the reference but reference must be flat over entire dataset
        trialBadChans=unique([trialBadChans; nonRefChans(find(~std(theFilteredData(nonRefChans,:)')))]);

        trialGoodChans=setdiff(EEGchans,[trialBadChans; globalBadChans]); %good channels in the trial
        
        if length(trialGoodChans) > 1
            %channel bad if too different from all other channels.  Check for difference from opposite polarity of all other channels too.
            for iChan=1:length(trialGoodChans)
                theChan=trialGoodChans(iChan);
                otherChans=setdiff(trialGoodChans,theChan);
                chanWave=theFilteredData(theChan,:)-repmat(theFilteredData(theChan,1),1,size(theData,2));
                chanDiffs1=min(max(abs(theFilteredData(otherChans,:)-repmat(chanWave,length(trialGoodChans)-1,1)-repmat(theFilteredData(otherChans,1),1,size(theData,2))),[],2));
                chanDiffs2=min(max(abs(theFilteredData(otherChans,:)+repmat(chanWave,length(trialGoodChans)-1,1)-repmat(theFilteredData(otherChans,1),1,size(theData,2))),[],2));
                chanDiffs(theChan)=min([chanDiffs1 chanDiffs2]);
            end;
            trialBadChans=unique([find(chanDiffs >= badDataCriteria.maxneighbor); trialBadChans]); %bad channels due to being different from all the other channels
        end;
        
        if strcmp(dataChunk.dataType,'average')
            badChanNum{iChunk}(theSubject,iTrial,trialBadChans)=NaN;
        else
            badChanNum{iChunk}(theSubject,iTrial,trialBadChans)=-1;
        end;
        
        %too many bad channels?
        badTrial=0;
        if badDataCriteria.badnum ~=0
            if length(unique([trialBadChans; globalBadChans]))>=(length(EEGchans)*(badDataCriteria.badnum/100))
                theData=0;
                badTrialNum{iChunk}(theSubject,iTrial)=1;
                badTrial=1;
                badTrialTooManyBadChanCount=badTrialTooManyBadChanCount+1;
            end;
        end;
        
        %if there are adjacent trial bad channels, trial is bad
        if badDataCriteria.neighbors && badDataCriteria.noadjacent && ~badTrial && ~isempty(trialBadChans) %neighboring bad channels?
            badNeighbors=0;
            for i=1:length(trialBadChans)
                chan=trialBadChans(i);
                if ~isempty(intersect(trialBadChans,neighbors(chan,:)))
                    badNeighbors=1;
                    break
                end;
            end;
            if badNeighbors
                badTrialNum{iChunk}(theSubject,iTrial)=1;
                badTrialBadNeighborsCount=badTrialBadNeighborsCount+1;
            end;
        end;
        
    end;
    totBadChanNum=[totBadChanNum badChanNum{iChunk}];
    totBadTrialNum=[totBadTrialNum badTrialNum{iChunk}];
    
    dataChunk.analysis.badChans=badChanNum{iChunk};
    if ~strcmp(dataChunk.dataType,'average')
        dataChunk.analysis.badTrials=badTrialNum{iChunk};
    end;
    eval (['save ''' deblank(inFile) '''-' num2str(iChunk) '.mat dataChunk;']);
end;


%are any channels bad in too many of the good trials?
if badDataCriteria.badtrials ~=0
    numTotTrials=length(totBadTrialNum);
    numTotGoodTrials=sum(totBadTrialNum(theSubject,:)==0);
    tooBadChans=find((sum(totBadChanNum(theSubject,find(totBadTrialNum(theSubject,:)==0),:),2)/numTotGoodTrials)>badDataCriteria.badtrials/100);
    tooBadChans=tooBadChans(:);
    
    totBadChanNum=[];
    for iChunk = startChunk:endChunk
        eval(['load ''' deblank(inFile) '''-' num2str(iChunk) '.mat']);
        badChanNum{iChunk}(theSubject,:,tooBadChans)=1;
        totBadChanNum=[totBadChanNum badChanNum{iChunk}];
    end;
end;

globalBadChans=unique([globalBadChans; tooBadChans]);

numBadTrials=sum(totBadTrialNum(theSubject,:)');
if strcmp(dataChunk.dataType,'average')
    meanBadChans=mean(sum(isnan(totBadChanNum(theSubject,find(dataChunk.avgNum(theSubject,:)~=-1),:)),3)); %calculate mean bad chans only for good trials
else
    meanBadChans=-mean(sum(totBadChanNum(theSubject,find(totBadTrialNum(theSubject,:)==0),:),3)); %calculate mean bad chans only for good trials
end;

if strcmp(dataChunk.dataType,'continuous')
    theSegment = 'one second epoch';
else
    theSegment = 'trial';
end;

if numBadTrials ==1
    outputLog{1}=['There was 1 bad ' theSegment '.'];
else
    outputLog{1}=['There were ' num2str(numBadTrials) ' bad ' theSegment 's.'];
end;
if numBadTrials > 0
    outputLog{end+1}=['   Originally there were ' num2str(badTrialOrigCount) ' bad ' theSegment 's.'];
    outputLog{end+1}=['   There were ' num2str(badTrialTooManyBadChanCount) ' bad ' theSegment 's due to too many bad channels.'];
    outputLog{end+1}=['   There were ' num2str(badTrialBadNeighborsCount) ' bad ' theSegment 's due to neighboring bad channels.'];
end;
outputLog{end+1}=['For good ' theSegment 's, there was an average of ' num2str(meanBadChans) ' bad channels per ' theSegment '.'];
