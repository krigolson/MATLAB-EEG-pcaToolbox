function cellNums=ep_segmentData(cellTable,inputFiles,importFormat,outputFormat,preview)
% cellNums=ep_segmentData(cellTable,inputFiles,preview) -
% segments continuous data files into single-trial segmented data.
%
%Input:
%  cellTable        : Cell array with segmenting parameters (name, prestim ms, poststim ms, delay, list of five specs with
%                       event spec (type, value, or key), relation, and value.  Empty when not used.
%  inputFiles       : List of session files including path.  If preview mode, then the EPdataset structured variable.
%  importFormat     : file format code for input data
%  outputFormat     : file format code for output data
%  preview          : indicates files should only be assessed for number of resulting epochs without
%                       generating output files.
%
%Output:
%   cellNums       : Array of number of epochs resulting for each cell (subject,cell).

%History
%  by Joseph Dien (11/16/13)
%  jdien07@mac.com
%
% bugfix 3/19/14 JD
% Fixed crash when segmenting or previewing.  Not sure why the syntax was working and suddenly was not.
%
% modified 6/18/14 JD
% Added starts, ends, contains, and -follows- keywords.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% modified 8/31/14 JD
% Added support for adding additional keys information to events in continuous data and trial specs in single-trial
% data.
%
% modified 9/4/14 JD
% Added delay field to segment function.
%
% modified 9/15/14 JD
% Added contents of type field to time-lock events list.
%
% modified 9/17/14 JD
% Adding event output for files that don't support events.
%
% modified 10/16/14 JD
% Passes screen size to ep_readData call.
%
% modified 10/30/14 JD
% Eliminated requirement that segmented epochs not overlap with each other.
%
% bugfix 8/14/15 JD
% Fixed criterion comparisons not being evaluated correctly for < and >.
%
% modified 8/17/15 JD
% Added ability to resegment (and thus reassign cells of) single-trial data.
%
% modified & bugfix 9/3/15 JD
% Adds next TRSP event after segmentation event without regard to whether
% the TRSP event fell within the segmentation period.
% The sectionTrial and SectionNumbers changed to start with 1 rather than 0.
% Added capability to resegment single-trial data to reassign cell
% conditions.
% Added capability to specify OR criteria for a condition by just giving
% them the same name.
% Trial spec names no longer continue to have TS- prefix when segmented
% data saved, which also resulted in trial specs not being recognized
% during segmentation.
%
% modified 12/15/15 JD
% When it runs into multiple matching specs, takes first one rather than dropping the segment.
%
% bugfix 12/18/15 JD
% Fixed crash when the criterion string is longer than the stimulus string for the "starts" and "ends" relations.
%
% modified 10/16/16 JD
% Added .stims field.
%
% modified 11/5/16 JD
% Added support for writing out subject spec text files.
%
% bugfix 11/8/16 JD
% Fixed segmentation ignoring the 6th criterion.
% Fixed section crit < and > relations being interpreted as the reverse direction.
%
% bugfix 3/16/17 JD
% Fixed error when batching multiple subject files and there are more than one row of criteria for the same cell name.
%
% bugfix 5/21/17 JD
% Fixed crash when resegmenting single-trial data.
% Fixed adding sec TrialSpec fields when already present, as in resegmenting data, resulting in crashes during averaging.
%
% modified & bugfix 6/21/17 JD
% Bad data edits passed on to segmented files.
% Flexible segments implemented.
% Fixed crash when segmenting time-frequency data.
% switch to amplitude scaling when adding freq data together.
% Fixed crash that can happen with continuous data.
%
% modified & bugfix 11/25/17 JD
% Added -precedes- crit option to the Segment function.
% Made fixes to -follows- function for relations other than = and ~=
% Added support for TS-EPoffset field.
%
% bugfix 3/1/18 JD
% Fixed unable to batch multiple files when the ced file contains BAD channels.
%
% bugfix 3/16/18 JD
% Now handles situation where the recording was aborted before the last trial spec was recorded.
%
% bugfix 4/29/18 JD
% No longer misses events when the .value field is a number rather than a string.
%
% bugfix 5/13/18 JD
% Now includes 'ep_segmentData' in .history field record.
%
% bugfix 6/14/18 JD
% Fixed using next available TRSP when a TRSP is missing instead of just dropping the trial.
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

global EPmain

cellNums=[];
numSpecs=6;
numOutCells=size(cellTable,1);
cellCount=zeros(length(inputFiles),numOutCells);
flexMode=strcmpi(cellTable{1,4}(1),'F');

%assume the first session file is representative of the rest.
for iFile=1:length(inputFiles)
    if ~preview
        disp(['Working on #' num2str(iFile) ': ' inputFiles{iFile} '.']);
        if iFile==1
            readArg{1}='format';
            readArg{2}=importFormat;
            readArg{3}='type';
            readArg{4}='continuous';
            readArg{5}='elecPrefs';
            readArg{6}=EPmain.preferences.general.rotateHead;
            readArg{7}='screenSize';
            readArg{8}=EPmain.scrsz;
            readArg{9}='FontSize';
            readArg{10}=EPmain.fontsize;
            
            SMIsuffix=EPmain.preferences.general.SMIsuffix;
            if ~isempty(SMIsuffix)
                readArg{end+1}='SMIsuffix';
                readArg{end+1}=SMIsuffix;
            end;
            specSuffix=EPmain.preferences.general.specSuffix;
            if ~isempty(specSuffix)
                readArg{end+1}='specSuffix';
                readArg{end+1}=specSuffix;
            end;
            subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
            if ~isempty(subjectSpecSuffix)
                readArg{end+1}='subjectSpecSuffix';
                readArg{end+1}=subjectSpecSuffix;
            end;
            BVheader=EPmain.preferences.general.BVheader;
            if ~isempty(BVheader)
                readArg{end+1}='BVheader';
                readArg{end+1}=BVheader;
            end;
            
            Name=deblank(inputFiles{1});
            thisReadArg=readArg;
            thisReadArg{end+1}='file';
            [pathstr, fileName, ext]=fileparts(Name);
            thisReadArg{end+1}=Name;
            inputData=ep_readData(thisReadArg);
            if isempty(inputData)
                return;
            end;
            readArg{end+1}='silent';
            readArg{end+1}='on';
            readArg{end+1}='ced';
            readArg{end+1}=inputData.ced; %assume all the files to be segmented will be using the same ced file.
%             readArg{end+1}='eloc';
%             readArg{end+1}=inputData.eloc; %assume all the files to be segmented will be using the same eloc info.
            readArg{end+1}='montage';
            readArg{end+1}=inputData.montage; %assume all the files to be segmented will be using the same montage.
        else
            thisReadArg=readArg;
            thisReadArg{end+1}='file';
            Name=deblank(inputFiles{iFile});
            thisReadArg{end+1}=Name;
            [pathstr, fileName, ext]=fileparts(Name);
            inputData=ep_readData(thisReadArg);
        end;
        if ~any(strcmp(inputData.dataType,{'continuous','single_trial'}))
            msg{1}=['Error: The file ' fileName ' is not a continuous or single-trial file.'];
            [msg]=ep_errorMsg(msg);
            continue
        end;
        if all(cellfun(@isempty,inputData.events))
            msg{1}=['Error: The file ' fileName ' has no events.'];
            [msg]=ep_errorMsg(msg);
            continue
        end;
        numChans=length(inputData.chanNames);
        samplesize=1000/inputData.Fs;
        outputData=inputData;
        outputData.cellNames=cell(0);
        outputData.cellTypes=cell(0);
        outputData.trialSpecs=cell(0);
        outputData.events=cell(0);
        outputData.trialNames=[];
        if strcmp(inputData.dataType,'continuous')
            outputData.dataType='single_trial';
            if flexMode
                numInterPoints=str2num(cellTable{1,4}(2:end));
                outputData.Fs=numInterPoints*10;
                outputData.timeNames=[0:numInterPoints-1]'*(1000/outputData.Fs);
                outputData.baseline=0;
                outputData.timeUnits='per';
            else
                outputData.baseline=-str2num(cellTable{1,3})/samplesize;
                outputData.timeNames=[str2num(cellTable{1,3}):samplesize:str2num(cellTable{1,4})-samplesize]';
            end;
            outputData.recTime=[];
        end;
        outputData.data=zeros(numChans,length(outputData.timeNames),0,1,1,max(1,length(outputData.freqNames)));
    else
        inputData=inputFiles;
        Name='preview';
        numChans=length(inputData.chanNames);
        samplesize=1000/inputData.Fs;
    end;
    
    if strcmp(inputData.dataType,'single_trial')
        %convert single_trial data so it is organized like continuous
            
        tempEvents=cell(1);
        numPoints=length(inputData.timeNames);
        if ~isempty(inputData.recTime)
            [~, trialOrder]=sort(inputData.recTime); %the trials are not necessarily organized chronologically
        else
            trialOrder=[1:length(inputData.cellNames)];
        end;
        if ~preview
            tempData=zeros(size(inputData.data,1),size(inputData.data,2)*size(inputData.data,3),1,size(inputData.data,4),size(inputData.data,5),size(inputData.data,6),size(inputData.data,7));
            for iTrial=1:length(trialOrder)
                theTrial=trialOrder(iTrial);
                tempData(:,(iTrial-1)*numPoints+1:iTrial*numPoints,1,:,:,:,:)=inputData.data(:,:,theTrial,:,:,:,:);
            end;
            inputData.data=tempData;
            inputData.timeNames=[1:size(tempData,2)];
        end
        
        recTime=zeros(length(trialOrder),1);
        for iTrial=1:length(trialOrder)
            theTrial=trialOrder(iTrial);
            trialEvents=inputData.events{1,theTrial};
            for iTrialEvent=1:length(trialEvents);
                trialEvents(iTrialEvent).sample=trialEvents(iTrialEvent).sample+(iTrial-1)*numPoints;
            end;
            tempEvents{1}=[tempEvents{1} trialEvents];
            recTime(iTrial)=inputData.recTime(theTrial);
        end;
        inputData.events=tempEvents;
        outputData.recTime=[];
        outputData.baseline=-str2num(cellTable{1,3})/samplesize;
        if flexMode
            numInterPoints=str2num(cellTable{1,4}(2:end));
            outputData.Fs=numInterPoints*10;
            outputData.timeNames=[0:numInterPoints-1]'*(1000/outputData.Fs);
            outputData.baseline=0;
            outputData.timeUnits='per';
        else
            outputData.timeNames=[str2num(cellTable{1,3}):samplesize:str2num(cellTable{1,4})-samplesize]';
        end;
        outputData.data=zeros(numChans,length(outputData.timeNames),0);
    end;
    
    %read in E-Prime text file if present and insert as TRSP information.
    if strcmp(dataType,'continuous')
        EPMtableData=[];
        [pathstr, name, fileSuffix] = fileparts(fileName);
        EPMfileName=[pathstr filesep name EPMsuffix];
        EPMmatchFileName=[pathstr filesep name '_EPMmatch.txt'];
        if exist(EPMfileName,'file')
            
        end
    end
            
            
            
            
            
            
    
    
    
    
    TRSPindex=find(strcmp('TRSP',{inputData.events{1}.value}));
    TRSPsamples=[inputData.events{1}(TRSPindex).sample];
    
    allEvents={inputData.events{1}.value};
    [allEvents{find(cellfun(@isempty,allEvents))}]=deal(' ');
    eventSamples=[inputData.events{1}.sample];
    [orderedEventSamples, eventOrder]=sort(eventSamples); %the events are not necessarily organized chronologically
    orderedEvents=inputData.events{1}(eventOrder);
    badTrials=cell(numOutCells,1);
    badChans=cell(numOutCells,1);
    
    for iCell=1:numOutCells       
        %generate list of potential epochs
        typeIndex=find(strcmp(cellTable{iCell,2},{inputData.events{1}.type}));
        valueIndex=find(strcmp(cellTable{iCell,2},cellfun(@num2str,{inputData.events{1}.value},'UniformOutput',false)));
        eventIndex=union(typeIndex,valueIndex);
        eventSamples=round([inputData.events{1}(eventIndex).sample]);
        delayTime=str2num(cellTable{iCell,5})/(1000/inputData.Fs);
        
        if ~flexMode
            prestim=str2num(cellTable{iCell,3})/(1000/inputData.Fs);
            poststim=str2num(cellTable{iCell,4})/(1000/inputData.Fs);
            epochLength=poststim-prestim;
            epochStart=eventSamples+prestim+delayTime;
            epochEnd=eventSamples+poststim+delayTime-1;
            
            if strcmp(inputData.dataType,'continuous')
                %check to see if full epoch within recording session
                goodEpochs=intersect(find(epochStart>0),find((epochStart+epochLength-1)<length(inputData.timeNames)+1));
                epochRecTime=epochStart(goodEpochs);
            elseif strcmp(inputData.dataType,'single_trial')
                %check to see if new epochs fall within bounds of original epochs
                eventEpoch=floor((eventSamples-1)/numPoints)+1;
                startEpoch=floor((epochStart-1)/numPoints)+1;
                endEpoch=floor((epochEnd-1)/numPoints)+1;
                goodEpochs=intersect(find(eventEpoch==startEpoch),find(eventEpoch==endEpoch));
                epochRecTime=recTime(eventEpoch);
            end;
        else
            %see if flex start event is followed by a flex end event without another intervening start event
            epochStart=eventSamples+delayTime;
            flexEndIndex=union(find(strcmp(cellTable{iCell,3},{inputData.events{1}.type})),find(strcmp(cellTable{iCell,3},{inputData.events{1}.value})));
            flexEndSamples=round([inputData.events{1}(flexEndIndex).sample]);
            epochEnd=[];
            goodEpochs=[];
            for iEvent=1:length(eventIndex)
                endSample=min(flexEndSamples(flexEndSamples>eventSamples(iEvent)));
                epochEnd(end+1)=endSample+delayTime;
                nextEventSample=min(eventSamples(eventSamples>eventSamples(iEvent)));
                if isempty(nextEventSample) || (endSample <= nextEventSample)
                    goodEpochs(end+1)=iEvent;
                end;
            end;
            epochRecTime=epochStart(goodEpochs);
        end;
        eventIndex=eventIndex(goodEpochs);
        epochStart=epochStart(goodEpochs);
        epochEnd=epochEnd(goodEpochs);
        
        if strcmp(inputData.dataType,'continuous')
            %check to see if boundary event falls within epoch
            boundaryIndex=find(strcmp('boundary',{inputData.events{1}.value}));
            boundarySamples=[inputData.events{1}(boundaryIndex).sample];
            badEpochs=[];
            for iBoundary=1:length(boundaryIndex)
                badEpochs=[badEpochs find((epochStart<=boundarySamples(iBoundary)) & ((epochStart+epochLength-1)>=boundarySamples(iBoundary)))];
            end;
            badEpochs=unique(badEpochs);
            eventIndex(badEpochs)=[];
            epochStart(badEpochs)=[];
        end;
        
        %determine which potential epochs meets trial spec criteria
        
        goodSpecs=ones(numSpecs,1);
        for iSpec=1:numSpecs
            if isempty(deblank(cellTable{iCell,3+iSpec*3}))
                %cellTable{iCell,3+iSpec*3}='none';
                goodSpecs(iSpec)=0;
            end;
            if isempty(deblank(cellTable{iCell,4+iSpec*3}))
                goodSpecs(iSpec)=0;
            elseif any(strcmp(deblank(cellTable{iCell,4+iSpec*3}),{'=','~='}))
                if isempty(deblank(cellTable{iCell,5+iSpec*3}))
                    %cellTable{iCell,5+iSpec*3}='none';
                    goodSpecs(iSpec)=0;
                end;
            elseif isempty(deblank(cellTable{iCell,4+iSpec*3}))
                %cellTable{iCell,4+iSpec*3}='none';
                goodSpecs(iSpec)=0;
            end;
        end;
        goodSpecs=find(goodSpecs);
        badEpochs=[];
        allTRSPspecValues=cell(length(eventIndex),1);
        TRSPspecNamesList=cell(0);
        sectionNumber=1;
        sectionTrial=1;
        sectionStart=0;
        sectionNumberList=[];
        sectionTrialList=[];
        for iEvent=1:length(eventIndex) %looping through events potentially defining trials
            theEvent=eventIndex(iEvent);
            stimSpecNames=cell(0); %list of spec names for this event
            stimSpecValues=cell(0);%list of the values for these specs
            stimSpecNames{end+1}='value';
            stimSpecValues{end+1}=inputData.events{1}(theEvent).value;
            for iKey=1:length(inputData.events{1}(theEvent).keys)
                stimSpecNames{end+1}=inputData.events{1}(theEvent).keys(iKey).code;
                stimSpecValues{end+1}=inputData.events{1}(theEvent).keys(iKey).data;
            end;
            
            TRSPspecNames=cell(0);
            TRSPspecValues=cell(0);
            if strcmp(inputData.dataType,'continuous')
                epochTRSP=TRSPindex(find((TRSPsamples>=epochStart(iEvent))));
                if ~isempty(epochTRSP)
                    epochTRSP=epochTRSP(1);
                    if (iEvent==length(eventIndex)) || (inputData.events{1}(epochTRSP).sample < inputData.events{1}(eventIndex(iEvent+1)).sample)
                        for iKey=1:length(inputData.events{1}(epochTRSP).keys)
                            TRSPspecNames{end+1,1}=['TS-' inputData.events{1}(epochTRSP).keys(iKey).code];
                            TRSPspecValues{end+1,1}=inputData.events{1}(epochTRSP).keys(iKey).data;
                        end;
                    end;
                end;
                allTRSPspecValues{iEvent}=TRSPspecValues;
            elseif strcmp(inputData.dataType,'single_trial')
                %theEpoch=floor((eventSamples(eventIndex(iEvent))-1)/length(inputData.timeNames))+1;
                for iKey=1:length(inputData.trialSpecNames)
                    TRSPspecNames{end+1,1}=['TS-' inputData.trialSpecNames{iKey}];
                    TRSPspecValues{end+1,1}=inputData.trialSpecs{iEvent,iKey};
                end;
                allTRSPspecValues{iEvent}=TRSPspecValues;                
            else
                msg{1}=['Error: The file ' fileName ' is not a continuous or single-trial file.'];
                [msg]=ep_errorMsg(msg);
                continue
            end;
            if ~isempty(TRSPspecNames)
                TRSPspecNamesList=TRSPspecNames;
            end;
            
            keepTrial=1;
            afterPrecede=0;
            firstPrecede=0; %event number (sorted) of the first follow criterion event
            precedeNeg=0; %next spec is following a ~= -precedes- criterion            
            afterFollow=0;
            firstFollow=0; %event number (sorted) of the first follow criterion event
            followNeg=0; %next spec is following a ~= -follows- criterion
            for iSpec=1:length(goodSpecs)
                theSpec=goodSpecs(iSpec);
                stimSpec=find(strcmp(cellTable{iCell,3+theSpec*3},stimSpecNames));
                if length(stimSpec) > 1
                    disp(['Error: multiple matching stimulus specs (' stimSpecNames{stimSpec(1)} ') for ' Name '.  Will use the first instance.']);
                    stimSpec=stimSpec(1);
                    %keepTrial=0;
                end;
                noStim=0;
                badStim=0;
                if afterFollow && ~any(ismember(theSpec-1,goodSpecs))
                    afterFollow=0;
                end;
                if afterPrecede && ~any(ismember(theSpec-1,goodSpecs))
                    afterPrecede=0;
                end;
                
                %if it is a '-precedes-' spec
                if strcmp('-precedes-',cellTable{iCell,3+theSpec*3})
                    afterPrecede=1;
                    thePrecedeName=cellTable{iCell,5+theSpec*3};
                    if firstPrecede %if there is already a precede in effect
                        precede1=min(find(orderedEventSamples==inputData.events{1}(theEvent).sample))+1; %ranked order
                        precede2=firstPrecede-1; %ranked order
                    else
                        precede1=min(find(orderedEventSamples==inputData.events{1}(theEvent).sample))+1; %ranked order
                        precede2=length(orderedEvents); %ranked order
                    end;
                    switch deblank(cellTable{iCell,4+theSpec*3})
                        case '='
                            if isempty(find(strcmp(thePrecedeName,{orderedEvents(precede1:precede2).value})))
                                badStim=1;
                            else
                                firstPrecede=max(find(strcmp(thePrecedeName,{orderedEvents(precede1:precede2).value}))); %if there is a further "precedes", it will check between the lock event and the first one.
                            end;
                        case '~='
                            if ~isempty(find(strcmp(thePrecedeName,allEvents(eventOrder(precede1:precede2)))))
                                if isempty(intersect(theSpec+1,goodSpecs)) %if bad for all specs
                                    badStim=1;
                                else
                                    precedeNeg=1;
                                end;
                            end;
                        case '<'
                            theEventsList={orderedEvents(precede1:precede2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                precedeList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}<thePrecedeName
                                            precedeList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                precedeEvent=min(precedeList);
                                if isempty(precedeEvent)
                                    badStim=1;
                                else
                                    firstPrecede=precedeEvent-1+precede1; %of the events in the range, the most immediately succeeding that is less in alphabetical order than the specified name.
                                end;
                            end;
                        case '>'
                            theEventsList={orderedEvents(precede1:precede2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                precedeList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}>thePrecedeName
                                            precedeList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                precedeEvent=min(precedeList);
                                if isempty(precedeEvent)
                                    badStim=1;
                                else
                                    firstPrecede=precedeEvent-1+precede1; %of the events in the range, the most immediately succeeding that is greater in alphabetical order than the specified name.
                                end;
                            end;
                        case '<='
                            theEventsList={orderedEvents(precede1:precede2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                precedeList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}<=thePrecedeName
                                            precedeList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                precedeEvent=min(precedeList);
                                if isempty(precedeEvent)
                                    badStim=1;
                                else
                                    firstPrecede=precedeEvent-1+precede1; %of the events in the range, the most immediately succeeding that is less than or equal in alphabetical order than the specified name.
                                end;
                            end;
                        case '>='
                            theEventsList={orderedEvents(precede1:precede2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                precedeList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}<thePrecedeName
                                            precedeList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                precedeEvent=min(precedeList);
                                if isempty(precedeEvent)
                                    badStim=1;
                                else
                                    firstPrecede=precedeEvent-1+precede1; %of the events in the range, the most immediately succeeding that is greater than or equal in alphabetical order than the specified name.
                                end;
                            end;
                        case 'starts'
                            eventCheckList=[];
                            for iEventCheck=precede1:precede2
                                if min(strfind(allEvents{eventOrder(iEventCheck)},thePrecedeName))==1
                                    eventCheckList=[eventCheckList; iEventCheck];
                                end;
                            end;
                            
                            if isempty(eventCheckList)
                                badStim=1;
                            else
                                firstPrecede=max(eventCheckList);
                            end;
                        case 'ends'
                            eventCheckList=[];
                            for iEventCheck=precede1:precede2
                                if max(strfind(allEvents{eventOrder(iEventCheck)},thePrecedeName))==(length(allEvents{iEventCheck})-length(thePrecedeName)+1)
                                    eventCheckList=[eventCheckList; iEventCheck];
                                end;
                            end;
                            
                            if isempty(eventCheckList)
                                badStim=1;
                            else
                                firstPrecede=max(eventCheckList);
                            end;
                        case 'contains'
                            eventCheckList=[];
                            for iEventCheck=precede1:precede2
                                if ~isempty(strfind(allEvents{eventOrder(iEventCheck)},thePrecedeName))
                                    eventCheckList=[eventCheckList; iEventCheck];
                                end;
                            end;
                            
                            if isempty(eventCheckList)
                                badStim=1;
                            else
                                firstPrecede=max(eventCheckList);
                            end;
                        otherwise
                            disp(['Error: cell table spec relationship (' cellTable{iCell,4+theSpec*3} ') not valid.']);
                            return
                    end;
                    if badStim
                        keepTrial=0;
                    end;
                    
                %if it is a '-follows-' spec
                elseif strcmp('-follows-',cellTable{iCell,3+theSpec*3})
                    afterFollow=1;
                    theFollowName=cellTable{iCell,5+theSpec*3};
                    if firstFollow %if there is already a follow in effect
                        follow1=firstFollow+1; %ranked order
                        follow2=min(find(orderedEventSamples==inputData.events{1}(theEvent).sample))-1; %ranked order
                    else
                        follow1=1; %ranked order
                        follow2=min(find(orderedEventSamples==inputData.events{1}(theEvent).sample))-1; %ranked order
                    end;
                    switch deblank(cellTable{iCell,4+theSpec*3})
                        case '='
                            if isempty(find(strcmp(theFollowName,{orderedEvents(follow1:follow2).value})))
                                badStim=1;
                            else
                                firstFollow=max(find(strcmp(theFollowName,{orderedEvents(follow1:follow2).value}))); %if there is a further "follow", it will check between the first one and the lock event.
                            end;
                        case '~='
                            if ~isempty(find(strcmp(theFollowName,allEvents(eventOrder(follow1:follow2)))))
                                if isempty(intersect(theSpec+1,goodSpecs)) %if bad for all specs
                                    badStim=1;
                                else
                                    followNeg=1;
                                end;
                            end;
                        case '<'
                            theEventsList={orderedEvents(follow1:follow2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                followsList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}<theFollowName
                                            followsList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                followEvent=max(followsList);
                                if isempty(followEvent)
                                    badStim=1;
                                else
                                    firstFollow=followEvent-1+follow1; %of the events in the range, the most immediately prior that is less in alphabetical order than the specified name.
                                end;
                            end;
                        case '>'
                            theEventsList={orderedEvents(follow1:follow2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                followsList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}>theFollowName
                                            followsList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                followEvent=max(followsList);
                                if isempty(followEvent)
                                    badStim=1;
                                else
                                    firstFollow=followEvent-1+follow1; %of the events in the range, the most immediately prior that is greater in alphabetical order than the specified name.
                                end;
                            end;
                        case '<='
                            theEventsList={orderedEvents(follow1:follow2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                followsList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}<=theFollowName
                                            followsList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                followEvent=max(followsList);
                                if isempty(followEvent)
                                    badStim=1;
                                else
                                    firstFollow=followEvent-1+follow1; %of the events in the range, the most immediately prior that is less than or equal in alphabetical order than the specified name.
                                end;
                            end;
                        case '>='
                            theEventsList={orderedEvents(follow1:follow2).value};
                            if isempty(theEventsList)
                                badStim=1;
                            else
                                followsList=[];
                                for iList=1:length(theEventsList)
                                    if ~isempty(theEventsList{iList})
                                        if theEventsList{iList}>theFollowName
                                            followsList(end+1)=iList;
                                        end;
                                    end;
                                end;
                                followEvent=max(followsList);
                                if isempty(followEvent)
                                    badStim=1;
                                else
                                    firstFollow=followEvent-1+follow1; %of the events in the range, the most immediately prior that is greater or equal in alphabetical order than the specified name.
                                end;
                            end;
                        case 'starts'
                            eventCheckList=[];
                            for iEventCheck=follow1:follow2
                                if min(strfind(allEvents{eventOrder(iEventCheck)},theFollowName))==1
                                    eventCheckList=[eventCheckList; iEventCheck];
                                end;
                            end;
                            
                            if isempty(eventCheckList)
                                badStim=1;
                            else
                                firstFollow=max(eventCheckList);
                            end;
                        case 'ends'
                            eventCheckList=[];
                            for iEventCheck=follow1:follow2
                                if max(strfind(allEvents{eventOrder(iEventCheck)},theFollowName))==(length(allEvents{iEventCheck})-length(theFollowName)+1)
                                    eventCheckList=[eventCheckList; iEventCheck];
                                end;
                            end;
                            
                            if isempty(eventCheckList)
                                badStim=1;
                            else
                                firstFollow=max(eventCheckList);
                            end;
                        case 'contains'
                            eventCheckList=[];
                            for iEventCheck=follow1:follow2
                                if ~isempty(strfind(allEvents{eventOrder(iEventCheck)},theFollowName))
                                    eventCheckList=[eventCheckList; iEventCheck];
                                end;
                            end;
                            
                            if isempty(eventCheckList)
                                badStim=1;
                            else
                                firstFollow=max(eventCheckList);
                            end;
                        otherwise
                            disp(['Error: cell table spec relationship (' cellTable{iCell,4+theSpec*3} ') not valid.']);
                            return
                    end;
                    if badStim
                        keepTrial=0;
                    else
                        if sectionStart ~= firstFollow
                            sectionStart=firstFollow;
                            sectionNumber=sectionNumber+1;
                            sectionTrial=1;
                        end;
                    end;
                    
                    %if it is a regular spec
                elseif ~any(strcmp(cellTable{iCell,3+theSpec*3},{'-secNum-','-secTrial-','-secTrialFromEnd-'}))
                    theCriterion=cellTable{iCell,5+theSpec*3};
                    if afterPrecede %match to -precedes- event specs rather than the lock event specs
                        thePrecedeEvent=firstPrecede;
                        precedeStimSpecNames=cell(0);
                        precedeSpecSpecValues=cell(0);
                        for iKey=1:length(inputData.events{1}(eventOrder(thePrecedeEvent)).keys)
                            theCode=inputData.events{1}(eventOrder(thePrecedeEvent)).keys(iKey).code;
                            theData=inputData.events{1}(eventOrder(thePrecedeEvent)).keys(iKey).data;
                            if ~isempty(theCode) && ~isempty(theData)
                                precedeStimSpecNames{end+1}=theCode;
                                precedeSpecSpecValues{end+1}=theData;
                            end;
                        end;
                        if ~isempty(followStimSpecNames) %if the -precedes- event has no specs and yet a spec condition was specified, then treated as no precede spec was present, resulting in a bad epoch.
                            stimSpec=find(strcmp(cellTable{iCell,3+theSpec*3},precedeStimSpecNames));
                            theStim=num2str(precedeSpecSpecValues{stimSpec});
                        end;
                    elseif afterFollow %match to -follows- event specs rather than the lock event specs
                        theFollowEvent=firstFollow;
                        followStimSpecNames=cell(0);
                        followSpecSpecValues=cell(0);
                        for iKey=1:length(inputData.events{1}(eventOrder(theFollowEvent)).keys)
                            theCode=inputData.events{1}(eventOrder(theFollowEvent)).keys(iKey).code;
                            theData=inputData.events{1}(eventOrder(theFollowEvent)).keys(iKey).data;
                            if ~isempty(theCode) && ~isempty(theData)
                                followStimSpecNames{end+1}=theCode;
                                followSpecSpecValues{end+1}=theData;
                            end;
                        end;
                        if ~isempty(followStimSpecNames) %if the -follows- event has no specs and yet a spec condition was specified, then treated as no follow spec was present, resulting in a bad epoch.
                            stimSpec=find(strcmp(cellTable{iCell,3+theSpec*3},followStimSpecNames));
                            theStim=num2str(followSpecSpecValues{stimSpec});
                        end;
                    else
                        if ~isempty(stimSpec)
                            theStim=num2str(stimSpecValues{stimSpec});
                        else
                            theStim=[];
                        end;
                    end;
                    
                    if isempty(stimSpec)
                        noStim=1;
                    else
                        cmp=ep_compareStrings(theStim,theCriterion);
                        if isempty(cmp)
                            badStim=1;
                        else
                            switch deblank(cellTable{iCell,4+theSpec*3})
                                case '='
                                    if cmp ~=0
                                        badStim=1;
                                    end;
                                case '~='
                                    if cmp ==0
                                        badStim=1;
                                    end;
                                case '<'
                                    if cmp > -1
                                        badStim=1;
                                    end;
                                case '>'
                                    if cmp < 1
                                        badStim=1;
                                    end;
                                case '<='
                                    if cmp > 0
                                        badStim=1;
                                    end;
                                case '>='
                                    if cmp < 0
                                        badStim=1;
                                    end;
                                case 'starts'
                                    if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(1:length(theCriterion)))
                                        badStim=1;
                                    end;
                                case 'ends'
                                    if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(end-length(theCriterion)+1:end))
                                        badStim=1;
                                    end;
                                case 'contains'
                                    if isempty(findstr(theCriterion,theStim))
                                        badStim=1;
                                    end;
                                otherwise
                                    disp(['Error: cell table spec relationship (' cellTable{iCell,4+theSpec*3} ') not valid.']);
                                    return
                            end;
                        end;
                        if (badStim && followNeg) || (badStim && precedeNeg)
                            badStim=0;
                        end;
                        followNeg=0;
                        precedeNeg=0;
                    end;
                    if badStim
                        keepTrial=0;
                    end;
                    
                    noTRSP=0;
                    if ~afterFollow && ~afterPrecede %TRSP crits do not succeed a -follows- or -precedes- criteria
                        TRSPspec=find(strcmp(cellTable{iCell,3+theSpec*3},TRSPspecNames));
                        if length(TRSPspec) > 1
                            disp(['Error: multiple matching TRSP specs (' TRSPspecNames{TRSPspec(1)} ') for ' Name '.']);
                            keepTrial=0;
                        end;
                        badTRSP=0;
                        if isempty(TRSPspec)
                            noTRSP=1;
                        else
                            cmp=ep_compareStrings(TRSPspecValues{TRSPspec},cellTable{iCell,5+theSpec*3});
                            if isempty(cmp)
                                badTRSP=1;
                            else
                                switch deblank(cellTable{iCell,4+theSpec*3})
                                    case '='
                                        if cmp ~=0
                                            badTRSP=1;
                                        end;
                                    case '~='
                                        if cmp ==0
                                            badTRSP=1;
                                        end;
                                    case '<'
                                        if cmp > -1
                                            badTRSP=1;
                                        end;
                                    case '>'
                                        if cmp < 1
                                            badTRSP=1;
                                        end;
                                    case '<='
                                        if cmp > 0
                                            badTRSP=1;
                                        end;
                                    case '>='
                                        if cmp < 0
                                            badTRSP=1;
                                        end;
                                    case 'starts'
                                        theStim=num2str(TRSPspecValues{TRSPspec});
                                        theCriterion=cellTable{iCell,5+theSpec*3};
                                        if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(1:length(theCriterion)))
                                            badTRSP=1;
                                        end;
                                    case 'ends'
                                        theStim=num2str(TRSPspecValues{TRSPspec});
                                        theCriterion=cellTable{iCell,5+theSpec*3};
                                        if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(end-length(theCriterion)+1:end))
                                            badTRSP=1;
                                        end;
                                    case 'contains'
                                        theStim=num2str(TRSPspecValues{TRSPspec});
                                        theCriterion=cellTable{iCell,5+theSpec*3};
                                        if isempty(findstr(theCriterion,theStim))
                                            badTRSP=1;
                                        end;
                                    otherwise
                                end;
                            end;
                            if badTRSP
                                keepTrial=0;
                            end;
                        end;
                    else
                        afterFollow=0;
                        afterPrecede=0;
                    end;
                    if noStim && noTRSP
                        keepTrial=0;
                    end;
                end;%follow or not -follows-
                if ~keepTrial
                    break %no need to check rest of specs
                end;
            end;%iSpec
            if ~keepTrial
                badEpochs=[badEpochs iEvent];
            else
                sectionNumberList=[sectionNumberList sectionNumber];
                sectionTrialList=[sectionTrialList sectionTrial];
                sectionTrial=sectionTrial+1;
            end;
        end;%iEvent
        eventIndex(badEpochs)=[];
        epochStart(badEpochs)=[];
        epochEnd(badEpochs)=[];
        allTRSPspecValues(badEpochs)=[];
        
        sectionNumHist=histc(sectionNumberList,unique(sectionNumberList));
        sectionNumTrials=[];
        for iSection=1:length(sectionNumHist)
            sectionNumTrials=[sectionNumTrials repmat(sectionNumHist(iSection),1,sectionNumHist(iSection))];
        end;
        sectionTrialFromEndList=sectionNumTrials-sectionTrialList+1;
        
        %Check for section crits now that the sections have been fully defined
        badEpochs=[];
        for iEvent=1:length(eventIndex) %looping through events potentially defining trials
            badStim=0;
            for iSpec=1:length(goodSpecs)
                theSpec=goodSpecs(iSpec);
                theCrit=cellTable{iCell,3+theSpec*3};
                theStim=cellTable{iCell,5+theSpec*3};
                if ischar(theStim)
                    theStim=str2num(theStim);
                end;
                if any(strcmp(theCrit,{'-secNum-','-secTrial-','-secTrialFromEnd-'}))
                    switch theCrit
                        case '-secNum-'
                            theValue=sectionNumberList(iEvent);
                        case '-secTrial-'
                            theValue=sectionTrialList(iEvent);
                        case '-secTrialFromEnd-'
                            theValue=sectionTrialFromEndList(iEvent);
                    end;
                    switch deblank(cellTable{iCell,4+theSpec*3})
                        case '='
                            if theValue ~= theStim
                                badStim=1;
                            end;
                        case '~='
                            if theValue == theStim
                                badStim=1;
                            end;
                        case '<'
                            if theValue >= theStim
                                badStim=1;
                            end;
                        case '>'
                            if theValue <= theStim
                                badStim=1;
                            end;
                        case '<='
                            if theValue > theStim
                                badStim=1;
                            end;
                        case '>='
                            if theValue < theStim
                                badStim=1;
                            end;
                    end;
                end;
            end;
            if badStim
                badEpochs=[badEpochs iEvent];
            end;
        end;
        eventIndex(badEpochs)=[];
        epochStart(badEpochs)=[];
        epochEnd(badEpochs)=[];
        allTRSPspecValues(badEpochs)=[];
        sectionNumberList(badEpochs)=[];
        sectionTrialList(badEpochs)=[];
        sectionTrialFromEndList(badEpochs)=[];
        
        %if there is a TS-EPoffset field, then adjust the epochs accordingly.
        if any(strcmp('TS-EPoffset',TRSPspecNamesList))
            offsetTRSP=find(strcmp('TS-EPoffset',TRSPspecNamesList));
            disp('Adjusting the epoch timepoints according to the TS-EPoffset values')
            badEpochs=[];
            for iEvent=1:length(eventIndex) %looping through events potentially defining trials
                theOffset=round(str2num(allTRSPspecValues{iEvent}{offsetTRSP})/(1000/inputData.Fs));
                if ~isempty(theOffset)
                    epochStart(iEvent)=epochStart(iEvent)+theOffset;
                    epochEnd(iEvent)=epochEnd(iEvent)+theOffset;
                    if (epochStart(iEvent) < 1) || (epochEnd(iEvent) > length(inputData.timeNames))
                        badEpochs(end+1)=iEvent;
                    end;
                end;
            end;
            eventIndex(badEpochs)=[];
            epochStart(badEpochs)=[];
            epochEnd(badEpochs)=[];
            allTRSPspecValues(badEpochs)=[];
            sectionNumberList(badEpochs)=[];
            sectionTrialList(badEpochs)=[];
            sectionTrialFromEndList(badEpochs)=[];
        end;
        
        cellCount(iFile,iCell)=length(eventIndex);
        if ~preview
            numCellTrials=length(eventIndex);
            if numCellTrials > 0
                for iSpec=1:length(TRSPspecNamesList)
                    TRSPspecNamesList{iSpec}=TRSPspecNamesList{iSpec}(4:end);
                end;
                if ~any(strcmp('secNum',TRSPspecNamesList))
                    TRSPspecNamesList{end+1}='secNum';
                end;
                if ~any(strcmp('secTrial',TRSPspecNamesList))
                    TRSPspecNamesList{end+1}='secTrial';
                end;
                if ~any(strcmp('secTrialFromEnd',TRSPspecNamesList))
                    TRSPspecNamesList{end+1}='secTrialFromEnd';
                end;
                secNumTRSP=find(strcmp('secNum',TRSPspecNamesList));
                secTrialTRSP=find(strcmp('secTrial',TRSPspecNamesList));
                secTrialFromEndTRSP=find(strcmp('secTrialFromEnd',TRSPspecNamesList));
                sectionNumTrials=histc(sectionNumberList,unique(sectionNumberList));
                for iEpoch=1:numCellTrials 
                    allTRSPspecValues{iEpoch}{secNumTRSP}=sectionNumberList(iEpoch);
                    allTRSPspecValues{iEpoch}{secTrialTRSP}=sectionTrialList(iEpoch);
                    allTRSPspecValues{iEpoch}{secTrialFromEndTRSP}=sectionTrialFromEndList(iEpoch);
                end;
                if iCell==1 %assume same trial specs for entire file
                    outputData.trialSpecNames=TRSPspecNamesList;
                end;
                allEventSamples=round([inputData.events{1}.sample]);
                trialCountStart=0;
                sameCellName=find(strcmp(cellTable{iCell,1},cellTable(1:iCell-1,1)));
                if ~isempty(sameCellName)
                    trialCountStart=sum(cellCount(iFile,sameCellName));
                end;
                if flexMode
                    badPoints=find(kron([inputData.analysis.badTrials(1,:)],ones(1,inputData.Fs)));
                    badChanPoints=cell(0);
                    if strcmp(inputData.dataType,'continuous')
                        for iChan=1:length(inputData.chanNames)
                            badChanPoints{iChan}=find(kron([inputData.analysis.badChans(1,:,iChan)==-1],ones(1,inputData.Fs)));
                        end;
                    else
                        for iChan=1:length(inputData.chanNames)
                            badChanPoints{iChan}=find(kron([inputData.analysis.badChans(1,:,iChan)==-1],ones(1,numPoints)));
                        end;
                    end;
                end;
                badTrials{iCell}=zeros(1,numCellTrials);
                badChans{iCell}=zeros(1,numCellTrials,length(inputData.chanNames));
                for iEpoch=1:numCellTrials
                    if flexMode
                        deltaT=floor((epochEnd(iEpoch)-epochStart(iEpoch)+1)/numInterPoints);
                        if strcmp(inputData.dataType,'continuous')
                            goodPoints=epochStart(iEpoch):epochEnd(iEpoch);
                            goodPoints=setdiff(goodPoints,badPoints);
                        else
                            theInTrial=(floor(epochStart(iEpoch)-1)/numPoints)+1;
                            if inputData.analysis.badTrials(1,theInTrial)
                                goodPoints=[];
                            else
                                goodPoints=epochStart(iEpoch):epochEnd(iEpoch);
                            end;
                        end;
                        outputData.data(:,:,end+1,1,1,1)=0;
                        if isempty(goodPoints)
                            badTrials{iCell}(1,iEpoch)=1;
                        else
                            for iChan=1:length(inputData.chanNames)
                                goodChanPoints=setdiff(goodPoints,badChanPoints{iChan});
                                for iPoint=1:numInterPoints
                                    epochPoints=intersect([(iPoint-1)*deltaT+epochStart(iEpoch):iPoint*deltaT+epochStart(iEpoch)-1],goodChanPoints);
                                    if ~isempty(epochPoints)
                                        if ~isempty(inputData.freqNames) && any(strcmp(inputData.chanTypes{iChan},{'EEG','REG'}))
                                            %switch to amplitude scaling when adding freq data together.
                                            outputData.data(iChan,iPoint,end,1,:,:)=squeeze(mean(abs(inputData.data(iChan,epochPoints,1,1,:,:)),2));
                                        else
                                            outputData.data(iChan,iPoint,end,1,:,:)=squeeze(mean(inputData.data(iChan,epochPoints,1,1,:,:),2));
                                        end;
                                    else
                                        outputData.data(:,iPoint,end,1,1,1)=0;
                                        badChans{iCell}(1,iEpoch,iChan)=-1;
                                    end;
                                end;
                            end;
                        end;
                        outputData.events{1,end+1}=inputData.events{1}(find((allEventSamples>=epochStart(iEpoch)) & (allEventSamples <= epochEnd(iEpoch))));
                        for iEvent=1:length(outputData.events{1,end})
                            outputData.events{1,end}(iEvent).sample=floor((outputData.events{1,end}(iEvent).sample-epochStart(iEpoch)+1-1)/numInterPoints)+1; %change event time to be relative to epoch
                        end;
                    else
                        outputData.data(:,:,end+1,1,:,:)=inputData.data(:,epochStart(iEpoch):epochEnd(iEpoch),1,1,:,:);
                        outputData.events{1,end+1}=inputData.events{1}(find((allEventSamples>=epochStart(iEpoch)) & (allEventSamples <= epochStart(iEpoch)+epochLength-1)));
                        for iEvent=1:length(outputData.events{1,end})
                            outputData.events{1,end}(iEvent).sample=outputData.events{1,end}(iEvent).sample-epochStart(iEpoch)+1; %change event time to be relative to epoch
                        end;
                        if strcmp(inputData.dataType,'continuous')
                            theInTrial=floor((epochStart(iEpoch)-1)/inputData.Fs)+1;
                        else
                            theInTrial=floor((epochStart(iEpoch)-1)/numPoints)+1;
                        end;
                        badTrials{iCell}(1,iEpoch)=inputData.analysis.badTrials(1,theInTrial);
                        badChans{iCell}(1,iEpoch,:)=inputData.analysis.badChans(1,theInTrial,:);
                    end;
                    outputData.cellNames{end+1,1}=cellTable{iCell,1};
                    outputData.cellTypes{end+1,1}='SGL';
                    outputData.trialNames(end+1,1)=iEpoch+trialCountStart;
                    if isempty(TRSPspecNamesList)
                        outputData.trialSpecs{end+1,1}=[];
                    else
                        for iSpec=1:length(TRSPspecNamesList)
                            if iSpec==1
                                if isempty(allTRSPspecValues{iEpoch})
                                    outputData.trialSpecs{end+1,iSpec}=[];
                                else
                                    outputData.trialSpecs{end+1,iSpec}=allTRSPspecValues{iEpoch}{iSpec};
                                end;
                            else
                                if isempty(allTRSPspecValues{iEpoch})
                                    outputData.trialSpecs{end,iSpec}=[];
                                else
                                    outputData.trialSpecs{end,iSpec}=allTRSPspecValues{iEpoch}{iSpec};
                                end;
                            end;
                        end;
                    end;
                    outputData.recTime(end+1,1)=epochRecTime(iEpoch);
                end;
            end;
        end;
    end;%iCell
    if ~preview
        numTrials=length(outputData.cellNames);
        if numTrials==0
            disp(['Error: Segmented file for ' Name ' did not result in any segments and therefore will not be saved.']);
        else
            outputData.avgNum=ones(1,numTrials);
            outputData.subNum=ones(1,numTrials);
            outputData.covNum=ones(1,numTrials);
            outputData.history{end+1}={'ep_segmentData',cellTable,inputFiles,importFormat,outputFormat,preview};
            outputData.analysis.blinkTrial=zeros(1,numTrials);
            outputData.analysis.saccadeTrial=zeros(1,numTrials);
            outputData.analysis.saccadeOnset=zeros(1,numTrials);
            outputData.analysis.moveTrial=zeros(1,numTrials);
            outputData.analysis.badTrials=zeros(1,numTrials);
            outputData.analysis.badChans=zeros(1,numTrials,numChans);
            count=0;
            for iCell=1:numOutCells
                outputData.analysis.badTrials(1,count+1:count+length(badTrials{iCell}))=badTrials{iCell};
                outputData.analysis.badChans(1,count+1:count+length(badTrials{iCell}),:)=badChans{iCell};
                count=count+length(badTrials{iCell});
            end;
            if isempty(outputData.trialSpecNames)
                outputData.trialSpecs=cell(numTrials,0);
            end;
            outputData.stims=struct('name',{},'image',{});
            for iStim=1:length(inputData.stims)
                stimFlag=0;
                for iCell=1:length(outputData.cellNames)
                    for iEvent=1:length(outputData.events{iCell})
                        if any(strcmp(inputData.stims(iStim).name,{outputData.events{iCell}(iEvent).keys.data}))
                            stimFlag=1;
                        end;
                    end;
                end;
                if stimFlag
                    outputData.stims(end+1)=inputData.stims(iStim); %keep only stim images whose events are still in the segmented data.
                end;
            end;
            
            [err]=ep_checkEPfile(outputData);
            if err
                disp(['Error: Segmented file for ' Name ' did not pass data integrity checks and therefore will not be saved.']);
            else
                [pathstr, name, ext] = fileparts(inputFiles{iFile});
                [fileSuffix,formatName]=ep_fileExtensions(outputFormat);
                
                sameName=1;
                theNumber=0;
                fileNameSuffix=[pathstr filesep name '_seg' fileSuffix];
                while sameName
                    sameName=0;
                    if exist(fileNameSuffix,'file')
                        sameName=1;
                    end;
                    if sameName
                        theNumber=theNumber+1;
                        fileNameSuffix=[pathstr filesep name '_seg-' num2str(theNumber) fileSuffix];
                    end;
                end;
                [err]=ep_writeData(outputData,fileNameSuffix,EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,outputFormat);
            end;
        end;
    end;
end;
cellNums=cellCount;





