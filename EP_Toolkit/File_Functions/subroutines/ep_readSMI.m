function [EPdataOut msgLog]=ep_readSMI(EPdataIn,fileName,matchTable);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EPdataOut=ep_readSMI(EPdataIn,fileName,matchTable);
% Reads in text output from SMI eye tracker and adds the sample-by-sample information to the data field and the events
% to the events field.  Assumes that the EEG is a continuous data file. Interpolates the SMI data to correct for uneven
% temporal sampling, using matching EEG-SMI events as the anchor points.
%
%Inputs
%   EPdataIn: Structured array with the data and accompanying information.  See readData.
%   fileName: name of the SMI text file
%   matchTable: Table with first column of unique EEG event values and a second table of the SMI events they correspond
%   to, if any.  If they do not correspond, then the value will be 'none'.
%
%Outputs
%   EPdataOut: Structured array with the data and accompanying information.  See readData.
%   msgLog: messages.
%
% History:
%
% by Joseph Dien (6/15/14)
% jdien07@mac.com
%
% bugfix 1/28/16 JD
% Fixed not detecting SMI blink events.
%
% bugfix 1/31/16 JD
% Fixed SMI_EEGoffset off by one sample.
%
% modified 11/2/16 JD
% Sets eye-tracker data to NaN for pauses over 50 ms and if it ends prior to the EEG.
% Eye-tracker data and ETsaccade and ETfixation events blanked out during blinks.
% Drop the "# Message: " from the stim file names.
% Added option to read in AOI information from SMI "word" files.
% Added mouse clicks.
%
% bugfix 12/5/16 JD
% Using POR rather than GVEC columns for SMI eye-tracker coordinates.
% Using last AOI of fixation and of saccade periods rather than first per SMI software.
% If first SMI sample of a section is a fixation, treat it as start of fixation regardless of whether there was no SMI info or a fixation prior to it (stimulus presentation serves as start of fixation event).
%
% bugfix 4/17/17 JD
% Fixed crash when there is only a single sample in the section.
% Fixed crash due to rounding error making it seem as if SMI started prior to EEG even if actually at the same time.
% Since something in the experiment procedure is messing up the clock synchs at the start of each run, just drop the SMI data prior to the first matchList event.


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

tic
EPdataOut=[];
msgLog=cell(0);

if isempty(matchTable)
    disp('Error: no match table specified.');
    return
end;

if size(matchTable,2) ~= 2
    disp('Error: match table defective.');
    return
end;

if isempty(EPdataIn)
    disp('Error: no EEG data file.');
    return
end;

if ~strcmp(EPdataIn.dataType,'continuous')
    disp('Error: EEG data file is not continuous.');
    return
end;

%read in the SMI text file

fid=fopen(fileName);
if fid==-1
    disp('Error: no such file.');
    return
end;

[pathstr, name, fileSuffix] = fileparts(fileName);

commentLine=1;
numComments=0;
while commentLine
    tempLine=fgetl(fid);
    if strcmp(tempLine(1:2),'##')
        numComments=numComments+1;
    else
        commentLine=0;
    end;
end;

delim='\t';
numcols=length(regexp(tempLine,delim))+1; %determine number of columns based on number of delimiters

frewind(fid);
for i=1:numComments
    tempLine=fgetl(fid);
end;

tempLine=fgetl(fid);
theHeaders=strsplit(tempLine,delim);

theData=textscan(fid, repmat('%s',1,numcols),'Delimiter',delim);
fclose(fid);

TimeCol=find(strcmp('Time',theHeaders));
TypeCol=find(strcmp('Type',theHeaders));
msgCol=find(strcmp('L Raw X [px]',theHeaders));
LeventCol=find(strcmp('L Event Info',theHeaders));
ReventCol=find(strcmp('R Event Info',theHeaders));
stimulusCol=find(strcmp('L AOI Hit',theHeaders));
LpupilCol=find(strcmp('L Pupil Diameter [mm]',theHeaders));
LXeyeCol=find(strcmp('L POR X [px]',theHeaders));
LYeyeCol=find(strcmp('L POR Y [px]',theHeaders));
RpupilCol=find(strcmp('R Pupil Diameter [mm]',theHeaders));
RXeyeCol=find(strcmp('R POR X [px]',theHeaders));
RYeyeCol=find(strcmp('R POR Y [px]',theHeaders));
LplaneCol=find(strcmp('L Plane',theHeaders));
RplaneCol=find(strcmp('R Plane',theHeaders));

if isempty(TimeCol) || isempty(TypeCol) || isempty(msgCol) || isempty(LeventCol) ||...
        isempty(ReventCol) || isempty(stimulusCol) || isempty(LpupilCol) || isempty(LXeyeCol) || isempty(LYeyeCol) ||...
        isempty(RpupilCol) || isempty(RXeyeCol) || isempty(RYeyeCol)
    msg=['Error: The file ' fileName ' did not contain the expected column headers.'];
    disp(msg);
    msgLog{end+1}=msg;
    return
end;

%load in word files with AOI data since SMI software has bug that precludes it being included in the full raw data output
if exist([pathstr filesep name(1:end-4) '_word01_smi.txt'],'file')
    disp('Reading in AOI information from SMI word files.')
    fprintf('%60s\n',' ' );
    fileCounter=1;
    wordFid=0;
    while wordFid ~= -1
        wordFile=[pathstr filesep name(1:end-4) '_word' sprintf('%02d',fileCounter) '_smi.txt'];
        if exist(wordFile,'file')
            wordFid=fopen(wordFile);
        else
            wordFid=-1;
        end;
        if wordFid ~= -1
            fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('Loading AOI file: %4d', fileCounter))
            commentLine=1;
            numComments=0;
            while commentLine
                tempLine=fgetl(wordFid);
                if strcmp(tempLine(1:2),'##')
                    numComments=numComments+1;
                else
                    commentLine=0;
                end;
            end;
            
            frewind(wordFid);
            for i=1:numComments+1 %comments plus the header line
                tempLine=fgetl(wordFid);
            end;
            
            tempData=textscan(fid, repmat('%s',1,numcols),'Delimiter',delim);
            fclose(wordFid);
            SMProws=strcmp(theData{TypeCol},'SMP');
            for iRow=1:length(tempData{1})
                if strcmp(tempData{TypeCol}(iRow),'SMP')
                    theRow=find(strcmp(tempData{TimeCol}(iRow),theData{TimeCol}) & SMProws);
                    theData{stimulusCol}(theRow)=tempData{stimulusCol}(iRow);
                end;
            end;
        end;
        fileCounter=fileCounter+1;
    end;
    fprintf('%60s\n',' ' );
end;

%match up the EEG and SMI events and add the stimulus info to the events
matchTable=matchTable(~strcmp('none',matchTable(:,2)),:);%drop EEG events for which there are no matching SMI events
EEGevents={EPdataIn.events{1}.value}';
MSGrows=find(strcmp(theData{TypeCol},'MSG'));

EEGlist=[]; %list of relevant events in events array
SMIlist=[]; %list of corresponding rows in SMI event text
badSection=[];
for iMSG=1:length(MSGrows) %since there are missing triggers in the Brainvision data, will treat the SMI data as the more reliable record of events
    theMSG=MSGrows(iMSG);
    if any(strcmp(theData{msgCol}{theMSG},matchTable(:,2)))
        matchEEG=find(strcmp(matchTable{find(strcmp(theData{msgCol}{theMSG},matchTable(:,2))),1},EEGevents));
        if isempty(EEGlist)
            theEEGevent=matchEEG(1);
        else
            theEEGevent=matchEEG(min(find(matchEEG>EEGlist(end))));
        end;
        %patch for missing brainvision triggers.  Assuming the .vhmrk contents and the SMI events are in chronological order.  also assuming the first trigger is present and will error out if not
        if ~isempty(EEGlist)
            lastEEG=EPdataIn.events{1}(EEGlist(end)).sample*(1000/EPdataIn.Fs);
            newEEG=EPdataIn.events{1}(theEEGevent).sample*(1000/EPdataIn.Fs);
            lastSMI=str2num(theData{TimeCol}{SMIlist(end)})/1000;
            newSMI=str2num(theData{TimeCol}{theMSG})/1000;
%             if (newEEG-lastEEG)-(newSMI-lastSMI) > 10000 %if more than 10 sec discrepancy, EEG trigger was missing
%                 newSamp=round((str2num(theData{TimeCol}{theMSG})-str2num(theData{TimeCol}{SMIlist(end)}))/(1000*(1000/EPdataIn.Fs))+EPdataIn.events{1}(EEGlist(end)).sample);
%                 newEvent=min(find([EPdataIn.events{1}.sample]>newSamp));
%                 disp(['Missing trigger detected.  Discrepancy between SMI and EEG timing is: ' sprintf('%05.2f',((newEEG-lastEEG)-(newSMI-lastSMI))/1000) ' seconds.']);
%                 theEvent.type='Stimulus';
%                 theEvent.sample=newSamp;
%                 theEvent.value=matchTable{find(strcmp(theData{msgCol}{theMSG},matchTable(:,2))),1};
%                 theEvent.duration=1;
%                 theEvent.keys=struct('code','','data','','datatype','','description','');
%                 EPdataIn.events{1}=[EPdataIn.events{1}(1:newEvent-1) theEvent EPdataIn.events{1}(newEvent:end)];
%                 EEGevents={EPdataIn.events{1}.value}';
%                 theEEGevent=newEvent;
%                 badSection(end+1)=0;
            if (newEEG-lastEEG)-(newSMI-lastSMI) > 1000
                disp(['Substantial timing discrepancy for line ' num2str(theMSG) ': ' theData{msgCol}{theMSG} ' of ' sprintf('%05.2f',((newEEG-lastEEG)-(newSMI-lastSMI))/1000) ' seconds.']);
                badSection(end+1)=1;
            else
                badSection(end+1)=0;
            end;
        end;
        if isempty(theEEGevent)
            msg='Warning: EEG-SMI event mismatch.  Have run out of SMI events to match up with.  Will integrate SMI events up to this point.  Be sure to verify the event matching.';
            msgLog{end+1}=msg;
            disp(msg);
            EEGlist=[];
            break
        else
            SMIlist(end+1)=theMSG;
            EEGlist(end+1)=theEEGevent;
        end;
    end;
end;

badSection(end+1)=0; %if there was a problem with the timesynch in the last section I'd have no way of knowing it so will just assume it's okay.  Probably won't be using the data anyway.

sampLength=(1000/EPdataIn.Fs);

if ~isempty(EEGlist)    
    if length(SMIlist) ~= length(EEGlist)
        disp('Programming error');
        keyboard
    end;
    
    %add the SMI event information to the EEG dataset
    EPadd=[];
    EPadd.chanNames={'pupil';'x-eye';'y-eye'};
    EPadd.chanTypes={'PPL';'XEY';'YEY'};
    [EPdataIn]=ep_addData(EPdataIn,EPadd,'channels');
    if isempty(EPdataIn)
        return;
    end;
    
    %default the SMI channels to NaN for time points where no SMI data were collected
    EPdataIn.data(end-2,:,:,:,:,:,:)=NaN;
    EPdataIn.data(end-1,:,:,:,:,:,:)=NaN;
    EPdataIn.data(end,:,:,:,:,:,:)=NaN;
    
    %add an initial "event" corresponding to the start of the EEG or SMI (whichever began first) to avoid losing the SMI data prior to the first event.
    SMI_EEGoffset=(str2num(theData{TimeCol}{SMIlist(1)})/1000)-((EPdataIn.events{1}(EEGlist(1)).sample-1)*(1000/EPdataIn.Fs)); %SMI time when first matchTable EEG event started in ms
    firstSMI=min(find(strcmp('SMP',theData{TypeCol})));
    firstEEG=EPdataIn.events{1}(EEGlist(1)).sample*(1000/EPdataIn.Fs)+SMI_EEGoffset;
    
%since something in the experiment procedure is messing up the clock synchs at the start of each run, just drop the SMI data prior to the first matchList event.
%     if round((str2num(theData{TimeCol}{firstSMI})/1000)-SMI_EEGoffset) < 0
%         %SMI events prior to EEG
%         SMItimes=theData{TimeCol}{find(strcmp('SMP',theData{TypeCol}))};
%         SMIlist=[min(find((SMItimes/1000) >= firstEEG)) SMIlist];
%         if isempty(min(find((SMItimes/1000) >= firstEEG)))
%             disp('Problem')
%             keyboard
%         end;
%         EEGlist=[0 EEGlist]; %list of events.  There isn't one for the start of the EEG so zero
%     else
%         %SMI events after EEG
%         if firstSMI < SMIlist(1)
%             SMIlist=[firstSMI SMIlist];
%             EEGlist=[0 EEGlist];
%         end;
%     end;
    
    for iEvent=1:length(EEGlist)
        if iEvent==1
            SMI_EEGoffset=(str2num(theData{TimeCol}{SMIlist(2)})/1000)-((EPdataIn.events{1}(EEGlist(2)).sample-1)*(1000/EPdataIn.Fs)); %offset between SMI and EEG in ms.
        else
            SMI_EEGoffset=(str2num(theData{TimeCol}{SMIlist(iEvent)})/1000)-((EPdataIn.events{1}(EEGlist(iEvent)).sample-1)*(1000/EPdataIn.Fs)); %SMI time when EEG started in ms, recalculate in case of drift
        end;
        %Treat as saccade or blink if either eye designated as such, to be conservative with judging onsets of new state.
        %also, treating blanks as being something other than fixation or saccades
        startSMIsample=SMIlist(iEvent);
        if iEvent == length(EEGlist)
            endSMIsample=length(theData{1});
        else
            endSMIsample=SMIlist(iEvent+1);
        end;
        startEEGsample=EEGlist(iEvent);
        if iEvent == length(EEGlist)
            endEEGsample=EEGlist(iEvent); %will figure out later in the code
        else
            endEEGsample=EEGlist(iEvent+1);
        end;
        
        %add in the stimulus file for each event.  Assumes will be in an immediately prior event if present and ends in .jpg.
        if (EEGlist(iEvent) > 1) && (SMIlist(iEvent)>1)
            backRow=1;
            doneFlag=0;
            oneStimFlag=0;
            while ((SMIlist(iEvent)-backRow) > 0) && ~doneFlag
                theMsg=theData{msgCol}(SMIlist(iEvent)-backRow);
                if ~strcmp('MSG',theData{TypeCol}(SMIlist(iEvent)-backRow))
                    if ~oneStimFlag %skip over an SMP event in between the MSG events
                        oneStimFlag=1;
                        backRow=backRow+1;
                    else
                        doneFlag=1;
                    end;
                else
                    if strcmp('.jpg',theMsg{1}(end-3:end))
                        doneFlag=1;
                        oneEvent.code='stim';
                        oneEvent.data=theData{msgCol}{SMIlist(iEvent)-backRow};
                        if strfind(oneEvent.data,'# Message: ') == 1
                            oneEvent.data=oneEvent.data(12:end);
                        end;
                        oneEvent.datatype='char';
                        oneEvent.description='';
                        if isempty(EPdataIn.events{1}(EEGlist(iEvent)).keys(end).code)
                            EPdataIn.events{1}(EEGlist(iEvent)).keys(end)=oneEvent;
                        else
                            EPdataIn.events{1}(EEGlist(iEvent)).keys(end+1)=oneEvent;
                        end;
                        if exist([pathstr filesep 'stimuli'],'dir')
                            if exist([pathstr filesep 'stimuli' filesep oneEvent.data],'file') && ~any(strcmp(oneEvent.data,{EPdataIn.stims.name}))
                                EPdataIn.stims(end+1).name=oneEvent.data;
                                EPdataIn.stims(end).image=imread([pathstr filesep 'stimuli' filesep oneEvent.data]);
                            end;
                        end;
                    else
                        backRow=backRow+1;
                    end;
                end;
            end;
        end;
        
        %determine mouse click samples
        clickSamples=find(strcmp('# Message: UE-mouseclick',theData{msgCol}(startSMIsample:endSMIsample)));
        
        for iClicks=1:length(clickSamples)
            EPdataIn.events{1}(end+1).type='eye-tracker';
            EPdataIn.events{1}(end).sample=((str2num(theData{TimeCol}{clickSamples(iClicks)+startSMIsample-1})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs);
            EPdataIn.events{1}(end).value='clickET';
            EPdataIn.events{1}(end).duration=1;
            EPdataIn.events{1}(end).keys(1).code='location';
            EPdataIn.events{1}(end).keys(1).datatype='short';
            EPdataIn.events{1}(end).keys(1).data=theData{msgCol}{clickSamples(iClicks)+startSMIsample-1}(strfind(theData{msgCol}{clickSamples(iClicks)+startSMIsample-1},'x'):end);
            EPdataIn.events{1}(end).keys(1).description='';
            EPdataIn.events{1}(end).keys(2).code='mouse';
            EPdataIn.events{1}(end).keys(2).datatype='char';
            EPdataIn.events{1}(end).keys(2).data=theData{msgCol}{clickSamples(iClicks)+startSMIsample-1}(26:26);
            EPdataIn.events{1}(end).keys(2).description='';    
        end;
        
        %determine blink samples so can drop saccades and fixation events that are situated in the middle of a blink period
        blinkSamples=union(find(strcmp('Blink',theData{LeventCol}(startSMIsample:endSMIsample))),find(strcmp('Blink',theData{ReventCol}(startSMIsample:endSMIsample))));

        %determine if there is a prior sample so that first sample isn't automatically dubbed a saccade start if something started prior to it.
        SMPsamples=find(find(strcmp('SMP',theData{TypeCol}))<startSMIsample);
        priorSample=[];
        if ~isempty(SMPsamples)
            if ((str2num(theData{TimeCol}{startSMIsample})/1000) - (str2num(theData{TimeCol}{SMPsamples(end)})/1000)) <= 100
                priorSample=SMPsamples(end); 
            end;
        end;
        
        saccadeSamples=union(find(strcmp('Saccade',theData{LeventCol}(startSMIsample:endSMIsample))),find(strcmp('Saccade',theData{ReventCol}(startSMIsample:endSMIsample))));
        saccadeStarts=saccadeSamples(find(diff([0;saccadeSamples])>1));
        theEnds=find(diff([saccadeSamples; inf])>1)-1;
        if ~isempty(theEnds) && theEnds(1)==0
            theEnds(1)=1;
        end;
        saccadeEnds=saccadeSamples(theEnds)+1;
        %count as saccades start run of samples where either eye registers as being saccades.
        
        if ~isempty(saccadeStarts) && (saccadeStarts(1)==1)
            if isempty(priorSample) || (strcmp('Saccade',theData{LeventCol}(priorSample)) || strcmp('Saccade',theData{ReventCol}(priorSample)))
                %do not count first SMI sample in this section as being saccade if there was either a saccade just prior to it or if there was no SMI information at all just before it
                saccadeStarts=saccadeStarts(2:end);
            end;
        end;
        
        if any(ismember(saccadeStarts,blinkSamples))
            saccadeStarts(find(ismember(saccadeStarts,blinkSamples)))=[]; %drop saccade starts situated in the middle of a blink
        end;

        SMIsamples=startSMIsample:endSMIsample;
        keepSamples=find(strcmp('SMP',theData{TypeCol}(startSMIsample:endSMIsample))); %keep only SMP events, dropping MSG events
        keepSamples=setdiff(keepSamples,find(strcmp('Blink',theData{LeventCol}(startSMIsample:endSMIsample)))); %drop left eye blinks as eye-tracking will be unreliable
        keepSamples=setdiff(keepSamples,find(strcmp('Blink',theData{ReventCol}(startSMIsample:endSMIsample)))); %drop right eye blinks as eye-tracking will be unreliable
        keepSamples=setdiff(keepSamples,find(~strcmp('1',theData{LplaneCol}(startSMIsample:endSMIsample)))); %drop bad SMI samples according to L Plane
        keepSamples=setdiff(keepSamples,find(~strcmp('1',theData{RplaneCol}(startSMIsample:endSMIsample)))); %drop bad SMI samples according to R Plane
        SMIsamples=SMIsamples(keepSamples);
        if iEvent==1
            EEGsamples=[round(((str2num(theData{TimeCol}{SMIlist(1)})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs))+1:EPdataIn.events{1}(endEEGsample).sample];
        else
            EEGsamples=[EPdataIn.events{1}(startEEGsample).sample:EPdataIn.events{1}(endEEGsample).sample];
        end;
        
        SMItimesRaw=cellfun(@str2num,theData{TimeCol}(SMIsamples));
        SMItimes=(SMItimesRaw-SMItimesRaw(1))/1000; %may need to add an offset if can determine what it should be
        
        if iEvent == length(EEGlist)
            EEGsamples=[0:(SMItimes(end)/(1000/EPdataIn.Fs))]+EEGsamples(1);
            EEGtimes=(EEGsamples-EEGsamples(1))*(1000/EPdataIn.Fs);
            if SMItimes(end) < EEGtimes(end)
                excessSamples=length(find(EEGtimes > SMItimes(end)));
                EEGsamples=EEGsamples(1:end-excessSamples);
                EEGtimes=EEGtimes(1:end-excessSamples);
            end;
            if EEGsamples(end)>length(EPdataIn.timeNames)
                excessSamples=length(find(EEGsamples > length(EPdataIn.timeNames)));
                EEGsamples=EEGsamples(1:end-excessSamples);
                EEGtimes=EEGtimes(1:end-excessSamples);
            end;
        else
            EEGtimes=(EEGsamples-EEGsamples(1))*(1000/EPdataIn.Fs);
        end;
        
        if badSection(iEvent)
            EPdataIn.data(end-2,EEGsamples,:,:,:,:,:)=NaN;
            EPdataIn.data(end-1,EEGsamples,:,:,:,:,:)=NaN;
            EPdataIn.data(end,EEGsamples,:,:,:,:,:)=NaN;
        else
            for iSaccade=1:length(saccadeStarts)
                %calculate saccade direction based on 50 ms following each saccade
                sacTime=str2num(theData{TimeCol}{saccadeStarts(iSaccade)+startSMIsample-1});
                sacSamples=find((SMItimesRaw <= sacTime+50000) & (SMItimesRaw >= sacTime));
                sacSMIsamples=SMIsamples(sacSamples);
                sacX=cellfun(@str2num,theData{LXeyeCol}(sacSMIsamples))+cellfun(@str2num,theData{RXeyeCol}(sacSMIsamples));
                sacY=cellfun(@str2num,theData{LYeyeCol}(sacSMIsamples))+cellfun(@str2num,theData{RYeyeCol}(sacSMIsamples));
                Bx=[ones(length(sacX),1),SMItimes(sacSamples)]\sacX; %SMI increases rightward
                By=[ones(length(sacY),1),SMItimes(sacSamples)]\(-sacY); %SMI increases downward
                sacAng=atand(By(2)/Bx(2)); %will be 360 degrees with zero upwards
                
                if (Bx(2)==0) || (By(2)==0)
                    if (Bx(2)==0) && (By(2)==0)
                        sacAng=NaN;
                    elseif (Bx(2)==0) && (By(2)~=0)
                        if By(2)>0
                            sacAng=0;
                        else
                            sacAng=180;
                        end;
                    elseif (Bx(2)~=0) && (By(2)==0)
                        if Bx(2)>0
                            sacAng=90;
                        else
                            sacAng=270;
                        end;
                    end;
                else
                    if (Bx(2)>0) && (By(2)>0)
                        sacAng=90-sacAng;
                    elseif (Bx(2)<0) && (By(2)>0)
                        sacAng=270-sacAng;
                    elseif (Bx(2)>0) && (By(2)<0)
                        sacAng=90-sacAng;
                    elseif (Bx(2)<0) && (By(2)<0)
                        sacAng=270-sacAng;
                    end;
                end;
                
                EPdataIn.events{1}(end+1).type='eye-tracker';
                EPdataIn.events{1}(end).sample=((str2num(theData{TimeCol}{saccadeStarts(iSaccade)+startSMIsample-1})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs);
                EPdataIn.events{1}(end).value='saccadeET';
                EPdataIn.events{1}(end).duration=1;
                EPdataIn.events{1}(end).keys(1).code='SMIevent';
                EPdataIn.events{1}(end).keys(1).datatype='short';
                EPdataIn.events{1}(end).keys(1).data=num2str(iEvent);
                EPdataIn.events{1}(end).keys(1).description='';
                EPdataIn.events{1}(end).keys(2).code='word';
                EPdataIn.events{1}(end).keys(2).datatype='char';
                EPdataIn.events{1}(end).keys(2).data=theData{stimulusCol}{saccadeEnds(iSaccade)+startSMIsample-1};
                EPdataIn.events{1}(end).keys(2).description='';
                EPdataIn.events{1}(end).keys(3).code='angle';
                EPdataIn.events{1}(end).keys(3).datatype='char';
                EPdataIn.events{1}(end).keys(3).data=num2str(round(sacAng));
                EPdataIn.events{1}(end).keys(3).description='';
            end;
            fixationSamples=union(find(strcmp('Fixation',theData{LeventCol}(startSMIsample:endSMIsample))),find(strcmp('Fixation',theData{ReventCol}(startSMIsample:endSMIsample))));
            %count as fixation start run of samples where either eye registers as being fixations.
            fixationStarts=fixationSamples(find(diff([0;fixationSamples])>1));
            fixationEnds=fixationSamples(find(diff([fixationSamples; inf])>1));
            
            %         if ~isempty(fixationStarts) && (fixationStarts(1)==1)
            %             if isempty(priorSample) || (strcmp('Fixation',theData{LeventCol}(priorSample)) || strcmp('Fixation',theData{ReventCol}(priorSample)))
            %                 %do not count first SMI sample in this section as being fixation if there was either a fixation just prior to it or if there was no SMI information at all just before it
            %                 fixationStarts=fixationStarts(2:end);
            %             end;
            %         end;
            
            if any(ismember(fixationStarts,blinkSamples))
                fixationStarts(find(ismember(fixationStarts,blinkSamples)))=[]; %drop fixation starts situated in the middle of a blink
            end;
            
            for iFixation=1:length(fixationStarts)
                EPdataIn.events{1}(end+1).type='eye-tracker';
                EPdataIn.events{1}(end).sample=((str2num(theData{TimeCol}{fixationStarts(iFixation)+startSMIsample-1})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs);
                EPdataIn.events{1}(end).value='fixationET';
                EPdataIn.events{1}(end).duration=1;
                EPdataIn.events{1}(end).keys(1).code='SMIevent';
                EPdataIn.events{1}(end).keys(1).datatype='short';
                EPdataIn.events{1}(end).keys(1).data=num2str(iEvent);
                EPdataIn.events{1}(end).keys(1).description='';
                EPdataIn.events{1}(end).keys(2).code='word';
                EPdataIn.events{1}(end).keys(2).datatype='char';
                EPdataIn.events{1}(end).keys(2).data=theData{stimulusCol}{fixationEnds(iFixation)+startSMIsample-1};
                EPdataIn.events{1}(end).keys(2).description='';
            end;
            %count as blink start run of samples where either eye registers as being blinks.
            blinkStarts=blinkSamples(find(diff([-inf;blinkSamples])>1));
            theEnds=find(diff([blinkSamples; inf])>1)-1;
            if ~isempty(theEnds) && theEnds(1)==0
                theEnds(1)=1;
            end;
            blinkEnds=blinkSamples(theEnds)+1;
            
            if ~isempty(blinkStarts) && (blinkStarts(1)==1)
                if isempty(priorSample) || (strcmp('Blink',theData{LeventCol}(priorSample)) || strcmp('Blink',theData{ReventCol}(priorSample)))
                    %do not count first SMI sample in this section as being blink if there was either a blink just prior to it or if there was no SMI information at all just before it
                    blinkStarts=blinkStarts(2:end);
                end;
            end;
            
            for iBlink=1:length(blinkStarts)
                EPdataIn.events{1}(end+1).type='eye-tracker';
                EPdataIn.events{1}(end).sample=((str2num(theData{TimeCol}{blinkStarts(iBlink)+startSMIsample-1})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs);
                EPdataIn.events{1}(end).value='blinkStartET';
                EPdataIn.events{1}(end).duration=1;
                EPdataIn.events{1}(end).keys(1).code='SMIevent';
                EPdataIn.events{1}(end).keys(1).datatype='short';
                EPdataIn.events{1}(end).keys(1).data=num2str(iEvent);
                EPdataIn.events{1}(end).keys(1).description='';
                EPdataIn.events{1}(end).keys(2).code='word';
                EPdataIn.events{1}(end).keys(2).datatype='char';
                EPdataIn.events{1}(end).keys(2).data=theData{stimulusCol}{blinkStarts(iBlink)+startSMIsample-1};
                EPdataIn.events{1}(end).keys(2).description='';
            end;
            for iBlink=1:length(blinkEnds)
                EPdataIn.events{1}(end+1).type='eye-tracker';
                EPdataIn.events{1}(end).sample=((str2num(theData{TimeCol}{blinkEnds(iBlink)+startSMIsample-1})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs);
                EPdataIn.events{1}(end).value='blinkEndET';
                EPdataIn.events{1}(end).duration=1;
                EPdataIn.events{1}(end).keys(1).code='SMIevent';
                EPdataIn.events{1}(end).keys(1).datatype='short';
                EPdataIn.events{1}(end).keys(1).data=num2str(iEvent);
                EPdataIn.events{1}(end).keys(1).description='';
                EPdataIn.events{1}(end).keys(2).code='word';
                EPdataIn.events{1}(end).keys(2).datatype='char';
                EPdataIn.events{1}(end).keys(2).data=theData{stimulusCol}{blinkEnds(iBlink)+startSMIsample-1};
                EPdataIn.events{1}(end).keys(2).description='';
            end;
            
            pupilData=(cellfun(@str2num,theData{LpupilCol}(SMIsamples))+cellfun(@str2num,theData{RpupilCol}(SMIsamples)))/2;
            
            %drop SMI samples where one eye has zeroes and then take the mean of the two eye measures
            X1=cellfun(@str2num,theData{LXeyeCol}(SMIsamples));
            X1samples=X1~=0;
            X2=cellfun(@str2num,theData{RXeyeCol}(SMIsamples));
            X2samples=X2~=0;
            XeyeData=(X1+X2)./(X1samples+X2samples);
            
            Y1=cellfun(@str2num,theData{LYeyeCol}(SMIsamples));
            Y1samples=Y1~=0;
            Y2=cellfun(@str2num,theData{RYeyeCol}(SMIsamples));
            Y2samples=Y2~=0;
            YeyeData=(Y1+Y2)./(Y1samples+Y2samples);
            
            %resample SMI data since its sampling rate is rather variable
            pupilData=pupilData(:);
            XeyeData=XeyeData(:);
            YeyeData=YeyeData(:);
            SMItimes=SMItimes(:);
            
            if length(EEGtimes) > 1
                if EEGtimes(end) > SMItimes(end)
                    theTimes=[SMItimes; EEGtimes(end)];
                    pupilData=interp1(theTimes,[pupilData; pupilData(end)],EEGtimes);
                    XeyeData=interp1(theTimes,[XeyeData; XeyeData(end)],EEGtimes);
                    YeyeData=interp1(theTimes,[YeyeData; YeyeData(end)],EEGtimes);
                else
                    theTimes=SMItimes;
                    pupilData=interp1(theTimes,pupilData,EEGtimes);
                    XeyeData=interp1(theTimes,XeyeData,EEGtimes);
                    YeyeData=interp1(theTimes,YeyeData,EEGtimes);
                end;
            end;
            
            blinkStarts=blinkSamples(find(diff([-inf;blinkSamples])>1)); %start of blinks in SMI samples relative to the section
            blinkEnds=blinkSamples(find(diff([blinkSamples;inf])>1)); %end of blinks in SMI samples relative to the section
            sampleOffset=round(((str2num(theData{TimeCol}{startSMIsample})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs)); %offset in samples between start of section and sample times to convert to section relative samples
            for iBlink=1:length(blinkStarts)
                EEGsampleStart=round(((str2num(theData{TimeCol}{blinkStarts(iBlink)+startSMIsample-1})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs)-sampleOffset+1);
                EEGsampleEnd=round(((str2num(theData{TimeCol}{blinkEnds(iBlink)+startSMIsample-1})/1000)-SMI_EEGoffset)/(1000/EPdataIn.Fs)-sampleOffset+1);
                if EEGsampleEnd > length(EEGsamples)
                    disp('blink sample too long for section.  Truncating.');
                    EEGsampleEnd=length(EEGsamples);
                end;
                pupilData(EEGsampleStart:EEGsampleEnd)=NaN;
                XeyeData(EEGsampleStart:EEGsampleEnd)=NaN;
                YeyeData(EEGsampleStart:EEGsampleEnd)=NaN;
            end;
            
            pausePoints=find(diff(SMItimes) > 50); %pauses of over 50 ms are treated as NaN gaps rather than interpolated
            for iPause=1:length(pausePoints)
                EEGsampleStart=round(SMItimes(pausePoints(iPause))/sampLength)+1;
                EEGsampleEnd=EEGsampleStart+round((SMItimes(pausePoints(iPause)+1)-SMItimes(pausePoints(iPause)))/sampLength)-2;
                if EEGsampleEnd > length(EEGsamples)
                    disp('pause sample too long for section.  Truncating.');
                    EEGsampleEnd=length(EEGsamples);
                end;
                pupilData(EEGsampleStart:EEGsampleEnd)=NaN;
                XeyeData(EEGsampleStart:EEGsampleEnd)=NaN;
                YeyeData(EEGsampleStart:EEGsampleEnd)=NaN;
            end;
            
            EPdataIn.data(end-2,EEGsamples,:,:,:,:,:)=pupilData;
            EPdataIn.data(end-1,EEGsamples,:,:,:,:,:)=XeyeData;
            EPdataIn.data(end,EEGsamples,:,:,:,:,:)=YeyeData;
        end;
    end;
    
    EPdataOut=EPdataIn;
    disp(['The subject took ' sprintf('%5.2f',toc/60) ' minutes to load in the SMI information.']);
end;