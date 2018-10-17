function [ALLEEG]=ep_ep2alleeg(EPdata)
% ep_ep2alleeg - [ALLEEG]=ep_ep2alleeg(EPdata) -
% Convert EPdata format to EEGlab's ALLEEG format.  Doesn't convert factor files.
%
%Inputs
%   EPdata        : Structured array with the data and accompanying information in EP file format.  See readData.
%
%Outputs:
%  ALLEEG         : EEGlab's ALLEG format.
%
%History
%  by Joseph Dien, with help by Grega Repovs (3/28/09)
%  jdien07@mac.com
%
%  EEGlab files don't really have set conventions but to the extent that there is one,
%  the events should have both a 'value' field to denote the generic type of event,
%  as in 'trigger', and a 'type' field to denote the nature of this generic event,
%  as in the condition of the experiment.
%  Note also that this is the reverse of the FieldTrip convention.


%
% bugfix 5/24/10 JD
% Now separates each cell into a separate dataset, not just each subject.
%
% bugfix 10/20/13 JD
% Fixed crash when saving single-trial data.
% Fixed including only first trial of each cell type.
% Fixed not including correct condition name.
% Fixed not setting up events field correctly.
% Fixed not filling out ref field.
% Fixed epoch field to correspond to epoch from entire session, not from just the one condition.
% Fixed not filling out rejected field based on bad channel and bad trial fields.
%
% bugfix 10/29/13 JD
% Fixed crash when translating single-trial EP file to EEGlab format and there are no events.
%
% modified 10/29/13 JD
% Added .keys field to events.
%
% bugfix 6/30/14 JD
% Fixed crash when file has an empty .keys field.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bugfix 7/21/14 JD
% Fixed labels for REG channels being blank.
% For single_trial data, the event latency values now conform to "for epoched datasets the event latencies are also
% encoded in sample points with respect to the beginning of the data (as if the data were continuous)"
% . http://sccn.ucsd.edu/wiki/Chapter_03:_Event_Processing rather than being in terms of the beginning of the original continuous data.
% Fixed adding 'trigger' events if they are already present.
% Fixed epoch event field referring to trial number from complete dataset rather than in terms of the EEG file (the one
% condition).
% Fixed urevent event field reflecting event numbering of full dataset rather than just the one condition in the EEG
% file.
% Fixed sometimes adding too many epoch entries when exporting single_trial .set files, resulting in aborted export process.
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

ALLEEG=[];

uniqueCells=unique(EPdata.cellNames);
numPoints=length(EPdata.timeNames);

for iEloc=1:length(EPdata.eloc)
    if isempty(EPdata.eloc(iEloc).labels)
        EPdata.eloc(iEloc).labels=EPdata.chanNames{iEloc};
    end;
end;

count=0;
for subject=1:length(EPdata.subNames)
    for theCell=1:length(uniqueCells)
        count=count+1;
        
        ALLEEG(count).setname=[EPdata.ename '-' EPdata.subNames{subject} '-' EPdata.cellNames{theCell}];
        ALLEEG(count).filename=[];
        ALLEEG(count).filepath=[];
        ALLEEG(count).subject=EPdata.subNames{subject};
        ALLEEG(count).group='';
        ALLEEG(count).condition=uniqueCells{theCell};
        ALLEEG(count).session=[];
        ALLEEG(count).comments=[];
        ALLEEG(count).nbchan=length(EPdata.chanNames);
        ALLEEG(count).trials=length(strmatch(uniqueCells{theCell},EPdata.cellNames,'exact'));
        ALLEEG(count).pnts=length(EPdata.timeNames);
        ALLEEG(count).srate=EPdata.Fs;
        ALLEEG(count).xmin=min(EPdata.timeNames)/1000;
        ALLEEG(count).xmax=max(EPdata.timeNames)/1000;
        ALLEEG(count).times=EPdata.timeNames;
        ALLEEG(count).data=squeeze(EPdata.data(:,:,strmatch(uniqueCells{theCell},EPdata.cellNames,'exact'),subject));
        ALLEEG(count).icaact=[];
        ALLEEG(count).icawinv=[];
        ALLEEG(count).icasphere=[];
        ALLEEG(count).icaweights=[];
        ALLEEG(count).icachansind=[];
        ALLEEG(count).chanlocs=EPdata.eloc;
        ALLEEG(count).urchanlocs=[];
        ALLEEG(count).chaninfo.icachansind=[];
        if strcmp(EPdata.reference.type,'REG')
            ALLEEG(count).ref='common';
        elseif strcmp(EPdata.reference.type,'AVG')
            ALLEEG(count).ref='averef';
        else
            ALLEEG(count).ref=[];
        end;
        
        trialList=strmatch(uniqueCells{theCell},EPdata.cellNames,'exact');
        clear events;
        if ~isempty(EPdata.events)
            for trial=1:length(trialList)
                theTrial=trialList(trial);
                for theEvent=1:length(EPdata.events{subject,theTrial})
                    EPdata.events{subject,theTrial}(theEvent).epoch=trial;
                    EPdata.events{subject,theTrial}(theEvent).sample=(trial-1)*ALLEEG(count).pnts+EPdata.events{subject,theTrial}(theEvent).sample;
                    if exist('events', 'var')
                        events(end+1)=EPdata.events{subject,theTrial}(theEvent);
                    else
                        events(1)=EPdata.events{subject,theTrial}(theEvent);
                    end;
                end;
            end;
            if exist('events', 'var')
                tempEvents=events;
                for i=1:length(events)
                    events(i).value=tempEvents(i).type;
                    events(i).type=tempEvents(i).value;
                    events(i).latency=events(i).sample;
                    events(i).urevent=i;
                    for iKey=1:length(events(i).keys)
                        keyName=events(i).keys(iKey).code;
                        if ~isempty(keyName)
                            if ~any(strcmp(keyName,{'epoch','urevent'}))
                                if ~isempty(strfind(keyName,'#'))
                                    keyName=strrep(keyName,'#','hash_');
                                end;
                                if ~isempty(strfind(keyName,'+'))
                                    keyName=strrep(keyName,'+','plus_');
                                end;
                                eval(['events(i).' keyName '=events(i).keys(iKey).data;']);
                            end;
                        end;
                    end;
                end;
                events = rmfield(events,'sample');
                events = rmfield(events,'keys');
            end;
        end;
        
        if ~exist('events', 'var')
            events(1)=struct('type',[],'value',[],'latency',[],'duration',[],'epoch',[]);
        end;
        
        %add trigger events to mark each trial if single-trial data.
        if strcmp(EPdata.dataType,'single_trial')
            for trial=1:length(trialList)
                theTrial=trialList(trial);
                if length(find(strcmp('trigger',{events.value}))) ~= length(trialList)
                    if isempty(events(end).value) && length(events)==1
                        events(end).value='trigger';
                    else
                        events(end+1).value='trigger';
                    end;
                    events(end).type=EPdata.cellNames{theTrial};
                    events(end).latency=(trial-1)*ALLEEG(count).pnts+EPdata.baseline+1;
                    events(end).duration=0;
                    events(end).epoch=trial;
                end;
                
                if ~isempty(EPdata.trialSpecNames)
                    events(end+1).value='TRSP';
                    events(end).type='';
                    events(end).latency=(trial-1)*ALLEEG(count).pnts+EPdata.baseline+1;
                    events(end).duration=0;
                    events(end).epoch=trial;
                    events(end).keys.names=EPdata.trialSpecNames;
                    events(end).keys.specs=EPdata.trialSpecs(theTrial,:);
                end;
            end;
        end;
        
        ALLEEG(count).event=events;
        ALLEEG(count).urevent=events;
        
        ALLEEG(count).eventdescription=[];
        
        clear epoch;
        epoch(1)=struct('event',[],'eventtype',[],'eventvalue',[],'eventlatency',[],'eventduration',[]);
        
        for trial=1:length(trialList)
            eventList=find([events.epoch]==trial);
            epoch(trial).event=eventList;
            epoch(trial).eventduration={events(eventList).duration};
            epoch(trial).eventlatency={events(eventList).latency};
            epoch(trial).eventtype={events(eventList).type};
            epoch(trial).eventvalue={events(eventList).value};
            if strcmp(EPdata.dataType,'single_trial')
                ALLEEG(count).reject.manualE=(squeeze(EPdata.analysis.badChans(subject,trialList,:))' ~=0);
                ALLEEG(count).reject.manual=EPdata.analysis.badTrials(subject,trialList)~=0;
            else
                ALLEEG(count).reject.manualE=false(size(squeeze(EPdata.analysis.badChans(subject,trialList,:))'));
                ALLEEG(count).reject.manual=false(EPdata.analysis.badTrials(subject,trialList));
            end;
        end;
        
        ALLEEG(count).epoch=epoch;
        ALLEEG(count).epochdescription=cell(0);
        ALLEEG(count).stats=[];
        ALLEEG(count).specdata=[];
        ALLEEG(count).specicaact=[];
        ALLEEG(count).splinefile='';
        ALLEEG(count).icasplinefile='';
        ALLEEG(count).dipfit=[];
        ALLEEG(count).history=[];
        ALLEEG(count).saved='yes';
        ALLEEG(count).etc=[];
        ALLEEG(count).datfile=[];
    end;
end;