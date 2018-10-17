function eventLines=ep_collateEventLines(theFunction)
%  eventLines=ep_collateEventLines(theFunction)
%
%generate list of event lines for View functions.
%
%Inputs:
%   theFunction : The calling function ('EPwaves' or 'EPtopos') to determine format of output.
%
%Outputs:
%   eventLines  : Cell array of samples of events to be displayed (4 colors)(cells/subs)

%History:
%  by Joseph Dien (6/7/18)
%  jdien07@mac.com
%
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

global EPdataset EPmain

eventLines=cell(4,1);
if strcmp('VLT',EPmain.view.dataTransform)
    if (isempty(EPmain.view.eventList) && EPmain.view.events==1) || (~isempty(EPmain.view.eventList) && (EPmain.view.events <= length(EPmain.view.eventList)+2))
        for iColor=1:4
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                subList=1;
                if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'continuous')
                    cellList=1;
                end;
                if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average')
                    if any(ismember(EPmain.view.allTrials(iColor),[5 6]))
                        cellList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames)];
                    else
                        cellList=EPmain.view.cell(iColor);
                    end;
                    if any(ismember(EPmain.view.allTrials(iColor),[1 2]))
                        subList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames)];
                    else
                        subList=EPmain.view.subject(iColor);
                    end;
                else  %if single_trial data
                    if any(ismember(EPmain.view.allTrials(iColor),[1 2]))
                        cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
                    else
                        cellList=intersect(find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames)),...
                            find(EPmain.view.trial(iColor)==EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames));
                    end;
                end;
                if strcmp('EPtopos',theFunction) || any(ismember(EPmain.view.allTrials(iColor),[2 4 6]))
                    eventLines{iColor}=cell(max(length(cellList),length(subList)),1);
                end;
                for iSub=1:length(subList)
                    theSub=subList(iSub);
                    for iCell=1:length(cellList)
                        theCell=cellList(iCell);
                        if (EPmain.view.events > length(EPmain.view.eventList)) && EPmain.view.RT
                            theRTspec=find(strcmp('RT',EPdataset.dataset(EPmain.view.dataset(iColor)).trialSpecNames));
                            if ~isempty(theRTspec)
                                theRT=EPdataset.dataset(EPmain.view.dataset(iColor)).trialSpecs{theCell,theRTspec,theSub};
                                if ~isempty(theRT)
                                    if ischar(theRT)
                                        theRT=str2double(theRT);
                                    end;
                                    theRT=ceil(theRT/(1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs))+EPdataset.dataset(EPmain.view.dataset(iColor)).baseline;
                                    if (theRT > 0) && (theRT <= length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames))
                                        if strcmp('EPtopos',theFunction) || any(ismember(EPmain.view.allTrials(iColor),[2 4 6]))
                                            eventLines{iColor}{max(iSub,iCell)}(end+1)=theRT;
                                        else
                                            eventLines{iColor}(end+1)=theRT;
                                        end;
                                    end;
                                end;
                            end;
                        end;
                        if (EPmain.view.events == length(EPmain.view.eventList)+2) || (EPmain.view.events <= length(EPmain.view.eventList))
                            theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{theSub,theCell};
                            for iEvent=1:length(theEvents)
                                if  (EPmain.view.events == length(EPmain.view.eventList)+2) || (strcmp([num2str(theEvents(iEvent).type) '-' num2str(theEvents(iEvent).value)],EPmain.view.eventList{EPmain.view.events}))
                                    if strcmp('EPtopos',theFunction) || any(ismember(EPmain.view.allTrials(iColor),[2 4 6]))
                                        eventLines{iColor}{max(iSub,iCell)}(end+1)=theEvents(iEvent).sample;
                                    else
                                        eventLines{iColor}(end+1)=theEvents(iEvent).sample;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
                %                     eventLines(iColor)=length(setdiff(eventNames,{'SESS','CELL','TRSP','bgin','boundary'}));
            end;
        end;
    end;
end;