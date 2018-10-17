function plotData=ep_makePlotData(figureHandle,displayPeriod,totalDisplayPeriod,decimateSamples,titleList,trialdata,EEGchans,theSubject) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plotData=makePlotData(figureHandle,displayPeriod,totalDisplayPeriod,decimateSamples,titleList,trialdata,EEGchans,theSubject) 
%
%	%decimate the data to be appropriate for the summary figure plots
%
%Inputs
%   figureHandle:  the handle for the output figure.
%	displayPeriod:     The undecimated length of the data to be displayed (for one subject if average data).
%	totalDisplayPeriod: The decimated length of the data to be displayed (for all the subjects if average data).
%   decimateSamples:   number of samples to skip when decimating the data
%   titleList:  THe title of the subplot to be updated.
%   trialdata:   The new data (chans, timepoints)
%   EEGchans:    Array of EEG channels in the input data.
%   theSubject: which subject of the file is being processed.
%
%Outputs
%	plotData: decimated data for plotting in the summary figure (channels, timepoints).
%
% History:
%
% by Joseph Dien (4/8/18)
% jdien07@mac.com
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

figure(figureHandle);
axesHandles=get(figureHandle,'children');
theAxes=0;
for iAxes=1:length(axesHandles)
    theTitle=get(axesHandles(iAxes),'title');
    if any(strcmp(theTitle.String,titleList))
        theAxes=iAxes;
    end;
end;
numChans=length(EEGchans);
plotData=zeros(numChans,length([1:decimateSamples:totalDisplayPeriod]));
plotPoints=[1:decimateSamples:totalDisplayPeriod];
if theAxes
    lineHandles=get(axesHandles(theAxes),'children');
    for iLine=1:length(lineHandles)
        plotData(iLine,:)=get(lineHandles(iLine),'YData');
    end;
end;
subjectPoints=find(ismember([1:displayPeriod]+(theSubject-1)*displayPeriod,plotPoints));
priorPoints=length(find(plotPoints<((theSubject-1)*displayPeriod+1)));
plotData(:,[1:length(subjectPoints)]+priorPoints)=trialdata(EEGchans,subjectPoints);

