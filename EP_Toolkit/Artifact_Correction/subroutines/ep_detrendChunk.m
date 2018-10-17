function [outputLog, graphCounter] = ep_detrendChunk(inFile, startChunk, endChunk, theSubject, doDetrend, baseline, butterflyFig, graphCounter, numGraphs);
% [outputLog] = ep_detrendChunk(inFile, startChunk, endChunk, theSubject);
%	Detrends and baseline corrects the data on a trialwise basis, reading in a chunk at a time.
%   Baseline correction performed after the detrending.
%
%Inputs
%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.
%	startChunk: starting chunk (usually 1)
%   endChunk:   ending chunk
%   theSubject: which subject of the file is being processed.
%   doDetrend:  whether to perform detrend (0=no, 1=yes)
%   baseline:   whether to perform baseline correction (empty=no, array of numbers=yes and what points to use)
%   butterflyFig:  the handle for the output figure.  Otherwise, will open a new figure.
%   graphCounter: the current subplot for the summary figure.
%   numGraphs: the total number of subgraphs in the summary figure.
%
%   The input chunks include: dataChunk in EP data format.
%
%Outputs
%   outputLog: output messages from bad data process.
%   graphCounter: the current subplot for the summary figure.
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
% modified 2/11/10 JD
% Will now work with subject average files with multiple subjects.
%
% modified 1/28/13 JD
% Added baseline correction.
%
% bugfix 5/9/13 JD
% Fixed detrend not working with a single chunk and crashing with multiple chunks.
%
% bugfix 10/18/17 JD
% Detrend and baseline correct are now applied per trial for single-trial and average data rather than applied as if it was a continuous dataset.
%
% modified 2/4/18 JD
% Added support for summary figure output.
%
% bugfix & modified 4/30/18 JD
% Fixed turning data to NaN when no baseline set and data are not continuous.
% Consolidated summary figure for average files so no longer one per subject.
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

if doDetrend && isempty(baseline)
    msg=['Detrending.'];
elseif ~doDetrend && ~isempty(baseline)
    msg=['Baseline Correcting.'];
elseif doDetrend && ~isempty(baseline)
    msg=['Detrending and Baseline Correcting.'];
else
    return
end;
disp(msg);
outputLog{1}=msg;

if ~exist('butterflyFig','var')
    butterflyFig=figure('Name','Artifact Correction','NumberTitle','off');
    colormap jet;
    standAlone=1;
else
    standAlone=0;
end;

for iChunk = startChunk:endChunk
    eval(['load ''' deblank(inFile) '''-' num2str(iChunk) '.mat']);
    numChans=length(dataChunk.chanNames);
    EEGchans=find(strcmp('EEG',dataChunk.chanTypes));
    numPoints = length(dataChunk.timeNames);
    if strcmp(dataChunk.dataType,'continuous')
        theData=squeeze(dataChunk.data(:,:,:,theSubject));
        if ~isempty(baseline)
            baseMeans = mean(theData(EEGchans,baseline),2);
            theData(EEGchans,:)=theData(EEGchans,:)-repmat(baseMeans,1,numPoints);
        end;
        if doDetrend
            theData=detrend(theData')';
        end;
        dataChunk.data(:,:,:,theSubject)=theData;
    else
        numTrials = length(dataChunk.cellNames);
        for iTrial = 1:numTrials
            theData=squeeze(dataChunk.data(:,:,iTrial,theSubject));
            if ~isempty(baseline)
                baseMeans = mean(theData(EEGchans,baseline),2);
                theData(EEGchans,:)=theData(EEGchans,:)-repmat(baseMeans,1,numPoints);
            end;
            if doDetrend
                theData=detrend(theData')';
            end;
            dataChunk.data(:,:,iTrial,theSubject)=theData;
        end;
    end;
    eval (['save ''' inFile '''-' num2str(iChunk) '.mat dataChunk']);
    if ~isempty(butterflyFig)
        trialdata=reshape(dataChunk.data(:,:,:,theSubject),numChans,[]);
        displayPeriod=size(trialdata,2);    %Number of timepoints to graph in display.
        totalDisplayPeriod=displayPeriod*size(dataChunk.data,4);
        decimateSamples=ceil(max(1,totalDisplayPeriod/10000));     
        figure(butterflyFig(iChunk));
%         axesHandles=get(butterflyFig(iChunk),'children');
%         theAxes=0;
%         for iAxes=1:length(axesHandles)
%             theTitle=get(axesHandles(iAxes),'title');
%             if any(strcmp(theTitle.String,{'Detrended','Baseline Corrected','Detrended and Baseline Corrected.'}))
%                 theAxes=iAxes;
%             end;
%         end;
%         plotData=zeros(numChans,length([1:decimateSamples:totalDisplayPeriod]));
%         plotPoints=[1:decimateSamples:totalDisplayPeriod];
%         if theAxes
%             lineHandles=get(axesHandles(theAxes),'children');
%             for iLine=1:length(lineHandles)
%                 plotData(iLine,:)=get(lineHandles(iLine),'YData');
%             end;
%         end;
%         subjectPoints=find(ismember([1:displayPeriod]+(theSubject-1)*displayPeriod,plotPoints));
%         priorPoints=length(find(plotPoints<((theSubject-1)*displayPeriod+1)));
        if doDetrend && isempty(baseline)
            theTitle='Detrended';
        elseif ~doDetrend && ~isempty(baseline)
            theTitle='Baseline Corrected';
        elseif doDetrend && ~isempty(baseline)
            theTitle='Detrended and Baseline Corrected.';
        else
            return
        end;
%         plotData(EEGchans,[1:length(subjectPoints)]+priorPoints)=trialdata(EEGchans,subjectPoints);
        plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,trialdata,EEGchans,theSubject);
        subplot(numGraphs,1,graphCounter), plot([1:decimateSamples:totalDisplayPeriod],plotData);
        title(theTitle,'Interpreter','none');
        axis([1 totalDisplayPeriod -200 200])
        set(gca,'XTickLabel','','XTick',[]);
        
        drawnow
    end;
    if standAlone
        try
            MATLABver=ver('MATLAB');
            [a b]=strtok(MATLABver.Version,'.');
            b=b(2:end);
            if ~isprop(butterflyFig,'Number')
                eval (['print -f' num2str(butterflyFig(iChunk)) ' -djpeg ''' inFile '''-' num2str(iChunk) 'detrend.jpg']);
            else
                eval (['print -f' num2str(butterflyFig(iChunk).Number) ' -djpeg ''' inFile '''-' num2str(iChunk) 'detrend.jpg']);
            end;
        catch
            disp('Couldn''t save a copy of the artifact correction figure.  Perhaps your version of Matlab is not current.');
        end;
        close(butterflyFig(iChunk));
    end;
end;
if ~isempty(butterflyFig)
    graphCounter=graphCounter+1;
end;