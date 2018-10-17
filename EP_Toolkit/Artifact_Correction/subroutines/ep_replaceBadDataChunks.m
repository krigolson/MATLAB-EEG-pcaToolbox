function [outputLog, graphCounter] = ep_replaceBadDataChunks(inFile, startChunk, endChunk, badChans, butterflyFig, graphCounter, numGraphs, theSubject);
% [outputLog, graphCounter] = ep_replaceBadDataChunks(inFile, startChunk, endChunk, badChans, butterflyFig, graphCounter, numGraphs, theSubject);
%	Interpolates bad channels.
%
%Inputs
%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.
%	startChunk: starting chunk (usually 1)
%   endChunk:   ending chunk
%   badChans:   list of globally bad channels.  Will be set to a flat line.
%   butterflyFig:  the handle for the output figure.  Otherwise, will open a new figure.
%   graphCounter: the current subplot for the summary figure.
%   numGraphs: the total number of subgraphs in the summary figure.
%   theSubject: which subject of the file is being processed.
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
% by Joseph Dien (4/9/09)
% jdien07@mac.com
%
% Based on markBadDataChunks.
%
% modified 4/17/09 JD
% Dropped eloc as separate input parameter (now part of data).
%
% modified 5/30/09 JD
% Ouputs to butterly figure.
%
% bugfix 7/14/09 JD
% Fixed occasional crash in bad channel replacement code.
%
%  bugfix 9/14/09 JD
%  Don't bother to go through the time points of a trial if none of the channels are bad.  Was slowing things down.
%
% modified 10/28/09 JD
% Added option to disable preprocessing figure for low memory situations.
%
%  bugfix 12/8/09 JD
%  Bad channels incorrectly interpolated (the x & y coordinates were reversed).
%
% modified 2/7/10 JD
% Changed bad channel field to negative numbers for still bad channels.
%
% modified 2/11/10 JD
% Will now work with subject average files with multiple subjects.
% analysis fields no longer optional.
% Gets bad channel and bad trial info from the data chunk rather than from the function call.
%
%  modified 10/12/10 JD
%  For continuous files, data now divided into one second epochs and can be artifact rejected in an epochwise fashion
%  in same fashion as segmented data.
% 
%  modified 10/16/10 JD
%  Added support for HEOG saccade correction.
% 
% modified 10/10/13 JD
% Restricted bad channel correction to EEG channels.
%
% bufix 3/11/14 JD
% Handles decimal sampling rates gracefully.
%
% bufix 4/16/14 JD
% Fixed crash when replacing bad channels in an average file and a cell is bad.
% 
% modified 4/24/14 JD
% Rereference an epoch after performing bad channel replacement
%
% bufix 6/6/14 JD
% Fixed crash when running replace bad channels in continuous data.
%
% bugfix 8/14/15 JD
% Fixed unable to save artifact correction summary figure starting with
% Matlab 2014b in stand-alone mode.
% 
% modified 9/29/16 JD
% eye-tracker channels set to NaN for bad trials.
%
% bugfix 10/5/17 JD
% Fixed crashes due to Matlab changing their graphics objects AGAIN in 2017b.
%
% bugfix 10/20/17 JD
% Eliminated x tick labels to address problem with subplots in summary
% artifact figure getting squeezed by formatting problem on PCs.
%
% modified 2/4/18 JD
% Made subplot specification for summary figure output more flexible.
%
% modified 4/8/18 JD
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

msg='Replacing bad channels.';
disp(msg);
outputLog{1}=msg;

if ~exist('butterflyFig','var')
    butterflyFig=figure('Name','Bad Channel Correction','NumberTitle','off');
    colormap jet;
    standAlone=1;
else
    standAlone=0;
end;

for iChunk = startChunk:endChunk
    disp([deblank(inFile) '-' num2str(iChunk)]);
    warning off MATLAB:griddata:DuplicateDataPoints
    fprintf('%60s\n',' ' );
    trialCount=1;
    eval(['load ''' deblank(inFile) '-' num2str(iChunk) '.mat''']);
    
    if length(dataChunk.facNames) > 1
        error('This function is not intended for application to factor data.');
    end;
    
    maxRad=0.5;
    GRID_SCALE=67;
    [y,x] = pol2cart(([dataChunk.eloc.theta]/360)*2*pi,[dataChunk.eloc.radius]);  % transform electrode locations from polar to cartesian coordinates
    y=-y; %flip y-coordinate so that nose is upwards.
    plotrad = min(1.0,max([dataChunk.eloc.radius])*1.02);            % default: just outside the outermost electrode location
    plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
    x = x*(maxRad/plotrad);
    y = y*(maxRad/plotrad);
    
    xmin = min(-maxRad,min(x));
    xmax = max(maxRad,max(x));
    ymin = min(-maxRad,min(y));
    ymax = max(maxRad,max(y));
    
    x=round(((x/(xmax-xmin))*GRID_SCALE)+ceil(GRID_SCALE/2));
    y=round(((y/(ymax-ymin))*GRID_SCALE)+ceil(GRID_SCALE/2));

    if strcmp(dataChunk.dataType,'continuous')
        numTrials=floor(size(dataChunk.data,2)/ceil(dataChunk.Fs)); %excess time points are tacked onto final epoch
        numSamples = min(ceil(dataChunk.Fs),size(dataChunk.data,2)); %one second epochs
    else
        numTrials = length(dataChunk.cellNames);
        numSamples=length(dataChunk.timeNames);
    end;
    numChans=length(dataChunk.chanNames);
    numSubs=length(dataChunk.subNames);
    EEGchans=find(strcmp('EEG',dataChunk.chanTypes));
    ETchans=find(ismember(dataChunk.chanTypes,{'PPL','XEY','YEY'}));
    if strcmp(dataChunk.dataType,'continuous')
        displayPeriod=size(dataChunk.data,2);    %Number of timepoints to graph in display.
    else
        displayPeriod=size(dataChunk.data,2)*size(dataChunk.data,3);
    end;
    decimateSamples=ceil(max(1,displayPeriod/10000));
    totalDisplayPeriod=displayPeriod*size(dataChunk.data,4);
    badData=zeros(numChans,displayPeriod);
    
    badChanNum=dataChunk.analysis.badChans;
    badTrialNum=dataChunk.analysis.badTrials;
    
    if standAlone
        trialdata=reshape(dataChunk.data(:,:,:,theSubject),numChans,[]);
        figure(butterflyFig(iChunk));
        subplot(3,1,1), plot([1:displayPeriod],trialdata(:,1:displayPeriod));
        axis([1 displayPeriod -200 200])
        set(gca,'XTickLabel','','XTick',[]);
        title([deblank(inFile) '-' num2str(iChunk)],'Interpreter','none');
    end;
    
    for trial=1:numTrials
        if strcmp(dataChunk.dataType,'continuous')
            realTrial=1;
        else
            realTrial=trial;
        end;
        if strcmp(dataChunk.dataType,'single_trial')
            fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d of %4d','Working on trial# ', trial+trialCount-1, numTrials))
        end;
        if strcmp(dataChunk.dataType,'continuous')
            fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d of %4d','Working on epoch# ', trial+trialCount-1,numTrials))
        end;
        if badTrialNum(theSubject,trial+trialCount-1)
            if strcmp(dataChunk.dataType,'continuous')
                if trial == numTrials %excess time points are tacked onto final epoch
                    badData(:,(trial-1)*numSamples+1:end)=dataChunk.data(:,(trial-1)*numSamples+1:end,1,theSubject);
                    dataChunk.data(:,(trial-1)*numSamples+1:end,1,theSubject)=0; %set bad trials to flat lines
                    if ~isempty(ETchans)
                        dataChunk.data(ETchans,(trial-1)*numSamples+1:end,1,theSubject)=NaN; %set eye tracker data to NaN
                    end;
                else
                    badData(:,(trial-1)*numSamples+1:trial*numSamples)=dataChunk.data(:,(trial-1)*numSamples+1:trial*numSamples,1,theSubject);
                    dataChunk.data(:,(trial-1)*numSamples+1:trial*numSamples,1,theSubject)=0; %set bad trials to flat lines
                    if ~isempty(ETchans)
                        dataChunk.data(ETchans,(trial-1)*numSamples+1:trial*numSamples,1,theSubject)=NaN; %set eye tracker data to NaN
                    end;
                end;
            else
                badData(:,(trial-1)*numSamples+1:trial*numSamples)=dataChunk.data(:,:,trial,theSubject);
                dataChunk.data(:,:,trial,theSubject)=0; %set bad trials to flat lines
                if ~isempty(ETchans)
                    dataChunk.data(ETchans,:,trial,theSubject)=NaN; %set eye tracker data to NaN
                end;
            end;
        else
            trialBadChans=find(badChanNum(theSubject,trialCount+trial-1,:));
            trialGoodChans=setdiff(EEGchans,trialBadChans);
            if ~isempty(trialBadChans) && ~isempty(trialGoodChans) && (dataChunk.avgNum(theSubject,realTrial)~=-1)  
                if strcmp(dataChunk.dataType,'continuous')
                    if trial == numTrials %excess time points are tacked onto final epoch
                        badData(trialBadChans,(trial-1)*numSamples+1:end)=dataChunk.data(trialBadChans,(trial-1)*numSamples+1:end,1,theSubject);
                        theData=squeeze(dataChunk.data(:,(trial-1)*numSamples+1:end,1,theSubject));
                    else
                        badData(trialBadChans,(trial-1)*numSamples+1:trial*numSamples)=dataChunk.data(trialBadChans,(trial-1)*numSamples+1:trial*numSamples,1,theSubject);
                        theData=squeeze(dataChunk.data(:,(trial-1)*numSamples+1:trial*numSamples,1,theSubject));
                    end;
                else
                    badData(trialBadChans,(trial-1)*numSamples+1:trial*numSamples)=dataChunk.data(trialBadChans,:,trial,theSubject);
                    theData=squeeze(dataChunk.data(:,:,trial,theSubject));
                end;
                for sample=1:numSamples
                    [Xi,Yi,Zi] = griddata(x(trialGoodChans),y(trialGoodChans),theData(trialGoodChans,sample),[1:GRID_SCALE]',[1:GRID_SCALE],'v4');
                    %v4 interpolates to the edge of the box.  With other interpolation options, if a bad channel is at the edge of the montage then it
                    %would just be NaN.
                    for theBadchan=1:length(trialBadChans)
                        badchan=trialBadChans(theBadchan);
                        theData(badchan,sample)=Zi(y(badchan),x(badchan));
                    end;
                end;
                %change badChan field to positive to indicate those channels have been corrected.
                dataChunk.analysis.badChans(theSubject,trialCount+trial-1,trialBadChans)=abs(dataChunk.analysis.badChans(theSubject,trialCount+trial-1,trialBadChans));
                
                %rereference the data.
                if strcmp('AVG',dataChunk.reference.type)
                    refChans=EEGchans;
                elseif strcmp('CSD',dataChunk.reference.type)
                    %do nothing
                else
                    refChans=dataChunk.reference.current;
                    if isempty(refChans)
                        refChans=dataChunk.reference.original;
                    end;
                end;
                if ~isempty(refChans)
                    for i = 1:numTrials
                        epoch=theData;
                        referenceData=mean(epoch(refChans,:),1);
                        for iChan=1:length(EEGchans)
                            theChan=EEGchans(iChan);
                            theData(theChan,:)=epoch(theChan,:)-referenceData;
                        end;
                    end;
                end;
                
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
    end;
    fprintf('%60s\n',' ' );
    
    trialdata=reshape(dataChunk.data(:,:,:,theSubject),numChans,[]);
    
    if ~isempty(butterflyFig)
        if standAlone
            figure(butterflyFig(iChunk));
            subplot(3,1,2), plot([1:decimateSamples:totalDisplayPeriod],badData(EEGchans,:));
            axis([1 displayPeriod -200 200])
            set(gca,'XTickLabel','','XTick',[]);
            title('bad data','Interpreter','none');
            subplot(3,1,3), plot([1:decimateSamples:totalDisplayPeriod],trialdata(EEGchans,:));
            axis([1 displayPeriod -200 200])
            set(gca,'XTickLabel','','XTick',[]);
            title('with bad channels replaced and bad trials zeroed','Interpreter','none');
        else
            figure(butterflyFig(iChunk));
            theTitle='bad data';
            plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,badData,EEGchans,theSubject);
            subplot(numGraphs,1,graphCounter), plot([1:decimateSamples:totalDisplayPeriod],plotData);
            axis([1 totalDisplayPeriod -200 200])
            set(gca,'XTickLabel','','XTick',[]);
            title(theTitle,'Interpreter','none');
            
            theTitle='with bad channels replaced and bad trials zeroed';
            plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,trialdata,EEGchans,theSubject);
            subplot(numGraphs,1,graphCounter+1), plot([1:decimateSamples:totalDisplayPeriod],plotData);
            axis([1 totalDisplayPeriod -200 200])
            set(gca,'XTickLabel','','XTick',[]);
            title(theTitle,'Interpreter','none');
        end;
    end;
    
    drawnow
    
    eval (['save ''' inFile '-' num2str(iChunk) '.mat'' dataChunk']);
    if standAlone
        try
            MATLABver=ver('MATLAB');
            [a b]=strtok(MATLABver.Version,'.');
            b=b(2:end);
            if ~isprop(butterflyFig,'Number')
                eval (['print -f' num2str(butterflyFig(iChunk)) ' -djpeg ''' inFile '''-' num2str(iChunk) 'badData.jpg']);
            else
                eval (['print -f' num2str(butterflyFig(iChunk).Number) ' -djpeg ''' inFile '''-' num2str(iChunk) 'badData.jpg']);
            end;
        catch
            disp('Couldn''t save a copy of the bad channel correction figure.  Perhaps your version of Matlab is not current.');
        end;
    end;
    
    trialCount=trialCount+numTrials;
end;
fprintf('%s\r',' ' );
warning on MATLAB:griddata:DuplicateDataPoints

if standAlone
    close(butterflyFig);
end;
if ~isempty(butterflyFig)
    graphCounter=graphCounter+2;
end;
