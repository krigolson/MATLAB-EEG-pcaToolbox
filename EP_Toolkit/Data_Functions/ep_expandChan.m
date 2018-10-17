function ep_expandChan(src,eventdata,theFunction)
% ep_expandChan - ep_expandChan -
% Expands a channel into a full window for show waves function.
%

%History
%  by Joseph Dien (11/22/13)
%  jdien07@mac.com
%
% bugfix 1/12/14 JD
% Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% modified 2/26/14 JD
% Added View function option to plot or erpimage all trials and all subjects.
%
% modified 4/20/14 JD
% Added coherence and phase-locking options, including support for complex numbers.
%
% modified 6/19/14 JD
% Added support for sample-by-sample t-tests, including STS chanType.
%
% modified 5/25/14 JD
% Set colormap to jet even for Matlab 2014b onwards.
%
% bugfix 6/12/17 JD
% Fixed crash when displaying TFT data and not all four colors are being used.
%
% modified 6/18/17 JD
% For TFT data, when set to display only one Hz bin, switches to waveform display rather than erpimage display.
% Fixed single-sample duration sampleTest results not displaying in View Waves.
% Now displays expanded channels for both Waves and Topos.
%
% modified 2/9/18 JD
% Added support for GFP plots and error bands.
%
% bugfix & modified 2/23/18 JD
% Fixed crash when expanding waveform in Topos.
% Added -all- and -erpimage- options for cells and factors.
%
% bugfix & modified 4/27/18 JD
% Fixed crash when expanding Topos figure of factors.
% Added support for View Topos listing all trials or cells or subjects.
%
% modified 6/5/18 JD
% Improved support for RT event marks.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global EPwaves EPtopos EPmain EPdataset

if isempty(EPmain.view)
    warndlg('When you leave the View pane, the data linked to waveform plots are no longer available.  You''ll need to generate a new waveform plot.');
    return
end;

scrsz = EPmain.scrsz;

switch theFunction
    case 'EPwaves'
        functionData=EPwaves;
    case 'EPtopos'
        functionData=EPtopos;
end;

if isfield(functionData.handles.waves,'hWave')
    chan=find([functionData.handles.waves.hWave] == src);
else
    chan=[];
end;

if isempty(chan)
    for iChan=1:length(functionData.handles.waves.hLines)
        for theLine=1:length(functionData.handles.waves.hLines{iChan})
            if src == functionData.handles.waves.hLines{iChan}(theLine)
                chan=iChan;
            end;
        end;
    end;
end;

if isempty(chan)
    disp('Could not find the channel for some reason.')
    return
end;

if strcmp('EPtopos',theFunction)
    row=chan;
    chan=EPtopos.chans(row);
    if EPtopos.complexData
        rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
    else
        rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
    end;
    if EPtopos.numRows > 1
        theFac=rowCounter;
    else
        theFac=1;
    end;
    if EPtopos.complexData && rem(row,2)
        theData2=EPtopos.totalImagData(:,:,theFac,:,:,EPtopos.plotColors);
    else
        theData2=EPtopos.totalData(:,:,theFac,:,:,EPtopos.plotColors);
    end
    if any(EPmain.view.rel(EPmain.view.dataset <= length(EPdataset.dataset)))
        theData2(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors)))=abs(theData2(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors))));
        theData2=mean(theData,5); %collapse over relations if any
    end;
    for iCol=1:length(EPtopos.plotColors)
        for iFreq=1:size(theData2,4)
            theData(:,:,iCol,iFreq)=theData2(:,:,:,iFreq,:,iCol);
        end;
    end;
    theChan=EPtopos.chans(rowCounter);
else
    theData=EPwaves.totalData;
    theChan=chan;
    rowCounter=1;
end;

% if isfield(functionData.handles.waves,'hExpandedFigure')
%     if ~isempty(functionData.handles.waves.hExpandedFigure)
%         if ishandle(functionData.handles.waves.hExpandedFigure)
%             close(functionData.handles.waves.hExpandedFigure)
%             functionData.handles.waves.hExpandedFigure=[];
%         end;
%     end;
% end;

expandChan.handles.waves.hExpandedFigure = figure('Name', [functionData.chanNames{theChan}], 'NumberTitle', 'off', 'Position',[scrsz(3)/2 scrsz(4)/2 600 400]);
colormap jet;

plotForm=EPmain.view.dataTransform;
if strcmp('TFT',EPmain.view.dataTransform) && (functionData.firstHz==functionData.lastHz)
    plotForm='VLT';
end;
if strcmp('TFT',EPmain.view.dataTransform) && (functionData.firstTime==functionData.lastTime)
    plotForm='FFT';
end;

if strcmp('FFT',plotForm)
    sampleSize=0;
    spacing=(functionData.lastHz-functionData.firstHz)/(functionData.numHz-1);
else
    sampleSize=functionData.sampleSize;
    spacing=functionData.sampleSize;
end;

waveList=find(~ismember(functionData.colorIndex,find(ismember(EPmain.view.allTrials,[2 4 6]))));
theMarker='none';
theMarkerSize=2;
plotDataList=find(EPmain.view.dataset <= length(EPdataset.dataset));

switch plotForm
    case 'VLT'
        numImages=length(find(ismember(EPmain.view.allTrials,[2 4 6])));
        if (length(functionData.plotColors)-numImages) > 0
            imageSpace=4;
        else
            imageSpace=numImages;
        end;
        expandChan.handles.waves.hExpandedAxes=[];
        if numImages %if any erpimage
            imageCount=0;
            for iColor=1:4
                if (EPmain.view.dataset(iColor) <= length(EPdataset.dataset)) && ismember(EPmain.view.allTrials(iColor),[2 4 6])
                    imageCount=imageCount+1;
                    trialList=find(functionData.colorIndex==iColor);
                    expandChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .05+(.9/imageSpace)*(imageSpace-imageCount) .9 (.9/imageSpace)]);
                    expandChan.handles.waves.hExpandedAxes(end+1,1) = imagesc(functionData.firstTime:functionData.lastTime,1:length(trialList),squeeze(theData(theChan,:,trialList,:))',[functionData.plotMVmin, functionData.plotMVmax]);
                    axis([functionData.firstTime functionData.lastTime 1 length(trialList)]);
                    line([0 0],[1 length(trialList)],'Color','black','LineWidth',1) %stimulus onset
                    if ~isempty(functionData.marker1)
                        line(repmat(functionData.marker1,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                    if ~isempty(functionData.marker2)
                        line(repmat(functionData.marker2,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                    %plot event lines
                    eventX=[];
                    eventY=[];
                    for iLine=1:length(functionData.eventLines{iColor})
                        if ~isempty(functionData.eventLines{iColor}{iLine})
                            eventX(end+1)=functionData.eventLines{iColor}{iLine};
                            eventY(end+1)=iLine;
                        end;
                    end;
                    if ~isempty(eventX)
                        theTimePoints=[functionData.firstTime:functionData.spacing:functionData.firstTime+(functionData.spacing*(functionData.numPoints-1))];
                        hold on
                        plot(theTimePoints(eventX),eventY)
                        hold off
                    end
                end;
            end;
        end
        if (length(functionData.plotColors)-numImages) > 0 %if there will be waveforms
            expandChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .05 .9 .9*((4-numImages)/4)]);
            hold on
            for iWave=1:length(waveList)
                theWave=waveList(iWave);
                theColor=functionData.colorIndex(theWave);
                if functionData.complexData
                    if strcmp(functionData.plotLineIndex{theWave},':')
                        theMarker='none';
                        theMarkerSize=2;
                    else
                        theMarker='none';
                        theMarkerSize=2;
                    end;
                end;
                if functionData.STSdata==functionData.colorIndex(iWave)
                    breakList=sort([find(diff([0 (squeeze(theData(theChan,:,theWave,:))>0) 0])<0)-1 find(diff([0 (squeeze(theData(theChan,:,theWave,:))>0) 0])>0)]);
                    if ~isempty(breakList)
                        theSTdata=squeeze(theData(theChan,:,theWave,:));
                        theData1=squeeze(theData(theChan,:,min(setdiff(plotDataList,functionData.STSdata)),:));
                        theData2=squeeze(theData(theChan,:,max(setdiff(plotDataList,functionData.STSdata)),:));
                        for iSigArea=1:length(breakList)/2
                            theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)];
                                if length(theTimePoints) == 1
                                    expandChan.handles.waves.hLines{chan}(iWave)=line(([theTimePoints theTimePoints]*spacing)+(functionData.firstTime-spacing),[theData1(theTimePoints) theData2(theTimePoints)],'LineWidth',1,'Color',[1 .5 .5]);
                                else
                                    expandChan.handles.waves.hLines{chan}(iWave)=patch(([theTimePoints flip(theTimePoints)]*spacing)+(functionData.firstTime-spacing),[theData1(theTimePoints) theData2(flip(theTimePoints))],functionData.plotColorIndex(iWave,:),'FaceColor',functionData.plotColorIndex(iWave,:),'EdgeColor','none','FaceAlpha',.25);
                                end;
                        end;
                        %                             [A maxPoint]=max(theSTdata);
                        %                             line([maxPoint maxPoint]*spacing-(functionData.firstTime+spacing),[theData1(maxPoint) theData2(maxPoint)],'Color',functionData.plotColorIndex(iWave,:),'LineWidth',1); %maximum sample test point
                        %                             patch(([theTimePoints flip(theTimePoints)]*spacing)-(functionData.firstTime+spacing),[theData1(theTimePoints) theData2(flip(theTimePoints))],functionData.plotColorIndex(iWave,:),'FaceColor',functionData.plotColorIndex(iWave,:),'EdgeColor','none','FaceAlpha',.5);
                        
                    end;
                else
                    if strcmp('EPtopos',theFunction)
                        theTimePoints=[EPtopos.firstTime:spacing:EPtopos.firstTime+(spacing*(EPtopos.numPoints-1))];
                        EPtopos.handles.waves.hLines{chan}(iWave)=plot(theTimePoints,squeeze(EPtopos.totalData(chan,:,theFac,:,:,theWave)),'LineStyle',EPtopos.plotLineIndex{theWave},'color',EPtopos.plotColorIndex(theWave,:),'Marker',theMarker,'MarkerSize',theMarkerSize);
                    else
                        theTimePoints=[EPwaves.firstTime:spacing:EPwaves.firstTime+(spacing*(EPwaves.numPoints-1))];
                        EPwaves.handles.waves.hLines{chan}(iWave)=plot(theTimePoints,squeeze(EPwaves.totalData(chan,:,theWave,:)),'LineStyle',EPwaves.plotLineIndex{theWave},'color',EPwaves.plotColorIndex(theWave,:),'Marker',theMarker,'MarkerSize',theMarkerSize);
                        if EPwaves.bandIndex(theColor) > 0
                            theBand1=squeeze(EPwaves.totalData(chan,:,theWave,:))+EPwaves.bandData(chan,:,1,:,theColor);
                            theBand2=squeeze(EPwaves.totalData(chan,:,theWave,:))-EPwaves.bandData(chan,:,1,:,theColor);
                            EPwaves.handles.waves.hLines{chan}(iWave+1)=patch([theTimePoints flip(theTimePoints)],[theBand1 flip(theBand2)],EPwaves.plotColorIndex(iWave,:),'FaceColor',EPwaves.plotColorIndex(iWave,:),'EdgeColor','none','FaceAlpha',.25);
                        end;
                    end;
                end;
            end;
            hold off
            axis([functionData.firstTime functionData.lastTime functionData.plotMVmin functionData.plotMVmax]);
            if functionData.direction ==2
                set(expandChan.handles.waves.hWave(chan),'YDir','reverse')
            end;
            line([functionData.firstTime functionData.lastTime-functionData.sampleSize],[0 0],'Color','black','LineWidth',1) % zero line
            line([0 0],[0 functionData.plotMVmax],'Color','black','LineWidth',1) %stimulus onset
            
            %plot event lines
             if (functionData.plotMVmin < 0) && (functionData.plotMVmax >= 0)
                for color=1:length(functionData.plotColors)
                    theColor=functionData.plotColors(color);
                    if strcmp('EPtopos',theFunction) && ~any(strcmp(functionData.type,{'time','freq'}))
                        theRow=rowCounter;
                    else
                        theRow=1;
                    end;
                    if ~isempty(functionData.eventWave{theColor}{theRow})
                        plotPoints=find(functionData.eventWave{theColor}{theRow}>min(functionData.eventWave{theColor}{theRow}));
                        plotTimes=[functionData.firstTime:spacing:functionData.firstTime+(spacing*(functionData.numPoints-1))];
                        if length(plotPoints)==1
                            line([plotTimes(plotPoints) plotTimes(plotPoints)],[functionData.plotMVmin functionData.eventWave{theColor}{theRow}(plotPoints)*(functionData.plotMVmin/2)],'Color',functionData.thePlotColors(theColor,:),'LineWidth',5) %event line
                        else
                            hold on
                            expandChan.handles.waves.eventLines{chan,theColor} = plot(plotTimes,(functionData.eventWave{theColor}{theRow}*(abs(functionData.plotMVmin/2)))+functionData.plotMVmin,'LineWidth',5,'Color',functionData.thePlotColors(theColor,:));
                            hold off
                        end;
                    end;
                end;
            end;
        end;
        
        text(.05,.8, functionData.chanNames(chan), 'Units','normalized');
        if ~isempty(functionData.marker1)
            try
                eval('line(repmat(functionData.marker1,2),[functionData.plotMVmin functionData.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1)');
            catch
            end
        end;
        if ~isempty(functionData.marker2)
            try
                eval('line(repmat(functionData.marker2,2),[functionData.plotMVmin functionData.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1);');
            catch
            end;
        end;
    case 'FFT'
        numImages=length(find(ismember(EPmain.view.allTrials,[2 4 6])));
        if (length(functionData.plotColors)-numImages) > 0
            imageSpace=4;
        else
            imageSpace=numImages;
        end;
        expandChan.handles.waves.hExpandedAxes=[];
        if numImages %if any erpimage
            imageCount=0;
            for iColor=1:4
                if (EPmain.view.dataset(iColor) <= length(EPdataset.dataset)) && ismember(EPmain.view.allTrials(iColor),[2 4 6])
                    imageCount=imageCount+1;
                    trialList=find(functionData.colorIndex==iColor);
                    expandChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .05+(.9/imageSpace)*(imageSpace-imageCount) .9 (.9/imageSpace)]);
                    expandChan.handles.waves.hExpandedAxes(end+1,1) = imagesc(functionData.firstHz+sampleSize:functionData.lastHz,1:length(trialList),squeeze(theData(theChan,:,trialList,:))',[functionData.plotMVmin, functionData.plotMVmax]);
                    axis([functionData.firstHz+sampleSize functionData.lastHz 1 length(trialList)]);
                    if ~isempty(functionData.marker1)
                        line(repmat(functionData.marker1,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                    if ~isempty(functionData.marker2)
                        line(repmat(functionData.marker2,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                end;
            end;
        end
        if (length(functionData.plotColors)-numImages) > 0 %if there will be waveforms
            expandChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .05 .9 .9*((4-numImages)/4)]);
            hold on
            for iWave=1:length(waveList)
                theWave=waveList(iWave);
                plot([functionData.firstHz+sampleSize:spacing:functionData.lastHz],squeeze(theData(theChan,:,theWave,:)),'LineStyle',functionData.plotLineIndex{theWave},'color',functionData.plotColorIndex(theWave,:));
            end;
            hold off
            axis([functionData.firstHz functionData.lastHz functionData.plotMVmin functionData.plotMVmax]);
            line([functionData.firstHz+sampleSize functionData.lastHz],[0 0],'Color','black','LineWidth',1) % zero line
        end;
        text(.1,.1, functionData.chanNames(theChan), 'Units','normalized');
        if ~isempty(functionData.marker1)
            try
                eval('line(repmat(functionData.marker1,2),[functionData.plotMVmin functionData.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1)');
            catch
            end
        end;
        if ~isempty(functionData.marker2)
            try
                eval('line(repmat(functionData.marker2,2),[functionData.plotMVmin functionData.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1);');
            catch
            end;
        end;
    case 'TFT'
        imageSpace=length(functionData.plotColors);
        imageCount=0;
        for iColor=1:4
            if (EPmain.view.dataset(iColor) <= length(EPdataset.dataset))
                imageCount=imageCount+1;
                expandChan.handles.waves.hExpandedAxes(iColor) = axes('position',[.05 .05+(.9/imageSpace)*(imageSpace-imageCount) .9 (.9/imageSpace)]);
                expandChan.handles.waves.hExpandedAxes(4+iColor) = imagesc(functionData.firstTime:functionData.lastTime,functionData.firstHz:functionData.lastHz,squeeze(theData(theChan,:,(functionData.plotColors==iColor),:))');
                axis([functionData.firstTime functionData.lastTime functionData.firstHz functionData.lastHz]);
                line([0 0],[functionData.firstHz functionData.lastHz],'Color','black','LineWidth',1) %stimulus onset
                if ~isempty(functionData.marker1)
                    line(repmat(functionData.marker1,2),[functionData.firstHz functionData.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
                end;
                if ~isempty(functionData.marker2)
                    line(repmat(functionData.marker2,2),[functionData.firstHz functionData.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
                end;
            end;
        end;
        %text(.1,.1, functionData.chanNames(theChan), 'Units','normalized');
end;



