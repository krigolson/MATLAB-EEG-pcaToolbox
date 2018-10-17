function ep_scanCont
% ep_scanCont - ep_scanCont -
% Presents continuous data for manual editing.
%
%Input:
%

%History
%  by Joseph Dien (1/24/16)
%  jdien07@mac.com
%
% modified 5/4/16 JD
% 2D plot now excludes bad channels.
%
% modified 9/6/16 JD
% Added temporally arranged event list for navigating to individual events.
% Event lists for setting boundaries of Redisplay are now alphabetic order.
% Global bad channels marked with faded line rather than red zone.
%
% bugfix 9/21/16 JD
% Fixed bad channels being displayed as still bad.
% Fixed crash when clicking on "- Channels"
%
% modified 9/22/16 JD
% Scaling of non-EEG channels adjusted so variations are visible.
%
% modified 11/4/16 JD
% Upgraded controls for navigating and examining events.
% Upgraded eye-position figures.
%
% bugfix & modified 11/26/16 JD
% Added support for .calibration field
% Fixed crash when saving edits in View>Scan for continuous data.
% Reversed + and - buttons to comply with norms for such buttons.
%
% bugfix 2/3/17 JD
% Accommodated the datasets having differing channels.
%
% modified 4/21/17 JD
% Added toggle for edit mode.
%
% bugfix 6/16/17 JD
% Fixed conversion to spectral density dividing by bin width rather than sqrt(bin width).
%
% bugfix 5/4/18 JD
% Fixed crash when data has channels without electrode coordinates.
% Fixed crash in View Scan when clicking for plot point without first advancing the plotting window.
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

global EPdataset EPscanCont EPmain

scrsz = EPmain.scrsz;

EPscanCont.edited=0;

EPmain.handles.scanContData = figure('Name', 'Scan Continuous Data', 'NumberTitle', 'off', 'Position',[270 1 scrsz(3)-270 scrsz(4)], 'MenuBar', 'none');
colormap jet;
drawnow

if ~isempty(EPdataset.dataset(EPmain.view.dataset(1)).freqNames)
    EPscanCont.numColors=1;
    EPscanCont.dataType='TFT';
    EPscanCont.startBins=min(find(EPmain.view.startHz <= EPdataset.dataset(EPmain.view.dataset(1)).freqNames));
    EPscanCont.lastBins=min(find(EPmain.view.endHz <= EPdataset.dataset(EPmain.view.dataset(1)).freqNames));
	EPscanCont.spacing=EPscanCont.lastBins-EPscanCont.startBins+1;
    EPscanCont.FFTunits=EPmain.view.FFTunits;
    if EPscanCont.FFTunits==1
        EPscanCont.spacing=EPscanCont.spacing*2; %complex numbers
    end;
    EPscanCont.freqList=cellstr(num2str(round(EPdataset.dataset(EPmain.view.dataset(1)).freqNames)));
    EPscanCont.freqList{end+1}='all';
    EPscanCont.freqShow=length(EPscanCont.freqList);
    EPscanCont.freqPos=1;
else
    EPscanCont.numColors=4;
    EPscanCont.dataType='VLT';
    EPscanCont.spacing=100; %+/-microvolt space for each waveform
    EPscanCont.freqList={'none'};
    EPscanCont.freqShow=length(EPscanCont.freqList);
    EPscanCont.freqPos=1;
end;

EPscanCont.chanNames=cell(0);
EPscanCont.chanTypes=cell(0);
EPscanCont.eloc=[];
EPscanCont.chanIndex=cell(EPscanCont.numColors,1);
for iColor=1:EPscanCont.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        EPdata=ep_loadEPdataset(EPdataset,EPmain.view.dataset(iColor));
        EPscanCont.pca{iColor}=EPdata.pca;
        EPdata=rmfield(EPdata,'pca'); %pca field contents too variable to merge into one structure so temporarily remove
        EPscanCont.EPdata(iColor)=EPdata;
        for iChan=1:length(EPdata.chanNames)
            theChan=find(strcmp(EPdata.chanNames{iChan},EPscanCont.chanNames));
            if isempty(theChan)
                EPscanCont.chanNames{end+1}=EPdata.chanNames{iChan};
                EPscanCont.chanTypes{end+1}=EPdata.chanTypes{iChan};
                EPscanCont.chanIndex{iColor}(iChan)=length(EPscanCont.chanNames);
                if ~isempty(EPdata.eloc)
                    if isempty(EPscanCont.eloc)
                        EPscanCont.eloc=EPdata.eloc(iChan);
                    else
                        EPscanCont.eloc(iChan)=EPdata.eloc(iChan);
                    end;
                end;
            else
                EPscanCont.chanIndex{iColor}(iChan)=theChan;
            end;
        end;
    end;
end;

EPscanCont.numChans=length(EPscanCont.chanNames);
EPscanCont.numPoints=length(EPscanCont.EPdata(1).timeNames);
EPscanCont.sampleSize=1000/EPscanCont.EPdata(1).Fs;
EPscanCont.xGraphMin=min(EPscanCont.numPoints,ceil(EPscanCont.EPdata(1).Fs)); %minimum number of time points to be displayable
EPscanCont.displayLength=min(EPscanCont.numPoints,ceil(EPscanCont.EPdata(1).Fs)*4); %start at displaying 4 seconds worth of data
EPscanCont.plotPoint=1; %start at displaying 4 seconds worth of data
EPscanCont.scaling=zeros(EPscanCont.numChans,1); %scaling for waveforms.  For EEG, scaling is one.
EPscanCont.displayChans=[1:EPscanCont.numChans]; %channels to be displayed
EPscanCont.displayPoints=[1:EPscanCont.displayLength]; %points to be displayed
EPscanCont.xGraphPos=0; %current x position of the graph from 0 to 1
EPscanCont.yGraphPos=0; %current y position of the graph from 0 to 1
EPscanCont.currentData=[]; %data to be displayed
EPscanCont.freqData=[]; %freq data to be displayed
EPscanCont.plotLineEpoch=0; %epoch of plot line (excess points counted as an epoch)
EPscanCont.freqBins=[1:40];
EPscanCont.freqLines=[4 8 12 20];
EPscanCont.eventList=cell(length(EPscanCont.EPdata(1).events{1}),1); %list of events
EPscanCont.firstEvent=1; %choice of event for start of new display (from alphabetically-ordered event list)
EPscanCont.lastEvent=1; %choice of event for end of new display (from alphabetically-ordered event list)
EPscanCont.currentEvent=[]; %choice of event to center on (from chronologically-ordered event list)
EPscanCont.starredSampleEventList=cell(0); %list of events in chronological order
EPscanCont.starredEventList=cell(0); %list of events in alphabetic order
EPscanCont.globalBadChans=zeros(EPscanCont.numChans,1);
EPscanCont.lineStyles=cell(EPscanCont.numChans,1);
EPscanCont.currentkey=1;

EPscanCont.chanList=EPscanCont.chanNames;
EPscanCont.chanList{end+1}='none';
EPscanCont.chanShow=length(EPscanCont.chanList);

%prepare to display stimulus screens in eye-plot figures and list of event names
sampleList=zeros(length(EPscanCont.EPdata(1).events{1}),1);
EPscanCont.stimList=[];
for iEvent=1:length(EPscanCont.EPdata(1).events{1})
    theEventValue=EPscanCont.EPdata(1).events{1}(iEvent).value;
    theEventSample=EPscanCont.EPdata(1).events{1}(iEvent).sample;
    sampleList(iEvent,1)=theEventSample;
    keyData=EPscanCont.EPdata(1).events{1}(iEvent).keys;
    keyStim=[];
    for iKey=1:length(keyData)
        if strcmp(keyData(iKey).code,'stim')
            keyStim=keyData(iKey).data;
            if any(strcmp(keyStim,{EPscanCont.EPdata(1).stims.name}))
                EPscanCont.stimList(end+1).sample=EPscanCont.EPdata(1).events{1}(iEvent).sample;
                EPscanCont.stimList(end).stim=find(strcmp(keyStim,{EPscanCont.EPdata(1).stims.name}));
            end;
        end;
    end;
    if ~isempty(theEventValue)
        eventNum=sprintf('%04d',(length(find([EPscanCont.EPdata(1).events{1}(strcmp(theEventValue,{EPscanCont.EPdata(1).events{1}.value})).sample]<=theEventSample))));
        EPscanCont.eventList{iEvent}=[theEventValue '(' eventNum ')' keyStim];
    else
        theEventType=EPscanCont.EPdata(1).events{1}(iEvent).type;
        if ~isempty(theEventType)
            eventNum=sprintf('%04d',(length(find([EPscanCont.EPdata(1).events{1}(strcmp(theEventType,{EPscanCont.EPdata(1).events{1}.value})).sample]<=theEventSample))));
        	EPscanCont.eventList{iEvent}=[theEventType '(' eventNum ')' keyStim];
        else
            EPscanCont.eventList{iEvent}='null_event';
        end;
    end;
end;

[B EPscanCont.alphabetEventOrder]=sort(EPscanCont.eventList);
[B EPscanCont.sampleEventOrder]=sort(sampleList);

if isempty(EPscanCont.eventList)
    EPscanCont.starredSampleEventList='none';
    EPscanCont.starredEventList='none';
else
    EPscanCont.starredSampleEventList=EPscanCont.eventList;
    EPscanCont.starredEventList=EPscanCont.eventList;
end;

EPscanCont.epochMarks=zeros(EPscanCont.numPoints,1); %provides the epoch boundary marks
EPscanCont.epochMarks(EPscanCont.EPdata(1).Fs:EPscanCont.EPdata(1).Fs:end,1)=1; %add boundary marks

EPscanCont.XEYchan=find(strcmp('XEY',EPscanCont.chanTypes));
if length(EPscanCont.XEYchan) > 1
    EPscanCont.XEYchan=EPscanCont.XEYchan(1);
end;
EPscanCont.YEYchan=find(strcmp('YEY',EPscanCont.chanTypes));
if length(EPscanCont.YEYchan) > 1
    EPscanCont.YEYchan=EPscanCont.YEYchan(1);
end;
if ~isempty(EPscanCont.XEYchan) && ~isempty(EPscanCont.YEYchan)
    EPscanCont.EYEsize=max([max(abs(EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.XEYchan),:))) max(abs(EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.YEYchan),:)))]);
    if isfield(EPscanCont.EPdata(1).calibration,'ET') && isfield(EPscanCont.EPdata(1).calibration.ET,'Xzero') && isfield(EPscanCont.EPdata(1).calibration.ET,'Yzero') && isfield(EPscanCont.EPdata(1).calibration.ET,'Xscale') && isfield(EPscanCont.EPdata(1).calibration.ET,'Yscale')
        EPscanCont.eyeCenterX=EPscanCont.EPdata(1).calibration.ET.Xzero;
        EPscanCont.eyeCenterY=EPscanCont.EPdata(1).calibration.ET.Yzero;
        EPscanCont.eyeScaleX=EPscanCont.EPdata(1).calibration.ET.Xscale;
        EPscanCont.eyeScaleY=EPscanCont.EPdata(1).calibration.ET.Yscale;
    else
        if ~isempty(EPscanCont.stimList) %SMI POR field is already in pixels of the stimuli jpegs.
            EPscanCont.eyeScaleX=max(size(EPscanCont.EPdata(1).stims(1).image));
            EPscanCont.eyeScaleY=-max(size(EPscanCont.EPdata(1).stims(1).image)); %need to flip SMI coordinates
            EPscanCont.eyeCenterX=round(size(EPscanCont.EPdata(1).stims(1).image,2)/2);
            EPscanCont.eyeCenterY=round(size(EPscanCont.EPdata(1).stims(1).image,1)/2);
        else
            nonNaNx=~isnan(EPscanCont.EPdata(1).data(EPscanCont.XEYchan,:));
            nonNaNy=~isnan(EPscanCont.EPdata(1).data(EPscanCont.YEYchan,:));
            EPscanCont.eyeCenterX=median(EPscanCont.EPdata(1).data(EPscanCont.XEYchan,nonNaNx));
            EPscanCont.eyeCenterY=median(EPscanCont.EPdata(1).data(EPscanCont.YEYchan,nonNaNy));
            EPscanCont.eyeScaleX=max(EPscanCont.EPdata(1).data(EPscanCont.XEYchan,nonNaNx)-EPscanCont.eyeCenterX)-min(EPscanCont.EPdata(1).data(EPscanCont.XEYchan,nonNaNx)-EPscanCont.eyeCenterX);
            EPscanCont.eyeScaleY=max(EPscanCont.EPdata(1).data(EPscanCont.YEYchan,nonNaNx)-EPscanCont.eyeCenterY)-min(EPscanCont.EPdata(1).data(EPscanCont.YEYchan,nonNaNx)-EPscanCont.eyeCenterY);
        end;
    end;
else
    EPscanCont.EYEsize=[];
end;

EPscanCont.HSACchan=find(strcmp('Hsaccade',EPscanCont.chanNames));
if length(EPscanCont.HSACchan) > 1
    EPscanCont.HSACchan=EPscanCont.HSACchan(1);
end;
EPscanCont.VSACchan=find(strcmp('Vsaccade',EPscanCont.chanNames));
if length(EPscanCont.VSACchan) > 1
    EPscanCont.VSACchan=EPscanCont.VSACchan(1);
end;
if ~isempty(EPscanCont.HSACchan) && ~isempty(EPscanCont.VSACchan)
    HSACchan=find(EPscanCont.chanIndex{1}==EPscanCont.HSACchan);
    VSACchan=find(EPscanCont.chanIndex{1}==EPscanCont.VSACchan);
    EPscanCont.SACCsize=max([max(abs(EPscanCont.EPdata(1).data(HSACchan,:))) max(abs(EPscanCont.EPdata(1).data(VSACchan,:)))]);
    if EPscanCont.SACCsize ==0
        EPscanCont.SACCsize=[];
    else
        if isfield(EPscanCont.EPdata(1).calibration,'SAC') && isfield(EPscanCont.EPdata(1).calibration.SAC,'Xzero') && isfield(EPscanCont.EPdata(1).calibration.SAC,'Yzero') && isfield(EPscanCont.EPdata(1).calibration.SAC,'Xscale') && isfield(EPscanCont.EPdata(1).calibration.SAC,'Yscale')
            EPscanCont.sacCenterX=EPscanCont.EPdata(1).calibration.SAC.Xzero;
            EPscanCont.sacCenterY=EPscanCont.EPdata(1).calibration.SAC.Yzero;
            EPscanCont.sacScaleX=EPscanCont.EPdata(1).calibration.SAC.Xscale;
            EPscanCont.sacScaleY=EPscanCont.EPdata(1).calibration.SAC.Yscale;
        else
            nonNaNx=~isnan(EPscanCont.EPdata(1).data(HSACchan,:));
            nonNaNy=~isnan(EPscanCont.EPdata(1).data(VSACchan,:));
            EPscanCont.sacCenterX=median(EPscanCont.EPdata(1).data(HSACchan,nonNaNx));
            EPscanCont.sacCenterY=median(EPscanCont.EPdata(1).data(VSACchan,nonNaNy));
            EPscanCont.sacScaleX=max(EPscanCont.EPdata(1).data(HSACchan,nonNaNx)-EPscanCont.sacCenterX)-min(EPscanCont.EPdata(1).data(HSACchan,nonNaNx)-EPscanCont.sacCenterX);
            EPscanCont.sacScaleY=max(EPscanCont.EPdata(1).data(VSACchan,nonNaNx)-EPscanCont.sacCenterY)-min(EPscanCont.EPdata(1).data(VSACchan,nonNaNx)-EPscanCont.sacCenterY);
        end;
    end;
else
    EPscanCont.SACCsize=[];
end;

for iChan=1:EPscanCont.numChans
    if any(strcmp(EPscanCont.chanTypes(iChan),{'EEG','REG'}))
        if strcmp(EPscanCont.dataType,'VLT') 
            EPscanCont.scaling(iChan)=1;
        elseif any(strcmp(EPscanCont.chanNames(iChan),{'Hsaccade','Vsaccade','blink','SacPot'}))
            EPscanCont.scaling(iChan)=1;
        else
            EPscanCont.scaling(iChan)=100;
        end;
    else        
        theMax=0;
        for iColor=1:4
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                theChan=find(EPscanCont.chanIndex{iColor}==iChan);
                if ~isempty(theChan)
                    theMax=max(theMax,max(abs(EPscanCont.EPdata(iColor).data(theChan,:))));
                end;
            end;
        end;
        EPscanCont.scaling(iChan)=EPscanCont.spacing/theMax;
    end;
end;

EPscanCont.EEGchans=find(strcmp('EEG',EPscanCont.chanTypes));
EPscanCont.hasLoc=[];
if ~isempty(EPscanCont.eloc)
    for iChan=1:length(EPscanCont.EEGchans)
        if ~isempty(EPscanCont.eloc(EPscanCont.EEGchans(iChan)).theta)
            EPscanCont.hasLoc(end+1)=EPscanCont.EEGchans(iChan);
        end;
    end;
end;

EPscanCont.waveColors={'blue','red','green','black'};
EPscanCont.graphCoords=[200 220 EPmain.scrsz(3)-700 EPmain.scrsz(4)-320];
EPscanCont.topoSize=180;

EPscanCont.playMode=false;
EPscanCont.eventLabels=true;
EPscanCont.editMode=true;

if ~isempty(EPscanCont.eloc)
    EPscanCont.gridSize=67;
    maxRad=0.5;
    [y,x] = pol2cart(([EPscanCont.eloc(EPscanCont.hasLoc).theta]/360)*2*pi,[EPscanCont.eloc(EPscanCont.hasLoc).radius]);  % transform electrode locations from polar to cartesian coordinates
    y=-y; %flip y-coordinate so that nose is upwards.
    plotrad = min(1.0,max([EPscanCont.eloc(EPscanCont.hasLoc).radius])*1.02);            % default: just outside the outermost electrode location
    plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
    x = x*(maxRad/plotrad);
    y = y*(maxRad/plotrad);
    
    xmin = min(-maxRad,min(x));
    xmax = max(maxRad,max(x));
    ymin = min(-maxRad,min(y));
    ymax = max(maxRad,max(y));
    
    EPscanCont.x=round(((x/(xmax-xmin))*EPscanCont.gridSize)+ceil(EPscanCont.gridSize/2));
    EPscanCont.y=round(((y/(ymax-ymin))*EPscanCont.gridSize)+ceil(EPscanCont.gridSize/2));
end;


for iChan=1:EPscanCont.numChans
    theChan=find(EPscanCont.chanIndex{1}==EPscanCont.displayChans(iChan));
    if all(EPscanCont.EPdata(1).analysis.badChans(1,:,theChan))
        EPscanCont.globalBadChans(theChan)=1;
        EPscanCont.lineStyles{iChan}=':';
    else
        EPscanCont.lineStyles{iChan}='-';
    end;
    EPscanCont.trialBadChans(theChan,:)=EPscanCont.EPdata(1).analysis.badChans(1,:,theChan);
end;
EPscanCont.badTrials=EPscanCont.EPdata(1).analysis.badTrials;

EPscanCont.handles.graphPlot = axes('Units','pixels','position',EPscanCont.graphCoords,'Tag','graphPlot');
if EPmain.preferences.view.positive ==2
    set(EPscanCont.handles.graphPlot,'YDir','reverse')
end;

axis([1 EPscanCont.displayLength -length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing])
set(EPscanCont.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
EPscanCont.handles.plotLine=line([EPscanCont.plotPoint EPscanCont.plotPoint],[-length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing],'Color','green','LineWidth',1);

EPscanCont.handles.freqPlot = axes('Units','pixels','position',[EPscanCont.graphCoords(1)+EPscanCont.graphCoords(3)+50 EPscanCont.graphCoords(2) 170 EPscanCont.graphCoords(4)],'Tag','freqPlot');
axis([1 length(EPscanCont.freqBins) -length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing])
set(EPscanCont.handles.freqPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

if ~isempty(EPscanCont.eloc)
    for iColor=1:EPscanCont.numColors
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            EPscanCont.handles.topos.topo(iColor) = axes('Units','pixel','position',[600+(iColor*EPscanCont.topoSize) 10 EPscanCont.topoSize EPscanCont.topoSize]);
            set(EPscanCont.handles.topos.topo(iColor),'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
        end;
    end;
end;

if ~isempty(EPscanCont.EYEsize)
    EPscanCont.handles.topos.eye = axes('Units','pixel','position',[600+(5*EPscanCont.topoSize)+50 10 EPscanCont.topoSize EPscanCont.topoSize]);
    set(EPscanCont.handles.topos.eye,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
    axis([-EPscanCont.EYEsize EPscanCont.EYEsize -EPscanCont.EYEsize EPscanCont.EYEsize])
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Eye Track','FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-20 50 20]);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Center X','FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-40 50 20]);
    EPscanCont.handles.topos.eyeCenterX = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.eyeCenterX),'FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-60 50 20],'Callback',@changeEyeSettings);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Center Y','FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-80 50 20]);
    EPscanCont.handles.topos.eyeCenterY = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.eyeCenterY),'FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-100 50 20],'Callback',@changeEyeSettings);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Scale X','FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-120 50 20]);
    EPscanCont.handles.topos.eyeScaleX = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.eyeScaleX),'FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-140 50 20],'Callback',@changeEyeSettings);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Scale Y','FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-160 50 20]);
    EPscanCont.handles.topos.eyeScaleY = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.eyeScaleY),'FontSize',EPmain.fontsize,...
        'Position',[600+(5*EPscanCont.topoSize) EPscanCont.topoSize+10-180 50 20],'Callback',@changeEyeSettings);
end;
if ~isempty(EPscanCont.SACCsize)
    EPscanCont.handles.topos.saccade = axes('Units','pixel','position',[600+(6*EPscanCont.topoSize)+100 10 EPscanCont.topoSize EPscanCont.topoSize]);
    set(EPscanCont.handles.topos.saccade,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
    axis([-EPscanCont.SACCsize EPscanCont.SACCsize -EPscanCont.SACCsize EPscanCont.SACCsize])
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Sacc Track','FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-20 50 20]);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Center X','FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-40 50 20]);
    EPscanCont.handles.topos.sacCenterX = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.sacCenterX),'FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-60 50 20],'Callback',@changeEyeSettings);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Center Y','FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-80 50 20]);
    EPscanCont.handles.topos.sacCenterY = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.sacCenterY),'FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-100 50 20],'Callback',@changeEyeSettings);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Scale X','FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-120 50 20]);
    EPscanCont.handles.topos.sacScaleX = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.sacScaleX),'FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-140 50 20],'Callback',@changeEyeSettings);
    uicontrol('Style','text','HorizontalAlignment','left','String', 'Scale Y','FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-160 50 20]);
    EPscanCont.handles.topos.sacScaleY = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.2f', EPscanCont.sacScaleY),'FontSize',EPmain.fontsize,...
        'Position',[600+(6*EPscanCont.topoSize)+50 EPscanCont.topoSize+10-180 50 20],'Callback',@changeEyeSettings);
end;

EPscanCont.handles.chanNames=[];

EPscanCont.handles.xSlider = uicontrol('Style', 'slider', 'Value',EPscanCont.xGraphPos,...
    'Position', [EPscanCont.graphCoords(1) EPscanCont.graphCoords(2)-20 EPscanCont.graphCoords(3) 20],'SliderStep',[EPscanCont.displayLength/EPscanCont.numPoints max(0.2,length(EPscanCont.displayLength)/EPscanCont.numPoints)], 'Min',0,'Max',1,'Callback', @xSlider);

if EPscanCont.numPoints <= EPscanCont.displayLength
    set(EPscanCont.handles.xSlider,'enable','off');
end;

EPscanCont.handles.ySlider = uicontrol('Style', 'slider', 'Value',EPscanCont.yGraphPos,...
    'Position', [EPscanCont.graphCoords(1)-20 EPscanCont.graphCoords(2) 20 EPscanCont.graphCoords(4)],'SliderStep',[length(EPscanCont.displayChans)/EPscanCont.numChans max(0.2,length(EPscanCont.displayChans)/EPscanCont.numChans)], 'Min',0,'Max',1,'Callback', @ySlider);

if EPscanCont.numChans == length(EPscanCont.displayChans)
    set(EPscanCont.handles.ySlider,'enable','off');
end;

EPscanCont.handles.epochMarks = axes('Units','pixel','Position', [EPscanCont.graphCoords(1) EPscanCont.graphCoords(2)-20 EPscanCont.graphCoords(3) 20]);
axis([1 EPscanCont.displayLength 0 1])
set(EPscanCont.handles.epochMarks,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

EPscanCont.handles.events=[];

EPscanCont.handles.editMode = uicontrol('Style', 'togglebutton', 'Value', EPscanCont.editMode,'String','edt','FontSize',EPmain.fontsize,...
    'Position', [100 scrsz(4)-100 30 30], 'Callback', @editMode);

EPscanCont.handles.eventLabels = uicontrol('Style', 'togglebutton', 'Value', EPscanCont.eventLabels,'String','evt','FontSize',EPmain.fontsize*2,...
    'Position', [130 scrsz(4)-100 30 30], 'Callback', @eventLabels);

EPscanCont.handles.displayStart=uicontrol('Style','text','HorizontalAlignment','left','String', ep_ms2time((EPscanCont.displayPoints(1)-1)*(1000/EPscanCont.EPdata(1).Fs)),'FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1) 160 70 20]);

EPscanCont.handles.displayEnd=uicontrol('Style','text','HorizontalAlignment','left','String', ep_ms2time(EPscanCont.displayPoints(end)*(1000/EPscanCont.EPdata(1).Fs)),'FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(3)+EPscanCont.graphCoords(1)-70 160 70 20]);

if ~isempty(EPscanCont.currentEvent)
    EPscanCont.handles.currentEvent = uicontrol('Style', 'popupmenu', 'String', EPscanCont.starredSampleEventList, 'FontSize', EPmain.fontsize,...
        'Callback', @changeCurrentEvent,...
        'Value', EPscanCont.currentEvent,'Position', [EPscanCont.graphCoords(1) 140 350 20]);
else
    EPscanCont.handles.currentEvent = uicontrol('Style', 'popupmenu', 'String','none', 'FontSize', EPmain.fontsize,...
        'Callback', @changeCurrentEvent,...
        'Value', 1,'Position', [EPscanCont.graphCoords(1) 140 350 20]);
end;

EPscanCont.handles.recenter = uicontrol('Style', 'pushbutton', 'String', 'Recenter','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1) 120 100 20], 'Callback', @recenter);

EPscanCont.handles.backEvent = uicontrol('Style', 'pushbutton', 'String', '<-','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+110 120 20 20], 'Callback', @shiftEvent);

EPscanCont.handles.nextEvent = uicontrol('Style', 'pushbutton', 'String', '->','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+130 120 20 20], 'Callback', @shiftEvent);

EPscanCont.handles.e1 = uicontrol('Style', 'pushbutton', 'String', 'e1','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+160 120 20 20], 'Callback', @e1e2);

EPscanCont.handles.e2 = uicontrol('Style', 'pushbutton', 'String', 'e2','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+180 120 20 20], 'Callback', @e1e2);

EPscanCont.handles.firstEvent = uicontrol('Style', 'popupmenu', 'String', EPscanCont.starredEventList, 'FontSize', EPmain.fontsize,...
    'Callback', @changeEvent,...
    'Value', EPscanCont.firstEvent,'Position', [EPscanCont.graphCoords(1) 100 350 20]);

EPscanCont.handles.lastEvent = uicontrol('Style', 'popupmenu', 'String', EPscanCont.starredEventList, 'FontSize', EPmain.fontsize,...
    'Callback', @changeEvent,...
    'Value', EPscanCont.lastEvent,'Position', [EPscanCont.graphCoords(1) 80 350 20]);

EPscanCont.handles.redisplay = uicontrol('Style', 'pushbutton', 'String', 'Redisplay','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1) 60 100 20], 'Callback', @redisplay);

EPscanCont.handles.addEvent = uicontrol('Style', 'pushbutton', 'String', '+ evt','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 170 50 20], 'Callback', @editEvent);

EPscanCont.handles.minusEvent = uicontrol('Style', 'pushbutton', 'String', '- evt','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+460 170 50 20], 'Callback', @editEvent);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Type','FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1)+360 150 50 20]);

EPscanCont.handles.eventType = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 150 150 20], 'Callback', @editEvent);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Value','FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1)+360 130 50 20]);

EPscanCont.handles.eventValue = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 130 150 20], 'Callback', @editEvent);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Sample','FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1)+360 110 50 20]);

EPscanCont.handles.eventSample = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 110 150 20], 'Callback', @editEvent);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Duration','FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1)+360 90 50 20]);

EPscanCont.handles.eventDuration = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 90 150 20], 'Callback', @editEvent);

EPscanCont.handles.eventKeys = uicontrol('Style', 'popupmenu', 'String', {''},'Value', EPscanCont.currentkey, 'FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+405 70 60 20], 'Callback', @editEvent);

EPscanCont.handles.addKey = uicontrol('Style', 'pushbutton', 'String', '+ key','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+460 70 50 20], 'Callback', @editEvent);

EPscanCont.handles.minusKey = uicontrol('Style', 'pushbutton', 'String', '- key','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+510 70 50 20], 'Callback', @editEvent);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Code','FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1)+360 50 50 20]);

EPscanCont.handles.keyCode = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 50 150 20], 'Callback', @editEvent);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Data','FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1)+360 30 50 20]);

EPscanCont.handles.keyData = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 30 150 20], 'Callback', @editEvent);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Type','FontSize',EPmain.fontsize,...
    'Position',[EPscanCont.graphCoords(1)+360 10 50 20]);

EPscanCont.handles.keyType = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+410 10 75 20], 'Callback', @editEvent);

EPscanCont.handles.keyDescrip = uicontrol('Style', 'edit', 'String', '','FontSize',EPmain.fontsize,...
    'Position', [EPscanCont.graphCoords(1)+485 10 75 20], 'Callback', @editEvent);

if isempty(EPscanCont.eventList)
    set(EPscanCont.handles.firstEvent,'enable','off');
    set(EPscanCont.handles.lastEvent,'enable','off');
    set(EPscanCont.handles.redisplay,'enable','off');
    set(EPscanCont.handles.currentEvent,'enable','off');
end;

if ~isempty(EPscanCont.EYEsize) || ~isempty(EPscanCont.SACCsize)
    EPscanCont.handles.freqList = uicontrol('Style', 'popupmenu', 'String', EPscanCont.freqList, 'FontSize', EPmain.fontsize,...
        'Callback', @changeFreq,...
        'Value', EPscanCont.freqShow,'Position', [92 140 100 20]);
    
    if strcmp(EPscanCont.dataType,'VLT')
        set(EPscanCont.handles.freqList,'enable','off');
    end;
    
    EPscanCont.handles.chanList = uicontrol('Style', 'popupmenu', 'String', EPscanCont.chanList, 'FontSize', EPmain.fontsize,...
        'Callback', @changeChan,...
        'Value', EPscanCont.chanShow,'Position', [92 120 100 20]);
end;

EPscanCont.handles.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
    'Position', [92 5 50 35], 'Callback', @cancel);

EPscanCont.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Keep','FontSize',EPmain.fontsize,...
    'Position', [152 5 50 35], 'Callback', @done);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Channels','FontSize',EPmain.fontsize,...
    'Position',[212 40 60 20]);

EPscanCont.handles.lessChans = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*2,...
    'Position', [212 10 30 30], 'Callback', @lessChans);

EPscanCont.handles.moreChans = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*2,...
    'Position', [242 10 30 30], 'Callback', @moreChans);

set(EPscanCont.handles.moreChans,'enable','off');

uicontrol('Style','text','HorizontalAlignment','left','String', 'Time','FontSize',EPmain.fontsize,...
    'Position',[292 40 60 20]);

EPscanCont.handles.lessTime = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*2,...
    'Position', [292 10 30 30], 'Callback', @lessTime);

EPscanCont.handles.moreTime = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*2,...
    'Position', [322 10 30 30], 'Callback', @moreTime);

EPscanCont.handles.play = uicontrol('Style', 'pushbutton', 'String', 'Play','FontSize',EPmain.fontsize*2,...
    'Position', [372 10 50 30], 'Callback', @startPlay);

EPscanCont.handles.stop = uicontrol('Style', 'pushbutton', 'String', 'Stop','FontSize',EPmain.fontsize*2,...
    'Position', [422 10 50 30], 'Callback', ['global EPscanCont;','EPscanCont.playMode=false;' ]);

EPscanCont.handles.movie = uicontrol('Style', 'pushbutton', 'String', 'Mov','FontSize',EPmain.fontsize*2,...
    'Position', [472 10 50 30], 'Callback', @startMovie);

set(EPmain.handles.scanContData,'WindowButtonDownFcn',@clickGraph);

updateDisplay(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cancel(~,~)
%quit from trim data window without saving edits

global EPmain EPscanCont 

close(EPmain.handles.scanContData);

EPscanCont=[];

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function done(~,~)
%quit from trim data window and save edits

global EPmain EPscanCont EPdataset

close(EPmain.handles.scanContData);

if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
    msg{1}='The work directory cannot be found.';
    [msg]=ep_errorMsg(msg);
    return
end;

if EPscanCont.edited
    for iChan=1:EPscanCont.numChans
        theChan=find(EPscanCont.chanIndex{1}==EPscanCont.displayChans(iChan));
        EPscanCont.EPdata(1).analysis.badChans(1,:,theChan)=EPscanCont.trialBadChans(iChan,:);
        if EPscanCont.globalBadChans(iChan)
            EPscanCont.EPdata(1).analysis.badChans(1,:,theChan)=1;
        end;
    end;
    EPscanCont.EPdata(1).analysis.badTrials=EPscanCont.badTrials;
    
    if ~isempty(EPscanCont.XEYchan) && ~isempty(EPscanCont.YEYchan)
        EPscanCont.EPdata(1).calibration.ET.Xzero=EPscanCont.eyeCenterX;
        EPscanCont.EPdata(1).calibration.ET.Yzero=EPscanCont.eyeCenterY;
        EPscanCont.EPdata(1).calibration.ET.Xscale=EPscanCont.eyeScaleX;
        EPscanCont.EPdata(1).calibration.ET.Yscale=EPscanCont.eyeScaleY;
    end;
    if ~isempty(EPscanCont.HSACchan) && ~isempty(EPscanCont.VSACchan)
        EPscanCont.EPdata(1).calibration.SAC.Xzero=EPscanCont.sacCenterX;
        EPscanCont.EPdata(1).calibration.SAC.Yzero=EPscanCont.sacCenterY;
        EPscanCont.EPdata(1).calibration.SAC.Xscale=EPscanCont.sacScaleX;
        EPscanCont.EPdata(1).calibration.SAC.Yscale=EPscanCont.sacScaleY;
    end;
    
    EPscanCont.EPdata(1).pca=EPscanCont.pca{1};

    delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(EPmain.view.dataset(1)).dataName '.mat']);
    EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],EPmain.view.dataset(1)));
    EPdataset=ep_saveEPdataset(EPdataset,EPscanCont.EPdata(1),EPmain.view.dataset(1),'no');
end;

EPscanCont=[];

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xSlider(~,~)
%update graph after xSlider moved

global EPscanCont EPmain EPdataset

EPscanCont.xGraphPos=get(EPscanCont.handles.xSlider,'Value');

EPscanCont.displayPoints=1+floor(EPscanCont.xGraphPos*EPscanCont.numPoints):EPscanCont.displayLength+floor(EPscanCont.xGraphPos*EPscanCont.numPoints);
if max(EPscanCont.displayPoints) > EPscanCont.numPoints
    EPscanCont.displayPoints=[max(1,EPscanCont.numPoints-EPscanCont.displayLength+1):EPscanCont.numPoints];
end;

updateDisplay(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ySlider(~,~)
%update graph after ySlider moved

global EPscanCont EPmain EPdataset

EPscanCont.yGraphPos=get(EPscanCont.handles.ySlider,'Value');

EPscanCont.displayChans=1+floor((1-EPscanCont.yGraphPos)*EPscanCont.numChans):length(EPscanCont.displayChans)+floor((1-EPscanCont.yGraphPos)*EPscanCont.numChans);
if max(EPscanCont.displayChans) > EPscanCont.numChans
    EPscanCont.displayChans=[max(1,EPscanCont.numChans-length(EPscanCont.displayChans)+1):EPscanCont.numChans];
end;

updateDisplay(0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lessChans(~,~)
%zoom in (less channels)

global EPscanCont EPmain EPdataset

if length(EPscanCont.displayChans) > 1
    EPscanCont.displayChans=[1:floor(length(EPscanCont.displayChans)/2)];
end

if EPscanCont.numChans > length(EPscanCont.displayChans)
    set(EPscanCont.handles.ySlider,'enable','on');
end;

EPscanCont.yGraphPos=1;
set(EPscanCont.handles.ySlider,'Value',EPscanCont.yGraphPos);
sliderStep=2*length(EPscanCont.displayChans)/EPscanCont.numChans;
set(EPscanCont.handles.ySlider,'SliderStep',[sliderStep/4 sliderStep]);

set(EPscanCont.handles.moreChans,'enable','on');
if length(EPscanCont.displayChans) == 1
    set(EPscanCont.handles.lessChans,'enable','off');
else
    set(EPscanCont.handles.lessChans,'enable','on');
end;

updateDisplay(0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function moreChans(~,~)
%zoom out (more channels)

global EPscanCont

EPscanCont.displayChans=[1:length(EPscanCont.displayChans)*2];
if max(EPscanCont.displayChans) > EPscanCont.numChans
    EPscanCont.displayChans=[1:EPscanCont.numChans];
end;

EPscanCont.yGraphPos=1;
set(EPscanCont.handles.ySlider,'Value',EPscanCont.yGraphPos);
sliderStep=2*length(EPscanCont.displayChans)/EPscanCont.numChans;
set(EPscanCont.handles.ySlider,'SliderStep',[sliderStep/4 sliderStep]);

if EPscanCont.numChans == length(EPscanCont.displayChans)
    set(EPscanCont.handles.ySlider,'enable','off');
end;

set(EPscanCont.handles.lessChans,'enable','on');
if length(EPscanCont.displayChans) == EPscanCont.numChans
    set(EPscanCont.handles.moreChans,'enable','off');
else
    set(EPscanCont.handles.moreChans,'enable','on');
end;

updateDisplay(0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lessTime(~,~)
%zoom in (less time)

global EPscanCont

if length(EPscanCont.displayPoints) > EPscanCont.xGraphMin
    EPscanCont.displayPoints=[EPscanCont.displayPoints(1):EPscanCont.displayPoints(1)+floor(length(EPscanCont.displayPoints)/2)-1];
end;

if length(EPscanCont.displayPoints) < EPscanCont.xGraphMin
    EPscanCont.displayPoints=[EPscanCont.displayPoints(1):EPscanCont.displayPoints(1)+EPscanCont.xGraphMin-1];
end;

EPscanCont.plotPoint=ceil(EPscanCont.plotPoint*2);
if EPscanCont.plotPoint > length(EPscanCont.displayPoints)
    EPscanCont.plotPoint = length(EPscanCont.displayPoints);
end;

EPscanCont.displayLength=length(EPscanCont.displayPoints);
set(EPscanCont.handles.xSlider,'SliderStep',[EPscanCont.displayLength/EPscanCont.numPoints max(0.2,length(EPscanCont.displayLength)/EPscanCont.numPoints)]);

if EPscanCont.numPoints > length(EPscanCont.displayPoints)
    set(EPscanCont.handles.xSlider,'enable','on');
end;

set(EPscanCont.handles.moreTime,'enable','on');
if length(EPscanCont.displayPoints) == EPscanCont.xGraphMin
    set(EPscanCont.handles.lessTime,'enable','off');
else
    set(EPscanCont.handles.lessTime,'enable','on');
end;

updateDisplay(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function moreTime(~,~)
%zoom out (more time)

global EPscanCont

EPscanCont.displayPoints=[EPscanCont.displayPoints(1):EPscanCont.displayPoints(1)+length(EPscanCont.displayPoints)*2-1];

if max(EPscanCont.displayPoints) > EPscanCont.numPoints
    EPscanCont.displayPoints=[EPscanCont.displayPoints(1):EPscanCont.numPoints];
end;

EPscanCont.plotPoint=ceil(EPscanCont.plotPoint/2);

EPscanCont.displayLength=length(EPscanCont.displayPoints);
set(EPscanCont.handles.xSlider,'SliderStep',[EPscanCont.displayLength/EPscanCont.numPoints max(0.2,length(EPscanCont.displayLength)/EPscanCont.numPoints)]);

if EPscanCont.numPoints == length(EPscanCont.displayPoints)
    set(EPscanCont.handles.xSlider,'enable','off');
end;

set(EPscanCont.handles.lessTime,'enable','on');
if length(EPscanCont.displayPoints) == EPscanCont.numPoints
    set(EPscanCont.handles.moreTime,'enable','off');
else
    set(EPscanCont.handles.moreTime,'enable','on');
end;

updateDisplay(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clickGraph(~,~)
%respond to click in graph area

global EPscanCont EPmain

pos = get(EPmain.handles.scanContData, 'CurrentPoint');
if (pos(1) < EPscanCont.graphCoords(1)) || (pos(1) > EPscanCont.graphCoords(1)+EPscanCont.graphCoords(3))
    return
end;
if (pos(2) < (EPscanCont.graphCoords(2))-20) || (pos(2) > EPscanCont.graphCoords(2)+EPscanCont.graphCoords(4))
    return
end;

clickX=ceil(((pos(1)-EPscanCont.graphCoords(1))/EPscanCont.graphCoords(3))*EPscanCont.displayLength); %convert to graph samples
if clickX < 1
    clickX=1;
end;
if clickX > EPscanCont.displayLength
    clickX=EPscanCont.displayLength;
end;

clickY=((EPmain.scrsz(4)-pos(2)-(EPmain.scrsz(4)-EPscanCont.graphCoords(2)-EPscanCont.graphCoords(4)))/EPscanCont.graphCoords(4))*(-(length(EPscanCont.displayChans)+1)*EPscanCont.spacing)+EPscanCont.spacing; %convert to graph microvolts

axes(EPscanCont.handles.graphPlot);
if ~strcmpi(get(EPmain.handles.scanContData,'selectiontype'),'normal')
    %right-click means move the green line
       
    delete(EPscanCont.handles.plotLine);
    
    EPscanCont.plotPoint=clickX;
    
    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
        EPscanCont.handles.plotLine=line([EPscanCont.plotPoint EPscanCont.plotPoint],[-length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing],'Color','green','LineWidth',3);
    else
        EPscanCont.handles.plotLine=line([EPscanCont.plotPoint EPscanCont.plotPoint],[1 EPscanCont.spacing*length(EPscanCont.displayChans)],'Color','green','LineWidth',3);
    end;
    
    drawPlots
    
    updateFreqGraph
else
    %left-click means edit
    theEpoch=ceil((EPscanCont.displayPoints(clickX))/EPscanCont.EPdata(1).Fs);
    graphEpoch=ceil(EPscanCont.displayPoints(clickX)/EPscanCont.EPdata(1).Fs)-floor(EPscanCont.displayPoints(1)/EPscanCont.EPdata(1).Fs);
    if theEpoch>length(EPscanCont.badTrials)
        %extra time points tacked onto last epoch
        theEpoch=theEpoch-1;
        graphEpoch=graphEpoch-1;
        endEpoch=EPscanCont.displayLength;
    elseif (theEpoch==length(EPscanCont.badTrials)) && rem(EPscanCont.numPoints,EPscanCont.EPdata(1).Fs)
        endEpoch=EPscanCont.displayLength;
    else
        endEpoch=graphEpoch*EPscanCont.EPdata(1).Fs-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
    end;
    startEpoch=((graphEpoch-1)*EPscanCont.EPdata(1).Fs)+1-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
    if (pos(2) < EPscanCont.graphCoords(2)) && (pos(2) >= (EPscanCont.graphCoords(2)-20))
        %bad trial edit
        if EPscanCont.badTrials(theEpoch)
            EPscanCont.badTrials(theEpoch)=0;
            delete(EPscanCont.handles.badTrials(theEpoch));
        else
            EPscanCont.badTrials(theEpoch)=1;
            EPscanCont.handles.badTrials(theEpoch)=patch([startEpoch startEpoch endEpoch endEpoch],...
                [-EPscanCont.spacing*length(EPscanCont.displayChans) EPscanCont.spacing*length(EPscanCont.displayChans) EPscanCont.spacing*length(EPscanCont.displayChans) -EPscanCont.spacing*length(EPscanCont.displayChans)],'red','EdgeColor','red');
            alpha(EPscanCont.handles.badTrials(theEpoch),.5);
        end;
    else
        %trial bad chan edit
        if EPscanCont.editMode
            theChan=ceil(-(clickY-25)/EPscanCont.spacing);
            if theChan < 1
                theChan=1;
            elseif theChan > length(EPscanCont.displayChans)
                theChan=length(EPscanCont.displayChans);
            end;
            theChan=EPscanCont.displayChans(theChan);
            if EPscanCont.trialBadChans(theChan,theEpoch) ==-1
                EPscanCont.trialBadChans(theChan,theEpoch)=0;
                delete(EPscanCont.handles.trialBadChans(theChan,theEpoch));
            else
                EPscanCont.trialBadChans(theChan,theEpoch)=-1;
                axes(EPscanCont.handles.graphPlot);
                EPscanCont.handles.trialBadChans(theChan,theEpoch)=patch([startEpoch startEpoch endEpoch endEpoch],...
                    [-EPscanCont.spacing*(theChan-.25) -EPscanCont.spacing*(theChan-1.25) -EPscanCont.spacing*(theChan-1.25) -EPscanCont.spacing*(theChan-.25)],'red','EdgeColor','red');
                alpha(EPscanCont.handles.trialBadChans(theChan,theEpoch),.5);
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drawPlots
%Redraw the plots

global EPscanCont EPmain EPdataset

if ~isempty(EPscanCont.plotPoint)
    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
        if ~isempty(EPscanCont.hasLoc)
            plotData=zeros(length(EPscanCont.hasLoc),EPscanCont.numColors);
            for iColor=1:EPscanCont.numColors
                if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                    plotData(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.hasLoc)),iColor)=EPscanCont.EPdata(iColor).data(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.hasLoc)),EPscanCont.displayPoints(EPscanCont.plotPoint),:,:,:,EPscanCont.freqPos);
                    plotData(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.hasLoc)),iColor)=plotData(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.hasLoc)),iColor)-EPscanCont.baseline(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.hasLoc)),iColor);
                end;
            end;
            
            if strcmp(EPscanCont.dataType,'TFT')
                plotData=abs(plotData); %convert complex number to real number
                plotData=plotData/mean(diff(EPscanCont.EPdata(1).freqNames)); %convert to spectral density
                plotData=plotData.^2; %convert amplitude to power
                plotData=log10(abs(plotData))*10; %convert to dB log scaling
                tempVar=plotData;
                tempVar(isinf(tempVar))=-flintmax;
                plotData=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
            end;
            theMin=inf;
            theMax=-inf;
            for iColor=1:EPscanCont.numColors
                if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                    theMin=min(theMin,min(plotData(:,iColor)));
                    theMax=max(theMax,max(plotData(:,iColor)));
                end;
            end;
            
            for iColor=1:EPscanCont.numColors
                if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                    theEpoch=ceil((EPscanCont.displayPoints(EPscanCont.plotPoint))/EPscanCont.EPdata(1).Fs);
                    if theEpoch > size(EPscanCont.EPdata(iColor).analysis.badChans,2)
                        theEpoch=size(EPscanCont.EPdata(iColor).analysis.badChans,2); %the stub is tacked onto the last epoch
                    end;
                    goodChans=find(~ismember(EPscanCont.hasLoc,EPscanCont.chanIndex{iColor}(find(EPscanCont.EPdata(iColor).analysis.badChans(1,theEpoch,:)==-1))));
                    if ~isempty(goodChans)
                        theData=plotData(EPscanCont.hasLoc(goodChans),iColor);
                        [Xi,Yi,Zi] = griddata(EPscanCont.x(goodChans),EPscanCont.y(goodChans),theData,[1:EPscanCont.gridSize]',[1:EPscanCont.gridSize],'linear');
                        
                        axes(EPscanCont.handles.topos.topo(iColor))
                        EPscanCont.handles.topoImage(iColor) = imagesc(Zi);
                        set(EPscanCont.handles.topos.topo(iColor),'XTickLabel','','YTickLabel','');
                        set(EPscanCont.handles.topoImage(iColor),'ButtonDownFcn',{@D3headPlot,iColor});
                    end;
                end;
            end;
        end;
    end;
    
    if ~isempty(EPscanCont.EYEsize) || ~isempty(EPscanCont.SACCsize)
        if EPscanCont.chanShow==length(EPscanCont.chanList)
            colorData=zeros(length(EPscanCont.displayPoints),3);
        else
            theData=squeeze(EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.chanShow),EPscanCont.displayPoints,:,:,:,EPscanCont.freqPos))';
            if strcmp(EPscanCont.dataType,'TFT')
                theData=abs(theData); %convert complex number to real number
                theData=theData/mean(diff(EPscanCont.EPdata(1).freqNames)); %convert to spectral density
                theData=theData.^2; %convert amplitude to power
                theData=log10(abs(theData))*10; %convert to dB log scaling
                tempVar=theData;
                tempVar(isinf(tempVar))=-flintmax;
                theData=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
            end;
            colorData(:,1)=(theData-mean(theData))/max(theData-mean(theData));
            colorData(:,3)=(theData-mean(theData))/min(theData-mean(theData));
            colorData(colorData<0)=0;
        end;
    end;
    
    stimPict=[];
    if ~isempty(EPscanCont.stimList)
        whichPict=max(find([EPscanCont.stimList.sample] <= EPscanCont.displayPoints(EPscanCont.plotPoint)));
        if ~isempty(whichPict)
            stimPict=EPscanCont.EPdata(1).stims(EPscanCont.stimList(whichPict).stim).image;
        end;
    end;
    
    if ~isempty(EPscanCont.EYEsize)
        Xdata=squeeze(EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.XEYchan),EPscanCont.displayPoints,:,:,:,1));
        Ydata=squeeze(EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.YEYchan),EPscanCont.displayPoints,:,:,:,1));
        if EPscanCont.eyeScaleX < 0
            Xdata=((Xdata-EPscanCont.eyeCenterX)/(-EPscanCont.eyeScaleX))+.5;
            Xdata=1-Xdata;
        else
            Xdata=((Xdata-EPscanCont.eyeCenterX)/EPscanCont.eyeScaleX)+.5;
        end;
        if EPscanCont.eyeScaleX < 0
            Ydata=((Ydata-EPscanCont.eyeCenterY)/(-EPscanCont.eyeScaleY))+.5;
            Ydata=1-Ydata;
        else
            Ydata=((Ydata-EPscanCont.eyeCenterY)/EPscanCont.eyeScaleY)+.5;
        end;
        cla(EPscanCont.handles.topos.eye);
        axes(EPscanCont.handles.topos.eye);
        axis([0 1 0 1])
        set(EPscanCont.handles.topos.eye,'Units','normalized')
        hold on
        if ~isempty(stimPict)
            imageSize=size(stimPict);
            imageMax=max(imageSize(1),imageSize(2));
            x=max(1,round((imageSize(1)-imageSize(2))/2))/imageMax;
            y=max(1,round((imageSize(2)-imageSize(1))/2))/imageMax;
            EPscanCont.handles.topos.eyeImage=image([x 1-x],[1-y y],stimPict);
            set(EPscanCont.handles.topos.eyeImage,'ButtonDownFcn',@ep_expandEyePlot);
        end;
        if EPscanCont.chanShow~=length(EPscanCont.chanList)
            for iPoint=1:length(EPscanCont.displayPoints)-1
                plot(Xdata(iPoint:iPoint+1),Ydata(iPoint:iPoint+1),'color',colorData(iPoint,:));
            end;
        else
            plot(Xdata,Ydata,'color','black');
        end;
        plot(Xdata(EPscanCont.plotPoint),Ydata(EPscanCont.plotPoint),'color','green','Marker','*','MarkerSize',2);
        hold off
        set(EPscanCont.handles.topos.eye,'ButtonDownFcn',@ep_expandEyePlot);
    end;
    
    if ~isempty(EPscanCont.SACCsize)
        Xdata=squeeze(EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.HSACchan),EPscanCont.displayPoints,:,:,:,1));
        Ydata=squeeze(EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.VSACchan),EPscanCont.displayPoints,:,:,:,1));
        Xdata=(-(Xdata-EPscanCont.sacCenterX)/EPscanCont.sacScaleX)+.5;
        Ydata=(-(Ydata-EPscanCont.sacCenterY)/EPscanCont.sacScaleY)+.5;
        cla(EPscanCont.handles.topos.saccade);
        axes(EPscanCont.handles.topos.saccade);
        axis([0 1 0 1])
        set(EPscanCont.handles.topos.saccade,'Units','normalized')
        hold on
        if ~isempty(stimPict)
            imageSize=size(stimPict);
            imageMax=max(imageSize(1),imageSize(2));
            x=max(1,round((imageSize(1)-imageSize(2))/2))/imageMax;
            y=max(1,round((imageSize(2)-imageSize(1))/2))/imageMax;
            EPscanCont.handles.topos.saccImage=image([x 1-x],[1-y y],stimPict);
            set(EPscanCont.handles.topos.saccImage,'ButtonDownFcn',@ep_expandSaccPlot);
        end;
        if EPscanCont.chanShow~=length(EPscanCont.chanList)
            for iPoint=1:length(EPscanCont.displayPoints)-1
                plot(Xdata(iPoint:iPoint+1),Ydata(iPoint:iPoint+1),'color',colorData(iPoint,:));
            end;
        else
            plot(Xdata,Ydata,'color','black');
        end;
        plot(Xdata(EPscanCont.plotPoint),Ydata(EPscanCont.plotPoint),'color','green','Marker','*','MarkerSize',2);
        hold off
        set(EPscanCont.handles.topos.saccade,'ButtonDownFcn',@ep_expandSaccPlot);
    end;
    
    %update the voltage read-outs
    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
        for iChan=1:length(EPscanCont.displayChans)
            set(EPscanCont.handles.plotVolt(iChan),'String',sprintf('%04.0f',EPscanCont.EPdata(1).data(find(EPscanCont.chanIndex{1}==EPscanCont.displayChans(iChan)),EPscanCont.displayPoints(EPscanCont.plotPoint))));
        end;
    end;
else
    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
        if ~isempty(EPscanCont.hasLoc)
            for iColor=1:EPscanCont.numColors
                if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                    cla(EPscanCont.handles.topos.topo(iColor))
                end;
            end;
        end;
    end;
    if ~isempty(EPscanCont.EYEsize)
        cla(EPscanCont.handles.topos.eye);
    end;
    if ~isempty(EPscanCont.SACCsize)
        cla(EPscanCont.handles.topos.saccade);
    end;
    %clear the voltage read-outs
    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
        for iChan=1:length(EPscanCont.displayChans)
            set(EPscanCont.handles.plotVolt(iChan),'String','');
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function startPlay(~,~)
%Start play mode

global EPscanCont EPdataset EPmain

EPscanCont.playMode=true;
playRate=(4/EPscanCont.sampleSize); %time points advanced per cycle
xChange=0;

while EPscanCont.playMode
    advancePoints=playRate;
    if EPscanCont.numPoints >= (EPscanCont.displayPoints(end)+playRate)
        EPscanCont.displayPoints=EPscanCont.displayPoints+advancePoints;
        xChange=1;
    elseif (EPscanCont.numPoints < (EPscanCont.displayPoints(end)+playRate)) && (EPscanCont.numPoints ~= EPscanCont.displayPoints(end))
        advancePoints=EPscanCont.numPoints-EPscanCont.displayPoints(end);
        EPscanCont.displayPoints=EPscanCont.displayPoints+advancePoints;
        xChange=1;
    elseif length(EPscanCont.displayPoints) >= EPscanCont.EPscanCont.plotPoint+advancePoints
        EPscanCont.plotPoint=EPscanCont.plotPoint+advancePoints;
    else
        EPscanCont.playMode=false;
    end;
    
    updateDisplay(xChange,0);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function startMovie(~,~)
%Start recording movie

global EPscanCont
[outFileName, pathname] = uiputfile('*.*','Save:','myVideo.avi');
tic
numPoints=length(EPscanCont.displayPoints);
%EPscanCont.movie.ET=struct('cdata', cell(numPoints,1), 'colormap', cell(numPoints,1));
writerObj=VideoWriter([pathname outFileName]);
writerObj.FrameRate=1;
open(writerObj);
for iPoint=1:numPoints
    if ~isempty(EPscanCont.EYEsize)
        ep_expandEyePlotFrame(iPoint);
%        EPscanCont.movie.ET(iPoint)=getframe(EPscanCont.handles.waves.frameFigure);
        theFrame=getframe(EPscanCont.handles.waves.frameFigure);
        writeVideo(writerObj,theFrame);
        close(EPscanCont.handles.waves.frameFigure)
    end;
end;
% for iPoint=1:numPoints
%     writeVideo(writerObj,EPscanCont.movie.ET(iPoint));
% end;
close(writerObj);
subjectTime = floor(toc/60);
disp(['The data took ' num2str(subjectTime) ' minutes to generate an ERP tracker movie.']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateDisplay(xChange,yChange)
%Redraw the display

global EPscanCont EPmain EPdataset

if yChange
    if ~isempty(EPscanCont.handles.chanNames)
        for iChan=1:length(EPscanCont.handles.chanNames)
            delete(EPscanCont.handles.chanNames(iChan));
            delete(EPscanCont.handles.plotVolt(iChan));
        end;
    end;
    
    theColor='black';
    EPscanCont.handles.chanNames=[];
    for iChan=1:length(EPscanCont.displayChans)
        thChan=EPscanCont.displayChans(iChan);
        EPscanCont.handles.chanNames(iChan)=uicontrol('Style','togglebutton','HorizontalAlignment','left','String', EPscanCont.chanNames{thChan},'FontSize',EPmain.fontsize,...
            'ForegroundColor',theColor,'Value', EPscanCont.globalBadChans(thChan),'Callback', {@globalBadChans,thChan},...
            'Position',[110 EPmain.scrsz(4)-120-EPscanCont.graphCoords(4)*(1/length(EPscanCont.displayChans))*(iChan-1) 50 20]);
        EPscanCont.handles.plotVolt(iChan)=uicontrol('Style','text','HorizontalAlignment','left','String', '0','FontSize',EPmain.fontsize,...
            'Position',[EPscanCont.graphCoords(3)+EPscanCont.graphCoords(1)+5 EPmain.scrsz(4)-120-EPscanCont.graphCoords(4)*(1/length(EPscanCont.displayChans))*(iChan-1) 40 20]);
    end;
end;

%update the current data matrix
if xChange || yChange
    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
        EPscanCont.currentData=zeros(length(EPscanCont.displayChans),length(EPscanCont.displayPoints),4);
    else
        EPscanCont.currentData=zeros(length(EPscanCont.displayChans)*EPscanCont.spacing,length(EPscanCont.displayPoints),4);
    end;
    for iColor=1:EPscanCont.numColors
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            for iChan=1:length(EPscanCont.displayChans)
                theChan=EPscanCont.displayChans(iChan);
                theDataChan=find(EPscanCont.chanIndex{iColor}==theChan);
                if isempty(theDataChan) 
                    EPscanCont.baseline(iChan,iColor)=0;
                    EPscanCont.currentData(iChan,:,iColor)=0;
                else
                    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
                        waveData=squeeze(EPscanCont.EPdata(iColor).data(theDataChan,EPscanCont.displayPoints,:,:,:,EPscanCont.freqPos));
                        goodPoints=~isnan(waveData);
                        EPscanCont.baseline(iChan,iColor)=mean(waveData(goodPoints));
                        EPscanCont.currentData(iChan,:,iColor)=((EPscanCont.EPdata(iColor).data(theDataChan,EPscanCont.displayPoints,:,:,:,EPscanCont.freqPos)-EPscanCont.baseline(iChan,iColor))*EPscanCont.scaling(theChan));
                    elseif EPscanCont.FFTunits==1 %complex numbers
                        if ismember(theChan,EPscanCont.EEGchans)
                            EPscanCont.currentData((iChan-1)*EPscanCont.spacing+1:iChan*EPscanCont.spacing-(EPscanCont.spacing/2),:,iColor)=real(squeeze((EPscanCont.EPdata(iColor).data(theDataChan,EPscanCont.displayPoints,:,:,:,EPscanCont.startBins:EPscanCont.lastBins))*EPscanCont.scaling(theChan))');
                            EPscanCont.currentData((iChan-1)*EPscanCont.spacing+1+(EPscanCont.spacing/2):iChan*EPscanCont.spacing,:,iColor)=imag(squeeze((EPscanCont.EPdata(iColor).data(theDataChan,EPscanCont.displayPoints,:,:,:,EPscanCont.startBins:EPscanCont.lastBins))*EPscanCont.scaling(theChan))');
                        else
                            EPscanCont.currentData((iChan-1)*EPscanCont.spacing+1:iChan*EPscanCont.spacing,:,iColor)=repmat(real(squeeze((EPscanCont.EPdata(iColor).data(theDataChan,EPscanCont.displayPoints,:,:,:,1))*EPscanCont.scaling(theChan))),EPscanCont.spacing,1);
                        end;
                    else
                        if ismember(theChan,EPscanCont.EEGchans)
                            EPscanCont.currentData((iChan-1)*EPscanCont.spacing+1:iChan*EPscanCont.spacing,:,iColor)=squeeze((EPscanCont.EPdata(iColor).data(theDataChan,EPscanCont.displayPoints,:,:,:,EPscanCont.startBins:EPscanCont.lastBins))*EPscanCont.scaling(theChan))';
                        else
                            EPscanCont.currentData((iChan-1)*EPscanCont.spacing+1:iChan*EPscanCont.spacing,:,iColor)=repmat(squeeze((EPscanCont.EPdata(iColor).data(theDataChan,EPscanCont.displayPoints,:,:,:,1))*EPscanCont.scaling(theChan)),EPscanCont.spacing,1);
                        end;
                    end;
                end;
            end;
        end;
    end;
    
    if strcmp(EPscanCont.dataType,'TFT')
        if EPscanCont.freqShow == length(EPscanCont.freqList)
            theRows=logical(kron(ismember(EPscanCont.displayChans,find(ismember(EPdataset.dataset(EPmain.view.dataset(1)).chanTypes,{'EEG','REG'})))',ones(EPscanCont.spacing,1)));
        else
            theRows=ismember(EPscanCont.displayChans,find(ismember(EPdataset.dataset(EPmain.view.dataset(1)).chanTypes,{'EEG','REG'})))';
        end;
        if (EPscanCont.FFTunits > 1)
            EPscanCont.currentData(theRows,:,1)=abs(EPscanCont.currentData(theRows,:,1)); %convert complex number to real number
        end;
        EPscanCont.currentData(theRows,:,1)=EPscanCont.currentData(theRows,:,1)/sqrt(mean(diff(EPscanCont.EPdata(1).freqNames))); %convert to spectral density
        if EPscanCont.FFTunits > 2
            EPscanCont.currentData(theRows,:,1)=EPscanCont.currentData(theRows,:,1).^2; %convert amplitude to power
        end;
        if (EPscanCont.FFTunits == 4)
            EPscanCont.currentData(theRows,:,1)=log10(abs(EPscanCont.currentData(theRows,:,1)))*10; %convert to dB log scaling
            tempVar=EPscanCont.currentData(theRows,:,1);
            tempVar(isinf(tempVar))=-flintmax;
            EPscanCont.currentData(theRows,:,1)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
        end;
    end;
    
    if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
        for iColor=1:EPscanCont.numColors
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                for iChan=1:length(EPscanCont.displayChans)
                    EPscanCont.currentData(iChan,:,iColor)=EPscanCont.currentData(iChan,:,iColor)-EPscanCont.spacing*(iChan-1);
                end;
            end;
        end;
    end;
end;

if xChange && EPscanCont.eventLabels
    for iEvent=1:length(EPscanCont.handles.events)
        if isgraphics(EPscanCont.handles.events(iEvent))
            delete(EPscanCont.handles.events(iEvent));
        end;
    end;
    EPscanCont.handles.events=[];
end;
EPscanCont.handles.eventLines=[];

if xChange
    %add epochmarks below the graph plot
    axes(EPscanCont.handles.epochMarks);
    EPscanCont.handles.epochMarksPlot=plot([1:EPscanCont.displayLength],EPscanCont.epochMarks(EPscanCont.displayPoints),'color','black');
    axis([1 EPscanCont.displayLength 0 1])
    set(EPscanCont.handles.epochMarks,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
end;

%update the current event to somewhere around the middle of the display area
if xChange
    if ~isempty(EPscanCont.eventList)
        EPscanCont.starredEventList=cell(length(EPscanCont.eventList),1);
        EPscanCont.starredSampleEventList=cell(length(EPscanCont.eventList),1);
        for iEvent=1:length(EPscanCont.eventList)
            if ismember(round(EPscanCont.EPdata(1).events{1}(iEvent).sample),EPscanCont.displayPoints)
                EPscanCont.starredEventList{find(iEvent==EPscanCont.alphabetEventOrder)}=['*' EPscanCont.eventList{iEvent}];
                sampleOrderNumber=find(iEvent==EPscanCont.sampleEventOrder);
                EPscanCont.starredSampleEventList{sampleOrderNumber}=['*' EPscanCont.eventList{iEvent}];
            else
                EPscanCont.starredEventList{find(iEvent==EPscanCont.alphabetEventOrder)}=EPscanCont.eventList{iEvent};
                EPscanCont.starredSampleEventList{find(iEvent==EPscanCont.sampleEventOrder)}=EPscanCont.eventList{iEvent};
            end;
        end;
        set(EPscanCont.handles.firstEvent,'String',EPscanCont.starredEventList);
        set(EPscanCont.handles.lastEvent,'String',EPscanCont.starredEventList);
        set(EPscanCont.handles.currentEvent,'String',EPscanCont.starredSampleEventList);
        
        if isempty(EPscanCont.plotPoint)
            [a currentEvent]=min(abs([EPscanCont.EPdata(1).events{1}.sample]-EPscanCont.displayPoints(round(EPscanCont.displayLength/2))));
            if ~isempty(currentEvent)
                EPscanCont.currentEvent=find(currentEvent==EPscanCont.sampleEventOrder);
                EPscanCont.currentkey=length(EPscanCont.EPdata(1).events{1}(currentEvent).keys);
                EPscanCont.plotPoint=find(round(EPscanCont.EPdata(1).events{1}(currentEvent).sample)==EPscanCont.displayPoints);
                set(EPscanCont.handles.currentEvent,'Value',EPscanCont.currentEvent);
            else
                EPscanCont.plotPoint=round(EPscanCont.displayLength/2);
            end;
        end;
    else
        EPscanCont.plotPoint=round(EPscanCont.displayLength/2);
    end;
end;

%update the current event info
if isempty(EPscanCont.currentEvent)
    set(EPscanCont.handles.addEvent,'enable','on');
    set(EPscanCont.handles.minusEvent,'enable','off');
    set(EPscanCont.handles.eventType,'enable','off');
    set(EPscanCont.handles.eventValue,'enable','off');
    set(EPscanCont.handles.eventSample,'enable','off');
    set(EPscanCont.handles.eventDuration,'enable','off');
    set(EPscanCont.handles.eventKeys,'enable','off');
    set(EPscanCont.handles.addKey,'enable','off');
    set(EPscanCont.handles.minusKey,'enable','off');
    set(EPscanCont.handles.keyCode,'enable','off');
    set(EPscanCont.handles.keyData,'enable','off');
    set(EPscanCont.handles.keyType,'enable','off');
    set(EPscanCont.handles.keyDescrip,'enable','off');
else
    set(EPscanCont.handles.addEvent,'enable','on');
    set(EPscanCont.handles.minusEvent,'enable','on');
    set(EPscanCont.handles.eventType,'enable','on');
    set(EPscanCont.handles.eventValue,'enable','on');
    set(EPscanCont.handles.eventSample,'enable','on');
    set(EPscanCont.handles.eventDuration,'enable','on');
    set(EPscanCont.handles.eventKeys,'enable','on');
    set(EPscanCont.handles.addKey,'enable','on');
    set(EPscanCont.handles.minusKey,'enable','on');
    set(EPscanCont.handles.keyCode,'enable','on');
    set(EPscanCont.handles.keyData,'enable','on');
    set(EPscanCont.handles.keyType,'enable','on');
    set(EPscanCont.handles.keyDescrip,'enable','on');
    
    theCurrentEvent=EPscanCont.sampleEventOrder(EPscanCont.currentEvent);
    set(EPscanCont.handles.eventType,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).type);
    set(EPscanCont.handles.eventValue,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).value);
    set(EPscanCont.handles.eventSample,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample);
    set(EPscanCont.handles.eventDuration,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).duration);
    
    theString=cell(length(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys),1);
    for iKey=1:length(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys)
        theString{iKey}=num2str(iKey);
    end;
    set(EPscanCont.handles.eventKeys,'String',theString);
    set(EPscanCont.handles.eventKeys,'Value',EPscanCont.currentkey);
    
    if EPscanCont.currentkey
        keyCode=EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).code;
        keyData=EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).data;
        keyType=EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).datatype;
        keyDescrip=EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).description;
    else
        keyCode='';
        keyData='';
        keyType='';
        keyDescrip='';
    end;
    set(EPscanCont.handles.keyCode,'String',keyCode);
    set(EPscanCont.handles.keyData,'String',keyData);
    set(EPscanCont.handles.keyType,'String',keyType);
    set(EPscanCont.handles.keyDescrip,'String',keyDescrip);
end;

%redraw the data graph
cla(EPscanCont.handles.graphPlot);
axes(EPscanCont.handles.graphPlot);
if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
    for iColor=1:EPscanCont.numColors
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            for iChan=1:length(EPscanCont.displayChans)
                hold on
                EPscanCont.handles.graphChans{iColor}=plot([1:EPscanCont.displayLength],squeeze(EPscanCont.currentData(iChan,:,iColor)),'color',EPscanCont.waveColors{iColor},'LineStyle',EPscanCont.lineStyles{EPscanCont.displayChans(iChan)});
                hold off
            end;
        end;
    end;
    axis([1 EPscanCont.displayLength -length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing])
    if ~isempty(EPscanCont.plotPoint)
        EPscanCont.handles.plotLine=line([EPscanCont.plotPoint EPscanCont.plotPoint],[-length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing],'Color','green','LineWidth',3);
    end;
else
    EPscanCont.handles.TFTgraph = imagesc([1 EPscanCont.displayLength],[1 EPscanCont.spacing*length(EPscanCont.displayChans)],squeeze(EPscanCont.currentData(:,:,1)));
    axis([1 EPscanCont.displayLength 1 EPscanCont.spacing*length(EPscanCont.displayChans)])
    if ~isempty(EPscanCont.plotPoint)
        EPscanCont.handles.plotLine=line([EPscanCont.plotPoint EPscanCont.plotPoint],[1 EPscanCont.spacing*length(EPscanCont.displayChans)],'Color','green','LineWidth',3);
    end;
end;

set(EPscanCont.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
if (EPmain.preferences.view.positive ==2) || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow == length(EPscanCont.freqList)))
    set(EPscanCont.handles.graphPlot,'YDir','reverse')
else
    set(EPscanCont.handles.graphPlot,'YDir','normal')
end;

%redraw the events
if ~isempty(EPscanCont.EPdata(1).events{1})
    eventList=find(ismember(round([EPscanCont.EPdata(1).events{1}.sample]),EPscanCont.displayPoints));
    if xChange
        EPscanCont.handles.events=gobjects(length(eventList),1);
    end;
    labelPos=EPscanCont.graphCoords(3)/EPscanCont.displayLength;
    for iEvent=1:length(eventList)
        theEventValue=EPscanCont.EPdata(1).events{1}(eventList(iEvent)).value;
        if ~isempty(theEventValue)
            if xChange
                if EPscanCont.eventLabels
                    if strcmp(EPscanCont.EPdata(1).events{1}(eventList(iEvent)).type,'eye-tracker')
                        yPos=80;
                    else
                        yPos=100;
                    end;
                    EPscanCont.handles.events(iEvent)=uicontrol('Style','text','HorizontalAlignment','left','String', theEventValue,'FontSize',EPmain.fontsize,...
                        'Position',[EPscanCont.graphCoords(1)+labelPos*(EPscanCont.EPdata(1).events{1}(eventList(iEvent)).sample-EPscanCont.displayPoints(1)+1) EPmain.scrsz(4)-yPos 60 20]);
                end;
            end;
            if strcmp(EPscanCont.EPdata(1).events{1}(eventList(iEvent)).value,'boundary')
                theColor='red';
            else
                theColor='black';
            end;
            if strcmp(EPscanCont.dataType,'VLT') || (strcmp(EPscanCont.dataType,'TFT') && (EPscanCont.freqShow ~= length(EPscanCont.freqList)))
                EPscanCont.handles.eventLines(iEvent)=line([EPscanCont.EPdata(1).events{1}(eventList(iEvent)).sample-EPscanCont.displayPoints(1)+1 EPscanCont.EPdata(1).events{1}(eventList(iEvent)).sample-EPscanCont.displayPoints(1)+1],[EPscanCont.displayLength -length(EPscanCont.displayChans)*EPscanCont.spacing],'Color',theColor,'LineWidth',1);
            else
                EPscanCont.handles.eventLines(iEvent)=line([EPscanCont.EPdata(1).events{1}(eventList(iEvent)).sample-EPscanCont.displayPoints(1)+1 EPscanCont.EPdata(1).events{1}(eventList(iEvent)).sample-EPscanCont.displayPoints(1)+1],[1 EPscanCont.spacing*length(EPscanCont.displayChans)],'Color',theColor,'LineWidth',1);
            end;
        end;
    end;
end;

%redraw the bad data markings

if EPscanCont.editMode
    badDataMarking
end;

set(EPscanCont.handles.displayStart,'String', ep_ms2time((EPscanCont.displayPoints(1)-1)*(1000/EPscanCont.EPdata(1).Fs)));
set(EPscanCont.handles.displayEnd,'String', ep_ms2time(EPscanCont.displayPoints(end)*(1000/EPscanCont.EPdata(1).Fs)));

drawPlots

updateFreqGraph

drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D3headPlot(src,eventdata,theCell)
global EPscanCont EPmain EPdataset

disp('Using EEGlab function headplot to perform 3D head display.');

theData=squeeze(EPscanCont.EPdata(theCell).data(find(ismember(EPscanCont.chanIndex{theCell},EPscanCont.hasLoc)),EPscanCont.displayPoints(EPscanCont.plotPoint))-EPscanCont.baseline(EPscanCont.hasLoc,theCell));

theMin=inf;
theMax=-inf;
for iColor=1:EPscanCont.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        theMin=min(theMin,min(EPscanCont.EPdata(iColor).data(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.hasLoc)),EPscanCont.displayPoints(EPscanCont.plotPoint))-EPscanCont.baseline(EPscanCont.hasLoc,theCell)));
        theMax=max(theMax,max(EPscanCont.EPdata(iColor).data(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.hasLoc)),EPscanCont.displayPoints(EPscanCont.plotPoint))-EPscanCont.baseline(EPscanCont.hasLoc,theCell)));
    end;
end;

%check to see if spline file needs to be generated
if ~isempty(EPscanCont.EPdata(theCell).ced) && ~any(strcmp(EPscanCont.EPdata(theCell).ced,{'none','internal'}))
    [pathstr, name, ext] = fileparts(EPscanCont.EPdata(theCell).ced);
    CEDloc=which(EPscanCont.EPdata(theCell).ced);
else
    name=EPscanCont.EPdata(theCell).dataName;
    CEDloc=[pwd filesep 'temp'];
end;
if isempty(which([name '.spl']))
    [pathstr2, name2, ext2] = fileparts(CEDloc);
    if max(abs([EPscanCont.EPdata(theCell).eloc(EPscanCont.hasLoc).X])) <= 1
        MNItransform=[ 0 -15 0 0.08 0 -1.571 102 93 100 ]; %assume eloc coordinates are from a .elp file
    else
        MNItransform=[ 0 -15 4 0.05 0 -1.571 10.2 12 12.2 ]; %assume eloc coordinates are from a .sfp file
    end;
    headplot('setup', EPscanCont.EPdata(1).eloc(EPscanCont.hasLoc), [pathstr2 filesep name '.spl'],'transform',MNItransform); %save spline file in same place as the ced file.
end;

figure
[hdaxis cbaraxis] = headplot(theData,which([name '.spl']),'cbar',0,'maplimits',[theMin theMax],'labels',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eventLabels(src,eventdata)
%toggles the event labels on and off
global EPscanCont EPmain

EPscanCont.eventLabels=get(EPscanCont.handles.eventLabels,'Value');

if ~EPscanCont.eventLabels
    for iEvent=1:length(EPscanCont.handles.events)
        if isgraphics(EPscanCont.handles.events(iEvent))
            delete(EPscanCont.handles.events(iEvent));
        end;
    end;
    EPscanCont.handles.events=[];
else
    if ~isempty(EPscanCont.EPdata(1).events)
        eventList=find(ismember(round([EPscanCont.EPdata(1).events{1}.sample]),EPscanCont.displayPoints));
        EPscanCont.handles.events=gobjects(length(eventList),1);
        labelPos=EPscanCont.graphCoords(3)/EPscanCont.displayLength;
        for iEvent=1:length(eventList)
            if ~isempty(EPscanCont.EPdata(1).events{1}(eventList(iEvent)).value)
                if strcmp(EPscanCont.EPdata(1).events{1}(eventList(iEvent)).type,'eye-tracker')
                    yPos=80;
                else
                    yPos=100;
                end;
                EPscanCont.handles.events(iEvent)=uicontrol('Style','text','HorizontalAlignment','left','String', EPscanCont.EPdata(1).events{1}(eventList(iEvent)).value,'FontSize',EPmain.fontsize,...
                    'Position',[EPscanCont.graphCoords(1)+labelPos*(EPscanCont.EPdata(1).events{1}(eventList(iEvent)).sample-EPscanCont.displayPoints(1)+1) EPmain.scrsz(4)-yPos 60 20]);
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editMode(src,eventdata)
%toggles the red edit markings on and off
global EPscanCont

EPscanCont.editMode=get(EPscanCont.handles.editMode,'Value');

if ~EPscanCont.editMode
    for iChan=1:size(EPscanCont.handles.trialBadChans,1)
        for iEpoch=1:size(EPscanCont.handles.trialBadChans,2)
            if isgraphics(EPscanCont.handles.trialBadChans(iChan,iEpoch))
                delete(EPscanCont.handles.trialBadChans(iChan,iEpoch));
            end;
        end;
    end;
    EPscanCont.handles.trialBadChans=[];
    for iEpoch=1:length(EPscanCont.handles.badTrials)
        if isgraphics(EPscanCont.handles.badTrials(iEpoch))
            delete(EPscanCont.handles.badTrials(iEpoch));
        end;
    end;
    EPscanCont.handles.badTrials=[];
else
    badDataMarking
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function globalBadChans(src,eventdata,theChan)
%toggles the global bad status of a channel
global EPscanCont EPmain

if ~EPscanCont.editMode
    set(EPscanCont.globalBadChans(theChan),1-get(EPscanCont.handles.chanNames(theChan),'Value'));
    drawnow
    return;
end;

EPscanCont.globalBadChans(theChan)=get(EPscanCont.handles.chanNames(theChan),'Value');

if ~EPscanCont.globalBadChans(theChan)
    EPscanCont.globalBadChans(theChan)=0;
    EPscanCont.lineStyles{theChan}='-';
    %delete(EPscanCont.handles.globalBadChans(theChan));
else
    EPscanCont.globalBadChans(theChan)=1;
    EPscanCont.lineStyles{theChan}=':';
    axes(EPscanCont.handles.graphPlot);
%     EPscanCont.handles.globalBadChans(theChan)=patch([1 1 EPscanCont.displayLength EPscanCont.displayLength],[-EPscanCont.spacing*(theChan-.25) -EPscanCont.spacing*(theChan-1.25) -EPscanCont.spacing*(theChan-1.25) -EPscanCont.spacing*(theChan-.25)],'red','EdgeColor','red');
%     alpha(EPscanCont.handles.globalBadChans(theChan),.5);
end;

EPscanCont.edited=1;

updateDisplay(0,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateFreqGraph
%update the contents of the spectral domain graph plots

global EPscanCont EPmain EPdataset

if ~isempty(EPscanCont.plotPoint)
    theEpoch=ceil((EPscanCont.displayPoints(EPscanCont.plotPoint))/EPscanCont.EPdata(1).Fs);
    if (theEpoch ~= EPscanCont.plotLineEpoch) || (theEpoch==(size(EPscanCont.trialBadChans,2)))
        if (theEpoch==(size(EPscanCont.trialBadChans,2)+1)) && rem(EPscanCont.numPoints,EPscanCont.EPdata(1).Fs)
            EPscanCont.freqData=ones(length(EPscanCont.displayChans),length(EPscanCont.freqBins),4);
            for iChan=1:length(EPscanCont.displayChans)
                EPscanCont.freqData(iChan,:,:)=EPscanCont.freqData(iChan,:,:)-EPscanCont.spacing*(iChan-1);
            end;
        else
            EPscanCont.plotLineEpoch=theEpoch;
            cfg=[];
            cfg.method='mtmfft';
            cfg.output='pow';
            cfg.taper='hanning';
            deltaT=1; %time window in seconds
            deltaF=1/deltaT; %frequency resolution
            cfg.foi=[deltaF:deltaF:EPscanCont.EPdata(1).Fs/2];
            cfg.pad='maxperlen';
            cfg.feedback='no';
            numPoints=EPscanCont.EPdata(1).Fs;
            numChans=size(EPscanCont.currentData,1);
            
            FTdata.hdr.Fs=EPscanCont.EPdata(1).Fs;
            FTdata.hdr.nChans=numChans;
            FTdata.hdr.label=EPscanCont.chanNames(EPscanCont.displayChans);
            FTdata.hdr.nTrials=1;
            FTdata.hdr.nSamplesPre=EPscanCont.EPdata(1).baseline/(1000/EPscanCont.EPdata(1).Fs);
            FTdata.hdr.nSamples=numPoints;
            FTdata.label=EPscanCont.chanNames(EPscanCont.displayChans);
            FTdata.time=cell(1,1);
            FTdata.trial=cell(1,1);
            FTdata.sampleinfo=zeros(1,2);
            FTdata.time{1,1}=([1:EPscanCont.sampleSize:1000]')/1000; %time in seconds
            
            FTdata.sampleinfo(1,1)=((1-1)*numPoints)+1;
            FTdata.sampleinfo(1,2)=1*numPoints;
            FTdata.dimord='{rpt}_chan_time';
            FTdata.fsample=EPscanCont.EPdata(1).Fs;
            
            EPscanCont.freqData=zeros(length(EPscanCont.displayChans),length(EPscanCont.freqBins),4);
            for iColor=1:EPscanCont.numColors
                if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                    if strcmp(EPscanCont.dataType,'VLT')
                        FTdata.trial{1,1}=zeros(numChans,numPoints);
                        FTdata.trial{1,1}(find(ismember(EPscanCont.displayChans,EPscanCont.chanIndex{iColor})),:)=EPscanCont.EPdata(iColor).data(find(ismember(EPscanCont.chanIndex{iColor},EPscanCont.displayChans)),((theEpoch-1)*EPscanCont.EPdata(1).Fs)+1:theEpoch*EPscanCont.EPdata(1).Fs);
                        evalc('[freq] = ft_freqanalysis(cfg, FTdata);'); %pow spectrum
                        EPscanCont.freqData(:,:,iColor)= log10(freq.powspctrm(:,EPscanCont.freqBins))*10; %dB scale
                    else
                        FTdata.hdr.nChans=1;
                        for iChan=1:length(EPscanCont.displayChans)
                            theChan=EPscanCont.displayChans(iChan);
                            FTdata.hdr.label=EPscanCont.chanNames(theChan);
                            FTdata.label=EPscanCont.chanNames(theChan);
                            cfg.channel=EPscanCont.chanNames(theChan);
                            if any(strcmp(EPscanCont.chanNames(theChan),{'Hsaccade','Vsaccade','blink','SacPot'}))
                                freqChan=zeros(length(EPscanCont.freqBins),1);
                            else
                                freqChan=squeeze(abs(EPscanCont.EPdata(iColor).data(find(EPscanCont.chanIndex{iColor}==theChan),((theEpoch-1)*EPscanCont.EPdata(1).Fs)+1:theEpoch*EPscanCont.EPdata(1).Fs,:,:,:,EPscanCont.freqBins))); %convert complex number to real number
                                freqChan=freqChan/mean(diff(EPscanCont.EPdata(1).freqNames)); %convert to spectral density
                                freqChan=freqChan.^2; %convert amplitude to power
                                freqChan=log10(abs(freqChan))*10; %convert to dB log scaling
                                freqChan(isinf(freqChan))=-flintmax;%log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                                freqChan=mean(freqChan);
                            end;
                            EPscanCont.freqData(iChan,:,iColor)= freqChan;
                        end;
                    end;
                    for iChan=1:length(EPscanCont.displayChans)
                        EPscanCont.freqData(iChan,:,iColor)=EPscanCont.freqData(iChan,:,iColor)-EPscanCont.spacing*(iChan-1);
                    end;
                end;
            end;
        end;
        axes(EPscanCont.handles.freqPlot);
        cla(EPscanCont.handles.freqPlot);
        for iColor=1:EPscanCont.numColors
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                hold on
                EPscanCont.handles.freqChans{iColor}=plot(EPscanCont.freqBins,squeeze(EPscanCont.freqData(:,:,iColor))','color',EPscanCont.waveColors{iColor});
                hold off
            end;
        end;
        for iLine=1:length(EPscanCont.freqLines)
            EPscanCont.handles.freqLines(iLine)=line([EPscanCont.freqLines(iLine) EPscanCont.freqLines(iLine)],[-length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing],'Color','black','LineWidth',1);
        end;
        axis([1 length(EPscanCont.freqBins) -length(EPscanCont.displayChans)*EPscanCont.spacing EPscanCont.spacing])
    end;
else
    cla(EPscanCont.handles.freqPlot);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeEvent(src,eventdata)
%change redisplay event controls.

global EPscanCont

EPscanCont.firstEvent=get(EPscanCont.handles.firstEvent,'Value');
EPscanCont.lastEvent=get(EPscanCont.handles.lastEvent,'Value');

if EPscanCont.lastEvent <= EPscanCont.firstEvent
    set(EPscanCont.handles.redisplay,'enable','off');
else
    set(EPscanCont.handles.redisplay,'enable','on');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function redisplay(src,eventdata)
%change display to show the period between the two chosen events.

global EPscanCont

if EPscanCont.lastEvent <= EPscanCont.firstEvent
    return
end;

oldDisplayPoints=EPscanCont.displayPoints;
theFirstEvent=EPscanCont.alphabetEventOrder(EPscanCont.firstEvent);
theLastEvent=EPscanCont.alphabetEventOrder(EPscanCont.lastEvent);
EPscanCont.displayPoints=[EPscanCont.EPdata(1).events{1}(theFirstEvent).sample:EPscanCont.EPdata(1).events{1}(theLastEvent).sample];
EPscanCont.displayLength=length(EPscanCont.displayPoints);
if ismember(oldDisplayPoints(EPscanCont.plotPoint),EPscanCont.displayPoints)
    EPscanCont.plotPoint=find(oldDisplayPoints(EPscanCont.plotPoint)==EPscanCont.displayPoints);
else
    EPscanCont.plotPoint=[];
end;

set(EPscanCont.handles.xSlider,'SliderStep',[EPscanCont.displayLength/EPscanCont.numPoints max(0.2,length(EPscanCont.displayLength)/EPscanCont.numPoints)]);
EPscanCont.xGraphPos=(EPscanCont.displayPoints(1)-1)/EPscanCont.numPoints;
set(EPscanCont.handles.xSlider,'Value',EPscanCont.xGraphPos);

updateDisplay(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeCurrentEvent(src,eventdata)
%change current event based on menu.

global EPscanCont

EPscanCont.currentEvent=get(EPscanCont.handles.currentEvent,'Value');
theCurrentEvent=EPscanCont.sampleEventOrder(EPscanCont.currentEvent);
EPscanCont.currentkey=length(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys);

if ~ismember(round(EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample),EPscanCont.displayPoints)
    EPscanCont.plotPoint=[];
else
    EPscanCont.plotPoint=find(round(EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample) == EPscanCont.displayPoints);
end;

updateDisplay(0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function recenter(src,eventdata)
%change display to center on the current event.

global EPscanCont

theCurrentEvent=EPscanCont.sampleEventOrder(EPscanCont.currentEvent);
EPscanCont.displayPoints=[round(EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample)-ceil(EPscanCont.displayLength/2)+1:round(EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample)+ceil(EPscanCont.displayLength/2)];
if EPscanCont.displayPoints(1) < 1
    EPscanCont.displayPoints=EPscanCont.displayPoints-EPscanCont.displayPoints(1)+1;
end;
if EPscanCont.displayPoints(end) > EPscanCont.numPoints
    EPscanCont.displayPoints=EPscanCont.displayPoints-(EPscanCont.displayPoints(end)-EPscanCont.numPoints);
end;
EPscanCont.displayLength=length(EPscanCont.displayPoints);
EPscanCont.plotPoint=find(round(EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample)==EPscanCont.displayPoints);

set(EPscanCont.handles.xSlider,'SliderStep',[EPscanCont.displayLength/EPscanCont.numPoints max(0.2,length(EPscanCont.displayLength)/EPscanCont.numPoints)]);
EPscanCont.xGraphPos=(EPscanCont.displayPoints(1)-1)/EPscanCont.numPoints;
set(EPscanCont.handles.xSlider,'Value',EPscanCont.xGraphPos);

updateDisplay(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeFreq(src,eventdata)
%change frquency displayed for TFT data.

global EPscanCont

EPscanCont.freqShow=get(EPscanCont.handles.freqList,'Value');

if EPscanCont.freqShow == length(EPscanCont.freqList)
    EPscanCont.freqPos=1;
	EPscanCont.spacing=EPscanCont.lastBins-EPscanCont.startBins+1;
    if EPscanCont.FFTunits==1
        EPscanCont.spacing=EPscanCont.spacing*2; %complex numbers
    end;
else
    EPscanCont.freqPos=EPscanCont.freqShow;
    EPscanCont.spacing=100;
end;

updateDisplay(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeChan(src,eventdata)
%change channel used to color the eye position plots

global EPscanCont

EPscanCont.chanShow=get(EPscanCont.handles.chanList,'Value');

drawPlots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeEyeSettings(src,eventdata)
%change the settings for the eye position plots

global EPscanCont

if ~isempty(EPscanCont.EYEsize)
    theVar=str2num(get(EPscanCont.handles.topos.eyeCenterX,'String'));
    if ~isempty(theVar)
        EPscanCont.eyeCenterX=theVar;
    end;
    
    theVar=str2num(get(EPscanCont.handles.topos.eyeCenterY,'String'));
    if ~isempty(theVar)
        EPscanCont.eyeCenterY=theVar;
    end;
    
    theVar=str2num(get(EPscanCont.handles.topos.eyeScaleX,'String'));
    if ~isempty(theVar)
        EPscanCont.eyeScaleX=theVar; %minus means flip
    end;
    
    theVar=str2num(get(EPscanCont.handles.topos.eyeScaleY,'String'));
    if ~isempty(theVar)
        EPscanCont.eyeScaleY=theVar; %minus means flip
    end;
end;

if ~isempty(EPscanCont.SACCsize)
    theVar=str2num(get(EPscanCont.handles.topos.sacCenterX,'String'));
    if ~isempty(theVar)
        EPscanCont.sacCenterX=theVar;
    end;
    
    theVar=str2num(get(EPscanCont.handles.topos.sacCenterY,'String'));
    if ~isempty(theVar)
        EPscanCont.sacCenterY=theVar;
    end;
    
    theVar=str2num(get(EPscanCont.handles.topos.sacScaleX,'String'));
    if ~isempty(theVar) && (theVar > 0)
        EPscanCont.sacScaleX=theVar;
    end;
    
    theVar=str2num(get(EPscanCont.handles.topos.sacScaleY,'String'));
    if ~isempty(theVar) && (theVar > 0)
        EPscanCont.sacScaleY=theVar;
    end;
end;

EPscanCont.edited=1;

drawPlots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editEvent(src,eventdata)
%edit the information of the current event

global EPscanCont

EPscanCont.edited=1;
changedEventFlag=0; %need to update events list
changedKeyFlag=0;  %need to update stims list
if ~isempty(EPscanCont.currentEvent)
    theCurrentEvent=EPscanCont.sampleEventOrder(EPscanCont.currentEvent);
end;
switch src
    case EPscanCont.handles.addEvent
        EPscanCont.EPdata(1).events{1}(end+1)=struct('type',' ','sample',EPscanCont.displayPoints(EPscanCont.plotPoint),'value',' ','duration',' ','keys',struct('code','','data','','datatype','','description',''));
        EPscanCont.currentEvent=find(EPscanCont.sampleEventOrder==length(EPscanCont.sampleEventOrder));
        theCurrentEvent=length(EPscanCont.sampleEventOrder);
        EPscanCont.currentkey=1;
        changedEventFlag=1;
    case EPscanCont.handles.minusEvent
        keyData=EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys;
        for iKey=1:length(keyData)
            if strcmp(keyData(iKey).code,'stim')
                changedKeyFlag=1;
            end;
        end;
        EPscanCont.EPdata(1).events{1}(theCurrentEvent)=[];
        changedEventFlag=1;
    case EPscanCont.handles.eventType
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).type=get(EPscanCont.handles.eventType,'String');
    case EPscanCont.handles.eventValue
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).value=get(EPscanCont.handles.eventValue,'String');
    case EPscanCont.handles.eventSample
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample=get(EPscanCont.handles.eventSample,'String');
    case EPscanCont.handles.eventDuration
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).duration=get(EPscanCont.handles.eventDuration,'String');
    case EPscanCont.handles.eventKeys
        EPscanCont.currentkey=get(EPscanCont.handles.eventKeys,'Value');
    case EPscanCont.handles.addKey
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(end+1)=struct('code','','data','','datatype','','description','');
        EPscanCont.currentkey=EPscanCont.currentkey+1;
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).datatype='char';
    case EPscanCont.handles.minusKey
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey)=[];
        if isempty(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys)
            EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(1)=struct('code','','data','','datatype','','description','');
        elseif EPscanCont.currentkey > length(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys)
            EPscanCont.currentkey=EPscanCont.currentkey-1;
        end;    
        changedKeyFlag=1;        
    case EPscanCont.handles.keyCode
        oldKey=EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).code;
        newKey=get(EPscanCont.handles.keyCode,'String');
        if strcmp(oldKey,'stim') || strcmp(newKey,'stim')
            changedKeyFlag=1;
        end;
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).code=newKey;
    case EPscanCont.handles.keyData
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).data=get(EPscanCont.handles.keyData,'String');
        if strcmp(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).code,'stim')
            changedKeyFlag=1;
        end;
    case EPscanCont.handles.keyType
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).datatype=get(EPscanCont.handles.keyType,'String');
    case EPscanCont.handles.keyDescrip
        EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).description=get(EPscanCont.handles.keyDescrip,'String');
end;

if changedEventFlag
    
    set(EPscanCont.handles.addEvent,'enable','off');
    set(EPscanCont.handles.minusEvent,'enable','off');
    set(EPscanCont.handles.eventType,'enable','off');
    set(EPscanCont.handles.eventValue,'enable','off');
    set(EPscanCont.handles.eventSample,'enable','off');
    set(EPscanCont.handles.eventDuration,'enable','off');
    set(EPscanCont.handles.eventKeys,'enable','off');
    set(EPscanCont.handles.addKey,'enable','off');
    set(EPscanCont.handles.minusKey,'enable','off');
    set(EPscanCont.handles.keyCode,'enable','off');
    set(EPscanCont.handles.keyData,'enable','off');
    set(EPscanCont.handles.keyType,'enable','off');
    set(EPscanCont.handles.keyDescrip,'enable','off');
    drawnow

    sampleList=zeros(length(EPscanCont.EPdata(1).events{1}),1);
    EPscanCont.eventList=cell(length(EPscanCont.EPdata(1).events{1}),1);
    for iEvent=1:length(EPscanCont.EPdata(1).events{1})
        theEventValue=EPscanCont.EPdata(1).events{1}(iEvent).value;
        theEventSample=EPscanCont.EPdata(1).events{1}(iEvent).sample;
        sampleList(iEvent,1)=theEventSample;
        if ~isempty(theEventValue)
            eventNum=sprintf('%04d',(length(find([EPscanCont.EPdata(1).events{1}(strcmp(theEventValue,{EPscanCont.EPdata(1).events{1}.value})).sample]<=theEventSample))));
            EPscanCont.eventList{iEvent}=[theEventValue '(' eventNum ')' keyStim];
        else
            theEventType=EPscanCont.EPdata(1).events{1}(iEvent).type;
            if ~isempty(theEventType)
                eventNum=sprintf('%04d',(length(find([EPscanCont.EPdata(1).events{1}(strcmp(theEventType,{EPscanCont.EPdata(1).events{1}.value})).sample]<=theEventSample))));
                EPscanCont.eventList{iEvent}=[theEventType '(' eventNum ')' keyStim];
            else
                EPscanCont.eventList{iEvent}='null_event';
            end;
        end;
    end;
    
    [B EPscanCont.alphabetEventOrder]=sort(EPscanCont.eventList);
    [B EPscanCont.sampleEventOrder]=sort(sampleList);
    if isempty(EPscanCont.eventList)
        EPscanCont.starredSampleEventList='none';
        EPscanCont.starredEventList='none';
    else
        EPscanCont.starredEventList=cell(length(EPscanCont.eventList),1);
        EPscanCont.starredSampleEventList=cell(length(EPscanCont.eventList),1);
        for iEvent=1:length(EPscanCont.eventList)
            if ismember(round(EPscanCont.EPdata(1).events{1}(iEvent).sample),EPscanCont.displayPoints)
                EPscanCont.starredEventList{find(iEvent==EPscanCont.alphabetEventOrder)}=['*' EPscanCont.eventList{iEvent}];
                sampleOrderNumber=find(iEvent==EPscanCont.sampleEventOrder);
                EPscanCont.starredSampleEventList{sampleOrderNumber}=['*' EPscanCont.eventList{iEvent}];
            else
                EPscanCont.starredEventList{find(iEvent==EPscanCont.alphabetEventOrder)}=EPscanCont.eventList{iEvent};
                EPscanCont.starredSampleEventList{find(iEvent==EPscanCont.sampleEventOrder)}=EPscanCont.eventList{iEvent};
            end;
        end;
    end;
    set(EPscanCont.handles.firstEvent,'String',EPscanCont.starredEventList);
    set(EPscanCont.handles.lastEvent,'String',EPscanCont.starredEventList);
    set(EPscanCont.handles.currentEvent,'String',EPscanCont.starredSampleEventList);
    
end;

if changedKeyFlag
    for iEvent=1:length(EPscanCont.EPdata(1).events{1})
        keyData=EPscanCont.EPdata(1).events{1}(iEvent).keys;
        keyStim=[];
        for iKey=1:length(keyData)
            if strcmp(keyData(iKey).code,'stim')
                keyStim=keyData(iKey).data;
                if any(strcmp(keyStim,{EPscanCont.EPdata(1).stims.name}))
                    EPscanCont.stimList(end+1).sample=EPscanCont.EPdata(1).events{1}(iEvent).sample;
                    EPscanCont.stimList(end).stim=find(strcmp(keyStim,{EPscanCont.EPdata(1).stims.name}));
                end;
            end;
        end;
    end;
end;

if ~isempty(EPscanCont.currentEvent)
    set(EPscanCont.handles.eventType,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).type);
    set(EPscanCont.handles.eventValue,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).value);
    set(EPscanCont.handles.eventSample,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample);
    set(EPscanCont.handles.eventDuration,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).duration);
    
    theString=cell(length(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys),1);
    for iKey=1:length(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys)
        theString{iKey}=num2str(iKey);
    end;
    set(EPscanCont.handles.eventKeys,'String',theString);
    set(EPscanCont.handles.eventKeys,'Value',EPscanCont.currentkey);
    
    set(EPscanCont.handles.keyCode,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).code);
    set(EPscanCont.handles.keyData,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).data);
    set(EPscanCont.handles.keyType,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).datatype);
    set(EPscanCont.handles.keyDescrip,'String',EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys(EPscanCont.currentkey).description);
end;
if changedEventFlag
    updateDisplay(1,0);
    set(EPscanCont.handles.addEvent,'enable','on');
    set(EPscanCont.handles.minusEvent,'enable','on');
    set(EPscanCont.handles.eventType,'enable','on');
    set(EPscanCont.handles.eventValue,'enable','on');
    set(EPscanCont.handles.eventSample,'enable','on');
    set(EPscanCont.handles.eventDuration,'enable','on');
    set(EPscanCont.handles.eventKeys,'enable','on');
    set(EPscanCont.handles.addKey,'enable','on');
    set(EPscanCont.handles.minusKey,'enable','on');
    set(EPscanCont.handles.keyCode,'enable','on');
    set(EPscanCont.handles.keyData,'enable','on');
    set(EPscanCont.handles.keyType,'enable','on');
    set(EPscanCont.handles.keyDescrip,'enable','on');
end;
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shiftEvent(src,eventdata)
%shift the current event left or right

global EPscanCont

eventSamples=round([EPscanCont.EPdata(1).events{1}(EPscanCont.sampleEventOrder).sample]');

if (src == EPscanCont.handles.backEvent) && (EPscanCont.currentEvent > 1)
    EPscanCont.currentEvent=max(find(eventSamples<EPscanCont.displayPoints(EPscanCont.plotPoint)));
elseif (src == EPscanCont.handles.nextEvent) && (EPscanCont.currentEvent < length(EPscanCont.eventList))
    EPscanCont.currentEvent=min(find(eventSamples>EPscanCont.displayPoints(EPscanCont.plotPoint)));
end;

set(EPscanCont.handles.currentEvent,'Value',EPscanCont.currentEvent);

theCurrentEvent=EPscanCont.sampleEventOrder(EPscanCont.currentEvent);
EPscanCont.currentkey=length(EPscanCont.EPdata(1).events{1}(theCurrentEvent).keys);

if ~ismember(round(EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample),EPscanCont.displayPoints)
    EPscanCont.plotPoint=[];
else
    EPscanCont.plotPoint=find(round(EPscanCont.EPdata(1).events{1}(theCurrentEvent).sample) == EPscanCont.displayPoints);
end;

updateDisplay(0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function e1e2(src,eventdata)
%change e1 menu or e2 menu to the current event

global EPscanCont

if (src == EPscanCont.handles.e1)
    set(EPscanCont.handles.firstEvent,'Value',find(EPscanCont.alphabetEventOrder == EPscanCont.sampleEventOrder(EPscanCont.currentEvent)));
elseif (src == EPscanCont.handles.e2)
    set(EPscanCont.handles.lastEvent,'Value',find(EPscanCont.alphabetEventOrder == EPscanCont.sampleEventOrder(EPscanCont.currentEvent)));
end;

changeEvent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function badDataMarking
%add bad data markings

global EPscanCont

axes(EPscanCont.handles.graphPlot);
for iChan=1:length(EPscanCont.displayChans)
    theChan=EPscanCont.displayChans(iChan);
    if any(EPscanCont.trialBadChans(find(EPscanCont.chanIndex{1}==theChan),:)==-1)
        for iEpoch=floor((EPscanCont.displayPoints(1)-1)/EPscanCont.EPdata(1).Fs)+1:floor((EPscanCont.displayPoints(end)-1)/EPscanCont.EPdata(1).Fs)+1
            graphEpoch=iEpoch-(floor((EPscanCont.displayPoints(1)-1)/EPscanCont.EPdata(1).Fs)+1)+1;
            excessEpoch=0;
            theEpoch=iEpoch;
            if iEpoch>size(EPscanCont.trialBadChans,2)
                continue
            elseif (iEpoch==size(EPscanCont.trialBadChans,2)) && rem(EPscanCont.numPoints,EPscanCont.EPdata(1).Fs)
                excessEpoch=1; %extra time points tacked onto last epoch
            end;
            if EPscanCont.trialBadChans(theChan,theEpoch)==-1
                if excessEpoch
                    startEpoch=((graphEpoch-1)*EPscanCont.EPdata(1).Fs)+1-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
                    endEpoch=EPscanCont.displayLength;
                else
                    startEpoch=((graphEpoch-1)*EPscanCont.EPdata(1).Fs)+1-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
                    endEpoch=graphEpoch*EPscanCont.EPdata(1).Fs-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
                end;
                EPscanCont.handles.trialBadChans(theChan,theEpoch)=patch([startEpoch startEpoch endEpoch endEpoch],...
                    [-EPscanCont.spacing*(iChan-.25) -EPscanCont.spacing*(iChan-1.25) -EPscanCont.spacing*(iChan-1.25) -EPscanCont.spacing*(iChan-.25)],'red','EdgeColor','red');
                alpha(EPscanCont.handles.trialBadChans(theChan,theEpoch),.5);
            end;
        end;
    end;
end;
for iEpoch=floor((EPscanCont.displayPoints(1)-1)/EPscanCont.EPdata(1).Fs)+1:floor((EPscanCont.displayPoints(end)-1)/EPscanCont.EPdata(1).Fs)+1
    graphEpoch=iEpoch-(floor((EPscanCont.displayPoints(1)-1)/EPscanCont.EPdata(1).Fs)+1)+1;
    theEpoch=iEpoch;
    excessEpoch=0;
    if iEpoch>length(EPscanCont.badTrials)
        continue
    elseif (iEpoch==size(EPscanCont.trialBadChans,2)) && rem(EPscanCont.numPoints,EPscanCont.EPdata(1).Fs)
        excessEpoch=1; %extra time points tacked onto last epoch
    end;
    if EPscanCont.badTrials(theEpoch)
        if excessEpoch
            startEpoch=((graphEpoch-1)*EPscanCont.EPdata(1).Fs)+1-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
            endEpoch=EPscanCont.displayLength;
        else
            startEpoch=((graphEpoch-1)*EPscanCont.EPdata(1).Fs)+1-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
            endEpoch=graphEpoch*EPscanCont.EPdata(1).Fs-rem(EPscanCont.displayPoints(1)-1,EPscanCont.EPdata(1).Fs);
        end;
        EPscanCont.handles.badTrials(theEpoch)=patch([startEpoch startEpoch endEpoch endEpoch],...
            [-EPscanCont.spacing*length(EPscanCont.displayChans) EPscanCont.spacing*length(EPscanCont.displayChans) EPscanCont.spacing*length(EPscanCont.displayChans) -EPscanCont.spacing*length(EPscanCont.displayChans)],'red','EdgeColor','red');
        alpha(EPscanCont.handles.badTrials(theEpoch),.5);
    end;
end;

