function ep_template
% ep_template
% Manual creation of blink and saccade templates.
%
%Input:
%  EPdataset (via global variable)
%
%Output: (via save file)
%   EPblink: structured file with information about blink template
%      .template   : the averaged voltage values for the template
%      .eloc       : the electrode coordinates
%      .num        : the number of blinks averaged into the template
%   EPsaccade: structured file with information about saccade template
%      .eloc       : the electrode coordinates (EEG chans)
%      .hSaccade.template   : the averaged voltage values for the template (EEG chans)
%      .hSaccade.num        : the number of saccades averaged into the template
%      .vSaccade.template   : the averaged voltage values for the template (EEG chans)
%      .vSaccade.num        : the number of saccades averaged into the template
%      .sacPot.template   : the averaged voltage values for the template (EEG chans)
%      .sacPot.num        : the number of saccade potentials averaged into the template

%History
%  by Joseph Dien (7/13/09)
%  jdien07@mac.com
%
%  bugfix 8/27/09 JD
%  Location of windows appearing partly off screen on some systems.
%
%  modified 8/30/09 JD
%  Added support for continuous data.
%
%  bugfix 9/5/09 JD
%  Crash when loading or saving blink template with spaces in the path name.
%  Maximum number of trials/segments not updating when switching between datasets.
%
%  bugfix 12/3/09 JD
%  Fixed crash when loading in blink template for a file or pathname with a space in it.
%
%  bugfix 12/8/09 JD
%  Bad channels incorrectly interpolated when adding a new blink to the template.
%
%  modified 10/28/10 JD
%  Added support for saccade template.  Added undo and clear buttons.
%
%  bugfix 6/1/12 JD
%  Fixed crash when there is an average file in the working set with multiple subjects.
%
%  bugfix 9/8/12 JD
%  Fixed blink templates being saved as text files.
%
%  bugfix 9/15/12 JD
%  Fixed blink and saccade waveforms not being displayed correctly.
% 
%  modified 10/10/13 JD
%  Supports presence of non-EEG channels.
%
%  bugfix 11/1/13 JD
%  Fixes font sizes on Windows.
%
%  bugfix 12/27/13 JD
%  Fixed crashes when there is a non EEG channel present.
%  Fixed crash when switching between datasets with different electrode montages or epoch length.
% 
%  modified 12/28/13 JD
%  Changed x-axis labels to milliseconds.
%
%  bugfix 1/12/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
%  bugfix 3/2/14 JD
%  Fix for not properly adding new blinks to the template for single trial data.
%
%  bufix 3/11/14 JD
%  Handles decimal sampling rates gracefully.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% bufix 5/21/14 JD
% Fixed crash when frequency or factor files present by excluding them entirely.
% Fixed template waveform plot not adjusting to new length when switching datasets.
% Fixed marker for blink and saccade not placed at correct latency after switching between datasets with the same number
% of samples but different baselines.
%
% bufix 7/28/14 JD
% Fixed crash when the dataset contains a boundary event.
%
% bufix 11/18/14 JD
% Fixed no time marker for blink and saccades when forming templates.
%
% modified 11/23/14 JD
% Added vertical saccade correction.
% Added click 2D topos to expand into 3D view.
% Added EOG channel controls.
%
% bufix 12/8/14 JD
% Fixed crash in template function when there are non-EEG channels present.
%
% modified 5/25/14 JD
% Set colormap to jet even for Matlab 2014b onwards.
%
% bufix 8/26/15 JD
% Fixed crash in template function when loading a template file and there
% are datasets with a different number of channels.
% Fixed number of trials in template not being reset when Clear button
% used.
%
% bufix 4/22/16 JD
% Fixed crash in template function when scrolling through continuous data and there is an event with an empty .value field.
%
% bufix 5/20/16 JD
% Fixed crash when clicking for 3D head of template head and there are bad channels present.
%
% bugfix and modified 9/8/16 JD
% Added saccade potential to set of templates.
% Can click and drag time point marker in waveform windows.
% Fixed crash in template function after performing a Clear Template for data with non-EEG channels.
%
% bugfix and modified 11/15/16 JD
% Added eye-tracker plot.
% Fixed waveform plots not updating when changing segment via popupmenu.
%
% bugfix and modified 5/4/17 JD
% Added saturation and badtrial preference settings to bad channel detection call.
% Various fixes to the behavior of the dataset and trial number controls.
% Made template undo specific to last change.
% Improved determination of horizontal saccade direction.
%
% modified 6/20/17 JD
% Improved dragging of red lines in Template Creation pane.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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

global EPmain EPdataset EPtemplate

EPtemplate.sessionFiles=[];
for theData=1:length(EPdataset.dataset)
    EPtemplate.EEGchans=find(strcmp('EEG',EPdataset.dataset(theData).chanTypes));
    if any(strcmp(EPdataset.dataset(theData).dataType,{'single_trial','continuous'})) && ~isempty(EPdataset.dataset(theData).eloc) && (length([EPdataset.dataset(theData).eloc(EPtemplate.EEGchans).theta]) == length(EPtemplate.EEGchans)) && (length([EPdataset.dataset(theData).eloc(EPtemplate.EEGchans).radius]) == length(EPtemplate.EEGchans)) && isempty(EPdataset.dataset(theData).freqNames) && isempty(EPdataset.dataset(theData).facNames)
        EPtemplate.sessionFiles=[EPtemplate.sessionFiles theData];
    end;
end;

if ~isempty(EPtemplate.sessionFiles)
    EPtemplate.dataset=length(EPtemplate.sessionFiles); %the dataset currently being examined
else
    msg{1}=['Error: There are no suitable single trial or continuous files with which to form a blink template.'];
    [msg]=ep_errorMsg(msg);
    return
end;

scrsz = EPmain.scrsz;

EPtemplate.windowHeight=700;

templateFigure=findobj('Name', 'Eye Artifact Templates Creation');

set(EPmain.handles.preprocess.template,'enable','off');
set(EPmain.handles.preprocess.preprocess,'enable','off');
set(EPmain.handles.preprocess.done,'enable','off');

if ~isempty(templateFigure)
    close(templateFigure)
end;
EPtemplate.handles.window = figure('Name', 'Eye Artifact Templates Creation', 'NumberTitle', 'off', 'Position',[201 scrsz(4)-EPtemplate.windowHeight 700 EPtemplate.windowHeight]);
colormap jet;

EPtemplate.handles.load = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
    'Position', [20 EPtemplate.windowHeight-100 60 30], 'Callback', @loadTemplate);

EPtemplate.handles.save = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
    'Position', [20 EPtemplate.windowHeight-130 60 30], 'Callback', @saveTemplate);

EPtemplate.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
    'Position', [20 EPtemplate.windowHeight-160 60 30], 'Callback', @clearTemplates);

EPtemplate.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Undo','FontSize',EPmain.fontsize,...
    'Position', [20 EPtemplate.windowHeight-190 60 30], 'Callback', @undo);

EPtemplate.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
    'Position', [20 EPtemplate.windowHeight-220 60 30], 'Callback', @done);


EPtemplate.trial=1; %the trial of the dataset currently being examined
EPtemplate.saved=1; %flag that the templates have just been saved

EPtemplate.EPdata=ep_loadEPdataset(EPdataset,EPtemplate.sessionFiles(EPtemplate.dataset));
if strcmp(EPtemplate.EPdata.dataType,'continuous')
    EPtemplate.segLength=min(ceil(EPtemplate.EPdata.Fs),length(EPtemplate.EPdata.timeNames));
    EPtemplate.maxSegs=floor(length(EPtemplate.EPdata.timeNames)/EPtemplate.segLength);
else
    EPtemplate.segLength=length(EPtemplate.EPdata.timeNames);
    EPtemplate.maxSegs=length(EPtemplate.EPdata.cellNames);
end;

EPtemplate.EEGchans=find(strcmp('EEG',EPtemplate.EPdata.chanTypes));

%remove bad channels
EPtemplate.criteria.neighbors=EPmain.preferences.preprocess.neighbors;
EPtemplate.criteria.badchan=EPmain.preferences.preprocess.badchan;
EPtemplate.criteria.saturation=[-EPmain.preferences.preprocess.saturation EPmain.preferences.preprocess.saturation];
EPtemplate.criteria.badtrials=EPmain.preferences.preprocess.badtrials;
EPtemplate.badChans=ep_detectBadChans(EPtemplate.EPdata, EPtemplate.criteria,1);

EPtemplate.minVolt=-400;
EPtemplate.maxVolt=400;
EPtemplate.minSacVolt=-400;
EPtemplate.maxSacVolt=400;
EPtemplate.gridSize=67;

EPtemplate.handles.events=[];
EPtemplate.screenCoords=[];
EPtemplate.clickGraph=0;

EPtemplate.lastChange=0;

initializeTemplate

EPtemplate.handles.dataset = uicontrol('Style','text',...
    'String','LUV','FontSize',EPmain.fontsize,'Position',[20 EPtemplate.windowHeight-250 30 20]);

EPtemplate.handles.LUVEOG = uicontrol('Style','popupmenu',...
    'String',EPtemplate.EPdata.chanNames(EPtemplate.EEGchans),'FontSize',EPmain.fontsize,...
    'Value', EPtemplate.VEOG(1),'Position',[50 EPtemplate.windowHeight-250 75 20],...
    'Callback', @changeEOG);

EPtemplate.handles.dataset = uicontrol('Style','text',...
    'String','RUV','FontSize',EPmain.fontsize,'Position',[20 EPtemplate.windowHeight-270 30 20]);

EPtemplate.handles.RUVEOG = uicontrol('Style','popupmenu',...
    'String',EPtemplate.EPdata.chanNames(EPtemplate.EEGchans),'FontSize',EPmain.fontsize,...
    'Value', EPtemplate.VEOG(2),'Position',[50 EPtemplate.windowHeight-270 75 20],...
    'Callback', @changeEOG);

EPtemplate.handles.dataset = uicontrol('Style','text',...
    'String','LLV','FontSize',EPmain.fontsize,'Position',[20 EPtemplate.windowHeight-290 30 20]);

EPtemplate.handles.LLVEOG = uicontrol('Style','popupmenu',...
    'String',EPtemplate.EPdata.chanNames(EPtemplate.EEGchans),'FontSize',EPmain.fontsize,...
    'Value', EPtemplate.VEOG(3),'Position',[50 EPtemplate.windowHeight-290 75 20],...
    'Callback', @changeEOG);

EPtemplate.handles.dataset = uicontrol('Style','text',...
    'String','RLV','FontSize',EPmain.fontsize,'Position',[20 EPtemplate.windowHeight-310 30 20]);

EPtemplate.handles.RLVEOG = uicontrol('Style','popupmenu',...
    'String',EPtemplate.EPdata.chanNames(EPtemplate.EEGchans),'FontSize',EPmain.fontsize,...
    'Value', EPtemplate.VEOG(4),'Position',[50 EPtemplate.windowHeight-310 75 20],...
    'Callback', @changeEOG);

EPtemplate.handles.dataset = uicontrol('Style','text',...
    'String','LH','FontSize',EPmain.fontsize,'Position',[20 EPtemplate.windowHeight-330 30 20]);

EPtemplate.handles.LHEOG = uicontrol('Style','popupmenu',...
    'String',EPtemplate.EPdata.chanNames(EPtemplate.EEGchans),'FontSize',EPmain.fontsize,...
    'Value', EPtemplate.HEOG(1),'Position',[50 EPtemplate.windowHeight-330 75 20],...
    'Callback', @changeEOG);

EPtemplate.handles.dataset = uicontrol('Style','text',...
    'String','RH','FontSize',EPmain.fontsize,'Position',[20 EPtemplate.windowHeight-350 30 20]);

EPtemplate.handles.RHEOG = uicontrol('Style','popupmenu',...
    'String',EPtemplate.EPdata.chanNames(EPtemplate.EEGchans),'FontSize',EPmain.fontsize,...
    'Value', EPtemplate.HEOG(2),'Position',[50 EPtemplate.windowHeight-350 75 20],...
    'Callback', @changeEOG);

EPtemplate.handles.leftTrial = uicontrol('Style', 'pushbutton', 'String', '<--','FontSize',EPmain.fontsize,...
    'Position', [120 EPtemplate.windowHeight-650 50 30], 'Callback', @leftTrial);
if EPtemplate.trial ==1
    set(EPtemplate.handles.leftTrial,'enable','off');
end;

EPtemplate.handles.rightTrial = uicontrol('Style', 'pushbutton', 'String', '-->','FontSize',EPmain.fontsize,...
    'Position', [330 EPtemplate.windowHeight-650 50 30], 'Callback', @rightTrial);
if EPtemplate.trial == EPtemplate.maxSegs
    set(EPtemplate.handles.rightTrial,'enable','off');
end;

EPtemplate.handles.trialNum=uicontrol('Style','popupmenu',...
    'String',num2cell([1:EPtemplate.maxSegs]),...
    'Value',EPtemplate.trial,'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[170 EPtemplate.windowHeight-655 160 30],'Callback', @menuTrial);

EPtemplate.handles.leftDataset = uicontrol('Style', 'pushbutton', 'String', '<--','FontSize',EPmain.fontsize,...
    'Position', [120 EPtemplate.windowHeight-680 50 30], 'Callback', @leftDataset);
if EPtemplate.dataset ==1
    set(EPtemplate.handles.leftDataset,'enable','off');
end;

EPtemplate.handles.rightDataset = uicontrol('Style', 'pushbutton', 'String', '-->','FontSize',EPmain.fontsize,...
    'Position', [330 EPtemplate.windowHeight-680 50 30], 'Callback', @rightDataset);
if EPtemplate.dataset == length(EPtemplate.sessionFiles)
    set(EPtemplate.handles.rightDataset,'enable','off');
end;

EPtemplate.handles.dataset = uicontrol('Style','popupmenu',...
    'String',{EPdataset.dataset(EPtemplate.sessionFiles).dataName},'FontSize',EPmain.fontsize,...
    'Value', EPtemplate.dataset,'Position',[170 EPtemplate.windowHeight-685 160 30],...
    'Callback', @menuDataset);

if strcmp(EPtemplate.EPdata.dataType,'continuous')
    EPtemplate.theTrial=squeeze(EPtemplate.EPdata.data(:,(EPtemplate.trial-1)*EPtemplate.segLength+1:EPtemplate.trial*EPtemplate.segLength,:,1,1));
    EPtemplate.startEpoch=(EPtemplate.trial-1)*EPtemplate.segLength+1;
    EPtemplate.endEpoch=EPtemplate.trial*EPtemplate.segLength;
else
    EPtemplate.theTrial=squeeze(EPtemplate.EPdata.data(:,:,EPtemplate.trial,1,1));
    EPtemplate.startEpoch=1;
    EPtemplate.endEpoch=length(EPtemplate.EPdata.timeNames);
end
EPtemplate.sampleLength=1000/EPtemplate.EPdata.Fs;

%baseline correct
baseline=max(EPtemplate.EPdata.baseline,1);
EPtemplate.theTrial(EPtemplate.EEGchans,:)=EPtemplate.theTrial(EPtemplate.EEGchans,:)-repmat(mean(EPtemplate.theTrial(EPtemplate.EEGchans,1:baseline),2),1,size(EPtemplate.theTrial,2));

%place markers at point most likely to be an artifact
[C EPtemplate.blinkMarker]=max(max(abs(EPtemplate.theTrial(EPtemplate.VEOG,:))));
[C EPtemplate.saccadeMarker]=max(abs(EPtemplate.theTrial(EPtemplate.HEOG(1),:)-EPtemplate.theTrial(EPtemplate.HEOG(2),:)));
[C EPtemplate.VsaccadeMarker]=max(abs(EPtemplate.theTrial(EPtemplate.VEOG(1),:)-EPtemplate.theTrial(EPtemplate.VEOG(3),:)+EPtemplate.theTrial(EPtemplate.VEOG(2),:)-EPtemplate.theTrial(EPtemplate.VEOG(4),:)));

%saccade potential
cannonicalSPthresh=25; %microvolt threshold for detecting an SP when constructing the automatic template
blinksign = [1 1 -1 -1]';
fourMS=ceil(4/(1000/(EPtemplate.EPdata.Fs))); %how many samples provides at least 4 ms
[theCzChan theOrder] = ep_findChan(EPtemplate.EPdata.eloc, EPtemplate.badChans, [0 -9.7989 107.9359]);
diffData=diff(EPtemplate.theTrial([EPtemplate.VEOG theCzChan],:),1,2);
theData=squeeze(diffData)-repmat(diffData(5,:),5,1); %rereference to Cz
theData=theData(1:4,:); %drop Cz, which is now a flat line
EPtemplate.sacPotMarker=1;
for iPoint=1:size(theData,2)-2
    if all(theData(:,iPoint) <= -cannonicalSPthresh) && all(theData(:,iPoint+(2*fourMS)) >= cannonicalSPthresh)
        if (mean(EPtemplate.theTrial(EPtemplate.VEOG(blinksign == 1),iPoint)) - mean(EPtemplate.theTrial(EPtemplate.VEOG(blinksign == -1),iPoint))) < 200
            EPtemplate.sacPotMarker=iPoint;
            break
        end;
    end;
end;

%Waveform plot of blink channels
EPtemplate.screenCoords(1,:)=[150 EPtemplate.windowHeight-120 200 100];
EPtemplate.handles.blinkPlot = axes('units','pixels','position',EPtemplate.screenCoords(1,:));
EPtemplate.handles.blinkWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.VEOG,:));
axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minVolt EPtemplate.maxVolt]);
EPtemplate.handles.blinkMarker=line(repmat(EPtemplate.blinkMarker,length([EPtemplate.minVolt:EPtemplate.maxVolt]),1),[EPtemplate.minVolt:EPtemplate.maxVolt],'Color','red','LineWidth',1); %marker
for i=1:length(EPtemplate.handles.blinkWaves)
    set(EPtemplate.handles.blinkWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
end;

%2D plot of blink channels

trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.blinkMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.blinkTopo = axes('units','pixels','position',[400 EPtemplate.windowHeight-120 100 100]);
EPtemplate.handles.blinkTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.blinkTopoImage,'ButtonDownFcn',@D3headPlot);

%2D plot of blink template
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.blinks.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.blink = axes('units','pixels','position',[550 EPtemplate.windowHeight-120 100 100]);
EPtemplate.handles.blinkImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.blinkImage,'ButtonDownFcn',@D3headPlot);

EPtemplate.handles.addBlink = uicontrol('Style', 'pushbutton', 'String', 'Add Blink','FontSize',EPmain.fontsize,...
    'Position', [550 EPtemplate.windowHeight-150 100 30], 'Callback', @addBlink);

EPtemplate.handles.blinkNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.blinks.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-145 30 20]);

%Waveform plot of H saccade channels
EPtemplate.screenCoords(2,:)=[150 EPtemplate.windowHeight-260 200 100];
EPtemplate.handles.saccadePlot = axes('units','pixels','position',EPtemplate.screenCoords(2,:));
EPtemplate.handles.saccadeWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.HEOG,:));
axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minSacVolt EPtemplate.maxSacVolt]);
EPtemplate.handles.saccadeMarker=line(repmat(EPtemplate.saccadeMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
for i=1:length(EPtemplate.handles.saccadeWaves)
    set(EPtemplate.handles.saccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.HEOG(i)) ',:)'])
end;

%2D plot of H saccade channels
trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.saccadeMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.saccadeTopo = axes('units','pixels','position',[400 EPtemplate.windowHeight-260 100 100]);
EPtemplate.handles.saccadeTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.saccadeTopoImage,'ButtonDownFcn',@D3headPlot);

%2D plot of H saccade template
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.blinks.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.saccade = axes('units','pixels','position',[550 EPtemplate.windowHeight-260 100 100]);
EPtemplate.handles.saccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.saccadeImage,'ButtonDownFcn',@D3headPlot);

EPtemplate.handles.addSaccade = uicontrol('Style', 'pushbutton', 'String', 'Add H Saccade','FontSize',EPmain.fontsize,...
    'Position', [550 EPtemplate.windowHeight-290 100 30], 'Callback', @addSaccade);

EPtemplate.handles.saccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.saccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-285 30 20]);

%Waveform plot of V saccade channels
EPtemplate.screenCoords(3,:)=[150 EPtemplate.windowHeight-400 200 100];
EPtemplate.handles.VsaccadePlot = axes('units','pixels','position',EPtemplate.screenCoords(3,:));
EPtemplate.handles.VsaccadeWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.VEOG,:));
axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minSacVolt EPtemplate.maxSacVolt]);
EPtemplate.handles.VsaccadeMarker=line(repmat(EPtemplate.VsaccadeMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
for i=1:length(EPtemplate.handles.VsaccadeWaves)
    set(EPtemplate.handles.VsaccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
end;

%2D plot of V saccade channels
trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.VsaccadeMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.VsaccadeTopo = axes('units','pixels','position',[400 EPtemplate.windowHeight-400 100 100]);
EPtemplate.handles.VsaccadeTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.VsaccadeTopoImage,'ButtonDownFcn',@D3headPlot);

%2D plot of V saccade template
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.Vsaccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.Vsaccade = axes('units','pixels','position',[550 EPtemplate.windowHeight-400 100 100]);
EPtemplate.handles.VsaccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.VsaccadeImage,'ButtonDownFcn',@D3headPlot);

EPtemplate.handles.addVSaccade = uicontrol('Style', 'pushbutton', 'String', 'Add V Saccade','FontSize',EPmain.fontsize,...
    'Position', [550 EPtemplate.windowHeight-430 100 30], 'Callback', @addVSaccade);

EPtemplate.handles.VsaccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.Vsaccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-425 30 20]);

%Waveform plot of saccade potential channels
EPtemplate.screenCoords(4,:)=[150 EPtemplate.windowHeight-540 200 100];
EPtemplate.handles.sacPotPlot = axes('units','pixels','position',EPtemplate.screenCoords(4,:));
EPtemplate.handles.sacPotWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.VEOG,:));
axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minSacVolt EPtemplate.maxSacVolt]);
EPtemplate.handles.sacPotMarker=line(repmat(EPtemplate.sacPotMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
for i=1:length(EPtemplate.handles.sacPotWaves)
    set(EPtemplate.handles.sacPotWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
end;

%2D plot of saccade potential channels
trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.sacPotMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.sacPotTopo = axes('units','pixels','position',[400 EPtemplate.windowHeight-540 100 100]);
EPtemplate.handles.sacPotTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.sacPotTopoImage,'ButtonDownFcn',@D3headPlot);

%2D plot of saccade potential template
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.sacPot.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

EPtemplate.handles.sacPot = axes('units','pixels','position',[550 EPtemplate.windowHeight-540 100 100]);
EPtemplate.handles.sacPotImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.sacPotImage,'ButtonDownFcn',@D3headPlot);

EPtemplate.handles.addSacPot = uicontrol('Style', 'pushbutton', 'String', 'Add Saccade Potentials','FontSize',EPmain.fontsize,...
    'Position', [550 EPtemplate.windowHeight-570 100 30], 'Callback', @addSacPot);

EPtemplate.handles.sacPotNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.sacPot.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-565 30 20]);

set(EPtemplate.handles.window,'WindowButtonDownFcn',@clickGraph,'WindowButtonMotionFcn',@dragGraph,'WindowButtonUpFcn',@unClickGraph);

initializeET

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function done(src,eventdata)
%quit from template window

global EPtemplate EPmain

if ~EPtemplate.saved
    choice = questdlg('Would you like to save the template first, or else lose your work?', ...
        'Save Template?', ...
        'Yes','No','Cancel','Yes');
    % Handle response
    switch choice
        case 'Yes'
            saveTemplate
        case 'No'
        case 'Cancel'
            return
    end
end;

close(EPtemplate.handles.window);

set(EPmain.handles.preprocess.template,'enable','on');
set(EPmain.handles.preprocess.preprocess,'enable','on');
set(EPmain.handles.preprocess.done,'enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leftTrial(src,eventdata)
%move to previous trial

global EPtemplate

if EPtemplate.trial > 1
    EPtemplate.trial=EPtemplate.trial-1;
end;

if EPtemplate.trial ==1
    set(EPtemplate.handles.leftTrial,'enable','off');
end;

set(EPtemplate.handles.rightTrial,'enable','on');
set(EPtemplate.handles.trialNum,'Value',EPtemplate.trial);

updateData
displayTrial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rightTrial(src,eventdata)
%move to next trial

global EPtemplate

if EPtemplate.trial < EPtemplate.maxSegs
    EPtemplate.trial=EPtemplate.trial+1;
end;

if EPtemplate.trial == EPtemplate.maxSegs
    set(EPtemplate.handles.rightTrial,'enable','off');
end;

set(EPtemplate.handles.leftTrial,'enable','on');
set(EPtemplate.handles.trialNum,'Value',EPtemplate.trial);

updateData
displayTrial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leftDataset(src,eventdata)
%move to previous dataset

global EPtemplate

EPtemplate.dataset=EPtemplate.dataset-1;
set(EPtemplate.handles.dataset,'Value',EPtemplate.dataset);

changeDataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rightDataset(src,eventdata)
%move to next dataset

global EPtemplate

EPtemplate.dataset=EPtemplate.dataset+1;
set(EPtemplate.handles.dataset,'Value',EPtemplate.dataset);

changeDataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function menuDataset(src,eventdata)
%move to new dataset using menu

global EPtemplate

EPtemplate.dataset=get(EPtemplate.handles.dataset,'Value');

changeDataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeDataset(src,eventdata)
%change the dataset.

global EPtemplate EPdataset EPmain

newData=ep_loadEPdataset(EPdataset,EPtemplate.sessionFiles(EPtemplate.dataset));
oldEloc=EPtemplate.EPdata.eloc;
newEloc=newData.eloc;

%check to see if new dataset has different electrode montage than prior dataset
newMontage=0;
if (length([oldEloc.theta]) ~= length([newEloc.theta])) || (length([oldEloc.radius]) ~= length([newEloc.radius]))
    newMontage=1;
else
    if any([newEloc.theta]-[oldEloc.theta]) || any([oldEloc.radius]-[newEloc.radius])
        newMontage=1;
    end;
end;

if strcmp(newData.dataType,'continuous')
    newLength=min(250,length(newData.timeNames));
else
    newLength=length(newData.timeNames);
end;
oldStart=EPtemplate.startEpoch;

if newMontage
    if (EPtemplate.blinks.num~=0) || (EPtemplate.saccades.num~=0)
        if ~EPtemplate.saved
            choice = questdlg('Switching to new electrode montage.  Would you like to save the template first, or else lose your work?', ...
                'Save Template?', ...
                'Yes','No','Cancel','Yes');
            % Handle response
            switch choice
                case 'Yes'
                    saveTemplate
                case 'No'
                case 'Cancel'
                    delete(EPtemplate.handles.dataset)
                    EPtemplate.handles.dataset=uicontrol('Style','text',...
                        'String',EPdataset.dataset(EPtemplate.sessionFiles(EPtemplate.dataset)).dataName,'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                        'Position',[220 EPtemplate.windowHeight-450 60 30]);
                    return
            end
        end;
    end;
    EPtemplate.EPdata=newData;
    EPtemplate.EEGchans=find(strcmp('EEG',EPtemplate.EPdata.chanTypes));
    initializeTemplate
    clearTemplates
    for i=1:length(EPtemplate.handles.blinkWaves)
        set(EPtemplate.handles.blinkWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
    end;
    for i=1:length(EPtemplate.handles.saccadeWaves)
        set(EPtemplate.handles.saccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.HEOG(i)) ',:)'])
    end;
    for i=1:length(EPtemplate.handles.VsaccadeWaves)
        set(EPtemplate.handles.VsaccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
    end;
    for i=1:length(EPtemplate.handles.sacPotWaves)
        set(EPtemplate.handles.sacPotWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
    end;
    set(EPtemplate.handles.LUVEOG,'String',newData.chanNames(EPtemplate.EEGchans));
    set(EPtemplate.handles.RUVEOG,'String',newData.chanNames(EPtemplate.EEGchans));
    set(EPtemplate.handles.LLVEOG,'String',newData.chanNames(EPtemplate.EEGchans));
    set(EPtemplate.handles.RLVEOG,'String',newData.chanNames(EPtemplate.EEGchans));
    set(EPtemplate.handles.LHEOG,'String',newData.chanNames(EPtemplate.EEGchans));
    set(EPtemplate.handles.RHEOG,'String',newData.chanNames(EPtemplate.EEGchans));
    
    set(EPtemplate.handles.LUVEOG,'Value',EPtemplate.VEOG(1));
    set(EPtemplate.handles.RUVEOG,'Value',EPtemplate.VEOG(2));
    set(EPtemplate.handles.LLVEOG,'Value',EPtemplate.VEOG(3));
    set(EPtemplate.handles.RLVEOG,'Value',EPtemplate.VEOG(4));
    set(EPtemplate.handles.LHEOG,'Value',EPtemplate.HEOG(1));
    set(EPtemplate.handles.RHEOG,'Value',EPtemplate.HEOG(2));

else
    EPtemplate.EPdata=newData;
end;

if EPtemplate.dataset == 1
    set(EPtemplate.handles.leftDataset,'enable','off');
else
    set(EPtemplate.handles.leftDataset,'enable','on');
end;

if EPtemplate.dataset == length(EPtemplate.sessionFiles)
    set(EPtemplate.handles.rightDataset,'enable','off');
else
    set(EPtemplate.handles.rightDataset,'enable','on');
end;

EPtemplate.trial=1;

if strcmp(EPtemplate.EPdata.dataType,'continuous')
    EPtemplate.theTrial=squeeze(EPtemplate.EPdata.data(:,(EPtemplate.trial-1)*newLength+1:EPtemplate.trial*newLength,:,1,1));
    EPtemplate.maxSegs=floor(length(EPtemplate.EPdata.timeNames)/newLength);
    EPtemplate.startEpoch=(EPtemplate.trial-1)*newLength+1;
    EPtemplate.endEpoch=EPtemplate.trial*newLength;
else
    EPtemplate.theTrial=squeeze(EPtemplate.EPdata.data(:,:,EPtemplate.trial,1,1));
    EPtemplate.maxSegs=length(EPtemplate.EPdata.cellNames);
    EPtemplate.startEpoch=1;
    EPtemplate.endEpoch=length(EPtemplate.EPdata.timeNames);
end;

set(EPtemplate.handles.trialNum,'Value',EPtemplate.trial);
set(EPtemplate.handles.trialNum,'String',num2cell([1:EPtemplate.maxSegs]));

%if epoch length of new dataset is different from the old
if (EPtemplate.segLength ~= newLength) || (oldStart ~= EPtemplate.startEpoch)
    EPtemplate.segLength = newLength;
    cla(EPtemplate.handles.blinkPlot);
    set(gcf,'CurrentAxes',EPtemplate.handles.blinkPlot);
    EPtemplate.handles.blinkWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.VEOG,:));
    axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minVolt EPtemplate.maxVolt]);
    EPtemplate.handles.blinkMarker=line(repmat(EPtemplate.blinkMarker,length([EPtemplate.minVolt:EPtemplate.maxVolt]),1),[EPtemplate.minVolt:EPtemplate.maxVolt],'Color','red','LineWidth',1); %marker
    for i=1:length(EPtemplate.handles.blinkWaves)
        set(EPtemplate.handles.blinkWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
    end;
    cla(EPtemplate.handles.saccadePlot);
    set(gcf,'CurrentAxes',EPtemplate.handles.saccadePlot)
    EPtemplate.handles.saccadeWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.HEOG,:));
    axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minSacVolt EPtemplate.maxSacVolt]);
    EPtemplate.handles.saccadeMarker=line(repmat(EPtemplate.saccadeMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
    for i=1:length(EPtemplate.handles.saccadeWaves)
        set(EPtemplate.handles.saccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.HEOG(i)) ',:)'])
    end;
    cla(EPtemplate.handles.VsaccadePlot);
    set(gcf,'CurrentAxes',EPtemplate.handles.VsaccadePlot)
    EPtemplate.handles.VsaccadeWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.VEOG,:));
    axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minSacVolt EPtemplate.maxSacVolt]);
    EPtemplate.handles.VsaccadeMarker=line(repmat(EPtemplate.VsaccadeMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
    for i=1:length(EPtemplate.handles.VsaccadeWaves)
        set(EPtemplate.handles.VsaccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
    end;
    cla(EPtemplate.handles.sacPotPlot);
    set(gcf,'CurrentAxes',EPtemplate.handles.sacPotPlot)
    EPtemplate.handles.sacPotWaves = plot([EPtemplate.startEpoch:EPtemplate.endEpoch],EPtemplate.theTrial(EPtemplate.VEOG,:));
    axis([EPtemplate.startEpoch EPtemplate.endEpoch EPtemplate.minSacVolt EPtemplate.maxSacVolt]);
    EPtemplate.handles.sacPotMarker=line(repmat(EPtemplate.sacPotMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
    for i=1:length(EPtemplate.handles.sacPotWaves)
        set(EPtemplate.handles.sacPotWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
    end;
end;

%remove bad channels
EPtemplate.criteria.neighbors=EPmain.preferences.preprocess.neighbors;
EPtemplate.criteria.badchan=EPmain.preferences.preprocess.badchan;
EPtemplate.criteria.saturation=[-EPmain.preferences.preprocess.saturation EPmain.preferences.preprocess.saturation];
EPtemplate.criteria.badtrials=EPmain.preferences.preprocess.badtrials;
EPtemplate.badChans=ep_detectBadChans(EPtemplate.EPdata, EPtemplate.criteria,1);

set(EPtemplate.handles.leftTrial,'enable','off');
set(EPtemplate.handles.rightTrial,'enable','on');

updateData
displayTrial

if ~isempty(EPtemplate.EYEsize)
    delete(EPtemplate.handles.topos.eye);
    delete(EPtemplate.handles.ETslider)
    delete(EPtemplate.handles.blinkMarkerET);
    delete(EPtemplate.handles.saccadeMarkerET);
    delete(EPtemplate.handles.VsaccadeMarkerET);
    delete(EPtemplate.handles.sacPotMarkerET);
    EPtemplate.XEYchan=[];
    EPtemplate.YEYchan=[];
    EPtemplate.EYEsize=[];
end;

initializeET

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateData
%update data to be displayed

global EPtemplate

if strcmp(EPtemplate.EPdata.dataType,'continuous')
    EPtemplate.theTrial=squeeze(EPtemplate.EPdata.data(:,(EPtemplate.trial-1)*EPtemplate.segLength+1:EPtemplate.trial*EPtemplate.segLength,:,1,1));
    EPtemplate.startEpoch=(EPtemplate.trial-1)*EPtemplate.segLength+1;
    EPtemplate.endEpoch=EPtemplate.trial*EPtemplate.segLength;
else
    EPtemplate.theTrial=squeeze(EPtemplate.EPdata.data(:,:,EPtemplate.trial,1,1));
    EPtemplate.startEpoch=1;
    EPtemplate.endEpoch=length(EPtemplate.EPdata.timeNames);
end
EPtemplate.sampleLength=1000/EPtemplate.EPdata.Fs;

%baseline correct
baseline=max(EPtemplate.EPdata.baseline,1);
EPtemplate.theTrial(EPtemplate.EEGchans,:)=EPtemplate.theTrial(EPtemplate.EEGchans,:)-repmat(mean(EPtemplate.theTrial(EPtemplate.EEGchans,1:baseline),2),1,size(EPtemplate.theTrial,2));
[C EPtemplate.blinkMarker]=max(max(abs(EPtemplate.theTrial(EPtemplate.VEOG,:))));
[C EPtemplate.saccadeMarker]=max(abs(EPtemplate.theTrial(EPtemplate.HEOG(1),:)-EPtemplate.theTrial(EPtemplate.HEOG(2),:)));
[C EPtemplate.VsaccadeMarker]=max(abs(EPtemplate.theTrial(EPtemplate.VEOG(1),:)-EPtemplate.theTrial(EPtemplate.VEOG(3),:)+EPtemplate.theTrial(EPtemplate.VEOG(2),:)-EPtemplate.theTrial(EPtemplate.VEOG(4),:)));

%saccade potential
cannonicalSPthresh=25; %microvolt threshold for detecting an SP when constructing the automatic template
blinksign = [1 1 -1 -1]';
fourMS=ceil(4/(1000/(EPtemplate.EPdata.Fs))); %how many samples provides at least 4 ms
[theCzChan theOrder] = ep_findChan(EPtemplate.EPdata.eloc, EPtemplate.badChans, [0 -9.7989 107.9359]);
diffData=diff(EPtemplate.theTrial([EPtemplate.VEOG theCzChan],:),1,2);
theData=squeeze(diffData)-repmat(diffData(5,:),5,1); %rereference to Cz
theData=theData(1:4,:); %drop Cz, which is now a flat line
EPtemplate.sacPotMarker=1;
for iPoint=1:size(theData,2)-2
    if all(theData(:,iPoint) <= -cannonicalSPthresh) && all(theData(:,iPoint+(2*fourMS)) >= cannonicalSPthresh)
        if (mean(EPtemplate.theTrial(EPtemplate.VEOG(blinksign == 1),iPoint)) - mean(EPtemplate.theTrial(EPtemplate.VEOG(blinksign == -1),iPoint))) < 200
            EPtemplate.sacPotMarker=iPoint;
            break
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displayTrial
%display the updated trial data

global EPtemplate EPmain

delete(EPtemplate.handles.blinkMarker);
delete(EPtemplate.handles.saccadeMarker);
delete(EPtemplate.handles.VsaccadeMarker);
delete(EPtemplate.handles.sacPotMarker);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.blinkPlot)
EPtemplate.handles.blinkMarker=line(repmat(EPtemplate.blinkMarker,length([EPtemplate.minVolt:EPtemplate.maxVolt]),1),[EPtemplate.minVolt:EPtemplate.maxVolt],'Color','red','LineWidth',1); %marker
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.saccadePlot)
EPtemplate.handles.saccadeMarker=line(repmat(EPtemplate.saccadeMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.VsaccadePlot)
EPtemplate.handles.VsaccadeMarker=line(repmat(EPtemplate.VsaccadeMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.sacPotPlot)
EPtemplate.handles.sacPotMarker=line(repmat(EPtemplate.sacPotMarker,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','red','LineWidth',1); %marker

for i=1:length(EPtemplate.handles.blinkWaves)
    refreshdata(EPtemplate.handles.blinkWaves(i))
end;
for i=1:length(EPtemplate.handles.saccadeWaves)
    refreshdata(EPtemplate.handles.saccadeWaves(i))
end;
for i=1:length(EPtemplate.handles.VsaccadeWaves)
    refreshdata(EPtemplate.handles.VsaccadeWaves(i))
end;
for i=1:length(EPtemplate.handles.sacPotWaves)
    refreshdata(EPtemplate.handles.sacPotWaves(i))
end;

%2D blink plot
trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.blinkMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.blinkTopoImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.blinkTopo)
EPtemplate.handles.blinkTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.blinkTopoImage,'ButtonDownFcn',@D3headPlot);

%2D saccade plot
trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.saccadeMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.saccadeTopoImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.saccadeTopo)
EPtemplate.handles.saccadeTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.saccadeTopoImage,'ButtonDownFcn',@D3headPlot);

%2D V saccade plot
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.VsaccadeMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.VsaccadeTopoImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.VsaccadeTopo)
EPtemplate.handles.VsaccadeTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.VsaccadeTopoImage,'ButtonDownFcn',@D3headPlot);

%2D saccade potential plot
[Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),EPtemplate.theTrial(trialGoodChans,EPtemplate.sacPotMarker),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.sacPotTopoImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.sacPotTopo)
EPtemplate.handles.sacPotTopoImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.sacPotTopoImage,'ButtonDownFcn',@D3headPlot);

if ~isempty(EPtemplate.EPdata.events)
    if strcmp(EPtemplate.EPdata.dataType,'continuous')
        theTrial=1;
        thePoints=[(EPtemplate.trial-1)*EPtemplate.segLength+1:EPtemplate.trial*EPtemplate.segLength];
    else 
        theTrial=EPtemplate.trial;
        thePoints=[1:EPtemplate.segLength];
    end
    eventList=find(ismember([EPtemplate.EPdata.events{theTrial}.sample],thePoints));
    for iLabel=1:length(EPtemplate.handles.events)
        if EPtemplate.handles.events(iLabel) ~= 0
            delete(EPtemplate.handles.events(iLabel));
        end;
    end;
    EPtemplate.handles.events=[];
    for iEvent=1:length(eventList)
        if ~isempty(EPtemplate.EPdata.events{theTrial}(eventList(iEvent)).value)
            EPtemplate.handles.events(iEvent)=uicontrol('Style','text','HorizontalAlignment','left','String', EPtemplate.EPdata.events{theTrial}(eventList(iEvent)).value,'FontSize',EPmain.fontsize,...
                'Position',[150+(200*(EPtemplate.EPdata.events{theTrial}(eventList(iEvent)).sample-thePoints(1)+1)/(thePoints(end)-thePoints(1)+1)) EPtemplate.windowHeight-20 50 15]);
        end;
    end;
end;

if ~isempty(EPtemplate.EYEsize)
    Xdata=EPtemplate.theTrial(EPtemplate.XEYchan,:);
    Ydata=EPtemplate.theTrial(EPtemplate.YEYchan,:);
    Xdata=((Xdata-EPtemplate.eyeCenterX)/EPtemplate.eyeScaleX)+.5;
    Ydata=((Ydata-EPtemplate.eyeCenterY)/EPtemplate.eyeScaleY)+.5;
    cla(EPtemplate.handles.topos.eye);
    axes(EPtemplate.handles.topos.eye);
    set(EPtemplate.handles.topos.eye,'Units','normalized')
    axis([0 1 0 1])
    set(EPtemplate.handles.topos.eye,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
    hold on
    plot(Xdata,Ydata,'color','black');
    plot(Xdata(EPtemplate.ETpoint),Ydata(EPtemplate.ETpoint),'color','green','Marker','*','MarkerSize',2);
    hold off
    axes(EPtemplate.handles.blinkPlot)
    delete(EPtemplate.handles.blinkMarkerET);
    EPtemplate.handles.blinkMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minVolt:EPtemplate.maxVolt]),1),[EPtemplate.minVolt:EPtemplate.maxVolt],'Color','green','LineWidth',1); %marker
    axes(EPtemplate.handles.saccadePlot)
    delete(EPtemplate.handles.saccadeMarkerET);
    EPtemplate.handles.saccadeMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','green','LineWidth',1); %marker
    axes(EPtemplate.handles.VsaccadePlot)
    delete(EPtemplate.handles.VsaccadeMarkerET);
    EPtemplate.handles.VsaccadeMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','green','LineWidth',1); %marker
    axes(EPtemplate.handles.sacPotPlot)
    delete(EPtemplate.handles.sacPotMarkerET);
    EPtemplate.handles.sacPotMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','green','LineWidth',1); %marker
end;

drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addBlink(src,eventdata)
%add data at marker to the blink template

global EPtemplate

if strcmp(EPtemplate.EPdata.dataType,'continuous')
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,(EPtemplate.trial-1)*EPtemplate.segLength+1:EPtemplate.trial*EPtemplate.segLength,:,1,1));
else
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,:,EPtemplate.trial,1,1));
end

%baseline correct
baseline=max(EPtemplate.EPdata.baseline,1);
theTrial=theTrial-diag(mean(theTrial(:,1:baseline),2))*ones(size(theTrial));

theBlink=theTrial(:,EPtemplate.blinkMarker);

%if there are bad channels, then first interpolate them
if ~isempty(EPtemplate.badChans)
    trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
    [Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),theBlink(trialGoodChans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'v4');
    for theBadchan=1:length(EPtemplate.badChans)
        badchan=EPtemplate.badChans(theBadchan);
        theBlink(badchan)=Zi(EPtemplate.y(badchan),EPtemplate.x(badchan));
    end;
end;

EPtemplate.undo.blinks.template=EPtemplate.blinks.template;
EPtemplate.undo.blinks.num=EPtemplate.blinks.num;

EPtemplate.blinks.template=EPtemplate.blinks.template*(EPtemplate.blinks.num/(EPtemplate.blinks.num+1))+theBlink/(EPtemplate.blinks.num+1);
EPtemplate.blinks.num=EPtemplate.blinks.num+1;

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.blinks.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.blinkImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.blink)
EPtemplate.handles.blinkImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.blinkImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.blinkNum);
EPtemplate.handles.blinkNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.blinks.num),'HorizontalAlignment','left',...
    'Position',[660 EPtemplate.windowHeight-145 30 20]);

EPtemplate.saved=0;
EPtemplate.lastChange=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addSaccade(src,eventdata)
%add data at marker to the saccade template

global EPtemplate

if strcmp(EPtemplate.EPdata.dataType,'continuous')
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,(EPtemplate.trial-1)*EPtemplate.segLength+1:EPtemplate.trial*EPtemplate.segLength,:,1,1));
else
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,:,EPtemplate.trial,1,1));
end

%baseline correct
baseline=max(EPtemplate.EPdata.baseline,1);
theTrial=theTrial-diag(mean(theTrial(:,1:baseline),2))*ones(size(theTrial));

theSaccade=theTrial(:,EPtemplate.saccadeMarker);

%if there are bad channels, then first interpolate them
if ~isempty(EPtemplate.badChans)
    trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
    [Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),theSaccade(trialGoodChans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'v4');
    for theBadchan=1:length(EPtemplate.badChans)
        badchan=EPtemplate.badChans(theBadchan);
        theSaccade(badchan)=Zi(EPtemplate.y(badchan),EPtemplate.x(badchan));
    end;
end;


if (theSaccade(EPtemplate.HEOG(2))+theSaccade(EPtemplate.VEOG(2))+theSaccade(EPtemplate.VEOG(4))-theSaccade(EPtemplate.HEOG(1))-theSaccade(EPtemplate.VEOG(1))-theSaccade(EPtemplate.VEOG(3))) > 0 %if rightward saccade
%     [Xi,Yi,Zi] = griddata(EPtemplate.x,EPtemplate.y,theSaccade,[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'v4');
%     for theChan=1:length(theSaccade)
%         leftSaccade(theChan,1)=Zi(EPtemplate.y(theChan),EPtemplate.gridSize-EPtemplate.x(theChan)+1); %flip the scalp topography
%     end;
    leftSaccade=-theSaccade;
else
    leftSaccade=theSaccade;
end

EPtemplate.undo.saccades.template=EPtemplate.saccades.template;
EPtemplate.undo.saccades.num=EPtemplate.saccades.num;

EPtemplate.saccades.template=EPtemplate.saccades.template*(EPtemplate.saccades.num/(EPtemplate.saccades.num+1))+leftSaccade/(EPtemplate.saccades.num+1);
EPtemplate.saccades.num=EPtemplate.saccades.num+1;

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.saccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.saccadeImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.saccade)
EPtemplate.handles.saccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.saccadeImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.saccadeNum);
EPtemplate.handles.saccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.saccades.num),'HorizontalAlignment','left',...
    'Position',[660 EPtemplate.windowHeight-285 30 20]);

EPtemplate.saved=0;
EPtemplate.lastChange=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addVSaccade(src,eventdata)
%add data at marker to the V saccade template

global EPtemplate

if strcmp(EPtemplate.EPdata.dataType,'continuous')
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,(EPtemplate.trial-1)*EPtemplate.segLength+1:EPtemplate.trial*EPtemplate.segLength,:,1,1));
else
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,:,EPtemplate.trial,1,1));
end

%baseline correct
baseline=max(EPtemplate.EPdata.baseline,1);
theTrial=theTrial-diag(mean(theTrial(:,1:baseline),2))*ones(size(theTrial));

theSaccade=theTrial(:,EPtemplate.VsaccadeMarker);

%if there are bad channels, then first interpolate them
if ~isempty(EPtemplate.badChans)
    trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
    [Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),theSaccade(trialGoodChans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'v4');
    for theBadchan=1:length(EPtemplate.badChans)
        badchan=EPtemplate.badChans(theBadchan);
        theSaccade(badchan)=Zi(EPtemplate.y(badchan),EPtemplate.x(badchan));
    end;
end;


if (theSaccade(EPtemplate.VEOG(1))+theSaccade(EPtemplate.VEOG(2))-theSaccade(EPtemplate.VEOG(3))-theSaccade(EPtemplate.VEOG(4))) > 0 %if upward saccade
%     [Xi,Yi,Zi] = griddata(EPtemplate.x,EPtemplate.y,theSaccade,[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'v4');
%     for theChan=1:length(theSaccade)
%         leftSaccade(theChan,1)=Zi(EPtemplate.y(theChan),EPtemplate.gridSize-EPtemplate.x(theChan)+1); %flip the scalp topography
%     end;
    downSaccade=-theSaccade;
else
    downSaccade=theSaccade;
end

EPtemplate.undo.Vsaccades.template=EPtemplate.Vsaccades.template;
EPtemplate.undo.Vsaccades.num=EPtemplate.Vsaccades.num;

EPtemplate.Vsaccades.template=EPtemplate.Vsaccades.template*(EPtemplate.Vsaccades.num/(EPtemplate.Vsaccades.num+1))+downSaccade/(EPtemplate.Vsaccades.num+1);
EPtemplate.Vsaccades.num=EPtemplate.Vsaccades.num+1;

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.Vsaccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.VsaccadeImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.Vsaccade)
EPtemplate.handles.VsaccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.VsaccadeImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.VsaccadeNum);
EPtemplate.handles.VsaccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.Vsaccades.num),'HorizontalAlignment','left',...
    'Position',[660 EPtemplate.windowHeight-425 30 20]);

EPtemplate.saved=0;
EPtemplate.lastChange=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addSacPot(src,eventdata)
%add data at marker to the Saccade Potential template

global EPtemplate

if strcmp(EPtemplate.EPdata.dataType,'continuous')
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,(EPtemplate.trial-1)*EPtemplate.segLength+1:EPtemplate.trial*EPtemplate.segLength,:,1,1));
else
    theTrial=squeeze(EPtemplate.EPdata.data(EPtemplate.EEGchans,:,EPtemplate.trial,1,1));
end

%baseline correct
baseline=max(EPtemplate.EPdata.baseline,1);
theTrial=theTrial-diag(mean(theTrial(:,1:baseline),2))*ones(size(theTrial));

theSacPot=theTrial(:,EPtemplate.sacPotMarker);

%if there are bad channels, then first interpolate them
if ~isempty(EPtemplate.badChans)
    trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
    [Xi,Yi,Zi] = griddata(EPtemplate.x(trialGoodChans),EPtemplate.y(trialGoodChans),theSacPot(trialGoodChans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'v4');
    for theBadchan=1:length(EPtemplate.badChans)
        badchan=EPtemplate.badChans(theBadchan);
        theSacPot(badchan)=Zi(EPtemplate.y(badchan),EPtemplate.x(badchan));
    end;
end;

EPtemplate.undo.sacPot.template=EPtemplate.sacPot.template;
EPtemplate.undo.sacPot.num=EPtemplate.sacPot.num;

EPtemplate.sacPot.template=EPtemplate.sacPot.template*(EPtemplate.sacPot.num/(EPtemplate.sacPot.num+1))+theSacPot/(EPtemplate.sacPot.num+1);
EPtemplate.sacPot.num=EPtemplate.sacPot.num+1;

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.sacPot.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.sacPotImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.sacPot)
EPtemplate.handles.sacPotImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.sacPotImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.sacPotNum);
EPtemplate.handles.sacPotNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.sacPot.num),'HorizontalAlignment','left',...
    'Position',[660 EPtemplate.windowHeight-565 30 20]);

EPtemplate.saved=0;
EPtemplate.lastChange=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadTemplate(src,eventdata)
%load artifact templates

global EPtemplate EPdataset EPmain

[FileName,PathName,FilterIndex] = uigetfile('*.mat','Load Blink Template','blinks.mat');
if FileName ~= 0
    eval(['load ''' PathName FileName '''']);
    if ~exist('EPblink','var')
        msg{1}='Not a blink template.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if length(EPblink.eloc) ~= length(EPtemplate.blinks.eloc)
        msg{1}='Number of blink electrodes different from current dataset.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if any([EPtemplate.blinks.eloc.theta]-[EPtemplate.EPdata.eloc.theta])
        msg{1}='Blink electrode locations not consistent with current data.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if ~exist('EPblink','var')
        EPblink.template=zeros(length(EPtemplate.EPdata.chanNames),1);
        EPblink.num=0;
    end;

    EPtemplate.undo.blinks.template=EPtemplate.blinks.template;
    EPtemplate.undo.blinks.num=EPtemplate.blinks.num;
    
    EPtemplate.blinks=EPblink;
end;

[FileName,PathName,FilterIndex] = uigetfile('*.mat','Load Saccade Template','saccades.mat');
if FileName ~= 0
    eval(['load ''' PathName FileName '''']);
    
    if ~exist('EPsaccade','var')
        msg{1}='Not a saccade template.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if length(EPsaccade.eloc) ~= length(EPtemplate.saccades.eloc)
        msg{1}=['Number of saccade electrodes (' num2str(length(EPsaccade.eloc)) ') different from current dataset (' num2str(length(EPtemplate.saccades.eloc)) ').'];
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if any([EPsaccade.eloc.theta]-[EPtemplate.saccades.eloc.theta])
        msg{1}='Saccade electrode locations not consistent with current data.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    EPtemplate.undo.saccades.template=EPtemplate.saccades.template;
    EPtemplate.undo.saccades.num=EPtemplate.saccades.num;
    if ~isfield(EPsaccade,'hSaccade') %backward compatibility
        EPsaccade.hSaccade.template=EPsaccade.template;
        EPsaccade.hSaccade.num=EPsaccade.num;
    end;
    EPtemplate.saccades=EPsaccade.hSaccade;
    EPtemplate.saccades.eloc=EPsaccade.eloc;
    
    if ~isfield(EPsaccade,'vSaccade')
        EPsaccade.vSaccade.template=zeros(length(EPtemplate.EPdata.chanNames),1);
        EPsaccade.vSaccade.num=0;
    end;
    
    EPtemplate.undo.Vsaccades.template=EPtemplate.Vsaccades.template;
    EPtemplate.undo.Vsaccades.num=EPtemplate.Vsaccades.num;
    EPtemplate.Vsaccades=EPsaccade.vSaccade;
    
    if ~isfield(EPsaccade,'sacPot')
        EPsaccade.sacPot.template=zeros(length(EPtemplate.EPdata.chanNames),1);
        EPsaccade.sacPot.num=0;
    end;

    EPtemplate.undo.sacPot.template=EPtemplate.sacPot.template;
    EPtemplate.undo.sacPot.num=EPtemplate.sacPot.num;
    EPtemplate.sacPot=EPsaccade.sacPot;
end;

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.blinks.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.blinkImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.blink)
EPtemplate.handles.blinkImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.blinkImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.blinkNum);
EPtemplate.handles.blinkNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.blinks.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-145 30 20]);

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.saccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.saccadeImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.saccade)
EPtemplate.handles.saccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.saccadeImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.saccadeNum);
EPtemplate.handles.saccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.saccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-285 30 20]);

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.Vsaccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.VsaccadeImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.Vsaccade)
EPtemplate.handles.VsaccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.VsaccadeImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.VsaccadeNum);
EPtemplate.handles.VsaccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.Vsaccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-425 30 20]);

[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.Vsaccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');

delete(EPtemplate.handles.sacPotImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.sacPot)
EPtemplate.handles.sacPotImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.sacPotImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.sacPotNum);
EPtemplate.handles.sacPotNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.sacPot.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-565 30 20]);

EPtemplate.saved=1;
EPtemplate.lastChange=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveTemplate(src,eventdata)
%save blink and saccade templates

global EPtemplate

if EPtemplate.blinks.num > 0
    [FileName,PathName,FilterIndex] = uiputfile('*.mat','Save Templates','blinks.mat');
    if FileName ==0
        
    else
        [pathstr, name, ext] = fileparts(FileName);
        if strcmp(ext,'.rtw') %fix for 2012a Matlab bug
            ext='.mat';
        end;
        FileName=[name ext];
        
        EPblink=EPtemplate.blinks;
        eval(['save ''' PathName FileName ''' EPblink']);
    end;
end;

if (EPtemplate.saccades.num > 0) || (EPtemplate.Vsaccades.num > 0) || (EPtemplate.sacPot.num > 0)
    [FileName,PathName,FilterIndex] = uiputfile('*.mat','Save Templates','saccades.mat');
    if FileName ==0
        
    else
        [pathstr, name, ext] = fileparts(FileName);
        if strcmp(ext,'.rtw') %fix for 2012a Matlab bug
            ext='.mat';
        end;
        FileName=[name ext];
        
        EPsaccade.hSaccade.template=EPtemplate.saccades.template;
        EPsaccade.hSaccade.num=EPtemplate.saccades.num;
        EPsaccade.vSaccade.template=EPtemplate.Vsaccades.template;
        EPsaccade.vSaccade.num=EPtemplate.Vsaccades.num;
        EPsaccade.sacPot.template=EPtemplate.sacPot.template;
        EPsaccade.sacPot.num=EPtemplate.sacPot.num;
        EPsaccade.eloc=EPtemplate.saccades.eloc;
        eval(['save ''' PathName FileName ''' EPsaccade']);
    end;
end;

EPtemplate.saved=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function undo(src,eventdata)
%undo last change to templates

global EPtemplate EPmain

if (EPtemplate.lastChange==1) || (EPtemplate.lastChange==100)
    
    tempVar.blinks.template=EPtemplate.blinks.template;
    tempVar.blinks.num=EPtemplate.blinks.num;
    EPtemplate.blinks.template=EPtemplate.undo.blinks.template;
    EPtemplate.blinks.num=EPtemplate.undo.blinks.num;
    EPtemplate.undo.blinks.template=tempVar.blinks.template;
    EPtemplate.undo.blinks.num=tempVar.blinks.num;
    delete(EPtemplate.handles.blinkImage);
    set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.blink)
    [Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.blinks.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
    EPtemplate.handles.blinkImage = imagesc(Zi);
    set(gca,'XTickLabel','','YTickLabel','');
    set(EPtemplate.handles.blinkImage,'ButtonDownFcn',@D3headPlot);
    delete(EPtemplate.handles.blinkNum);
    EPtemplate.handles.blinkNum=uicontrol('Style','text',...
        'String',sprintf('%d',EPtemplate.blinks.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'Position',[660 EPtemplate.windowHeight-145 30 20]);
    
end;

if (EPtemplate.lastChange==2) || (EPtemplate.lastChange==100)
    
    tempVar.saccades.template=EPtemplate.saccades.template;
    tempVar.saccades.num=EPtemplate.saccades.num;
    EPtemplate.saccades.template=EPtemplate.undo.saccades.template;
    EPtemplate.saccades.num=EPtemplate.undo.saccades.num;
    EPtemplate.undo.saccades.template=tempVar.saccades.template;
    EPtemplate.undo.saccades.num=tempVar.saccades.num;
    delete(EPtemplate.handles.saccadeImage);
    set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.saccade)
    [Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.saccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
    EPtemplate.handles.saccadeImage = imagesc(Zi);
    set(gca,'XTickLabel','','YTickLabel','');
    set(EPtemplate.handles.saccadeImage,'ButtonDownFcn',@D3headPlot);
    
    delete(EPtemplate.handles.saccadeNum);
    EPtemplate.handles.saccadeNum=uicontrol('Style','text',...
        'String',sprintf('%d',EPtemplate.saccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'Position',[660 EPtemplate.windowHeight-285 30 20]);
    
end;

if (EPtemplate.lastChange==3) || (EPtemplate.lastChange==100)
    
    tempVar.Vsaccades.template=EPtemplate.Vsaccades.template;
    tempVar.Vsaccades.num=EPtemplate.Vsaccades.num;
    EPtemplate.Vsaccades.template=EPtemplate.undo.Vsaccades.template;
    EPtemplate.Vsaccades.num=EPtemplate.undo.Vsaccades.num;
    EPtemplate.undo.Vsaccades.template=tempVar.Vsaccades.template;
    EPtemplate.undo.Vsaccades.num=tempVar.Vsaccades.num;
    delete(EPtemplate.handles.VsaccadeImage);
    set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.Vsaccade)
    [Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.Vsaccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
    EPtemplate.handles.VsaccadeImage = imagesc(Zi);
    set(gca,'XTickLabel','','YTickLabel','');
    set(EPtemplate.handles.VsaccadeImage,'ButtonDownFcn',@D3headPlot);
    
    delete(EPtemplate.handles.VsaccadeNum);
    EPtemplate.handles.VsaccadeNum=uicontrol('Style','text',...
        'String',sprintf('%d',EPtemplate.Vsaccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'Position',[660 EPtemplate.windowHeight-425 30 20]);
    
end;

if (EPtemplate.lastChange==4) || (EPtemplate.lastChange==100)
    
    tempVar.sacPot.template=EPtemplate.sacPot.template;
    tempVar.sacPot.num=EPtemplate.sacPot.num;
    EPtemplate.sacPot.template=EPtemplate.undo.sacPot.template;
    EPtemplate.sacPot.num=EPtemplate.undo.sacPot.num;
    EPtemplate.undo.sacPot.template=tempVar.sacPot.template;
    EPtemplate.undo.sacPot.num=tempVar.sacPot.num;
    delete(EPtemplate.handles.sacPotImage);
    set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.sacPot)
    [Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.sacPot.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
    EPtemplate.handles.sacPotImage = imagesc(Zi);
    set(gca,'XTickLabel','','YTickLabel','');
    set(EPtemplate.handles.sacPotImage,'ButtonDownFcn',@D3headPlot);
    
    delete(EPtemplate.handles.sacPotNum);
    EPtemplate.handles.sacPotNum=uicontrol('Style','text',...
        'String',sprintf('%d',EPtemplate.sacPot.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'Position',[660 EPtemplate.windowHeight-425 30 20]);
    
end;

EPtemplate.saved=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clearTemplates(src,eventdata)
%clear templates

global EPtemplate EPmain

EPtemplate.undo.blinks.template=EPtemplate.blinks.template;
EPtemplate.undo.blinks.num=EPtemplate.blinks.num;
EPtemplate.undo.saccades.template=EPtemplate.saccades.template;
EPtemplate.undo.saccades.num=EPtemplate.saccades.num;
EPtemplate.undo.Vsaccades.template=EPtemplate.Vsaccades.template;
EPtemplate.undo.Vsaccades.num=EPtemplate.Vsaccades.num;
EPtemplate.undo.sacPot.template=EPtemplate.sacPot.template;
EPtemplate.undo.sacPot.num=EPtemplate.sacPot.num;

EPtemplate.blinks.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.blinks.num=0;
EPtemplate.saccades.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.saccades.num=0;
EPtemplate.Vsaccades.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.Vsaccades.num=0;
EPtemplate.sacPot.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.sacPot.num=0;

delete(EPtemplate.handles.blinkImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.blink)
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.blinks.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
EPtemplate.handles.blinkImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.blinkImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.blinkNum);
EPtemplate.handles.blinkNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.blinks.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-145 30 20]);

delete(EPtemplate.handles.saccadeImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.saccade)
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.saccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
EPtemplate.handles.saccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.saccadeImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.saccadeNum);
EPtemplate.handles.saccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.saccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-285 30 20]);

delete(EPtemplate.handles.VsaccadeImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.Vsaccade)
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.Vsaccades.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
EPtemplate.handles.VsaccadeImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.VsaccadeImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.VsaccadeNum);
EPtemplate.handles.VsaccadeNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.Vsaccades.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-425 30 20]);

delete(EPtemplate.handles.sacPotImage);
set(EPtemplate.handles.window,'CurrentAxes',EPtemplate.handles.sacPot)
[Xi,Yi,Zi] = griddata(EPtemplate.x(EPtemplate.EEGchans),EPtemplate.y(EPtemplate.EEGchans),EPtemplate.sacPot.template(EPtemplate.EEGchans),[1:EPtemplate.gridSize]',[1:EPtemplate.gridSize],'linear');
EPtemplate.handles.sacPotImage = imagesc(Zi);
set(gca,'XTickLabel','','YTickLabel','');
set(EPtemplate.handles.sacPotImage,'ButtonDownFcn',@D3headPlot);

delete(EPtemplate.handles.sacPotNum);
EPtemplate.handles.sacPotNum=uicontrol('Style','text',...
    'String',sprintf('%d',EPtemplate.sacPot.num),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[660 EPtemplate.windowHeight-565 30 20]);

EPtemplate.saved=0;
EPtemplate.lastChange=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initializeTemplate
%set up new template

global EPtemplate EPmain

if (length(EPmain.preferences.preprocess.EOGchans) ~=6) || (max(EPmain.preferences.preprocess.EOGchans) > length(EPtemplate.EPdata.eloc))
    eog=ep_findEOGchans(EPtemplate.EPdata.eloc);
else
    eog.LUVEOG=EPmain.preferences.preprocess.EOGchans(1);
    eog.RUVEOG=EPmain.preferences.preprocess.EOGchans(2);
    eog.LLVEOG=EPmain.preferences.preprocess.EOGchans(3);
    eog.RLVEOG=EPmain.preferences.preprocess.EOGchans(4);
    eog.LHEOG=EPmain.preferences.preprocess.EOGchans(5);
    eog.RHEOG=EPmain.preferences.preprocess.EOGchans(6);
end;

EPtemplate.VEOG=[eog.LUVEOG eog.RUVEOG eog.LLVEOG eog.RLVEOG];
EPtemplate.HEOG=[eog.LHEOG eog.RHEOG];

EPtemplate.blinks.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.undo.blinks.template=EPtemplate.blinks.template;
EPtemplate.blinks.eloc=EPtemplate.EPdata.eloc(EPtemplate.EEGchans);
EPtemplate.blinks.num=0;
EPtemplate.undo.blinks.num=EPtemplate.blinks.num;

EPtemplate.saccades.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.undo.saccades.template=EPtemplate.saccades.template;
EPtemplate.saccades.eloc=EPtemplate.EPdata.eloc(EPtemplate.EEGchans);
EPtemplate.saccades.num=0;
EPtemplate.undo.saccades.num=EPtemplate.saccades.num;

EPtemplate.Vsaccades.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.undo.Vsaccades.template=EPtemplate.Vsaccades.template;
EPtemplate.Vsaccades.eloc=EPtemplate.EPdata.eloc(EPtemplate.EEGchans);
EPtemplate.Vsaccades.num=0;
EPtemplate.undo.Vsaccades.num=EPtemplate.Vsaccades.num;

EPtemplate.sacPot.template=zeros(length(EPtemplate.EEGchans),1);
EPtemplate.undo.sacPot.template=EPtemplate.sacPot.template;
EPtemplate.sacPot.eloc=EPtemplate.EPdata.eloc(EPtemplate.EEGchans);
EPtemplate.sacPot.num=0;
EPtemplate.undo.sacPot.num=EPtemplate.sacPot.num;

maxRad=0.5;
[y,x] = pol2cart(([EPtemplate.EPdata.eloc.theta]/360)*2*pi,[EPtemplate.EPdata.eloc.radius]);  % transform electrode locations from polar to cartesian coordinates
y=-y; %flip y-coordinate so that nose is upwards.
plotrad = min(1.0,max([EPtemplate.EPdata.eloc.radius])*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
x = x*(maxRad/plotrad);
y = y*(maxRad/plotrad);

xmin = min(-maxRad,min(x));
xmax = max(maxRad,max(x));
ymin = min(-maxRad,min(y));
ymax = max(maxRad,max(y));

EPtemplate.x=round(((x/(xmax-xmin))*EPtemplate.gridSize)+ceil(EPtemplate.gridSize/2));
EPtemplate.y=round(((y/(ymax-ymin))*EPtemplate.gridSize)+ceil(EPtemplate.gridSize/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D3headPlot(src,eventdata)
global EPtemplate

disp('Using EEGlab function headplot to perform 3D head display.');

trialGoodChans=setdiff(EPtemplate.EEGchans,EPtemplate.badChans);
switch src
    case EPtemplate.handles.blinkTopoImage
        theData=EPtemplate.theTrial(trialGoodChans,EPtemplate.blinkMarker);
    case EPtemplate.handles.blinkImage
        theData=EPtemplate.blinks.template(EPtemplate.EEGchans);
    case EPtemplate.handles.saccadeTopoImage
        theData=EPtemplate.theTrial(trialGoodChans,EPtemplate.saccadeMarker);
    case EPtemplate.handles.saccadeImage
        theData=EPtemplate.blinks.template(EPtemplate.EEGchans);
    case EPtemplate.handles.VsaccadeTopoImage
        theData=EPtemplate.theTrial(trialGoodChans,EPtemplate.VsaccadeMarker);
    case EPtemplate.handles.VsaccadeImage
        theData=EPtemplate.Vsaccades.template(EPtemplate.EEGchans);
    case EPtemplate.handles.sacPotTopoImage
        theData=EPtemplate.theTrial(trialGoodChans,EPtemplate.sacPotMarker);
    case EPtemplate.handles.sacPotImage
        theData=EPtemplate.sacPot.template(EPtemplate.EEGchans);
end;

%check to see if spline file needs to be generated
tempFlag=0;
if any(EPtemplate.badChans)
    name='badChanSpline';
    CEDloc=[pwd filesep 'temp'];
    tempFlag=1;
    disp('Since there are bad channels, will need to generate a new spline file');
elseif ~isempty(EPtemplate.EPdata.ced) && ~any(strcmp(EPtemplate.EPdata.ced,{'none','internal'}))
    [pathstr, name, ext] = fileparts(EPtemplate.EPdata.ced);
    CEDloc=which(EPtemplate.EPdata.ced);
else
    name=EPtemplate.EPdata.dataName;
    CEDloc=[pwd filesep 'temp'];
end;
if isempty(which([name '.spl'])) || tempFlag
    [pathstr2, name2, ext2] = fileparts(CEDloc);
    if max(abs([EPtemplate.EPdata.eloc.X])) <= 1
        MNItransform=[ 0 -15 0 0.08 0 -1.571 102 93 100 ]; %assume eloc coordinates are from a .elp file
    else
        MNItransform=[ 0 -15 4 0.05 0 -1.571 10.2 12 12.2 ]; %assume eloc coordinates are from a .sfp file
    end;
    headplot('setup', EPtemplate.EPdata.eloc(trialGoodChans), [pathstr2 filesep name '.spl'],'transform',MNItransform); %save spline file in same place as the ced file.
end;

figure
%[hdaxis cbaraxis] = headplot(theData,which([name '.spl']),'cbar',0,'maplimits',[EPtemplate.minVolt EPtemplate.maxVolt],'labels',2);
[hdaxis cbaraxis] = headplot(theData,which([name '.spl']),'cbar',0,'maplimits',[-max(abs(theData)) max(abs(theData))],'labels',2);

if tempFlag
    delete([pwd filesep 'badChanSpline.spl']);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeEOG(src,eventdata)
global EPtemplate

EPtemplate.VEOG(1)=get(EPtemplate.handles.LUVEOG,'Value');
EPtemplate.VEOG(2)=get(EPtemplate.handles.RUVEOG,'Value');
EPtemplate.VEOG(3)=get(EPtemplate.handles.LLVEOG,'Value');
EPtemplate.VEOG(4)=get(EPtemplate.handles.RLVEOG,'Value');
EPtemplate.HEOG(1)=get(EPtemplate.handles.LHEOG,'Value');
EPtemplate.HEOG(2)=get(EPtemplate.handles.RHEOG,'Value');

for i=1:length(EPtemplate.handles.blinkWaves)
    set(EPtemplate.handles.blinkWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
end;
for i=1:length(EPtemplate.handles.saccadeWaves)
    set(EPtemplate.handles.saccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.HEOG(i)) ',:)'])
end;
for i=1:length(EPtemplate.handles.VsaccadeWaves)
    set(EPtemplate.handles.VsaccadeWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
end;
for i=1:length(EPtemplate.handles.sacPotWaves)
    set(EPtemplate.handles.sacPotWaves(i),'YDataSource',['EPtemplate.theTrial(' num2str(EPtemplate.VEOG(i)) ',:)'])
end;
displayTrial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function menuTrial(src,eventdata)
%move to new dataset using menu

global EPtemplate

EPtemplate.trial=get(EPtemplate.handles.trialNum,'Value');

set(EPtemplate.handles.leftTrial,'enable','on');
set(EPtemplate.handles.rightTrial,'enable','on');
if EPtemplate.trial ==1
    set(EPtemplate.handles.leftTrial,'enable','off');
end;
if EPtemplate.trial == EPtemplate.maxSegs
    set(EPtemplate.handles.rightTrial,'enable','off');
end;

updateData
displayTrial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clickGraph(~,~)
%respond to click in graph area

global EPtemplate

if ~strcmpi(get(EPtemplate.handles.window,'selectiontype'),'normal')
    return %ignore if not left-click
end;

pos = get(EPtemplate.handles.window, 'CurrentPoint');
for iGraph=1:size(EPtemplate.screenCoords,1)
    if (pos(2) >= EPtemplate.screenCoords(iGraph,2)) && (pos(2) <= EPtemplate.screenCoords(iGraph,2)+EPtemplate.screenCoords(iGraph,4)) && ...
            (pos(1) >= EPtemplate.screenCoords(iGraph,1)) && (pos(1) <= EPtemplate.screenCoords(iGraph,1)+EPtemplate.screenCoords(iGraph,3))
        EPtemplate.clickGraph=iGraph;
        dragGraph
        break
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dragGraph(~,~)
%respond to drag in graph area

global EPtemplate

if ~EPtemplate.clickGraph
    return
end;

pos=get(EPtemplate.handles.window, 'CurrentPoint');

clickX=ceil(((pos(1)-EPtemplate.screenCoords(EPtemplate.clickGraph,1))/EPtemplate.screenCoords(EPtemplate.clickGraph,3))*EPtemplate.segLength);
if clickX < 1
    clickX=1;
end;
if clickX > EPtemplate.segLength
    clickX=EPtemplate.segLength;
end;

switch EPtemplate.clickGraph
    case 1 %blink
        EPtemplate.blinkMarker=clickX;
    case 2 %H saccade
        EPtemplate.saccadeMarker=clickX;
    case 3 %V saccade
        EPtemplate.VsaccadeMarker=clickX;
    case 4 %saccade potential
        EPtemplate.sacPotMarker=clickX;
end;

displayTrial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function unClickGraph(~,~)
%respond to unclick in graph area

global EPtemplate

EPtemplate.clickGraph=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initializeET
%initialize eye-tracker figure

global EPtemplate

EPtemplate.XEYchan=find(strcmp('XEY',EPtemplate.EPdata.chanTypes));
if length(EPtemplate.XEYchan) > 1
    EPtemplate.XEYchan=EPtemplate.XEYchan(1);
end;
EPtemplate.YEYchan=find(strcmp('YEY',EPtemplate.EPdata.chanTypes));
if length(EPtemplate.YEYchan) > 1
    EPtemplate.YEYchan=EPtemplate.YEYchan(1);
end;
if ~isempty(EPtemplate.XEYchan) && ~isempty(EPtemplate.YEYchan)
    EPtemplate.EYEsize=max([max(abs(EPtemplate.EPdata.data(EPtemplate.XEYchan,:))) max(abs(EPtemplate.EPdata.data(EPtemplate.YEYchan,:)))]);
    if isfield(EPtemplate.EPdata.calibration,'ET') && isfield(EPtemplate.EPdata.calibration.ET,'Xzero') && isfield(EPtemplate.EPdata.calibration.ET,'Yzero') && isfield(EPtemplate.EPdata.calibration.ET,'Xscale') && isfield(EPtemplate.EPdata.calibration.ET,'Yscale')
        EPtemplate.eyeCenterX=EPtemplate.EPdata.calibration.ET.Xzero;
        EPtemplate.eyeCenterY=EPtemplate.EPdata.calibration.ET.Yzero;
        EPtemplate.eyeScaleX=EPtemplate.EPdata.calibration.ET.Xscale;
        EPtemplate.eyeScaleY=EPtemplate.EPdata.calibration.ET.Yscale;
    else
        nonNaNx=~isnan(EPtemplate.EPdata.data(EPtemplate.XEYchan,:));
        nonNaNy=~isnan(EPtemplate.EPdata.data(EPtemplate.YEYchan,:));
        EPtemplate.eyeCenterX=median(EPtemplate.EPdata.data(EPtemplate.XEYchan,nonNaNx));
        EPtemplate.eyeCenterY=median(EPtemplate.EPdata.data(EPtemplate.YEYchan,nonNaNy));
        EPtemplate.eyeScaleX=max(EPtemplate.EPdata.data(EPtemplate.XEYchan,nonNaNx)-EPtemplate.eyeCenterX)-min(EPtemplate.EPdata.data(EPtemplate.XEYchan,nonNaNx)-EPtemplate.eyeCenterX);
        EPtemplate.eyeScaleY=max(EPtemplate.EPdata.data(EPtemplate.YEYchan,nonNaNx)-EPtemplate.eyeCenterY)-min(EPtemplate.EPdata.data(EPtemplate.YEYchan,nonNaNx)-EPtemplate.eyeCenterY);
    end;
    EPtemplate.ETpoint=round(EPtemplate.segLength/2);
    EPtemplate.ETslider=.5;
    EPtemplate.handles.topos.eye = axes('Units','pixel','position',[400 30 100 100]);
    EPtemplate.handles.ETslider = uicontrol('Style', 'slider', 'Value',EPtemplate.ETslider,...
    'Position', [400 10 100 20],'SliderStep',[(1/EPtemplate.segLength) 10*(1/EPtemplate.segLength)], 'Min',0,'Max',1,'Callback', @ETslider);
else
    EPtemplate.EYEsize=[];
end;

if ~isempty(EPtemplate.EYEsize)
    EPtemplate.handles.blinkMarkerET=[];
    EPtemplate.handles.saccadeMarkerET=[];
    EPtemplate.handles.VsaccadeMarkerET=[];
    EPtemplate.handles.sacPotMarkerET=[];
    ETupdate
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ETslider(~,~)
%respond to eye-tracker figure slider

global EPtemplate

EPtemplate.ETslider=get(EPtemplate.handles.ETslider,'Value');
EPtemplate.ETpoint=round(EPtemplate.ETslider*(EPtemplate.segLength-1))+1;

ETupdate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ETupdate(~,~)
%update ET displays

global EPtemplate

Xdata=EPtemplate.theTrial(EPtemplate.XEYchan,:);
Ydata=EPtemplate.theTrial(EPtemplate.YEYchan,:);
Xdata=((Xdata-EPtemplate.eyeCenterX)/EPtemplate.eyeScaleX)+.5;
Ydata=((Ydata-EPtemplate.eyeCenterY)/EPtemplate.eyeScaleY)+.5;
cla(EPtemplate.handles.topos.eye);
axes(EPtemplate.handles.topos.eye);
set(EPtemplate.handles.topos.eye,'Units','normalized')
axis([0 1 0 1])
set(EPtemplate.handles.topos.eye,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
hold on
plot(Xdata,Ydata,'color','black');
plot(Xdata(EPtemplate.ETpoint),Ydata(EPtemplate.ETpoint),'color','green','Marker','*','MarkerSize',2);
hold off
axes(EPtemplate.handles.blinkPlot)
delete(EPtemplate.handles.blinkMarkerET);
EPtemplate.handles.blinkMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minVolt:EPtemplate.maxVolt]),1),[EPtemplate.minVolt:EPtemplate.maxVolt],'Color','green','LineWidth',1); %marker
axes(EPtemplate.handles.saccadePlot)
delete(EPtemplate.handles.saccadeMarkerET);
EPtemplate.handles.saccadeMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','green','LineWidth',1); %marker
axes(EPtemplate.handles.VsaccadePlot)
delete(EPtemplate.handles.VsaccadeMarkerET);
EPtemplate.handles.VsaccadeMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','green','LineWidth',1); %marker
axes(EPtemplate.handles.sacPotPlot)
delete(EPtemplate.handles.sacPotMarkerET);
EPtemplate.handles.sacPotMarkerET=line(repmat(EPtemplate.ETpoint,length([EPtemplate.minSacVolt:EPtemplate.maxSacVolt]),1),[EPtemplate.minSacVolt:EPtemplate.maxSacVolt],'Color','green','LineWidth',1); %marker
set(EPtemplate.handles.ETslider,'Value',EPtemplate.ETslider);
