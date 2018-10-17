function ep_trimData
% ep_trimData - ep_trimData -
% Provides graph view of data for manually trimming continuous data.
%
%Input:
%

%History
%  by Joseph Dien (3/12/14)
%  jdien07@mac.com
%
% bufix 7/28/14 JD
% Fixed edits not being saved when clicking on "keep".
%
%  modified 5/25/14 JD
%  Set colormap to jet even for Matlab 2014b onwards.

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

global EPdataset EPtrimData EPmain

try
    set(EPmain.handles.segment.trim,'enable','off');
    set(EPmain.handles.segment.preview,'enable','off');
    set(EPmain.handles.segment.segment,'enable','off');
    set(EPmain.handles.segment.main,'enable','off');
catch
    ep('start')
end;

scrsz = EPmain.scrsz;

EPtrimData.edited=0;

EPmain.handles.trimData = figure('Name', 'Trim Data', 'NumberTitle', 'off', 'Position',[270 1 scrsz(3)-270 scrsz(4)], 'MenuBar', 'none');
colormap jet;
drawnow

EPtrimData.EPdata=ep_loadEPdataset(EPdataset,EPmain.segment.dataset);

EPtrimData.numChans=length(EPtrimData.EPdata.chanNames);
EPtrimData.numPoints=length(EPtrimData.EPdata.timeNames);
EPtrimData.xGraphMin=min(EPtrimData.numPoints,ceil(EPtrimData.EPdata.Fs)); %minimum number of time points to be displayable
EPtrimData.displayLength=min(EPtrimData.numPoints,ceil(EPtrimData.EPdata.Fs)*4); %start at displaying 4 seconds worth of data
EPtrimData.spacing=100; %+/-microvolt space for each waveform
EPtrimData.displayChans=[1:EPtrimData.numChans]; %channels to be displayed
EPtrimData.displayPoints=[1:EPtrimData.displayLength]; %points to be displayed
EPtrimData.xGraphPos=0; %current x position of the graph from 0 to 1
EPtrimData.yGraphPos=0; %current y position of the graph from 0 to 1
EPtrimData.currentData=[]; %data to be displayed

EPtrimData.ColorOrder=[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; .75 0 .75; .75 .75 0; .25 .25 .25];
EPtrimData.graphCoords=[200 200 EPmain.scrsz(3)-600 EPmain.scrsz(4)-300];

EPtrimData.handles.graphPlot = axes('Units','pixels','position',EPtrimData.graphCoords,'Tag','graphPlot');
if EPmain.preferences.view.positive ==2
    set(EPtrimData.handles.graphPlot,'YDir','reverse')
end;

EPtrimData.handles.deletionZone=[];
EPtrimData.deletionZone.x=[];
EPtrimData.deletionZone.y=[];
EPtrimData.deletionZone.status=0;
EPtrimData.deletionZone.dataSamples=[];

EPtrimData.handles.deleteMenu = uicontextmenu;
hcb1 = [@deleteZone];
item1 = uimenu(EPtrimData.handles.deleteMenu, 'Label', 'delete', 'Callback', hcb1);

for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,1:EPtrimData.displayLength)-EPtrimData.EPdata.data(theChan,1,:,:,:)-EPtrimData.spacing*(iChan-1);
end;
EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(:,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])
set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

EPtrimData.handles.chanNames=[];
for iChan=1:length(EPtrimData.displayChans)
    theColor=EPtrimData.ColorOrder(rem(EPtrimData.displayChans(iChan)-1,size(EPtrimData.ColorOrder,1))+1,:);
    EPtrimData.handles.chanNames(iChan)=uicontrol('Style','text','HorizontalAlignment','left','String', EPtrimData.EPdata.chanNames{EPtrimData.displayChans(iChan)},'FontSize',EPmain.fontsize,'ForegroundColor',theColor,...
    'Position',[140 EPmain.scrsz(4)-120-((EPmain.scrsz(4)-300)*(length(EPtrimData.displayChans)/(length(EPtrimData.displayChans)+1))/length(EPtrimData.displayChans))*iChan 40 20]);
end;

EPtrimData.handles.xSlider = uicontrol('Style', 'slider', 'Value',EPtrimData.xGraphPos,...
    'Position', [200 180 scrsz(3)-600 20],'SliderStep',[EPtrimData.displayLength/EPtrimData.numPoints max(0.2,length(EPtrimData.displayLength)/EPtrimData.numPoints)], 'Min',0,'Max',1,'Callback', @xSlider);

if EPtrimData.numPoints <= EPtrimData.displayLength
    set(EPtrimData.handles.xSlider,'enable','off');
end;

EPtrimData.handles.ySlider = uicontrol('Style', 'slider', 'Value',EPtrimData.yGraphPos,...
    'Position', [180 200 20 scrsz(4)-300],'SliderStep',[length(EPtrimData.displayChans)/EPtrimData.numChans max(0.2,length(EPtrimData.displayChans)/EPtrimData.numChans)], 'Min',0,'Max',1,'Callback', @ySlider);

if EPtrimData.numChans == length(EPtrimData.displayChans)
    set(EPtrimData.handles.ySlider,'enable','off');
end;

if ~isempty(EPtrimData.EPdata.events)
    eventList=find(ismember([EPtrimData.EPdata.events{1}.sample],EPtrimData.displayPoints));
    for iEvent=1:length(eventList)
        if ~isempty(EPtrimData.EPdata.events{1}(eventList(iEvent)).value)
            EPtrimData.handles.events(iEvent)=uicontrol('Style','text','HorizontalAlignment','left','String', EPtrimData.EPdata.events{1}(eventList(iEvent)).value,'FontSize',EPmain.fontsize,...
                'Position',[200+((EPmain.scrsz(3)-600)*(EPtrimData.displayLength/(EPtrimData.displayLength+1))/EPtrimData.displayLength)*EPtrimData.EPdata.events{1}(iEvent).sample EPmain.scrsz(4)-100 50 20]);
            if strcmp(EPtrimData.EPdata.events{1}(eventList(iEvent)).value,'boundary')
                theColor='red';
            else
                theColor='black';
            end;
            
            EPtrimData.handles.eventLines(iEvent)=line([EPtrimData.EPdata.events{1}(eventList(iEvent)).sample-EPtrimData.displayPoints(1)+1 EPtrimData.EPdata.events{1}(eventList(iEvent)).sample-EPtrimData.displayPoints(1)+1],[EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing],'Color',theColor,'LineWidth',1);
        end;
    end;
end;

EPtrimData.handles.displayStart=uicontrol('Style','text','HorizontalAlignment','left','String', ep_ms2time((EPtrimData.displayPoints(1)-1)*(1000/EPtrimData.EPdata.Fs)),'FontSize',EPmain.fontsize,...
    'Position',[200 160 70 20]);

EPtrimData.handles.displayEnd=uicontrol('Style','text','HorizontalAlignment','left','String', ep_ms2time(EPtrimData.displayPoints(end)*(1000/EPtrimData.EPdata.Fs)),'FontSize',EPmain.fontsize,...
    'Position',[EPmain.scrsz(3)-470 160 70 20]);

EPtrimData.handles.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
    'Position', [92 5 50 35], 'Callback', @cancel);

EPtrimData.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Keep','FontSize',EPmain.fontsize,...
    'Position', [152 5 50 35], 'Callback', @done);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Channels','FontSize',EPmain.fontsize,...
    'Position',[212 40 60 20]);

EPtrimData.handles.lessChans = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*2,...
    'Position', [212 10 30 30], 'Callback', @lessChans);

EPtrimData.handles.moreChans = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*2,...
    'Position', [242 10 30 30], 'Callback', @moreChans);

set(EPtrimData.handles.moreChans,'enable','off');

uicontrol('Style','text','HorizontalAlignment','left','String', 'Time','FontSize',EPmain.fontsize,...
    'Position',[292 40 60 20]);

EPtrimData.handles.lessTime = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*2,...
    'Position', [292 10 30 30], 'Callback', @lessTime);

EPtrimData.handles.moreTime = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*2,...
    'Position', [322 10 30 30], 'Callback', @moreTime);

uicontrol('Style','text','HorizontalAlignment','left','String', 'Deletion Zone','FontSize',EPmain.fontsize,...
    'Position',[372 40 100 20]);

EPtrimData.handles.preEvent = uicontrol('Style', 'pushbutton', 'String', '|<-','FontSize',EPmain.fontsize*2,...
    'Position', [372 10 30 30], 'Callback', @preEvent);

EPtrimData.handles.preSample = uicontrol('Style', 'pushbutton', 'String', '<-','FontSize',EPmain.fontsize*2,...
    'Position', [402 10 30 30], 'Callback', @preSample);

EPtrimData.handles.deletionStart=uicontrol('Style','text','HorizontalAlignment','left','String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)),'FontSize',EPmain.fontsize,...
    'Position',[432 10 70 20]);

EPtrimData.handles.deletionEnd=uicontrol('Style','text','HorizontalAlignment','left','String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)),'FontSize',EPmain.fontsize,...
    'Position',[502 10 70 20]);

EPtrimData.handles.postSample = uicontrol('Style', 'pushbutton', 'String', '->','FontSize',EPmain.fontsize*2,...
    'Position', [572 10 30 30], 'Callback', @postSample);

EPtrimData.handles.postEvent = uicontrol('Style', 'pushbutton', 'String', '->|','FontSize',EPmain.fontsize*2,...
    'Position', [602 10 30 30], 'Callback', @postEvent);

set(EPtrimData.handles.preEvent,'enable','off');
set(EPtrimData.handles.preSample,'enable','off');
set(EPtrimData.handles.postSample,'enable','off');
set(EPtrimData.handles.postEvent,'enable','off');

set(EPmain.handles.trimData,'WindowButtonDownFcn',@clickGraph,'WindowButtonMotionFcn',@dragGraph,'WindowButtonUpFcn',@unClickGraph);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cancel(~,~)
%quit from trim data window without saving edits

global EPmain EPtrimData 

close(EPmain.handles.trimData);

EPtrimData=[];

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function done(~,~)
%quit from trim data window and save edits

global EPmain EPtrimData EPdataset

close(EPmain.handles.trimData);

if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
    msg{1}='The work directory cannot be found.';
    [msg]=ep_errorMsg(msg);
    return
end;

if EPtrimData.edited
    delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(EPmain.segment.dataset).dataName '.mat']);
    EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],EPmain.segment.dataset));
    EPdataset=ep_saveEPdataset(EPdataset,EPtrimData.EPdata,EPmain.segment.dataset,'no');
end;

EPtrimData=[];

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xSlider(~,~)
%update graph after xSlider moved

global EPtrimData EPmain

EPtrimData.xGraphPos=get(EPtrimData.handles.xSlider,'Value');

oldDisplayPoints=EPtrimData.displayPoints;
EPtrimData.displayPoints=1+floor(EPtrimData.xGraphPos*EPtrimData.numPoints):EPtrimData.displayLength+floor(EPtrimData.xGraphPos*EPtrimData.numPoints);
if max(EPtrimData.displayPoints) > EPtrimData.numPoints
    EPtrimData.displayPoints=[max(1,EPtrimData.numPoints-EPtrimData.displayLength+1):EPtrimData.numPoints];
end;

for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints)-EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints(1),:,:,:)-EPtrimData.spacing*(iChan-1);
end;

EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(EPtrimData.displayChans,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])
set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

drawEvents(1);

set(EPtrimData.handles.displayStart,'String', ep_ms2time((EPtrimData.displayPoints(1)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.displayEnd,'String', ep_ms2time(EPtrimData.displayPoints(end)*(1000/EPtrimData.EPdata.Fs)));

if ~isempty(EPtrimData.deletionZone.x) && ~isempty(EPtrimData.deletionZone.y)
    EPtrimData.deletionZone.x=EPtrimData.deletionZone.x+oldDisplayPoints(1)-EPtrimData.displayPoints(1);
    EPtrimData.handles.deletionZone=patch(EPtrimData.deletionZone.x,EPtrimData.deletionZone.y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
    alpha(EPtrimData.handles.deletionZone,.5);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ySlider(~,~)
%update graph after ySlider moved

global EPtrimData EPmain

EPtrimData.yGraphPos=get(EPtrimData.handles.ySlider,'Value');

EPtrimData.displayChans=1+floor((1-EPtrimData.yGraphPos)*EPtrimData.numChans):length(EPtrimData.displayChans)+floor((1-EPtrimData.yGraphPos)*EPtrimData.numChans);
if max(EPtrimData.displayChans) > EPtrimData.numChans
    EPtrimData.displayChans=[max(1,EPtrimData.numChans-length(EPtrimData.displayChans)+1):EPtrimData.numChans];
end;

theFirstColor=rem(EPtrimData.displayChans(1)-1,size(EPtrimData.ColorOrder,1))+1;
theColors=[EPtrimData.ColorOrder(theFirstColor:end,:); EPtrimData.ColorOrder(1:theFirstColor-1,:)];

for iChan=1:length(EPtrimData.handles.chanNames)
    delete(EPtrimData.handles.chanNames(iChan));
end;

EPtrimData.handles.chanNames=[];

for iChan=1:length(EPtrimData.displayChans)
    theColor=EPtrimData.ColorOrder(rem(EPtrimData.displayChans(iChan)-1,size(EPtrimData.ColorOrder,1))+1,:);
    EPtrimData.handles.chanNames(iChan)=uicontrol('Style','text','HorizontalAlignment','left','String', EPtrimData.EPdata.chanNames{EPtrimData.displayChans(iChan)},'FontSize',EPmain.fontsize,'ForegroundColor',theColor,...
    'Position',[140 EPmain.scrsz(4)-120-((EPmain.scrsz(4)-300)*(length(EPtrimData.displayChans)/(length(EPtrimData.displayChans)+1))/length(EPtrimData.displayChans))*iChan 40 20]);
end;

for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints)-EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints(1),:,:,:)-EPtrimData.spacing*(iChan-1);
end;

set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[],'ColorOrder',theColors,'NextPlot','replacechildren');
EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(EPtrimData.displayChans,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])

drawEvents(0);

if ~isempty(EPtrimData.deletionZone.x) && ~isempty(EPtrimData.deletionZone.y)
    bottomSide=-length(EPtrimData.displayChans)*EPtrimData.spacing;
    height=(length(EPtrimData.displayChans)+1)*EPtrimData.spacing;
    EPtrimData.deletionZone.y=[bottomSide+height;bottomSide;bottomSide;bottomSide+height];
    EPtrimData.handles.deletionZone=patch(EPtrimData.deletionZone.x,EPtrimData.deletionZone.y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
    alpha(EPtrimData.handles.deletionZone,.5);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lessChans(~,~)
%zoom in (less channels)

global EPtrimData EPmain

if length(EPtrimData.displayChans) > 1
    EPtrimData.displayChans=[EPtrimData.displayChans(1):EPtrimData.displayChans(1)+floor(length(EPtrimData.displayChans)/2)-1];
end

theFirstColor=rem(EPtrimData.displayChans(1)-1,size(EPtrimData.ColorOrder,1))+1;
theColors=[EPtrimData.ColorOrder(theFirstColor:end,:); EPtrimData.ColorOrder(1:theFirstColor-1,:)];

for iChan=1:length(EPtrimData.handles.chanNames)
    delete(EPtrimData.handles.chanNames(iChan));
end;

EPtrimData.handles.chanNames=[];

for iChan=1:length(EPtrimData.displayChans)
    theColor=EPtrimData.ColorOrder(rem(EPtrimData.displayChans(iChan)-1,size(EPtrimData.ColorOrder,1))+1,:);    EPtrimData.handles.chanNames(iChan)=uicontrol('Style','text','HorizontalAlignment','left','String', EPtrimData.EPdata.chanNames{EPtrimData.displayChans(iChan)},'FontSize',EPmain.fontsize,'ForegroundColor',theColor,...
    'Position',[140 EPmain.scrsz(4)-120-((EPmain.scrsz(4)-300)*(length(EPtrimData.displayChans)/(length(EPtrimData.displayChans)+1))/length(EPtrimData.displayChans))*iChan 40 20]);
end;

if EPtrimData.numChans > length(EPtrimData.displayChans)
    set(EPtrimData.handles.ySlider,'enable','on');
end;

set(EPtrimData.handles.ySlider,'SliderStep',[length(EPtrimData.displayChans)/EPtrimData.numChans max(0.2,length(EPtrimData.displayChans)/EPtrimData.numChans)],'Value',1-(EPtrimData.displayChans(1)/EPtrimData.numChans));

set(EPtrimData.handles.moreChans,'enable','on');
if length(EPtrimData.displayChans) == 1
    set(EPtrimData.handles.lessChans,'enable','off');
else
    set(EPtrimData.handles.lessChans,'enable','on');
end;

for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints)-EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints(1),:,:,:)-EPtrimData.spacing*(iChan-1);
end;

EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(EPtrimData.displayChans,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])
set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

drawEvents(0);

if ~isempty(EPtrimData.deletionZone.x) && ~isempty(EPtrimData.deletionZone.x)
    bottomSide=-length(EPtrimData.displayChans)*EPtrimData.spacing;
    height=(length(EPtrimData.displayChans)+1)*EPtrimData.spacing;
    EPtrimData.deletionZone.y=[bottomSide+height;bottomSide;bottomSide;bottomSide+height];
    EPtrimData.handles.deletionZone=patch(EPtrimData.deletionZone.x,EPtrimData.deletionZone.y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
    alpha(EPtrimData.handles.deletionZone,.5);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function moreChans(~,~)
%zoom out (more channels)

global EPtrimData EPmain

EPtrimData.displayChans=[EPtrimData.displayChans(1):EPtrimData.displayChans(1)+length(EPtrimData.displayChans)*2-1];
if max(EPtrimData.displayChans) > EPtrimData.numChans
    EPtrimData.displayChans=[EPtrimData.displayChans(1):EPtrimData.numChans];
end;

theFirstColor=rem(EPtrimData.displayChans(1)-1,size(EPtrimData.ColorOrder,1))+1;
theColors=[EPtrimData.ColorOrder(theFirstColor:end,:); EPtrimData.ColorOrder(1:theFirstColor-1,:)];

for iChan=1:length(EPtrimData.handles.chanNames)
    delete(EPtrimData.handles.chanNames(iChan));
end;

EPtrimData.handles.chanNames=[];

for iChan=1:length(EPtrimData.displayChans)
    theColor=EPtrimData.ColorOrder(rem(EPtrimData.displayChans(iChan)-1,size(EPtrimData.ColorOrder,1))+1,:);
    EPtrimData.handles.chanNames(iChan)=uicontrol('Style','text','HorizontalAlignment','left','String', EPtrimData.EPdata.chanNames{EPtrimData.displayChans(iChan)},'FontSize',EPmain.fontsize,'ForegroundColor',theColor,...
    'Position',[140 EPmain.scrsz(4)-120-((EPmain.scrsz(4)-300)*(length(EPtrimData.displayChans)/(length(EPtrimData.displayChans)+1))/length(EPtrimData.displayChans))*iChan 40 20]);
end;

set(EPtrimData.handles.ySlider,'SliderStep',[length(EPtrimData.displayChans)/EPtrimData.numChans max(0.2,length(EPtrimData.displayChans)/EPtrimData.numChans)],'Value',1-(EPtrimData.displayChans(1)/EPtrimData.numChans));

if EPtrimData.numChans == length(EPtrimData.displayChans)
    set(EPtrimData.handles.ySlider,'enable','off');
end;

set(EPtrimData.handles.lessChans,'enable','on');
if length(EPtrimData.displayChans) == EPtrimData.numChans
    set(EPtrimData.handles.moreChans,'enable','off');
else
    set(EPtrimData.handles.moreChans,'enable','on');
end;

for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints)-EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints(1),:,:,:)-EPtrimData.spacing*(iChan-1);
end;

EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(EPtrimData.displayChans,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])
set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

drawEvents(0);

if ~isempty(EPtrimData.deletionZone.x) && ~isempty(EPtrimData.deletionZone.x)
    bottomSide=-length(EPtrimData.displayChans)*EPtrimData.spacing;
    height=(length(EPtrimData.displayChans)+1)*EPtrimData.spacing;
    EPtrimData.deletionZone.y=[bottomSide+height;bottomSide;bottomSide;bottomSide+height];
    EPtrimData.handles.deletionZone=patch(EPtrimData.deletionZone.x,EPtrimData.deletionZone.y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
    alpha(EPtrimData.handles.deletionZone,.5);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lessTime(~,~)
%zoom in (less time)

global EPtrimData EPmain

if length(EPtrimData.displayPoints) > EPtrimData.xGraphMin
    EPtrimData.displayPoints=[EPtrimData.displayPoints(1):EPtrimData.displayPoints(1)+floor(length(EPtrimData.displayPoints)/2)-1];
end;

if length(EPtrimData.displayPoints) < EPtrimData.xGraphMin
    EPtrimData.displayPoints=[EPtrimData.displayPoints(1):EPtrimData.displayPoints(1)+EPtrimData.xGraphMin-1];
end;

EPtrimData.displayLength=length(EPtrimData.displayPoints);
set(EPtrimData.handles.xSlider,'SliderStep',[EPtrimData.displayLength/EPtrimData.numPoints max(0.2,length(EPtrimData.displayLength)/EPtrimData.numPoints)]);

if EPtrimData.numPoints > length(EPtrimData.displayPoints)
    set(EPtrimData.handles.xSlider,'enable','on');
end;

set(EPtrimData.handles.moreTime,'enable','on');
if length(EPtrimData.displayPoints) == EPtrimData.xGraphMin
    set(EPtrimData.handles.lessTime,'enable','off');
else
    set(EPtrimData.handles.lessTime,'enable','on');
end;

EPtrimData.currentData=zeros(EPtrimData.numChans,EPtrimData.displayLength);
for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints)-EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints(1),:,:,:)-EPtrimData.spacing*(iChan-1);
end;

EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(EPtrimData.displayChans,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])
set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

drawEvents(1);

set(EPtrimData.handles.displayStart,'String', ep_ms2time((EPtrimData.displayPoints(1)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.displayEnd,'String', ep_ms2time(EPtrimData.displayPoints(end)*(1000/EPtrimData.EPdata.Fs)));

if ~isempty(EPtrimData.deletionZone.x) && ~isempty(EPtrimData.deletionZone.x)
    EPtrimData.handles.deletionZone=patch(EPtrimData.deletionZone.x,EPtrimData.deletionZone.y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
    alpha(EPtrimData.handles.deletionZone,.5);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function moreTime(~,~)
%zoom out (more time)

global EPtrimData EPmain

EPtrimData.displayPoints=[EPtrimData.displayPoints(1):EPtrimData.displayPoints(1)+length(EPtrimData.displayPoints)*2-1];

if max(EPtrimData.displayPoints) > EPtrimData.numPoints
    EPtrimData.displayPoints=[EPtrimData.displayPoints(1):EPtrimData.numPoints];
end;

EPtrimData.displayLength=length(EPtrimData.displayPoints);
set(EPtrimData.handles.xSlider,'SliderStep',[EPtrimData.displayLength/EPtrimData.numPoints max(0.2,length(EPtrimData.displayLength)/EPtrimData.numPoints)]);

if EPtrimData.numPoints == length(EPtrimData.displayPoints)
    set(EPtrimData.handles.xSlider,'enable','off');
end;

set(EPtrimData.handles.lessTime,'enable','on');
if length(EPtrimData.displayPoints) == EPtrimData.numPoints
    set(EPtrimData.handles.moreTime,'enable','off');
else
    set(EPtrimData.handles.moreTime,'enable','on');
end;

EPtrimData.currentData=zeros(EPtrimData.numChans,EPtrimData.displayLength);
for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints)-EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints(1),:,:,:)-EPtrimData.spacing*(iChan-1);
end;

EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(EPtrimData.displayChans,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])
set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

drawEvents(1);

set(EPtrimData.handles.displayStart,'String', ep_ms2time((EPtrimData.displayPoints(1)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.displayEnd,'String', ep_ms2time(EPtrimData.displayPoints(end)*(1000/EPtrimData.EPdata.Fs)));

EPtrimData.handles.deletionZone=patch(EPtrimData.deletionZone.x,EPtrimData.deletionZone.y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
alpha(EPtrimData.handles.deletionZone,.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clickGraph(~,~)
%respond to click in graph area

global EPtrimData EPmain

if ~strcmpi(get(EPmain.handles.trimData,'selectiontype'),'normal')
    return %ignore if not left-click
end;

pos = get(EPmain.handles.trimData, 'CurrentPoint');
if (pos(2) < EPtrimData.graphCoords(2)) || (pos(2) > EPtrimData.graphCoords(2)+EPtrimData.graphCoords(4))
    return
end;
if ~isempty(EPtrimData.deletionZone.x) && ~isempty(EPtrimData.deletionZone.x)
    minDeletionZone=((min(EPtrimData.deletionZone.x)/EPtrimData.displayLength)*EPtrimData.graphCoords(3))+EPtrimData.graphCoords(1); %convert to figure pixels
    maxDeletionZone=((max(EPtrimData.deletionZone.x)/EPtrimData.displayLength)*EPtrimData.graphCoords(3))+EPtrimData.graphCoords(1); %convert to figure pixels
    if (pos(1) >= minDeletionZone) && (pos(1) <= maxDeletionZone)
        %clicking within existing deletion zone
        delete(EPtrimData.handles.deletionZone);
        EPtrimData.handles.deletionZone=[];
        EPtrimData.deletionZone.x=[];
        EPtrimData.deletionZone.y=[];
        EPtrimData.deletionZone.status=0;
        EPtrimData.deletionZone.dataSamples=[];
        set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
        set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));
        set(EPtrimData.handles.preEvent,'enable','off');
        set(EPtrimData.handles.preSample,'enable','off');
        set(EPtrimData.handles.postSample,'enable','off');
        set(EPtrimData.handles.postEvent,'enable','off');
        return
    end;
end;

clickX=ceil(((pos(1)-EPtrimData.graphCoords(1))/EPtrimData.graphCoords(3))*EPtrimData.displayLength); %convert to graph samples
if clickX < 1
    clickX=1;
end;
if clickX > EPtrimData.displayLength
    clickX=EPtrimData.displayLength;
end;

bottomSide=-length(EPtrimData.displayChans)*EPtrimData.spacing;
height=(length(EPtrimData.displayChans)+1)*EPtrimData.spacing;
x=[clickX;clickX;clickX;clickX];
y=[bottomSide+height;bottomSide;bottomSide;bottomSide+height];
if ~isempty(EPtrimData.handles.deletionZone)
    delete(EPtrimData.handles.deletionZone);
    EPtrimData.handles.deletionZone=[];
end;
EPtrimData.handles.deletionZone=patch(x,y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
alpha(EPtrimData.handles.deletionZone,.5);

EPtrimData.deletionZone.x=x;
EPtrimData.deletionZone.y=y;
EPtrimData.deletionZone.status=1;
EPtrimData.deletionZone.dataSamples=[EPtrimData.displayPoints(clickX) EPtrimData.displayPoints(clickX)];


set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));

set(EPtrimData.handles.preEvent,'enable','on');
set(EPtrimData.handles.preSample,'enable','on');
set(EPtrimData.handles.postSample,'enable','on');
set(EPtrimData.handles.postEvent,'enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dragGraph(~,~)
%respond to click in graph area

global EPtrimData EPmain

if ~EPtrimData.deletionZone.status
    return
end;

pos=get(EPmain.handles.trimData, 'CurrentPoint');
clickX=ceil(((pos(1)-EPtrimData.graphCoords(1))/EPtrimData.graphCoords(3))*EPtrimData.displayLength);
if clickX < 1
    clickX=1;
end;
if clickX > EPtrimData.displayLength
    clickX=EPtrimData.displayLength;
end;

x=EPtrimData.deletionZone.x;
y=EPtrimData.deletionZone.y;

x(3)=clickX;
x(4)=clickX;

if ~isempty(EPtrimData.handles.deletionZone)
    delete(EPtrimData.handles.deletionZone);
    EPtrimData.handles.deletionZone=[];
end;

EPtrimData.handles.deletionZone=patch(x,y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
alpha(EPtrimData.handles.deletionZone,.5);

EPtrimData.deletionZone.x=x;
EPtrimData.deletionZone.y=y;
EPtrimData.deletionZone.dataSamples(2)=EPtrimData.displayPoints(clickX);

set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));

if min(EPtrimData.deletionZone.dataSamples)==1
    set(EPtrimData.handles.preEvent,'enable','off');
    set(EPtrimData.handles.preSample,'enable','off');
end;

if max(EPtrimData.deletionZone.dataSamples)==EPtrimData.numPoints
    set(EPtrimData.handles.postSample,'enable','off');
    set(EPtrimData.handles.postEvent,'enable','off');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function unClickGraph(~,~)
%respond to unclick in graph area

global EPtrimData EPmain

EPtrimData.deletionZone.status=0;
set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function preEvent(~,~)
%grow deletion zone by one earlier event

global EPtrimData

minX=min(EPtrimData.deletionZone.dataSamples);

if minX==1
    return
end;

if ~isempty(EPtrimData.EPdata.events)
    eventList=find([EPtrimData.EPdata.events{1}.sample] < min(EPtrimData.deletionZone.dataSamples));
else
    eventList=[];
end;

if isempty(eventList)
    newStart=1;
else
    newStart=max([EPtrimData.EPdata.events{1}(eventList).sample]);
end;

x=EPtrimData.deletionZone.x;
y=EPtrimData.deletionZone.y;

x(find(x == min(x)))=x(find(x == min(x)))-(min(EPtrimData.deletionZone.dataSamples)-newStart);

if ~isempty(EPtrimData.handles.deletionZone)
    delete(EPtrimData.handles.deletionZone);
    EPtrimData.handles.deletionZone=[];
end;

EPtrimData.handles.deletionZone=patch(x,y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
alpha(EPtrimData.handles.deletionZone,.5);

EPtrimData.deletionZone.x=x;
EPtrimData.deletionZone.y=y;
EPtrimData.deletionZone.dataSamples(find(EPtrimData.deletionZone.dataSamples == min(EPtrimData.deletionZone.dataSamples)))=newStart;

set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));

if min(EPtrimData.deletionZone.dataSamples)==1
    set(EPtrimData.handles.preEvent,'enable','off');
    set(EPtrimData.handles.preSample,'enable','off');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function preSample(~,~)
%grow deletion zone by one earlier sample

global EPtrimData

minX=min(EPtrimData.deletionZone.dataSamples);

if minX==1
    return
end;

x=EPtrimData.deletionZone.x;
y=EPtrimData.deletionZone.y;

x(find(x == min(x)))=x(find(x == min(x)))-1;

if ~isempty(EPtrimData.handles.deletionZone)
    delete(EPtrimData.handles.deletionZone);
    EPtrimData.handles.deletionZone=[];
end;

EPtrimData.handles.deletionZone=patch(x,y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
alpha(EPtrimData.handles.deletionZone,.5);

EPtrimData.deletionZone.x=x;
EPtrimData.deletionZone.y=y;
EPtrimData.deletionZone.dataSamples(find(EPtrimData.deletionZone.dataSamples == min(EPtrimData.deletionZone.dataSamples)))=EPtrimData.deletionZone.dataSamples(find(EPtrimData.deletionZone.dataSamples == min(EPtrimData.deletionZone.dataSamples)))-1;

set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));

if min(EPtrimData.deletionZone.dataSamples)==1
    set(EPtrimData.handles.preEvent,'enable','off');
    set(EPtrimData.handles.preSample,'enable','off');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function postSample(~,~)
%grow deletion zone by one later sample

global EPtrimData

maxX=max(EPtrimData.deletionZone.dataSamples);

if maxX==EPtrimData.numPoints
    return
end;

x=EPtrimData.deletionZone.x;
y=EPtrimData.deletionZone.y;

x(find(x == max(x)))=x(find(x == max(x)))+1;

if ~isempty(EPtrimData.handles.deletionZone)
    delete(EPtrimData.handles.deletionZone);
    EPtrimData.handles.deletionZone=[];
end;

EPtrimData.handles.deletionZone=patch(x,y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
alpha(EPtrimData.handles.deletionZone,.5);

EPtrimData.deletionZone.x=x;
EPtrimData.deletionZone.y=y;
EPtrimData.deletionZone.dataSamples(find(EPtrimData.deletionZone.dataSamples == max(EPtrimData.deletionZone.dataSamples)))=EPtrimData.deletionZone.dataSamples(find(EPtrimData.deletionZone.dataSamples == max(EPtrimData.deletionZone.dataSamples)))+1;

set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));

if max(EPtrimData.deletionZone.dataSamples)==EPtrimData.numPoints
    set(EPtrimData.handles.postSample,'enable','off');
    set(EPtrimData.handles.postEvent,'enable','off');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function postEvent(~,~)
%grow deletion zone by just shy of one later event

global EPtrimData

maxX=max(EPtrimData.deletionZone.dataSamples);

if maxX==EPtrimData.numPoints
    return
end;

if ~isempty(EPtrimData.EPdata.events)
    eventList=find([EPtrimData.EPdata.events{1}.sample] > max(EPtrimData.deletionZone.dataSamples)+1);
else
    eventList=[];
end;

if isempty(eventList)
    newStart=EPtrimData.numPoints;
else
    newStart=min([EPtrimData.EPdata.events{1}(eventList).sample])-1;
end;

x=EPtrimData.deletionZone.x;
y=EPtrimData.deletionZone.y;

x(find(x == max(x)))=x(find(x == max(x)))+(newStart-max(EPtrimData.deletionZone.dataSamples));

if ~isempty(EPtrimData.handles.deletionZone)
    delete(EPtrimData.handles.deletionZone);
    EPtrimData.handles.deletionZone=[];
end;

EPtrimData.handles.deletionZone=patch(x,y,'red','EdgeColor','red','UIContextMenu',EPtrimData.handles.deleteMenu);
alpha(EPtrimData.handles.deletionZone,.5);

EPtrimData.deletionZone.x=x;
EPtrimData.deletionZone.y=y;
EPtrimData.deletionZone.dataSamples(find(EPtrimData.deletionZone.dataSamples == max(EPtrimData.deletionZone.dataSamples)))=newStart;

set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));

if max(EPtrimData.deletionZone.dataSamples)==EPtrimData.numPoints
    set(EPtrimData.handles.postSample,'enable','off');
    set(EPtrimData.handles.postEvent,'enable','off');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deleteZone(~,~)
%Delete the deletion zone

global EPtrimData EPmain

if abs(diff(EPtrimData.deletionZone.dataSamples))+1 == EPtrimData.numPoints
    disp('Cannot delete the complete set of data');
    return
end;

[EPtrimData.EPdata]=ep_selectData(EPtrimData.EPdata,{[],[1:min(EPtrimData.deletionZone.dataSamples)-1 max(EPtrimData.deletionZone.dataSamples)+1:EPtrimData.numPoints],[],[],[],[]});

EPtrimData.numPoints=length(EPtrimData.EPdata.timeNames);
EPtrimData.displayPoints=[EPtrimData.displayPoints(1):EPtrimData.displayPoints(1)+EPtrimData.displayLength-1];

if max(EPtrimData.displayPoints) > EPtrimData.numPoints
    EPtrimData.displayPoints=[max(1,EPtrimData.numPoints-EPtrimData.displayLength+1):EPtrimData.numPoints];
    if min(EPtrimData.displayPoints) < 1
        EPtrimData.displayPoints=[1:EPtrimData.numPoints];
    end;
end;

EPtrimData.displayLength=length(EPtrimData.displayPoints);
set(EPtrimData.handles.xSlider,'SliderStep',[EPtrimData.displayLength/EPtrimData.numPoints max(0.2,length(EPtrimData.displayLength)/EPtrimData.numPoints)]);

if EPtrimData.numPoints == length(EPtrimData.displayPoints)
    set(EPtrimData.handles.xSlider,'enable','off');
end;

set(EPtrimData.handles.lessTime,'enable','on');
if length(EPtrimData.displayPoints) == EPtrimData.numPoints
    set(EPtrimData.handles.moreTime,'enable','off');
else
    set(EPtrimData.handles.moreTime,'enable','on');
end;

EPtrimData.currentData=zeros(EPtrimData.numChans,EPtrimData.displayLength);
for iChan=1:length(EPtrimData.displayChans)
    theChan=EPtrimData.displayChans(iChan);
    EPtrimData.currentData(theChan,:)=EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints)-EPtrimData.EPdata.data(theChan,EPtrimData.displayPoints(1),:,:,:)-EPtrimData.spacing*(iChan-1);
end;

EPtrimData.handles.graphChans=plot([1:EPtrimData.displayLength],EPtrimData.currentData(EPtrimData.displayChans,:));
axis([1 EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing EPtrimData.spacing])
set(EPtrimData.handles.graphPlot,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);

drawEvents(1);

set(EPtrimData.handles.displayStart,'String', ep_ms2time((EPtrimData.displayPoints(1)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.displayEnd,'String', ep_ms2time(EPtrimData.displayPoints(end)*(1000/EPtrimData.EPdata.Fs)));

EPtrimData.edited=1;
EPtrimData.handles.deletionZone=[];
EPtrimData.deletionZone.x=[];
EPtrimData.deletionZone.y=[];
EPtrimData.deletionZone.status=0;
EPtrimData.deletionZone.dataSamples=[];
set(EPtrimData.handles.deletionStart,'String', ep_ms2time((min(EPtrimData.deletionZone.x)-1)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.deletionEnd,'String', ep_ms2time(max(EPtrimData.deletionZone.x)*(1000/EPtrimData.EPdata.Fs)));
set(EPtrimData.handles.preEvent,'enable','off');
set(EPtrimData.handles.preSample,'enable','off');
set(EPtrimData.handles.postSample,'enable','off');
set(EPtrimData.handles.postEvent,'enable','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drawEvents(xChange)
%Draw the event lines

global EPtrimData EPmain

if xChange
    for iEvent=1:length(EPtrimData.handles.events)
        if isgraphics(EPtrimData.handles.events(iEvent))
            delete(EPtrimData.handles.events(iEvent));
        end;
    end;
    EPtrimData.handles.events=[];
end;
EPtrimData.handles.eventLines=[];

if ~isempty(EPtrimData.EPdata.events)
    eventList=find(ismember(round([EPtrimData.EPdata.events{1}.sample]),EPtrimData.displayPoints));
    for iEvent=1:length(eventList)
        if ~isempty(EPtrimData.EPdata.events{1}(eventList(iEvent)).value)
            if xChange
                EPtrimData.handles.events(iEvent)=uicontrol('Style','text','HorizontalAlignment','left','String', EPtrimData.EPdata.events{1}(eventList(iEvent)).value,'FontSize',EPmain.fontsize,...
                    'Position',[200+((EPmain.scrsz(3)-600)*(EPtrimData.displayLength/(EPtrimData.displayLength+1))/EPtrimData.displayLength)*(EPtrimData.EPdata.events{1}(eventList(iEvent)).sample-EPtrimData.displayPoints(1)+1) EPmain.scrsz(4)-100 50 20]);
            end;
            if strcmp(EPtrimData.EPdata.events{1}(eventList(iEvent)).value,'boundary')
                theColor='red';
            else
                theColor='black';
            end;
            
            EPtrimData.handles.eventLines(iEvent)=line([EPtrimData.EPdata.events{1}(eventList(iEvent)).sample-EPtrimData.displayPoints(1)+1 EPtrimData.EPdata.events{1}(eventList(iEvent)).sample-EPtrimData.displayPoints(1)+1],[EPtrimData.displayLength -length(EPtrimData.displayChans)*EPtrimData.spacing],'Color',theColor,'LineWidth',1);
        end;
    end;
end;