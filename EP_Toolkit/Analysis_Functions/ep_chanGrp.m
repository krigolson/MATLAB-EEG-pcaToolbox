function ep_chanGrp(EPdata,factorData)
% ep_chanGrp(EPdata,factorData) -
% Defines channel groups for windowing data.
% Strips EPdata of adds (combined channels, cells, subjects, and factors) first.
%
%Input:
%  EPdata         : Structured array with the data and accompanying information in EP file format.  See readData.
%  factorData       : Structured array with information about PCA datasets.
%       .name       : Name of the dataset.
%       .FacPat     : Factor pattern matrix.
%       .facNames   : Labels for each of the factors.
%
%Output: (via global variable)
%   EPchanGrp       : Structured array with information about channel groups
%       .

%History
%  by Joseph Dien (6/13/09)
%  jdien07@mac.com
%
%  bugfix 9/5/09 JD
%  Location of windows appearing partly off screen on some systems due to offset from having Dock (on Mac) or Task Bar
%  (on PC) at the bottom of the screen.
%
%  bugfix 1/17/10 JD
%  Workaround for Matlab bug that appears to prevent buttons from being colored when pressed on Macs using version 2009b
%  under at least some conditions.
%
%  bugfix 2/27/10 JD
%  Applied workaround to PCA guided option too.
%
% bugfix 7/2/10 JD
% In order fix crash due to Matlab bug on Windows, put in brief pause after spec name requestor diaglog box.
%
% bugfix 1/22/11 JD
% Fixed channel controls not working on Matlab versions prior to 7.8.
% Fixed some controls going off screen on a laptop screen.
%
% bugfix 8/14/11 JD
% Fixed crash when loading or saving channel group files where there is a space in the name or the path.
%
% bugfix 9/9/12 JD
% Fixed channel numbers appearing black on black background on windows computers.
%
% bugfix 2/5/13 JD
% Fixed crash when all channels are along midline or center line.
%
% bugfix 3/28/13 JD
% Fixed crash when applying factor loadings and there is more than one PCA dataset in the working set.
%
% modfied 4/1/13 JD
% Ensures there cannot be multiple channgle groups loaded with the same name.
% The active channel group in the Window pane is set to the last one being worked on in the channel group when the Done
% button is pressed.
% When saving name of channel group, change it that of the new name.
% Cancel only drops new changes rather than entirely erasing the channel group list.
%
% bugfix 5/23/13 JD
% Fixed crash when loading in a channel group file and there currently isn't any channel group defined.
% Fixed crash when going into channel group subpane, cancel out without defining a channel group, and then return to the
% subpane.
% Fixed crash during windowing when factors were used to define the areas of a channel group.
%
% modified 9/16/13 JD
% Modified display to work with 640 pixels high screen.
%
% bugfix 9/17/13 JD
% Fixed not making factor data available for windowing when first in analysis set had different number of channels.
% Fixed not asking for channel group name in windowing function when creating a new one via factor loadings and is currently the initial blank one.
% Fixed not updating the display of the area names when using factor loadings to define the channel groups.
% Fixed when using factor loadings to define channel groups, area names not reflecting factors defining them.
% Fixed when using factor loadings to define channel groups, wrong channels could be chosen.
%
%  bugfix 11/1/13 JD
%  Fixes font sizes on Windows.
%
%  modified 12/23/13 JD
%  Added windowing function option to window 'all' the channels.
%
%  bugfix 1/12/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
%  modified 5/12/14 JD
%  When using factors to set areas, can now specify whether to use largest absolute, negative, or positive loadings.
%  Will no longer keep resetting the factor loading threshold back to the original number.
%
%  modified 5/25/14 JD
%  Set colormap to jet even for Matlab 2014b onwards.
%
%  Copyright (C) 1999-2018  Joseph Dien
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

global EPchanGrp EPmain

scrsz = EPmain.scrsz;
if scrsz(4) > 900
    windowHeight=900;
else
    windowHeight=scrsz(4)-40;
end;

chanGrpFigure=findobj('Name', 'Channel Group Window');

if nargin > 0
    [EPdata]=ep_stripAdds(EPdata);
    if isempty(EPdata.data)
        msg{1}='Error: The file had no data left after additions were removed.';
        [msg]=ep_errorMsg(msg);
        return
    end;

    if ~isempty(EPchanGrp)  
        eloc=EPdata.eloc;
        
        if length(eloc) ~= length(EPchanGrp.eloc)
            EPchanGrp =[]; %if the new dataset has a different number of channels than the old one, reinitialize this window.
        else
            EPchanGrp.factorData=[];
            for i=1:length(factorData)
                if length(eloc) ==size(factorData(i).FacPat,1)
                    if isempty(EPchanGrp.factorData)
                        EPchanGrp.factorData=factorData(i);
                    else
                        EPchanGrp.factorData(end+1)=factorData(i);
                    end;
                end;
            end;
            EPchanGrp.chanNames=EPdata.chanNames;
            EPchanGrp.eloc=eloc;
        end;
    end;
end;

if isempty(EPchanGrp)    
    if ~isempty(chanGrpFigure)
        close(chanGrpFigure)
    end;
    EPchanGrp.handles.window = figure('Name', 'Channel Group Window', 'NumberTitle', 'off', 'Position',[201 scrsz(4)-windowHeight 700 windowHeight]);
    colormap jet;
    
    EPchanGrp.activeArea=1;
    EPchanGrp.activeGrp=1;
    
    EPchanGrp.factorData=[];
    for i=1:length(factorData)
        if length(EPdata.eloc) ==size(factorData(i).FacPat,1)
            if isempty(EPchanGrp.factorData)
                EPchanGrp.factorData=factorData(i);
            else
                EPchanGrp.factorData(end+1)=factorData(i);
            end;
        end;
    end;
    EPchanGrp.chanNames=EPdata.chanNames;
    EPchanGrp.eloc=EPdata.eloc;
    
    EPchanGrp.numAreas=12;
    EPchanGrp.leftMargin=.15;
    EPchanGrp.rightMargin=.1;    
    EPchanGrp.topMargin=.1;    
    EPchanGrp.bottomMargin=.5; 
    EPchanGrp.work.group(EPchanGrp.activeGrp).channel(1:length(EPchanGrp.chanNames))=EPchanGrp.numAreas+1;
    EPchanGrp.areaColors=[0 0 1;0 1 0;1 0 0;0 1 1;1 0 1;1 1 0;0 0 .5;0 .5 0;.5 0 0;0 .5 .5;.5 0 .5;.5 .5 0;0 0 0];
    
%     grpName = inputdlg('Name of new electrode grouping?');
%     pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
%     if iscell(grpName) && (length(grpName) == 1)
%         grpName=grpName{1};
%     end;
%     if isempty(grpName)
%         close(EPchanGrp.handles.window);
%         EPchanGrp=[];
%         return
%     end;
    grpName=[];
    EPchanGrp.work.group(EPchanGrp.activeGrp).name=char(grpName);
    EPchanGrp.work.group(EPchanGrp.activeGrp).threshold=[];
    for area =1:12
        EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{area}=['area' num2str(area)];
    end;
    EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{13}='none';
    EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac=0;
    EPchanGrp.work.group(EPchanGrp.activeGrp).threshold=.6;
    EPchanGrp.work.group(EPchanGrp.activeGrp).factorSelection=[];
    EPchanGrp.work.group(EPchanGrp.activeGrp).factorData=[];
    
    EPchanGrp.group=EPchanGrp.work.group;
else
    if isempty(chanGrpFigure)
        EPchanGrp.handles.window = figure('Name', 'Channel Group Window', 'NumberTitle', 'off', 'Position',[201 windowHeight 700 windowHeight]);
        EPchanGrp.work.group=EPchanGrp.group;
        if EPmain.window.chanGrp > length(EPchanGrp.group)
            EPchanGrp.activeGrp=1;
        else
            EPchanGrp.activeGrp=EPmain.window.chanGrp;
        end;
    else
        figure(EPchanGrp.handles.window);
    end;
end;

numChans=length(EPchanGrp.chanNames);

%compute screen locations of channels
EPchanGrp.X=zeros(numChans,1);
EPchanGrp.Y=zeros(numChans,1);
hasLoc=zeros(numChans,1);
% for chan=1:numChans    
%     theta=((EPchanGrp.eloc(chan).theta+90)/180)*pi;
%     radius=EPchanGrp.eloc(chan).radius;
%     if isempty(theta)
%         EPchanGrp.X(chan)=[];
%         EPchanGrp.Y(chan)=[];
%     else
%         [EPchanGrp.X(chan),EPchanGrp.Y(chan)] = pol2cart(theta,radius);
%         hasLoc(chan)=1;
%     end;    
% end;

if ~isempty(EPchanGrp.eloc) %if there are electrode coordinates
    nonLoc=0;
    for chan=1:numChans
        theta=((EPchanGrp.eloc(chan).theta+90)/180)*pi;
        radius=EPchanGrp.eloc(chan).radius;
        if isempty(theta)
            nonLoc=nonLoc+1;
            EPchanGrp.X(chan)=.02;
            EPchanGrp.Y(chan)=.02+((.05)*nonLoc);
        else
            [EPchanGrp.X(chan),EPchanGrp.Y(chan)] = pol2cart(theta,radius);
            hasLoc(chan)=1;
        end;
    end;
else %if there are no electrode coordinates
    graphSize=ceil(sqrt(numChans));
    for chan=1:numChans
        EPchanGrp.X(chan)=1-(mod(chan-1,graphSize)+1)*(.05);
        EPchanGrp.Y(chan)=1-(floor((chan-1)/graphSize)+1)*(.05);
        hasLoc(chan)=1;
    end;
end;

EPchanGrp.X=(1-EPchanGrp.X)/2;
EPchanGrp.Y=(EPchanGrp.Y+1)/2;

if any(diff(EPchanGrp.X)) %don't adjust if all along the midline
    EPchanGrp.X=((EPchanGrp.X-min(EPchanGrp.X))/(max(EPchanGrp.X-min(EPchanGrp.X)))*(1-EPchanGrp.leftMargin-EPchanGrp.rightMargin))+EPchanGrp.leftMargin;
end;
if any(diff(EPchanGrp.Y)) %don't adjust if all along the center line
    EPchanGrp.Y=((EPchanGrp.Y-min(EPchanGrp.Y))/(max(EPchanGrp.Y-min(EPchanGrp.Y)))*(1-EPchanGrp.topMargin-EPchanGrp.bottomMargin))+EPchanGrp.bottomMargin;
end;

EPchanGrp.handles.chans=zeros(numChans,1);

MATLABver=ver('MATLAB');
[a b]=strtok(MATLABver.Version,'.');
b=b(2:end);

if (str2num(a) < 7) || ((str2num(a) == 7) && (str2num(b) > 8)) && ismac %workaround for Matlab bug
    
    for chan=1:numChans
        if hasLoc(chan)
            EPchanGrp.handles.chans(chan)=uicontrol('units','normalized','position',[EPchanGrp.X(chan) EPchanGrp.Y(chan) .06 .02],'FontSize',EPmain.fontsize,...
                'Style', 'pushbutton', 'String', EPchanGrp.chanNames(chan),'ForegroundColor',EPchanGrp.areaColors(EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan),:),...
                'BackgroundColor',[1 1 1],'Callback', {@clickChannel,chan});
        end;
    end;
    
else
    
    
    for chan=1:numChans
        if hasLoc(chan)
            EPchanGrp.handles.chans(chan)=uicontrol('units','normalized','position',[EPchanGrp.X(chan) EPchanGrp.Y(chan) .06 .02],'FontSize',EPmain.fontsize,...
                'Style', 'pushbutton', 'String', EPchanGrp.chanNames(chan),'ForegroundColor',EPchanGrp.areaColors(EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan),:),...
                'Callback', {@clickChannel,chan});
        end;
    end;
    
end;

if scrsz(4) > 900
    for area =1:EPchanGrp.numAreas+1
        EPchanGrp.handles.group(EPchanGrp.activeGrp).areaName(area) = uicontrol('Style','edit','ForegroundColor', EPchanGrp.areaColors(area,:),'FontSize',EPmain.fontsize,...
            'String',EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{area},'Position',[20 windowHeight-500-area*25 100 25],...
            'Callback', ['global EPchanGrp;','EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{' num2str(area) '}=get(EPchanGrp.handles.group(EPchanGrp.activeGrp).areaName(' num2str(area) '),''String'');']);
        EPchanGrp.handles.area.channel(area)=uicontrol('position',[130 windowHeight-500-area*25 20 20],'FontSize',EPmain.fontsize,...
            'Min',0,'Max',1,...
            'Style', 'radiobutton', 'ForegroundColor', EPchanGrp.areaColors(area,:),...
            'Callback', ['global EPchanGrp;','set(EPchanGrp.handles.area.channel(EPchanGrp.activeArea),''Value'',0);','EPchanGrp.activeArea=' num2str(area) ';','set(EPchanGrp.handles.area.channel(EPchanGrp.activeArea),''Value'',1);']);
    end;
else
    for area =1:EPchanGrp.numAreas+1
        EPchanGrp.handles.group(EPchanGrp.activeGrp).areaName(area) = uicontrol('Style','edit','ForegroundColor', EPchanGrp.areaColors(area,:),'FontSize',EPmain.fontsize,...
            'String',EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{area},'Position',[20+floor((area-1)/4)*150 windowHeight-500-(rem((area-1),4)+1)*25 100 25],...
            'Callback', ['global EPchanGrp;','EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{' num2str(area) '}=get(EPchanGrp.handles.group(EPchanGrp.activeGrp).areaName(' num2str(area) '),''String'');']);
        EPchanGrp.handles.area.channel(area)=uicontrol('position',[130+floor((area-1)/4)*150 windowHeight-500-(rem((area-1),4)+1)*25 20 20],'FontSize',EPmain.fontsize,...
            'Min',0,'Max',1,...
            'Style', 'radiobutton', 'ForegroundColor', EPchanGrp.areaColors(area,:),...
            'Callback', ['global EPchanGrp;','set(EPchanGrp.handles.area.channel(EPchanGrp.activeArea),''Value'',0);','EPchanGrp.activeArea=' num2str(area) ';','set(EPchanGrp.handles.area.channel(EPchanGrp.activeArea),''Value'',1);']);
    end;
end;

set(EPchanGrp.handles.area.channel(EPchanGrp.activeArea),'Value',1);


set(EPchanGrp.handles.group(EPchanGrp.activeGrp).areaName(13),'Style','text');

EPchanGrp.handles.new = uicontrol('Style', 'pushbutton', 'String', 'New','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-110 80 35], 'Callback', @newGroup);

EPchanGrp.handles.delete = uicontrol('Style', 'pushbutton', 'String', 'Delete','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-150 80 35], 'Callback', @deleteGroup);

EPchanGrp.handles.load = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-190 80 35], 'Callback', @loadGroup);

EPchanGrp.handles.save = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-230 80 35], 'Callback', @saveGroup);

EPchanGrp.handles.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-270 80 35], 'Callback', ['global EPchanGrp; close(EPchanGrp.handles.window); EPchanGrp.work.group=EPchanGrp.group; ep(''start'');']);

EPchanGrp.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-310 80 35], 'Callback', ['global EPchanGrp EPmain; EPchanGrp.group=EPchanGrp.work.group; EPmain.window.chanGrp=EPchanGrp.activeGrp; close(EPchanGrp.handles.window); ep(''start''); ']);

h = uicontrol('Style','text',...
    'String','Channel Group','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[105 windowHeight-95 90 20]);

EPchanGrp.handles.chanGrp = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',{EPchanGrp.work.group.name},...
    'Value',EPchanGrp.activeGrp,'Position',[100 windowHeight-115 100 20], 'Callback', ['global EPchanGrp;','EPchanGrp.activeGrp=get(EPchanGrp.handles.chanGrp,''Value''); ep_chanGrp;']);



EPchanGrp.handles.factor = uicontrol('Style', 'pushbutton', 'String', 'Factor','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-350 80 35], 'Callback', ['global EPchanGrp; EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac=1; ep_chanGrp; ']);

if isempty(EPchanGrp.factorData)
    set(EPchanGrp.handles.factor,'enable','off');
end;


if EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac
    h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
        'String','Factor File','HorizontalAlignment','left',...
        'Position',[20 windowHeight-370 90 20]);
    
    EPchanGrp.handles.facFile = uicontrol('Style','popupmenu',...
        'String',{EPchanGrp.factorData.name},'FontSize',EPmain.fontsize,...
        'Value',EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac,'Position',[20 windowHeight-390 100 20], 'Callback', ['global EPchanGrp;','EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac=get(EPchanGrp.handles.facFile,''Value'');']);

    h = uicontrol('Style','text',...
        'String','Threshold','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'Position',[20 windowHeight-410 90 20]);
    
    EPchanGrp.handles.threshold = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
        'String',EPchanGrp.work.group(EPchanGrp.activeGrp).threshold,'Position',[20 windowHeight-430 100 20], 'Callback', ['global EPchanGrp;','EPchanGrp.work.group(EPchanGrp.activeGrp).threshold=str2num(get(EPchanGrp.handles.threshold,''String''));']);
    
    EPchanGrp.handles.minusPlus = uicontrol('Style', 'pushbutton', 'String', '-/+','FontSize',EPmain.fontsize*1.5,...
        'Position', [20 windowHeight-470 30 35], 'Callback', @applyFactor);

    EPchanGrp.handles.plus = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*1.5,...
        'Position', [55 windowHeight-470 30 35], 'Callback', @applyFactor);

    EPchanGrp.handles.minus = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*1.5,...
        'Position', [90 windowHeight-470 30 35], 'Callback', @applyFactor);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadGroup(src,eventdata) %load an electrode group
global EPchanGrp

[FileName,PathName,FilterIndex] = uigetfile('','Load Channel Group','Group1');
if FileName ~= 0
    eval(['load ''' PathName FileName '''']);
    
    %backward compatibility checks
    if ~isfield(EPgroup,'factorSelection')
        EPgroup.factorSelection=[];
    end;
    if ~isfield(EPgroup,'factorData')
        EPgroup.factorData=[];
    end;
    
    if length(EPgroup.channel) == length(EPchanGrp.eloc)
        if ~isempty(EPchanGrp.work.group(EPchanGrp.activeGrp).name)
            sameName=1;
            suffix=0;
            grpNameSuffix=EPgroup.name;
            while sameName
                sameName=0;
                for i=1:length(EPchanGrp.work.group)
                    if strcmp(EPchanGrp.work.group(i).name,grpNameSuffix)
                        sameName=1;
                    end;
                end;
                if sameName
                    suffix=suffix+1;
                    grpNameSuffix=[EPgroup.name '-' num2str(suffix)];
                end;
            end;
            EPgroup.name=grpNameSuffix;
            EPchanGrp.work.group(end+1)=EPgroup;
            EPchanGrp.activeGrp=length(EPchanGrp.work.group);
        else
            EPchanGrp.work.group(EPchanGrp.activeGrp)=EPgroup;
        end;
    else
        warndlg('Number of electrodes different from current dataset.');
    end;
    ep_chanGrp;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveGroup(src,eventdata) %save an electrode group
global EPchanGrp

[fileName,PathName,FilterIndex] = uiputfile('*.mat','Save Channel Group',EPchanGrp.work.group(EPchanGrp.activeGrp).name);
[pathstr, grpName, ext] = fileparts(fileName);

EPchanGrp.work.group(EPchanGrp.activeGrp).name=grpName;

EPgroup=EPchanGrp.work.group(EPchanGrp.activeGrp);
eval(['save ''' PathName fileName ''' EPgroup']);
ep_chanGrp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newGroup(src,eventdata) %Add a new electrode group
global EPchanGrp

grpName = inputdlg('Name of new electrode grouping?');

sameName=1;
suffix=0;
grpNameSuffix=grpName{1};
while sameName
    sameName=0;
    for i=1:length(EPchanGrp.work.group)
        if strcmp(EPchanGrp.work.group(i).name,grpNameSuffix)
            sameName=1;
        end;
    end;
    if sameName
        suffix=suffix+1;
        grpNameSuffix=[grpName{1} '-' num2str(suffix)];
    end;
end;

pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
if ~isempty(EPchanGrp.work.group(1).name)
    EPchanGrp.work.group(end+1).name=char(grpNameSuffix);
else
    EPchanGrp.work.group(1).name=char(grpNameSuffix);
end;
EPchanGrp.activeGrp=length(EPchanGrp.work.group);
EPchanGrp.work.group(EPchanGrp.activeGrp).channel(1:length(EPchanGrp.eloc))=EPchanGrp.numAreas+1;
for area =1:EPchanGrp.numAreas
    EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{area}=['area' num2str(area)];
end;
EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{EPchanGrp.numAreas+1}='none';
EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac=0;
EPchanGrp.work.group(EPchanGrp.activeGrp).threshold=.6;

clf(EPchanGrp.handles.window)
ep_chanGrp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clickChannel(src,eventdata,chan) %Click on electrode
global EPchanGrp

if EPchanGrp.activeArea == EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan)
    EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan)=EPchanGrp.numAreas+1;
else
    EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan)=EPchanGrp.activeArea;
end;

set(EPchanGrp.handles.chans(chan),'ForegroundColor',EPchanGrp.areaColors(EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan),:));

if isempty(EPchanGrp.work.group(EPchanGrp.activeGrp).name)
    grpName = inputdlg('Name of new electrode grouping?');
    pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
    if iscell(grpName) && (length(grpName) == 1)
        grpName=grpName{1};
    end;
    if isempty(grpName)
        close(EPchanGrp.handles.window);
        EPchanGrp=[];
        return
    end;
    EPchanGrp.work.group(EPchanGrp.activeGrp).name=char(grpName);
    ep_chanGrp;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function applyFactor(src,eventdata) %Use factor results to form channel areas
global EPchanGrp

switch src
    case EPchanGrp.handles.minus
        polarity='minus';
    case EPchanGrp.handles.plus
        polarity='plus';
    case EPchanGrp.handles.minusPlus
        polarity='minusPlus';
    otherwise
        disp('Programmer error');
        return
end;

numChans=length(EPchanGrp.eloc);
threshold=abs(str2num(get(EPchanGrp.handles.threshold,'string')));

%select factors to apply
[factorSelection,ok] = listdlg('PromptString',['Choose up to' num2str(EPchanGrp.numAreas) ' factors'],'ListString',EPchanGrp.factorData(EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac).facNames);
numSelectedFacs=length(factorSelection);
if ok && ~isempty(factorSelection) && numSelectedFacs <=EPchanGrp.numAreas
    EPchanGrp.work.group(EPchanGrp.activeGrp).factorSelection=factorSelection;
    EPchanGrp.work.group(EPchanGrp.activeGrp).factorData=EPchanGrp.factorData(EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac);
    for area =1:numSelectedFacs
        EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{area}=EPchanGrp.factorData(1).facNames{factorSelection(area)};
    end;
    %assign each channel to factor with loading that is higher than the threshold
    facSigns=zeros(numChans,numSelectedFacs);
    for chan=1:numChans
        switch polarity
            case 'minus'
                [theLoading biggestFac]= min(EPchanGrp.work.group(EPchanGrp.activeGrp).factorData(1).FacPat(chan,factorSelection));
                if theLoading < 0
                    theLoading=abs(theLoading);
                else
                    theLoading =0;
                end;
            case 'plus'
                [theLoading biggestFac]= max(EPchanGrp.work.group(EPchanGrp.activeGrp).factorData(1).FacPat(chan,factorSelection));
                if theLoading < 0
                    theLoading = 0;
                end;
            case 'minusPlus'
                [theLoading biggestFac]= max(abs(EPchanGrp.work.group(EPchanGrp.activeGrp).factorData(1).FacPat(chan,factorSelection)));
        end;
        if theLoading >= threshold
            facSigns(chan,biggestFac)=sign(EPchanGrp.work.group(EPchanGrp.activeGrp).factorData(1).FacPat(chan,factorSelection(biggestFac)));
            EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan)=biggestFac;
        end;
    end;
    if strcmp(polarity,'minusPlus')
        %if a factor has both positive and negative channels then keep only channels with the same sign factor loading
        for factor =1:numSelectedFacs
            if any(facSigns(:,factor)>0) && any(facSigns(:,factor)<0)
                if max(EPchanGrp.work.group(EPchanGrp.activeGrp).factorData(1).FacPat(find(facSigns(:,factor)>0),factorSelection(factor))) >= max(abs(EPchanGrp.work.group(EPchanGrp.activeGrp).factorData(1).FacPat(find(facSigns(:,factor)<0),factorSelection(factor))))
                    EPchanGrp.work.group(EPchanGrp.activeGrp).channel(find(facSigns(:,factor)<0))=EPchanGrp.numAreas+1;
                else
                    EPchanGrp.work.group(EPchanGrp.activeGrp).channel(find(facSigns(:,factor)>0))=EPchanGrp.numAreas+1;
                end;
            end;
        end;
    end;
    for chan=1:numChans
        set(EPchanGrp.handles.chans(chan),'ForegroundColor',EPchanGrp.areaColors(EPchanGrp.work.group(EPchanGrp.activeGrp).channel(chan),:));
    end;
    drawnow
end;

if isempty(EPchanGrp.work.group(EPchanGrp.activeGrp).name)
    grpName = inputdlg('Name of new electrode grouping?');
    pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
    if iscell(grpName) && (length(grpName) == 1)
        grpName=grpName{1};
    end;
    if isempty(grpName)
        close(EPchanGrp.handles.window);
        EPchanGrp=[];
        return
    end;
    EPchanGrp.work.group(EPchanGrp.activeGrp).name=char(grpName);
    ep_chanGrp;
end;

ep_chanGrp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deleteGroup(src,eventdata) %Delete current electrode group
global EPchanGrp

numGroups=length(EPchanGrp.work.group);
if  numGroups == 1 %if only one group, then need to replace it with a new one
    grpName = inputdlg('Name of new electrode grouping?');
    pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
    if ~isempty(grpName)
        EPchanGrp.work.group(1).name=char(grpName);
        EPchanGrp.activeGrp=1;
        EPchanGrp.work.group(EPchanGrp.activeGrp).channel(1:length(EPchanGrp.eloc))=EPchanGrp.numAreas+1;
        for area =1:EPchanGrp.numAreas
            EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{area}=['area' num2str(area)];
        end;
        EPchanGrp.work.group(EPchanGrp.activeGrp).areaName{EPchanGrp.numAreas+1}='none';
        EPchanGrp.work.group(EPchanGrp.activeGrp).activeFac=0;
        EPchanGrp.work.group(EPchanGrp.activeGrp).threshold=.6;
    else
        return
    end;
else
    EPchanGrp.work.group=EPchanGrp.work.group(setdiff([1:numGroups],EPchanGrp.activeGrp));
    EPchanGrp.activeGrp=max(1,EPchanGrp.activeGrp-1);
end;

clf(EPchanGrp.handles.window)
ep_chanGrp;




