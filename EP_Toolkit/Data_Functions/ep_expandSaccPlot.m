function ep_expandSaccPlot(src,eventdata)
% ep_expandSaccPlot - ep_expandEyePlot -
% Expands a saccade plot into a full window for scan waves function.
%

%History
%  by Joseph Dien (5/20/16)
%  jdien07@mac.com
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

global EPscanCont EPmain

if isempty(EPmain.view)
    warndlg('When you leave the View pane, the data linked to waveform plots are no longer available.  You''ll need to generate a new waveform plot.');
    return
end;

scrsz = EPmain.scrsz;

EPscanCont.handles.waves.ExpandedSaccFigure = figure('Name', 'Saccade-Tracker Plot', 'NumberTitle', 'off', 'Position',[scrsz(3)/2 scrsz(4)/2 800 800]);
colormap jet;

colorData=zeros(length(EPscanCont.displayPoints),3);
if EPscanCont.chanShow~=length(EPscanCont.chanList)
    theData=squeeze(EPscanCont.EPdata(1).data(EPscanCont.chanShow,EPscanCont.displayPoints,:,:,:,EPscanCont.freqPos))';
    if strcmp(EPscanCont.dataType,'TFT')
        theData=abs(theData); %convert complex number to real number
        theData=theData/mean(diff(EPscanCont.EPdata(1).freqNames)); %convert to spectral density
        theData=theData.^2; %convert amplitude to power
        theData=log10(abs(theData))*10; %convert to dB log scaling
        tempVar=theData;
        tempVar(isinf(tempVar))=-flintmax;
        theData=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
    end;
    colorData(:,1)=((theData-mean(theData))/max(theData-mean(theData)))>.5;
    %colorData(:,3)=(theData-mean(theData))/min(theData-mean(theData));
    %colorData(colorData<0)=0;
end;

stimPict=[];
if ~isempty(EPscanCont.stimList)
    whichPict=max(find([EPscanCont.stimList.sample] <= EPscanCont.displayPoints(EPscanCont.plotPoint)));
    if ~isempty(whichPict)
        stimPict=EPscanCont.EPdata(1).stims(EPscanCont.stimList(whichPict).stim).image;
    end;
end;

Xdata=squeeze(EPscanCont.EPdata(1).data(EPscanCont.HSACchan,EPscanCont.displayPoints,:,:,:,1));
Ydata=squeeze(EPscanCont.EPdata(1).data(EPscanCont.VSACchan,EPscanCont.displayPoints,:,:,:,1));
Xdata=(-(Xdata-EPscanCont.sacCenterX)/EPscanCont.sacScaleX)+.5;
Ydata=(-(Ydata-EPscanCont.sacCenterY)/EPscanCont.sacScaleY)+.5;
EPscanCont.handles.topos.expandSacc=axes(EPscanCont.handles.waves.ExpandedSaccFigure);
axis([0 1 0 1])
set(EPscanCont.handles.topos.expandSacc,'Units','normalized')
set(EPscanCont.handles.topos.expandSacc,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
hold on
if ~isempty(stimPict)
    imageSize=size(stimPict);
    imageMax=max(imageSize(1),imageSize(2));
    x=max(1,round((imageSize(1)-imageSize(2))/2))/imageMax;
    y=max(1,round((imageSize(2)-imageSize(1))/2))/imageMax;
    EPscanCont.handles.topos.saccImage=image([x 1-x],[1-y y],stimPict);
end;
if EPscanCont.chanShow~=length(EPscanCont.chanList)
    for iPoint=1:length(EPscanCont.displayPoints)-1
%        if any(colorData(iPoint,:)>.5)
            %plot(Xdata(iPoint:iPoint+1),Ydata(iPoint:iPoint+1),'color',colorData(iPoint,:));
            scatter(Xdata(iPoint),Ydata(iPoint),30,'filled','MarkerFaceAlpha',1/8,'MarkerFaceColor',colorData(iPoint,:));
%        end;
    end;
else
    plot(Xdata,Ydata,'color','black');
end;
plot(Xdata(EPscanCont.plotPoint),Ydata(EPscanCont.plotPoint),'color','green','Marker','*','MarkerSize',2);
hold off
