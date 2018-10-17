function ep_showTopos(startTime,endTime,startHz,endHz,marker1,marker2,FFTunits,minVolt,maxVolt,eventLines)
% ep_showTopos - ep_showTopos(startTime,endTime,startHz,endHz,FFTunits,minVolt,maxVolt) -
% Displays the scalp topographies of the data and performs source analyses.
% Will choose the maximum channel, time point, and herz for display.  If multiple datasets are being displayed,
% will use the one in the leftmost column not set to "none" as the index dataset and only judge based on columns
% corresponding to that one dataset.  That way there is no problem of apples and orange comparison for max values,
% as in correlations for coherence versus power for FFT.
% Regional channels are not included in the peak calculations.
%
%Input:
%  startTime         : Start of epoch
%  endTime           : End of epoch in ms (right side unless both startTime and endTime are the same number)
%  startHz           : Start of Hz band
%  endHz             : End of Hz band
%  marker1           : Optional marker for plots
%  marker2           : Second optional marker for plots
%  FFTunits          : spectral data units (1=complex, 3=asd, 3=psd, 4=dB)
%  minVolt           : Minimum data value in raw units.
%  maxVolt           : Maximum data value in raw units.
%  eventLines        : Cell array of samples of events to be displayed (4 colors)(cells/subs)
%
%  Time is counted for an epoch as, say, -200 to 800 ms.  In this case the samples have some length to them depending on
% the digitization rate, such as 4ms for a 250Hz sampling rate.  The points prior to the baseline are counted from the
% left side of the baseline sample and the points following the baseline are counted from the right side of the baseline
% sample.  Thus, when the actual samples are calculated, one first adds the baseline to the numbers, yielding 0 to
% 1000ms.  One then must add the length of the baseline to the first number to obtain 4 to 1000ms.  Then one would
% divide the numbers by the sample size, resulting in 1 to 250 samples, which are the correct samples to use.
% For continuous data or other such data where there is no baseline, one must by this convention start with 0 ms,
% as in 0-1000ms and 1000-2000ms.

%History
%  by Joseph Dien (4/20/10)
%  jdien07@mac.com

% bugfix 5/9/10 JD
% Fixed crash when using Topo button of View EEG pane and not all four colors are being used.
% Fixed can only change channel and latency settings for number of rows equal to number of pages of factors.
% Fixed sometimes crashes when first color is set to none.
%
% modified 5/15/10 JD
% Added white marker to topos for electrode corresponding to the waveform figures.
% Small black dots indicate electrode locations in topographical plots.
% May click on electrode dots to move the waveform plot channel.
% May right-click on topographical plot to obtain expanded 2D plot.
% May right-click on topographical plot to obtain expanded 3D plot.
% May right-click on topographical plot to obtain basic dipole analysis.
% Eliminated sorting of channel names.
%
% bugfix 6/6/10 JD
% Fixed crash when trying to display 3D plot or dipole source using .ced file generated from a .elp file.
%
% modified 6/15/10 JD
% Added dipole analysis of jack-knifed PCA results to the Topos function of the View Pane.
% Marks when a file isn't saved yet.
%
% bugfix 7/23/10 JD
% Fixed crash when doing dipole or 3D function and ced is either empty or "none".
% Fixed 3D, dipole, and jackknife results incorrect for topos not on the first page (when there are multiple pages of
% topos).
%
% bugfix 12/26/10 JD
% Fixed bottom of Topos window being cut off on laptop screens.
%
% bugfix 1/9/11 JD
% Fixed channels not being unstandardized correctly in jack-knife, resulting in some inaccuracy in dipole results.
% Fixed crash when trying to display data where some channels are missing electrode coordinates.
%
% bugfix 2/12/11 JD
% Fixed contents of topos window getting shifted upwards off window when OS X Dock is at bottom of screen.
%
% bugfix 2/14/11 JD
% Fixed crash when displaying data with only one timepoint.
%
% bugfix 3/24/11 JD
% Fixed crash when displaying data with only two timepoints.
%
% modified 3/15/12 JD
% Added support for plotting FFT data, including psd and dB scaling.
%
% bugfix 6/3/12 JD
% Fixed crash when conducting dipole analysis on data with regional channels.
%
% bugfix 6/4/12 JD
% Fixed crash after clicking on Done when the final colors are not PCA data.
%
% bugfix 9/1/12 JD
% Fixed minimum voltage scaling being set to -1000 whenever minimum voltages are below 1000.
%
% bugfix 11/4/12 JD
% Fixed crash when invoking 2D plots.
%
% modified 1/16/13 JD
% Handles situation where FFT data has negative values (e.g., due to imprecision in PCA results) and transforming to dB
% would result in complex numbers, by taking absolute value first.
%
% modified 1/24/13 JD
% Added option to contextual menu to rescale figures according to selected topo map.
% Changed peak Hz calculation for FFT data to operate on absolute values to accommodate sometimes negative numbers from
% PCA results.
%
% modified 1/25/13 JD
% Added markers and expanding window to waveform figures.
%
% bugfix 1/31/13 JD
% Fixed a variety of issues with the automatic scaling in the topomap figures, especially for spectral data.
%
% bugfix 2/6/13 JD
% Fixed peak channels not being identified correctly.
%
% bugfix 2/19/13 JD
% Fixed crash when displaying data with a regional channel.
%
% bugfix 2/25/13 JD
% Fixed crash when displaying frequency data in dB scaling and the maximum value is negative.
%
% bugfix 3/24/13 JD
% Fixed peak samples of ERP data not being identified by absolute amplitude.
% Fixed Topos not allowing two datasets to be shown in parallel when they have different regional channels.
% Fixed Topos crashing when trying to display factor and non-factor data side by side.
%
% modified 3/24/13 JD
% Added display of topos at every 50 ms for ERP and TFT data and every Hz for FFT data for non-factor data.
%
% modified 4/2/13 JD
% Added peak point/Hz line to expanded waveform windows.
%
% bugfix 4/2/13 JD
% Markers in waveform plots can be set at zero ms.
%
% bugfix 4/12/13 JD
% Scaling of topos now obeys the values on the View pane.  Also, fixed manual changes to plotting range no longer
% working.
%
% modified 9/26/13 JD
% Restricted peak channels to EEG channels.
%
% bugfix 10/10/13 JD
% Fixed topoplot not showing the correct peak channel for non-factor data.
%
% bugfix 11/1/13 JD
% Fixed font sizes on Windows.
%
% bugfix 11/22/13 JD
% Fixed crash when trying to change scaling of frequency-domain data.
%
% bugfix 1/12/14 JD
% Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% bugfix 1/15/14 JD
% Fixed crash when trying to plot frequency data where electrode coordinate information is present but all coordinates are missing.
%
% bugfix 2/25/14 JD
% Fixed crash when trying to plot in dB data where the power equals zero.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% modified 4/23/14 JD
% Added coherence and phase-locking options, including support for complex numbers.
% For peak channels and points and Hz, uses only the first dataset as the index for this.
%
% bufix 4/28/14 JD
% Fixed crash when changing pages in Topos view with frequency data.
% Fixed TFT color scale for expanded plots.
% Fixed 2D expanded head plots for Topos view for frequency data.
% Fixed jack-knife test possibly conducted on wrong factor or just crashing in Topos view.
%
% bufix 5/29/14 JD
% Fixed crash when there are electrodes without coordinates.
%
% modified 6/19/14 JD
% Added support for sample-by-sample t-tests, including STS chanType.
%
% bufix 8/5/14 JD
% Fixed not displaying frequency-domain data.
%
% modified 8/12/14 JD
% Added explore rereferencing option to contextual menus.
%
% bufix 12/31/14 JD
% Fixed crash with 3D head when original CED file not available.
%
% bufix 1/9/15 JD
% Fixed axis labels being shown for waveform plots under Matlab 2014b due
% to Matlab bug.
%
% bufix 3/22/15 JD
% Fixed peak chans sometimes not being computed correctly for voltage data.
%
%  modified 5/25/14 JD
%  Set colormap to jet even for Matlab 2014b onwards.
%
%  modified 5/29/14 JD
%  Changed min and max scale to be set by plotted data unless overriden.
%
%  bugfix 7/4/15 JD
%  Fixed crash when using rereference function to change the displayed
%  referencing.
%  Fixed crash when plotting two datasets with different numbers of
%  electrodes and the larger one is the first one.
%
% modified 7/5/15 JD
% Performs average reference prior to performing PARE correction.
%
% bugfix 8/30/15 JD
% Fixed amplitude spectral density calculated as divided by Hz rather than
% by square root of Hz.
% Fixed dB of amplitude data not converted to power first.
% Fixed dB units should be labeled as dBV since it is a ratio and therefore
% has no units.
%
% bugfix 1/15/16 JD
% Now allows power scaled data to be displayed as amplitudes.
% Now handles complex FFT numbers.
%
% modified 1/21/16 JD
% Consolidated spectral unit controls so just four options (cm, am, pw, dB).
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% bugfix 10/9/16 JD
% Fixed FFT data on amplitude scaling not converted to real numbers.
%
% bugfix & modified 6/20/17 JD
% Fixed displaying imaginary component of FFT data even if FFT units not specified as being complex.
% Improved auto scaling for dB unit FFT data.
% Fixed conversion to spectral density dividing by bin width rather than sqrt(bin width).
% Fixed sometimes crashing when frequency band is just one Hz.
% Now uses ep_expandChan for expanded view as well.
% Added Compare Channels option to the View Topos pop-up menu.
%
% modified & bugfix 5/23/18 JD
% Fixed sometimes one too many pages when displaying factors, resulting in a crash when one visits the erroneous one.
% Fixed color bar missing.
% Added support for listing all trials or cells or subjects.
%
% modified 6/5/18 JD
% Added option to add marks for RT and selected events to View figures.
% Added Synch checkbox to Topos view so that when the channel or the time point is changed for one row, it is changed for all of them.
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

global EPdataset EPmain EPtopos

scrsz = EPmain.scrsz;

windowHeight=scrsz(4)-100;
textColors={'blue','red',[0 .5 0],'black'};

if EPtopos.page==0 %set up topoplot for the first time
    
    EPtopos.RGBcolors=[0 0 1;1 0 0;0 1 0;0 0 0]; %colors for the waveforms
    
    EPtopos.handles.topos.topoWindow = figure('Name', 'TopoWindow', 'NumberTitle', 'off', 'Position',[201 101 800 windowHeight]);
    colormap jet;
    EPtopos.plotWidth=140;
    EPtopos.plotHeight=100;
    EPtopos.topoSize=100;
    EPtopos.margin=.1;
    
    EPtopos.perPage=min(floor(windowHeight/120)-1,10); %how many rows of figures per page
    EPtopos.page=1; %current page
    EPtopos.chans=[];
    EPtopos.points=[];
    EPtopos.synch=0;
    
    EPtopos.FFTunits=FFTunits;
    EPtopos.CSD=NaN;
    
    EPtopos.marker1=marker1;
    EPtopos.marker2=marker2;
    
    EPtopos.twoChan=[];
    EPtopos.handles.twoChan.figure=gobjects(0);
    
    EPtopos.eventLines=eventLines;
    
    %make sure the plotted datasets are compatible
    numChans=0;
    EPtopos.numRows=1;
    for color=1:4
        if EPmain.view.dataset(color) <= length(EPdataset.dataset) %if not set to "none"
            if numChans == 0
                EPtopos.flexMode=strcmp(EPdataset.dataset(EPmain.view.dataset(color)).timeUnits,'per');
                %windowName=[EPdataset.dataset(EPmain.view.dataset(color)).dataName '-' EPdataset.dataset(EPmain.view.dataset(color)).subNames{EPmain.view.subject(color)}];
                windowName=EPdataset.dataset(EPmain.view.dataset(color)).dataName;
                set(EPtopos.handles.topos.topoWindow,'Name',windowName);
                EPtopos.firstTime=min([EPdataset.dataset(EPmain.view.dataset(color)).timeNames]);
                EPtopos.lastTime=max([EPdataset.dataset(EPmain.view.dataset(color)).timeNames]);
                EPtopos.firstHz=min([EPdataset.dataset(EPmain.view.dataset(color)).freqNames]);
                EPtopos.lastHz=max([EPdataset.dataset(EPmain.view.dataset(color)).freqNames]);
                EPtopos.sampleSize=mean(diff(EPdataset.dataset(EPmain.view.dataset(color)).timeNames));
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(color)).freqNames)
                    EPtopos.binSize=mean(diff(EPdataset.dataset(EPmain.view.dataset(color)).freqNames));
                end;
                implicit=EPdataset.dataset(EPmain.view.dataset(color)).implicit; %use implicits of the first dataset
                EPtopos.chanNames=EPdataset.dataset(EPmain.view.dataset(color)).chanNames; %use channel names of the first dataset
                EPtopos.nonRegChans=find(~ismember(EPdataset.dataset(EPmain.view.dataset(color)).chanTypes,{'ANS','ECG','REG'}));
                
                EPtopos.eloc=EPdataset.dataset(EPmain.view.dataset(color)).eloc;
                if isempty(EPtopos.eloc) || all(isempty([EPtopos.eloc.radius])) || all(isempty([EPtopos.eloc.theta]))
                    msg{1}=['Error: The dataset ' EPdataset.dataset(EPmain.view.dataset(color)).dataName ' has no electrode coordinates.'];
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end;
                hasLoc=[];
                for iChan=1:length(EPtopos.nonRegChans)
                    if ~isempty(EPtopos.eloc(EPtopos.nonRegChans(iChan)).theta)
                        hasLoc(end+1)=EPtopos.nonRegChans(iChan);
                    end;
                end;  
                EPtopos.nonRegChans=hasLoc;
                
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(color)).relNames)
                    EPtopos.nonRegRelChans=EPtopos.nonRegChans;
                else
                    EPtopos.nonRegRelChans=1;
                end;
                EPtopos.ced=EPdataset.dataset(EPmain.view.dataset(color)).ced;
                EPtopos.dataName=EPdataset.dataset(EPmain.view.dataset(color)).dataName;
                if ~isempty(implicit)
                    implicit=implicit(find(~strcmp('FID',{implicit.type})));
                    EPtopos.chanNames=[EPtopos.chanNames {implicit.labels}]; %use channel names of the first dataset
                    EPtopos.eloc=[EPtopos.eloc implicit];
                end;
                numChans=length(EPtopos.nonRegChans);
                if numChans==0
                    msg{1}='Error: No channels with known locations.';
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end;
                
                if any(EPmain.view.allTrials==1)
                    if strcmp(EPdataset.dataset(EPmain.view.dataset(color)).dataType,'average')
                        EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(color)).subNames);
                        EPtopos.type='subject';
                    else %single-trial
                        EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(color)).trialNames);
                        EPtopos.type='trial';
                    end;
                elseif any(EPmain.view.allTrials==3)
                    EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(color)).facNames);
                    EPtopos.type='factor';
                elseif any(EPmain.view.allTrials==5)
                    EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(color)).cellNames);
                    EPtopos.type='cell';
                elseif any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
                    EPtopos.type='time';
                elseif strcmp(EPmain.view.dataTransform,'FFT')
                    EPtopos.type='freq';
                else
                    error('Programming Error - aborting.');
                end
                
                if isfield(EPtopos.eloc,'labels')
                    eLabels={EPtopos.eloc.labels};
                else
                    eLabels=cell(numChans,1);
                end;
                for iChan=1:length(eLabels)
                    if isempty(eLabels{iChan})
                        eLabels{iChan}=EPtopos.chanNames{iChan};
                    end;
                end;
                %                 [B,chanIX{color}] = sort(eLabels);
                %                 if ~isempty(EPtopos.eloc)
                %                     sortedEloc=EPtopos.eloc(chanIX{color});
                %                 end;
                %EPtopos.chanNames=EPtopos.chanNames(chanIX{color});
                EPtopos.theFirstColor=color; %the leftmost color that is not set to "none" and hence used as the index color
                EPtopos.rowList=[];
            else
                newImplicit=EPdataset.dataset(EPmain.view.dataset(color)).implicit;
                if ~isempty(newImplicit)
                    newImplicit=newImplicit(find(~strcmp('FID',{newImplicit.type})));
                end;
                newNonRegChans=find(~ismember(EPdataset.dataset(EPmain.view.dataset(color)).chanTypes,{'ANS','ECG','REG'}));
                hasLoc=[];
                for iChan=1:length(newNonRegChans)
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(color)).eloc(newNonRegChans(iChan)).theta)
                        hasLoc(end+1)=newNonRegChans(iChan);
                    end;
                end;  
                newNonRegChans=hasLoc;

                if length(EPtopos.nonRegChans) ~= length(newNonRegChans)
                    msg{1}='Error: The number of channels are not compatible.';
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end;
                EPtopos.firstTime =max(min([EPdataset.dataset(EPmain.view.dataset(color)).timeNames]),EPtopos.firstTime);
                EPtopos.lastTime = min(max([EPdataset.dataset(EPmain.view.dataset(color)).timeNames]),EPtopos.lastTime);
                EPtopos.firstHz =max(min([EPdataset.dataset(EPmain.view.dataset(color)).freqNames]),EPtopos.firstHz);
                EPtopos.lastHz = min(max([EPdataset.dataset(EPmain.view.dataset(color)).freqNames]),EPtopos.lastHz);
                if ~isnan(EPtopos.sampleSize) && EPtopos.sampleSize ~= mean(diff(EPdataset.dataset(EPmain.view.dataset(color)).timeNames))
                    msg{1}='Error: The sampling rates are not compatible.';
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end;
                
                switch EPtopos.type
                    case 'subject'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(color)).subNames));
                    case 'trial'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(color)).trialNames));
                    case 'factor'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(color)).facNames));
                    case 'cell'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(color)).cellNames));
                end;
                        
                newEloc=[EPdataset.dataset(EPmain.view.dataset(color)).eloc newImplicit];
                if ~isempty(EPtopos.eloc) && ~isempty(newEloc)
                    eLabels={newEloc.labels};
                    for chan=1:length(eLabels)
                        if isempty(eLabels{chan})
                            eLabels{chan}=EPdataset.dataset(EPmain.view.dataset(color)).chanNames{chan};
                        end;
                    end;
                    %[B,chanIX{color}] = sort(eLabels);
                    %newEloc=newEloc(chanIX{color});
                    if any([EPtopos.eloc.theta]-[newEloc.theta]) || any([EPtopos.eloc.radius]-[newEloc.radius])
                        msg{1}='Error: The electrode coordinates are not compatible.';
                        [msg]=ep_errorMsg(msg);
                        done
                        return
                    end;
                    %                 elseif xor(isempty(EPtopos.eloc),isempty(newEloc)) %if one but not the other has no electrode coordinates
                    %                     msg{1}='Error: The electrode coordinates are not compatible.';
                    %                     [msg]=ep_errorMsg(msg);
                    %                     done
                    %                     return
                else
                    if isempty(newEloc)
                        msg{1}=['Error: The dataset ' EPdataset.dataset(EPmain.view.dataset(color)).dataName ' has no electrode coordinates.'];
                        [msg]=ep_errorMsg(msg);
                        done
                        return
                    end;
                end;
            end;
        end;
    end;
    
    if startTime==endTime
        EPtopos.firstTime=max(EPtopos.firstTime,startTime);
        EPtopos.lastTime=min(EPtopos.lastTime,endTime); %assume already  left side of sample
    else
        EPtopos.firstTime=max(EPtopos.firstTime,startTime);
        EPtopos.lastTime=min(EPtopos.lastTime,endTime-EPtopos.sampleSize); %shift endTime to left side of sample
    end;
    EPtopos.firstHz=max(EPtopos.firstHz,startHz);
    EPtopos.lastHz=min(EPtopos.lastHz,endHz);
    
    %compute epoch samples for each color
    EPtopos.numPoints=1;
    EPtopos.numHz=1;
    startBins=ones(1,4);
    lastBins=ones(1,4);
    startSamp=ones(1,4);
    lastSamp=ones(1,4);
    
    if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
        for color=1:4
            if EPmain.view.dataset(color) <= length(EPdataset.dataset)
                [M startBins(color)]=min(abs(EPdataset.dataset(EPmain.view.dataset(color)).freqNames-EPtopos.firstHz));
                [M lastBins(color)]=min(abs(EPdataset.dataset(EPmain.view.dataset(color)).freqNames-EPtopos.lastHz));
                lastBins(color)=max(lastBins(color),startBins(color));
                EPtopos.numHz=lastBins(color)-startBins(color)+1;
            end;
        end
        EPtopos.freqNameList=cell(0);
        plotDataList=find(EPmain.view.dataset <= length(EPdataset.dataset));
        for iBin=min(startBins(plotDataList)):max(lastBins(plotDataList))
            EPtopos.freqNameList{end+1}=sprintf('%5.1f',round(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).freqNames(iBin)*10)/10);
        end;
    end;
    if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
        for color=1:4
            if EPmain.view.dataset(color) <= length(EPdataset.dataset)
                [M startSamp(color)]=min(abs(EPdataset.dataset(EPmain.view.dataset(color)).timeNames-EPtopos.firstTime));
                [M lastSamp(color)]=min(abs(EPdataset.dataset(EPmain.view.dataset(color)).timeNames-EPtopos.lastTime));
                lastSamp(color)=max(lastSamp(color),startSamp(color));
                EPtopos.numPoints=lastSamp(color)-startSamp(color)+1;
            end;
        end
    end;
    
    %organize the data for plotting
    EPtopos.colorIndex=[];
    EPtopos.plotColors=[];
    EPtopos.thePlotColors=[];
    EPtopos.plotLineIndex=cell(0);
    EPtopos.plotColorIndex=[];
    theDataset=0;
    EPtopos.totalData=zeros(numChans,EPtopos.numPoints,EPtopos.numRows,EPtopos.numHz,length(EPtopos.nonRegRelChans),4); %default is zero voltage if, for example, two factor datasets with different numbers of factors
    EPtopos.totalImagData=zeros(numChans,EPtopos.numPoints,EPtopos.numRows,EPtopos.numHz,length(EPtopos.nonRegRelChans),4); %default is zero voltage if, for example, two factor datasets with different numbers of factors
    EPtopos.colsForMax=[];
    EPtopos.complexData=0;
    EPtopos.plotForm=EPmain.view.dataTransform;
    if strcmp('TFT',EPmain.view.dataTransform) && (EPtopos.firstHz==EPtopos.lastHz)
        EPtopos.plotForm='VLT';
        EPtopos.type='time';
    end;
    if strcmp('TFT',EPmain.view.dataTransform) && (EPtopos.firstTime==EPtopos.lastTime)
        EPtopos.plotForm='FFT';
        EPtopos.type='freq';
    end;
    EPtopos.eventWave=cell(4,1);
    EPtopos.boundary=cell(4,1);
    theMax=0;
    for iColor=1:4
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            EPtopos.colorIndex(end+1)=iColor;
            if isempty(EPtopos.rowList)
                switch EPtopos.type
                    case 'subject'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'trial'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'factor'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'cell'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'time'
                        if EPtopos.flexMode
                            EPtopos.rowList=ceil([1:(5/EPtopos.sampleSize):EPtopos.numPoints]);
                        else
                            EPtopos.rowList=ceil([1:(50/EPtopos.sampleSize):EPtopos.numPoints]);
                        end;
                    case 'freq'
                        EPtopos.rowList=[1:EPtopos.numHz];
                    otherwise
                        error('Programming Error - aborting.');
                end;
            end;
            
            if theDataset ~= EPmain.view.dataset(iColor) %load in the dataset if different from the one already loaded in.
                EPdata=ep_loadEPdataset(EPdataset,EPmain.view.dataset(iColor));
                theDataset=EPmain.view.dataset(iColor);
            end;
            
            if EPmain.view.dataset(EPtopos.theFirstColor) == EPmain.view.dataset(iColor)
                EPtopos.colsForMax=[EPtopos.colsForMax; iColor];
            end;
            
            if any(strcmp(EPtopos.type,{'subject','time','freq','factor'})) || ~EPmain.view.allTrials(iColor)
                if isempty(EPdata.trialNames) %averaged data
                    theCell=EPmain.view.cell(iColor);
                else  %if single_trial data
                    theCell=intersect(find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdata.cellNames)),...
                        find(EPmain.view.trial(iColor)==EPdata.trialNames));
                end;
            else
                theCell=1;
            end;
            
            refChans=EPdata.reference.current;
            if isempty(refChans) && ~any(ismember({'AVG','CSD'},EPdata.reference.type))
                refChans=EPdata.reference.original;
            end;
            EPtopos.reference(iColor)=EPdata.reference;
            
            if any(EPmain.view.rel(EPmain.view.dataset <= length(EPdataset.dataset)))
                if strcmp(EPdata.dataType,'average')
                    EPtopos.goodRelChans{iColor}=find(squeeze(any(~isnan(EPdata.analysis.badChans(EPmain.view.subject(iColor),theCell,:)),2)));
                else
                    EPtopos.goodRelChans{iColor}=find(squeeze(any((EPdata.analysis.badChans(EPmain.view.subject(iColor),theCell,:)~=-1),2)));
                end;
                
                if length(refChans)==1
                    EPtopos.goodRelChans{iColor}=setdiff(EPtopos.goodRelChans{iColor},refChans); %coherence with a single reference channel is NaN.
                end;
            end;
            
            if strcmp(EPdata.reference.type,'CSD')
                if EPtopos.CSD==0
                    disp('Warning: Some data are CSD and some are not.  Labeling will be according to the last dataset.');
                end;
                EPtopos.CSD=1;
            else
                if EPtopos.CSD==1
                    disp('Warning: Some data are CSD and some are not.  Labeling will be according to the last dataset.');
                end;
                EPtopos.CSD=0;
            end;
            EPtopos.plotColors=[EPtopos.plotColors iColor];
            EPtopos.thePlotColors=[EPtopos.thePlotColors; EPtopos.RGBcolors(iColor,:)];
            EPtopos.plotLineIndex{end+1}='-';
            EPtopos.plotColorIndex=[EPtopos.plotColorIndex; EPtopos.RGBcolors(iColor,:)];
            
            switch EPtopos.type
                case 'subject'
                    if ~isempty(EPdata.facNames)
                        theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),theCell,[],EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                    else
                        theData=EPdata.data(EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),theCell,:,EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                    end;
                case {'trial','cell'}
                    if EPmain.view.allTrials(iColor)
                        if ~isempty(EPdata.facNames)
                            theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),[],EPmain.view.subject(iColor),EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                        else
                            theData=EPdata.data(EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),:,EPmain.view.subject(iColor),EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                        end;
                    else
                        if ~isempty(EPdata.facNames)
                            theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                        else
                            theData=EPdata.data(EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                        end;
                    end;
                case 'factor'
                    theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),theCell,EPmain.view.subject(iColor),[],startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                    if isfield(EPdata.pca,'PCAmode')
                        if isfield(EPdata.pca,'jack')
                            if isfield(EPdata.pca.jack,'FacPat')
                                for theJK=1:size(EPdata.pca.jack.FacPat,3)
                                    for theFac=1:size(EPdata.pca.jack.FacPat,2)
                                        EPtopos.jack(iColor).FacPat(:,theFac,theJK)=EPdata.pca.jack.FacPat(:,theFac,theJK).*EPdata.pca.jack.varSD(:,theJK);
                                        EPtopos.jack(iColor).PCAmode=EPdata.pca.PCAmode;
                                    end;
                                end;
                            elseif isfield(EPdata.pca.jack,'FacPatST')
                                for theJK=1:size(EPdata.pca.jack.FacPatST,3)
                                    for theFac=1:size(EPdata.pca.jack.FacPatST,2)
                                        EPtopos.jack(iColor).FacPat(:,theFac,theJK)=EPdata.pca.jack.FacPatST(:,theFac,theJK).*EPdata.pca.jack.varSDST(floor((theFac-1)/EPdata.pca.numFacs2)+1,:,theJK)';
                                        EPtopos.jack(iColor).PCAmode=EPdata.pca.PCAmode2;
                                    end;
                                end;
                            elseif isfield(EPdata.pca.jack,'FacPat3')
                                for theJK=1:size(EPdata.pca.jack.FacPat3,3)
                                    for theFac=1:size(EPdata.pca.jack.FacPat3,2)
                                        EPtopos.jack(iColor).FacPat(:,theFac,theJK)=EPdata.pca.jack.FacPat3(:,theFac,theJK).*EPdata.pca.jack.varSD3(floor((theFac-1)/EPdata.pca.numFacs3)+1,:,theJK)';
                                        EPtopos.jack(iColor).PCAmode=EPdata.pca.PCAmode3;
                                    end;
                                end;
                            end;
                        end;
                    end;
                case {'time','freq'}
                    if ~isempty(EPdata.facNames)
                        theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                    else
                        theData=EPdata.data(EPtopos.nonRegChans,startSamp(iColor):lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),startBins(iColor):lastBins(iColor),EPtopos.nonRegRelChans);
                    end;
                otherwise
                    error('Programming Error - aborting.');
            end;
            
            if ~isreal(theData) && (FFTunits ==1) %if the data has an imaginary component, as in spectral data
                theDataImag=imag(theData);
                theData=real(theData);
                EPtopos.complexData=1;
                EPtopos.plotLineIndex{end+1}=':';
                EPtopos.plotColorIndex=[EPtopos.plotColorIndex; EPtopos.RGBcolors(iColor,:)];
            else
                theDataImag=[];
            end;
            
            if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
                for iRow=1:EPtopos.numRows
                    switch EPtopos.type
                        case 'subject'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,4)
                                    theRowData=squeeze(theData(:,:,:,iRow,:,:,:));
                                else
                                    theRowData=[];
                                end;
                            else
                                theRowData=theData;
                            end;
                        case 'trial'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,3)
                                    theRowData=squeeze(theData(:,:,iRow,:,:,:,:));
                                else
                                    theRowData=[];
                                end;
                            else
                                theRowData=theData;
                            end;
                        case 'factor'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,5)
                                    theRowData=squeeze(theData(:,:,:,:,iRow,:,:));
                                else
                                    theRowData=[];
                                end;
                            else
                                theRowData=theData;
                            end;
                        case 'cell'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,3)
                                    theRowData=squeeze(theData(:,:,iRow,:,:,:,:));
                                else
                                    theRowData=[];
                                end;
                            else
                                theRowData=theData;
                            end;
                    end;
                    if ~isempty(theDataImag)
                        switch EPtopos.type
                            case 'subject'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,4)
                                        theRowDataImag=squeeze(theDataImag(:,:,:,iRow,:,:,:));
                                    else
                                        theRowDataImag=[];
                                    end;
                                else
                                    theRowData=theData;
                                end;
                            case 'trial'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,3)
                                        theRowDataImag=squeeze(theDataImag(:,:,iRow,:,:,:,:));
                                    else
                                        theRowDataImag=[];
                                    end;
                                else
                                    theRowData=theData;
                                end;
                            case 'factor'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,5)
                                        theRowDataImag=squeeze(theDataImag(:,:,:,:,iRow,:,:));
                                    else
                                        theRowDataImag=[];
                                    end;
                                else
                                    theRowData=theData;
                                end;
                            case 'cell'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,3)
                                        theRowDataImag=squeeze(theDataImag(:,:,iRow,:,:,:,:));
                                    else
                                        theRowDataImag=[];
                                    end;
                                else
                                    theRowDataImag=[];
                                end;
                        end;
                    end;
                    if ~isempty(theRowData)
                        EPtopos.totalData(1:size(theData,1),1:size(theData,2),iRow,1:size(theData,6),1:size(theData,7),iColor)=theRowData; %allows for combining datasets with differing dimensions
                    end;
                    if ~isempty(theDataImag) && ~isempty(theRowDataImag)
                        EPtopos.totalImagData(1:size(theDataImag,1),1:size(theDataImag,2),iRow,1:size(theDataImag,6),1:size(theDataImag,7),iColor)=theRowDataImag; %allows for combining datasets with differing dimensions
                    end;
                end;
            else
                EPtopos.totalData(1:size(theData,1),1:size(theData,2),1,1:size(theData,6),1:size(theData,7),iColor)=theData; %allows for combining datasets with differing dimensions
                if ~isempty(theDataImag)
                    EPtopos.totalImagData(1:size(theDataImag,1),1:size(theDataImag,2),1,1:size(theDataImag,6),1:size(theDataImag,7),iColor)=theDataImag; %allows for combining datasets with differing dimensions
                end;
            end;
            
            %event marks
            EPtopos.eventWave{iColor}=cell(EPtopos.numRows,1);
            EPtopos.boundary{iColor}=cell(EPtopos.numRows,1);
            for iRow=1:EPtopos.numRows
                theRow=EPtopos.rowList(iRow);
                tempEvents=[];
                if strcmp(EPdata.dataType,'continuous')
                    if strcmp('VLT',EPmain.view.dataTransform)
                        for theChan=1:size(EPtopos.totalData,1)
                            EPtopos.totalData(theChan,:,:,:,:)=EPtopos.totalData(theChan,:,:,:,:)-mean(EPtopos.totalData(theChan,:,:,:,:),2); %center the waveforms if continuous
                        end;
                    end;
                    if ~isempty(EPdata.events{EPmain.view.subject(iColor),1})
                        tempEvents=EPdata.events{EPmain.view.subject(iColor),1}(([EPdata.events{1}.sample]>=EPtopos.startSamp(iColor)) & ([EPdata.events{1}.sample]<=EPtopos.lastSamp(iColor)));
                    end;
                else
                    switch EPtopos.type
                        case 'subject'
                            tempEvents=EPdata.events{theRow,EPmain.view.cell(iColor)};
                            whichEvent=iRow;
                        case 'cell'
                            tempEvents=EPdata.events{EPmain.view.subject(iColor),theRow};
                            whichEvent=iRow;
                        otherwise
                            tempEvents=EPdata.events{EPmain.view.subject(iColor),EPmain.view.cell(iColor)};
                            whichEvent=1;
                    end;
                end;
                if ~isempty(EPtopos.eventLines{iColor}) && ~isempty(EPtopos.eventLines{iColor}{:}) && ~isempty(EPtopos.eventLines{iColor}{whichEvent})
                    EPtopos.eventWave{iColor}{iRow}=histc(EPtopos.eventLines{iColor}{whichEvent},[1:size(EPtopos.totalData,2)]);
                    theMax=max([theMax,max(EPtopos.eventWave{iColor}{whichEvent})]);
                end;
                if ~isempty(tempEvents)
                    boundaryEvents=find(strcmp('boundary',{tempEvents.value}));
                    if ~isempty(boundaryEvents)
                        EPtopos.boundary{iColor}{iRow}=tempEvents(boundaryEvents).sample;
                    end;
                end;
            end;
        end;
    end;
    for iColor=1:4
        for iRow=1:EPtopos.numRows
            if ~isempty(EPtopos.eventWave{iColor})
                if length(EPtopos.eventWave{iColor})>0
                    if ~isempty(EPtopos.eventWave{iColor}{iRow})
                        EPtopos.eventWave{iColor}{iRow}=EPtopos.eventWave{iColor}{iRow}/theMax; %rescale event information to between zero and one.
                    end;
                end;
            end;
        end;
    end;
    
    for color=1:4
        if EPmain.view.dataset(color) <= length(EPdataset.dataset) && ~EPmain.view.correl(color)
            if (EPmain.view.cell(color) > length(EPdataset.dataset(EPmain.view.dataset(color)).cellTypes)) || ~strcmp('STS',EPdataset.dataset(EPmain.view.dataset(color)).cellTypes(EPmain.view.cell(color)))
                %if the data are not correlations or STS output then rescale to the chosen units if FFT data
                if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
                    if (FFTunits > 1)
                        EPtopos.totalData(:,:,:,:,:,color)=abs(EPtopos.totalData(:,:,:,:,:,color)); %convert complex number to real number
                    end;
                    EPtopos.totalData(:,:,:,:,:,color)=EPtopos.totalData(:,:,:,:,:,color)/sqrt(mean(diff(EPdataset.dataset(EPmain.view.dataset(color)).freqNames))); %convert to spectral density
                    if FFTunits > 2
                        EPtopos.totalData(:,:,:,:,:,color)=EPtopos.totalData(:,:,:,:,:,color).^2; %convert amplitude to power
                    end;
                    if (FFTunits == 4)
                        if ~all(EPtopos.totalData(:,:,:,:,:,color) >=0)
                            disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
                        end;
                        EPtopos.totalData(:,:,:,:,:,color)=log10(abs(EPtopos.totalData(:,:,:,:,:,color)))*10; %convert to dB log scaling
                        tempVar=EPtopos.totalData(:,:,:,:,:,color);
                        tempVar(isinf(tempVar))=-flintmax;
                        EPtopos.totalData(:,:,:,:,:,color)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                    end;
                end;
            end;
        end;
    end;
end;

%check for sample test output
plotDataList=find(EPmain.view.dataset <= length(EPdataset.dataset));
EPtopos.STSdata=[];
for iColor=1:length(plotDataList)
    theColor=plotDataList(iColor);
    if (EPmain.view.cell(theColor) <= length(EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames))
        theCell=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(theColor)}(EPmain.view.cell(theColor)),EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames));
        theCell=theCell(1);
        if strcmp('STS',EPdataset.dataset(EPmain.view.dataset(theColor)).cellTypes(theCell))
            EPtopos.STSdata=[EPtopos.STSdata; theColor];
        end;
    end;
end;
if (length(EPtopos.STSdata)~=1)  || (length(plotDataList)~=3) || (length(EPtopos.plotColors)~=3)
    EPtopos.STmode=0;
else
    EPtopos.STmode=1; %present STS results as red zone between the two waves.
end;

if EPtopos.complexData
    EPtopos.perPage=floor(EPtopos.perPage/2)*2; %if complex, ensure even number of rows per page and reduce by half so can present both imaginary and real
end;

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    totalRows=length(EPtopos.rowList);
else
    totalRows=length(EPtopos.rowList)+1;
end;

if EPtopos.complexData
    EPtopos.numPages=ceil((totalRows*2)/EPtopos.perPage);
else
    EPtopos.numPages=ceil((totalRows)/EPtopos.perPage);
end;

if isempty(EPtopos.chans) %if first time through
%     %find max channels and latencies and Hz based on the selected data.
    maxData=EPtopos.totalData;    
    if any(EPmain.view.rel(EPtopos.colsForMax))
        maxData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.colsForMax)))=abs(maxData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.colsForMax))));
        maxData2=zeros([size(maxData,1),size(maxData,2),size(maxData,3),size(maxData,4),1,size(maxData,6)]);
        for iColor=1:length(EPtopos.colsForMax) %collapse over relations if any
            maxData2(:,:,:,:,:,iColor)=mean(maxData(:,:,:,:,EPtopos.goodRelChans{EPtopos.colsForMax(iColor)},iColor),5); %skip bad channels including the NaN for reference channel
        end;
        maxData=maxData2;
    end;
    
    if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
        theMax=-5;
    else
        theMax=0;
    end;
    switch EPtopos.type
        case {'subject','trial','cell','factor'}
            EPtopos.chans=ones(length(EPtopos.rowList),1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList),1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList),1); %default freq is one
            for theFactor=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1,[])));
                if length(EPtopos.colsForMax)==1
                    [C EPtopos.points(theFactor)]=max(max(reshape(shiftdim(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
                    [C EPtopos.chans(theFactor)]=max(max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])'));
                    [C EPtopos.freqs(theFactor)]=max(max(reshape(shiftdim(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
                else
                    [C EPtopos.points(theFactor)]=max(max(reshape(shiftdim(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
                    [C EPtopos.chans(theFactor)]=max(max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])'));
                    [C EPtopos.freqs(theFactor)]=max(max(reshape(shiftdim(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
                end;
            end;
        case 'time'
            EPtopos.chans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList)+1,1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
            [C EPtopos.points(1)]=max(max(reshape(shiftdim(maxData(:,:,1,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
            EPtopos.points(2:end)=EPtopos.rowList;
            EPtopos.rowList=[0 EPtopos.rowList];
            for thePoint=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),1,[])));
                [C EPtopos.chans(thePoint)]=max(max(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                [C EPtopos.freqs(thePoint)]=max(max(reshape(shiftdim(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
            end;
        case 'freq'
            EPtopos.chans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList)+1,1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
            theMax=max(theMax,max(reshape(maxData(:,:,1,:,:,EPtopos.colsForMax),1,[])));
            [C EPtopos.freqs(1)]=max(max(reshape(shiftdim(maxData(:,:,1,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
            EPtopos.freqs(2:end)=EPtopos.rowList;
            EPtopos.rowList=[0 EPtopos.rowList];
            for theFreq=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),1,[])));
                [C EPtopos.points(theFreq)]=max(max(reshape(shiftdim(maxData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),1),EPtopos.numPoints,[])',[],1));
                [C EPtopos.chans(theFreq)]=max(max(reshape(maxData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
            end;
        case 'chan'
            EPtopos.chans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList)+1,1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
            [C EPtopos.chans(1)]=max(max(reshape(maxData(:,:,1,:,EPtopos.colsForMax),length(EPtopos.nonRegRelChans),[])'));
            EPtopos.chans(2:end)=EPtopos.rowList;
            EPtopos.rowList=[0 EPtopos.rowList];
            for theChan=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),1,[])));
                [C EPtopos.points(theChan)]=max(max(reshape(maxData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),EPtopos.numPoints,[])',[],1));
                [C EPtopos.freqs(theChan)]=max(max(reshape(shiftdim(maxData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
            end;
        otherwise
            error('Programming Error - aborting.');
            return
    end;
    if EPtopos.complexData
        maxData=EPtopos.totalImagData;
        if (EPtopos.FFTunits ==4) && any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'})) && ~all(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
            if any(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
                %if mixing coherence (correlations) and regular FFT under dB scaling, then undo the rescaling to microvolts,
                %take absolute value of the correlations, then redo the rescaling
                maxData(:,:,:,:,:,find(EPmain.view.correl))=(maxData(:,:,:,:,:,find(EPmain.view.correl))-EPtopos.plotMVmin)/(EPtopos.plotMVmax-EPtopos.plotMVmin);
                maxData(:,:,:,:,:,find(EPmain.view.correl))=abs(maxData(:,:,:,:,:,find(EPmain.view.correl)));
                maxData(:,:,:,:,:,find(EPmain.view.correl))=(maxData(:,:,:,:,:,find(EPmain.view.correl))*(EPtopos.plotMVmax-EPtopos.plotMVmin))+EPtopos.plotMVmin;
            end;
        elseif (EPtopos.FFTunits < 4) || ~any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
            maxData=abs(maxData);
        elseif any(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
            maxData(:,:,:,:,:,find(EPmain.view.correl))=abs(maxData(:,:,:,:,:,find(EPmain.view.correl)));
        end;
        maxData=mean(maxData,5); %collapse over relations if any
        
        switch EPtopos.type
            case {'subject','trial','cell','factor'}
                EPtopos.iChans=ones(length(EPtopos.rowList),1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList),1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList),1); %default freq is one
                for theFactor=1:length(EPtopos.rowList)
                    [C EPtopos.iPoints(theFactor)]=max(max(reshape(shiftdim(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
                    [C EPtopos.iChans(theFactor)]=max(max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])'));
                    [C EPtopos.iFreqs(theFactor)]=max(max(reshape(shiftdim(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
                end;
            case 'time'
                EPtopos.rowList=EPtopos.rowList(2:end); %undo the addition of the first row so it doesn't end up getting added twice.
                EPtopos.iChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList)+1,1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
                [C EPtopos.points(1)]=max(max(reshape(shiftdim(maxData(:,:,1,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
                
                EPtopos.iPoints(2:end)=EPtopos.rowList;
                EPtopos.rowList=[0 EPtopos.rowList];
                for thePoint=1:length(EPtopos.rowList)
                    [C EPtopos.iChans(thePoint)]=max(max(reshape(maxData(:,EPtopos.iPoints(thePoint),:,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                    [C EPtopos.iFreqs(thePoint)]=max(max(reshape(shiftdim(maxData(:,EPtopos.iPoints(thePoint),:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
                end;
            case 'freq'
                EPtopos.rowList=EPtopos.rowList(2:end); %undo the addition of the first row so it doesn't end up getting added twice.
                EPtopos.iChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList)+1,1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
                [C EPtopos.iFreqs(1)]=max(max(reshape(shiftdim(maxData(:,:,1,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
                
                EPtopos.iFreqs(2:end)=EPtopos.rowList;
                EPtopos.rowList=[0 EPtopos.rowList];
                for theFreq=1:length(EPtopos.rowList)
                    [C EPtopos.iPoints(theFreq)]=max(max(reshape(shiftdim(maxData(:,:,:,EPtopos.iFreqs(theFreq),:,EPtopos.colsForMax),1),EPtopos.numPoints,[])',[],1));
                    [C EPtopos.iChans(theFreq)]=max(max(reshape(maxData(:,:,:,EPtopos.iFreqs(theFreq),:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                end;
            case 'chan'
                EPtopos.rowList=EPtopos.rowList(2:end); %undo the addition of the first row so it doesn't end up getting added twice.
                EPtopos.iChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList)+1,1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
                [C EPtopos.iChans(1)]=max(max(reshape(maxData(:,:,1,:,:,EPtopos.colsForMax),length(EPtopos.nonRegRelChans),[])'));
                
                EPtopos.iChans(2:end)=EPtopos.rowList;
                EPtopos.rowList=[0 EPtopos.rowList];
                for theChan=1:length(EPtopos.rowList)
                    [C EPtopos.iPoints(theChan)]=max(max(reshape(shiftdim(maxData(EPtopos.iChans(theChan),:,:,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])',[],1));
                    [C EPtopos.iFreqs(theChan)]=max(max(reshape(shiftdim(maxData(EPtopos.iChans(theChan),:,:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
                end;
            otherwise
                error('Programming Error - aborting.');
                return
        end;
    end;
    if ~isempty(minVolt)
        EPtopos.plotMVmin=minVolt;
    else
        if theMax < 0
            EPtopos.plotMVmin=theMax*2;
        else
            EPtopos.plotMVmin=-theMax;
        end;
    end;
    if ~isempty(maxVolt)
        EPtopos.plotMVmax=maxVolt;
    else
        EPtopos.plotMVmax=theMax;
    end;
    if ~all(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset))) && any(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset))) 
        %if not all the waveforms are correlations but some are, then rescale the correlations to match the plots
        EPtopos.totalData(:,:,:,:,:,find(EPmain.view.correl))=(EPtopos.totalData(:,:,:,:,:,find(EPmain.view.correl))*(EPtopos.plotMVmax-EPtopos.plotMVmin))+EPtopos.plotMVmin;
        if EPtopos.complexData
            EPtopos.totalImagData(:,:,:,:,:,find(EPmain.view.correl))=(EPtopos.totalImagData(:,:,:,:,:,find(EPmain.view.correl))*(EPtopos.plotMVmax-EPtopos.plotMVmin))+EPtopos.plotMVmin;
        end;
    end;
end;

if EPtopos.page < EPtopos.numPages
    numRows=EPtopos.perPage;
else
    if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
        if EPtopos.complexData
            numRows=mod((EPtopos.numRows*2)-1,EPtopos.perPage)+1;
        else
            numRows=mod(EPtopos.numRows-1,EPtopos.perPage)+1;
        end;
    else
        if EPtopos.complexData
            numRows=mod((length(EPtopos.rowList)*2)-1,EPtopos.perPage)+1;
        else
            numRows=mod(length(EPtopos.rowList)-1,EPtopos.perPage)+1;
        end;
    end;
end;

clf(EPtopos.handles.topos.topoWindow)

EPtopos.handles.synch = uicontrol('Style', 'checkbox', 'String', 'Synch','Value',EPtopos.synch,'FontSize',EPmain.fontsize,'Position', [20 windowHeight-80 60 30],...
    'TooltipString','When checked, changes in time or channel to one row changes all of them.',...
    'Callback', ['global EPtopos;','EPtopos.synch=get(EPtopos.handles.synch,''Value'');']);

EPtopos.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-110 50 30], 'Callback', @done);

EPtopos.handles.plotMVmin = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
    'String',num2str(EPtopos.plotMVmin),'HorizontalAlignment','left',...
    'Position',[120 windowHeight-110 50 20],'Callback',@changeMV);

EPtopos.handles.plotMVmax = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
    'String',num2str(EPtopos.plotMVmax),'HorizontalAlignment','left',...
    'Position',[180 windowHeight-110 50 20],'Callback',@changeMV);

theLabel='';

if all(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
    theLabel='Correlation';
elseif strcmp(EPmain.view.dataTransform,'VLT')
    if strcmp(EPtopos.CSD,'CSD')
        theLabel='V/m^2';
    else
        theLabel='µv';
    end;
elseif any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
    if (EPtopos.FFTunits > 2)
        if (EPtopos.FFTunits == 4)
        	theLabel='dBV';
        else
            if strcmp(EPtopos.CSD,'CSD')
                theLabel='(V^2)/(Hz*m^4)';
            else
                theLabel='(µv^2)/Hz';
            end;
        end;
    else
        if strcmp(EPtopos.CSD,'CSD')
            theLabel='V/(sqrt(Hz)*m^2)';
        else
            theLabel='µv/sqrt(Hz)';
        end;
    end;
end;

uicontrol('Style','text','HorizontalAlignment','left','String',theLabel,'FontSize',EPmain.fontsize,...
    'Position',[130 windowHeight-65 50 20]);

EPtopos.handles.topos.colorbar = axes('Units','pixel','position',[120 windowHeight-95 109 20],'Visible','off');
colorbar('location','southoutside','Visible','on', 'XTick', [])

for color=1:4
    if EPmain.view.dataset(color) <= length(EPdataset.dataset)
        uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
            'String',EPdataset.dataset(EPmain.view.dataset(color)).dataName,'HorizontalAlignment','left',...
            'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-70 80 20]);
    end;
    if EPmain.view.dataset(color) <= length(EPdataset.dataset)
        if EPmain.view.allTrials(color)==5
            uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                'String','-all cells-','HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-90 80 20]);
        else
            uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                'String',EPdataset.dataset(EPmain.view.dataset(color)).cellNames{EPmain.view.cell(color)},'HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-90 80 20]);
        end;
    end;
    if EPmain.view.dataset(color) <= length(EPdataset.dataset)
        if strcmp(EPdataset.dataset(EPmain.view.dataset(color)).dataType,'single-trial')
            if EPmain.view.allTrials(color)==1
                uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                    'String','-all trials-','HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-110 80 20]);
            else
                uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                    'String',EPdataset.dataset(EPmain.view.dataset(color)).trialNames{EPmain.view.trial(color)},'HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-110 80 20]);
            end;
        else
            if EPmain.view.allTrials(color)==1
                uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                    'String','-all subs-','HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-110 80 20]);
            else
                uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                    'String',EPdataset.dataset(EPmain.view.dataset(color)).subNames{EPmain.view.subject(color)},'HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-110 80 20]);
            end;
        end;
    end;
    if (EPmain.view.dataset(color) <= length(EPdataset.dataset)) && ~isempty(EPdataset.dataset(EPmain.view.dataset(color)).facNames)
        if EPmain.view.allTrials(color)==3
            uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                'String','-all facs-','HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-130 80 20]);
        else
            uicontrol('Style','text','ForegroundColor',textColors{color},'FontSize',EPmain.fontsize,...
                'String',EPdataset.dataset(EPmain.view.dataset(color)).facNames{EPmain.view.factor(color)},'HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(color-1)*(EPtopos.topoSize+20) windowHeight-130 80 20]);
        end;
    end;
end;

if EPtopos.numPages > 1
    %page buttons
    if EPtopos.numPages <= 19
        for i=1:EPtopos.numPages
            EPtopos.handles.pageButton(i) = uicontrol('Style', 'pushbutton', 'String', num2str(i*EPtopos.perPage),'FontSize',EPmain.fontsize,...
                'Position', [(20 + (i-1)*40) windowHeight-40 30 20], 'Callback', ['global EPtopos;','EPtopos.page=' num2str(i) ';','ep_showTopos(EPtopos.firstTime,EPtopos.lastTime)']);
        end;
        set(EPtopos.handles.pageButton(EPtopos.page),'ForegroundColor','blue');
    else
        EPtopos.handles.pageMenu = uicontrol('Style', 'popupmenu', 'String', num2str([ceil((1:EPtopos.numPages)*EPtopos.perPage)]'),'Value',EPtopos.page,'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-40 100 20],...
            'Callback', ['global EPtopos;','tempVar=get(EPtopos.handles.pageMenu,''Value'');','if tempVar ~=0,EPtopos.page=tempVar;end;','ep_showTopos(EPtopos.firstTime,EPtopos.lastTime)']);
    end;
end;

%start drawing the figures
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesXTickLabel',[])
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesYTickLabel',[])
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesXTick',[])
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesYTick',[])
set(EPtopos.handles.topos.topoWindow,'DefaultAxesColorOrder',EPtopos.thePlotColors)

numChans=length(EPtopos.eloc);
%start adapted code from EEGlab's topoplot

locChans=true(numChans,1);
Rd=[];
Th=[];
for chan=1:numChans
    theRd=EPtopos.eloc(chan).radius;
    if isempty(theRd)
        locChans(chan)=false;
    end;
    theTheta=EPtopos.eloc(chan).theta;
    if isempty(theTheta)
        locChans(chan)=false;
    end;
    if locChans(chan)
        Rd=[Rd theRd];
        Th=[Th theTheta];
    end;
end;
EPtopos.locChans=find(locChans);

Th = pi/180*Th;
[x,y] = pol2cart(Th,Rd);
plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
squeezefac = .5/plotrad;
x    = x*squeezefac;
y    = y*squeezefac;
%end code from EEGlab's topoplot

if EPtopos.page ==1 && ~any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    if EPtopos.complexData
        uicontrol('Style','frame',...
            'Position',[0 windowHeight-345 800 1]);
    else
        uicontrol('Style','frame',...
            'Position',[0 windowHeight-225 800 1]);
    end;
end;

for iRow=1:numRows
    if EPtopos.complexData
        rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(iRow/2);
    else
        rowCounter=(EPtopos.page-1)*EPtopos.perPage+iRow;
    end;
    theRow=1;
    if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
        theRow=rowCounter;
    end;
    
    switch EPtopos.type
        case 'subject'
            uicontrol('Style','text',...
                'String',EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).subNames{rowCounter},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
        case 'trial'
            uicontrol('Style','text',...
                'String',num2str(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).trialNames{rowCounter}),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
        case 'factor'
            uicontrol('Style','text',...
                'String',EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).facNames{rowCounter},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
        case 'cell'
            uicontrol('Style','text',...
                'String',EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).cellNames{rowCounter},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
    end;

    peakPoint=NaN;
    peakHz=NaN;
    if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
        EPtopos.sampleSize=EPtopos.sampleSize;
        EPtopos.spacing=EPtopos.sampleSize;
        if EPtopos.complexData && rem(iRow,2)
            peakPoint=EPtopos.firstTime+(EPtopos.iPoints(rowCounter)-1)*EPtopos.sampleSize;
        else
            peakPoint=EPtopos.firstTime+(EPtopos.points(rowCounter)-1)*EPtopos.sampleSize;
        end
    end;
    if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
        if EPtopos.complexData && rem(iRow,2)
            peakHz=EPtopos.iFreqs(rowCounter);
        else
            peakHz=EPtopos.freqs(rowCounter);
        end;
    else
        EPtopos.freqNameList='none';
        peakHz=1;
    end;
    if EPtopos.complexData && rem(iRow,2)
        peakChan=EPtopos.iChans(rowCounter);
    else
        peakChan=EPtopos.chans(rowCounter);
    end;
    EPtopos.handles.Hz(iRow)=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPtopos.freqNameList,'HorizontalAlignment','left',...
        'Value',peakHz,...
        'Position',[15 windowHeight-(160+(iRow-1)*120) 80 20],...
        'Callback',{@changeHz,iRow});
    EPtopos.handles.HzLabel(iRow)=uicontrol('Style','text','HorizontalAlignment','left','String', 'Hz','FontSize',EPmain.fontsize,...
        'Position',[90 windowHeight-(165+(iRow-1)*120) 20 20]);
    if  strcmp('VLT',EPmain.view.dataTransform)
        set(EPtopos.handles.Hz(iRow),'enable','off');
        set(EPtopos.handles.HzLabel(iRow),'enable','off');
    end
    EPtopos.handles.latency(iRow)=uicontrol('Style','edit','FontSize',EPmain.fontsize,...
        'String',num2str(round(peakPoint)),'HorizontalAlignment','left',...
        'Position',[20 windowHeight-(180+(iRow-1)*120) 50 20],...
        'Callback',{@changeLatency,iRow});
    if EPtopos.flexMode
        theUnit='%';
    else
        theUnit='ms';
    end;
     EPtopos.handles.latencyLabel(iRow)=uicontrol('Style','text','HorizontalAlignment','left','String',theUnit,'FontSize',EPmain.fontsize,...
        'Position',[70 windowHeight-(185+(iRow-1)*120) 20 20]);
    if  strcmp('FFT',EPmain.view.dataTransform)
        set(EPtopos.handles.latency(iRow),'enable','off');
        set(EPtopos.handles.latencyLabel(iRow),'enable','off');
    end
    EPtopos.handles.channel(iRow)=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPtopos.chanNames(EPtopos.nonRegChans)',...
        'Value',peakChan,...
        'Position',[15 windowHeight-(200+(iRow-1)*120) 80 20],...
        'Callback',{@changeChannel,iRow});
    
    %draw waveforms
    
    if strcmp('FFT',EPtopos.plotForm)
        EPtopos.sampleSize=0;
        EPtopos.spacing=(EPtopos.lastHz-EPtopos.firstHz)/(EPtopos.numHz-1);
    else
        EPtopos.sampleSize=EPtopos.sampleSize;
        EPtopos.spacing=EPtopos.sampleSize;
    end;
    
    switch EPtopos.plotForm
        case 'VLT'
            if EPtopos.numPoints > 1
                EPtopos.handles.waves.hWave(iRow) = axes('Units','pixel','position',[120 windowHeight-(220+(iRow-1)*120) EPtopos.plotWidth EPtopos.plotHeight],'XTickMode','manual','YTickMode','manual');
                
                if EPtopos.complexData && rem(iRow,2)
                    theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                else
                    theData=EPtopos.totalData(EPtopos.chans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                end
                if any(EPmain.view.rel(EPmain.view.dataset <= length(EPdataset.dataset)))
                    theData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors)))=abs(theData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors))));
                    theData=mean(theData,5); %collapse over relations if any
                end;
                theData=squeeze(theData);
                if length(EPtopos.plotColors)==1
                    theData=theData(:);
                end;
                
                hold on
                for iWave=1:size(theData,2)
                    if EPtopos.STmode && (EPtopos.STSdata==iWave)
                        breakList=sort([find(diff([0 (theData(:,iWave)'>0) 0])<0)-1 find(diff([0 (theData(:,iWave)'>0) 0])>0)]);
                        if ~isempty(breakList)
                            theSTdata=theData(:,iWave);
                            theData1=theData(:,min(setdiff(EPtopos.plotColors,EPtopos.STSdata)));
                            theData2=theData(:,max(setdiff(EPtopos.plotColors,EPtopos.STSdata)));
                            for iSigArea=1:length(breakList)/2
                                theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)];
                                patch(([theTimePoints flip(theTimePoints)]*EPtopos.spacing)+(EPtopos.firstTime-EPtopos.spacing),[theData1(theTimePoints)' theData2(flip(theTimePoints))'],EPtopos.thePlotColors(iWave,:),'FaceColor',EPtopos.thePlotColors(iWave,:),'EdgeColor','none','FaceAlpha',.25);
                            end;
                        end;
                    else
                        EPtopos.handles.waves.hLines{iRow} = plot([EPtopos.firstTime:EPtopos.spacing:EPtopos.lastTime],theData(:,iWave),'color',EPtopos.thePlotColors(iWave,:));
                    end;
                end;
                hold off
                
                axis([EPtopos.firstTime EPtopos.lastTime EPtopos.plotMVmin EPtopos.plotMVmax]); %left side of first sample to left side of last sample
                if EPtopos.direction ==2
                    set(EPtopos.handles.waves.hWave(iRow),'YDir','reverse')
                end;
                line([EPtopos.firstTime EPtopos.lastTime],[0 0],'Color','black','LineWidth',1) % zero line
                line([0 0],[0 EPtopos.plotMVmax],'Color','black','LineWidth',1) %stimulus onset
                %text(EPtopos.firstTime+(EPtopos.lastTime-EPtopos.firstTime)/20, EPtopos.plotMVmin+(EPtopos.plotMVmax-EPtopos.plotMVmin)/10, EPtopos.chanNames(EPtopos.chans(rowCounter)));
                line(repmat(peakPoint,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
                if ~isempty(EPtopos.marker1)
                    line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
                end;
                if ~isempty(EPtopos.marker2)
                    line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
                end;
                set(EPtopos.handles.waves.hWave(iRow),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                for i=1:length(EPtopos.handles.waves.hLines{iRow})
                    set(EPtopos.handles.waves.hLines{iRow}(i),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                end;

                %plot event lines
                for iColor=1:length(EPtopos.plotColors)
                    theColor=EPtopos.plotColors(iColor);
                    switch EPtopos.type
                        case 'subject'
                            theEventWave=rowCounter;
                        case 'cell'
                            theEventWave=rowCounter;
                        otherwise
                            theEventWave=1;
                    end;
                    if ~isempty(EPtopos.eventWave{theColor}{theEventWave}) && (EPtopos.plotMVmin < 0)
                        plotPoints=find(EPtopos.eventWave{theColor}{theEventWave}>min(EPtopos.eventWave{theColor}{theEventWave}));
                        plotTimes=[EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))];
                        if (EPtopos.plotMVmin < 0) && (EPtopos.plotMVmax >= 0)
                            if length(plotPoints)==1
                                line([plotTimes(plotPoints) plotTimes(plotPoints)],[EPtopos.plotMVmin EPtopos.eventWave{theColor}{theEventWave}(plotPoints)*(EPtopos.plotMVmin/2)],'Color',EPtopos.thePlotColors(theColor,:),'LineWidth',2) %event line
                            else
                                hold on
                                EPtopos.handles.waves.eventLines{iRow,theColor} = plot([EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))],(EPtopos.eventWave{theColor}{theEventWave}*(abs(EPtopos.plotMVmin/2)))+EPtopos.plotMVmin,'LineWidth',5,'Color',EPtopos.thePlotColors(theColor,:));
                                hold off
                                for iLine=1:length(EPtopos.handles.waves.eventLines{iRow,theColor})
                                    set(EPtopos.handles.waves.eventLines{iRow,theColor}(iLine),'YDataSource',['EPtopos.eventWave{' num2str(theColor) '}(' num2str(iLine) ',:)']);
                                end;
                            end;
                        end;
                    end;
                    if ~isempty(EPtopos.boundary{theColor})
                        hold on
                        for iBoundary=1:length(EPtopos.boundary{theColor}{theEventWave})
                            theSample=EPtopos.boundary{theColor}{theEventWave}(iBoundary);
                            EPtopos.handles.waves.boundary{iRow,theColor} = line([theSample theSample],[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',EPtopos.thePlotColors(theColor,:),'LineWidth',1);
                        end;
                        hold off
                    end;
                end;
            end;
        case 'FFT'
            if EPtopos.numHz > 1
                EPtopos.handles.waves.hWave(iRow) = axes('Units','pixel','position',[120 windowHeight-(220+(iRow-1)*120) EPtopos.plotWidth EPtopos.plotHeight],'XTickMode','manual','YTickMode','manual');
                
                if EPtopos.complexData && rem(iRow,2)
                    theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                    theLineStyle=':';
                else
                    theData=EPtopos.totalData(EPtopos.chans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                    theLineStyle='-';
                end
                
                if any(EPmain.view.rel(EPmain.view.dataset <= length(EPdataset.dataset)))
                    theData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors)))=abs(theData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors))));
                    theData2=zeros([size(theData,1),size(theData,2),size(theData,3),size(theData,4),1,size(theData,6)]);
                    for iColor=1:length(EPtopos.plotColors)
                        theData2(:,:,:,:,:,iColor)=mean(theData(:,:,:,:,EPtopos.goodRelChans{EPtopos.plotColors(iColor)},iColor),5); %collapse over relations if any
                    end;
                    theData=theData2;
                end;
                theData=squeeze(theData);
                if length(EPtopos.plotColors)==1
                    theData=theData(:);
                end;
                
                EPtopos.handles.waves.hLines{iRow} = plot([EPtopos.firstHz:EPtopos.spacing:EPtopos.lastHz],theData,'LineStyle',theLineStyle);
                axis([EPtopos.firstHz EPtopos.lastHz EPtopos.plotMVmin EPtopos.plotMVmax]);
                
                line([EPtopos.firstHz EPtopos.lastHz],[0 0],'Color','black','LineWidth',1) % zero line
                line(repmat(peakHz,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
                if ~isempty(EPtopos.marker1)
                    line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
                end;
                if ~isempty(EPtopos.marker2)
                    line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
                end;
                set(EPtopos.handles.waves.hWave(iRow),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                for i=1:length(EPtopos.handles.waves.hLines{iRow})
                    set(EPtopos.handles.waves.hLines{iRow}(i),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                end;
            end;
        case 'TFT'
            if EPtopos.numPoints > 1
                imageSpace=length(EPtopos.plotColors);
                imageCount=0;
                for i=1:4
                    if EPmain.view.dataset(i) <= length(EPdataset.dataset)
                        imageCount=imageCount+1;
                        EPtopos.handles.waves.hLines{iRow}(i) = axes('Units','pixel','position',[120 windowHeight-(220+(iRow-1)*120)+(EPtopos.plotHeight/imageSpace)*(imageSpace-imageCount) EPtopos.plotWidth (EPtopos.plotHeight/imageSpace)]);
                        
                        if EPtopos.complexData && rem(iRow,2)
                            theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),:,theRow,:,:,i);
                        else
                            theData=EPtopos.totalData(EPtopos.chans(rowCounter),:,theRow,:,:,i);
                        end
                        if EPmain.view.rel(i)
                            theData=abs(theData);
                            theData=mean(theData(:,:,:,:,EPtopos.goodRelChans{i},1),5); %collapse over relations if any
                        end;
                        theData=squeeze(theData)';
                        
                        EPtopos.handles.waves.hLines{iRow}(4+i) = imagesc(EPtopos.firstTime:EPtopos.lastTime+EPtopos.sampleSize,EPtopos.firstHz:EPtopos.lastHz,theData,[EPtopos.plotMVmin EPtopos.plotMVmax]);
                        axis([EPtopos.firstTime EPtopos.lastTime+EPtopos.sampleSize EPtopos.firstHz EPtopos.lastHz]);
                        line([0 0],[EPtopos.firstHz EPtopos.lastHz],'Color','black','LineWidth',1) %stimulus onset
                        line([EPtopos.firstTime EPtopos.lastTime],[peakHz peakHz],'Color','white','LineWidth',1) % Hz line
                        line([peakPoint peakPoint],[EPtopos.firstHz EPtopos.lastHz],'Color',[.5 .5 .5],'LineWidth',1); %ms line
                        if ~isempty(EPtopos.marker1)
                            line(repmat(EPtopos.marker1,2),[EPtopos.firstHz EPtopos.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
                        end;
                        if ~isempty(EPtopos.marker2)
                            line(repmat(EPtopos.marker2,2),[EPtopos.firstHz EPtopos.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
                        end;
                        set(EPtopos.handles.waves.hLines{iRow}(i),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                        set(EPtopos.handles.waves.hLines{iRow}(4+i),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                    end;
                end;
            end;
    end;
    
    %draw topos
    for iCol=1:4
        if EPmain.view.dataset(iCol) <= length(EPdataset.dataset)
            EPtopos.handles.topos.topo(iRow,iCol) = axes('Units','pixel','position',[100+EPtopos.plotWidth+50+(iCol-1)*(EPtopos.topoSize+20) windowHeight-(220+(iRow-1)*120) EPtopos.topoSize EPtopos.topoSize]);
            if EPmain.view.rel(iCol) %if relational data
                if EPtopos.complexData && rem(iRow,2)
                    theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),EPtopos.iChans(rowCounter),iCol);
                else
                    theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),EPtopos.chans(rowCounter),iCol);
                end
            else
                if EPtopos.complexData && rem(iRow,2)
                    theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),1,iCol);
                else
                    theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),1,iCol);
                end
            end;
            theData=squeeze(theData);
            theData(isnan(theData))=0;
            el_topoplot(theData, EPtopos.eloc(EPtopos.nonRegChans),'maplimits',[EPtopos.plotMVmin EPtopos.plotMVmax]);
            if EPtopos.complexData && rem(iRow,2)
                line([y(EPtopos.iChans(rowCounter)) y(EPtopos.iChans(rowCounter))],[x(EPtopos.iChans(rowCounter)) x(EPtopos.iChans(rowCounter))],'Marker','o','MarkerFaceColor','white'); %add marker for the electrode of the waveform plot
            else
                line([y(EPtopos.chans(rowCounter)) y(EPtopos.chans(rowCounter))],[x(EPtopos.chans(rowCounter)) x(EPtopos.chans(rowCounter))],'Marker','o','MarkerFaceColor','white'); %add marker for the electrode of the waveform plot
            end;
            for i=1:length(EPtopos.locChans)
                EPtopos.handles.topos.line(iRow,iCol,i)=line([y(i) y(i)],[x(i) x(i)],'color','black','ButtonDownFcn',['global EPtopos;','sel_typ = get(gcbf,''SelectionType'');','if strcmp(sel_typ,''normal'')','EPtopos.chans(' num2str(rowCounter) ')=' num2str(EPtopos.locChans(i)) ';','ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);','end;']);
            end;
            
            %         EPtopos.handles.topos.d2Button(iRow,iCol) = uicontrol('Style', 'pushbutton', 'String', '2D',...
            %             'Position', [100+EPtopos.plotWidth+50+(iCol-1)*(EPtopos.topoSize+50) windowHeight-(200+(iRow-1)*100)-20 30 20],...
            %             'Callback', ['global EPtopos;','figure;','el_topoplot(EPtopos.totalData(:,EPtopos.points(' num2str(rowCounter) '),' num2str(rowCounter) ',' num2str(iCol) '), EPtopos.eloc,''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]);']);
            
            hcmenu = uicontextmenu;
            % Define callbacks for context menu items that change linestyle
            
            if EPtopos.numRows > 1
                theRow=rowCounter;
            else
                theRow=1;
            end;
            
            if EPmain.view.rel(iCol) %if relational data
                if EPtopos.complexData && rem(iRow,2)
                    hcb1 = ['global EPtopos;','figure;','el_topoplot(EPtopos.totalImagData(EPtopos.iChans(' num2str(rowCounter) '),EPtopos.iPoints(' num2str(rowCounter) '),' num2str(theRow) ',EPtopos.iFreqs(' num2str(rowCounter) '),:,' num2str(iCol) '), EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                else
                    hcb1 = ['global EPtopos;','figure;','el_topoplot(EPtopos.totalData(EPtopos.chans(' num2str(rowCounter) '),EPtopos.points(' num2str(rowCounter) '),' num2str(theRow) ',EPtopos.freqs(' num2str(rowCounter) '),:,' num2str(iCol) '), EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                end
            else
                if EPtopos.complexData && rem(iRow,2)
                    hcb1 = ['global EPtopos;','figure;','el_topoplot(EPtopos.totalImagData(:,EPtopos.iPoints(' num2str(rowCounter) '),' num2str(theRow) ',EPtopos.iFreqs(' num2str(rowCounter) '),1,' num2str(iCol) '), EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                else
                    hcb1 = ['global EPtopos;','figure;','el_topoplot(EPtopos.totalData(:,EPtopos.points(' num2str(rowCounter) '),' num2str(theRow) ',EPtopos.freqs(' num2str(rowCounter) '),1,' num2str(iCol) '), EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                end
            end;
            hcb3 = ['set(gco, ''LineStyle'', ''-'')'];
            % Define the context menu items and install their callbacks
            item1 = uimenu(hcmenu, 'Label', '2D', 'Callback', hcb1);
            item2 = uimenu(hcmenu, 'Label', '3D', 'Callback', @D3head);
            item3 = uimenu(hcmenu, 'Label', 'Rescale',  'Callback', @rescaleFigures);
            item4 = uimenu(hcmenu, 'Label', '2-Channels',  'Callback', @twoChan);
            if strcmp('VLT',EPmain.view.dataTransform) && ~strcmp(EPtopos.CSD,'CSD')
                item5 = uimenu(hcmenu, 'Label', 'Rereference',  'Callback', @referenceChan);
            end;
            if strcmp('VLT',EPmain.view.dataTransform) || ((strcmp('FFT',EPmain.view.dataTransform) && EPmain.view.rel(iCol))) && ~strcmp(EPtopos.CSD,'CSD')
                item6 = uimenu(hcmenu, 'Label', 'Dipole',  'Callback', @dipoles);
                if isfield(EPtopos,'jack')
                    if length(EPtopos.jack) >=iCol
                        if strcmp(EPtopos.jack(iCol).PCAmode,'spat')
                            item6 = uimenu(hcmenu, 'Label', 'Jack-Knife',  'Callback', @jackknife);
                        end;
                    end;
                end;
            end;
            % Attach the context
            set(EPtopos.handles.topos.topo(iRow,iCol),'UIContextMenu',hcmenu)
            theChildren=get(EPtopos.handles.topos.topo(iRow,iCol),'Children');
            for i=1:length(theChildren)
                set(theChildren(i),'UIContextMenu',hcmenu)
            end;
            for i=1:length(EPtopos.locChans)
                set(EPtopos.handles.topos.line(iRow,iCol,i),'UIContextMenu',hcmenu)
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function done(src,eventdata)
global EPtopos EPmain EPdataset

%save results of jack-knife PCA dipole analyses
if isfield(EPtopos,'jack')
    for col=1:length(EPtopos.jack)
        if EPmain.view.dataset(col) <= length(EPdataset.dataset) %if not set to "none"
            if ~isempty(EPtopos.jack(col))
                if isfield(EPtopos.jack(col),'sources')
                    if ~isempty(EPtopos.jack(col).sources)
                        EPdata=ep_loadEPdataset(EPdataset,EPmain.view.dataset(col));
                        EPdata.pca.jack.sources=EPtopos.jack(col).sources;
                        
                        try
                            EPver=ver('EP_Toolkit');
                        catch
                            EPver='unavailable'; %workaround for bug in earlier version of Matlab
                        end;
                        EPdata.EPver=EPver;
                        EPdata.ver=ver;
                        EPdata.date=date;
                        EPdata.history(end+1)={'PCA jackknife source analysis'};
                        
                        [err]=ep_checkEPfile(EPdata);
                        
                        if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
                            warndlg('The work directory cannot be found.')
                            return
                        end;
                        
                        delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(EPmain.view.dataset(col)).dataName '.mat']);
                        EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],EPmain.view.dataset(col)));
                        EPdataset=ep_saveEPdataset(EPdataset,EPdata,EPmain.view.dataset(col),'no');
                    end;
                end;
            end;
        end;
    end;
end;

set(EPmain.handles.view.waves,'enable','on');
set(EPmain.handles.view.topos,'enable','on');
set(EPmain.handles.view.done,'enable','on');
set(EPmain.handles.view.marker1,'enable','on');
set(EPmain.handles.view.marker2,'enable','on');
set(EPmain.handles.view.startSamp,'enable','on');
set(EPmain.handles.view.endSamp,'enable','on');
set(EPmain.handles.view.startHz,'enable','on');
set(EPmain.handles.view.endHz,'enable','on');
set(EPmain.handles.view.topVolt,'enable','on');
set(EPmain.handles.view.bottomVolt,'enable','on');
for color=1:4
    set(EPmain.handles.view.dataset(color),'enable','on');
    if EPmain.view.dataset(color) <= length(EPdataset.dataset)
        set(EPmain.handles.view.cell(color),'enable','on');
        set(EPmain.handles.view.subject(color),'enable','on');
        set(EPmain.handles.view.trial(color),'enable','on');
        set(EPmain.handles.view.factor(color),'enable','on');
    end;
end;

close(EPtopos.handles.topos.topoWindow);
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeMV(src,eventdata)
global EPtopos

plotMVmin=str2num(get(EPtopos.handles.plotMVmin,'string'));
plotMVmax=str2num(get(EPtopos.handles.plotMVmax,'string'));

if plotMVmin < plotMVmax
    EPtopos.plotMVmin=plotMVmin;
    EPtopos.plotMVmax=plotMVmax;
end;

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeLatency(src,eventdata,theRow)
global EPtopos EPmain EPdataset

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(theRow/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+theRow;
end;
theLatency=str2num(get(EPtopos.handles.latency(theRow),'string')); %the input number in ms
if EPtopos.complexData && rem(theRow,2)
    oldTime=EPtopos.firstTime+(EPtopos.iPoints(rowCounter)-1)*EPtopos.sampleSize; %existing number number in ms
else
    oldTime=EPtopos.firstTime+(EPtopos.points(rowCounter)-1)*EPtopos.sampleSize; %existing number number in ms
end
newSample=round((theLatency-EPtopos.firstTime)/EPtopos.sampleSize)+1;
if theLatency ~= oldTime
    if (theLatency >= EPtopos.firstTime) && (theLatency <= EPtopos.lastTime)
        if (newSample > 0) && (newSample <= EPtopos.numPoints)
            if EPtopos.synch
                EPtopos.points(:)=newSample; %update all rows to new number (samples or bins)
                EPtopos.iPoints(:)=newSample; %update all rows to new number (samples or bins)
            else
                if EPtopos.complexData && rem(theRow,2)
                    EPtopos.iPoints(rowCounter)=newSample; %update to new number (samples or bins)
                else
                    EPtopos.points(rowCounter)=newSample; %update to new number (samples or bins)
                end
            end;
        end;
    end;
end;

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function changeHz(src,eventdata,theRow)
global EPtopos EPmain EPdataset

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(theRow/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+theRow;
end;
theHz=get(EPtopos.handles.Hz(theRow),'value');
if EPtopos.synch
    EPtopos.freqs(:)=theHz; %update all rows to new bins
    EPtopos.iFreqs(:)=theHz; %update all rows to new bins
else
    if EPtopos.complexData && rem(theRow,2)
        EPtopos.iFreqs(rowCounter)=theHz; %update to new bins
    else
        EPtopos.freqs(rowCounter)=theHz; %update to new bins
    end
end;

ep_showTopos(EPtopos.firstHz,EPtopos.lastHz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeChannel(src,eventdata,theRow)
global EPtopos

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(theRow/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+theRow;
end;
theChannel=EPtopos.nonRegChans(get(EPtopos.handles.channel(theRow),'value'));
if EPtopos.synch
    EPtopos.chans(:)=theChannel; %update all rows to new chans
    EPtopos.iChans(:)=theChannel; %update all rows to new chans
else
    if EPtopos.complexData && rem(theRow,2)
        EPtopos.iChans(rowCounter)=theChannel;
    else
        EPtopos.chans(rowCounter)=theChannel;
    end
end;

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D3head(src,eventdata)
global EPtopos EPmain

disp('Using EEGlab function headplot to perform 3D head display.');

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end;

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end;

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end;

%check to see if spline file needs to be generated
if ~isempty(EPtopos.ced) && ~any(strcmp(EPtopos.ced,{'none','internal'}))
    [pathstr, name, ext] = fileparts(EPtopos.ced);
    CEDloc=which(EPtopos.ced);
    if isempty(CEDloc)
        name=EPtopos.dataName;
        CEDloc=[pwd filesep 'temp'];
    end;
else
    name=EPtopos.dataName;
    CEDloc=[pwd filesep 'temp'];
end;
if isempty(which([name '.spl']))
    [pathstr2, name2, ext2] = fileparts(CEDloc);
    if max(abs([EPtopos.eloc.X])) <= 1
        MNItransform=[ 0 -15 0 0.08 0 -1.571 102 93 100 ]; %assume eloc coordinates are from a .elp file
    else
        MNItransform=[ 0 -15 4 0.05 0 -1.571 10.2 12 12.2 ]; %assume eloc coordinates are from a .sfp file
    end;
    headplot('setup', EPtopos.eloc(EPtopos.nonRegChans), [pathstr2 filesep name '.spl'],'transform',MNItransform); %save spline file in same place as the ced file.
end;

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end;

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),:,col);
    else
        theData=EPtopos.totalData(EPtopos.chans(rowCounter),EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),:,col);
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),1,col);
    else
        theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),1,col);
    end
end;

theData(isnan(theData))=0;

figure
[hdaxis cbaraxis] = headplot(theData,which([name '.spl']),'cbar',0,'maplimits',[EPtopos.plotMVmin EPtopos.plotMVmax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sources=computeDipoles(topo,eloc,sampleSize)

disp('Using FieldTrip function ft_dipolefitting to perform dipole analysis.');

data=[];
cfg=[];

%adapted from EEGlab's eeglab2fieldtrip under GPL
data.elec.pnt   = zeros(length( eloc ), 3);
for ind = 1:length( eloc )
    data.elec.label{ind} = eloc(ind).labels;
    if ~isempty(eloc(ind).X)
        data.elec.pnt(ind,1) = eloc(ind).X;
        data.elec.pnt(ind,2) = eloc(ind).Y;
        data.elec.pnt(ind,3) = eloc(ind).Z;
    else
        data.elec.pnt(ind,:) = [0 0 0];
    end;
end;

if max(abs([eloc.X])) <= 1
    MNItransform=[ 0 -15 0 0.08 0 -1.571 102 93 100 ]; %assume eloc coordinates are from a .elp file
else
    MNItransform=[ 0 -15 4 0.05 0 -1.571 10.2 12 12.2 ]; %assume eloc coordinates are from a .sfp file
end;
transfmat = traditionaldipfit(MNItransform);
data.elec.pnt = transfmat * [ data.elec.pnt ones(size(data.elec.pnt,1),1) ]';
data.elec.pnt = data.elec.pnt(1:3,:)';
%end adapted section

data.fsample = 1000/(sampleSize);
data.avg  = topo;
%data.trial{1}  = topo;
data.var  = ones(length(topo),1);
data.time(1) = 1;
data.label = { eloc.labels };
data.dimord='rpt_chan_time';
%data.unmixing = diag(ones(length(data.label),1)); % workaround for FieldTrip bug

cfg.numdipoles  = 2;
cfg.symmetry    = 'x';
cfg.channel     = 'all';
cfg.gridsearch  = 'yes';
cfg.nonlinear   = 'yes';
cfg.hdmfile     = which('standard_BESA.mat');
cfg.grid.xgrid  = 'auto';
cfg.grid.ygrid  = 'auto';
cfg.grid.zgrid  = 'auto';
cfg.grid.resolution = 50;
cfg.elec=data.elec;

[grid, cfg] = ft_prepare_sourcemodel(cfg);

[source] = ft_dipolefitting(cfg, data);

sources.momxyz(1,1:3)=source.dip.mom(1:3);
sources.momxyz(2,1:3)=source.dip.mom(4:6);

if (abs(source.dip.pos(1,1)) < .1) && (abs(source.dip.pos(2,1)) < .1)
    disp('Restarting dipole analysis with one dipole - two dipole solution blew up.');
    cfg.numdipoles  = 1;
    cfg.symmetry    = [];
    [source] = ft_dipolefitting(cfg, data);
    sources.momxyz=[];
    sources.momxyz(1,1:3)=source.dip.mom(1:3);
end;

if isfield(source.dip,'rv')
    sources.posxyz=source.dip.pos;
    sources.rv=source.dip.rv;
else
    msg{1}='Error: Source analysis has failed.';
    [msg]=ep_errorMsg(msg);
    return
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dipoles(src,eventdata)
global EPtopos EPmain

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end;

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end;

if strcmp('FFT',EPmain.view.dataTransform) && ~EPmain.view.rel(col)
    msg{1}='Error: Dipole analysis only works with voltage data.';
    [msg]=ep_errorMsg(msg);
    return
end;

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end;

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end;

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),EPtopos.iPoints(rowCounter),theRow,EPtopos.freqs(rowCounter),:,col);
    else
        theData=EPtopos.totalData(EPtopos.chans(rowCounter),EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),:,col);
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.freqs(rowCounter),1,col);
    else
        theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),1,col);
    end
end;

theData=squeeze(theData);
theData(isnan(theData))=0;

eloc=EPtopos.eloc(EPtopos.nonRegChans);

goodChans=find(~isnan(theData));

sources=computeDipoles(theData(goodChans),eloc(goodChans),EPtopos.sampleSize);

dipplot(sources,'mri',which('avg152t1.mat'),'color', { 'g' 'b' },'coordformat','MNI','verbose','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jackknife(src,eventdata)
global EPtopos

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end;

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end;

row=row+EPtopos.perPage*(EPtopos.page-1);

if ~isfield(EPtopos,'jack')
    msg{1}='Error: Jack-knife dipole analysis only works with PCA dataset.';
    [msg]=ep_errorMsg(msg);
    return
end;

if ~strcmp(EPtopos.jack(col).PCAmode,'spat')
    msg{1}='Error: Jack-knife dipole analysis only works with spatial or temporo-spatial PCA dataset.';
    [msg]=ep_errorMsg(msg);
    return
end;

eloc=EPtopos.eloc;
topo=squeeze(EPtopos.totalData(:,EPtopos.points(row),row,col));
sources=computeDipoles(topo,eloc,EPtopos.sampleSize);
dipColors{1}='r';
for theJK=1:size(EPtopos.jack(col).FacPat,3)
    topo=EPtopos.jack(col).FacPat(:,row,theJK);
    sources(end+1)=computeDipoles(topo,eloc,EPtopos.sampleSize);
    dipColors{end+1}='b';
end;

EPtopos.jack(col).sources=sources;

dipplot(sources,'mri',which('avg152t1.mat'),'color', dipColors,'coordformat','MNI','verbose','off','spheres','on','dipolelength',0,'summary','3d');

%jack-knife t-test for hemispheric main effect
%using statistical test presented by: Miller, J., Patterson, T., & Ulrich, R. (1998). Jackknife-based method for measuring LRP onset latency differences. Psychophysiology, 35(1), 99-115.

numJK=length(sources);
DipAmp=zeros(numJK,2);

numJKpairs=0;
for JK=1:numJK
    momxyz=EPtopos.jack(col).sources(JK).momxyz;
    if size(momxyz,1) == 2
        numJKpairs=numJKpairs+1;
        DipAmp(numJKpairs,1)=sqrt(momxyz(1,1)^2+momxyz(1,2)^2+momxyz(1,3)^2);
        DipAmp(numJKpairs,2)=sqrt(momxyz(2,1)^2+momxyz(2,2)^2+momxyz(2,3)^2);
    end;
end;

JD=DipAmp(1:numJKpairs,1)-DipAmp(1:numJKpairs,2);
JM=mean(JD);

SDjk=sqrt(((numJKpairs-1)/numJKpairs)*sum((JD-JM).^2));

if numJKpairs > 1
    disp(['The t-value for the amplitude of the two hemispheric dipoles is ' num2str(sum(JD)/SDjk) ' with ' num2str(numJKpairs-1) ' degrees of freedom.']);
    if sum(JD) > 0
        disp('Right hemisphere larger than left.');
    elseif sum(JD) < 0
        disp('Left hemisphere larger than right.');
    else
        disp('No difference between hemispheres.');
    end;
    if ft_hastoolbox('STATS', 0, 1)
        disp(['It has a two-sided p-value of ' num2str(2*tcdf(-abs(sum(JD)/SDjk),numJKpairs-1)) '.']);
    end;
else
    disp('Too few solutions with hemispheric dipoles to calculate hemispheric t-test.');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rescaleFigures(src,eventdata)
global EPtopos EPmain

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end;

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end;

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end;

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end;

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),:,col);
    else
        theData=EPtopos.totalData(EPtopos.chans(rowCounter),EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),:,col);
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),1,col);
    else
        theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),1,col);
    end
end;

if ~any(theData)
    disp('There is no data to rescale to in this topoplot.')
    return
end;

EPtopos.plotMVmin=min(min(min(theData)));
EPtopos.plotMVmax=max(max(max(theData)));

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%explore effects of reference scheme on the data
function referenceChan(src,eventdata)
global EPtopos EPmain

scrsz = EPmain.scrsz;

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end;

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end;

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end;

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end;

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),:,col);
    else
        theData=EPtopos.totalData(EPtopos.chans(rowCounter),EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),:,col);
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),1,col);
    else
        theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),1,col);
    end
end;

if ~any(theData)
    disp('There is no data to depict to in this topoplot.')
    return
end;

EPtopos.referenceFigure.refchan1=1;
EPtopos.referenceFigure.refchan2=1;
switch EPtopos.reference(col).type
    case 'REG'
        switch length(EPtopos.reference(col).current)
            case 0
                EPtopos.referenceFigure.reference=5;
            case 1
                EPtopos.referenceFigure.reference=1;
                EPtopos.referenceFigure.refchan1=EPtopos.reference(col).current(1);
                EPtopos.referenceFigure.refchan2=1;
            case 2
                EPtopos.referenceFigure.reference=2;
                EPtopos.referenceFigure.refchan1=EPtopos.reference(col).current(1);
                EPtopos.referenceFigure.refchan2=EPtopos.reference(col).current(2);
            otherwise
                EPtopos.referenceFigure.reference=5;
        end;
    case 'AVG'
        EPtopos.referenceFigure.reference=3;
    case 'PAR'
        EPtopos.referenceFigure.reference=4;
    otherwise
        EPtopos.referenceFigure.reference=5;
end;

EPtopos.handles.referenceFigure.figure = figure('Name', 'Explore effect of rereferencing', 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure', 'Position',[scrsz(3)/2 scrsz(4)/2 600 400]);
colormap jet;

uicontrol('Style','text',...
    'String','Reference Scheme','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[5 70 105 20]);

EPtopos.handles.referenceFigure.reference= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',{'One Chan','Two Chans','Average','PARE','unknown'},...
    'CallBack',@toposRereference,...
    'Value',EPtopos.referenceFigure.reference,'Position',[10 50 100 20]);

EPtopos.handles.referenceFigure.refchan1= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',EPtopos.chanNames(EPtopos.nonRegChans),...
    'CallBack',@toposRereference,...
    'Value',EPtopos.referenceFigure.refchan1,'Position',[10 30 100 20]);

EPtopos.handles.referenceFigure.refchan2= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',EPtopos.chanNames(EPtopos.nonRegChans),...
    'CallBack',@toposRereference,...
    'Value',EPtopos.referenceFigure.refchan2,'Position',[10 10 100 20]);

if ~any(ismember([1 2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan1,'enable','off');
end;
if ~any(ismember([2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan2,'enable','off');
end;

EPtopos.handles.referenceFigure.axes=axes('Units','pixel','position',[200 20 360 360]);

nSides=20;
[sphereCoords, elecInSphere]=ep_sphereHead(nSides, EPtopos.eloc);
sphereValues=ep_interpolateHead(theData, elecInSphere, sphereCoords);
X=reshape(sphereCoords(:,1),nSides+1,nSides+1);
Y=reshape(sphereCoords(:,2),nSides+1,nSides+1);
Z=reshape(sphereCoords(:,3),nSides+1,nSides+1);
sphereValues=reshape(sphereValues,nSides+1,nSides+1);
EPtopos.handles.referenceFigure.sphere=surf(X,Y,Z,sphereValues);
set(EPtopos.handles.referenceFigure.sphere,'CDataSource','EPtopos.referenceFigure.sphereValues');
set(EPtopos.handles.referenceFigure.sphere,'FaceColor','interp');
theMax=max(abs([min(theData) max(theData)]));
set(EPtopos.handles.referenceFigure.axes,'CLim',[-theMax theMax]);
cmap=colormap(EPtopos.handles.referenceFigure.figure);
cmap(29:36,:)=1;
set(EPtopos.handles.referenceFigure.figure,'colormap',cmap);

EPtopos.referenceFigure.elecValues=theData;
EPtopos.referenceFigure.sphereValues=sphereValues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%change the reference on the Topos rereference panel
function toposRereference(src,eventdata)
global EPtopos

EPtopos.referenceFigure.reference=get(EPtopos.handles.referenceFigure.reference,'Value');
EPtopos.referenceFigure.refchan1=get(EPtopos.handles.referenceFigure.refchan1,'Value');
EPtopos.referenceFigure.refchan2=get(EPtopos.handles.referenceFigure.refchan2,'Value');

if ~any(ismember([1 2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan1,'enable','off');
else
    set(EPtopos.handles.referenceFigure.refchan1,'enable','on');
end;
if ~any(ismember([2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan2,'enable','off');
else
    set(EPtopos.handles.referenceFigure.refchan2,'enable','on');
end;

switch EPtopos.referenceFigure.reference
    case 1
        refValue=EPtopos.referenceFigure.elecValues(EPtopos.referenceFigure.refchan1);
    case 2
        refValue=(EPtopos.referenceFigure.elecValues(EPtopos.referenceFigure.refchan1)+EPtopos.referenceFigure.elecValues(EPtopos.referenceFigure.refchan2))/2;
    case 3
        refValue=mean(EPtopos.referenceFigure.elecValues);
    case 4
        %first average reference
        refValue=mean(EPtopos.referenceFigure.elecValues);
        EPtopos.referenceFigure.elecValues=EPtopos.referenceFigure.elecValues-refValue;
        
        %then interpolate surface
        theData=EPtopos.referenceFigure.elecValues;
        nSides=20;
        [sphereCoords, elecInSphere]=ep_sphereHead(nSides, EPtopos.eloc);
        sphereValues=ep_interpolateHead(theData, elecInSphere, sphereCoords);
        X=reshape(sphereCoords(:,1),nSides+1,nSides+1);
        Y=reshape(sphereCoords(:,2),nSides+1,nSides+1);
        Z=reshape(sphereCoords(:,3),nSides+1,nSides+1);
        sphereValues=reshape(sphereValues,nSides+1,nSides+1);
        EPtopos.handles.referenceFigure.sphere=surf(X,Y,Z,sphereValues);
        set(EPtopos.handles.referenceFigure.sphere,'CDataSource','EPtopos.referenceFigure.sphereValues');
        set(EPtopos.handles.referenceFigure.sphere,'FaceColor','interp');
        theMax=max(abs([min(theData) max(theData)]));
        set(EPtopos.handles.referenceFigure.axes,'CLim',[-theMax theMax]);
        cmap=colormap(EPtopos.handles.referenceFigure.figure);
        cmap(29:36,:)=1;
        set(EPtopos.handles.referenceFigure.figure,'colormap',cmap);
        EPtopos.referenceFigure.sphereValues=sphereValues;
        
        %then sum up the surface values
        sphereValues=EPtopos.referenceFigure.sphereValues(2:end-1,2:end-1);
        sphereValues=[sphereValues(:); EPtopos.referenceFigure.sphereValues(1,1); EPtopos.referenceFigure.sphereValues(end,end)];
        refValue=mean(sphereValues);
    case 5
        refValue=0;
    otherwise
        refValue=0;
end;
EPtopos.referenceFigure.elecValues=EPtopos.referenceFigure.elecValues-refValue;
EPtopos.referenceFigure.sphereValues=EPtopos.referenceFigure.sphereValues-refValue;
refreshdata(EPtopos.handles.referenceFigure.sphere,'caller');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display two channels and plot sample-by-sample t-test comparison
function twoChan(src,eventdata)
global EPtopos EPmain EPdataset

if isempty(EPtopos.handles.twoChan.figure) || ~isgraphics(EPtopos.handles.twoChan.figure)
    scrsz = EPmain.scrsz;
    
    EPtopos.twoChan.sigTest=0;
    
    [EPtopos.twoChan.row,EPtopos.twoChan.col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));
    
    if isempty(EPtopos.twoChan.row)
        [EPtopos.twoChan.row,EPtopos.twoChan.col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
    end;
    
    if isempty(EPtopos.twoChan.row)
        tempVar=get(gco);
        [EPtopos.twoChan.row,EPtopos.twoChan.col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
    end;
    
    if EPtopos.complexData
        EPtopos.twoChan.rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(EPtopos.twoChan.row/2);
    else
        EPtopos.twoChan.rowCounter=(EPtopos.page-1)*EPtopos.perPage+EPtopos.twoChan.row;
    end;
    
    if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
        theRow=EPtopos.twoChan.rowCounter;
    else
        theRow=1;
    end;
    
    if EPmain.view.rel(EPtopos.twoChan.col) %if relational data
        if EPtopos.complexData && rem(EPtopos.twoChan.row,2)
            EPtopos.twoChan.theData=EPtopos.totalImagData(EPtopos.iChans(EPtopos.twoChan.rowCounter),EPtopos.iPoints(EPtopos.twoChan.rowCounter),theRow,EPtopos.iFreqs(EPtopos.twoChan.rowCounter),:,EPtopos.twoChan.col);
        else
            EPtopos.twoChan.theData=EPtopos.totalData(EPtopos.chans(EPtopos.twoChan.rowCounter),EPtopos.points(EPtopos.twoChan.rowCounter),theRow,EPtopos.freqs(EPtopos.twoChan.rowCounter),:,EPtopos.twoChan.col);
        end
    else
        if EPtopos.complexData && rem(EPtopos.twoChan.row,2)
            if isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).timeNames) && ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
                EPtopos.twoChan.theData=EPtopos.totalImagData(:,EPtopos.iPoints(EPtopos.twoChan.rowCounter),theRow,:,1,EPtopos.twoChan.col);
            else
                EPtopos.twoChan.theData=EPtopos.totalImagData(:,EPtopos.iPoints(EPtopos.twoChan.rowCounter),theRow,EPtopos.iFreqs(EPtopos.twoChan.rowCounter),1,EPtopos.twoChan.col);
            end;
        else
            if isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).timeNames) && ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
                EPtopos.twoChan.theData=EPtopos.totalData(:,:,theRow,:,1,EPtopos.twoChan.col);
            else
                EPtopos.twoChan.theData=EPtopos.totalData(:,:,theRow,EPtopos.freqs(EPtopos.twoChan.rowCounter),1,EPtopos.twoChan.col);
            end;
        end
    end;
    
    if ~any(EPtopos.twoChan.theData)
        disp('There is no data to depict to in this topoplot.')
        return
    end;
    
    EPtopos.handles.twoChan.figure = figure('Name', 'Compare two channels', 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure', 'Position',[scrsz(3)/2 scrsz(4)/2 600 400]);
    
    EPtopos.twoChan.chan1=EPtopos.chans(EPtopos.twoChan.rowCounter);
    eDist=zeros(length(EPtopos.nonRegChans),1);
    for iChan=1:length(EPtopos.nonRegChans)
        theChan=EPtopos.nonRegChans(iChan);
        if theChan==EPtopos.twoChan.chan1
            eDist(iChan)=inf;
        else
            %Y is left-right
            eDist(iChan)=norm([EPtopos.eloc(EPtopos.twoChan.chan1).X-EPtopos.eloc(theChan).X -EPtopos.eloc(EPtopos.twoChan.chan1).Y-EPtopos.eloc(theChan).Y EPtopos.eloc(EPtopos.twoChan.chan1).Z-EPtopos.eloc(theChan).Z]);
        end;
        
    end;
    [X I]=min(eDist);
    EPtopos.twoChan.chan2=EPtopos.nonRegChans(I);
    
    if isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).timeNames) && ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
        EPtopos.twoChan.plotForm='FFT';
    else
        EPtopos.twoChan.plotForm='VLT';
    end;
    
    uicontrol('Style','text',...
        'String','Channel 1','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'ForegroundColor',EPtopos.RGBcolors(1,:),'Position',[15 380 50 20]);
    
    EPtopos.handles.twoChan.chan1=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPtopos.chanNames(EPtopos.nonRegChans),...
        'Value',EPtopos.twoChan.chan1,...
        'Position',[70 380 100 20],...
        'CallBack',@twoChan);

    
    uicontrol('Style','text',...
        'String','Channel 2','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'ForegroundColor',EPtopos.RGBcolors(3,:),'Position',[175 380 50 20]);
    
    EPtopos.handles.twoChan.chan2=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPtopos.chanNames(EPtopos.nonRegChans),...
        'Value',EPtopos.twoChan.chan2,...
        'Position',[230 380 100 20],...
        'CallBack',@twoChan);
    
    EPtopos.handles.twoChan.sigTest=uicontrol('Style','checkbox',...
        'String','Sample-by-sample t-test','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'Value',EPtopos.twoChan.sigTest,'Position',[335 380 150 20],'CallBack',@twoChan);
    
    if strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'continuous') ||...
        (strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'single_trial') && length(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).trialNames)==1) ||...
        (strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'average') && length(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).subNames)==1)
            set(EPtopos.handles.twoChan.sigTest,'enable','off');
    end;
end;

EPtopos.twoChan.chan1=EPtopos.nonRegChans(get(EPtopos.handles.twoChan.chan1,'Value'));
EPtopos.twoChan.chan2=EPtopos.nonRegChans(get(EPtopos.handles.twoChan.chan2,'Value'));

EPtopos.twoChan.sigTest=get(EPtopos.handles.twoChan.sigTest,'Value');
if EPtopos.twoChan.sigTest
    
    disp('Will need to load the dataset to perform the sample-by-sample statistical test so there will be a delay.');
    
    cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(EPtopos.twoChan.col)}(EPmain.view.cell(EPtopos.twoChan.col)),EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).cellNames));
    if strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'average')
        subList=find(strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).subTypes,'RAW'));
    else
        subList=EPmain.view.subject(EPtopos.twoChan.col);
    end;
    if strcmp(EPtopos.plotForm,'FFT')
        freqList=[];
    elseif ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
        freqList=EPtopos.freqs(EPtopos.twoChan.col);
    else
        freqList=[];
    end;
    if strcmp(EPtopos.plotForm,'factor')
        facList=EPtopos.twoChan.rowCounter;
    else
        facList=EPmain.view.factor(EPtopos.twoChan.col);
    end;
    
    EPdata=ep_loadEPdataset(EPdataset,EPmain.view.dataset(EPtopos.twoChan.col));
    EPdataIn=ep_selectData(ep_stripAdds(EPdata),{[EPtopos.twoChan.chan1 EPtopos.twoChan.chan2],[],cellList,subList,facList,freqList});
    [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdataIn,'sample',.05,4,[],[],[],'t-test',[],[]);
    outputData=squeeze(outputData);
    outputData=outputData(:);
end;

%single-trial, average, FFT, TFT

waveSize=.8;

switch EPtopos.twoChan.plotForm
    case 'VLT'
        numImages=length(find(EPmain.view.allTrials == 2));
        if (length(EPtopos.plotColors)-numImages) > 0
            imageSpace=4;
        else
            imageSpace=numImages;
        end;
        EPtopos.twoChan.handles.waves.hExpandedAxes=[];
        if numImages %if any erpimage
            imageCount=0;
            for i=1:4
                if (EPmain.view.dataset(i) <= length(EPdataset.dataset)) && (EPmain.view.allTrials(i) == 2)
                    imageCount=imageCount+1;
                    trialList=find(EPtopos.colorIndex==i);
                    EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .10+(waveSize/imageSpace)*(imageSpace-imageCount) waveSize (waveSize/imageSpace)]);
                    EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = imagesc(EPtopos.firstTime:EPtopos.lastTime,1:length(trialList),squeeze(EPtopos.twoChan.theData([EPtopos.twoChan.chan1 EPtopos.twoChan.chan2],:,trialList,:))',[EPtopos.plotMVmin, EPtopos.plotMVmax]);
                    axis([EPtopos.firstTime EPtopos.lastTime 1 length(trialList)]);
                    line([0 0],[1 length(trialList)],'Color','black','LineWidth',1) %stimulus onset
                    if ~isempty(EPtopos.marker1)
                        line(repmat(EPtopos.marker1,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                    if ~isempty(EPtopos.marker2)
                        line(repmat(EPtopos.marker2,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                end;
            end;
        end
        if (length(EPtopos.plotColors)-numImages) > 0 %if there will be waveforms
            EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .10 waveSize waveSize*((4-numImages)/4)]);
            hold on
            if EPtopos.twoChan.sigTest
                breakList=sort([find(diff([0; (outputData>0); 0])<0)-1; find(diff([0; (outputData>0); 0])>0)]);
                if ~isempty(breakList)
                    theData1=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,:,:));
                    theData2=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,:,:));
                    theData1=theData1(:);
                    theData2=theData2(:);
                    for iSigArea=1:length(breakList)/2
                        theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)]';
                        if length(theTimePoints) == 1
                            EPtopos.twoChan.handles.waves.hLines=line(([theTimePoints theTimePoints]*EPtopos.spacing)+(EPtopos.firstTime-EPtopos.spacing),[theData1(theTimePoints) theData2(theTimePoints)],'LineWidth',1,'Color',[1 .5 .5]);
                        else
                            EPtopos.twoChan.handles.waves.hLines=patch(([theTimePoints; flip(theTimePoints)]*EPtopos.spacing)+(EPtopos.firstTime-EPtopos.spacing),[theData1(theTimePoints); theData2(flip(theTimePoints))],EPtopos.plotColorIndex(EPtopos.twoChan.col,:),'FaceColor','red','EdgeColor','none','FaceAlpha',.25);
                        end;
                    end;
                end;
            end;
            plot([EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPtopos.RGBcolors(1,:));
            plot([EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPtopos.RGBcolors(3,:));
            hold off
            axis([EPtopos.firstTime EPtopos.lastTime EPtopos.plotMVmin EPtopos.plotMVmax]);
            if EPtopos.direction ==2
                set(EPtopos.twoChan.handles.waves.hWave([EPtopos.twoChan.chan1 EPtopos.twoChan.chan2]),'YDir','reverse')
            end;
            line([EPtopos.firstTime EPtopos.lastTime-EPtopos.sampleSize],[0 0],'Color','black','LineWidth',1) % zero line
            line([0 0],[0 EPtopos.plotMVmax],'Color','black','LineWidth',1) %stimulus onset
        end;
        
        if ~isempty(EPtopos.marker1)
            try
                eval('line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1)');
            catch
            end
        end;
        if ~isempty(EPtopos.marker2)
            try
                eval('line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1);');
            catch
            end;
        end;
    case 'FFT'
        EPtopos.twoChan.handles.waves.hExpandedAxes=[];
        EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .05 waveSize waveSize]);
        hold on
        if EPtopos.twoChan.sigTest
            breakList=sort([find(diff([0; (outputData>0); 0])<0)-1; find(diff([0; (outputData>0); 0])>0)]);
            if ~isempty(breakList)
                theData1=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,:,:));
                theData2=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,:,:));
                theData1=theData1(:);
                theData2=theData2(:);
                for iSigArea=1:length(breakList)/2
                    theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)]';
                    if length(theTimePoints) == 1
                        EPtopos.twoChan.handles.waves.hLines=line(([theTimePoints theTimePoints]*EPtopos.spacing)+(EPtopos.firstHz-EPtopos.spacing),[theData1(theTimePoints) theData2(theTimePoints)],'LineWidth',1,'Color',[1 .5 .5]);
                    else
                        EPtopos.twoChan.handles.waves.hLines=patch(([theTimePoints; flip(theTimePoints)]*EPtopos.spacing)+(EPtopos.firstHz-EPtopos.spacing),[theData1(theTimePoints); theData2(flip(theTimePoints))],EPtopos.plotColorIndex(EPtopos.twoChan.col,:),'FaceColor','red','EdgeColor','none','FaceAlpha',.25);
                    end;
                end;
            end;
        end;
        plot([EPtopos.firstHz+EPtopos.sampleSize:EPtopos.spacing:EPtopos.lastHz],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPtopos.RGBcolors(1,:));
        plot([EPtopos.firstHz+EPtopos.sampleSize:EPtopos.spacing:EPtopos.lastHz],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPtopos.RGBcolors(3,:));
        hold off
        axis([EPtopos.firstHz EPtopos.lastHz EPtopos.plotMVmin EPtopos.plotMVmax]);
        line([EPtopos.firstHz+EPtopos.sampleSize EPtopos.lastHz],[0 0],'Color','black','LineWidth',1) % zero line
        if ~isempty(EPtopos.marker1)
            try
                eval('line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1)');
            catch
            end
        end;
        if ~isempty(EPtopos.marker2)
            try
                eval('line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',1);');
            catch
            end;
        end;
end;

