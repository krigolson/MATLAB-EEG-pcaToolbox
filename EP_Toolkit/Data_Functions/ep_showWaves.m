function err=ep_showWaves(plotMVmin,plotMVmax,startTime,endTime,startHz,endHz,marker1,marker2,FFTunits,eventLines,showLegend)
% ep_showWaves - ep_showWaves(plotMVmin,plotMVmax,startTime,endTime,startHz,endHz,marker1,marker2,FFTunits,eventLines,showLegend) -
% Displays the waveforms of the data.
%
%Input:
%  plotMVmin         : Minimum of voltage scale.  For spectral data, already converted to psd and dB as necessary.
%  plotMVmax         : Maximum of voltage scale  For spectral data, already converted to psd and dB as necessary.
%  startTime         : Start of epoch in ms
%  endTime           : End of epoch in ms (right side unless both startTime and endTime are the same number)
%  startHz           : Start of Hz band
%  endHz             : End of Hz band
%  marker1           : Optional marker for plots
%  marker2           : Second optional marker for plots
%  FFTunits          : spectral data units (1=comp, 2=asd, 3=psd, 4=dB)
%  eventLines        : Cell array of samples of events to be displayed (4 colors)
%  showLegend        : Show the legend box (1='yes',0='no')
%
%Output:
%  err             : Returns 0 if no problems and 1 if there was an error.

%  Time is counted for an epoch as, say, -200 to 800 ms.  In this case the samples have some length to them depending on
% the digitization rate, such as 4ms for a 250Hz sampling rate.  The points prior to the baseline are counted from the
% left side of the baseline sample and the points following the baseline are counted from the right side of the baseline
% sample.  Thus, when the actual samples are calculated, one first adds the baseline to the numbers, yielding 0 to
% 1000ms.  One then must add the length of the samples to the first number to obtain 4 to 1000ms.  Then one would
% divide the numbers by the sample size, resulting in 1 to 250 samples, which are the correct samples to use.
% For continuous data or other such data where there is no baseline, one must by this convention start with 0 ms,
% as in 0-1000ms and 1000-2000ms.
%
% For error bars:
% Payton, M. E., Greenstone, M. H., & Schenker, N. (2003). Overlapping confidence intervals or standard error intervals: what do they mean in terms of statistical significance. Journal of Insect Science, 3(1), 34. 

%History
%  by Joseph Dien (5/15/09)
%  jdien07@mac.com
%
%  bugfix 8/27/09 JD
%  Location of windows appearing partly off screen on some systems.  Fixed.
%
%  modified 1/16/10 JD
%  Rounds off ms labels in the scale graph.
%
%  modified 1/26/10 JD
%  Sped up waveform drawing by relying more on EPwork cache.
%  Added ability to plot data from continuous data files by showing only one second at a time.
%  Epoch ms now displayed as from beginning of first sample to end of last sample.
%
%  modified 2/16/10 JD
%  Can now handle displaying waves from multiple datasets where some of the channels are implicit in one but not the
%  other dataset.
%  Eliminated chantype field for implicit channels.
%
%  bugfix 2/27/10 JD
%  Fixed crash when viewing data with regional channels.
%
%  bugfix 3/16/10 JD
%  Fixed crash when viewing data with no implicit channels.
%
%  bugfix 3/23/10 JD
%  Fixed crash when viewing data with no electrode coordinate information.
%
%  modified 3/28/10 JD
%  Reduced memory use of figures by odd workaround of dividing the data by an arbitrary scaling factor.
%
%  bugfix 3/30/10 JD
%  Fixed baseline bar not appearing in expanded waveform figures.
%  Fixed stimulus onset bar not appearing in waveforms figures.
%  Fixed baseline bar and stimulus onset bar not appearing in expanded waveform figures when click lands on a waveform.
%  Fixed real reason figures using much more memory than needed.
%
%  modified 4/23/10 JD
%  Scale of waveforms no longer constrained to be at least +/- 1 microvolt.
%
% modified 3/15/12 JD
% Added support for plotting FFT data, including psd and dB scaling.
%
%  bugfix 1/10/13 JD
%  Fixed missing markers when data range is less than one.
%
%  bugfix 1/13/13 JD
%  Fixed error message when expanding channel for waveform plot for which the data are no longer available.
%
% modified 1/16/13 JD
% Handles situation where FFT data has negative values (e.g., due to imprecision in PCA results) and transforming to dB
% would result in complex numbers, by taking absolute value first.
%
% bugfix 2/5/13 JD
% Fixed crash when all channels are along midline or center line.
%
% modified 4/1/13 JD
% Wave plots can accommodate datasets with different sets of channels.
%
% bugfix 4/2/13 JD
% Markers in waveform plots can be set at zero ms.
%
% bugfix 5/9/13 JD
% Fixed View crashing if the scaling is from zero to zero, as with a bad trial.
%
% bugfix 8/4/13 JD
% Fixed View crashing if none of the data have electrode coordinates available.
%
% bugfix 10/15/13 JD
% Fixed epoch one sample longer than intended for continuous data.
%
% modified 11/3/13 JD
% For continuous data, baseline correct each channel by entire one second epoch so that waves will be visible.
% Channel numbers are now positioned in relative units so not moved when scaling changes.
%
% bugfix 11/6/13 JD
% Fixed crash in view waves function when superimposing two datasets where one is missing channel coordinates.
%
% modified 11/6/13 JD
% Scan and Waves functions can now present event markings.
%
% bugfix 1/12/14 JD
% Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% bugfix 2/25/14 JD
% Fixed crash for viewing TFT data.
%
% modified 2/26/14 JD
% Added View function option to plot or erpimage all trials and all subjects.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% modified 3/20/14 JD
% The events 'SESS','CELL','TRSP','bgin' are not excluded from being displayed for continuous files.
%
% bugfix 4/8/14 JD
% Fixed crash when overplotting two sets of data where one has regional channels not present in the other.
%
% modified 4/9/14 JD
% -all- and -erpimage- options in View leave out adds if lower levels are present (e.g., for subjects, leaves out grand averages if individual subjects are present).
%
% modified 4/24/14 JD
% Added coherence and phase-locking options, including support for complex numbers.
%
% modified 5/11/14 JD
% Added legend to wave plots.
%
% modified 6/19/14 JD
% Added support for sample-by-sample t-tests, including STS chanType.
%
% bugfix 12/2/14 JD
% Fixed channels left out when the last dataset has fewer channels
% (including regional) than earlier datasets.
%
% bufix 1/9/15 JD
% Fixed axis labels being shown for waveform plots under Matlab 2014b due
% to Matlab bug.
%
%  modified 5/25/14 JD
%  Set colormap to jet even for Matlab 2014b onwards.
%
% bugfix 8/30/15 JD
% Fixed amplitude spectral density calculated as divided by Hz rather than
% by square root of Hz.
% Fixed dB of amplitude data not converted to power first.
% Fixed dB units should be labeled as dBV since it is a ratio and therefore
% has no units.
%
% bugfix 1/14/16 JD
% Fixed crash when minimum voltage equals the maximum voltage.
%
% bugfix 1/15/16 JD
% Fixed when clicking on Views channels, can only get expanded channel window if one clicks along the very top edge for FFT data.
% Now allows power scaled data to be displayed as amplitudes.
% Now handles complex FFT numbers.
%
% bugfix & modified 1/22/16 JD
% Fixed Matlab crash when only erpimages are being displayed.
% Consolidated spectral unit controls so just four options (cmp, asd, psd, dB).
% Fixed crash when unable to generate unique labels for the four datasets.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% bugfix 10/9/16 JD
% Fixed FFT data on amplitude scaling not converted to real numbers.
%
% bugfix & modified 6/19/17 JD
% Fixed displaying separate imaginary component of spectral data when the FFT units are not set for complex units.
% For TFT data, when set to display only one Hz bin, switches to waveform display rather than erpimage display.
% Supports display of flexible segments.
% Only EEG chans are subjected to FFT unit conversions.
% Fixed single-sample duration sampleTest results not displaying in View Waves.
% Fixed conversion to spectral density dividing by bin width rather than sqrt(bin width).
% Now allows just one sample of TFT data to be chosen for display using the View function.
% Improved default auto-scaling for View Waves and View Topos.
%
% modified & bugfix 2/9/18 JD
% Added support for GFP plots and error bands.
% Fixed crash when displaying -all- trials.
%
% modified & bugfix 2/23/18 JD
% Fixed not using subNum for grand averages.
% Implemented Cousineau-Morey confidence interval.
% Added -all- and -erpimage- options for cells in View Waves function.
% Fixed crash when only erpimages are being presented.
%
% modified 6/5/18 JD
% Added option to add marks for RT and selected events to View figures.
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

err=0;

global EPdataset EPmain EPwaves

EPwaves.marker1=marker1;
EPwaves.marker2=marker2;
EPwaves.eventLines=eventLines;

scrsz = EPmain.scrsz;

EPwaves.RGBcolors=[0 0 1;1 0 0;0 1 0;0 0 0]; %colors for the waveforms

EPwaves.handles.waves=[];

EPwaves.handles.waves.hWaveWindow = figure('Name', 'WaveWindow', 'NumberTitle', 'off', 'Position',[270 1 scrsz(3)-270 scrsz(4)]);
colormap jet;
%, 'MenuBar', 'none'
EPwaves.plotWidth=.05;
EPwaves.plotHeight=.04;
EPwaves.margin=.1;
noStatsWarn=0;

%make sure the plotted datasets are compatible
numChans=0;
for iColor=1:4
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset) %if not set to "none"
        if numChans == 0
            if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average')
                if EPmain.view.allTrials(iColor)==1
                    subName='all';
                elseif EPmain.view.allTrials(iColor)==2
                    subName='erpimage';
                else
                    subName=EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)};
                end;
            else
                subName=EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)};
            end;
            
            windowName=[EPdataset.dataset(EPmain.view.dataset(iColor)).dataName '-' subName];
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                if EPmain.view.allTrials(iColor)==3
                    facName='all';
                elseif EPmain.view.allTrials(iColor)==4
                    facName='erpimage';
                else
                    facName=EPdataset.dataset(EPmain.view.dataset(iColor)).facNames{EPmain.view.factor(iColor)};
                end;
                windowName=[windowName '-' facName];
            end;
            set(EPwaves.handles.waves.hWaveWindow,'Name',windowName);
            EPwaves.firstTime=min([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]);
            EPwaves.lastTime=max([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]);
            EPwaves.firstHz=min([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]);
            EPwaves.lastHz=max([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]);
            if EPdataset.dataset(EPmain.view.dataset(iColor)).Fs==1
                EPwaves.sampleSize=1;
            else
                EPwaves.sampleSize=1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs;
            end;
            
            implicit=EPdataset.dataset(EPmain.view.dataset(iColor)).implicit; %use implicits of the first dataset
            EPwaves.chanNames=EPdataset.dataset(EPmain.view.dataset(iColor)).chanNames; %use channel names of the first dataset
            eloc=EPdataset.dataset(EPmain.view.dataset(iColor)).eloc;
            if ~isempty(implicit)
                implicit=implicit(find(~strcmp('FID',{implicit.type})));
                EPwaves.chanNames=[EPwaves.chanNames {implicit.labels}]; %use channel names of the first dataset
                eloc=[eloc implicit];
            end;
            numChans=length(EPwaves.chanNames);
            EEGchans=find(ismember(EPdataset.dataset(EPmain.view.dataset(iColor)).chanTypes,{'EEG','REG'}));
            EPwaves.chanIX{iColor}=[1:numChans];
            baseline=EPdataset.dataset(EPmain.view.dataset(iColor)).baseline;
        else
            if baseline ~=EPdataset.dataset(EPmain.view.dataset(iColor)).baseline
                disp('Note: Baseline differs between the datasets.  If a baseline is incorrectly specified, it could result in misalignment of the waveforms.');
            end;
            
            newImplicit=EPdataset.dataset(EPmain.view.dataset(iColor)).implicit;
            if ~isempty(newImplicit)
                newImplicit=newImplicit(find(~strcmp('FID',{newImplicit.type})));
            end;
            newChans=setdiff(EPdataset.dataset(EPmain.view.dataset(iColor)).chanNames,EPwaves.chanNames);
            newEloc=[EPdataset.dataset(EPmain.view.dataset(iColor)).eloc newImplicit];
            if ~isempty(newChans)
                EPwaves.chanNames=[EPwaves.chanNames; newChans]; %add new channels to the list
                for i=1:length(newChans)
                    eloc=[eloc, EPdataset.dataset(EPmain.view.dataset(iColor)).eloc(find(strcmp(newChans(i),EPdataset.dataset(EPmain.view.dataset(iColor)).chanNames)))]; %add new channels to the electrode coordinate list
                end;
                numChans=length(EPwaves.chanNames);
            end;
            EPwaves.firstTime =max(min([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]),EPwaves.firstTime);
            EPwaves.lastTime = min(max([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]),EPwaves.lastTime);
            EPwaves.firstHz =max(min([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]),EPwaves.firstHz);
            EPwaves.lastHz = min(max([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]),EPwaves.lastHz);
            
            if (EPwaves.sampleSize ~= 1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs) && ~((EPwaves.sampleSize ==1) && (EPdataset.dataset(EPmain.view.dataset(iColor)).Fs==1))
                msg{1}='Error: The sampling rates are not compatible.';
                [msg]=ep_errorMsg(msg);
                close(EPwaves.handles.waves.hWaveWindow);
                err=1;
                return
            end;
            if ~isempty(eloc) && ~isempty(newEloc)
                eLabels=EPdataset.dataset(EPmain.view.dataset(iColor)).chanNames;
                EPwaves.chanIX{iColor}=zeros(length(eLabels),1);
                for i=1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).chanNames)
                    EPwaves.chanIX{iColor}(i)=find(strcmp(eLabels{i},EPwaves.chanNames)); %which of the master chan list each current channel corresponds to
                end;
                %[B,EPwaves.chanIX{iColor}] = sort(eLabels);
                oldEloc=eloc(EPwaves.chanIX{iColor});
                newEloc=EPdataset.dataset(EPmain.view.dataset(iColor)).eloc;
%                 numChans=length(newEloc);
                if (length([oldEloc.theta]) ~= length([newEloc.theta])) || (length([oldEloc.radius]) ~= length([newEloc.radius]))
                    msg{1}='Error: The electrode coordinates are not compatible.';
                    [msg]=ep_errorMsg(msg);
                    close(EPwaves.handles.waves.hWaveWindow);
                    err=1;
                    return
                end;
                if any([newEloc.theta]-[oldEloc.theta]) || any([oldEloc.radius]-[newEloc.radius])
                    msg{1}='Error: The electrode coordinates are not compatible.';
                    [msg]=ep_errorMsg(msg);
                    close(EPwaves.handles.waves.hWaveWindow);
                    err=1;
                    return
                end;
            elseif xor(isempty(eloc),isempty(newEloc)) %if one but not the other has no electrode coordinates
                msg{1}='Error: The electrode coordinates are not compatible.';
                [msg]=ep_errorMsg(msg);
                close(EPwaves.handles.waves.hWaveWindow);
                err=1;
                return
            elseif length(eloc) ~= length(newEloc)
                msg{1}='Error: The electrode coordinates are not compatible.';
                [msg]=ep_errorMsg(msg);
                close(EPwaves.handles.waves.hWaveWindow);
                err=1;
                return
            else
                %if both have no electrode coordinates but have the same number of channels, assume they are compatible
                EPwaves.chanIX{iColor}=[1:numChans];
            end;
        end;
    end;
end;

if startTime==endTime
    EPwaves.firstTime=max(EPwaves.firstTime,startTime);
    EPwaves.lastTime=min(EPwaves.lastTime,endTime); %assume already  left side of sample
else
    EPwaves.firstTime=max(EPwaves.firstTime,startTime);
    EPwaves.lastTime=min(EPwaves.lastTime,endTime-EPwaves.sampleSize); %shift endTime to left side of sample
end;
EPwaves.firstHz=max(EPwaves.firstHz,startHz);
EPwaves.lastHz=min(EPwaves.lastHz,endHz);

if strcmp('FFT',EPmain.view.dataTransform)
    if ~isempty(EPwaves.marker1)
        if (EPwaves.marker1 > EPwaves.lastHz) || (EPwaves.marker1 < EPwaves.firstHz)
            EPwaves.marker1=0;
        end;
    end;
    if ~isempty(EPwaves.marker2)
        if (EPwaves.marker2 > EPwaves.lastHz) || (EPwaves.marker2 < EPwaves.firstHz)
            EPwaves.marker2=0;
        end;
    end;
else
    if ~isempty(EPwaves.marker1)
        if (EPwaves.marker1 > EPwaves.lastTime) || (EPwaves.marker1 < EPwaves.firstTime)
            EPwaves.marker1=0;
        end;
    end;
    if ~isempty(EPwaves.marker2)
        if (EPwaves.marker2 > EPwaves.lastTime) || (EPwaves.marker2 < EPwaves.firstTime)
            EPwaves.marker2=0;
        end;
    end;
end;

%compute epoch samples for each color
EPwaves.numPoints=1;
EPwaves.numHz=1;
EPwaves.startBins=ones(1,4);
EPwaves.lastBins=ones(1,4);
EPwaves.startSamp=ones(1,4);
EPwaves.lastSamp=ones(1,4);

if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
    for iColor=1:4
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                [M EPwaves.startBins(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames-EPwaves.firstHz));
                [M EPwaves.lastBins(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames-EPwaves.lastHz));
                EPwaves.lastBins(iColor)=max(EPwaves.lastBins(iColor),EPwaves.startBins(iColor));
            EPwaves.numHz=EPwaves.lastBins(iColor)-EPwaves.startBins(iColor)+1;
        end;
    end
end;
if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
    for iColor=1:4
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                [M EPwaves.startSamp(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames-EPwaves.firstTime));
                [M EPwaves.lastSamp(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames-EPwaves.lastTime));
                EPwaves.lastSamp(iColor)=max(EPwaves.lastSamp(iColor),EPwaves.startSamp(iColor));
            EPwaves.numPoints=EPwaves.lastSamp(iColor)-EPwaves.startSamp(iColor)+1;
        end;
    end
end;
    
%organize the data for plotting
EPwaves.plotColors=[];
EPwaves.thePlotColors=[];
EPwaves.totalData=zeros(numChans,EPwaves.numPoints,0,EPwaves.numHz); %the 4 dimensions of totalData are chans, points, waves(trials/cells/subjects), and frequencies
EPwaves.bandData=zeros(numChans,EPwaves.numPoints,0,EPwaves.numHz,4); %the 5 dimensions of bandData are chans, points, waves(trials/cells/subjects), frequencies, and band colors
EPwaves.colorIndex=[];
EPwaves.bandIndex=zeros(4,1);
EPwaves.eventWave=cell(4,1);
EPwaves.boundary=cell(4,1);
EPwaves.plotColorIndex=[];
EPwaves.plotLineIndex=cell(0);
EPwaves.correlPlotScaleIndex=[];
EPwaves.STScells=[];
theDataset=0;
theMax=0;
EPwaves.complexData=0;

for iColor=1:4
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if theDataset ~= EPmain.view.dataset(iColor) %load in the dataset if different from the one already loaded in.
            EPdata=ep_loadEPdataset(EPdataset,EPmain.view.dataset(iColor));
            theDataset=EPmain.view.dataset(iColor);
        end;
        theSub=EPmain.view.subject(iColor);
        if strcmp(EPdata.dataType,'average')
            theCell=EPmain.view.cell(iColor);
            if ismember(EPmain.view.allTrials(iColor),[5 6])
                theCell=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames)];
            end;
            if ismember(EPmain.view.allTrials(iColor),[1 2])
                theSub=find(strcmp('AVG',EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes));
                if isempty(theSub)
                    theSub=find(strcmp('GAV',EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes));
                end;
            end;
        else  %if single_trial data or continuous
            cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
            if ismember(EPmain.view.allTrials(iColor),[1 2 5 6]) || strcmp(EPdata.dataType,'continuous')
                theCell=cellList;
            else
                theCell=cellList(EPmain.view.trial(iColor));
            end;
        end;
        theFactor=EPmain.view.factor(iColor);
        if ismember(EPmain.view.allTrials(iColor),[3 4])
            theFactor=find(strcmp('SGL',EPdataset.dataset(EPmain.view.dataset(iColor)).facTypes));
            if isempty(theFactor)
                theFactor=find(strcmp('CMB',EPdataset.dataset(EPmain.view.dataset(iColor)).facTypes));
            end; 
        end;
        tempData=ep_expandFacs(EPdata,[],EPwaves.startSamp(iColor):EPwaves.lastSamp(iColor),theCell,theSub,theFactor,EPwaves.startBins(iColor):EPwaves.lastBins(iColor));
        
        if strcmp(EPdata.dataType,'average') && (length(theCell)==1)
            EPwaves.bandIndex(iColor)=0;
            switch EPmain.view.trial(iColor)
                case 1 %nothing
                case 2 %GFP
                    GFP=std(tempData(EEGchans,:,:,:,:,:,:),1);
                    tempData(EEGchans,:,:,:,:,:,:)=repmat(GFP,length(EEGchans),1);
                case 3 %noise
                    EPwaves.bandData(:,:,1,:,iColor)=EPdata.noise(:,:,theCell,theSub,theFactor);
                    EPwaves.bandIndex(iColor)=1;
                case 4 %StDev
                    if ~isempty(EPdata.std)
                        EPwaves.bandData(:,:,1,:,iColor)=EPdata.std(:,:,theCell,theSub,theFactor);
                        EPwaves.bandIndex(iColor)=2;
                    else
                        disp('No std information');
                    end;
                case 5 %CI
                    if ~isempty(EPdata.std)
                        if EPdata.subNum(theSub,theCell) > 1 %if grand average
                            numCIobs=EPdata.subNum(theSub,theCell);
                            theStd=EPdata.stdCM(:,:,theCell,theSub,theFactor);
                        else
                            numCIobs=EPdata.avgNum(theSub,theCell);
                            theStd=EPdata.std(:,:,theCell,theSub,theFactor);
                        end;
                        if ft_hastoolbox('STATS', 0, 1)
                            ts=tinv(.92,numCIobs);
                        else
                            if ~noStatsWarn
                                disp('Since the Matlab Statistics toolbox is not installed, will use a value appropriate for n=30 as a rough estimate for the standard error.')
                                ts=1.44;
                                %1.44 provides 84% CI for sample size of 30, which is about right to provide an alpha of .05 for overlapping error bar test when standard errors are about equal.
                            end;
                        end;
                        EPwaves.bandData(:,:,1,:,iColor)=ts*theStd/sqrt(numCIobs);
                        EPwaves.bandIndex(iColor)=3;
                    else
                        disp('No std information');
                    end;
            end;
        end;
        
        tempEvents=[];
        if strcmp(EPdata.dataType,'continuous')
            if strcmp('VLT',EPmain.view.dataTransform)
                for theChan=1:size(tempData,1)
                    tempData(theChan,:,:,:,:)=tempData(theChan,:,:,:,:)-mean(tempData(theChan,:,:,:,:),2); %center the waveforms if continuous
                end;
            end;
            if ~isempty(EPdata.events{EPmain.view.subject(iColor),1})
                tempEvents=EPdata.events{EPmain.view.subject(iColor),1}(([EPdata.events{1}.sample]>=EPwaves.startSamp(iColor)) & ([EPdata.events{1}.sample]<=EPwaves.lastSamp(iColor)));
            end;
        else
            for iCell=1:length(theCell)
                for iSub=1:length(theSub)
                    tempEvents=[tempEvents EPdata.events{theSub(iSub),theCell(iCell)}];
                end;
            end;
        end;
        if ~ismember(EPmain.view.allTrials(iColor),[2 4 6]) %if not erpimage
            EPwaves.eventWave{iColor}=cell(1);
            EPwaves.eventWave{iColor}{1}=zeros(1,size(tempData,2));
            EPwaves.boundary{iColor}=[];
            if ~isempty(EPwaves.eventLines{iColor})
                EPwaves.eventWave{iColor}{1}=histc(EPwaves.eventLines{iColor},[1:size(tempData,2)]);
                theMax=max([theMax EPwaves.eventWave{iColor}{1}]);
            end;
            if ~isempty(tempEvents)
                boundaryEvents=find(strcmp('boundary',{tempEvents.value}));
                if ~isempty(boundaryEvents)
                    EPwaves.boundary{iColor}=tempEvents(boundaryEvents).sample;
                end;
            end;
            if ~any(any(EPwaves.eventWave{iColor}{1}))
                EPwaves.eventWave{iColor}{1}=[];
            end;
        end;
        numWaves=max([size(tempData,3) size(tempData,4) size(tempData,5)]); %if any of them have been set to -all- or -erpimage-
        if ~isreal(tempData) && (FFTunits==1) %if the data has an imaginary component, as in spectral data
            tempDataImag=imag(tempData);
            tempData=real(tempData);
        else
            tempDataImag=[];
        end;
        if ~isempty(EPdata.relNames) %collapse over the relations dimension if any to provide a mean
            if strcmp(EPdata.dataType,'average')
                goodRelChans=find(squeeze(any(any(~isnan(EPdata.analysis.badChans(theSub,theCell,:)),2),1)));
            else
                goodRelChans=find(squeeze(any(any((EPdata.analysis.badChans(theSub,theCell,:)~=-1),2),1)));
            end;
            refChans=EPdata.reference.current;
            if isempty(refChans) && ~any(ismember({'AVG','CSD'},EPdata.reference.type))
                refChans=EPdata.reference.original;
            end;
            if length(refChans)==1
                goodRelChans=setdiff(goodRelChans,refChans); %coherence with a single reference channel is NaN.
            end;
            tempData=mean(abs(tempData(:,:,:,:,:,:,goodRelChans)),7);
            if ~isempty(tempDataImag)
                tempDataImag=mean(abs(tempDataImag(:,:,:,:,:,:,goodRelChans)),7); 
            end;
        end;
        EPwaves.totalData(EPwaves.chanIX{iColor},:,end+1:end+numWaves,:)=squeeze(tempData); %rearrange order of channels to be consistent with other datasets
        EPwaves.colorIndex(end+1:end+numWaves)=iColor; %which of the four data colors does the waveform belong to
        EPwaves.plotColorIndex=[EPwaves.plotColorIndex; repmat(EPwaves.RGBcolors(iColor,:),numWaves,1)]; %the plotting colors of the lines for each waveform
        EPwaves.plotColors=[EPwaves.plotColors iColor]; %the data colors that will actually be plotted
        EPwaves.thePlotColors=[EPwaves.thePlotColors; EPwaves.RGBcolors(iColor,:)]; %the plotting colors of the lines for each of these data colors
        [EPwaves.plotLineIndex{end+1:end+numWaves}]=deal('-'); %the plotting colors of the lines for each waveform
        EPwaves.correlPlotScaleIndex=[EPwaves.correlPlotScaleIndex; repmat(EPmain.view.correl(iColor),numWaves,1)]; %is each waveform a correlation?
        if length(theCell)==1
            EPwaves.STScells=[EPwaves.STScells; repmat(strcmp('STS',{EPdataset.dataset(EPmain.view.dataset(iColor)).cellTypes{theCell}}),numWaves,1)]; %is each waveform an STS output?
        else
            EPwaves.STScells=[EPwaves.STScells; strcmp('STS',{EPdataset.dataset(EPmain.view.dataset(iColor)).cellTypes{theCell}})']; %is each waveform an STS output?
        end;
        if ~isempty(tempDataImag)
            EPwaves.totalData(EPwaves.chanIX{iColor},:,end+1:end+numWaves,:)=squeeze(tempDataImag); %rearrange order of channels to be consistent with other datasets
            EPwaves.colorIndex(end+1:end+numWaves)=iColor; %which of the four data colors does the waveform belong to
            EPwaves.plotColorIndex=[EPwaves.plotColorIndex; repmat(EPwaves.RGBcolors(iColor,:),numWaves,1)]; %the plotting colors of the lines for each waveform
            [EPwaves.plotLineIndex{end+1:end+numWaves}]=deal(':'); %the plotting colors of the lines for each waveform
            EPwaves.correlPlotScaleIndex=[EPwaves.correlPlotScaleIndex; repmat(EPmain.view.correl(iColor),numWaves,1)]; %is each waveform a correlation?
            EPwaves.complexData=1;
        end;
    end;
end;


for iColor=1:4
    if ~isempty(EPwaves.eventLines{iColor}) && ~ismember(EPmain.view.allTrials(iColor),[2 4 6]) %if not erpimage
        if ~isempty(EPwaves.eventLines{iColor})
            EPwaves.eventWave{iColor}{1}=EPwaves.eventWave{iColor}{1}/theMax; %rescale event information to between zero and one.
        end;
    end;
end;

%check for sample test output
plotDataList=find(EPmain.view.dataset <= length(EPdataset.dataset));
EPwaves.STSdata=[];
for iColor=1:length(plotDataList)
    theColor=plotDataList(iColor);
    if ismember(EPmain.view.allTrials(iColor),[1 2 5 6])
        theCell=1;
    else
        theCell=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(theColor)}(EPmain.view.cell(theColor)),EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames));
        theCell=theCell(1);
    end;
    if strcmp('STS',EPdataset.dataset(EPmain.view.dataset(theColor)).cellTypes(theCell))
        EPwaves.STSdata=[EPwaves.STSdata; theColor];
    end;
end;
if (length(EPwaves.STSdata)~=1)  || (length(plotDataList)~=3) || (length(EPwaves.plotColors)~=3)
    EPwaves.STmode=0;
else
    EPwaves.STmode=1; %present STS results as red zone between the two waves.
end;

%scale waveforms
if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
    nonCorrWaves=setdiff(find(~EPwaves.correlPlotScaleIndex),EPwaves.STSdata);
    if (FFTunits > 1)
        EPwaves.totalData(EEGchans,:,nonCorrWaves,:)=abs(EPwaves.totalData(EEGchans,:,nonCorrWaves,:)); %convert complex number to real number
    end;
    EPwaves.totalData(EEGchans,:,nonCorrWaves,:)=EPwaves.totalData(EEGchans,:,nonCorrWaves,:)/sqrt(mean(diff(EPdata.freqNames))); %convert to spectral density
    if FFTunits > 2
        EPwaves.totalData(EEGchans,:,nonCorrWaves,:)=EPwaves.totalData(EEGchans,:,nonCorrWaves,:).^2; %convert amplitude to power
    end;
    if (FFTunits == 4)
        if ~all(EPwaves.totalData(EEGchans,:,nonCorrWaves,:) >=0)
            disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
        end;
        EPwaves.totalData(EEGchans,:,nonCorrWaves,:)=log10(abs(EPwaves.totalData(EEGchans,:,nonCorrWaves,:)))*10; %convert to dB log scaling
        tempVar=EPwaves.totalData(EEGchans,:,nonCorrWaves,:);
        tempVar(isinf(tempVar))=-flintmax;
        EPwaves.totalData(EEGchans,:,nonCorrWaves,:)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
    end;
end;

theMax=max(max(max(max(EPwaves.totalData(EEGchans,:,intersect(find(~EPwaves.correlPlotScaleIndex),find(~EPwaves.STScells)),:)))));

if ~isempty(plotMVmin)
    EPwaves.plotMVmin=plotMVmin;
else
    EPwaves.plotMVmin=-theMax;
end;
if ~isempty(plotMVmax)
    EPwaves.plotMVmax=plotMVmax;
else
    EPwaves.plotMVmax=theMax;
end;

if ~all(EPwaves.correlPlotScaleIndex) && any(EPwaves.correlPlotScaleIndex) %if not all the waveforms are correlations but some are, then rescale the correlations to match the plots
    EPwaves.totalData(EEGchans,:,find(EPwaves.correlPlotScaleIndex),:)=(EPwaves.totalData(EEGchans,:,find(EPwaves.correlPlotScaleIndex),:)*(EPwaves.plotMVmax-EPwaves.plotMVmin))+EPwaves.plotMVmin;
end;

%compute screen locations of channels
EPwaves.X=zeros(numChans,1);
EPwaves.Y=zeros(numChans,1);
hasLoc=zeros(numChans,1);

if ~isempty(eloc) %if there are electrode coordinates
    nonLoc=0;
    for chan=1:numChans
        theta=((eloc(chan).theta+90)/180)*pi;
        radius=eloc(chan).radius;
        if isempty(theta)
            nonLoc=nonLoc+1;
            EPwaves.X(chan)=.02;
            EPwaves.Y(chan)=.02+((EPwaves.plotHeight+.01)*nonLoc);
        else
            [EPwaves.X(chan),EPwaves.Y(chan)] = pol2cart(theta,radius);
            hasLoc(chan)=1;
        end;
    end;
    
    for chan=1:length(implicit)
        theta=((implicit(chan).theta+90)/180)*pi;
        radius=implicit(chan).radius;
        if isempty(theta)
            EPwaves.X(end+1)=0;
            EPwaves.Y(end+1)=0;
            hasLoc(end+1)=0;
        else
            [EPwaves.X(end+1),EPwaves.Y(end+1)] = pol2cart(theta,radius);
            hasLoc(end+1)=1;
        end;
    end;
else %if there are no electrode coordinates
    graphSize=ceil(sqrt(numChans));
    for chan=1:numChans
        EPwaves.X(chan)=1-(mod(chan-1,graphSize)+1)*(EPwaves.plotWidth+.01);
        EPwaves.Y(chan)=1-(floor((chan-1)/graphSize)+1)*(EPwaves.plotHeight+.01);
        hasLoc(chan)=1;
    end;
    for chan=1:length(implicit)
        if strcmp(implicit(chan).type,'REF')
            EPwaves.X(end+1)=1-(mod(length(EPwaves.X),graphSize)+1)*(EPwaves.plotWidth+.01);
            EPwaves.Y(end+1)=1-(floor(length(EPwaves.Y)/graphSize)+1)*(EPwaves.plotHeight+.01);
            hasLoc(end+1)=1;
        end;
    end;
end;

locList=find(hasLoc);
EPwaves.X(locList)=(1-EPwaves.X(locList))/2;
EPwaves.Y(locList)=(EPwaves.Y(locList)+1)/2;

if any(diff(EPwaves.X)) %don't adjust if all along the midline
    EPwaves.X(locList)=((EPwaves.X(locList)-min(EPwaves.X(locList)))/(max(EPwaves.X(locList)-min(EPwaves.X(locList)))*(1+2*EPwaves.margin)))+EPwaves.margin;
end;
if any(diff(EPwaves.Y)) %don't adjust if all along the center line
    EPwaves.Y(locList)=((EPwaves.Y(locList)-min(EPwaves.Y(locList)))/(max(EPwaves.Y(locList)-min(EPwaves.Y(locList)))*(1+2*EPwaves.margin)))+EPwaves.margin;
end;

set(EPwaves.handles.waves.hWaveWindow ,'DefaultAxesXTickLabel',[])
set(EPwaves.handles.waves.hWaveWindow ,'DefaultAxesYTickLabel',[])
set(EPwaves.handles.waves.hWaveWindow ,'DefaultAxesXTick',[])
set(EPwaves.handles.waves.hWaveWindow ,'DefaultAxesYTick',[])
nonERPimages=find(~ismember(EPwaves.colorIndex,find(ismember(EPmain.view.allTrials,[2 4 6]))));
if ~isempty(nonERPimages)
    set(EPwaves.handles.waves.hWaveWindow,'DefaultAxesColorOrder',EPwaves.plotColorIndex(nonERPimages,:)); %do not plot erpimages
end;

%plot waveforms

EPwaves.handles.waves.hWave=zeros(numChans,1);

EPwaves.plotForm=EPmain.view.dataTransform;
if strcmp('TFT',EPmain.view.dataTransform) && (EPwaves.firstHz==EPwaves.lastHz)
    EPwaves.plotForm='VLT';
end;
if strcmp('TFT',EPmain.view.dataTransform) && (EPwaves.firstTime==EPwaves.lastTime)
    EPwaves.plotForm='FFT';
end;

if strcmp('FFT',EPwaves.plotForm)
    sampleSize=0;
    EPwaves.spacing=(EPwaves.lastHz-EPwaves.firstHz)/(EPwaves.numHz-1);
else
    sampleSize=EPwaves.sampleSize;
    EPwaves.spacing=EPwaves.sampleSize;
end;
EPwaves.handles.waves.hLines=cell(numChans,1);

waveList=find(~ismember(EPwaves.colorIndex,find(ismember(EPmain.view.allTrials,[2 4 6]))));
theMarker='none';
theMarkerSize=2;
switch EPwaves.plotForm
    case 'VLT'
        for iChan=1:numChans
            numImages=length(find(ismember(EPmain.view.allTrials,[2 4 6])));
            if (length(EPwaves.plotColors)-numImages) > 0
                imageSpace=4;
            else
                imageSpace=numImages;
            end;
            tempHandles=[];
            if numImages %if any erpimage
                imageCount=0;
                for iColor=1:4
                    if (EPmain.view.dataset(iColor) <= length(EPdataset.dataset)) && ismember(EPmain.view.allTrials(iColor),[2 4 6])
                        imageCount=imageCount+1;
                        trialList=find(EPwaves.colorIndex==iColor);
                        tempHandles(1,end+1) = axes('position',[EPwaves.X(iChan) EPwaves.Y(iChan)+(EPwaves.plotHeight/imageSpace)*(imageSpace-imageCount) EPwaves.plotWidth (EPwaves.plotHeight/imageSpace)],'XTickMode','manual','YTickMode','manual');
                        tempHandles(1,end+1) = imagesc(EPwaves.firstTime:EPwaves.lastTime,1:length(trialList),squeeze(EPwaves.totalData(iChan,:,trialList,:))',[EPwaves.plotMVmin, EPwaves.plotMVmax]);
                        axis([EPwaves.firstTime EPwaves.lastTime 1 length(trialList)]);
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])
                        line([0 0],[1 length(trialList)],'Color','black','LineWidth',1) %stimulus onset
                        if ~isempty(EPwaves.marker1)
                            line(repmat(EPwaves.marker1,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                        end
                        if ~isempty(EPwaves.marker2)
                            line(repmat(EPwaves.marker2,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                        end
                        %plot event lines
                        eventX=[];
                        eventY=[];
                        for iLine=1:length(EPwaves.eventLines{iColor})
                            if ~isempty(EPwaves.eventLines{iColor}{iLine})
                                eventX(end+1)=EPwaves.eventLines{iColor}{iLine};
                                eventY(end+1)=iLine;
                            end;
                        end;
                        if ~isempty(eventX)
                            theTimePoints=[EPwaves.firstTime:EPwaves.spacing:EPwaves.firstTime+(EPwaves.spacing*(EPwaves.numPoints-1))];
                            hold on
                            plot(theTimePoints(eventX),eventY)
                            hold off
                        end
                        set(tempHandles(end-1),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                        set(tempHandles(end),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                    end
                end
            end
            if (length(EPwaves.plotColors)-numImages) > 0 %if there will be waveforms
                EPwaves.handles.waves.hWave(iChan) = axes('position',[EPwaves.X(iChan) EPwaves.Y(iChan) EPwaves.plotWidth EPwaves.plotHeight*((4-numImages)/4)],'XTickMode','manual','YTickMode','manual');
                hold on
                for iWave=1:length(waveList)
                    theWave=waveList(iWave);
                    theColor=EPwaves.colorIndex(theWave);
                    if EPwaves.complexData
                        if strcmp(EPwaves.plotLineIndex{theWave},':')
                            theMarker='none';
                            theMarkerSize=2;
                        else
                            theMarker='none';
                            theMarkerSize=2;
                        end
                    end
                    if EPwaves.STmode && EPwaves.STSdata==EPwaves.colorIndex(iWave)
                        breakList=sort([find(diff([0 (squeeze(EPwaves.totalData(iChan,:,theWave,:))>0) 0])<0)-1 find(diff([0 (squeeze(EPwaves.totalData(iChan,:,theWave,:))>0) 0])>0)]);
                        if ~isempty(breakList)
                            theSTdata=squeeze(EPwaves.totalData(iChan,:,theWave,:));
                            theData1=squeeze(EPwaves.totalData(iChan,:,min(setdiff(plotDataList,EPwaves.STSdata)),:));
                            theData2=squeeze(EPwaves.totalData(iChan,:,max(setdiff(plotDataList,EPwaves.STSdata)),:));
                            for iSigArea=1:length(breakList)/2
                                theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)];
                                if length(theTimePoints) == 1
                                    EPwaves.handles.waves.hLines{iChan}(iWave)=plot(([theTimePoints theTimePoints]*EPwaves.spacing)+(EPwaves.firstTime-EPwaves.spacing),[theData1(theTimePoints) theData2(theTimePoints)],'LineWidth',1,'Color',[1 .5 .5]);
                                else
                                    EPwaves.handles.waves.hLines{iChan}(iWave)=patch(([theTimePoints flip(theTimePoints)]*EPwaves.spacing)+(EPwaves.firstTime-EPwaves.spacing),[theData1(theTimePoints) theData2(flip(theTimePoints))],EPwaves.plotColorIndex(iWave,:),'FaceColor',EPwaves.plotColorIndex(iWave,:),'EdgeColor','none','FaceAlpha',.25);
                                end;
                            end;
                        end;
                    else
                        theTimePoints=[EPwaves.firstTime:EPwaves.spacing:EPwaves.firstTime+(EPwaves.spacing*(EPwaves.numPoints-1))];
                        EPwaves.handles.waves.hLines{iChan}(iWave)=plot(theTimePoints,squeeze(EPwaves.totalData(iChan,:,theWave,:)),'LineStyle',EPwaves.plotLineIndex{theWave},'color',EPwaves.plotColorIndex(theWave,:),'Marker',theMarker,'MarkerSize',theMarkerSize);
                        if EPwaves.bandIndex(theColor) > 0
                            theBand1=squeeze(EPwaves.totalData(iChan,:,theWave,:))+EPwaves.bandData(iChan,:,1,:,theColor);
                            theBand2=squeeze(EPwaves.totalData(iChan,:,theWave,:))-EPwaves.bandData(iChan,:,1,:,theColor);
                            EPwaves.handles.waves.hLines{iChan}(iWave+length(waveList))=patch([theTimePoints flip(theTimePoints)],[theBand1 flip(theBand2)],EPwaves.plotColorIndex(iWave,:),'FaceColor',EPwaves.plotColorIndex(iWave,:),'EdgeColor',EPwaves.plotColorIndex(iWave,:),'FaceAlpha',.25);
                        end;
                    end;
                end;
                hold off
                axis([EPwaves.firstTime EPwaves.lastTime EPwaves.plotMVmin EPwaves.plotMVmax]);
                if EPwaves.direction ==2
                    set(EPwaves.handles.waves.hWave(iChan),'YDir','reverse')
                end;
                EPwaves.handles.waves.zero(iChan)=line([EPwaves.firstTime EPwaves.lastTime-EPwaves.sampleSize],[0 0],'Color','black','LineWidth',1); % zero line
                EPwaves.handles.waves.onset(iChan)=line([0 0],[0 EPwaves.plotMVmax],'Color','black','LineWidth',1); %stimulus onset
                set(EPwaves.handles.waves.zero(iChan),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                set(EPwaves.handles.waves.onset(iChan),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                
                %plot event lines
                for iColor=1:length(EPwaves.plotColors)
                    theColor=EPwaves.plotColors(iColor);
                    if ~isempty(EPwaves.eventWave{theColor}{1})
                        plotPoints=find(EPwaves.eventWave{theColor}{1}>min(EPwaves.eventWave{theColor}{1}));
                        plotTimes=[EPwaves.firstTime:EPwaves.spacing:EPwaves.firstTime+(EPwaves.spacing*(EPwaves.numPoints-1))];
                        if length(plotPoints)==1
                            EPwaves.handles.waves.eventLines{iChan,theColor} = line([plotTimes(plotPoints) plotTimes(plotPoints)],[EPwaves.plotMVmin EPwaves.eventWave{theColor}{1}(plotPoints)*(EPwaves.plotMVmin/2)],'Color',EPwaves.thePlotColors(theColor,:),'LineWidth',2); %event line
                        else
                            hold on
                            EPwaves.handles.waves.eventLines{iChan,theColor} = plot([EPwaves.firstTime:EPwaves.spacing:EPwaves.firstTime+(EPwaves.spacing*(EPwaves.numPoints-1))],(EPwaves.eventWave{theColor}{1}*(abs(EPwaves.plotMVmin/2)))+EPwaves.plotMVmin,'LineWidth',5,'Color',EPwaves.thePlotColors(theColor,:));
                            hold off
                            for iLine=1:length(EPwaves.handles.waves.eventLines{iChan,theColor})
                                set(EPwaves.handles.waves.eventLines{iChan,theColor}(iLine),'YDataSource',['EPwaves.eventWave{' num2str(theColor) '}{1}(' num2str(iLine) ')']);
                            end;
                        end;
                    end;
                    if ~isempty(EPwaves.boundary{theColor})
                        hold on
                        for iBoundary=1:length(EPwaves.boundary{theColor})
                            theSample=EPwaves.boundary{theColor}(iBoundary);
                            EPwaves.handles.waves.boundary{iChan,theColor} = line([theSample theSample],[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',EPwaves.thePlotColors(theColor,:),'LineWidth',1);
                        end;
                        hold off
                    end;
                end;
            end;
            EPwaves.handles.waves.hLines{iChan}=[EPwaves.handles.waves.hLines{iChan}, tempHandles];
            text(.05,.8, EPwaves.chanNames(iChan), 'Units','normalized');
            if ~isempty(EPwaves.marker1)
                line(repmat(EPwaves.marker1,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
            end;
            if ~isempty(EPwaves.marker2)
                line(repmat(EPwaves.marker2,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
            end;
            if (length(EPwaves.plotColors)-numImages) > 0
                set(EPwaves.handles.waves.hWave(iChan),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
            end;
            for iLine=1:length(EPwaves.handles.waves.hLines{iChan})
                if isgraphics(EPwaves.handles.waves.hLines{iChan}(iLine))
                    set(EPwaves.handles.waves.hLines{iChan}(iLine),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                end;
            end;
            for iLine=1:length(EPwaves.handles.waves.hLines{iChan})-(numImages*2)
                if isprop(EPwaves.handles.waves.hLines{iChan}(iLine),'YDataSource')
                    set(EPwaves.handles.waves.hLines{iChan}(iLine),'YDataSource',['squeeze(EPwaves.totalData(' num2str(iChan) ',:,' num2str(iLine) ',:,:))']);
                end;
            end;
        end;
        
    case 'FFT'
        for iChan=1:numChans
            numImages=length(find(ismember(EPmain.view.allTrials,[2 4 6])));
            if (length(EPwaves.plotColors)-numImages) > 0
                imageSpace=4;
            else
                imageSpace=numImages;
            end;
            tempHandles=[];
            if numImages %if any erpimage
                imageCount=0;
                for i=1:4
                    if (EPmain.view.dataset(i) <= length(EPdataset.dataset)) && ismember(EPmain.view.allTrials(i),[2 4 6])
                        imageCount=imageCount+1;
                        trialList=find(EPwaves.colorIndex==i);
                        tempHandles(1,end+1) = axes('position',[EPwaves.X(iChan) EPwaves.Y(iChan)+(EPwaves.plotHeight/imageSpace)*(imageSpace-imageCount) EPwaves.plotWidth (EPwaves.plotHeight/imageSpace)],'XTickMode','manual','YTickMode','manual');
                        tempHandles(1,end+1) = imagesc(EPwaves.firstHz+EPwaves.sampleSize:EPwaves.lastHz,1:length(trialList),squeeze(EPwaves.totalData(iChan,:,trialList,:))',[EPwaves.plotMVmin, EPwaves.plotMVmax]);
                        axis([EPwaves.firstHz+EPwaves.sampleSize EPwaves.lastHz 1 length(trialList)]);
                        if ~isempty(EPwaves.marker1)
                            line(repmat(EPwaves.marker1,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                        end;
                        if ~isempty(EPwaves.marker2)
                            line(repmat(EPwaves.marker2,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',1);
                        end;
                        set(tempHandles(end-1),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                        set(tempHandles(end),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                    end;
                end;
            end
            if (length(EPwaves.plotColors)-numImages) > 0 %if there will be waveforms
                EPwaves.handles.waves.hWave(iChan) = axes('position',[EPwaves.X(iChan) EPwaves.Y(iChan) EPwaves.plotWidth EPwaves.plotHeight*((4-numImages)/4)],'XTickMode','manual','YTickMode','manual');
                hold on
                for iWave=1:length(waveList)
                    theWave=waveList(iWave);
                    if EPwaves.complexData
                        if strcmp(EPwaves.plotLineIndex{theWave},':')
                            theMarker='none';
                            theMarkerSize=2;
                        else
                            theMarker='*';
                            theMarkerSize=2;
                        end;
                    end;
                    EPwaves.handles.waves.hLines{iChan}(iWave)=plot([EPwaves.firstHz+sampleSize:EPwaves.spacing:EPwaves.lastHz],squeeze(EPwaves.totalData(iChan,:,theWave,:)),'LineStyle',EPwaves.plotLineIndex{theWave},'color',EPwaves.plotColorIndex(theWave,:),'Marker',theMarker,'MarkerSize',theMarkerSize);
                end;
                hold off
                axis([EPwaves.firstHz+sampleSize EPwaves.lastHz EPwaves.plotMVmin EPwaves.plotMVmax]);
                set(EPwaves.handles.waves.hWave(iChan),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                EPwaves.handles.waves.zero(iChan)=line([EPwaves.firstHz+EPwaves.sampleSize EPwaves.lastHz],[0 0],'Color','black','LineWidth',1); % zero line
                set(EPwaves.handles.waves.zero(iChan),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
            end;
            EPwaves.handles.waves.hLines{iChan}=[EPwaves.handles.waves.hLines{iChan}, tempHandles];
            text(.05,.2, EPwaves.chanNames(iChan), 'Units','normalized');
            if ~isempty(EPwaves.marker1)
                line(repmat(EPwaves.marker1,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
            end;
            if ~isempty(EPwaves.marker2)
                line(repmat(EPwaves.marker2,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
            end;
            for i=1:length(EPwaves.handles.waves.hLines{iChan})
                set(EPwaves.handles.waves.hLines{iChan}(i),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
            end;
            for i=1:length(EPwaves.handles.waves.hLines{iChan})-(numImages*2)
                if ~EPwaves.STmode || EPwaves.STSdata~=EPwaves.plotColors(i)
                    set(EPwaves.handles.waves.hLines{iChan}(i),'YDataSource',['squeeze(EPwaves.totalData(' num2str(iChan) ',:,' num2str(i) ',:,:))']);
                end;
            end;
        end;
    case 'TFT'
        imageSpace=length(EPwaves.plotColors);
        for iChan=1:numChans
            imageCount=0;
            for i=1:4
                if EPmain.view.dataset(i) <= length(EPdataset.dataset)
                    imageCount=imageCount+1;
                    EPwaves.handles.waves.hLines{iChan}(i) = axes('position',[EPwaves.X(iChan) EPwaves.Y(iChan)+(EPwaves.plotHeight/imageSpace)*(imageSpace-imageCount) EPwaves.plotWidth (EPwaves.plotHeight/imageSpace)]);
                    EPwaves.handles.waves.hLines{iChan}(4+i) = imagesc(EPwaves.firstTime:EPwaves.lastTime,EPwaves.firstHz:EPwaves.lastHz,squeeze(EPwaves.totalData(iChan,:,(EPwaves.plotColors==i),:))');
                    axis([EPwaves.firstTime EPwaves.lastTime EPwaves.firstHz EPwaves.lastHz]);
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    line([0 0],[EPwaves.firstHz EPwaves.lastHz],'Color','black','LineWidth',1) %stimulus onset
                    if ~isempty(EPwaves.marker1)
                        line(repmat(EPwaves.marker1,2),[EPwaves.firstHz EPwaves.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                    if ~isempty(EPwaves.marker2)
                        line(repmat(EPwaves.marker2,2),[EPwaves.firstHz EPwaves.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
                    end;
                    set(EPwaves.handles.waves.hLines{iChan}(i),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                    set(EPwaves.handles.waves.hLines{iChan}(4+i),'ButtonDownFcn',{@ep_expandChan,'EPwaves'});
                end;
            end;
            text(.05,-.1, EPwaves.chanNames(iChan), 'Units','normalized');
        end;
end;

if any(strcmp(EPwaves.plotForm,{'FFT','VLT'})) && showLegend
    %set up the legend
    numCols=length(find(EPmain.view.dataset <= length(EPdataset.dataset)));
    colData=cell(0);
    colCell=cell(0);
    colSub=cell(0);
    colTrial=cell(0);
    colFactor=cell(0);
    
    EPwaves.legendNames=cell(0);
    EPwaves.legendNames{1}='';
    for iColor=1:4
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            
            colData{end+1}=EPdataset.dataset(EPmain.view.dataset(iColor)).dataName;
            
            if EPmain.view.allTrials(iColor)==5
                colCell{end+1}='all';
            elseif EPmain.view.allTrials(iColor)==6
                colCell{end+1}='erpimage';
            else
                colCell(end+1)=EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor));
            end;
            
            if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average')
                if EPmain.view.allTrials(iColor)==1
                    subName='all';
                elseif EPmain.view.allTrials(iColor)==2
                    subName='erpimage';
                else
                    subName=EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)};
                end;
            else
                subName=EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)};
            end;
            colSub{end+1}=subName;
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                if ismember(EPmain.view.allTrials(iColor),[1 2])
                    colTrial{end+1}='all trials';
                else
                    colTrial{end+1}=num2str(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(EPmain.view.trial(iColor)));
                end;
            elseif strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average')
                switch EPmain.view.trial(iColor)
                    case 1 %nothing
                        colTrial{end+1}='data';
                    case 2 %GFP
                        colTrial{end+1}='GFP';
                    case 3 %noise
                        colTrial{end+1}='noise';
                    case 4 %StDev
                        colTrial{end+1}='StDev';
                    case 5 %95CI
                        colTrial{end+1}='95%CI';
                end;
            end;
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                if EPmain.view.allTrials(iColor)==3
                    facName='all';
                elseif EPmain.view.allTrials(iColor)==4
                    facName='erpimage';
                else
                    facName=EPdataset.dataset(EPmain.view.dataset(iColor)).facNames{EPmain.view.factor(iColor)};
                end;
                colFactor{end+1}=facName;
            end;
            
        end;
    end;
    if length(unique(colData))==numCols
        EPwaves.legendNames=colData;
    elseif length(unique(colCell))==numCols
        EPwaves.legendNames=colCell;
    elseif length(unique(colSub))==numCols
        EPwaves.legendNames=colSub;
    elseif (length(unique(colTrial))==numCols) && ~isempty(colTrial)
        EPwaves.legendNames=colTrial;
    elseif (length(unique(colFactor))==numCols) && ~isempty(colFactor)
        EPwaves.legendNames=colFactor;
    else
        nameCounter=1;
        EPwaves.legendNames=cell(1,numCols);
        for iColor=1:4
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                if length(unique(colData))~=1
                    EPwaves.legendNames{nameCounter}=[colData{nameCounter}];
                end;
                if length(unique(colCell))~=1
                    if ~isempty(EPwaves.legendNames{nameCounter})
                        EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} '-'];
                    end;
                    EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} colCell{nameCounter}];
                end;
                if length(unique(colSub))~=1
                    if ~isempty(EPwaves.legendNames{nameCounter})
                        EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} '-'];
                    end;
                    EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} colSub{nameCounter}];
                end;
                if ~isempty(colTrial) && length(unique(colTrial))~=1
                    if ~isempty(EPwaves.legendNames{nameCounter})
                        EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} '-'];
                    end;
                    EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} colTrial{nameCounter}];
                end;
                if ~isempty(colFactor) && length(unique(colFactor))~=1
                    if ~isempty(EPwaves.legendNames{nameCounter})
                        EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} '-'];
                    end;
                    EPwaves.legendNames{nameCounter}=[EPwaves.legendNames{nameCounter} colFactor{nameCounter}];
                end;
                if isempty(EPwaves.legendNames{nameCounter})
                    EPwaves.legendNames{nameCounter}=sprintf('Dataset %d',iColor);
                end;
                nameCounter=nameCounter+1;
            end;
        end;
    end;
    set(EPwaves.handles.waves.hWaveWindow,'DefaultAxesColorOrder',EPwaves.RGBcolors(find(EPmain.view.dataset <= length(EPdataset.dataset)),:));
    legendAxes = axes('position',[1-EPwaves.plotWidth 1-EPwaves.plotHeight EPwaves.plotWidth EPwaves.plotHeight],'XTickMode','manual','YTickMode','manual');
    hold on
    plot(legendAxes,[1:2],ones(2,length(find(EPmain.view.dataset <= length(EPdataset.dataset)))));
    hleg1 = legend('string',EPwaves.legendNames,'Location',[1-EPwaves.plotWidth 1-EPwaves.plotHeight EPwaves.plotWidth EPwaves.plotHeight],'Interpreter','none');
end;

theLabel='';
theLabel2='';

if all(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
    theLabel='Correlation';
elseif strcmp(EPmain.view.dataTransform,'VLT')
    if strcmp(EPdata.reference.type,'CSD')
        theLabel='V/m^2';
    else
        theLabel='µv';
    end;
elseif any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
    if (FFTunits > 2)
        if (FFTunits == 4)
        	theLabel='dBV';
        else
            if strcmp(EPdata.reference.type,'CSD')
                theLabel='(V^2)/(Hz*m^4)';
            else
                theLabel='(µv^2)/Hz';
            end;
        end;
    else
        if strcmp(EPdata.reference.type,'CSD')
            theLabel='V/(sqrt(Hz)*m^2)';
        else
            theLabel='µv/sqrt(Hz)';
        end;
    end;
end;

%scale graph
switch EPwaves.plotForm
    case 'VLT'
        EPwaves.handles.waves.graph = axes('position',[.02 .02 EPwaves.plotWidth EPwaves.plotHeight],'XTickMode','manual','YTickMode','manual');
        axis([EPwaves.firstTime EPwaves.lastTime EPwaves.plotMVmin EPwaves.plotMVmax]);
        if EPwaves.direction ==2
            set(EPwaves.handles.waves.graph,'YDir','reverse')
        end;
        line([EPwaves.firstTime EPwaves.lastTime],[0 0],'Color','black','LineWidth',1) % zero line
        if strcmp('VLT',EPmain.view.dataTransform)
            line([0 0],[0 EPwaves.plotMVmax],'Color','black','LineWidth',1) %stimulus onset
        end;
        text(EPwaves.firstTime+(EPwaves.lastTime-EPwaves.firstTime)/20, EPwaves.plotMVmin-(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+d',round(EPwaves.firstTime)));
        text(EPwaves.lastTime-(EPwaves.lastTime-EPwaves.firstTime)/5, EPwaves.plotMVmin-(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+d',round(EPwaves.lastTime)));
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.plotMVmax-(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+2.2f',EPwaves.plotMVmax));
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.plotMVmin+(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+2.2f',EPwaves.plotMVmin));
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.plotMVmax-(EPwaves.plotMVmax-EPwaves.plotMVmin)/4, theLabel);
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.plotMVmax-(EPwaves.plotMVmax-EPwaves.plotMVmin)/2, theLabel2);        
        
        if ~isempty(EPwaves.marker1)
            line(repmat(EPwaves.marker1,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
        end;
        if ~isempty(EPwaves.marker2)
            line(repmat(EPwaves.marker2,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
        end;
    case 'FFT'
        EPwaves.handles.waves.graph = axes('position',[.02 .02 EPwaves.plotWidth EPwaves.plotHeight],'XTickMode','manual','YTickMode','manual');
        axis([EPwaves.firstHz+sampleSize EPwaves.lastHz EPwaves.plotMVmin EPwaves.plotMVmax]);
        if EPwaves.direction ==2
            set(EPwaves.handles.waves.graph,'YDir','reverse')
        end;
        line([EPwaves.firstHz+sampleSize EPwaves.lastHz],[0 0],'Color','black','LineWidth',1) % zero line
        if strcmp('VLT',EPmain.view.dataTransform)
            line([0 0],[0 EPwaves.plotMVmax],'Color','black','LineWidth',1) %stimulus onset
        end;
        text(EPwaves.firstHz+(EPwaves.lastHz-EPwaves.firstHz)/20, EPwaves.plotMVmin-(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+d',round(EPwaves.firstHz)));
        text(EPwaves.lastHz-(EPwaves.lastHz-EPwaves.firstHz)/5, EPwaves.plotMVmin-(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+d',round(EPwaves.lastHz)));
        text(EPwaves.firstHz-(EPwaves.lastHz-EPwaves.firstHz)/3, EPwaves.plotMVmax-(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+2.2f',EPwaves.plotMVmax));
        text(EPwaves.firstHz-(EPwaves.lastHz-EPwaves.firstHz)/3, EPwaves.plotMVmin+(EPwaves.plotMVmax-EPwaves.plotMVmin)/8, sprintf('%+2.2f',EPwaves.plotMVmin));
        text(EPwaves.firstHz-(EPwaves.lastHz-EPwaves.firstHz)/3, EPwaves.plotMVmax+(EPwaves.plotMVmax-EPwaves.plotMVmin)/2, theLabel);
        text(EPwaves.firstHz-(EPwaves.lastHz-EPwaves.firstHz)/3, EPwaves.plotMVmax+(EPwaves.plotMVmax-EPwaves.plotMVmin)/4, theLabel2);

        if ~isempty(EPwaves.marker1)
            line(repmat(EPwaves.marker1,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
        end;
        if ~isempty(EPwaves.marker2)
            line(repmat(EPwaves.marker2,2),[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',[.5 .5 .5],'LineWidth',1);
        end;
    case 'TFT'
        EPwaves.handles.waves.graph = axes('position',[.02 .02 EPwaves.plotWidth EPwaves.plotHeight],'XTickMode','manual','YTickMode','manual');
        axis([EPwaves.firstTime EPwaves.lastTime EPwaves.firstHz EPwaves.lastHz]);
        line([EPwaves.firstTime EPwaves.lastTime],[0 0],'Color','black','LineWidth',1) % zero line
        if strcmp('VLT',EPmain.view.dataTransform)
            line([0 0],[0 EPwaves.lastHz],'Color','black','LineWidth',1) %stimulus onset
        end;
        text(EPwaves.firstTime+(EPwaves.lastTime-EPwaves.firstTime)/20, EPwaves.firstHz-(EPwaves.lastHz-EPwaves.firstHz)/8, sprintf('%+d',round(EPwaves.firstTime)));
        text(EPwaves.lastTime-(EPwaves.lastTime-EPwaves.firstTime)/5, EPwaves.firstHz-(EPwaves.lastHz-EPwaves.firstHz)/8, sprintf('%+d',round(EPwaves.lastTime)));
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.lastHz-(EPwaves.lastHz-EPwaves.firstHz)/8, sprintf('%2.2f',EPwaves.firstHz));
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.firstHz+(EPwaves.lastHz-EPwaves.firstHz)/8, sprintf('%2.2f',EPwaves.lastHz));
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.lastHz+(EPwaves.lastHz-EPwaves.firstHz)/2, theLabel);
        text(EPwaves.firstTime-(EPwaves.lastTime-EPwaves.firstTime)/3, EPwaves.lastHz+(EPwaves.lastHz-EPwaves.firstHz)/4, theLabel2);

        if ~isempty(EPwaves.marker1)
            line(repmat(EPwaves.marker1,2),[EPwaves.firstHz EPwaves.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
        end;
        if ~isempty(EPwaves.marker2)
            line(repmat(EPwaves.marker2,2),[EPwaves.firstHz EPwaves.lastHz],'Color',[.5 .5 .5],'LineWidth',1);
        end;
end;

set(EPwaves.handles.waves.hWaveWindow ,'Color',[.5 .5 .5])


