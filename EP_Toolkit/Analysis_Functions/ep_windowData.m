function outputdata=ep_windowData(EPdata, chanGrp, inputCells, outCellNames, subjectSpecs, startwin, endwin, measure, fname, factor, startfreq, endfreq, FFTunits, chanColl, missingData, peakWidth, numberType);

%  outputdata=ep_windowData(EPdata, chanGrp, inputCells, outCellNames, subjectSpecs, startwin, endwin, measure, fname, factor, startfreq, endfreq, FFTunits, missingData, peakWidth,numberType);
%
%  Outputs text file suitable for ANOVA.  Each line is a subject.  As for columns, cells will vary slowest while chan groups will
%  vary fastest.  Chan groups will be generated for spatial as well as temporal PCAs even though chan effects will not occur for spatial PCAs.
%  Results will be in microvolts (or msec) for a specific channel (or mean of group of channels) for a specific window to make results exactly comparable
%  to conventional analyses.  Strips EPdata of adds (combined channels, cells, subjects, and factors) first.
%  Starts with seven header lines: name of dataset, window, type of measure, channel group, factor if any, cell names, and channel area
%  names.
%  For single-trial data, each line will be a trial, the columns will be cells, and no specs will be included.
%
%Inputs:
%   EPdata         : Structured array with the data and accompanying information.  See readData.
%
% chanGrp   : Cell array of channels going into each output channel grouping.
%   .name       : Name of the channel grouping
%   .threshold  : Threshold used for factor groups (not used)
%   .areaname   : Cell array of the channel area names
%   .activeFac  : Active factor (not used)
%   .channel    : Vector of area each channel is assigned to
%
% inputCells : Cell array of conditions included in analysis (output cells, each containing array of input cells)
% outCellNames: Cell array of the output cell names
% subjectSpecs : Logical vector array of which specs to include with the output measures
% startwin   : Start of time window in samples (1-based).
% endwin	 : End of time window in samples.
% measure	 : Type of measure for windows: mean, minpeak, maxpeak, minlatency, maxlatency, mincentroid, maxcentroid, behavioral.
% fname		 : name of output file
% factor     : Factor if a PCA dataset, otherwise set to 1.
% startfreq   : Start of frequency window in bins (1-based).
% endfreq	 : End of frequency window in bins.
%  FFTunits  : spectral data units (1=complex, 3=asd, 3=psd, 4=dB)
% chanColl:   : Collapse channels then meausure (1) or measure then collapse (2)
% missingData : Number used to denote missing data.
% peakWidth   : Number of samples on either side of peak for peak measures (e.g., "2" equals 5 total samples).
% numberType  : Type of data to output: 'imag', 'real', or 'normal' for EEG data or 'RT' or 'ACC' for behavioral data.
%
%Outputs:
%   outputdata  : The resulting windowed measures.  Also written out to disk.
%
%  Centroid measure of latency:
%  Dien, J., Spencer, K. M., & Donchin, E. (2004). Parsing the "Late Positive Complex": Mental chronometry and the ERP
%  components that inhabit the neighborhood of the P300. Psychophysiology, 41(5), 665-678.
%  The present equation is corrected from that in the published article.
%
%  General introductory treatment of ANOVA of ERP data:
%  Dien, J., & Santuzzi, A. M. (2005). Application of repeated measures ANOVA to high-density ERP datasets: A review and
%  tutorial. In T. Handy (Ed.), Event-Related Potentials: A Methods Handbook. Cambridge, Mass: MIT Press.

%History:
%  by Joseph Dien (7/17/09)
%  jdien07@mac.com
%
% bugfix 10/10/10 JD
% Fixed centroid measures too large due to addition of the left side of the window (e.g., +200 ms for a window of
% 200-300 ms).
%
% modified 4/13/12 JD
% Added support for 6th dimension of frequency, including scaling to psd and dB.
%
% modified 4/23/12 JD
% Changed default for channel group windowing to collapsing channels first
% then taking measures.  Added preference option so can use original
% approach too of measuring first.
%
% modified 5/24/12 JD
% Added support for missing data when windowing data.
%
% bugfix 8/4/12 JD
% Fixed bad channel handling referring only to first channel in channel area.  Fixed missing data numbers being
% transformed to psd and dB for frequency data.
%
% modified 1/16/13 JD
% Handles situation where FFT data has negative values (e.g., due to imprecision in PCA results) and transforming to dB
% would result in complex numbers, by taking absolute value first.
%
% modified 1/29/13 JD
% Added frequency-domain peak measures.
%
% bugfix 2/4/13 JD
% Changed ms window information provided in header so it ranges from onset of sample to offset of sample rather than
% onset to onset (e.g., 0-4 ms rather than 0-0 ms for first sample).
%
% bugfix 3/27/13 JD
% Fixed windowed files were labeled as being in dB even when voltage data.
%
% modified 5/7/13 JD
% Added option to average together samples around a peak to minpeak and maxpeak measures.
% Implemented Luck (2005) suggestion to only count as a dip/peak a sample where both neighboring samples are higher/lower.
% When there is missing data (more likely now with new peak measure code), the NaNs are converted to missing data code
% specified in preferences setting.
%
% bugfix 5/23/13 JD
% Fixed crash when windowing using the minpeak or maxpeak measures with the adjoining samples option and the peak latency was at the upper end of the window.
% Fixed crash when there are multiple channels in the area of a channel group.
% Fixed minpeak and maxpeak measures yielding missing data numbers when the window size was less than three samples, as
% in the autoPCA mode.
%
% modified 3/18/14 JD
% Added option to window single-trial data.
%
% modified 4/24/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure, including support for complex numbers.
%
% bugfix 8/30/15 JD
% Fixed amplitude spectral density calculated as divided by Hz rather than
% by square root of Hz.
% Fixed dB of amplitude data not converted to power first.
%
% bugfix 11/11/15 JD
% Fixed crash when windowing with the "measure then collapse" option.
%
% bugfix 1/15/16 JD
% Now allows power scaled data to be output as amplitudes.
% Now handles complex FFT numbers.
%
% modified 1/21/16 JD
% Consolidated spectral unit controls so just four options (cmp, asd, psd, dB).
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% bugfix 6/18/17 JD
% Fixed crash when windowing frequency-domain data.
% Fixed Window function text file headers not mentioning dB scaling for spectral data.
% Switch to amplitude scaling when adding freq data together other than channels.
%
% modified 4/8/18 JD
% Added support for outputing behavioral data.
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

outputdata=[];

if all(chanGrp.channel == length(chanGrp.areaName))
    msg{1}='Error: No channels specified.';
    [msg]=ep_errorMsg(msg);
    return
end;

[EPdata]=ep_stripAdds(EPdata);
if isempty(EPdata.data)
    msg{1}='Error: The file had no data left after additions were removed.';
    [msg]=ep_errorMsg(msg);
    return
end;

switch numberType
    case {'normal','RT','ACC'}
        %do nothing
    case 'real'
        EPdata.data=real(EPdata.data);
    case 'imag'
        EPdata.data=imag(EPdata.data);
    otherwise
        disp('oops, programming mistake');
        return
end;

numPoints=length(EPdata.timeNames);
numChans=length(EPdata.chanNames);
numCells=length(unique(EPdata.cellNames));
numWaves=length(EPdata.cellNames);
numSubs=length(EPdata.subNames);
numFacs=length(EPdata.facNames);
numFreqs=length(EPdata.freqNames);
cellNames=unique(EPdata.cellNames);
numRels=length(EPdata.relNames);

if (endwin < startwin) && (numPoints > 0)
    msg{1}='Endwindow sample is smaller than startwindow sample.';
    [msg]=ep_errorMsg(msg);
    return
end;

if (endfreq < startfreq) && (numFreqs > 0)
    msg{1}='Endfrequency bin is smaller than startfrequency bin.';
    [msg]=ep_errorMsg(msg);
    return
end;

if (endwin > numPoints) && (numPoints > 0)
    msg{1}='Endwindow sample is larger than number of time points.';
    [msg]=ep_errorMsg(msg);
    return
end;

if factor == 0
    factor =1;
end;

if ~numPoints
    startwin=1;
    endwin=1;
end;

if ~numFreqs
    startfreq=1;
    endfreq=1;
end;

if strcmp(measure,'behavioral')
    numAreas=1;
    chanArea=1;
    chanAreaName{1}=numberType;
else
    %determine number of channel areas and which channels are assigned to them.
    numAreas=0;
    for theArea=1:length(chanGrp.areaName)-1
        chanArray=find(chanGrp.channel == theArea);
        if ~isempty(chanArray)
            numAreas=numAreas+1;
            chanArea{numAreas}=chanArray;
            chanAreaName{numAreas}=chanGrp.areaName{theArea};
        end;
    end;
end;

%construct array of windowed measures
if strcmp(EPdata.dataType,'single_trial')
    outputdata = zeros(sum(cellfun(@length,inputCells)),numAreas);
    trialCount=0;
    trialLabels=cell(size(outputdata,1),2);
else
    outputdata = zeros(numSubs,(numAreas*length(inputCells)));
end;
multChansFlag=0;
for theSub = 1:numSubs
    for theOutCell = 1:length(inputCells)
        for theArea = 1:numAreas
            theData = 0;
            theTrialData=ones(length(inputCells{theOutCell}),1)*missingData; %if single_trial data
            count = 0;
            theOkayCells=inputCells{theOutCell}(find(EPdata.avgNum(theSub,inputCells{theOutCell})>=0));
            if ~isempty(theOkayCells)
                for theInCell = 1:length(theOkayCells)
                    theCell=theOkayCells(theInCell);
                    if strcmp(measure,'behavioral')
                        if strcmp(numberType,'RT')
                            theData=EPdata.trialSpecs{theCell,find(strcmp('RT',EPdata.trialSpecNames)),theSub};
                        else
                            theData=EPdata.trialSpecs{theCell,find(strcmp('ACC',EPdata.trialSpecNames)),theSub};
                        end;
                        count=1;
                    else
                        if strcmp(EPdata.dataType,'average')
                            goodChans=find(~isnan(EPdata.analysis.badChans(theSub,theCell,:)));
                        else
                            goodChans=find(EPdata.analysis.badChans(theSub,theCell,:) >= 0);
                        end;
                        chanList=intersect(chanArea{theArea},goodChans);
                        if ~isempty(chanList)
                            chanLength=1;
                            numChansInArea=length(chanList);
                            if numChansInArea > 1
                                multChansFlag=1;
                            end;
                            if chanColl ==2 %collapse channels together after measuring
                                chanLength=numChansInArea;
                            end;
                            if numRels
                                if multChansFlag
                                    relChans=chanList; %take mean across the coherences between the channels in the area group, if any
                                else
                                    relChans=[]; %take mean across the coherences of the single channel with all the other channels.
                                end;
                            else
                                relChans=[]; %no coherences
                            end;
                            
                            if strcmp(EPdata.dataType,'average')
                                goodRelChans=find(squeeze(~isnan(EPdata.analysis.badChans(theSub,theCell,:))));
                            else
                                goodRelChans=find(squeeze((EPdata.analysis.badChans(theSub,theCell,:)~=-1)));
                            end;
                            refChans=EPdata.reference.current;
                            if isempty(refChans) && ~any(ismember({'AVG','CSD'},EPdata.reference.type))
                                refChans=EPdata.reference.original;
                            end;
                            if length(refChans)==1
                                goodRelChans=setdiff(goodRelChans,refChans); %coherence with a single reference channel is NaN.
                            end;
                            
                            for theChan = 1:chanLength
                                if chanColl ==1 %collapse channels together prior to measuring
                                    epochData=zeros(1,length(startwin:endwin),1,1,1,length(startfreq:endfreq));
                                    for theChanColl=chanList(:)'
                                        chanEpoch=ep_expandFacs(EPdata,theChanColl,startwin:endwin,theCell,theSub,factor,startfreq:endfreq,relChans);
                                        if numRels
                                            if multChansFlag
                                                chanEpoch=mean(chanEpoch(:,:,:,:,:,:,goodRelChans),7);
                                            else
                                                chanEpoch=mean(abs(chanEpoch(:,:,:,:,:,:,goodRelChans)),7); %for coherence, measure absolute value when taking the mean of all channels (since they would more or less sum to zero).
                                            end;
                                        end;
                                        epochData=epochData+chanEpoch;
                                    end;
                                    epochData=epochData./numChansInArea;
                                else
                                    chanEpoch=ep_expandFacs(EPdata,chanList(theChan),startwin:endwin,theCell,theSub,factor,startfreq:endfreq,relChans);
                                    if numRels
                                        if multChansFlag
                                            chanEpoch=mean(chanEpoch(:,:,:,:,:,:,goodRelChans),7);
                                        else
                                            chanEpoch=mean(abs(chanEpoch(:,:,:,:,:,:,goodRelChans)),7); %for coherence, measure absolute value when taking the mean of all channels (since they would more or less sum to zero).
                                        end;
                                    end;
                                    epochData=chanEpoch;
                                end;
                                if ~isempty(EPdata.freqNames) && any(strcmp(EPdata.chanTypes{chanList},{'EEG','REG'}))
                                    epochData=abs(epochData); %convert to amplitude scaling when adding freq data together.
                                end;
                                epoch=squeeze(mean(epochData,6)); %take mean across the frequency window, if any
                                epoch=epoch(:)';
                                HzEpoch=epochData;
                                HzEpoch=squeeze(mean(HzEpoch,2)); %take mean across the time window, if any
                                HzEpoch=HzEpoch(:)';
                                switch measure
                                    case 'mean', theData = theData + mean(epoch,2);
                                    case 'minpeak',
                                        if length(epoch) > 2
                                            dips=find(diff(sign(diff(epoch)))==2)+1;
                                            dips=dips(find(dips > peakWidth));
                                            dips=dips(find(dips <= length(epoch)-peakWidth));
                                        else
                                            dips=[1:length(epoch)];
                                        end;
                                        [dummy, minpeak]=min(epoch(dips));
                                        peakLatency=dips(minpeak);
                                        peakMean=mean(epoch(peakLatency-peakWidth:peakLatency+peakWidth),2);
                                        theData = theData + peakMean;
                                    case 'maxpeak',
                                        if length(epoch) > 2
                                            peaks=find(diff(sign(diff(epoch)))==-2)+1;
                                            peaks=peaks(find(peaks > peakWidth));
                                            peaks=peaks(find(peaks <= length(epoch)-peakWidth));
                                        else
                                            dips=[1:length(epoch)];
                                        end;
                                        [dummy, maxpeak]=max(epoch(peaks));
                                        peakLatency=peaks(maxpeak);
                                        peakMean=mean(epoch(peakLatency-peakWidth:peakLatency+peakWidth),2);
                                        theData = theData + peakMean;
                                    case 'minlatency',
                                        dips=find(diff(sign(diff(epoch)))==2)+1;
                                        [dummy, minpeak]=min(epoch(dips));
                                        peakLatency=dips(minpeak);
                                        theData = theData + (peakLatency+startwin-1-EPdata.baseline)*(1000/EPdata.Fs);
                                    case 'maxlatency',
                                        peaks=find(diff(sign(diff(epoch)))==-2)+1;
                                        [dummy, maxpeak]=max(epoch(peaks));
                                        peakLatency=peaks(maxpeak);
                                        theData = theData + (peakLatency+startwin-1-EPdata.baseline)*(1000/EPdata.Fs);
                                    case 'mincentroid',
                                        maxvalue=max(epoch,[],2);
                                        theData = theData + (mean(([startwin:endwin].*(epoch-maxvalue))/mean(epoch-maxvalue))- EPdata.baseline)*(1000/EPdata.Fs);
                                    case 'maxcentroid',
                                        minvalue=min(epoch,[],2);
                                        theData = theData + (mean(([startwin:endwin].*(epoch-minvalue))/mean(epoch-minvalue))- EPdata.baseline)*(1000/EPdata.Fs);
                                    case 'minHzPeak', theData = theData + EPdata.freqNames(round(min(HzEpoch,[],2)));
                                    case 'maxHzPeak', theData = theData + EPdata.freqNames(round(max(HzEpoch,[],2)));
                                    case 'minHzLatency', [dummy, latency] = min(HzEpoch,[],2);
                                        theData = theData + EPdata.freqNames(latency);
                                    case 'maxHzLatency', [dummy, latency] = max(HzEpoch,[],2);
                                        theData = theData + EPdata.freqNames(latency);
                                    case 'minHzCentroid',
                                        maxvalue=max(HzEpoch,[],2);
                                        theData = theData + EPdata.freqNames(round(mean(([startfreq:endfreq].*(HzEpoch-maxvalue))/mean(HzEpoch-maxvalue))));
                                    case 'maxHzCentroid',
                                        minvalue=min(HzEpoch,[],2);
                                        theData = theData + EPdata.freqNames(round(mean(([startfreq:endfreq].*(HzEpoch-minvalue))/mean(HzEpoch-minvalue))));
                                    otherwise,
                                        msg{1}='Measure should be one of the following options: mean, minpeak, maxpeak, minlatency, maxlatency, mincentroid, maxcentroid, minHzPeak, maxHzPeak, minHzLatency, maxHzLatency, minHzCentroid, maxHzCentroid';
                                        [msg]=ep_errorMsg(msg);
                                        return
                                end;
                                count = count +1;
                            end;
                        end;
                    end;
                    if strcmp(EPdata.dataType,'single_trial')
                        theTrialData(theInCell)=theData/count;
                        theData=0;
                        count=0;
                    end;
                end;
            end;
            if (count == 0 && ~strcmp(EPdata.dataType,'single_trial'))
                outputdata(theSub,(theOutCell-1)*numAreas+theArea) = missingData;
            else
                if strcmp(EPdata.dataType,'single_trial')
                    theData=theTrialData;
                else
                    theData = theData /count;
                end;
                if any(strcmp(measure,{'mean','minpeak','maxpeak'})) && numFreqs && ~numRels
                    if length(EPdata.freqNames) > 1
                        if (FFTunits > 1)
                            theData=abs(theData); %convert complex number to real number
                        end;
                        theData=theData/sqrt(mean(diff(EPdata.freqNames))); %convert to spectral density
                        if FFTunits > 2
                            theData=theData.^2; %convert amplitude to power
                        end;
                        if (FFTunits == 4)
                            if ~all(theData >=0)
                                disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
                            end;
                            theData=log10(abs(theData))*10; %convert to dB log scaling
                            theData(isinf(theData))=-flintmax; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                        end;
                    end;
                end;
                if strcmp(EPdata.dataType,'single_trial')
                    if isempty(theData)
                        outputdata(trialCount+1:trialCount+length(inputCells{theOutCell}),theArea) = missingData;
                    else
                        outputdata(trialCount+1:trialCount+length(inputCells{theOutCell}),theArea) = theData;
                    end;
                    [trialLabels{trialCount+1:trialCount+length(inputCells{theOutCell}),1}]=deal(outCellNames{theOutCell});
                    [trialLabels(trialCount+1:trialCount+length(inputCells{theOutCell}),2)]=deal(num2cell([1:length(inputCells{theOutCell})]));
                else
                    if isempty(theData)
                        outputdata(theSub,(theOutCell-1)*numAreas+theArea) = missingData;
                    else
                        outputdata(theSub,(theOutCell-1)*numAreas+theArea) = theData;
                    end;
                end;
            end;
        end;
        if strcmp(EPdata.dataType,'single_trial')
            trialCount=trialCount+length(inputCells{theOutCell});
        end;
    end;
end;

outputdata(isnan(outputdata))=missingData; %change NaNs to missing data code number


if any(strcmp(measure,{'mean','minpeak','maxpeak'})) && numFreqs && ~strcmp(measure,'behavioral')
    if length(EPdata.freqNames) == 1
        disp('Cannot determine frequency bin size for spectral density conversion.');
    end;
end;

switch numberType
    case 'normal'
        %do nothing
    case 'real'
        fname=[fname(1:end-4) '_real.txt'];
    case 'imag'
        fname=[fname(1:end-4) '_imag.txt'];
    case 'RT'
        fname=[fname(1:end-4) '_RT.txt'];
    case 'ACC'
        fname=[fname(1:end-4) '_ACC.txt'];
    otherwise
        disp('oops, programming mistake');
        return
end;

theName=fname;
sameName=1;
theNumber=0;
[pathstr, name, ext] = fileparts(theName);
fileName=name;
while sameName
    sameName=0;
    if exist([pathstr filesep fileName ext],'file')
        sameName=1;
    end;
    if sameName
        theNumber=theNumber+1;
        fileName=[name '-' num2str(theNumber)];
    end;
end;

outFID=fopen([pathstr filesep fileName ext],'w');
if (outFID == -1)
    msg{1}='Error creating output file!';
    [msg]=ep_errorMsg(msg);
    return
end;

%file header
if strcmp(measure,'behavioral')
    fprintf(outFID,'%s\r','behavioral');
    fprintf(outFID,'\r');
    fprintf(outFID,'\r');
    fprintf(outFID,'\r');
    fprintf(outFID,'\r');
else
    fprintf(outFID,'%s\r',EPdata.dataName);
    
    if numPoints
        fprintf(outFID,'%d-%d ms',EPdata.timeNames(startwin),EPdata.timeNames(endwin)+(1000/EPdata.Fs));
        if numFreqs
            fprintf(outFID,' ');
        end;
    end;
    if numFreqs
        fprintf(outFID,'%5.1f-%5.1f Hz',EPdata.freqNames(startfreq),EPdata.freqNames(endfreq));
    end;
    fprintf(outFID,'\r');
    if isempty(EPdata.freqNames)
        dataType='voltage';
    else
        if isempty(EPdata.timeNames)
            dataType='FFT';
        else
            dataType='TFT';
        end;
        if ~isempty(EPdata.relNames)
            switch numberType
                case 'normal'
                    dataType=['coherence ' dataType];
                case 'real'
                    dataType=['real coherence ' dataType];
                case 'imag'
                    dataType=['imaginary coherence ' dataType];
            end;
        end;
        switch FFTunits
            case 1
                dataType=['complex amplitude spectral density ' dataType];
            case 2
                dataType=['amplitude spectral density ' dataType];
            case 3
                dataType=['internally calculated as amplitude with final units as power spectral density ' dataType];
            case 4
                dataType=['internally calculated as amplitude with final units as power spectral density ' dataType];
        end;
    end;
    
    if (FFTunits == 4) && ~isempty(EPdata.freqNames)
        dataType=[dataType ' in dB'];
    end;
    
    if strcmp(EPdata.reference.type,'CSD')
        dataType=['CSD of ' dataType];
    end;
    
    if multChansFlag
        if chanColl ==1
            dataType=[dataType '. Channels were collapsed together prior to taking the measure.'];
        elseif chanColl ==2
            dataType=[dataType '. Channels were measured and then the measures were collapsed together.'];
        else
            disp('Programming error.');
        end;
    else
        dataType=[dataType '.'];
    end;
    
    if (startwin == endwin) && (startfreq == endfreq) && strcmp(measure,'mean')
        if isempty(EPdata.freqNames)
            fprintf(outFID,'%s\r',['One sample of ' dataType]);
        else
            if isempty(EPdata.timeNames)
                fprintf(outFID,'%s\r',['One frequency bin of ' dataType]);
            else
                fprintf(outFID,'%s\r',['One frequency wavelet of ' dataType]);
            end;
        end;
    else
        fprintf(outFID,'%s\r',[measure ' of ' dataType]);
    end;
    
    fprintf(outFID,'%s\r',chanGrp.name);
    if isempty(EPdata.facNames)
        fprintf(outFID,'\r');
    else
        fprintf(outFID,'Factor: %d\r',factor);
    end;
end;

if strcmp(EPdata.dataType,'single_trial')
    fprintf(outFID,'\r');
    for theArea = 1:length(chanArea)
        fprintf(outFID,'%s\t',chanAreaName{theArea});
    end
    fprintf(outFID,'\tcell\ttrial\r');
    
    %file data
    for iTrial=1:size(outputdata,1)
        for iChan = 1:size(outputdata,2)
            fprintf(outFID,'%f\t',outputdata(iTrial,iChan));
        end;
        fprintf(outFID,'\t%s', trialLabels{iTrial,1});
        fprintf(outFID,'\t%4.0f\r',trialLabels{iTrial,2});
    end;
else
    for theCell = 1:length(outCellNames)
        for theArea = 1:length(chanArea)
            fprintf(outFID,'%s\t',outCellNames{theCell});
        end
    end;
    
    for theSpec = 1:length(subjectSpecs)
        fprintf(outFID,'%s',EPdata.subjectSpecNames{subjectSpecs(theSpec)});
        if theSpec < length(subjectSpecs)
            fprintf(outFID,'\t');
        end;
    end;
    fprintf(outFID,'\r');
    for theCell = 1:length(outCellNames)
        for theArea = 1:length(chanArea)
            fprintf(outFID,'%s\t',chanAreaName{theArea});
        end
    end;
    
    for theSpec = 1:length(subjectSpecs)
        fprintf(outFID,'spec\t');
    end;
    fprintf(outFID,'\r');
    
    %file data
    for theSub=1:numSubs
        for theCell = 1:size(outputdata,2)
            fprintf(outFID,'%f\t',outputdata(theSub,theCell));
        end;
        for theSpec = 1:length(subjectSpecs)
            theSpec=EPdata.subjectSpecs{theSub,subjectSpecs(theSpec)};
            if isnumeric(theSpec)
                theSpec=num2str(theSpec);
            end;
            fprintf(outFID,'%s\t',theSpec);
        end;
        fprintf(outFID,[EPdata.subNames{theSub} '\r']);
    end;
end;



fclose(outFID);
