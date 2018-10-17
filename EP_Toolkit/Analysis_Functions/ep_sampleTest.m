function [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdataIn,method,alpha,contiguous,scale,templateWaveform,templateTopo,test,thresh,channel);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdata,method,alpha,contiguous,scale,templateWaveform,templateTopo,test,thresh,channel);% Performs sample-by-sample, CWT, or Woody Filter (non-iterative, covariance) using templates tests.% Uses either non-parametric test or jack-knife onset measures.% The resulting significance waveform is then added to the EPdata format for subsequent inspection.% It is assumed that the input data will already have only two cells selected to be present.% They will be dubbed cell1 and cell2 by alphabetical order.% Additionally, continuous data can be processed, in which case no significance testing is performed but a raw Woody filter waveform is added to the output dataset.%%Inputs%   EPdataIn         : Structured array with the data and accompanying information.  See readData.%   method: type of sample test ('sample','CWT','Template Woody','PCA Woody').  Both Woody options are treated the same.%   alpha: Alpha threshold for determining significance%   contiguous: Minimum number of contiguous significant samples to be deemed especially significant (a crude MCP control method)%   scale: The width in ms of the Mexican hat wavelet used for the t-CWT method%   templateWaveform: The microvolt-scaled waveform for the Woody options.%   templateTopo: The microvolt-scaled topographies for the Woody options.%   test: Non-parametric 't-test' or 'jack-knife'%   thresh: threshold for onset measure (jack-knife only).  If threshold is negative then will look for below it.%   channel: channel for onset measure if not full-head template (jack-knife or continuous data)%%Outputs%   outputData         : Array of sample-by-sample results (chans,samples)%                        0=non-significant and .25=significant and .5=significant and meets contiguity criterion and%                        .75=the most significant regardless of contiguity.%   sampLat1           : ms latencies for sample 1 if using Woody option with templateTopo%   sampLat2           : ms latencies for sample 2 if using Woody option with templateTopo%   sampAmp1           : max amplitude for sample 1 if using Woody option with templateTopo%   sampAmp2           : max amplitude for sample 2 if using Woody option with templateTopo%%  Bostanov, V., & Kotchoubey, B. (2006). The t-CWT: a new ERP detection and quantification method based on %  the continuous wavelet transform and Student?s t-statistics. Clin Neurophysiol, 117(12), 2627-2644.%%  Woody, C. D. (1967). Characterization of an adaptive filter for the analysis %  of variable latency neuroelectric signals. Medical and Biological Engineering, 5, 539-553.%%  Miller, J., Patterson, T., & Ulrich, R. (1998). Jackknife-based method for measuring LRP onset latency differences. %  Psychophysiology, 35(1), 99-115.%% History:%% by Joseph Dien (7/13/14)% jdien07@mac.com%% bugfix 3/21/16 JD% Fixed PCA Woody option of SampleTest function crashing with spatial PCA templates.%% modified 11/3/16 JD% Added support for template Woody and for continuous data.%% modified 6/20/17 JD% Added support for using sampleTest function with TFT data.% Added support for contrast between channels.% Added support for contrasting in the frequency domain rather than the time domain.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Copyright (C) 1999-2018  Joseph Dien%%     This program is free software: you can redistribute it and/or modify%     it under the terms of the GNU General Public License as published by%     the Free Software Foundation, either version 3 of the License, or%     (at your option) any later version.%%     This program is distributed in the hope that it will be useful,%     but WITHOUT ANY WARRANTY; without even the implied warranty of%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the%     GNU General Public License for more details.%%     You should have received a copy of the GNU General Public License%     along with this program.  If not, see <http://www.gnu.org/licenses/>.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%outputData=[];sampLat1=[];sampLat2=[];sampAmp1=[];sampAmp2=[];if ~isempty(EPdataIn.facVecT)    msg{1}=['Error: Temporal PCA data cannot undergo a sample-by-sample test.'];    [msg]=ep_errorMsg(msg);    returnend;[EPdataIn]=ep_stripAdds(EPdataIn);EEGchans=find(strcmp('EEG',EPdataIn.chanTypes));numChans=length(EEGchans);numPoints=length(EPdataIn.timeNames);if ~isempty(templateWaveform)    numTemplatePoints=length(templateWaveform);else    numTemplatePoints=1;end;numSubs=length(EPdataIn.subNames);numFreqs=max(1,length(EPdataIn.freqNames));numCells=length(EPdataIn.cellNames);sampleSize=(1000/EPdataIn.Fs);if numPoints==0    theDim='freq';else    theDim='time';end;if ~strcmp(EPdataIn.dataType,'continuous')        cellNames=unique(EPdataIn.cellNames);        if length(cellNames) == 2        theMode='cell';    elseif (numChans == 2) && (length(unique(cellNames)) == 1)        theMode='chan';    else        msg{1}=['Error: There must be two cells or two channels present to conduct a sample-by-sample test.'];        [msg]=ep_errorMsg(msg);        return    end;        if ~strcmp(EPdataIn.dataType,'average') && strcmp(test,'jack-knife')        msg{1}=['Error: Jack-knife analyses currently limited to averaged data.'];        [msg]=ep_errorMsg(msg);        return    end;        sigPoint=.25;    contSigPoint=.50;    maxSigPoint=.75;        if numPoints==0        freqNum=1;    else        freqNum=numFreqs;    end;        for iFreq=1:freqNum        if freqNum==1            freqList=[];        else            freqList=iFreq;        end;        switch theMode            case 'cell'                sample1=squeeze(ep_expandFacs(EPdataIn,EEGchans,[],strcmp(cellNames{1},EPdataIn.cellNames),[],[],freqList,[]));                sample2=squeeze(ep_expandFacs(EPdataIn,EEGchans,[],strcmp(cellNames{2},EPdataIn.cellNames),[],[],freqList,[]));                if strcmp(theDim,'freq')                    sample1=permute(sample1,[1 3 2]);                    sample2=permute(sample2,[1 3 2]);                end;                if strcmp(EPdataIn.dataType,'single_trial')                    sample1=sample1(:,:,EPdataIn.analysis.badTrials(1,strcmp(cellNames{1},EPdataIn.cellNames))==0);                    sample2=sample2(:,:,EPdataIn.analysis.badTrials(1,strcmp(cellNames{2},EPdataIn.cellNames))==0);                else %average data                    goodCells=(EPdataIn.avgNum(:,strcmp(cellNames{1},EPdataIn.cellNames))~=-1) & (EPdataIn.avgNum(:,strcmp(cellNames{2},EPdataIn.cellNames))~=-1);                    sample1=sample1(:,:,goodCells);                    sample2=sample2(:,:,goodCells);                end;            case 'chan'                sample1(1,:,:)=squeeze(ep_expandFacs(EPdataIn,1,[],[],[],[],freqList,[]));                sample2(1,:,:)=squeeze(ep_expandFacs(EPdataIn,2,[],[],[],[],freqList,[]));                if strcmp(theDim,'freq')                    sample1=permute(sample1,[1 3 2]);                    sample2=permute(sample2,[1 3 2]);                end;                if strcmp(EPdataIn.dataType,'single_trial')                    sample1=sample1(1,:,EPdataIn.analysis.badTrials(1,:)==0);                    sample2=sample2(1,:,EPdataIn.analysis.badTrials(1,:)==0);                else %average data                    goodCells=(EPdataIn.avgNum(:,:)~=-1);                    sample1=sample1(1,:,goodCells);                    sample2=sample2(1,:,goodCells);                end;        end;                if (size(sample1,3)<2) || (size(sample2,3)<2)            msg{1}=['Error: Not enough good data to perform the analysis.'];            [msg]=ep_errorMsg(msg);            return        end;                switch method            case 'sample'                chanNum=numChans;            case 'CWT'                chanNum=numChans;                %will use ms as the common metric                sampSizeMs=1000/EPdataIn.Fs;                tMs=[0:sampSizeMs:(numPoints*sampSizeMs)-sampSizeMs]; %the time range in ms, counting from first sample of epoch rather than from stimulus-onset                                fprintf('%60s\n',' ' );                for iTau=1:length(tMs)                    fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d of %4d',['Applying CWT to sample # '], iTau, numPoints))                    theTau=tMs(iTau);                    psi=(1-16.*(((tMs-theTau)/(scale/1000))/1000).^2).*exp(-8.*(((tMs-theTau)/(scale/1000))/1000).^2); %The Mexican Hat wavelet                    for iChan=1:numChans                        for iTrial=1:size(sample1,3)                            theEpoch=squeeze(sample1(iChan,:,iTrial));                            W1(iChan,iTau,iTrial)=(1/sqrt(scale/1000)).*theEpoch*psi';                        end;                        for iTrial=1:size(sample2,3)                            theEpoch=squeeze(sample2(iChan,:,iTrial));                            W2(iChan,iTau,iTrial)=(1/sqrt(scale/1000)).*theEpoch*psi';                        end;                    end;                end;                fprintf('%60s\n',' ' );                sample1=W1;                sample2=W2;            case {'Template Woody','PCA Woody'}                                if ~isempty(templateWaveform)                    theTemplate=templateWaveform';                    [a peakPoint]=max(abs(templateWaveform));                else                    theTemplate=1;                    peakPoint=1;                end;                if ~isempty(templateTopo)                    chanNum=1;                    theTemplate=templateTopo*theTemplate;                else                    chanNum=numChans;                end;                                %make sure templates are same size as the data                                fprintf('%60s\n',' ' );                for iPoint=1:numPoints                    fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d of %4d',['Applying Woody filter to sample # '], iPoint, numPoints));                    pointTemplate=zeros(size(theTemplate,1),numPoints);                    if ~isempty(templateWaveform)                        pointTemplate(:,max(iPoint-peakPoint+1,1):min(numPoints-peakPoint+iPoint,numPoints))=theTemplate(:,max(peakPoint+1-iPoint,1):min(numPoints-iPoint+peakPoint,numPoints));                    else                        pointTemplate(:,iPoint)=theTemplate;                    end;                    if ~isempty(templateTopo)                        for iTrial=1:size(sample1,3)                            theEpoch=squeeze(sample1(:,:,iTrial));                            W1(:,iPoint,iTrial)=mean(mean((theEpoch.*pointTemplate(EEGchans,:))'));                        end;                        for iTrial=1:size(sample2,3)                            theEpoch=squeeze(sample2(:,:,iTrial));                            W2(:,iPoint,iTrial)=mean(mean((theEpoch.*pointTemplate(EEGchans,:))'));                        end;                    else                        for iChan=1:numChans                            for iTrial=1:size(sample1,3)                                theEpoch=squeeze(sample1(iChan,:,iTrial));                                W1(iChan,iPoint,iTrial)=mean(abs(mean((theEpoch.*pointTemplate)')));                            end;                            for iTrial=1:size(sample2,3)                                theEpoch=squeeze(sample2(iChan,:,iTrial));                                W2(iChan,iPoint,iTrial)=mean(abs(mean((theEpoch.*pointTemplate)')));                            end;                        end;                    end;                end;                fprintf('%60s\n',' ' );                if isempty(templateTopo)                    for iChan=1:numChans                        sample1(iChan,:,:)=W1(1,:,:);                        sample2(iChan,:,:)=W2(1,:,:);                    end;                else                    sample1=W1;                    sample2=W2;                end;            otherwise                msg{1}=['Programming Error: Selected method not recognized.'];                [msg]=ep_errorMsg(msg);                return        end;        if strcmp(theMode,'chan')            chanNum=1;        end;                %         fprintf('%160s\n',' ' );        outputData=zeros(chanNum,numPoints,numFreqs);        sigResults=zeros(chanNum,numPoints,numFreqs);                switch test            case 't-test'                if strcmp(EPdataIn.dataType,'single_trial')                    for iChan=1:chanNum                        %fprintf('%s%-160s',repmat(sprintf('\b'),1,160),sprintf('%s%4d of %4d',['Applying sample-by-sample Mann-Whitney U-test to channel# '], iChan, numChans))                        if strcmp(theDim,'time')                            for iPoint=1:numPoints                                [sigResults(iChan,iPoint,iFreq), h]=ranksum(squeeze(sample1(iChan,iPoint,:)),squeeze(sample2(iChan,iPoint,:)));                            end;                        else                            for iPoint=1:numFreqs                                [sigResults(iChan,1,iPoint), h]=ranksum(squeeze(sample1(iChan,1,:)),squeeze(sample2(iChan,iPoint,:)));                            end;                        end;                    end;                else %averaged data                    for iChan=1:chanNum                        %                         fprintf('%s%-160s',repmat(sprintf('\b'),1,160),sprintf('%s%4d of %4d',['Applying sample-by-sample Wilcoxon sign rank t-test to channel# '], iChan, numChans))                        if strcmp(theDim,'time')                            for iPoint=1:numPoints                                [sigResults(iChan,iPoint,iFreq), h]=signrank(squeeze(sample1(iChan,iPoint,:)),squeeze(sample2(iChan,iPoint,:)));                            end;                        else                            for iPoint=1:numFreqs                                [sigResults(iChan,1,iPoint), h]=signrank(squeeze(sample1(iChan,1,:)),squeeze(sample2(iChan,iPoint,:)));                            end;                        end;                    end;                end;                %fprintf('%160s\n',' ' );            case 'jack-knife'                if strcmp(EPdataIn.dataType,'single_trial')                    msg{1}=['Error: Jack-knife analyses currently limited to averaged data.'];                    [msg]=ep_errorMsg(msg);                    outputData=[];                    return                else %averaged data                    jack=zeros(numSubs,1);                    badSubs=zeros(numSubs,1);                    for iSub=1:numSubs                        fprintf('%s%-160s',repmat(sprintf('\b'),1,160),sprintf('%s%4d of %4d',['Applying jack-knife test to subsample# '], iSub, numSubs))                        subjectList=setdiff([1:numSubs],iSub);                        jackGAV1=squeeze(mean(sample1(channel,:,subjectList),3));                        jackGAV2=squeeze(mean(sample2(channel,:,subjectList),3));                        if thresh >= 0                            onsetSample=min(find(jackGAV1 >= thresh));                        else                            onsetSample=min(find(jackGAV1 <= thresh));                        end;                        if isempty(onsetSample)                            badSubs(iSub)=1;                            continue                        end;                        if onsetSample==1                            jack(iSub,1)=EPdataIn.timeNames(1);                        else                            jack(iSub,1)=EPdataIn.timeNames(onsetSample)+(sampleSize*(thresh-jackGAV1(onsetSample-1))/(jackGAV1(onsetSample)-jackGAV1(onsetSample-1)));                        end;                        if thresh >= 0                            onsetSample=min(find(jackGAV2 >= thresh));                        else                            onsetSample=min(find(jackGAV2 <= thresh));                        end;                        if isempty(onsetSample)                            badSubs(iSub)=1;                            continue                        end;                        if onsetSample==1                            jack(iSub,2)=EPdataIn.timeNames(1);                        else                            jack(iSub,2)=EPdataIn.timeNames(onsetSample)+(sampleSize*(thresh-jackGAV2(onsetSample-1))/(jackGAV2(onsetSample)-jackGAV2(onsetSample-1)));                        end;                    end;                    fprintf('%160s\n',' ' );                    jack(find(badSubs))=[];                    if size(jack,1) < 2                        msg{1}=['Error: Too few averages met the threshold to compute a jack-knife statistic.'];                        [msg]=ep_errorMsg(msg);                        outputData=[];                        return                    end;                    J=mean(jack(:,2)-jack(:,1));                    sD=0;                    numSubSamples=size(jack,1);                    for iSub=1:numSubSamples                        sD=sD+(jack(iSub,2)-jack(iSub,1)-J)^2;                    end;                    sD=sqrt(((numSubSamples-1)/numSubSamples)*sD);                                        %calculate onset difference in grand average                    jackGAV1=squeeze(mean(sample1(channel,:,find(~badSubs)),3));                    if thresh >= 0                        onsetSample=min(find(jackGAV1 >= thresh));                    else                        onsetSample=min(find(jackGAV1 <= thresh));                    end;                    if isempty(onsetSample)                        msg{1}=['Error: First condition does not meet the threshold.'];                        [msg]=ep_errorMsg(msg);                        outputData=[];                        return                    end;                    if onsetSample==1                        onset1=EPdataIn.timeNames(1);                    else                        onset1=EPdataIn.timeNames(onsetSample)+(sampleSize*(thresh-jackGAV1(onsetSample-1))/(jackGAV1(onsetSample)-jackGAV1(onsetSample-1)));                    end;                    jackGAV2=squeeze(mean(sample2(channel,:,find(~badSubs)),3));                    if thresh >= 0                        onsetSample=min(find(jackGAV2 >= thresh));                    else                        onsetSample=min(find(jackGAV2 <= thresh));                    end;                    if isempty(onsetSample)                        msg{1}=['Error: Second condition does not meet the threshold.'];                        [msg]=ep_errorMsg(msg);                        outputData=[];                        return                    end;                    if onsetSample==1                        onset2=EPdataIn.timeNames(1);                    else                        onset2=EPdataIn.timeNames(onsetSample)+(sampleSize*(thresh-jackGAV2(onsetSample-1))/(jackGAV2(onsetSample)-jackGAV2(onsetSample-1)));                    end;                    D=onset2-onset1;                    tJ=D/sD;                    disp(['The jack-knife onset is a difference of ' num2str(D) ' with a t-value of ' num2str(tJ) ', ' num2str(numSubSamples-1) ' degrees of freedom.']);                    if ft_hastoolbox('STATS', 0, 1)                        disp(['It has a two-sided p-value of ' num2str(2*tcdf(-abs(tJ),numSubSamples-1)) '.']);                    end;                end;            otherwise                msg{1}=['Programming Error: Selected test not recognized.'];                [msg]=ep_errorMsg(msg);                outputData=[];                return        end;                %contiguity criterion        for iChan=1:chanNum            contigFlag=0;            if strcmp(theDim,'freq')                for iPoint=1:numFreqs                    if sigResults(iChan,1,iPoint) <= alpha                        outputData(iChan,1,iPoint)=sigPoint;                        contigFlag=contigFlag+1;                        if contigFlag >= contiguous                            outputData(iChan,1,iPoint)=contSigPoint;                            if contigFlag == contiguous %if just passed threshold for continguity, update prior samples too                                outputData(iChan,1,iPoint-contigFlag+1:iPoint)=contSigPoint;                            end;                        end;                    else                        contigFlag=0;                    end;                    [sigLevel,Index] = min(sigResults(iChan,1,iPoint));                    if sigLevel <= alpha                        outputData(iChan,1,Index)=maxSigPoint;                    end;                end;            else                for iPoint=1:numPoints                    if sigResults(iChan,iPoint,iFreq) <= alpha                        outputData(iChan,iPoint,iFreq)=sigPoint;                        contigFlag=contigFlag+1;                        if contigFlag >= contiguous                            outputData(iChan,iPoint,iFreq)=contSigPoint;                            if contigFlag == contiguous %if just passed threshold for continguity, update prior samples too                                outputData(iChan,iPoint-contigFlag+1:iPoint,iFreq)=contSigPoint;                            end;                        end;                    else                        contigFlag=0;                    end;                end;                [sigLevel,Index] = min(sigResults(iChan,:,iFreq));                if sigLevel <= alpha                    outputData(iChan,Index,iFreq)=maxSigPoint;                end;            end;        end;                if ~isempty(templateTopo) && any(strcmp(method,{'Template Woody','PCA Woody'}))            outputData=repmat(outputData(1,:,iFreq),numChans,1); %if spatial factor, then apply results to all the channels            [sampAmp1(iFreq), sampLat1(iFreq)]=max(abs(squeeze(W1(1,:,:)))); %extract latency and amplitude of maximum absolute cross-product for Sample 1            [sampAmp2(iFreq), sampLat2(iFreq)]=max(abs(squeeze(W2(1,:,:)))); %extract latency and amplitude of maximum absolute cross-product for Sample 2            sampLat1(iFreq)=(sampLat1-EPdataIn.baseline)*(1000/EPdataIn.Fs);            sampLat2(iFreq)=(sampLat2-EPdataIn.baseline)*(1000/EPdataIn.Fs);        end;    end;else %continuous    if ~isempty(templateWaveform)        theTemplate=templateWaveform';        [a peakPoint]=max(abs(templateWaveform));    else        theTemplate=1;        peakPoint=1;    end;    if ~isempty(templateTopo)        theTemplate=templateTopo*theTemplate;    end;        outputData=zeros(1,numPoints);        fprintf('%60s\n',' ' );    for iFreq=1:numFreqs        for iPoint=1:numPoints            fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d of %4d',['Applying Woody filter to sample # '], iPoint, numPoints));            %pointTemplate(:,max(iPoint-peakPoint+1,1):min(numTemplatePoints-peakPoint+iPoint,numPoints))=theTemplate(:,max(peakPoint+1-iPoint,1):min(numPoints-iPoint+peakPoint,numTemplatePoints));            d1=max(iPoint-peakPoint+1,1);            d2=min(numTemplatePoints-peakPoint+iPoint,numPoints);            t1=max(peakPoint+1-iPoint,1);            t2=min(numPoints-iPoint+peakPoint,numTemplatePoints);            if ~isempty(templateTopo)                outputData(1,iPoint)=mean(mean((EPdataIn.data(EEGchans,d1:d2).*theTemplate(EEGchans,t1:t2))'));            else                outputData(1,iPoint)=mean(abs(mean((EPdataIn.data(channel,d1:d2).*theTemplate(1,t1:t2))')));            end;        end;    end;    fprintf('%60s\n',' ' );end;disp('Finished sample-by-sample test.');