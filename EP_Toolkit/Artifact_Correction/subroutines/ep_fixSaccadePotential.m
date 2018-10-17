function [totsaccadeTrialNum, templates, outputLog, graphCounter] = ep_fixSaccadePotential(inFile, outFile, startChunk, endChunk, badDataCriteria, badChans, eog, templateSource, saccadeFile, refChan, excludePoints, butterflyFig, graphCounter, numGraphs, theSubject)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [totsaccadeTrialNum, templates, outputLog, graphCounter] = ep_fixSaccadePotential(inFile, outFile, startChunk, endChunk, badDataCriteria, badChans, eog, templateSource, saccadeFile, refChan, excludePoints, butterflyFig, graphCounter, numGraphs, theSubject)%%	Reads in file chunks generated by chunkInputFile function and corrects and marks saccade potentials.%   A canonical saccade potential is used to automatically construct a custom template based on all the chunks which is then used%   to perform a full correction run.%%Inputs%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.%	outFile:    filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.%	startChunk: starting chunk (usually 1)%   endChunk:   ending chunk%   badDataCriteria:  Criteria for detecting bad data.%       .window:    moving average window for smoothing%       .minmax:    difference from minimum to maximum for bad channel%       .trialminmax:  difference from minimum to maximum for bad trial%       .badnum:    percent of bad channels exceeded to declare bad trial, rounding down%       .hminmax:   difference from minimum to maximum for bad horizontal EOG%       .neighbors: number of electrodes considered to be neighbors%       .badchan:   maximum microvolt difference allowed from best matching neighbor%       .maxneighbor:   maximum microvolt difference allowed from best matching neighbor%       .blink:     threshold correlation with blink template, 0 to 1%       .saccade:     threshold correlation with saccade template, 0 to 1%       .saccademin:  �v Saccade Fac is the minimum HEOG voltage difference required to constitute a possible saccade.%       .detrend:   1 to detrend%       .badtrials: percentage of good trials chan is bad to declare a channel globally bad%       .replace:   1 to interpolate bad channels from neighbors.%       .noadjacent:1 to not allow adjacent bad channels (trial or subject declared bad)%       .movefacs  : number of factors to retain during movement correction.%       .channelMode: 'replace' to interpolate bad channels, 'mark' to mark them with a spike, and 'none' to do nothing.%       .trialMode: 'fix' to fix bad trial data and 'none' to do nothing.%       .saturation: followed by range of acceptable data values.  Time points with a channel outside this range will be excluded.%   badChans:   list of global bad channels to exclude from saccade potential detection process.%   eog:        EOG channels.%   templateSource:   source of saccade templates (fileTemplate: load file.  autoTemplate: automatically generate saccade potential template.%                  bothTemplate: check automatic template and then manual template).%   saccadeFile:  file with saccade templates.  Assumed to be in the same directory as the data file.%   refChan:    Array of current reference channels.%   excludePoints: time points to exclude, numbering end-to-end for single-trial data (cell array of number of chunks)%   butterflyFig:  the handle for the output figure.  Otherwise, will open a new figure.%   graphCounter: the current subplot for the summary figure.%   numGraphs: the total number of subgraphs in the summary figure.%   theSubject: which subject of the file is being processed.%%   The input chunks are EP format data files.%%Outputs%	Saves files with saccade potential removed, replacing the original chunked files.%   totsaccadeTrialNum: Total list of saccade trials.%   templates: the templates used for the correction, both file and auto if both used.  Returns empty if error occurs.%       .sacPot%           .manual: manual template from file (number of EEG channels,1)%           .auto: auto template (number of EEG channels,1)%   outputLog: output messages from saccade potential fixing process%   graphCounter: the current subplot for the summary figure.%% History:%% by Joseph Dien (2/5/17)% jdien07@mac.com%% modified 4/18/17 JD% Excludes time points beyond a certain range from global bad channel detection, blink, and saccade routines.%% bugfix 10/4/17 JD% Median corrects data prior to saturation check to ensure channels with merely high offsets are not treated as bad data.% Saccade Potential correction works with sampling rates other than 250 Hz.%   `% bugfix 10/5/17 JD% Fixed crashes due to Matlab changing their graphics objects AGAIN in 2017b.%% bugfix 10/20/17 JD% Eliminated x tick labels to address problem with subplots in summary% artifact figure getting squeezed by formatting problem on PCs.%% bugfix 10/5/17 JD% Fixed crash when preprocessing data containing impedance values.% Fixed crash when preprocessing multi-subject average files.%% bugfix 12/13/17 JD% Fixed saccade potential preprocessing crashing when channels marked as being missing in preprocessing preferences with a -1.% Fixed bugs causing saccade potential to not work as effectively for single-trial data.%% modified 2/4/18 JD% Made subplot specification for summary figure output more flexible.%% modified 4/8/18 JD% Consolidated summary figure for average files so no longer one per subject.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Copyright (C) 1999-2018  Joseph Dien%%     This program is free software: you can redistribute it and/or modify%     it under the terms of the GNU General Public License as published by%     the Free Software Foundation, either version 3 of the License, or%     (at your option) any later version.%%     This program is distributed in the hope that it will be useful,%     but WITHOUT ANY WARRANTY; without even the implied warranty of%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the%     GNU General Public License for more details.%%     You should have received a copy of the GNU General Public License%     along with this program.  If not, see <http://www.gnu.org/licenses/>.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%templates=[];msg='Starting saccade potentials routine.';disp(msg);outputLog{1}=msg;totsaccadeTrialNum=[];badChans=badChans(:); %make sure it's a column vectortemplateSPthresh=badDataCriteria.sacpot;SPminSep=100; %ms minimum separation between saccade potentialsif ~exist('butterflyFig','var')    butterflyFig=figure('Name','Artifact Correction','NumberTitle','off');    colormap jet;    standAlone=1;else    standAlone=0;end;if ~isempty(butterflyFig)    graphCounter=graphCounter+2;end;theEOG(1) = eog.LUVEOG;theEOG(2) = eog.RUVEOG;theEOG(3) = eog.LLVEOG;theEOG(4) = eog.RLVEOG;theEOG(5) = eog.LHEOG;theEOG(6) = eog.RHEOG;VEOG=[];goodVEOG=[];for i=1:4    if isempty(intersect(theEOG(i),badChans)) && (theEOG(i) ~= -1)        VEOG=[VEOG theEOG(i)];        goodVEOG(end+1)=i;    end;end;blinksign = [1 1 -1 -1]';blinksign=blinksign(goodVEOG,1);HEOG=[];for i=5:6    if isempty(intersect(theEOG(i),badChans)) && (theEOG(i) ~= -1)        HEOG=[HEOG theEOG(i)];    end;end;eval(['load ''' deblank(inFile) '''-' num2str(startChunk) '.mat']);numChans=length(dataChunk.chanNames);numSubs=length(dataChunk.subNames);EEGchans=find(strcmp('EEG',dataChunk.chanTypes));chans = setdiff(EEGchans,badChans);if dataChunk.Fs < 200    msg='The sampling rate needs to be at least 200 Hz to correct saccade potentials.';    outputLog{end+1}=msg;    disp(' ');    disp('**************************************************************');    disp(msg);    disp('**************************************************************');    disp(' ');    returnend;eloc=dataChunk.eloc;[theChan, theOrder] = ep_findChan(eloc, badChans, [0 -9.7989 107.9359]);if theChan ==0    msg='No good chans so aborting saccade potential correction.';    outputLog{end+1}=msg;    disp(' ');    disp('**************************************************************');    disp(msg);    disp('**************************************************************');    disp(' ');    returnend;VEOG(end+1)=theChan;blinksign(end+1)=0;numVEOG=length(VEOG);if theOrder > 1    msg=['Cz is a bad channel so instead using channel ' eloc(VEOG(end)).labels '.'];elseif strcmpi(eloc(VEOG(end)).labels,'cz')    msg=['Channel Cz identified.'];else    msg=['Channel ' eloc(VEOG(end)).labels ' is assumed to be Cz.'];end;outputLog{end+1}=msg;disp(msg);cannonicalSPthresh=25; %microvolt threshold for detecting an SP when constructing the automatic templatefourMS=ceil(4/(1000/(dataChunk.Fs))); %how many samples provides at least 4 msoutputLog=[];SPmanualTemplate=[];SPautoTemplateL=[];SPautoTemplateR=[];if any(strcmp(templateSource, {'fileTemplate','bothTemplate'}))    [fileDir, name, ext] = fileparts(inFile);    eval(['load ''' saccadeFile '''']);        if ~exist('EPsaccade','var')        errMsg{1}='Not a saccade template.';        [msg]=ep_errorMsg(errMsg);        return    end;        if ~isfield(EPsaccade,'sacPot')        errMsg{1}='No saccade potential template in saccades template file.';        [msg]=ep_errorMsg(errMsg);        return    end;        SPmanualTemplate=EPsaccade.sacPot.template;    numTemplateEEG=length(EPsaccade.sacPot.template);    numDataEEG=length(EEGchans);    if numTemplateEEG ~= numDataEEG        msg=['Number of saccade potential template EEG electrodes (' num2str(numTemplateEEG) ') different from the data (' num2str(numDataEEG) ').'];        outputLog{end+1}=msg;        disp(' ');        disp('**************************************************************');        disp(msg);        disp('**************************************************************');        disp(' ');        return;    end;        if any([EPsaccade.eloc.theta]-[dataChunk.eloc(EEGchans).theta])        msg='Saccade potential template electrode locations not consistent with the data.';        outputLog{end+1}=msg;        disp(' ');        disp('**************************************************************');        disp(msg);        disp('**************************************************************');        disp(' ');        return;    end;end;if any(strcmp(templateSource, {'autoTemplate','bothTemplate','eyeTracker'}))    priorPoints=0;    for iChunk = startChunk:endChunk        if endChunk > startChunk            msg=[deblank(inFile) '-' num2str(iChunk)];            disp(msg);            outputLog{end+1}=msg;        end;        if iChunk > startChunk            eval(['load ''' deblank(inFile) '''-' num2str(iChunk) '.mat']);        end;                if strcmp(dataChunk.dataType,'continuous')            theSegment = 'one second epoch';        else            theSegment = 'trial';        end;                if length(dataChunk.facNames) > 1            msg='This function is not intended for application to factor data.';            outputLog{end+1}=msg;            disp(' ');            disp('**************************************************************');            disp(msg);            disp('**************************************************************');            disp(' ');            return;        end;                trialdata=dataChunk.data(:,:,:,theSubject);        contData=reshape(dataChunk.data(:,:,:,theSubject),numChans,[]);                if strcmp(dataChunk.dataType,'continuous')            displayPeriod=size(trialdata,2);    %Number of timepoints to graph in display.        else            displayPeriod=size(trialdata,2)*size(trialdata,3);        end;        decimateSamples=ceil(max(1,displayPeriod/10000));        totalDisplayPeriod=displayPeriod*size(dataChunk.data,4);                if displayPeriod == 1            msg='There is only one time point and so the data cannot be saccade potential corrected.';            outputLog{end+1}=msg;            disp(' ');            disp('**************************************************************');            disp(msg);            disp('**************************************************************');            disp(' ');            return;        end;                if iChunk ==1            autoCountL=0;            autoCountR=0;            SPautoTemplateL=zeros(length(chans),1);            SPautoTemplateR=zeros(length(chans),1);        end;                if strcmp(dataChunk.dataType,'continuous')            numTrials=1; %excess time points are tacked onto final epoch            trialSize = length(dataChunk.timeNames);        else            trialSize = length(dataChunk.timeNames);            numTrials = length(dataChunk.cellNames);        end;                goodPoints = find((max(contData(chans,:)-repmat(median(contData(chans,:)')',1,size(contData(chans,:),2))) < badDataCriteria.saturation(2)) & (min(contData(chans,:)-repmat(median(contData(chans,:)')',1,size(contData(chans,:),2))) > badDataCriteria.saturation(1)));                goodPoints=setdiff(goodPoints,excludePoints{iChunk});                %automatically construct a custom template        if strcmp(templateSource, 'eyeTracker')            saccadeSamps=round([dataChunk.events{1}(find(strcmp('saccadeET',{dataChunk.events{1}.value}))).sample]);            saccadeSamps=saccadeSamps-priorPoints;            saccadeSamps=saccadeSamps((saccadeSamps > 0) & (saccadeSamps <= size(contData,2))); %only saccadeET events in the current data chunk.            SPautoTemplate=mean(contData(chans,saccadeSamps+1*fourMS),2)-mean(contData(chans,saccadeSamps),2);            autoCount=1;        elseif length(HEOG) ~= 2            disp('Need both HEOG channels to form custom template');        else            for iTrial=1:size(dataChunk.data,3)                boundaryPoints=[dataChunk.events{iTrial}(find(strcmp('boundary',{dataChunk.events{iTrial}.type}))).sample];                diffData=diff(dataChunk.data(VEOG,:,iTrial,theSubject),1,2);                theData=squeeze(diffData(:,:,1,1))-repmat(diffData(numVEOG,:,1,1),numVEOG,1); %rereference to Cz                theData=theData(1:numVEOG-1,:); %drop Cz, which is now a flat line                pointList=[];                trialGoodPoints=find((goodPoints >= (((iTrial-1)*trialSize)+1)) & (goodPoints <= (iTrial*trialSize)));                for iPoint=1:size(theData,2)-2                    if ismember(iPoint,trialGoodPoints) && ismember(iPoint+(2*fourMS),trialGoodPoints) && ~any(ismember([iPoint:iPoint+(2*fourMS)],boundaryPoints))                        if all(theData(:,iPoint) <= -cannonicalSPthresh) && all(theData(:,iPoint+(2*fourMS)) >= cannonicalSPthresh) %in the VEOG channels, look for a negative spike where iPoint is just as it starts and iPoint+2 is the peak (diff = Xt1 - Xt)                            if isempty(pointList) || (pointList(end) < (iPoint-(fourMS)*2))                                putativeSP=(trialdata(chans,iPoint+fourMS,iTrial,theSubject)-trialdata(VEOG(end),iPoint+fourMS,iTrial,theSubject))-(trialdata(chans,iPoint,iTrial,theSubject)-trialdata(VEOG(end),iPoint,iTrial,theSubject)); %Cz referenced                                if all(putativeSP < 100) %must be less than 100 uv to exclude blinks and other artifacts                                    EOGx=min([dataChunk.eloc(goodVEOG(1:end-1)).X]);                                    if max(putativeSP([goodVEOG(1:end-1)])) < min(putativeSP(find(ismember(chans,find(EOGx>[dataChunk.eloc.X])))))                                        pointList(end+1)=iPoint;                                        thePoints=[4*fourMS:12*fourMS]+iPoint;                                        SPdir=[ones(length(thePoints),1),thePoints']\(trialdata(eog.LHEOG,thePoints,iTrial)-trialdata(eog.RHEOG,thePoints,iTrial))'; %CRD is positive towards direction of gaze                                        if SPdir(2) > 0                                            autoCountL=autoCountL+1;                                            SPautoTemplateL=SPautoTemplateL+(trialdata(chans,iPoint+fourMS,iTrial,theSubject)-trialdata(chans,iPoint,iTrial,theSubject));                                        else                                            autoCountR=autoCountR+1;                                            SPautoTemplateR=SPautoTemplateR+(trialdata(chans,iPoint+fourMS,iTrial,theSubject)-trialdata(chans,iPoint,iTrial,theSubject));                                        end;                                    end;                                end;                            end;                        end;                    end;                end;            end;        end;        priorPoints=priorPoints+size(trialdata,2);    end;    if (autoCountL == 0) && (autoCountR == 0)        msg='No automatic template could be formed.';        outputLog{end+1}=msg;        disp(' ');        disp('**************************************************************');        disp(msg);        disp('**************************************************************');        disp(' ');        if any(strcmp(templateSource, {'autoTemplate','eyeTracker'}))            if ~isempty(butterflyFig)                if numSubs > 1                    theTitle='subtracted saccade potentials';                else                    theTitle='no saccade potentials to subtract that were detected';                end;                plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,zeros(size(contData)),EEGchans,theSubject);                figure(butterflyFig(iChunk));                subplot(numGraphs,1,graphCounter-2), plot([1:decimateSamples:totalDisplayPeriod],plotData);                title(theTitle,'Interpreter','none');                axis([1 totalDisplayPeriod -200 200])                set(gca,'XTickLabel','','XTick',[]);                                if numSubs > 1                    theTitle='with saccade potentials subtracted';                else                    theTitle='no saccade potentials subtracted';                end;                plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,contData,EEGchans,theSubject);                subplot(numGraphs,1,graphCounter-1), plot([1:decimateSamples:totalDisplayPeriod],plotData);                title(theTitle,'Interpreter','none');                axis([1 totalDisplayPeriod -200 200])                set(gca,'XTickLabel','','XTick',[]);            end;            return;        elseif strcmp(templateSource, 'bothTemplate')            templateSource='fileTemplate';        end;    else        if autoCountL == 0            msg='No left automatic template could be formed so will use flipped right template.';            outputLog{end+1}=msg;                        autoCountL=autoCountR;            oldElocs=dataChunk.eloc(chans);            newElocs=oldElocs;            for iEloc=1:length(newElocs)                newElocs(iEloc).theta=-newElocs(iEloc).theta; %flip L/R            end;            [SPautoTemplateL]=ep_interpChans(SPautoTemplateR, oldElocs, newElocs);        end;        if autoCountR == 0            msg='No right automatic template could be formed so will use flipped left template.';            outputLog{end+1}=msg;                        autoCountR=autoCountL;            oldElocs=dataChunk.eloc(chans);            newElocs=oldElocs;            for iEloc=1:length(newElocs)                newElocs(iEloc).theta=-newElocs(iEloc).theta; %flip L/R            end;            [SPautoTemplateR]=ep_interpChans(SPautoTemplateL, oldElocs, newElocs);        end;        SPautoTemplateL=SPautoTemplateL./autoCountL;        SPautoTemplateR=SPautoTemplateR./autoCountR;        SPautoTemplateL=SPautoTemplateL-SPautoTemplateL(VEOG(end));        SPautoTemplateR=SPautoTemplateR-SPautoTemplateR(VEOG(end));    end;end;%use the template to better detect saccade potentialssubtractedSacPotTopo=zeros(length(EEGchans),1);for iChunk = startChunk:endChunk        if startChunk ~= endChunk        eval(['load ''' deblank(inFile) '''-' num2str(iChunk) '.mat']);    end;        if length(dataChunk.facNames) > 1        msg='This function is not intended for application to factor data.';        outputLog{end+1}=msg;        disp(' ');        disp('**************************************************************');        disp(msg);        disp('**************************************************************');        disp(' ');        return;    end;        trialdata=dataChunk.data(:,:,:,theSubject);    contData=reshape(dataChunk.data(:,:,:,theSubject),numChans,[]);        if strcmp(dataChunk.dataType,'continuous')        displayPeriod=size(trialdata,2);    %Number of timepoints to graph in display.    else        displayPeriod=size(trialdata,2)*size(trialdata,3);    end;    totalDisplayPeriod=displayPeriod*size(dataChunk.data,4);    decimateSamples=ceil(max(1,totalDisplayPeriod/10000));    if displayPeriod == 1        msg='There is only one time point and so the data cannot be saccade potential corrected.';        outputLog{end+1}=msg;        disp(' ');        disp('**************************************************************');        disp(msg);        disp('**************************************************************');        disp(' ');        return;    end;        if strcmp(dataChunk.dataType,'continuous')        numTrials=floor(size(trialdata,2)/ceil(dataChunk.Fs)); %excess time points are tacked onto final epoch        trialSize = min(ceil(dataChunk.Fs),size(trialdata,2)); %one second epochs    else        trialSize = length(dataChunk.timeNames);        numTrials = length(dataChunk.cellNames);    end;        goodPoints = find((max(contData(chans,:)-repmat(median(contData(chans,:)')',1,size(contData(chans,:),2))) < badDataCriteria.saturation(2)) & (min(contData(chans,:)-repmat(median(contData(chans,:)')',1,size(contData(chans,:),2))) > badDataCriteria.saturation(1)));        goodPoints=setdiff(goodPoints,excludePoints{iChunk});        count=0;    SPwaveAL=zeros(1,size(trialdata,2),size(trialdata,3));    SPwaveAR=zeros(1,size(trialdata,2),size(trialdata,3));    SPwaveM=zeros(1,size(trialdata,2),size(trialdata,3));    subtractedSPs =zeros(size(trialdata));    sacPotEOG=zeros(1,size(trialdata,2),size(trialdata,3));    EEGchansNoHEOG=setdiff(chans,[eog.LHEOG eog.RHEOG]);    saccadeTrialNum{iChunk}=dataChunk.analysis.saccadeTrial;    saccadeOnsetNum{iChunk}=dataChunk.analysis.saccadeOnset;    if any(strcmp(templateSource, {'autoTemplate','bothTemplate','eyeTracker'}))        saccadesScaleAL=abs(SPautoTemplateL(VEOG(end))-mean(SPautoTemplateL(VEOG(blinksign == -1))));        saccadesScaleAR=abs(SPautoTemplateR(VEOG(end))-mean(SPautoTemplateR(VEOG(blinksign == -1))));    end;    if any(strcmp(templateSource,{'fileTemplate','bothTemplate'}))        saccadesScaleM=abs(SPmanualTemplate(VEOG(end))-mean(SPmanualTemplate(VEOG(blinksign == -1))));    end;    for iTrial=1:size(dataChunk.data,3)        if any(strcmp(templateSource, {'autoTemplate','bothTemplate','eyeTracker'}))            SPwaveAL(1,1:end-1,iTrial)=(-SPautoTemplateL(EEGchansNoHEOG)'*squeeze(diff(dataChunk.data(EEGchansNoHEOG,:,1,theSubject)-repmat(dataChunk.data(VEOG(end),:,1,theSubject),length(EEGchansNoHEOG),1),1,2)))/(length(EEGchansNoHEOG)*saccadesScaleAL);  %drop HEOG channels from template as they are too variable            SPwaveAR(1,1:end-1,iTrial)=(-SPautoTemplateR(EEGchansNoHEOG)'*squeeze(diff(dataChunk.data(EEGchansNoHEOG,:,1,theSubject)-repmat(dataChunk.data(VEOG(end),:,1,theSubject),length(EEGchansNoHEOG),1),1,2)))/(length(EEGchansNoHEOG)*saccadesScaleAR);  %drop HEOG channels from template as they are too variable        end;        if any(strcmp(templateSource,{'fileTemplate','bothTemplate'}))            SPwaveM(1,1:end-1,iTrial)=(-SPmanualTemplate(EEGchansNoHEOG)'*squeeze(diff(dataChunk.data(EEGchansNoHEOG,:,1,theSubject),1,2)))/(length(EEGchansNoHEOG)*saccadesScaleM);  %drop HEOG channels from template as they are too variable        end;        boundaryPoints=[dataChunk.events{iTrial}(find(strcmp('boundary',{dataChunk.events{iTrial}.type}))).sample];        pointList2=[];        for iPoint=1:size(trialdata,2)-(6*fourMS)-1            if all(ismember(iPoint:iPoint+(3*fourMS),goodPoints)) && ~any(ismember(iPoint:iPoint+(3*fourMS),boundaryPoints))                SPdur=0; %duration of the SP in 4ms increments                if any(strcmp(templateSource, {'autoTemplate','bothTemplate','eyeTracker'}))                    if SPwaveAL(iPoint) < SPwaveAR(iPoint)                        SPwaveA=SPwaveAL;                        SPTemplate=SPautoTemplateL;                    else                        SPwaveA=SPwaveAR;                        SPTemplate=SPautoTemplateR;                    end;                    if (SPwaveA(iPoint) <= -(templateSPthresh/2)) && (SPwaveA(iPoint+(2*fourMS)) >= templateSPthresh) %the second time points has a larger threshold because the onset of the saccade makes the voltages more variable so more conservative threshold needed.                        SPdur=2;                    end;                    if (SPwaveA(iPoint) <= -(templateSPthresh/2)) && (SPwaveA(iPoint+(1*fourMS)) >= templateSPthresh)                        SPdur=1;                    end;                end;                if any(strcmp(templateSource,{'fileTemplate','bothTemplate'})) && (SPdur == 0)                    if (SPwaveM(iPoint) <= -(templateSPthresh/2)) && (SPwaveM(iPoint+(2*fourMS)) >= templateSPthresh)                        SPdur=2;                    end;                    if (SPwaveM(iPoint) <= -(templateSPthresh/2)) && (SPwaveM(iPoint+(1*fourMS)) >= templateSPthresh)                        SPdur=1;                    end;                    SPTemplate=SPmanualTemplate;                end;                if SPdur > 0                    B=[ones(1,3*fourMS+1);1:3*fourMS+1]'\[mean(trialdata(VEOG(blinksign == 1),iPoint+1+(3*fourMS):iPoint+1+(6*fourMS),iTrial),1)]'; %correction factor for effects of vertical saccades on the saccade potential                    %ensure the upper VEOG channels actually have the appropriate temporal morphology                    if mean(mean(trialdata(VEOG(blinksign == 1),[iPoint iPoint+1+(SPdur*fourMS)],iTrial))) - mean(mean(trialdata(VEOG(blinksign == 1),iPoint+1:iPoint+(SPdur*fourMS),iTrial))) -B(2)/2  >= (templateSPthresh/10)                        if isempty(pointList2) || (pointList2(end) < (iPoint-(fourMS*SPminSep/4))) %enforce minimum spacing for SacPot events                            pointList2(end+1)=iPoint;                            count=count+1;                            saccadeTrialNum{iChunk}(theSubject,iTrial)=1;                            if ~saccadeOnsetNum{iChunk}(theSubject,iTrial)                                saccadeOnsetNum{iChunk}(theSubject,iTrial)=iPoint+1; %assumes no more than one saccade per trial.  If more than only the first one will be recorded.                            end;                            if strcmp(dataChunk.dataType,'continuous')                                dataChunk.events{theSubject,iTrial}(end+1).sample=iPoint+1;                                dataChunk.events{theSubject,iTrial}(end).type='artifact';                                dataChunk.events{theSubject,iTrial}(end).value='SacPot';                                dataChunk.events{theSubject,iTrial}(end).duration=0;                                dataChunk.events{theSubject,iTrial}(end).keys=struct('code','','data','','datatype','','description','');                            else                                dataChunk.events{theSubject,iTrial}(end+1).sample=iPoint+1;                                dataChunk.events{theSubject,iTrial}(end).type='artifact';                                dataChunk.events{theSubject,iTrial}(end).value='SacPot';                                dataChunk.events{theSubject,iTrial}(end).duration=0;                                dataChunk.events{theSubject,iTrial}(end).keys=struct('code','','data','','datatype','','description','');                            end;                            for iDur=1:SPdur*fourMS                                %SPamp=mean(diff(trialdata(VEOG(blinksign == 1),[iPoint iPoint+iDur]),1,2))/mean(SPTemplate(VEOG(blinksign == 1)));                                SPamp=[ones(length(EEGchansNoHEOG),1),SPTemplate(find(EEGchansNoHEOG))]\[diff(trialdata(EEGchansNoHEOG,[iPoint iPoint+iDur],iTrial),1,2)];                                trialdata(chans,iPoint+iDur,iTrial)=trialdata(chans,iPoint+iDur,iTrial)-(SPTemplate*SPamp(2)+SPamp(1));                                sacPotEOG(1,iPoint+iDur,iTrial)=SPamp(2);                                subtractedSPs(chans,iPoint+iDur,iTrial)=SPTemplate*SPamp(2)+SPamp(1);                            end;                        end;                    end;                end;            end;        end;    end;        if ~isempty(pointList2)        if count == 1            msg=['One saccade potential corrected.'];        else            msg=[num2str(count) ' saccade potentials corrected.'];        end        disp(msg);        outputLog{end+1}=msg;                if ~isempty(butterflyFig)            figure(butterflyFig(iChunk));            theTitle='subtracted saccade potentials';            plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,subtractedSPs,EEGchans,theSubject);            subplot(numGraphs,1,graphCounter-2), plot([1:decimateSamples:totalDisplayPeriod],plotData);            title(theTitle,'Interpreter','none');            axis([1 totalDisplayPeriod -200 200])            set(gca,'XTickLabel','','XTick',[]);                        theTitle='with saccade potentials subtracted';            plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,contData,EEGchans,theSubject);            subplot(numGraphs,1,graphCounter-1), plot([1:decimateSamples:totalDisplayPeriod],plotData);            title('with saccade potentials subtracted','Interpreter','none');            axis([1 totalDisplayPeriod -200 200])            set(gca,'XTickLabel','','XTick',[]);        end;    else        msg='No components match saccade potential template so no correction performed.';        disp(msg);        outputLog{end+1}=msg;                if ~isempty(butterflyFig)            figure(butterflyFig(iChunk));            if numSubs > 1                theTitle='subtracted saccade potentials';            else                theTitle='no saccade potentials to subtract that were detected';            end;            plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,zeros(size(contData)),EEGchans,theSubject);            subplot(numGraphs,1,graphCounter-2), plot([1:decimateSamples:totalDisplayPeriod],plotData);            axis([1 totalDisplayPeriod -200 200])            set(gca,'XTickLabel','','XTick',[]);            title('no saccade potentials to subtract that were detected','Interpreter','none');                        if numSubs > 1                theTitle='with saccade potentials subtracted';            else                theTitle='no saccade potentials subtracted';            end;            plotData=ep_makePlotData(butterflyFig(iChunk),displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,contData,EEGchans,theSubject);            subplot(numGraphs,1,graphCounter11), plot([1:decimateSamples:totalDisplayPeriod],plotData);            title('no saccade potentials subtracted','Interpreter','none');            axis([1 totalDisplayPeriod -200 200])            set(gca,'XTickLabel','','XTick',[]);        end;    end;        drawnow        goodPoints=[1:size(contData,2)];        %rereference the data.    if ~isempty(refChan)        for iTrial = 1:numTrials            epoch=goodPoints(find((goodPoints>((iTrial-1)*trialSize)) & (goodPoints<=(iTrial*trialSize))));            if ~isempty(epoch)                referenceData=mean(trialdata(refChan,epoch),1);                for iChan=1:length(EEGchans)                    theChan=EEGchans(iChan);                    trialdata(theChan,epoch)=trialdata(theChan,epoch)-referenceData;                end;            end;        end;    end;        dataChunk.analysis.saccadeTrial(theSubject,:)=saccadeTrialNum{iChunk};	dataChunk.analysis.saccadeOnset(theSubject,:)=saccadeOnsetNum{iChunk};    dataChunk.data(:,:,:,theSubject)=trialdata;        SacPotChan=find(strcmp('SacPot',dataChunk.chanNames));    if isempty(SacPotChan)        EPadd.chanNames{1}='SacPot';        EPadd.chanTypes{1}='REG';        EPadd.data=zeros(length(EPadd.chanNames),length(dataChunk.timeNames),length(dataChunk.cellNames),length(dataChunk.subNames));        EPadd.data(:,:,:,theSubject)=sacPotEOG;        if isfield(dataChunk,'interpChans')            interpChans=dataChunk.interpChans;            dataChunk=rmfield(dataChunk,'interpChans'); %remove interpChans temporarily so it doesn't trigger error by checkEPfile.            [dataChunk]=ep_addData(dataChunk,EPadd,'channels');            dataChunk.interpChans=interpChans;        else            [dataChunk]=ep_addData(dataChunk,EPadd,'channels');        end;    else        dataChunk.data(saccChan,:,:,theSubject)=sacPotEOG;    end;        if isempty(dataChunk)        disp('Warning: No file saved due to program error.');    end;    eval (['save ''' outFile '''-' num2str(iChunk) '.mat dataChunk;']);    if standAlone        try            MATLABver=ver('MATLAB');            [a b]=strtok(MATLABver.Version,'.');            b=b(2:end);            if ~isprop(butterflyFig,'Number')                eval (['print -f' num2str(butterflyFig(iChunk)) ' -djpeg ''' inFile '-' num2str(iChunk) 'EOG.jpg''']);            else                eval (['print -f' num2str(butterflyFig(iChunk).Number) ' -djpeg ''' inFile '-' num2str(iChunk) 'EOG.jpg''']);            end;        catch            disp('Couldn''t save a copy of the artifact correction figure.  Perhaps your version of Matlab is not current.');        end;    end;    subtractedSacPotTopo=subtractedSacPotTopo+mean(subtractedSPs(EEGchans,:),2);end;for iChunk = startChunk:endChunk    totsaccadeTrialNum=[totsaccadeTrialNum saccadeTrialNum{iChunk}];end;if any(strcmp(templateSource, {'fileTemplate','bothTemplate'}))    templates.sacPot.manual=SPmanualTemplate;else    templates.sacPot.manual=[];end;if any(strcmp(templateSource, {'autoTemplate','bothTemplate','eyeTracker'}))    templates.sacPot.autoL(chans,1)=SPautoTemplateL;    templates.sacPot.autoL=templates.sacPot.autoL(EEGchans);    templates.sacPot.autoR(chans,1)=SPautoTemplateR;    templates.sacPot.autoR=templates.sacPot.autoR(EEGchans);else    templates.sacPot.autoL=[];    templates.sacPot.autoR=[];end;templates.sacPot.sacPotTopo=subtractedSacPotTopo;