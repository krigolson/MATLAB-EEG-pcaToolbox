function [EPdataOut]=ep_combineData(EPdataIn,dataDimension,combineData,combineWeights,combineName,trialWeights);
%  [EPdataOut]=ep_combineData(EPdataIn,dataDimension,combineData,combineWeights,combineName);
%       Combines data levels (like a subset of subjects) and adds the result to the data.
%
%Inputs:
%  EPdataIn       : Structured array with the input data and accompanying information in EP file format.  See readData.
%  dataDimension  : Which data dimension to combine (i.e., 'channels', 'cells', 'subjects', 'factors')
%  combineData    : Which levels of the data dimension to combine (e.g., [1 3 4])
%  combineWeights : Weights to use when combining the levels of the data (e.g., [.5 .25 .25]).
%                   If totals more than one, weights will be rescaled so that they do total to one.
%                   If weights total zero (as in a difference wave) then no rescaling.
%                   For factors, no rescaling is applied so simple addition rather than averages are computed.
%  combineName    : Name of the new combined level (e.g., 'grand average')
%  trialWeights   : Weight by number of trials in average if available, for cells (0=no, 1=yes)
%
%Outputs:
%  EPdataOut      : Structured array with the output data and accompanying information in EP file format.  See readData.

%History:
%  by Joseph Dien (5/24/12)
%  jdien07@mac.com
%
%  modified 7/17/12 JD
%  Added option to weight cell combinations by number of trials in averages.
%
%  bugfix 8/6/12 JD
%  Fixed edit's add cells trial weighting option not working correctly when the cells are not a consecutive series starting with the first.
%
%  bugfix 11/4/12 JD
%  Fixed crash when using combining cells or chans for spatial PCA data.
%
%  bugfix 1/10/13 JD
%  Fixed crash when combining channels under certain circumstances.
%
%  bugfix 1/30/13 JD
%  Fixed crash when combining channels for factor data under certain circumstances (presence of facData due to adds).
%
%  bugfix 3/25/13 JD
%  Fixed combining of subjects and chans not correct when weights not the same (as in difference wave).
%  Fixed weighting of difference waves for cells and chans and subjects incorrect (waves too small).
%
%  bugfix 7/14/13 JD
%  Fixed combining channels results in flat waveform.
%
% modified 10/9/13 JD
% Added recTime field.
%
%  bugfix 11/25/13 JD
%  Fixed noise and std fields set equal to the data when combining cells or subjects.
%
%  bugfix 11/28/13 JD
%  If new add fails data check then output empty matrix.
%
% modified 3/19/14 JD
% Added combining factors
%
% bugfix 3/20/14 JD
% Makes sure that additions to the name fields are added as a column vector.
% Fixed bad subtype when adding a single subject.
%
%  bugfix 3/22/14 JD
%  Fixed all but one channel is flat for grand average combined factors, as in the "all" factor from PCAs.
%  Fixed all but one channel is flat for combined cells if one already has a combined factor, as when one uses the Edit
%  function to combined cells on a factor cell containing an "all" cell.
%
% modified 3/24/14 JD
% Added .cov field.
%
% bugfix 4/2/14 JD
% Fixed crash when performing combination of subjects with file containing .cov information.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% bugfix 4/2/14 JD
% Fixed weighting not correct when computing difference waves that do not sum to zero.
% Fixed calculation of the SubNum field (number of subjects going into averages) when combining cells.
%
% modified 4/24/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% bugfix 6/1/14 JD
% Fixed crash when combining trials.
%
% modified 6/1/14 JD
% Added simple averaging of trials in cell subpane for single-trial data.
%
% bugfix 6/29/14 JD
% Fixed cov.Nq field not being formed correctly when combining subjects, resulting in crashes later on.
%
% bugfix 4/27/15 JD
% Fixed error when adding cells to single-trial data.
%
% modified 9/4/15 JD
% Added trial specs for average files.
%
% bugfix 10/24/15 JD
% Fixed crash when combining cells or subjects and there are trial specs with numbers.
%
% bugfix 6/18/17 JD
% Fixed not combining numeric trial specs correctly over subjects, resulting in crashes down the line in ANOVA function.
% Switch to amplitude scaling when adding freq data together other than channels.
%
% bugfix 7/2/17 JD
% Fixed crash when combining subjects for spatial PCA data.
%
% modified 11/22/17 JD
% Added support for impedances field.
%
% bugfix 11/25/17 JD
% Fixed crash when conducting PCA on data where there are trial spec fields that are identical across all the trials/cells and are characters.
%
% bugfix 12/9/17 JD
% Fixed crash when combining cells and the data is frequency-domain.
%
% modified 2/11/18 JD
% Changed std field of combined subjects (as via the Edit function) to be std of the newly generated grand average data rather than a combination of their std values.
% Changed noise field of combined subjects (as via the Edit function) to be noise of the newly generated grand average data rather than a combination of their noise values.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% bugfix 4/8/18 JD
% Fixed bug when combining subjects with trial specs where they are all the same string value across the subjects, resulting in failure to combine.
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

EPdataOut=EPdataIn;

if nargin < 6
    trialWeights=0;
end

numChans=length(EPdataIn.chanNames);
numPoints=length(EPdataIn.timeNames);
numCells=length(EPdataIn.cellNames);
numSubs=length(EPdataIn.subNames);
numFacs=length(EPdataIn.facNames);
numFreqs=length(EPdataIn.freqNames);
numRels=length(EPdataIn.relNames);
if numFacs==0
    numFacs=1;
end;
if ~isempty(EPdataIn.facData)
    numCMBfacs=size(EPdataIn.facData,5);
else
    numCMBfacs=0;
end;
numSGLfacs=numFacs-numCMBfacs;

combineData=combineData(:);
combineWeights=combineWeights(:);

switch dataDimension
    case 'channels'
        chanList=combineData;
        newChanData=zeros(size(EPdataIn.data(1,:,:,:,:,:,:)));
        if ~isempty(EPdataIn.noise)
            newNoiseData=zeros(size(EPdataIn.noise(1,:,:,:,:,:,:)));
        end;
        if ~isempty(EPdataIn.std)
            newStdData=zeros(size(EPdataIn.std(1,:,:,:,:,:,:)));
        end;
        if ~isempty(EPdataIn.stdCM)
            newStdCMData=zeros(size(EPdataIn.stdCM(1,:,:,:,:,:,:)));
        end;
        
        if ~isempty(EPdataIn.facVecT)
            newChanData=newChanData(:,1,:,:,:,:,:);
        end;
        if ~isempty(EPdataIn.facVecF)
            newChanData=newChanData(:,:,:,:,:,1,:);
        end;
        newFacVecSData=zeros(1,numSGLfacs);
        newBadChanData=zeros(numSubs,numCells,1);
        newFacData=zeros(1,max(numPoints,1),numCells,numSubs,numCMBfacs,max(numFreqs,1),max(numRels,1));
        
        %no need to skip bad cells and subjects since the EP Toolkit will ignore them anyway and maybe they'll still be of interest
        for iSub=1:numSubs
            for theCell=1:numCells
                if strcmp(EPdataIn.dataType,'average')
                    goodChans=find(~isnan(EPdataIn.analysis.badChans(iSub,theCell,:)));
                else
                    goodChans=find(EPdataIn.analysis.badChans(iSub,theCell,:) >= 0);
                end;
                chanList=intersect(chanList,goodChans);
                totalWeight=sum(combineWeights(find(ismember(combineData,chanList))));
                if totalWeight==0
                    totalWeight=1;
                else
                    totalWeight=sum(abs(combineWeights(find(ismember(combineData,chanList)))));
                end;
                if any(ismember(combineData,chanList))
                    for chanLoop=1:length(chanList)
                        theChan=chanList(chanLoop);
                        if goodChans
                            theWeight=combineWeights(find(ismember(combineData,theChan)));
                            if ~isempty(EPdataIn.facVecS)
                                newFacVecSData=newFacVecSData+(theWeight/totalWeight)*EPdataIn.facVecS(theChan,:);
                            else
                                newChanData(1,:,theCell,iSub,:,:,:)=newChanData(1,:,theCell,iSub,:,:,:)+(theWeight/totalWeight)*EPdataIn.data(theChan,:,theCell,iSub,:,:,:);
                            end;
                            if ~isempty(EPdataIn.facData)
                                newFacData(1,:,theCell,iSub,:,:,:)=newFacData(1,:,theCell,iSub,:,:,:)+(theWeight/totalWeight)*EPdataIn.facData(theChan,:,theCell,iSub,:,:,:);
                            end;
                            if ~isempty(EPdataIn.noise)
                                newNoiseData(1,:,theCell,iSub,:,:)=newNoiseData(1,:,theCell,iSub,:,:)+(theWeight/totalWeight)*EPdataIn.noise(theChan,:,theCell,iSub,:,:);
                            end;
                            newBadChanData(iSub,theCell,1)=newBadChanData(iSub,theCell,1)+(theWeight/totalWeight)*EPdataIn.analysis.badChans(iSub,theCell,theChan);
                        end;
                    end;
                else %no good data available for this regional channel
                    if strcmp(EPdataIn.dataType,'average')
                        newBadChanData(iSub,theCell,1)=NaN;
                    else
                        newBadChanData(iSub,theCell,1)=-1;
                    end;
                end;
            end;
        end;
        if ~isempty(EPdataIn.facVecS)
            EPdataOut.facVecS(end+1,:)=newFacVecSData;
        else
            EPdataOut.data(end+1,:,:,:,:,:,:)=newChanData;
        end;
        
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData(end+1,:,:,:,:,:,:)=newFacData;
        end;
        
        EPdataOut.analysis.badChans(:,:,end+1)=newBadChanData;
        
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise(end+1,:,:,:,:)=newNoiseData;
        end;
        if ~isempty(EPdataIn.std)
            EPdataOut.std(end+1,:,:,:,:,:)=newStdData;
        end;
        if ~isempty(EPdataIn.stdCM)
            EPdataOut.stdCM(end+1,:,:,:,:,:)=newStdCMData;
        end;
        if ~isempty(EPdataIn.cov)
            EPdataOut.cov.covMatrix(:,end+1,end+1)=NaN;
        end;
        numChans=numChans+1;
        EPdataOut.chanNames{end+1,1}=combineName;
        EPdataOut.chanTypes{end+1,1}='REG';
        if ~isempty(EPdataIn.eloc)
            EPdataOut.eloc(end+1)=cell2struct(cell(1,length(fieldnames(EPdataIn.eloc))),fieldnames(EPdataIn.eloc),2);
        end;
        
        if ~isempty(EPdataIn.relNames)
            numRels=numRels+1; %the new channel
            EPdataOut.relNames{end+1,1}=combineName;
            EPdataOut.data(1:end-1,:,:,:,:,:,end+1)=permute(newChanData,[7 2 3 4 5 6 1]);
            EPdataOut.data(end,:,:,:,:,:,end)=1; %coherence of channel with itself is real number one.
        end;
        
        if ~isempty(EPdataIn.impedances.channels)
        	EPdataOut.impedances.channels(end+1,:)=NaN;
        end;
        
    case 'cells'
        newCellData=zeros(size(EPdataIn.data(:,:,1,:,:,:,:)));
        if ~isempty(EPdataIn.noise)
            newCellNoise=zeros(size(EPdataIn.noise(:,:,1,:,:,:,:)));
        end;
        if ~isempty(EPdataIn.std)
            newCellStd=zeros(size(EPdataIn.std(:,:,1,:,:,:,:)));
        end;
        if ~isempty(EPdataIn.stdCM)
            newCellStdCM=zeros(size(EPdataIn.stdCM(:,:,1,:,:,:,:)));
        end;
        newFacData=zeros(size(EPdataIn.facData(:,:,1,:,:,:,:)));
        if ~isempty(EPdataIn.facVecT)
            newCellData=newCellData(:,1,:,:,:,:,:);
        end;
        if ~isempty(EPdataIn.facVecF)
            newCellData=newCellData(:,:,:,:,:,1,:);
        end;
        if ~isempty(EPdataIn.facVecS)
            newCellData=newCellData(1,:,:,:,:,:,:);
        end;
        EPdataOut.analysis.badChans(:,numCells+1,:)= 0;
        for iSub=1:numSubs
            theOkayCells=combineData(find(EPdataOut.avgNum(iSub,combineData)>=0));
            if isempty(theOkayCells)
                EPdataOut.avgNum(iSub,numCells+1)=-1;
                EPdataOut.covNum(iSub,numCells+1)=-1;
                EPdataOut.subNum(iSub,numCells+1)=-1;
                EPdataOut.analysis.blinkTrial(iSub,numCells+1)= 0;
                EPdataOut.analysis.saccadeTrial(iSub,numCells+1)= 0;
                EPdataOut.analysis.saccadeOnset(iSub,numCells+1)= 0;
                EPdataOut.analysis.moveTrial(iSub,numCells+1)= 0;
                EPdataOut.analysis.badTrials(iSub,numCells+1)= 0;
                EPdataOut.analysis.badChans(iSub,numCells+1,:)= 0;
                EPdataOut.recTime(numCells+1)=1;
                
            else
                EPdataOut.avgNum(iSub,numCells+1)=sum(EPdataIn.avgNum(iSub,theOkayCells));
                if ~isempty(EPdataOut.covNum)
                    %assume cov matrix is the same so just need to figure out the new effective sample size for the new combination
                    %per p.128 of the 3.7.2 MNE manual, 1/Leff=Sigma weight-squared/L
                    covWeights=combineWeights(find(ismember(combineData,theOkayCells)));
                    EPdataOut.covNum(iSub,numCells+1)=sum([EPdataIn.covNum(iSub,theOkayCells)'.^-1].*covWeights.^2)^-1;
                end;
                EPdataOut.subNum(iSub,numCells+1)=max(EPdataIn.subNum(iSub,theOkayCells));
                EPdataOut.analysis.blinkTrial(iSub,numCells+1)= sum(EPdataIn.analysis.blinkTrial(iSub,theOkayCells));
                EPdataOut.analysis.saccadeTrial(iSub,numCells+1)= sum(EPdataIn.analysis.saccadeTrial(iSub,theOkayCells));
                EPdataOut.analysis.saccadeOnset(iSub,numCells+1)= sum(EPdataIn.analysis.saccadeOnset(iSub,theOkayCells));
                EPdataOut.analysis.moveTrial(iSub,numCells+1)= sum(EPdataIn.analysis.moveTrial(iSub,theOkayCells));
                EPdataOut.analysis.badTrials(iSub,numCells+1)= sum(EPdataIn.analysis.badTrials(iSub,theOkayCells));
                EPdataOut.recTime(numCells+1)=min(EPdataIn.recTime(theOkayCells));
            end;
            if trialWeights && ~any(EPdataIn.avgNum(iSub,theOkayCells) == 0) && strcmp(EPdataIn.dataType,'average')
                useTrialWeights=1;
            else
                useTrialWeights=0;
            end;
            for iChan=1:numChans
                if strcmp(EPdataIn.dataType,'average')
                    goodCellList=theOkayCells(find(~isnan(EPdataIn.analysis.badChans(iSub,theOkayCells,iChan))));
                else
                    goodCellList=theOkayCells(find(EPdataIn.analysis.badChans(iSub,theOkayCells,iChan) >= 0));
                end;
                if ~isempty(goodCellList)
                    if useTrialWeights
                        totalWeight=sum(combineWeights(find(ismember(combineData,goodCellList)))'.*EPdataIn.avgNum(iSub,combineData(find(ismember(combineData,goodCellList)))));
                        if totalWeight==0
                            totalWeight=1;
                        else
                            totalWeight=sum(abs(combineWeights(find(ismember(combineData,goodCellList))))'.*EPdataIn.avgNum(iSub,combineData(find(ismember(combineData,goodCellList)))));
                        end;
                    else
                        totalWeight=sum(combineWeights(find(ismember(combineData,goodCellList))));
                        if totalWeight==0
                            totalWeight=1;
                        else
                            totalWeight=sum(abs(combineWeights(find(ismember(combineData,goodCellList)))));
                        end;
                    end;
                    for iCell=1:length(goodCellList)
                        theCell=goodCellList(iCell);
                        if useTrialWeights
                            theWeight=combineWeights(find(ismember(combineData,goodCellList(iCell))))*EPdataIn.avgNum(iSub,combineData(find(ismember(combineData,goodCellList(iCell)))));
                        else
                            theWeight=combineWeights(find(ismember(combineData,goodCellList(iCell))));
                        end;
                        if (iChan==1) || isempty(EPdataIn.facVecS)
                            theData=EPdataIn.data(iChan,:,theCell,iSub,:,:,:);
                            if ~isempty(EPdataIn.freqNames) && any(strcmp(EPdataIn.chanTypes{iChan},{'EEG','REG'}))
                                theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                            end;
                            newCellData(iChan,:,1,iSub,:,:,:)=newCellData(iChan,:,1,iSub,:,:,:)+(theWeight/totalWeight)*theData;
                        end;
                        if ~isempty(EPdataIn.facData)
                            theData=EPdataIn.facData(iChan,:,theCell,iSub,:,:,:);
                            if ~isempty(EPdataIn.freqNames) && any(strcmp(EPdataIn.chanTypes{iChan},{'EEG','REG'}))
                                theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                            end;
                            newFacData(iChan,:,1,iSub,:,:,:)=newFacData(iChan,:,1,iSub,:,:,:)+(theWeight/totalWeight)*theData;
                        end;
                        if ~isempty(EPdataIn.noise)
                            newCellNoise(iChan,:,1,iSub,:,:)=newCellNoise(iChan,:,1,iSub,:,:)+(theWeight/totalWeight)*EPdataIn.noise(iChan,:,theCell,iSub,:,:);
                        end;
                        EPdataOut.analysis.badChans(iSub,numCells+1,iChan)= EPdataOut.analysis.badChans(iSub,numCells+1,iChan)+EPdataIn.analysis.badChans(iSub,theCell,iChan);
                    end;
                else %none of the cells have good data for this channel
                    if strcmp(EPdataIn.dataType,'average')
                        EPdataOut.analysis.badChans(iSub,numCells+1,iChan)=NaN;
                    else
                        EPdataOut.analysis.badChans(iSub,numCells+1,iChan)=-1;
                    end;
                end;
            end;
        end;
        EPdataOut.data(:,:,end+1,:,:,:,:)=newCellData;
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData(:,:,end+1,:,:,:,:)=newFacData;
        end;
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise(:,:,end+1,:,:,:)=newCellNoise;
        end;
        if ~isempty(EPdataIn.std)
            EPdataOut.std(:,:,end+1,:,:,:)=newCellStd;
        end;
        if ~isempty(EPdataIn.stdCM)
            EPdataOut.stdCM(:,:,end+1,:,:,:)=newCellStdCM;
        end;
        
        if length(combineData) > 1
            EPdataOut.cellTypes{end+1,1}='CMB';
        else
            EPdataOut.cellTypes{end+1}=EPdataIn.cellTypes(combineData);
        end;
        
        if ~isempty(EPdataIn.trialNames)
            EPdataOut.trialNames(end+1,1)=1;
        end;
        
        numCells=numCells+1;
        EPdataOut.cellNames{end+1,1}=combineName;
        
        if ~isempty(EPdataIn.trialSpecs)
            %generate the averaged trial specs
            for iSpec=1:length(EPdataIn.trialSpecNames)
                for iSub=1:length(EPdataIn.subNames)
                    if ischar(EPdataIn.trialSpecs{combineData(1),iSpec,iSub}) && all(strcmp(EPdataIn.trialSpecs{combineData(1),iSpec,iSub},EPdataIn.trialSpecs(combineData,iSpec,iSub)))
                        EPdataOut.trialSpecs{numCells,iSpec,iSub}=EPdataIn.trialSpecs{1,iSpec,iSub};
                    else
                        theSpecs=EPdataIn.trialSpecs(combineData,iSpec,iSub);
                        theStrSpecs=find(cellfun(@ischar,theSpecs));
                        for iStrSpec=1:length(theStrSpecs)
                            theStrSpec=theStrSpecs(iStrSpec);
                            theSpecs{theStrSpec}=str2double(theSpecs{theStrSpec});
                        end;
                        EPdataOut.trialSpecs{numCells,iSpec,iSub}=mean(cell2mat(theSpecs));
                    end;
                end;
            end;
        else
            EPdataOut.trialSpecs=cell(numCells,0,numSubs);
        end;
        
        EPdataOut.events(:,end+1)=cell(size(EPdataIn.events,1),1);
    case 'subjects'
        newCellData=zeros(size(EPdataIn.data(:,:,:,1,:,:,:)));
        if isempty(EPdataIn.noise)
            EPdataOut.noise=zeros(max(numChans,1),max(numPoints,1),numCells,numSubs,numFacs,max(numFreqs,1),max(numRels,1));
        end;        
        if isempty(EPdataIn.std)
            EPdataOut.std=zeros(max(numChans,1),max(numPoints,1),numCells,numSubs,numFacs,max(numFreqs,1),max(numRels,1));
        end;
        if isempty(EPdataIn.stdCM)
            EPdataOut.stdCM=zeros(max(numChans,1),max(numPoints,1),numCells,numSubs,numFacs,max(numFreqs,1),max(numRels,1));
        end;
        newCellNoise=zeros(size(EPdataOut.noise(:,:,:,1,:,:,:)));
        newCellStd=zeros(size(EPdataOut.std(:,:,:,1,:,:,:)));
        newCellStdCM=zeros(size(EPdataOut.stdCM(:,:,:,1,:,:,:)));
        newFacData=zeros(size(EPdataIn.facData(:,:,:,1,:,:,:)));
        newCellCov=zeros(1,numChans,numChans);
        if ~isempty(EPdataIn.facVecT)
            newCellData=newCellData(:,1,:,:,:,:,:);
        end;
        if ~isempty(EPdataIn.facVecF)
            newCellData=newCellData(:,:,:,:,:,1,:);
        end;
        if ~isempty(EPdataIn.facVecS)
            newCellData=newCellData(1,:,:,:,:,:,:);
        end;
        EPdataOut.analysis.badChans(numSubs+1,:,:)= 0;
        EPdataOut.subjectSpecs(end+1,:)=cell(1,size(EPdataIn.subjectSpecs,2));
        
        for iCell=1:numCells
            theOkaySubs=combineData(find(EPdataOut.avgNum(combineData,iCell)>=0));
            if isempty(theOkaySubs)
                EPdataOut.avgNum(numSubs+1,iCell)=-1;
                EPdataOut.covNum(numSubs+1,iCell)=-1;
                EPdataOut.subNum(numSubs+1,iCell)=-1;
                EPdataOut.analysis.blinkTrial(numSubs+1,iCell)= 0;
                EPdataOut.analysis.saccadeTrial(numSubs+1,iCell)= 0;
                EPdataOut.analysis.saccadeOnset(numSubs+1,iCell)= 0;
                EPdataOut.analysis.moveTrial(numSubs+1,iCell)= 0;
                EPdataOut.analysis.badTrials(numSubs+1,iCell)= 0;
                EPdataOut.analysis.badChans(numSubs+1,iCell,:)= 0;
                
            else
                EPdataOut.avgNum(numSubs+1,iCell)=sum(EPdataIn.avgNum(theOkaySubs,iCell));
                
                if ~isempty(EPdataOut.covNum)
                    %assume cov matrix is different so need to figure out new covariance matrix as well as the new effective sample size for the new combination
                    %per p.128 of the 3.7.2 MNE manual, 1/Leff=Sigma weight-squared/L
                    covWeights=combineWeights(find(ismember(combineData,theOkaySubs)));
                    EPdataOut.covNum(numSubs+1,iCell)=sum([EPdataIn.covNum(theOkaySubs,iCell).^-1].*covWeights.^2)^-1;
                end;
                
                EPdataOut.subNum(numSubs+1,iCell)=sum(EPdataIn.subNum(theOkaySubs,iCell));
                EPdataOut.analysis.blinkTrial(numSubs+1,iCell)= sum(EPdataIn.analysis.blinkTrial(theOkaySubs,iCell));
                EPdataOut.analysis.saccadeTrial(numSubs+1,iCell)=sum(EPdataIn.analysis.saccadeTrial(theOkaySubs,iCell));
                EPdataOut.analysis.saccadeOnset(numSubs+1,iCell)=sum(EPdataIn.analysis.saccadeOnset(theOkaySubs,iCell));
                EPdataOut.analysis.moveTrial(numSubs+1,iCell)=sum(EPdataIn.analysis.moveTrial(theOkaySubs,iCell));
                EPdataOut.analysis.badTrials(numSubs+1,iCell)=sum(EPdataIn.analysis.badTrials(theOkaySubs,iCell));
                for iChan=1:numChans
                    if strcmp(EPdataIn.dataType,'average')
                        goodSubList=theOkaySubs(find(~isnan(EPdataIn.analysis.badChans(theOkaySubs,iCell,iChan))));
                    else
                        goodSubList=theOkaySubs(find((EPdataIn.analysis.badChans(theOkaySubs,iCell,iChan) >= 0)));
                    end;
                    if ~isempty(goodSubList)
                        totalWeight=sum(combineWeights(find(ismember(combineData,goodSubList))));
                        if totalWeight==0
                            totalWeight=1;
                        else
                            totalWeight=sum(abs(combineWeights(find(ismember(combineData,goodSubList)))));
                        end;
                        if isempty(EPdataIn.facVecS)
                            theoldChan=iChan;
                        else
                            theoldChan=1;
                        end;
                        for iFactor=1:numSGLfacs;
                            newCellStd(iChan,:,iCell,1,iFactor,:)=std(squeeze(EPdataIn.data(theoldChan,:,iCell,goodSubList,iFactor,:,:))');
                            numGoodSubs=length(goodSubList);
                            theSigns=ones(numGoodSubs,1);
                            if numGoodSubs > 1
                                theSigns(2:2:end)=-1;
                            end;
                            for iSub=1:numGoodSubs
                                newCellNoise(iChan,:,iCell,1,iFactor,:)=newCellNoise(iChan,:,iCell,1,iFactor,:)+mean(theSigns(iSub)*EPdataIn.data(theoldChan,:,iCell,goodSubList(iSub),iFactor,:,:),4);
                            end;
                        end;
                        newCellNoise=newCellNoise/numGoodSubs;
                        for iFactor=1:numCMBfacs;
                            newCellStd(iChan,:,iCell,1,numSGLfacs+iFactor,:)=std(squeeze(EPdataIn.facData(iChan,:,iCell,goodSubList,iFactor,:,:))');
                            numGoodSubs=length(goodSubList);
                            theSigns=ones(numGoodSubs,1);
                            if numGoodSubs > 1
                                theSigns(2:2:end)=-1;
                            end;
                            for iSub=1:numGoodSubs
                                newCellNoise(iChan,:,iCell,1,numSGLfacs+iFactor,:)=newCellNoise(iChan,:,iCell,1,numSGLfacs+iFactor,:)+mean(theSigns(iSub)*EPdataIn.facData(iChan,:,iCell,goodSubList(iSub),iFactor,:,:),4);
                            end;
                        end;
                        newCellNoise=newCellNoise/numGoodSubs;
                        for iSub=1:length(goodSubList)
                            theSub=goodSubList(iSub);
                            theWeight=combineWeights(find(ismember(combineData,goodSubList(iSub))));
                            if (iChan==1) || isempty(EPdataIn.facVecS)
                                theData=EPdataIn.data(iChan,:,iCell,theSub,:,:,:);
                                if ~isempty(EPdataIn.freqNames) && any(strcmp(EPdataIn.chanTypes{iChan},{'EEG','REG'}))
                                    theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                                end;
                                newCellData(iChan,:,iCell,1,:,:,:)=newCellData(iChan,:,iCell,1,:,:,:)+(theWeight/totalWeight)*theData;
                            end;
                            if ~isempty(EPdataIn.facData)
                                theData=EPdataIn.facData(iChan,:,iCell,theSub,:,:,:);
                                if ~isempty(EPdataIn.freqNames) && any(strcmp(EPdataIn.chanTypes{iChan},{'EEG','REG'}))
                                    theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                                end;
                                newFacData(iChan,:,iCell,1,:,:,:)=newFacData(iChan,:,iCell,1,:,:,:)+(theWeight/totalWeight)*theData;
                            end;
                            
                            if ~isempty(EPdataIn.cov)
                                newCellCov(1,iChan,:)=newCellCov(1,iChan,:)+(EPdataIn.cov.covMatrix(theSub,iChan,:)*EPdataIn.cov.Nq(theSub)); %inefficient but whatever.  Not weighted as assumed equally good estimates of noise, except for sample size differences.
                            end;
                            EPdataOut.analysis.badChans(numSubs+1,iCell,iChan)= EPdataOut.analysis.badChans(numSubs+1,iCell,iChan)+EPdataIn.analysis.badChans(theSub,iCell,iChan);
                        end;
                    else %none of the subjects have good data for this channel
                        if strcmp(EPdataIn.dataType,'average')
                            EPdataOut.analysis.badChans(numSubs+1,iCell,iChan)=NaN;
                        else
                            EPdataOut.analysis.badChans(numSubs+1,iCell,iChan)=-1;
                        end;
                    end;
                end;
            end;
        end;
        EPdataOut.data(:,:,:,end+1,:,:,:)=newCellData;
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData(:,:,:,end+1,:,:,:)=newFacData;
        end;
        EPdataOut.noise(:,:,:,end+1,:,:)=newCellNoise;
        EPdataOut.std(:,:,:,end+1,:,:)=newCellStd;
        %generate stdCM
        for iChan=1:numChans
            if isempty(EPdataIn.facVecS)
                theoldChan=iChan;
            else
                theoldChan=1;
            end;
            for iSample = 1:numPoints
                if isempty(EPdataIn.facVecT)
                    theoldSamp=iSample;
                else
                    theoldSamp=1;
                end;
                for iFactor = 1:numFacs
                    for iFreq = 1:max(numFreqs,1)
                        if isempty(EPdataIn.facVecF)
                            theoldFreq=iFreq;
                        else
                            theoldFreq=1;
                        end;
                        for iRel = 1:max(numRels,1)
                            grandMean=0;
                            for iSub=1:numSubs
                                for iCell=1:numCells
                                    if (EPdataIn.avgNum(iSub,iCell)>-1)
                                        if iFactor > numSGLfacs
                                            grandMean=grandMean+EPdataIn.facData(iChan,iSample,iCell,iSub,iFactor-numSGLfacs,iFreq,iRel);
                                        else
                                            grandMean=grandMean+EPdataIn.data(theoldChan,theoldSamp,iCell,iSub,iFactor,theoldFreq,iRel);
                                        end;
                                    end;
                                end;
                            end;
                            CMdata=nan(numSubs,numCells);
                            for iSub=1:numSubs
                                goodCells=find(EPdataIn.avgNum(iSub,:)>-1);
                                if iFactor > numSGLfacs
                                    CMdata(iSub,goodCells)=squeeze((EPdataIn.facData(iChan,iSample,goodCells,iSub,iFactor-numSGLfacs,iFreq,iRel))-mean(squeeze(EPdataIn.facData(iChan,iSample,goodCells,iSub,iFactor-numSGLfacs,iFreq,iRel)))+grandMean);
                                else
                                    CMdata(iSub,goodCells)=squeeze((EPdataIn.data(theoldChan,theoldSamp,goodCells,iSub,iFactor,theoldFreq,iRel))-mean(squeeze(EPdataIn.data(theoldChan,theoldSamp,goodCells,iSub,iFactor,theoldFreq,iRel)))+grandMean);
                                end;
                            end;
                            for iCell=1:numCells
                                CMcol=CMdata(:,iCell);
                                numCMcol=length(find(~isnan(CMcol)));
                                newCellStdCM(iChan,iSample,iCell,1,iFactor,iFreq,iRel)=std(CMcol(~isnan(CMcol)))*sqrt(numCMcol/(numCMcol-1));
                            end;
                        end;
                    end;
                end;
            end;
        end;
        EPdataOut.stdCM(:,:,:,end+1,:,:)=newCellStdCM;
        if ~isempty(EPdataIn.cov)
            %add together cov matrices weighted by their sample size, per MNE manual 2.7.3 p.90.
            Nq=sum(EPdataIn.cov.Nq(ismember(combineData,theOkaySubs)));
            EPdataOut.cov.covMatrix(end+1,:,:)=newCellCov/Nq;
            EPdataOut.cov.Nq(end+1)=Nq;
        end;
        if length(combineData) > 1
            EPdataOut.subTypes{end+1,1}='GAV';
        else
            EPdataOut.subTypes{end+1,1}=EPdataIn.subTypes{combineData};
        end;
        
        numSubs=numSubs+1;
        EPdataOut.subNames{end+1,1}=combineName;
        
        if ~isempty(EPdataIn.trialSpecs)
            %generate the averaged trial specs
            for iSpec=1:length(EPdataIn.trialSpecNames)
                for iCell=1:length(EPdataIn.cellNames)
                    if ischar(EPdataIn.trialSpecs{iCell,iSpec,combineData(1)}) && all(strcmp(EPdataIn.trialSpecs{iCell,iSpec,combineData(1)},EPdataIn.trialSpecs(iCell,iSpec,combineData)))
                        EPdataOut.trialSpecs{iCell,iSpec,numSubs}=EPdataIn.trialSpecs{iCell,iSpec,1};
                    else
                        theSpecs=EPdataIn.trialSpecs(iCell,iSpec,combineData);
                        theSpecs=squeeze(theSpecs);
                        theSpecs=theSpecs(:);
                        theStrSpecs=find(cellfun(@ischar,theSpecs));
                        for iStrSpec=1:length(theStrSpecs)
                            theStrSpec=theStrSpecs(iStrSpec);
                            theSpecs{theStrSpec}=str2num(theSpecs{theStrSpec});
                            if ~isnumeric(theSpecs{theStrSpec})
                                theSpecs{theStrSpec}=[]; %dealing with weird Matlab bug where 'cab' was being converted into a driver handle
                            end;
                        end;
                        EPdataOut.trialSpecs{iCell,iSpec,numSubs}=mean(cell2mat(theSpecs));
                    end;
                end;
            end;
        else
            EPdataOut.trialSpecs=cell(numCells,0,numSubs);
        end;

        EPdataOut.events(end+1,:)=cell(1,size(EPdataIn.events,2));
        
        if ~isempty(EPdataIn.impedances.channels)
        	EPdataOut.impedances.channels(:,end+1)=NaN;
        end;
        if ~isempty(EPdataIn.impedances.ground)
        	EPdataOut.impedances.ground(end+1)=NaN;
        end;

    case 'factors'
        %totalWeight=sum(abs(combineWeights));
        newFacData=zeros(numChans,max(numPoints,1),numCells,numSubs,1,max(numFreqs,1),max(numRels,1));
        for i=find(combineWeights)'
            theWeight=combineWeights(i);
            newFacData=newFacData+(theWeight)*ep_expandFacs(EPdataIn,[],[],[],[],i,[]);
        end;
        
        if isempty(EPdataIn.facData)
            EPdataOut.facData=newFacData;
        else
            EPdataOut.facData(:,:,:,:,end+1,:,:)=newFacData;
        end;
        
        numFacs=numFacs+1;
        EPdataOut.facNames{end+1,1}=combineName;
        EPdataOut.facTypes{end+1,1}='CMB';
    otherwise
        disp('Data dimension not recognized.');
        EPdataOut=[];
        return
end;

[err]=ep_checkEPfile(EPdataOut);
if err
    EPdataOut=[];
end;
