function [EPdataOut]=ep_addData(EPdataIn,EPadd,dataDimension);
%  [EPdataOut]=ep_addData(EPdataIn,EPadd,dataDimension);
%       Adds a level of data to one of the data dimensions, using a partial EPdata variable as the source.
%       Defaults added where not present in the partial EPdata variable.
%       The corresponding names field is required.
%
%Inputs:
%  EPdataIn       : Structured array with the input data and accompanying information in EP file format.  See readData.
%  EPadd          : Structured array with the information to be added.
%  dataDimension  : Which data dimension to add to (i.e., 'channels', 'cells', 'subjects', 'factors')
%
%Outputs:
%  EPdataOut      : Structured array with the output data and accompanying information in EP file format.  See readData.

%History:
%  by Joseph Dien (6/4/14)
%  jdien07@mac.com
%
% modified 10/2/14 JD
% Can add datasets with different trial specs.
%
% bugfix 12/2/14 JD
% If adding two sets of data with cells with the same trialNames, then make
% sure they no longer overlap.
%
% bugfix 5/19/15 JD
% Fixed crash when adding non single-trial datasets.
%
% modified 9/4/15 JD
% Added trial specs for average files.
%
% modified 10/12/15 JD
% Fixed crash when adding an average file with trial specs.
%
% bugfix 10/26/15 JD
% Fixed crash when adding to a data file with an empty trialspecs field.
%
% modified 10/16/16 JD
% Added .stims field.
%
% bugfix 6/1/17 JD
% Fixed trialnames not being correctly numbered when adding cells to a single-trial file.
%
% modified 12/22/17 JD
% Added support for impedances field.
%
% bugfix 12/6/17 JD
% Fixed crash when input file has impedances field and channels are being added but none are being added to the impedances field.
%
% modified & bugfix 2/11/18 JD
% Added support for stdCM field.
% No longer tries to combine std values.  Instead sets to zero.
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

switch dataDimension
    case 'channels'
        if ~isfield(EPadd,'chanNames')
            disp('Error: No channel names specified.');
            return
        end;
        numAdded=length(EPadd.chanNames);
        if isempty(EPdataIn.facVecS)
            if isfield(EPadd,'data')
                EPdataOut.data(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.data;
            else
                EPdataOut.data(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.noise)
            if isfield(EPadd,'noise')
                EPdataOut.noise(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.noise;
            else
                EPdataOut.noise(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.std)
            if isfield(EPadd,'std')
                EPdataOut.std(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.std;
            else
                EPdataOut.std(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.stdCM)
            if isfield(EPadd,'stdCM')
                EPdataOut.stdCM(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.stdCM;
            else
                EPdataOut.stdCM(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.facVecS)
            if isfield(EPadd,'facVecS')
                EPdataOut.facVecS(end+1:end+numAdded,:)=EPadd.facVecS;
            else
                EPdataOut.facVecS(end+1:end+numAdded,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData')
                EPdataOut.facData(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.facData;
            else
                EPdataOut.facData(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badChans')
                EPdataOut.analysis.badChans(:,:,end+1:end+numAdded)=EPadd.analysis.badChans;
            else
                EPdataOut.analysis.badChans(:,:,end+1:end+numAdded)=0;
            end;
        else
            EPdataOut.analysis.badChans(:,:,end+1:end+numAdded)=0;
        end;
        if ~isempty(EPdataIn.cov)
            if isfield(EPadd,'cov')
                EPdataOut.cov.covMatrix(:,end+1:end+numAdded,end+1:end+numAdded)=EPadd.cov.covMatrix;
            else
                EPdataOut.cov.covMatrix(:,end+1:end+numAdded,end+1:end+numAdded)=NaN;
            end;
        end;
        EPdataOut.chanNames(end+1:end+numAdded,1)=EPadd.chanNames;
        if isfield(EPadd,'chanTypes')
            EPdataOut.chanTypes(end+1:end+numAdded,1)=EPadd.chanTypes;
        else
            [EPdataOut.chanTypes{end+1:end+numAdded,1}]=deal('SGL');
        end;
        
        if ~isempty(EPdataIn.eloc)
            if isfield(EPadd,'eloc')
                EPdataOut.eloc(end+1:end+numAdded)=EPadd.eloc;
            else
                EPdataOut.eloc(end+1:end+numAdded)=cell2struct(cell(numAdded,length(fieldnames(EPdataIn.eloc))),fieldnames(EPdataIn.eloc),2);
            end;
        end;
        
        if ~isempty(EPdataIn.relNames)
            EPdataOut.relNames{end+1:end+numAdded,1}=EPadd.relNames;
            if ~isfield(EPadd,'data')
                EPdataOut.data(1:end-1,:,:,:,:,:,end+1)=0;
                EPdataOut.data(end,:,:,:,:,:,end)=1; %coherence of channel with itself is real number one.
            end;
        end;
        
        if ~isempty(EPdataIn.impedances.channels)
            if isfield(EPadd,'impedances') && isfield(EPadd.impedances,'channels') && ~isempty(EPadd.impedances.channels)
                EPdataOut.impedances.channels(end+1:end+numAdded,:)=EPadd.impedances.channels;
            else
                EPdataOut.impedances.channels(end+1:end+numAdded,:)=NaN;
            end;
        end;
        
    case 'cells'
        
        if ~isfield(EPadd,'cellNames')
            disp('Error: No cell names specified.');
            return
        end;
        if ~isempty(EPdataIn.trialNames)
            if ~isfield(EPadd,'trialNames')
                disp('Error: No trial names specified.');
                return
            end;
        end;
        
        numAdded=length(EPadd.cellNames);
        if ~isempty(EPdataIn.avgNum)
            if isfield(EPadd,'avgNum')
                EPdataOut.avgNum(:,end+1:end+numAdded)=EPadd.avgNum;
            else
                EPdataOut.avgNum(:,end+1:end+numAdded)=-1;
            end;
        end;
        if ~isempty(EPdataIn.covNum)
            if isfield(EPadd,'covNum')
                EPdataOut.covNum(:,end+1:end+numAdded)=EPadd.covNum;
            else
                EPdataOut.covNum(:,end+1:end+numAdded)=-1;
            end;
        end;
        if ~isempty(EPdataIn.subNum)
            if isfield(EPadd,'subNum')
                EPdataOut.subNum(:,end+1:end+numAdded)=EPadd.subNum;
            else
                EPdataOut.subNum(:,end+1:end+numAdded)=-1;
            end;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'blinkTrial')
                EPdataOut.analysis.blinkTrial(:,end+1:end+numAdded)=EPadd.analysis.blinkTrial;
            else
                EPdataOut.analysis.blinkTrial(:,end+1:end+numAdded)=0;
            end;
        else
            EPdataOut.analysis.blinkTrial(:,end+1:end+numAdded)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeTrial')
                EPdataOut.analysis.saccadeTrial(:,end+1:end+numAdded)=EPadd.analysis.saccadeTrial;
            else
                EPdataOut.analysis.saccadeTrial(:,end+1:end+numAdded)=0;
            end;
        else
            EPdataOut.analysis.saccadeTrial(:,end+1:end+numAdded)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeOnset')
                EPdataOut.analysis.saccadeOnset(:,end+1:end+numAdded)=EPadd.analysis.saccadeOnset;
            else
                EPdataOut.analysis.saccadeOnset(:,end+1:end+numAdded)=0;
            end;
        else
            EPdataOut.analysis.saccadeOnset(:,end+1:end+numAdded)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'moveTrial')
                EPdataOut.analysis.moveTrial(:,end+1:end+numAdded)=EPadd.analysis.moveTrial;
            else
                EPdataOut.analysis.moveTrial(:,end+1:end+numAdded)=0;
            end;
        else
            EPdataOut.analysis.moveTrial(:,end+1:end+numAdded)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badTrials')
                EPdataOut.analysis.badTrials(:,end+1:end+numAdded)=EPadd.analysis.badTrials;
            else
                EPdataOut.analysis.badTrials(:,end+1:end+numAdded)=0;
            end;
        else
            EPdataOut.analysis.badTrials(:,end+1:end+numAdded)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badChans')
                EPdataOut.analysis.badChans(:,end+1:end+numAdded,:)=EPadd.analysis.badChans;
            else
                EPdataOut.analysis.badChans(:,end+1:end+numAdded,:)=0;
            end;
        else
            EPdataOut.analysis.badChans(:,end+1:end+numAdded,:)=0;
        end;
        if isfield(EPadd,'recTime')
            EPdataOut.recTime(end+1:end+numAdded)=EPadd.recTime;
        else
            EPdataOut.recTime(end+1:end+numAdded)=1;
        end;
        
        if isfield(EPadd,'data')
            EPdataOut.data(:,:,end+1:end+numAdded,:,:,:,:)=EPadd.data;
        else
            EPdataOut.data(:,:,end+1:end+numAdded,:,:,:,:)=0;
        end;
        
        if ~isempty(EPdataIn.noise)
            if isfield(EPadd,'noise')
                EPdataOut.noise(:,:,end+1:end+numAdded,:,:,:,:)=EPadd.noise;
            else
                EPdataOut.noise(:,:,end+1:end+numAdded,:,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.std)
            if isfield(EPadd,'std')
                EPdataOut.std(:,:,end+1:end+numAdded,:,:,:,:)=EPadd.std;
            else
                EPdataOut.std(:,:,end+1:end+numAdded,:,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.stdCM)
            if isfield(EPadd,'stdCM')
                EPdataOut.stdCM(:,:,end+1:end+numAdded,:,:,:,:)=EPadd.stdCM;
            else
                EPdataOut.stdCM(:,:,end+1:end+numAdded,:,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData')
                EPdataOut.facData(:,:,end+1:end+numAdded,:,:,:,:)=EPadd.facData;
            else
                EPdataOut.facData(:,:,end+1:end+numAdded,:,:,:,:)=0;
            end;
        end;
        EPdataOut.cellNames(end+1:end+numAdded,1)=EPadd.cellNames;
        if isfield(EPadd,'cellTypes')
            EPdataOut.cellTypes(end+1:end+numAdded,1)=EPadd.cellTypes;
        else
            [EPdataOut.cellTypes{end+1:end+numAdded,1}]=deal('SGL');
        end;
        if ~isempty(EPadd.trialNames)
            cellNameList=unique(EPadd.cellNames);
            newTrials=zeros(numAdded,1);
            for iCell=1:length(cellNameList)
                theCell=cellNameList{iCell};
                whichCells=find(strcmp(theCell,EPadd.cellNames));
                theAddTrials=EPadd.trialNames(whichCells);
                theInTrials=EPdataIn.trialNames(find(strcmp(theCell,EPdataIn.cellNames)));
                if ~isempty(intersect(theAddTrials,theInTrials))
                    theAddTrials=theAddTrials+max(theInTrials);
                end;
                newTrials(whichCells)=theAddTrials;
            end;
            EPdataOut.trialNames(end+1:end+numAdded,1)=newTrials;
        end;
        
        if isfield(EPadd,'trialSpecs') && ~isempty(EPadd.trialSpecs)
            if ~iscell(EPdataOut.trialSpecs) && isempty(EPdataOut.trialSpecs)
                EPdataOut.trialSpecs=cell(0);
            end;
            EPdataOut.trialSpecs(end+1:end+numAdded,:,:)=cell(numAdded,size(EPdataIn.trialSpecs,2),size(EPdataIn.trialSpecs,3));
            for iSpec=1:length(EPadd.trialSpecNames)
                oldSpec=find(strcmp(EPadd.trialSpecNames{iSpec},EPdataOut.trialSpecNames));
                if ~isempty(oldSpec)
                    EPdataOut.trialSpecs(end-numAdded+1:end,oldSpec,:)=EPadd.trialSpecs(:,iSpec,:);
                else
                    EPdataOut.trialSpecNames(end+1)=EPadd.trialSpecNames(iSpec);
                    EPdataOut.trialSpecs(end-numAdded:end,end+1,:)=EPadd.trialSpecs(:,iSpec,:);
                end;
            end;
        else
            if ~isempty(EPdataIn.trialSpecs)
                EPdataOut.trialSpecs(end+1:end+numAdded,:,:)=cell(numAdded,size(EPdataIn.trialSpecs,2),size(EPdataIn.trialSpecs,3));
            end;
        end;
        
        if ~isempty(EPdataIn.events)
            if isfield(EPadd,'events')
                EPdataOut.events(:,end+1:end+numAdded)=EPadd.events;
            else
                EPdataOut.events(:,end+1:end+numAdded)=cell(size(EPdataIn.events,1),numAdded);
            end;
            
            if ~isempty(EPdataIn.stims)
                for iStim=1:length(EPdataIn.stims)
                    if ~any(strcmp(EPdataIn.stims(iStim).name,{EPdataOut.stims.name}))
                        EPdataOut.stims(end+1)=EPdataIn.stims(iStim); %keep only stim images whose events are still in the segmented data.
                    end;
                end;
            end;
        end;
        
    case 'subjects'
        
        if ~isfield(EPadd,'subNames')
            disp('Error: No subject names specified.');
            return
        end;
        numAdded=length(EPadd.subNames);
        
        if ~isempty(EPdataIn.subjectSpecs)
            if isfield(EPadd,'subjectSpecs')
                EPdataOut.subjectSpecs(end+1:end+numAdded,:)=EPadd.subjectSpecs;
            else
                EPdataOut.subjectSpecs(end+1:end+numAdded,:)=cell(numAdded,size(EPdataIn.subjectSpecs,2));
            end;
        end;
        
        if ~isempty(EPdataIn.avgNum)
            if isfield(EPadd,'avgNum')
                EPdataOut.avgNum(end+1:end+numAdded,:)=EPadd.avgNum;
            else
                EPdataOut.avgNum(end+1:end+numAdded,:)=-1;
            end;
        end;
        if ~isempty(EPdataIn.covNum)
            if isfield(EPadd,'covNum')
                EPdataOut.covNum(end+1:end+numAdded,:)=EPadd.covNum;
            else
                EPdataOut.covNum(end+1:end+numAdded,:)=-1;
            end;
        end;
        if ~isempty(EPdataIn.subNum)
            if isfield(EPadd,'subNum')
                EPdataOut.subNum(end+1:end+numAdded,:)=EPadd.subNum;
            else
                EPdataOut.subNum(end+1:end+numAdded,:)=-1;
            end;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'blinkTrial')
                EPdataOut.analysis.blinkTrial(end+1:end+numAdded,:)=EPadd.analysis.blinkTrial;
            else
                EPdataOut.analysis.blinkTrial(end+1:end+numAdded,:)=0;
            end;
        else
            EPdataOut.analysis.blinkTrial(end+1:end+numAdded,:)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeTrial')
                EPdataOut.analysis.saccadeTrial(end+1:end+numAdded,:)=EPadd.analysis.saccadeTrial;
            else
                EPdataOut.analysis.saccadeTrial(end+1:end+numAdded,:)=0;
            end;
        else
            EPdataOut.analysis.saccadeTrial(end+1:end+numAdded,:)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeOnset')
                EPdataOut.analysis.saccadeOnset(end+1:end+numAdded,:)=EPadd.analysis.saccadeOnset;
            else
                EPdataOut.analysis.saccadeOnset(end+1:end+numAdded,:)=0;
            end;
        else
            EPdataOut.analysis.saccadeOnset(end+1:end+numAdded,:)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'moveTrial')
                EPdataOut.analysis.moveTrial(end+1:end+numAdded,:)=EPadd.analysis.moveTrial;
            else
                EPdataOut.analysis.moveTrial(end+1:end+numAdded,:)=0;
            end;
        else
            EPdataOut.analysis.moveTrial(end+1:end+numAdded,:)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badTrials')
                EPdataOut.analysis.badTrials(end+1:end+numAdded,:)=EPadd.analysis.badTrials;
            else
                EPdataOut.analysis.badTrials(end+1:end+numAdded,:)=0;
            end;
        else
            EPdataOut.analysis.badTrials(end+1:end+numAdded,:)=0;
        end;
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badChans')
                EPdataOut.analysis.badChans(end+1:end+numAdded,:,:)=EPadd.analysis.badChans;
            else
                EPdataOut.analysis.badChans(end+1:end+numAdded,:,:)=0;
            end;
        else
            EPdataOut.analysis.badChans(end+1:end+numAdded,:,:)=0;
        end;        
        
        if isfield(EPadd,'data')
            EPdataOut.data(:,:,:,end+1:end+numAdded,:,:,:)=EPadd.data;
        else
            EPdataOut.data(:,:,:,end+1:end+numAdded,:,:,:)=0;
        end;
        
        if ~isempty(EPdataIn.noise)
            if isfield(EPadd,'noise')
                EPdataOut.noise(:,:,:,end+1:end+numAdded,:,:,:)=EPadd.noise;
            else
                EPdataOut.noise(:,:,:,end+1:end+numAdded,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.std)
            if isfield(EPadd,'std')
                EPdataOut.std(:,:,:,end+1:end+numAdded,:,:,:)=EPadd.std;
            else
                EPdataOut.std(:,:,:,end+1:end+numAdded,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.stdCM)
            if isfield(EPadd,'stdCM')
                EPdataOut.stdCM(:,:,:,end+1:end+numAdded,:,:,:)=EPadd.stdCM;
            else
                EPdataOut.stdCM(:,:,:,end+1:end+numAdded,:,:,:)=0;
            end;
        end;
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData')
                EPdataOut.facData(:,:,:,end+1:end+numAdded,:,:,:)=EPadd.facData;
            else
                EPdataOut.facData(:,:,:,end+1:end+numAdded,:,:,:)=0;
            end;
        end;
        EPdataOut.subNames(end+1:end+numAdded,1)=EPadd.subNames;
        if isfield(EPadd,'subTypes')
            EPdataOut.subTypes(end+1:end+numAdded,1)=EPadd.subTypes;
        else
            if strcmp(EPdataIn,'average')
                [EPdataOut.subTypes{end+1:end+numAdded,1}]=deal('AVG');
            else
                [EPdataOut.subTypes{end+1:end+numAdded,1}]=deal('RAW');
            end;
        end;        
        if ~isempty(EPdataIn.cov)
            if isfield(EPadd,'cov')
                EPdataOut.cov.covMatrix(end+1:end+numAdded,:,:)=EPadd.cov.covMatrix;
                EPdataOut.cov.Nq(end+1:end+numAdded)=EPadd.cov.Nq;
            else
                EPdataOut.cov.covMatrix(end+1:end+numAdded,:,:)=NaN;
                EPdataOut.cov.Nq(end+1:end+numAdded)=NaN;
            end;
        end;        
        if ~isempty(EPdataIn.events)
            if isfield(EPadd,'events')
                EPdataOut.events(end+1:end+numAdded,:)=EPadd.events;
            else
                EPdataOut.events(end+1:end+numAdded,:)=cell(numAdded,size(EPdataIn.events,2));
            end;
            
            if ~isempty(EPdataIn.stims)
                for iStim=1:length(EPdataIn.stims)
                    if ~any(strcmp(EPdataIn.stims(iStim).name,{EPdataOut.stims.name}))
                        EPdataOut.stims(end+1)=EPdataIn.stims(iStim); %keep only stim images whose events are still in the segmented data.
                    end;
                end;
            end;
        end;        

        if isfield(EPadd,'trialSpecs') && ~isempty(EPadd.trialSpecs)
            if ~iscell(EPdataOut.trialSpecs) && isempty(EPdataOut.trialSpecs)
                EPdataOut.trialSpecs=cell(0);
            end;
            EPdataOut.trialSpecs(:,:,end+1:end+numAdded)=cell(size(EPdataIn.trialSpecs,1),size(EPdataIn.trialSpecs,2),numAdded);
            for iSpec=1:length(EPadd.trialSpecNames)
                oldSpec=find(strcmp(EPadd.trialSpecNames{iSpec},EPdataOut.trialSpecNames));
                if ~isempty(oldSpec)
                    EPdataOut.trialSpecs(:,oldSpec,end-numAdded+1:end)=EPadd.trialSpecs(:,iSpec,:);
                else
                    EPdataOut.trialSpecNames(end+1)=EPadd.trialSpecNames(iSpec);
                    EPdataOut.trialSpecs(:,end+1,end-numAdded+1:end)=EPadd.trialSpecs(:,iSpec,:);
                end;
            end;
        else
            if ~isempty(EPdataIn.trialSpecs)
                EPdataOut.trialSpecs(:,:,end+1:end+numAdded)=cell(size(EPdataIn.trialSpecs,1),size(EPdataIn.trialSpecs,2),numAdded);
            end;
        end;
        
        if ~isempty(EPdataIn.impedances.channels)
            if isfield(EPadd,'impedances') && isfield(EPadd.impedances,'channels') && ~isempty(EPadd.impedances.channels)
                EPdataOut.impedances.channels(:,end+1:end+numAdded)=EPadd.impedances.channels;
            else
                EPdataOut.impedances.channels(:,end+1:end+numAdded)=NaN;
            end;
        end;
        
        if ~isempty(EPdataIn.impedances.ground)
            if isfield(EPadd,'impedances') && isfield(EPadd.impedances,'ground') && ~isempty(EPadd.impedances.ground)
                EPdataOut.impedances.ground(end+1:end+numAdded)=EPadd.impedances.ground;
            else
                EPdataOut.impedances.ground(end+1:end+numAdded)=NaN;
            end;
        end;
        
    case 'factors'
        
        if ~isfield(EPadd,'facNames')
            disp('Error: No factor names specified.');
            return
        end;
        numAdded=length(EPadd.facNames);
        
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData')
                EPdataOut.facData(:,:,:,:,end+1:end+numAdded,:,:)=EPadd.facData;
            else
                EPdataOut.facData(:,:,:,:,end+1:end+numAdded,:,:)=0;
            end;
        end;        
        
        EPdataOut.facNames(end+1:end+numAdded,1)=EPadd.facNames;
        if isfield(EPadd,'facTypes')
            EPdataOut.facTypes(end+1:end+numAdded,1)=EPadd.facTypes;
        else
            [EPdataOut.facTypes{end+1:end+numAdded,1}]=deal('SGL');
        end;         

    otherwise
        disp('Data dimension not recognized.');
        EPdataOut=[];
        return
end;

[err]=ep_checkEPfile(EPdataOut);
if err
    EPdataOut=[];
end;
