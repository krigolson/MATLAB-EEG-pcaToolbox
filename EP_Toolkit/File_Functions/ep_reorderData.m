function [EPdataOut]=ep_reorderData(EPdataIn,dataDimension,sortOrder);
%  [EPdataOut]=ep_reorderData(EPdataIn,dataDimension,sortSpec);
%       Reorders data according to specified order and dimension. 
%
%Inputs:
%  EPdataIn       : Structured array with the input data and accompanying information in EP file format.  See readData.
%  dataDimension  : The data dimension ('cells', 'subjects', 'channels', 'factors').
%  sortOrder       : The new order of the selected data dimension.
%
%Outputs:
%  EPdataOut      : Structured array with the output data and accompanying information in EP file format.  See readData.

%History:
%  by Joseph Dien (6/12/18)
%  jdien07@mac.com
%
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

switch dataDimension
    case 'cells'
        EPdataOut.trialSpecs=EPdataIn.trialSpecs(sortOrder,:,:);
        EPdataOut.avgNum=EPdataIn.avgNum(:,sortOrder);
        EPdataOut.covNum=EPdataIn.covNum(:,sortOrder);
        EPdataOut.subNum=EPdataIn.subNum(:,sortOrder);
        EPdataOut.analysis.blinkTrial=EPdataIn.analysis.blinkTrial(:,sortOrder);
        EPdataOut.analysis.saccadeTrial=EPdataIn.analysis.saccadeTrial(:,sortOrder);
        EPdataOut.analysis.saccadeOnset=EPdataIn.analysis.saccadeOnset(:,sortOrder);
        EPdataOut.analysis.moveTrial=EPdataIn.analysis.moveTrial(:,sortOrder);
        EPdataOut.analysis.badTrials=EPdataIn.analysis.badTrials(:,sortOrder);
        EPdataOut.analysis.badChans=EPdataIn.analysis.badChans(:,sortOrder,:);
        EPdataOut.recTime=EPdataIn.recTime(sortOrder);
        EPdataOut.data=EPdataIn.data(:,:,sortOrder,:,:,:,:);
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(:,:,sortOrder,:,:,:,:);
        end;
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise=EPdataIn.noise(:,:,sortOrder,:,:,:,:);
        end;
        if ~isempty(EPdataIn.std)
            EPdataOut.std=EPdataIn.std(:,:,sortOrder,:,:,:,:);
        end;
        if ~isempty(EPdataIn.stdCM)
            EPdataOut.stdCM=EPdataIn.stdCM(:,:,sortOrder,:,:,:,:);
        end;
        EPdataOut.cellNames=EPdataIn.cellNames(sortOrder);
        EPdataOut.cellTypes=EPdataIn.cellTypes(sortOrder); 
        if ~isempty(EPdataIn.trialNames)
            EPdataOut.trialNames=EPdataIn.trialNames(sortOrder);
        end;
        EPdataOut.events=EPdataIn.events(:,sortOrder);  

    case 'subjects'
        
        EPdataOut.subjectSpecs=EPdataIn.subjectSpecs(sortOrder,:);
        EPdataOut.avgNum=EPdataIn.avgNum(sortOrder,:);
        EPdataOut.covNum=EPdataIn.covNum(sortOrder,:);
        EPdataOut.subNum=EPdataIn.subNum(sortOrder,:);
        EPdataOut.analysis.blinkTrial=EPdataIn.analysis.blinkTrial(sortOrder,:);
        EPdataOut.analysis.saccadeTrial=EPdataIn.analysis.saccadeTrial(sortOrder,:);
        EPdataOut.analysis.saccadeOnset=EPdataIn.analysis.saccadeOnset(sortOrder,:);
        EPdataOut.analysis.moveTrial=EPdataIn.analysis.moveTrial(sortOrder,:);
        EPdataOut.analysis.badTrials=EPdataIn.analysis.badTrials(sortOrder,:);
        EPdataOut.analysis.badChans=EPdataIn.analysis.badChans(sortOrder,:,:);
        EPdataOut.data=EPdataIn.data(:,:,:,sortOrder,:,:,:);
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(:,:,:,sortOrder,:,:,:);
        end;
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise=EPdataIn.noise(:,:,:,sortOrder,:,:,:);
        end;
        if ~isempty(EPdataIn.std)
            EPdataOut.std=EPdataIn.std(:,:,:,sortOrder,:,:,:);
        end;
        if ~isempty(EPdataIn.stdCM)
            EPdataOut.stdCM=EPdataIn.stdCM(:,:,:,sortOrder,:,:,:);
        end;
        EPdataOut.subNames=EPdataIn.subNames(sortOrder);
        EPdataOut.subTypes=EPdataIn.subTypes(sortOrder);
        EPdataOut.events=EPdataIn.events(sortOrder,:);
        if ~isempty(EPdataIn.impedances.channels)
            EPdataOut.impedances.channels=EPdataIn.impedances.channels(:,sortOrder);
        end;
        if ~isempty(EPdataIn.impedances.ground)
            EPdataOut.impedances.ground=EPdataIn.impedances.ground(sortOrder);
        end;
        
    case 'channels'
        
        if ~isempty(EPdataIn.facVecS)
            EPdataOut.facVecS=EPdataIn.facVecS(sortOrder,:);
        else
            EPdataOut.data=EPdataIn.data(sortOrder,:,:,:,:,:,:);
        end;
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise=EPdataIn.noise(sortOrder,:,:,:,:);
        end;
        if ~isempty(EPdataIn.cov)
            EPdataOut.cov.covMatrix=EPdataIn.cov.covMatrix(:,sortOrder,sortOrder);
        end;
        if ~isempty(EPdataIn.std)
            EPdataOut.std=EPdataIn.std(sortOrder,:,:,:,:,:);
        end;
        if ~isempty(EPdataIn.stdCM)
            EPdataOut.stdCM=EPdataIn.stdCM(sortOrder,:,:,:,:,:);
        end;
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(sortOrder,:,:,:,:,:,:);
        end
        
        reference.type=EPdataIn.reference.type;
        reference.original=[];
        reference.current=[];
        
        for i=1:length(EPdataIn.reference.original)
            reference.original(i)=find(EPdataIn.reference.original(i) == sortOrder);
        end;
        for i=1:length(EPdataIn.reference.current)
            reference.current(i)=find(EPdataIn.reference.current(i) == sortOrder);
        end;
        EPdataOut.reference=reference;
        
        EPdataOut.chanNames=EPdataIn.chanNames(sortOrder);
        EPdataOut.chanTypes=EPdataIn.chanTypes(sortOrder);
        EPdataOut.eloc=EPdataIn.eloc(sortOrder);
        EPdataOut.analysis.badChans=EPdataIn.analysis.badChans(:,:,sortOrder);
        
    case 'factors'
        EPdataOut.facNames=EPdataIn.facNames(sortOrder);
        EPdataOut.facTypes=EPdataIn.facTypes(sortOrder);
        
        singleFacs=find(ismember(EPdataIn.facTypes,'SGL')); %combined factors always last and so not reordered.
        combFacs=find(ismember(EPdataIn.facTypes,'CMB'))-length(singleFacs);
        EPdataOut.data=EPdataIn.data(:,:,:,:,sortOrder(singleFacs),:,:);
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(:,:,:,:,sortOrder(combFacs),:,:);
        end;
        if ~isempty(EPdataIn.facVecT)
            EPdataOut.facVecT=EPdataIn.facVecT(:,sortOrder(singleFacs));
        end;
        if ~isempty(EPdataIn.facVecS)
            EPdataOut.facVecS=EPdataIn.facVecS(:,sortOrder(singleFacs));
        end;
        if ~isempty(EPdataIn.facVar)
            EPdataOut.facVar=EPdataIn.facVar(:,sortOrder(singleFacs));
        end;
        if ~isempty(EPdataIn.facVarQ)
            EPdataOut.facVarQ=EPdataIn.facVarQ(:,sortOrder(singleFacs));
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
