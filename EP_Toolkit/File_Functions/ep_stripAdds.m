function [EPdata]=ep_stripAdds(EPdata,types);
%  [err]=ep_checkEPfile(data,types);
%       Strips out adds (combined channels, cells, subjects, and factors) to EP data.
%       If no types are specified, all chan, cell, and factor adds are stripped out and lowest subject level is kept.
%
%Inputs:
%  EPdata         : Structured array with the data and accompanying information in EP file format.  See readData.
%  types          : Optional cell array of types of information to keep:
%                   SGLchan (single channels, including EEG, MGM, MGA, MGP, ANS, ECG, PPL, XEY, YEY, and REF)
%                   REGchan (regional channels)
%                   SGLcell (single cells)
%                   CMBcell (combined cells)
%                   STScell (sample test statistics output)
%                   RAW (single trial data)
%                   AVG (averaged data)
%                   GAV (grand average data)
%                   SGLfac (single factors)
%                   CMBfac (combined factors)
%
%Outputs:
%  EPdata        : Structured array with the data and accompanying information in EP file format.  See readData.

%History:
%  by Joseph Dien (7/27/09)
%  jdien07@mac.com
%
% modified 2/11/10 JD
% analysis fields no longer optional.
%
%  bugfix 2/15/10 JD
%  No longer crashes when subject adds are being stripped from data with no subject specs.
%
% bugfix 4/24/10 JD
% Fixed CMB factors not dropped from data when present in .data rather than .facData.
% 
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
%
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% modified 9/22/13 JD
% Added support for ECG channel type.
%
% modified 10/9/13 JD
% Added recTime field.
%
% bugfix 1/12/14 JD
% Fixed crash when keeping CMB factors.
%
% bugfix 2/27/14 JD
% Fixed crash when the file has different types of adds in a category (e.g., GAV and AVG in subjects)
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% modified 3/24/14 JD
% Added .cov field.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% modified 5/28/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 6/9/14 JD
% Added STS cell type and PPL, XEY, and YEY channel types.
%
% modified 9/4/15 JD
% Added trial specs for average files.
%
% modified 10/16/16 JD
% Added .stims field.
%
% bugfix 2/27/17 JD
% Fixed crash when averaging data containing impedance values.
%
% modified 2/11/18 JD
% Added support for stdCM field.
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

if (nargin < 1) || isempty(EPdata)
    msg{1}='No data provided.';
    [msg]=ep_errorMsg(msg);
    return
end;

if nargin < 2
    types{1}='SGLchan';
end;

if isempty(find(ismember(types,{'SGLchan','REGchan'}), 1))
    types{end+1}='SGLchan'; %at least some type of channel needs to be kept
end;

if isempty(find(ismember(types,{'SGLcell','CMBcell'}), 1))
    types{end+1}='SGLcell'; %at least some type of cell needs to be kept
end;

if isempty(find(ismember(types,{'RAW','AVG','GAV'}), 1))
    if ~isempty(find(strcmp('RAW',EPdata.subTypes), 1))
        types{end+1}='RAW';
    elseif ~isempty(find(strcmp('AVG',EPdata.subTypes), 1))
        types{end+1}='AVG';
    else
        types{end+1}='GAV';
    end;
end;

if isempty(find(ismember(types,{'SGLfac','CMBfac'}), 1))
    types{end+1}='SGLfac'; %at least some type of factor needs to be kept
end;

chanAdds=[];
if isempty(find(strcmp('REGchan',types), 1)) %if regional channels not listed as type to keep, then add to list to strip out.
    chanAdds=[chanAdds; find(strcmp('REG',EPdata.chanTypes))];
end;
if isempty(find(strcmp('SGLchan',types), 1)) %if single channels not listed as type to keep, then add to list to strip out.
    chanAdds=[chanAdds; find(ismember(EPdata.chanTypes,{'EEG','MGM','MGA','MGP','ANS','ECG','PPL','XEY','YEY'}))];
end;
if ~isempty(chanAdds)
    keepChans=setdiff([1:length(EPdata.chanNames)],chanAdds);
    if ~isempty(EPdata.relNames)
        keepRelChans=keepChans;
    else
        keepRelChans=1;
    end;
    if ~isempty(EPdata.facVecS)
        EPdata.facVecS=EPdata.facVecS(keepChans,:);
    else
        EPdata.data=EPdata.data(keepChans,:,:,:,:,:,keepRelChans);
    end;
    if ~isempty(EPdata.facData)
        EPdata.facData=EPdata.facData(keepChans,:,:,:,:,:,keepRelChans);
    end;
    if ~isempty(EPdata.noise)
        EPdata.noise=EPdata.noise(keepChans,:,:,:,:,:);
    end;
    if ~isempty(EPdata.std)
        EPdata.std=EPdata.std(keepChans,:,:,:,:,:);
    end;
    if ~isempty(EPdata.stdCM)
        EPdata.stdCM=EPdata.stdCM(keepChans,:,:,:,:,:);
    end;
    if ~isempty(EPdata.cov)
        EPdata.cov.covMatrix=EPdata.cov.covMatrix(:,keepChans,keepChans);
    end;

    EPdata.chanNames=EPdata.chanNames(keepChans);
    EPdata.chanTypes=EPdata.chanTypes(keepChans);
    if ~isempty(EPdata.eloc)
        EPdata.eloc=EPdata.eloc(keepChans);
    end;

    EPdata.analysis.badChans=EPdata.analysis.badChans(:,:,keepChans);
    if ~isempty(EPdata.impedances.channels)
        EPdata.impedances.channels=EPdata.impedances.channels(keepChans,:);
    end;
end;

cellAdds=[];
if isempty(find(strcmp('CMBcell',types), 1))
    cellAdds=[cellAdds; find(strcmp('CMB',EPdata.cellTypes))];
end;
if isempty(find(strcmp('SGLcell',types), 1))
    cellAdds=[cellAdds; find(strcmp('SGL',EPdata.cellTypes))];
end;
if isempty(find(strcmp('STScell',types), 1))
    cellAdds=[cellAdds; find(strcmp('STS',EPdata.cellTypes))];
end;
if ~isempty(cellAdds)
    keepCells=setdiff([1:length(EPdata.cellNames)],cellAdds);

    EPdata.avgNum=EPdata.avgNum(:,keepCells);
    EPdata.covNum=EPdata.covNum(:,keepCells);
    EPdata.subNum=EPdata.subNum(:,keepCells);
    EPdata.data=EPdata.data(:,:,keepCells,:,:,:,:);
    if ~isempty(EPdata.facData)
        EPdata.facData=EPdata.facData(:,:,keepCells,:,:,:,:);
    end;
    if ~isempty(EPdata.noise)
        EPdata.noise=EPdata.noise(:,:,keepCells,:,:,:);
    end;
    if ~isempty(EPdata.std)
        EPdata.std=EPdata.std(:,:,keepCells,:,:,:);
    end;
    if ~isempty(EPdata.stdCM)
        EPdata.stdCM=EPdata.stdCM(:,:,keepCells,:,:,:);
    end;
    EPdata.cellNames=EPdata.cellNames(keepCells);
    EPdata.cellTypes=EPdata.cellTypes(keepCells);
    EPdata.events=EPdata.events(:,keepCells);
    if ~isempty(EPdata.trialSpecNames)
        EPdata.trialSpecs=EPdata.trialSpecs(keepCells,:,:);
    end;
    if strcmp(EPdata.dataType,'single_trial')
        EPdata.trialNames=EPdata.trialNames(keepCells);
    end;
    
    EPdata.analysis.blinkTrial= EPdata.analysis.blinkTrial(:,keepCells);
    EPdata.analysis.saccadeTrial= EPdata.analysis.saccadeTrial(:,keepCells);
    EPdata.analysis.saccadeOnset= EPdata.analysis.saccadeOnset(:,keepCells);
    EPdata.analysis.moveTrial= EPdata.analysis.moveTrial(:,keepCells);
    EPdata.analysis.badTrials= EPdata.analysis.badTrials(:,keepCells);
    EPdata.analysis.badChans= EPdata.analysis.badChans(:,keepCells,:);
    
    EPdata.recTime= EPdata.recTime(keepCells);
end;

subAdds=[];
if isempty(find(strcmp('GAV',types), 1))
    subAdds=[subAdds; find(strcmp('GAV',EPdata.subTypes))];
end;
if isempty(find(strcmp('AVG',types), 1))
    subAdds=[subAdds; find(strcmp('AVG',EPdata.subTypes))];
end;
if isempty(find(strcmp('RAW',types), 1))
    subAdds=[subAdds; find(strcmp('RAW',EPdata.subTypes))];
end;
if ~isempty(subAdds)
    keepSubjects=setdiff([1:length(EPdata.subNames)],subAdds);

    EPdata.avgNum=EPdata.avgNum(keepSubjects,:);
    EPdata.covNum=EPdata.covNum(keepSubjects,:);
    EPdata.subNum=EPdata.subNum(keepSubjects,:);
    EPdata.data=EPdata.data(:,:,:,keepSubjects,:,:,:);
    if ~isempty(EPdata.facData)
        EPdata.facData=EPdata.facData(:,:,:,keepSubjects,:,:,:);
    end;
    if ~isempty(EPdata.noise)
        EPdata.noise=EPdata.noise(:,:,:,keepSubjects,:,:);
    end;
    if ~isempty(EPdata.std)
        EPdata.std=EPdata.std(:,:,:,keepSubjects,:,:);
    end;
    if ~isempty(EPdata.stdCM)
        EPdata.stdCM=EPdata.stdCM(:,:,:,keepSubjects,:,:);
    end;
    if ~isempty(EPdata.cov)
        EPdata.cov.covMatrix=EPdata.cov.covMatrix(keepSubjects,:,:);
        EPdata.cov.Nq=EPdata.cov.Nq(keepSubjects);
    end;
    EPdata.subNames=EPdata.subNames(keepSubjects);
    EPdata.subTypes=EPdata.subTypes(keepSubjects);
    if ~isempty(EPdata.subjectSpecs)
        EPdata.subjectSpecs=EPdata.subjectSpecs(keepSubjects,:);
    end;
    EPdata.events=EPdata.events(keepSubjects,:);
    if ~isempty(EPdata.trialSpecNames)
        EPdata.trialSpecs=EPdata.trialSpecs(:,:,keepSubjects);
    end;
    
    EPdata.analysis.blinkTrial= EPdata.analysis.blinkTrial(keepSubjects,:);
    EPdata.analysis.saccadeTrial= EPdata.analysis.saccadeTrial(keepSubjects,:);
    EPdata.analysis.saccadeOnset= EPdata.analysis.saccadeOnset(keepSubjects,:);
    EPdata.analysis.moveTrial= EPdata.analysis.moveTrial(keepSubjects,:);
    EPdata.analysis.badTrials= EPdata.analysis.badTrials(keepSubjects,:);
    EPdata.analysis.badChans= EPdata.analysis.badChans(keepSubjects,:,:);

    if ~isempty(EPdata.impedances.channels)
        EPdata.impedances.channels=EPdata.impedances.channels(:,keepSubjects);
    end;
    if ~isempty(EPdata.impedances.ground)
        EPdata.impedances.ground=EPdata.impedances.ground(keepSubjects);
    end;
end;

facAdds=[]; %facs to drop
CMBfacAdds=[]; %facData facs to drop
if isempty(find(strcmp('CMBfac',types)))
    CMBfacAdds=find(strcmp('CMB',EPdata.facTypes))-size(EPdata.data,5);
    facAdds=[facAdds find(strcmp('CMB',EPdata.facTypes))];
end
SGLfacAdds=[];
if isempty(find(strcmp('SGLfac',types), 1))
    SGLfacAdds=find(strcmp('SGL',EPdata.facTypes));
    facAdds=[facAdds SGLfacAdds];
end;
if ~isempty(facAdds) && ~isempty(EPdata.facNames) %don't strip out factors unless they are present in the first place
    keepFacs=setdiff([1:length(EPdata.facNames)],facAdds);
    SGLkeepFacs=setdiff([1:size(EPdata.data,5)],SGLfacAdds);
    
    EPdata.data=EPdata.data(:,:,:,:,SGLkeepFacs,:,:);
    if ~isempty(EPdata.facVecS)
        EPdata.facVecS=EPdata.facVecS(:,SGLkeepFacs);
    end;
    if ~isempty(EPdata.facVecT)
        EPdata.facVecT=EPdata.facVecT(:,SGLkeepFacs);
    end;
    if ~isempty(EPdata.facData)
        CMBkeepFacs=setdiff([1:size(EPdata.facData,5)],CMBfacAdds);
        EPdata.facData=EPdata.facData(:,:,:,:,CMBkeepFacs,:,:);
    end;
    if ~isempty(EPdata.noise)
        EPdata.noise=EPdata.noise(:,:,:,:,keepFacs,:);
    end;
    if ~isempty(EPdata.std)
        EPdata.std=EPdata.std(:,:,:,:,keepFacs,:);
    end;
    if ~isempty(EPdata.stdCM)
        EPdata.stdCM=EPdata.stdCM(:,:,:,:,keepFacs,:);
    end;

    EPdata.facNames=EPdata.facNames(keepFacs);
    EPdata.facTypes=EPdata.facTypes(keepFacs);
end;

if ~isempty(subAdds) || ~isempty(cellAdds)
    stims=EPdata.stims;
    EPdata.stims=struct('name',{},'image',{});
    for iStim=1:length(stims)
        for iSub=1:length(EPdata.subNames)
            for iCell=1:length(EPdata.cellNames)
                for iEvent=1:length(EPdata.events{iSub,iCell})
                    if any(strcmp(stims(iStim).name,{EPdata.events{iSub,iCell}(iEvent).keys.data}))
                        EPdata.stims(end+1)=stims(iStim); %keep only stim images whose events are still in the data.
                    end;
                end;
            end;
        end;
    end;
end;

[err]=ep_checkEPfile(EPdata);

