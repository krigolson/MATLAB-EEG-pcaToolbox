function EPdata=ep_loadEPdataset(EPdataset,theDataset);
%  EPdata=ep_loadEPdataset(EPdataset,theDataset);
%       Loads a dataset from the work directory.
%
%Inputs:
%  EPdataset      : Structured array with the list of files in the work directory
%     .EPwork     : The path of the work directory.
%     .dataName   : The list of dataset names.
%  theDataset     : Which dataset to load in.
%
%Outputs:
%  EPdata         : Structured array with the data and accompanying information.  See readData.

%History:
%  by Joseph Dien (7/11/09)
%  jdien07@mac.com
%
%  bugfix 8/3/09 JD
%  File names with multiple periods would confuse Matlab about the file type,
%  resulting in it trying to load the file as an ASCII file rather than as a .mat file.
%
%  modified 2/24/10 JD
%  Backward compatibility for changes in EPdata format.
%
%  modified 10/12/10 JD
%  Backward compatibility for shift to epochwise format for analysis fields for continuous data files.
% 
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
%
% modified 4/19/11 JD
% Added support for transforms field
%
% modified 2/5/12 JD
% Eliminated transforms field and added instead freqNames and facVecF fields
%
% modified 2/22/12 JD
% Backward compatibility for noise and std no longer optional.
%
% modified 1/11/13 JD
% Added option to do internal calculations of frequency data in either amplitude or power form.
%
% modified 10/13/13 JD
% Added recTime field.
%
% bugfix 10/21/13 JD
% Ensures that power comes after analysis field so order of fields is always the same.
%
% modified 10/29/13 JD
% Added .keys field to events.
%
% bugfix 1/7/14 JD
% Added check for .recTime field in .pca field.
%
% modified 2/26/14 JD
% Made .pca field obligatory.  Reorders fields into a standard order.
%
% bufix 3/12/14 JD
% Handles decimal sampling rates gracefully.
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% bufix 3/19/14 JD
% Fixed recTime field not including space for the 'all' cell in PCA output, resulting in crashes when edited.
%
% modified 3/24/14 JD
% Added cov field.
%
% bufix 4/8/14 JD
% Fixed not putting factor variance information in correct location when loading PCA .ept files from older versions of
% the EP Toolkit, resulting in "There are 0 factors that meet the minimum variance criterion" error messages when trying
% to autoPCA them.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 4/17/14 JD
% Added .covNum field.
%
% bufix 4/25/14 JD
% Added conversion of REF channel type to EEG for older EP files.
%
% bufix 6/12/14 JD
% Fixed blank keys field of events being produced without .key (e.g., .keys.keyCod instead of .keys.key.keyCode)
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bufix 8/27/14 JD
% Fixed crash when accessing a dataset in the working set with a different number of fields
% than that of the current EP format.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% modified 10/16/16 JD
% Added .stims field.
%
% modified 11/13/16 JD
% Added .calibration field.
%
% modified 6/13/17 JD
% Added .timeUnits field.
%
% modified 11/22/17 JD
% Added support for impedances field.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% bugfix 3/18/18 JD
% Fixed crash when .ept file has event with empty key.
%
% bugfix 6/12/18 JD
% Fixed not ensuring subject spec names and trial spec names are column vectors.
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

EPdata=[];

if theDataset > length(EPdataset.dataset)
    msg{1}='Error: there is no dataset corresponding to the one requested.';
    [msg]=ep_errorMsg(msg);
    return
end;

if isempty(theDataset)
    msg{1}='Error: which dataset to load in was not specified.';
    [msg]=ep_errorMsg(msg);
    return
end;

if length(theDataset) > 1
    msg{1}='Error: more than one dataset specified.';
    [msg]=ep_errorMsg(msg);
    return
end;

try
    eval(['load ''' EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(theDataset).dataName '.mat'';']);
catch
    msg{1}='Error: Working directory has been corrupted.  Delete the contents of the EPwork directory, except for the EPprefs file, and restart the EP Toolkit.';
    [msg]=ep_errorMsg(msg);
    beep
    return
end;

%backward compatibility
if isfield(EPdata,'implicit')
    if isfield(EPdata.implicit,'chantype')
        EPdata.implicit=rmfield(EPdata.implicit,'chantype');
    end;
end;

if ~isfield(EPdata,'analysis')
    EPdata.analysis=[];
end;

numSubs=length(EPdata.subNames);
numWaves=length(EPdata.cellNames);
numChan=length(EPdata.chanNames);

numEpochs=numWaves;
if strcmp(EPdata.dataType,'continuous')
    numEpochs=floor(size(EPdata.data,2)/ceil(EPdata.Fs)); %excess time points are tacked onto final epoch
    if numEpochs == 0
        numEpochs =1;
    end;
end;

if ~isfield(EPdata.analysis,'blinkTrial')
    EPdata.analysis.blinkTrial=zeros(numSubs,numEpochs);
end;
if ~isfield(EPdata.analysis,'saccadeTrial')
    EPdata.analysis.saccadeTrial=zeros(numSubs,numEpochs);
end;
if ~isfield(EPdata.analysis,'saccadeOnset')
    EPdata.analysis.saccadeOnset=zeros(numSubs,numEpochs);
end;
if ~isfield(EPdata.analysis,'moveTrial')
    EPdata.analysis.moveTrial=zeros(numSubs,numEpochs);
end;
if ~isfield(EPdata.analysis,'badTrials')
    EPdata.analysis.badTrials=zeros(numSubs,numEpochs);
end;
if ~isfield(EPdata.analysis,'badChans')
    EPdata.analysis.badChans=zeros(numSubs,numEpochs,numChan);
end;

if ~isfield(EPdata,'freqNames')
    EPdata.freqNames=[];
end;

if ~isfield(EPdata,'timeUnits')
    EPdata.timeUnits='ms';
end;

if ~isfield(EPdata,'cov')
    EPdata.cov=[]; 
end;

if ~isfield(EPdata,'facVecF')
    EPdata.facVecF=[];
end;

if ~isfield(EPdata,'facVar')
    if isfield(EPdata,'pca')
        if isfield(EPdata.pca,'facVar')
            EPdata.facVar=EPdata.pca.facVar;
        else
            EPdata.facVar=[];
        end;
    else
        EPdata.facVar=[];
    end;
end;

if ~isfield(EPdata,'facVarQ')
    if isfield(EPdata,'pca')
        if isfield(EPdata.pca,'facVarQ')
            EPdata.facVarQ=EPdata.pca.facVarQ;
        else
            EPdata.facVarQ=[];
        end;
    else
        EPdata.facVarQ=[];
    end;
end;

if ~isfield(EPdata,'noise')
    EPdata.noise=[];
end;

if ~isfield(EPdata,'std')
    EPdata.std=[];
end;

if ~isfield(EPdata,'stdCM')
    EPdata.stdCM=[];
end;

if ~isfield(EPdata,'cov')
    EPdata.cov=[]; 
end;

if ~isfield(EPdata,'reference')
    EPdata.reference=[];
end;

if ~isfield(EPdata.reference,'original')
    EPdata.reference.original=[];
end;

if ~isfield(EPdata.reference,'current')
    EPdata.reference.current=[];
end;

if ~isfield(EPdata.reference,'type')
    EPdata.reference.type=[];
end;

if ~isfield(EPdata,'pca')
    EPdata.pca=[]; 
end;

if isfield(EPdata,'pca')
    if ~isfield(EPdata.pca,'freqNames')
        EPdata.pca.freqNames=[];
    end;
    if ~isfield(EPdata.pca,'numFreqs')
        EPdata.pca.numFreqs=length(EPdata.pca.freqNames);
    end;
    if ~isfield(EPdata.pca,'recTime')
        EPdata.pca.recTime=[];
    end;
end;

if ~isfield(EPdata,'recTime')
    EPdata.recTime=[1:EPdata.Fs:EPdata.Fs*(length(EPdata.cellNames)-1)+1]';
elseif ~isempty(EPdata.facNames)
    if length(EPdata.recTime) == (length(EPdata.cellNames)-1)
        EPdata.recTime(end+1)=0; %repair effects of bug
    end;
end;

if ~isfield(EPdata,'covNum') || isempty(EPdata.covNum)
    EPdata.covNum=EPdata.avgNum;
end;

if isfield(EPdata,'power')
    if strcmp(EPdata.power,'Power')
        EPdata=sqrt(EPdata.data);
    end;
    EPdata=rmfield(EPdata,'power');
end;

if ~isfield(EPdata,'stims')
    EPdata.stims=struct('name',{},'image',{});
end;

if ~isfield(EPdata,'calibration')
    EPdata.calibration=[];
end;

if ~isfield(EPdata,'impedances')
    EPdata.impedances=[];
end;
if ~isfield(EPdata.impedances,'channels')
    EPdata.impedances.channels=[];
end;
if ~isfield(EPdata.impedances,'ground')
    EPdata.impedances.ground=[];
end;

if xor(isempty(EPdata.trialSpecNames),isempty(EPdata.trialSpecs))
    EPdata.trialSpecNames=[];
    EPdata.trialSpecs=[];
end;

if (size(EPdata.data,2) >= EPdata.Fs*2) && size(EPdata.analysis.blinkTrial,2) ==1 && strcmp(EPdata.dataType,'continuous')
    EPdata.analysis.blinkTrial=zeros(numSubs,numEpochs);
    EPdata.analysis.saccadeTrial=zeros(numSubs,numEpochs);
    EPdata.analysis.saccadeOnset=zeros(numSubs,numEpochs);
    EPdata.analysis.moveTrial=zeros(numSubs,numEpochs);
    EPdata.analysis.badTrials=zeros(numSubs,numEpochs);
    EPdata.analysis.badChans=zeros(numSubs,numEpochs,numChan);
    disp('Resetting artifact correction fields to zero to update to epochwise approach to continuous files.');
end;

if ~isempty(EPdata.events) %backward compatibility conversion
    for i=1:size(EPdata.events,1)
        for k=1:size(EPdata.events,2)
            if ~isempty(EPdata.events{i,k})
                if ~isfield(EPdata.events{i,k},'keys')
                    EPdata.events{i,k}(1).keys=struct('code','','data','','datatype','','description','');
                elseif isfield(EPdata.events{i,k}(1).keys,'key')
                    if isfield(EPdata.events{i,k}(1).keys(1).key,'keyCode')
                        for iEvent=1:length(EPdata.events{i,k})
                            newKeys=[];
                            for iKey=1:length(EPdata.events{i,k}(iEvent).keys)
                                if ~isempty(EPdata.events{i,k}(iEvent).keys(iKey).key)
                                    newKeys(iKey).code=EPdata.events{i,k}(iEvent).keys(iKey).key.keyCode;
                                    newKeys(iKey).data=EPdata.events{i,k}(iEvent).keys(iKey).key.data.data;
                                    newKeys(iKey).datatype=EPdata.events{i,k}(iEvent).keys(iKey).key.data.dataType;
                                    newKeys(iKey).description='';
                                end;
                            end;
                            EPdata.events{i,k}(iEvent).keys=newKeys;
                        end;
                    end;
                elseif (length(EPdata.events{i,k}(1).keys)>0) && isfield(EPdata.events{i,k}(1).keys(1),'keyCode')
                    for iEvent=1:length(EPdata.events{i,k})
                        newKeys=[];
                        for iKey=1:length(EPdata.events{i,k}(iEvent).keys)
                            newKeys(iKey).code=EPdata.events{i,k}(iEvent).keys(iKey).keyCode;
                            newKeys(iKey).data=EPdata.events{i,k}(iEvent).keys(iKey).data.data;
                            if isfield(EPdata.events{i,k}(iEvent).keys(iKey).data,'dataType')
                                newKeys(iKey).datatype=EPdata.events{i,k}(iEvent).keys(iKey).data.dataType;
                            else
                                newKeys(iKey).datatype='';
                            end;
                            if isfield(EPdata.events{i,k}(iEvent).keys(iKey).data,'description')
                                newKeys(iKey).description=EPdata.events{i,k}(iEvent).keys(iKey).data.description;
                            else
                                newKeys(iKey).description='';
                            end;
                        end;
                        EPdata.events{i,k}(iEvent).keys=newKeys;
                    end;
                end;
            end;
        end;
    end;
end;

MEGchans=find(strcmp('MEG',EPdata.chanTypes));
if ~isempty(MEGchans)
    chanTypes(MEGchans)='MGA';
    disp('Converting MEG channel types to MGA (axial gradiometer MEG).  If any of the MEG channels are actually planar or magnetometers, you will need to use the Edit function to correct them.');
end;

refChans=find(strcmp('REF',EPdata.chanTypes)); %assume that REF normally indicates original reference channels
if ~isempty(refChans)
    EPdata.chanTypes{refChans}='EEG'; %assume all REF channels are EEG channels
end;

%Ensure that one-dimensional string arrays are column vectors
EPdata.chanNames=EPdata.chanNames(:);
EPdata.timeNames=EPdata.timeNames(:);
EPdata.subNames=EPdata.subNames(:);
EPdata.cellNames=EPdata.cellNames(:);
EPdata.trialNames=EPdata.trialNames(:);
EPdata.facNames=EPdata.facNames(:);
EPdata.freqNames=EPdata.freqNames(:);
EPdata.relNames=EPdata.relNames(:);
EPdata.chanTypes=EPdata.chanTypes(:);
EPdata.subTypes=EPdata.subTypes(:);
EPdata.cellTypes=EPdata.cellTypes(:);
EPdata.facTypes=EPdata.facTypes(:);
EPdata.recTime=EPdata.recTime(:);
EPdata.trialSpecNames=EPdata.trialSpecNames(:);
EPdata.subjectSpecNames=EPdata.subjectSpecNames(:);

%ensure fields are in standard order.
[EPfieldNames]=ep_fieldNames;

modelEPdata=[];
for i=1:length(EPfieldNames)
    modelEPdata.(EPfieldNames{i})=[];
end;

dataFieldNames=fieldnames(EPdata);
if ~isequal(EPfieldNames,dataFieldNames)
    disp(['Modifying fields of ' EPdataset.dataset(theDataset).dataName ' to conform to current EPdata structure.']);
    if length(EPfieldNames) < length(dataFieldNames)
        for iField = 1:length(dataFieldNames)
            if ~any(strcmp(dataFieldNames{iField},EPfieldNames))
                EPdata=rmfield(EPdata,dataFieldNames{iField});
            end;
        end;
    end;
    if length(EPfieldNames) > length(dataFieldNames)
        for iField = 1:length(EPfieldNames)
            if ~any(strcmp(EPfieldNames{iField},dataFieldNames))
                eval(['EPdata.' EPfieldNames{iField} '=[];']);
            end;
        end;
    end;
    EPdata = orderfields(EPdata, modelEPdata);
end;









