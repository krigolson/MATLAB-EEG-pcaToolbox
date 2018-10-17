function [outData eloc]=ep_mergeEPfiles(files,mergeName)
% ep_mergeEPfiles - [outData]=ep_mergeEPfiles(files,mergeName) -
% Merge list of EP data files so that they are combined together.
% They must be either subject average or combined subject average files or single subject files from the same subject.
% Factor files are not supported.
% They must have the same set of channels and timepoints.
% If for a single subject file, must have the same trial spec names.
% The result is a combined set of cells and a combined number of subjects.
% Any subjects without the full set of cells will have zeroed data for the missing cells.
% The information about the experiment and so forth are taken from the first file in the list.
%
%Inputs
% Optional
%   files: A 2D cell array with file names to merge (including path names) in the first column.
%          the remaining columns are input arguments for the readData function.
%   mergeName: the name of the resulting merged file.
%
%Outputs:
%  EPdata         : Structured array with the data and accompanying information.  See readData.
%  eloc:    The raw eloc prior to any editing.
%
%History
%  by Joseph Dien (3/4/09)
%  jdien07@mac.com
%
% modified 5/30/09 JD
% Added cellType, subType, chanType, facName, facType, and analysis fields.  Reads in files directly rather than .mat files.
%
% modified 6/26/09 JD
% added mergeName input and dataName field.
%
%  bugfix 10/28/09 JD
%  Crash due to not encoding cell and subject labels correctly.
%  Not filling in cellType and trialSpec fields correctly for single-trial data.
%
%  bugfix 1/18/10 JD
%  Fixed crash in EP when merge operation aborted due to problem with files.
%
%  bugfix 1/23/10 JD
%  Fixed contents of analysis fields being not arranged correctly.
%
% modified 2/11/10 JD
% analysis fields no longer optional.
%
%  bugfix 4/30/10 JD
%  Fixed crash when reading single_trial data using single file mode.
%
%  bugfix 5/23/10 JD
%  Fixed crash when merging average files.
%  Added 'exact' keyword to strmatch commands to ensure that the cell and subject names are matched exactly rather than
%  based on the beginning of the names.
%
% modified 5/26/10 JD
% Added eloc to parameters passed to files to handle cases where ced is non-empty but eloc is empty so doesn't have to
% keep accessing .ced file on every pass.
% Fixed not successfully merging average files together.
%
%  bugfix 6/2/10 JD
% Only adds trial names for single trial data files.
% 
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
% 
% modified 1/24/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
% 
% modified 3/25/12 JD
% Added support for frequency files.
%
% bugfix 10/18/12 JD
% Fixed subNames, subTypes, cellNames, and cellTypes fields not necessarily being column vectors, resulting in crashes in other parts of the Toolkit.
%
% bugfix 2/6/13 JD
% Fixed power field not being included in merged files, causing crashes elsewhere.
%
% bugfix 5/20/13 JD
% Fixed single file mode single-trial data files not being merged successfully.
%
% modified 5/20/13 JD
% Channels not present in the initial file, as in channels identified as being BAD in the CED file, simply dropped from succeeding files.
% Added original eloc output.
%
% modified 10/9/13 JD
% Added recTime field.
%
% modified 2/27/14 JD
% Fields output in standard order.
%
% modified 3/24/14 JD
% Added .cov field.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% modified 6/1/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure, including support for complex numbers.
%
% modified 9/4/15 JD
% Added trial specs for average files.
%
% bugfix 11/4/15 JD
% Fixed crash when merging files that do not have subject specs.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% modified 10/16/16 JD
% Added .stims field.
%
% modified 11/14/16 JD
% No longer assumes subject specs are the same for all merged files. 
% Added .calibration field.
%
% bugfix 5/19/17 JD
% Fixed crash when merging multiple subjects.
%
% modified 6/13/17 JD
% Added .timeUnits field.
%
% modified 11/22/17 JD
% Added support for impedances field.
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


outData=[];
eloc=[];

if nargin < 1
    files=[];
end

if isempty(files)
    [files, pathname] = uigetfile('*.*','Open:','MultiSelect','on');
    activeDirectory=pathname;
    if ~iscell(files)
        temp=files;
        files=[];
        files{1}=temp;
    end;
    if files{1}==0
        msg{1}='No filenames selected. You have to click on a name';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if ~iscell(files)
        temp=files;
        files=[];
        files{1}=temp;
    end;
else
    for theFile=1:size(files,1)
        if ~exist(char(files{theFile,1}),'file')
            msg{1}=['Error: The file ' files{theFile,2} ' is not in the directory.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
    end;
end;

numFiles = size(files,1);

if numFiles==1
    msg{1}='Only one file selected.';
    [msg]=ep_errorMsg(msg);
    return
end

totalSubjectSpecs=cell(0);
totalSubjectNames=cell(0);
totalCellNames=cell(0);
uniqueCellNames=cell(0);
uniqueFreqNames=[];
uniqueRelNames=[];
ced=[];
eloc=[];

%compile data from each file
for theFile=1:numFiles
    thisFile=files{theFile,1};
    temp=[];
    temp{1}='file';
    temp(2:size(files,2)+1)=files(theFile,:);
    if ~isempty(ced) && ~any(strcmp('ced',files(theFile,:)))
        temp{end+1}='ced';
        temp{end+1}=ced;
    end;
    if ~isempty(eloc)
        temp{end+1}='eloc';
        temp{end+1}=eloc;
    end;
    
    if theFile >1
        temp{end+1}='silent';
        temp{end+1}='on';
    end;
    
    [data eloc]=ep_readData(temp);
    
    if isempty(data)
        outData=[];
        msg{1}='File not read successfully.  Aborting effort to merge data files.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if strcmp(data.dataType,'single_trial')
        if theFile==1
            theSingleSubject=data.subNames{1};
        else
            if ~strcmp(data.subNames{1},theSingleSubject)
                outData=[];
                msg{1}=['Error: When merging single trial data, all must be from the same subject (' theSingleSubject ').'];
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
    elseif ~any(strcmp(data.dataType,{'average'}))
        outData=[];
        msg{1}=['Error: The file ' thisFile ' is the ' data.dataType ' data type.'];
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    %information that is assumed to be identical between all the files and is culled from the first such file
    if theFile ==1
        outData.dataType=data.dataType;
        outData.fileFormat = data.fileFormat;
        outData.montage = data.montage;
        outData.chanNames = data.chanNames;
        outData.timeNames = data.timeNames;
        outData.subNames = data.subNames;
        outData.chanTypes= data.chanTypes;
        outData.timeUnits= data.timeUnits;
        outData.cellTypes = data.cellTypes;
        outData.subjectSpecNames = data.subjectSpecNames;
        outData.facNames= data.facNames;
        outData.facTypes = data.facTypes;
        outData.freqNames = data.freqNames;
        outData.relNames = data.relNames;
        try
            EPver=ver('EP_Toolkit');
        catch
            EPver='unavailable'; %workaround for bug in earlier version of Matlab
        end;
        outData.EPver = EPver;
        outData.ver = ver;
        outData.date = date;
        outData.Fs = data.Fs;
        outData.baseline = data.baseline;
        outData.ename = data.ename;
        outData.dataName = mergeName;
        outData.trialSpecNames = data.trialSpecNames;
        outData.trialNames = [];
        outData.trialSpecs = [];
        outData.fileName = data.fileName;
        outData.ced = data.ced;
        outData.eloc = data.eloc;
        outData.reference = data.reference;
        outData.implicit = data.implicit;
        outData.facVecS = data.facVecS;
        outData.facVecT = data.facVecT;
        outData.facVecF = data.facVecF;
        outData.facVar = data.facVar;
        outData.facVarQ = data.facVarQ;
        outData.pca = data.pca;
        outData.stims = data.stims;
        outData.calibration = data.calibration;
        outData.impedances = data.impedances;
        ced = data.ced;
        if isempty(ced)
            ced='none';
        end;
        numChans = size(data.data,1);
        numPoints = size(data.data,2);
        numEvents = length(data.events);
        
        newBlinkTrial=[];
        newMoveTrial=[];
        newBadTrials=[];
        newBadChans=[];
    else
        
        if ~isempty(setxor(outData.chanNames,data.chanNames))
            if isempty(setdiff(outData.chanNames,data.chanNames))
                keepChans=[];
                for i=1:length(data.chanNames)
                    if any(strcmp(data.chanNames{i},outData.chanNames))
                        keepChans(end+1)=i;
                    end;
                end;
                data=ep_selectData(data,{keepChans,[],[],[],[],[]});
                disp(['Dropping channels in ' thisFile ' that are not in the initial file, including those marked as BAD in the CED file.']);
            else
                outData=[];
                msg{1}=['Error: The file ' thisFile ' has different channels from the first file.'];
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        if ~isempty(setxor(outData.timeNames,data.timeNames))
            outData=[];
            msg{1}=['Error: The file ' thisFile ' has different time points from the first file.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
        if ~isempty(setxor(outData.trialSpecNames,data.trialSpecNames)) && strcmp(data.dataType,'single_trial')
            outData=[];
            msg{1}=['Error: The file ' thisFile ' has different trial specs from the first file.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
        if ~isempty(data.stims)
            for iStim=1:length(data.stims)
                if ~any(strcmp(data.stims(iStim).name,{outData.stims.name}))
                    outData.stims(end+1)=data.stims(iStim);
                end;
            end;
        end;
    end;
    
    newTrialNames{theFile}= data.trialNames;
    newTrialSpecs{theFile}= data.trialSpecs;
    
    newCellNames{theFile}= data.cellNames;
    newData{theFile}= data.data;
    newNoise{theFile}= data.noise;
    newStd{theFile}= data.std;
    newStdCM{theFile}= data.stdCM;
    newCov{theFile}= data.cov;
    newSubNames{theFile}= data.subNames;
    newSubTypes{theFile}= data.subTypes;
    newCellTypes{theFile}= data.cellTypes;
    newSubjectSpecs{theFile}= data.subjectSpecs;
    newEvents{theFile}= data.events;
    newAvgNum{theFile}= data.avgNum;
    newCovNum{theFile}= data.covNum;
    newSubNum{theFile}= data.subNum;
    newFreqNames{theFile}=data.freqNames;
    newRelNames{theFile}=data.relNames;
    
    newBlinkTrial{theFile}= data.analysis.blinkTrial;
    newSaccadeTrial{theFile}= data.analysis.saccadeTrial;
    newSaccadeOnset{theFile}= data.analysis.saccadeOnset;
    newMoveTrial{theFile}= data.analysis.moveTrial;
    newBadTrials{theFile}= data.analysis.badTrials;
    newBadChans{theFile}= data.analysis.badChans;
    
    for i=1:length(data.cellNames)
        totalCellNames(end+1)=data.cellNames(i);
        if ~any(strcmp(data.cellNames{i},uniqueCellNames))
            uniqueCellNames{end+1}=char(data.cellNames{i});
        end;
    end;
    
    for i=1:length(data.freqNames)
        if ~isempty(setdiff(data.freqNames(i),uniqueFreqNames))
            uniqueFreqNames(end+1)=data.freqNames(i);
        end;
    end;

    for i=1:length(data.relNames)
        if ~isempty(setdiff(data.relNames(i),uniqueRelNames))
            uniqueRelNames(end+1)=data.relNames(i);
        end;
    end;
    
    numNewSubs=length(data.subNames);
    lastOldSub=length(outData.subNames);
    if isempty(lastOldSub)
        lastOldSub=0;
    end;
    for iSpec=1:length(data.subjectSpecs)
        theSpec=find(strcmp(data.subjectSpecNames{iSpec},outData.subjectSpecNames));
        if ~isempty(theSpec)
            outData.subjectSpecs{lastOldSub+1:lastOldSub+numNewSubs,theSpec}=data.subjectSpecs{:,iSpec};
        else
            outData.subjectSpecNames{end+1}=data.subjectSpecNames{iSpec};
            outData.subjectSpecs{lastOldSub+1:lastOldSub+numNewSubs,end+1}=data.subjectSpecs{:,iSpec};
        end;
    end;
    
    for i=1:length(data.subNames)
        if ~any(strcmp(data.subNames{i},totalSubjectNames))
            totalSubjectNames(end+1,1)=data.subNames(i);
            if ~isempty(outData.subjectSpecNames)
                totalSubjectSpecs{end+1,1}=[];
                if ~isempty(data.subjectSpecs)
                    totalSubjectSpecs(end,1)=data.subjectSpecs(i);
                end;
            end;
        end;
    end;
end;

numCells=length(uniqueCellNames);
numSubs=length(totalSubjectNames);
numFreqs=length(uniqueFreqNames);
numRels=length(uniqueRelNames);

if strcmp(outData.dataType,'average')
    numWaves=numCells;
elseif strcmp(data.dataType,'single_trial')
    numWaves=length(totalCellNames);
end;

outData.data=zeros(numChans,numPoints,numWaves,numSubs,1,max(1,numFreqs),max(1,numRels));
outData.noise=zeros(numChans,numPoints,numWaves,numSubs,1);
outData.std=zeros(numChans,numPoints,numWaves,numSubs,1,max(1,numFreqs));
outData.stdCM=zeros(numChans,numPoints,numWaves,numSubs,1,max(1,numFreqs));
outData.cov.covMatrix=NaN(numSubs,numChans,numChans);
outData.cov.Nq=NaN;
outData.facData=[];
outData.avgNum=zeros(numSubs,numCells);
outData.covNum=zeros(numSubs,numCells);
outData.subNum=zeros(numSubs,numCells);
outData.events=cell(numSubs,numCells);
outData.analysis.blinkTrial=zeros(numSubs,numWaves);
outData.analysis.saccadeTrial=zeros(numSubs,numWaves);
outData.analysis.saccadeOnset=zeros(numSubs,numWaves);
outData.analysis.moveTrial=zeros(numSubs,numWaves);
outData.analysis.badTrials=zeros(numSubs,numWaves);
outData.analysis.badChans=zeros(numSubs,numWaves,numChans);

outData.subNames=totalSubjectNames;
outData.subjectSpecs=totalSubjectSpecs;

outData.recTime=[1:outData.Fs:outData.Fs*(numWaves-1)+1];

%merge the data from the files together
trialCount=zeros(numCells,1);
outCount=0;
for theFile=1:numFiles
    for theSub=1:length(newSubNames{theFile})
        newSub=strmatch(newSubNames{theFile}{theSub},totalSubjectNames,'exact');
        outData.subTypes(newSub,1)=newSubTypes{theFile}(theSub);
        for theCell=1:length(newCellNames{theFile})
            newCell=strmatch(newCellNames{theFile}{theCell},uniqueCellNames,'exact');
            newWave=newCell;
            if strcmp(data.dataType,'single_trial')
                newWave=outCount+theCell;
            end;
            for theSpec=1:size(newTrialSpecs{theFile},2)
                outData.trialSpecs{newWave,theSpec,newSub}=newTrialSpecs{theFile}{theCell,theSpec,theSub};
            end;
            trialCount(newCell)=trialCount(newCell)+1;
            if strcmp(data.dataType,'single_trial')
                outData.trialNames(newWave,1)=trialCount(newCell);
            end;
            outData.analysis.blinkTrial(newSub,newWave)=newBlinkTrial{theFile}(theSub,theCell);
            outData.analysis.saccadeTrial(newSub,newWave)=newSaccadeTrial{theFile}(theSub,theCell);
            outData.analysis.saccadeOnset(newSub,newWave)=newSaccadeOnset{theFile}(theSub,theCell);
            outData.analysis.moveTrial(newSub,newWave)=newMoveTrial{theFile}(theSub,theCell);
            outData.analysis.badTrials(newSub,newWave)=newBadTrials{theFile}(theSub,theCell);
            outData.analysis.badChans(newSub,newWave,:)=newBadChans{theFile}(theSub,theCell,:);
            outData.cellNames(newWave,1)=uniqueCellNames(newCell);
            outData.cellTypes(newWave,1)=newCellTypes{theFile}(theCell);
            
            for theFreq=1:max(1,length(newFreqNames{theFile}))
                if ~isempty(newFreqNames{theFile})
                    newFreq=find(newFreqNames{theFile}(theFreq) == uniqueFreqNames);
                    outData.freqNames(newFreq,1)=uniqueFreqNames(newFreq);
                else
                    newFreq=1;
                end;
                
                for theRel=1:max(1,length(newRelNames{theFile}))
                    if ~isempty(newRelNames{theFile})
                        newRel=find(newRelNames{theFile}(theFreq) == uniqueRelNames);
                        outData.relNames(newRel,1)=uniqueRelNames(newRel);
                    else
                        newRel=1;
                    end;

                    if any(any(squeeze(outData.data(:,:,newWave,newSub,1,newFreq,newRel))))
                        disp(' ');
                        disp('**************************************************************');
                        disp(['Warning: The file ' thisFile ' duplicated the data for subject ' totalSubjectNames{newSub} ' for cell ' totalCellNames{newCell} '.']);
                        disp('**************************************************************');
                        disp(' ');
                    end;
                    outData.data(:,:,newWave,newSub,1,newFreq)=newData{theFile}(:,:,theCell,theSub,1,theFreq,newRel);
                end;
                if ~isempty(newNoise{theFile})
                    outData.noise(:,:,newWave,newSub,1)=newNoise{theFile}(:,:,theCell,theSub,1);
                end;
                if ~isempty(newStd{theFile})
                    outData.std(:,:,newWave,newSub,1,newFreq)=newStd{theFile}(:,:,theCell,theSub,1,theFreq);
                end;
                if ~isempty(newStdCM{theFile})
                    outData.stdCM(:,:,newWave,newSub,1,newFreq)=newStdCM{theFile}(:,:,theCell,theSub,1,theFreq);
                end;
                if ~isempty(newCov{theFile})
                    outData.cov.covMatrix(newSub,:,:)=newCov{theFile}.covMatrix(theSub,:,:);
                    outData.cov.Nq(Nq)=newCov{theFile}.Nq(theSub);
                end;
            end;
            outData.avgNum(newSub,newWave)=newAvgNum{theFile}(theSub,theCell);
            outData.covNum(newSub,newWave)=newCovNum{theFile}(theSub,theCell);
            outData.subNum(newSub,newWave)=newSubNum{theFile}(theSub,theCell);
            outData.events(newSub,newWave)=newEvents{theFile}(theSub,theCell);
        end;
        outCount=outCount+length(newCellNames{theFile});
    end;
end;

try
    EPver=ver('EP_Toolkit');
catch
    EPver='unavailable'; %workaround for bug in earlier version of Matlab
end;
outData.EPver=EPver;
outData.ver=ver;
outData.date=date;
outData.history={'mergeEPfiles',files};
outData.facData=[]; %factor files not supported so assumed no such data.
if isempty(outData.trialSpecs)
    outData.trialSpecs=cell(size(outData.data,3),0,size(outData.data,4));
end;
outData.pca=[];

if ~any(any(any(any(any(any(outData.std))))))
    outData.std=[];
end;

if ~any(any(any(any(any(any(outData.stdCM))))))
    outData.stdCM=[];
end;

if ~any(any(any(any(any(outData.noise)))))
    outData.noise=[];
end;

if all(all(all(isnan(outData.cov.covMatrix))))
    outData.cov=[];
end;

%ensure fields are in standard order.
[EPfieldNames]=ep_fieldNames;

modelEPdata=[];
for i=1:length(EPfieldNames)
    modelEPdata.(EPfieldNames{i})=[];
end;

if ~isequal(EPfieldNames,fieldnames(outData))
    outData = orderfields(outData, modelEPdata);
end;

[err]=ep_checkEPfile(outData);

if err
    outData=[];
    msg{1}='Defective file will not be loaded.  File can be fixed by loading manually into Matlab and editing, as in load(''filename.ept'',''-mat''); and then save(''filename.ept'',''-mat'');';
    [msg]=ep_errorMsg(msg);
    return
end;
