function [err]=ep_checkEPfile(data);
%  [err]=ep_checkEPfile(data);
%       Checks EP data file for problems.
%
%Inputs:
%  data         : Structured array with the data and accompanying information.  See readData.
%
%Outputs:
%  err        : Array of detected error codes.  If no errors, returns zero.

%History:
%  by Joseph Dien (3/4/09)
%  jdien07@mac.com
%
%
% bugfix 10/2/13 JD
% Fixed rejecting files if sampling rate and time names different past three decimals due to rounding errors.
%
% bugfix 1/12/14 JD
% Fixed giving errors when .facdata is not empty but .data is.
%
% bugfix 10/25/16 JD
% Fixed giving errors with averaged data when there are no trial specs.
%
% bugfix 3/15/18 JD
% Fixed giving errors with averaged data when there are more subjects than channels.
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

err=[];

%largest is 130

if ~isfield(data,'data')
    err = [err ;1];
    disp('The .data field is missing.');
    data.data=[];
end

if ~isfield(data,'montage')
    err = [err ;2];
    disp('The .montage field is missing.');
end;

if ~isfield(data,'ced')
    err = [err ;50];
    disp('The .ced field is missing.');
end;

if ~isfield(data,'fileFormat')
    err = [err ;3];
    disp('The .fileFormat field is missing.');
end;

if ~isfield(data,'dataType')
    err = [err ;4];
    disp('The .dataType field is missing.');
    data.dataType=[];
end;

if ~isfield(data,'chanTypes')
    err = [err ;46];
    disp('The .chanTypes field is missing.');
else
    if size(data.chanTypes,2) > 1
        err = [err ;96];
        disp('The .chanTypes field is not a column vector.');
    end;
end;

if ~isfield(data,'subTypes')
    err = [err ;47];
    disp('The .subTypes field is missing.');
else
    if size(data.subTypes,2) > 1
        err = [err ;95];
        disp('The .subTypes field is not a column vector.');
    end;
    if length(data.subTypes) ~= size(data.data,4)
        err = [err ;41];
        disp('Number of subject types does not match the data.');
    end;
end;

if ~isfield(data,'cellTypes')
    err = [err ;48];
    disp('The .cellTypes field is missing.');
else
    if size(data.cellTypes,2) > 1
        err = [err ;94];
        disp('The .cellTypes field is not a column vector.');
    end;
    if length(data.cellTypes) ~= size(data.data,3)
        err = [err ;42];
        disp('Number of cell types does not match the data.');
    end;
end;

if ~isfield(data,'facTypes')
    err = [err ;49];
    disp('The .facTypes field is missing.');
else
    if size(data.facTypes,2) > 1
        err = [err ;93];
        disp('The .facTypes field is not a column vector.');
    end;
end;

if ~isfield(data,'chanNames')
    err = [err ;5];
    disp('The .chanNames field is missing.');
else
    if size(data.chanNames,2) > 1
        err = [err ;92];
        disp('The .chanNames field is not a column vector.');
    end;
    if any(cellfun('isempty',data.chanNames))
        err = [err ;77];
        disp('Empty channel name present.');
    end;

end;

if ~isfield(data,'timeNames')
    err = [err ;6];
    disp('The .timeNames field is missing.');
else
    if size(data.timeNames,2) > 1
        err = [err ;91];
        disp('The .timeNames field is not a column vector.');
    end;
end;

if ~isfield(data,'subNames')
    err = [err ;7];
    disp('The .subNames field is missing.');
else
    if size(data.subNames,2) > 1
        err = [err ;87];
        disp('The .subNames field is not a column vector.');
    end;
    if length(data.subNames) ~= size(data.data,4)
        err = [err ;53];
        disp('Number of subject names does not match the data.');
    end;
    if any(cellfun('isempty',data.subNames))
        err = [err ;76];
        disp('Empty subject name present.');
    end;
end;

if ~isfield(data,'cellNames')
    err = [err ;8];
    disp('The .cellNames field is missing.');
else
    if size(data.cellNames,2) > 1
        err = [err ;88];
        disp('The .cellNames field is not a column vector.');
    end;
    if length(data.cellNames) ~= size(data.data,3)
        err = [err ;52];
        disp('Number of cell names does not match the data.');
    end;
    if any(cellfun('isempty',data.cellNames))
        err = [err ;75];
        disp('Empty cell name present.');
    end;
end;

if ~isfield(data,'facNames')
    err = [err ;45];
    disp('The .facNames field is missing.');
else
    if size(data.facNames,2) > 1
        err = [err ;89];
        disp('The .facNames field is not a column vector.');
    end;
end;

if ~isfield(data,'trialNames')
    err = [err ;9];
    disp('The .trialNames field is missing.');
else
    if size(data.trialNames,2) > 1
        err = [err ;90];
        disp('The .trialNames field is not a column vector.');
    end;
end;

if ~isfield(data,'relNames')
    err = [err ;114];
    disp('The .relNames field is missing.');
else
    if size(data.relNames,2) > 1
        err = [err ;92];
        disp('The .relNames field is not a column vector.');
    end;
    if length(data.relNames) ~= length(data.relNames)
        err = [err ;117];
        disp('The .relNames field is not the same size as the chanNames field.');
    end;
end;

if ~isfield(data,'EPver')
    err = [err ;10];
    disp('The .EPver field is missing.');
end;

if ~isfield(data,'ver')
    err = [err ;11];
    disp('The .ver field is missing.');
end;

if ~isfield(data,'date')
    err = [err ;12];
    disp('The .date field is missing.');
end;

if ~isfield(data,'Fs')
    err = [err ;13];
    disp('The .Fs field is missing.');
    data.Fs=[];
end;

if ~isfield(data,'baseline')
    err = [err ;14];
    disp('The .baseline field is missing.');
end;

if ~isfield(data,'ename')
    err = [err ;15];
    disp('The .ename field is missing.');
end;

if ~isfield(data,'dataName')
    err = [err ;68];
    disp('The .dataName field is missing.');
end;

if ~isfield(data,'trialSpecs')
    err = [err ;16];
    disp('The .trialSpecs field is missing.');
end;

if ~isfield(data,'trialSpecNames')
    err = [err ;17];
    disp('The .trialSpecNames field is missing.');
end;

if ~isfield(data,'subjectSpecs')
    err = [err ;18];
    disp('The .subjectSpecs field is missing.');
end;

if ~isfield(data,'subjectSpecNames')
    err = [err ;19];
    disp('The .subjectSpecNames field is missing.');
end;

if ~isfield(data,'events')
    err = [err ;20];
    disp('The .events field is missing.');
elseif ~isempty(data.events)
    for i=1:size(data.events,1)
        for k=1:size(data.events,2)
            if ~isempty(data.events{i,k})
                if ~isfield(data.events{i,k},'type')
                    err = [err ;100];
                    disp('The .events.type field is missing.');
                end;
                if ~isfield(data.events{i,k},'sample')
                    err = [err ;101];
                    disp('The .events.sample field is missing.');
                end;
                if ~isfield(data.events{i,k},'value')
                    err = [err ;102];
                    disp('The .events.value field is missing.');
                end;
                if ~isfield(data.events{i,k},'duration')
                    err = [err ;103];
                    disp('The .events.duration field is missing.');
                end;
                if ~isfield(data.events{i,k},'keys')
                    err = [err ;104];
                    disp('The .events.keys field is missing.');
                end;
            end;
        end;
    end;
end;

if ~isfield(data,'avgNum')
    err = [err ;21];
    disp('The .avgNum field is missing.');
end;

if ~isfield(data,'covNum')
    err = [err ;115];
    disp('The .covNum field is missing.');
end;

if ~isfield(data,'subNum')
    err = [err ;39];
    disp('The .subNum field is missing.');
else
    if size(data.subNum,1) ~= size(data.data,4) || size(data.subNum,2) ~= size(data.data,3)
        err = [err ;40];
        disp('Number of subjects going into averages does not match the data.');
    end;
end;

if ~isfield(data,'fileName')
    err = [err ;31];
    disp('The .fileName field is missing.');
end;

if ~isfield(data,'history')
    err = [err ;32];
    disp('The .history field is missing.');
end;

if ~isfield(data,'facData')
    err = [err ;62];
    disp('The .facData field is missing.');
end;

if ~isfield(data,'reference')
    err = [err ;73];
    disp('The .reference field is missing.');
else
    if ~isfield(data.reference,'original')
        err = [err ;78];
        disp('The .reference.original field is missing.');
    end;
    if ~isfield(data.reference,'current')
        err = [err ;79];
        disp('The .reference.current field is missing.');
    end;
    if ~isfield(data.reference,'type')
        err = [err ;80];
        disp('The .reference.type field is missing.');
    end;
end;

if ~isfield(data,'noise')
    err = [err ;85];
    disp('The .noise field is missing.');
end;

if ~isfield(data,'std')
    err = [err ;86];
    disp('The .std field is missing.');
end;

if ~isfield(data,'stdCM')
    err = [err ;129];
    disp('The .stdCM field is missing.');
end;

numEpochs=size(data.data,3);
if strcmp(data.dataType,'continuous')
    numEpochs=floor(size(data.data,2)/ceil(data.Fs)); %excess time points are tacked onto final epoch
    if numEpochs == 0
        numEpochs =1;
    end;
end;

if ~isfield(data,'analysis')
    err = [err ;63];
    disp('The .analysis field is missing.');
else
    if ~isfield(data.analysis,'blinkTrial')
        err = [err ;64];
        disp('The .analysis.blinkTrial field is missing.');
    else
        if size(data.analysis.blinkTrial,2) ~= numEpochs || size(data.analysis.blinkTrial,1) ~= size(data.data,4)
            err = [err ;44];
            disp('Number of blink trials does not match the data.');
        end;
        
    end;
    if ~isfield(data.analysis,'saccadeTrial')
        err = [err ;69];
        disp('The .analysis.saccadeTrial field is missing.');
    else
        if size(data.analysis.saccadeTrial,2) ~= numEpochs || size(data.analysis.saccadeTrial,1) ~= size(data.data,4)
            err = [err ;71];
            disp('Number of saccade trials does not match the data.');
        end;
    end;
    if ~isfield(data.analysis,'saccadeOnset')
        err = [err ;70];
        disp('The .analysis.saccadeOnset field is missing.');
    else
        if size(data.analysis.saccadeOnset,2) ~= numEpochs || size(data.analysis.saccadeOnset,1) ~= size(data.data,4)
            err = [err ;72];
            disp('Number of saccade onset trials does not match the data.');
        end;
    end;
    if ~isfield(data.analysis,'moveTrial')
        err = [err ;65];
        disp('The .analysis.moveTrial field is missing.');
    else
        if size(data.analysis.moveTrial,2) ~= numEpochs || size(data.analysis.moveTrial,1) ~= size(data.data,4)
            err = [err ;45];
            disp('Number of movement trials does not match the data.');
        end;
    end;
    if ~isfield(data.analysis,'badTrials')
        err = [err ;66];
        disp('The .analysis.badTrials field is missing.');
    else
        if size(data.analysis.badTrials,2) ~= numEpochs || size(data.analysis.badTrials,1) ~= size(data.data,4)
            err = [err ;46];
            disp('Number of bad trials does not match the data.');
        end;
    end;
    if ~isfield(data.analysis,'badChans')
        err = [err ;67];
        disp('The .analysis.badChans field is missing.');
    else
        if ((size(data.analysis.badChans,3) ~= size(data.data,1)) && (size(data.analysis.badChans,3) ~= size(data.facVecS,1))) || (size(data.analysis.badChans,2) ~= numEpochs) || size(data.analysis.badChans,1) ~= size(data.data,4)
            err = [err ;47];
            disp('Number of bad channels does not match the data.');
        end;
    end;
end;


if ~isfield(data,'facVecS')
    err = [err ;56];
    disp('The .facVecS field is missing.');
elseif ~isempty(data.facVecS)
    if (size(data.facVecS,1) ~= length(data.chanNames))
        err = [err ;22];
        disp('Number of channel names does not match the data.');
    end;
    if ~isempty(data.facVecS)
        if size(data.facVecS,2) ~= size(data.data,5)
            err = [err ;54];
            disp('Number of factors does not match spatial factor vector.');
        end;
    end;
else
    if (length(data.chanNames) ~= size(data.data,1)) && ~isempty(data.data)
        err = [err ;22];
        disp('Number of channel names does not match the data.');
    end;
end;

numFacData=0;
if ~isfield(data,'facData')
    err = [err ;108];
    disp('The .facData field is missing.');
elseif ~isempty(data.facData)
    if length(data.chanNames) ~= size(data.facData,1)
        err = [err ;58];
        disp('Number of channel names does not match the combined factor data.');
    end;
    if length(data.timeNames) ~= size(data.facData,2)
        if ~(isempty(data.timeNames) && (size(data.data,2) == 1)) %frequency data
            err = [err ;59];
            disp('Number of sample names does not match the combined factor data.');
        end;
    end;
    if length(data.cellNames) ~= size(data.facData,3)
        err = [err ;60];
        disp('Number of cell names does not match the combined factor data.');
    end;
    
    if length(data.subNames) ~= size(data.facData,4)
        err = [err ;61];
        disp('Number of subject names does not match the combined factor data.');
    end;
    if ~isempty(data.facData)
        numFacData=size(data.facData,5);
    end;
end;

if ~isfield(data,'facVecT')
    err = [err ;57];
    disp('The .facVecT field is missing.');
elseif ~isempty(data.facVecT)
    if (size(data.facVecT,1) ~= length(data.timeNames)) || (size(data.data,2) > 1)
        err = [err ;51];
        disp('Number of sample names does not match the data.');
    end;
    if ~isempty(data.facVecT)
        if size(data.facVecT,2) ~= size(data.data,5)
            err = [err ;55];
            disp('Number of factors does not match temporal factor vector.');
        end;
    end;
else
    if isempty(data.timeNames)
        if size(data.data,2) ~= 1 %no timeNames and second dimension of data array is size one means frequency data
            err = [err ;106];
            disp('Number of sample names does not match the data.');
        end;
    elseif length(data.timeNames) ~= size(data.data,2) && ~isempty(data.data)
        err = [err ;107];
        disp('Number of sample names does not match the data.');
    end;
end;

if ~isfield(data,'facVecF')
    err = [err ;81];
    disp('The .facVecF field is missing.');
elseif ~isempty(data.facVecF)
    if (size(data.facVecF,1) ~= length(data.freqNames)) || (size(data.data,6) > 1)
        err = [err ;82];
        disp('Number of frequency names does not match the data.');
    end;
else
    if isempty(data.freqNames)
        if size(data.data,6) ~= 1 %no freqNames and sixth dimension of data array is size one means not frequency data
            err = [err ;83];
            disp('Number of frequency names does not match the data.');
        end;
    elseif length(data.freqNames) ~= size(data.data,6) && ~isempty(data.data)
        err = [err ;84];
        disp('Number of frequency names does not match the data.');
    end;
end;

if ~isfield(data,'pca')
    err = [err ;110];
    disp('The .pca field is missing.');
end;

if ~isfield(data,'facVar')
    err = [err ;111];
    disp('The .facVar field is missing.');
end;

if ~isfield(data,'facVarQ')
    err = [err ;112];
    disp('The .facVarQ field is missing.');
end;

if ~isfield(data,'recTime')
    err = [err ;99];
    disp('The .recTime field is missing.');
else
    if length(data.recTime) ~= size(data.data,3)
        err = [err ;105];
        disp('Number of recording times does not match the data.');
    end;
end;

if ~isfield(data,'stims')
    err = [err ;121];
    disp('The .stims field is missing.');
end;

if ~isfield(data,'calibration')
    err = [err ;122];
    disp('The .calibration field is missing.');
end;

if ~isfield(data,'timeUnits')
    err = [err ;124];
    disp('The .timeUnits field is missing.');
end;

if isfield(data,'facNames')
    if length(data.facNames) ~= (size(data.data,5)+numFacData) && not (isempty(data.facNames) && size(data.data,5)==1)
        err = [err ;44];
        disp('Number of factor names does not match the data.');
    end;
end;

if ~isfield(data,'impedances')
    err = [err ;125];
    disp('The .impedances field is missing.');
else
    if ~isfield(data.impedances,'channels')
        err = [err ;126];
        disp('The .impedances.channels field is missing.');
    else
        if ~isempty(data.impedances.channels)
            if length(data.chanNames) ~= size(data.impedances.channels,1)
                err = [err ;128];
                disp('The number of channel names does not match the number of channel impedances.');
            end;
            if length(data.subNames) ~= size(data.impedances.channels,2)
                err = [err ;130];
                disp('The number of subject names does not match the number of channel impedances.');
            end;
        end;
    end;
    if ~isfield(data.impedances,'ground')
        err = [err ;127];
        disp('The .impedances.ground field is missing.');
    end;
end;

if isfield(data,'subjectSpecs')
    if (size(data.subjectSpecs,1) ~= size(data.data,4)) && ~isempty(data.subjectSpecs)
        err = [err ;27];
        disp('Number of subject specs does not match the data.');
    end;
    
    if size(data.subjectSpecs,2) ~= length(data.subjectSpecNames)
        err = [err ;28];
        disp('Number of subject spec names does not match the subject specs.');
    end;
end;

if isfield(data,'events')
    if size(data.events,1) ~= size(data.data,4) || size(data.events,2) ~= size(data.data,3)
        err = [err ;29];
        disp('Number of events does not match the data.');
    end;
end;

if isfield(data,'avgNum')
    if size(data.avgNum,1) ~= size(data.data,4) || size(data.avgNum,2) ~= size(data.data,3)
        err = [err ;30];
        disp('Number of trials going into averages does not match the data.');
    end;
end;

if isfield(data,'covNum')
    if size(data.covNum,1) ~= size(data.data,4) || size(data.covNum,2) ~= size(data.data,3)
        err = [err ;116];
        disp('Effective number of trials going into averages does not match the data.');
    end;
end;

if ~isfield(data,'eloc')
    err = [err ;38];
    disp('The .eloc field is missing.');
elseif ~isempty(data.eloc)
    if ~isempty(data.facVecS)
        if (size(data.facVecS,1) ~= length(data.eloc))
            err = [err ;35];
            disp('Number of electrode coordinates does not match the data.');
        end;
    else
        if length(data.eloc) ~= size(data.data,1) && ~isempty(data.data)
            err = [err ;35];
            disp('Number of electrode coordinates does not match the data.');
        end;
    end;
    if isfield(data.eloc,'sph')
        beep
        disp('debugging note:  sph field has appeared in eloc structure.')
    end;
end;

implicitRef=0;
if ~isfield(data,'implicit')
    err = [err ;109];
    disp('The .implicit field is missing.');
elseif ~isempty(data.implicit)
    implicitRef=sum(strcmp('REF',{data.implicit.type}));
end;

if isfield(data,'reference')
    if (length(data.reference.current)+ implicitRef) > 2
        err = [err ;36];
        disp('There cannot be more than two reference channels.');
    end;
end;

if isfield(data,'facVecT')
    if ~isempty(data.facVecS)
        if (size(data.facVecS,1) ~= length(data.chanTypes))
            err = [err ;37];
            disp('Number of channel types does not match the data.');
        end;
    end;
else
    if isfield(data,'chanTypes')
        if length(data.chanTypes) ~= size(data.data,1) && ~isempty(data.data)
            err = [err ;37];
            disp('Number of channel types does not match the data.');
        end;
    end;
end;

if ~isfield(data,'facVecT')
    if isfield(data,'chanTypes')
        if length(data.chanTypes) ~= size(data.data,1) && ~isempty(data.data)
            err = [err ;37];
            disp('Number of channel types does not match the data.');
        end;
    end;
end;

if isfield(data,'facTypes')
    if (length(data.facTypes) ~= (size(data.data,5)+numFacData)) && not (isempty(data.facTypes) && size(data.data,5)==1)
        err = [err ;43];
        disp('Number of factor types does not match the data.');
    end;
end;

if isfield(data,'timeUnits')
    if ~any(strcmp(data.timeUnits,{'ms','per'}))
        err = [err ;123];
        disp('Time units field is not a valid value.');
    end;
end;

if ~isempty(data.dataType)
    switch data.dataType
        case 'continuous'
            if length(data.subNames) ~= 1
                err = [err ;69];
                disp('Continuous data should have only one subject waveform.');
            end;
        case 'single_trial'
            if length(data.trialNames) ~= size(data.data,3)
                err = [err ;23];
                disp('Number of trial names does not match the data.');
            end;
            if ~isempty(data.trialSpecNames)
                if size(data.trialSpecs,1) ~= size(data.data,3)
                    err = [err ;25];
                    disp('Number of trial specs does not match the data.');
                end;
                if size(data.trialSpecs,2) ~= length(data.trialSpecNames)
                    err = [err ;26];
                    disp('Number of trial spec names does not match the trial specs.');
                end;
            end;
            errFlag=0;
            cellNameList=unique(data.cellNames);
            for iCell=1:length(cellNameList)
                theCell=cellNameList{iCell};
                theTrials=data.trialNames(find(strcmp(theCell,data.cellNames)));
                if length(unique(theTrials)) ~= length(theTrials)
                    errFlag=1;
                end;
            end;
            if errFlag
                err = [err ;120];
                disp('Trial names for a given cell are not unique.');
            end;
        case 'average'
            if ~isempty(data.trialNames)
                err = [err ;68];
                disp('Average files should not have trial names.');
            end;
            if ~isempty(data.trialSpecNames)
                if size(data.trialSpecs,1) ~= size(data.data,3)
                    err = [err ;25];
                    disp('Number of trial specs does not match the number of cells.');
                end;
                if size(data.trialSpecs,2) ~= length(data.trialSpecNames)
                    err = [err ;26];
                    disp('Number of trial spec names does not match the trial specs.');
                end;
                if size(data.trialSpecs,3) ~= size(data.data,4)
                    err = [err ;120];
                    disp('Number of trial specs does not match the number of subjects.');
                end;
            end;
        otherwise
            err = [err ;24];
            disp(['Data type ' data.dataType ' is not a recognized type.']);
    end;
end;

if ~isfield(data,'cov')
    err = [err ;113];
    disp('The .cov field is missing.');
else
    if ~isempty(data.cov)
        if (size(data.cov.covMatrix,1) ~= length(data.subNames)) || (size(data.cov.covMatrix,2) ~= length(data.chanNames)) || (size(data.cov.covMatrix,3) ~= length(data.chanNames))
            err = [err ;113];
            disp('The size of the cov field does not match the data.');
        end;
    end;
    if ~isempty(data.cov)
        if (length(data.cov.Nq) ~= length(data.subNames))
            err = [err ;119];
            disp('The size of the cov.Nq field does not match the data.');
        end;
    end;
end;

[EPfieldNames]=ep_fieldNames;
theFieldNames=fieldnames(data);

if ~isequal(EPfieldNames,theFieldNames)
    err = [err ;110];
    disp(['Fields are not correct or at least not in the correct order:']);
    for iField=1:length(EPfieldNames)
        if iField > length(theFieldNames)
            disp([' -' EPfieldNames{iField}]);
        else
            if ~strcmp(theFieldNames{iField},EPfieldNames{iField})
                disp([theFieldNames{iField} '-' EPfieldNames{iField}]);
            end;
        end;
    end;
end;

if isempty(err)
    err=0;
else
    beep;
end;

