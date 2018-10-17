function [eloc,chanNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM,badChans,reference,facVecS,covMatrix2,impedances2]=ep_addEloc(ced,eloc,fileFormat,dataType,chanNames,timeNames,cellNames,subNames,freqNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM,covMatrix,badChans,reference,facVecS,hdr,impedances);
% [eloc,chanNames,chanTypes,sevenDdata,noise,stanDev,badChans,reference,facVecS,covMatrix2,impedances2]=ep_addEloc(ced,eloc,fileFormat,dataType,chanNames,timeNames,cellNames,subNames,freqNames,chanTypes,sevenDdata,noise,stanDev,covMatrix,badChans,reference,facVecS,hdr,impedances) -
% Adds electrode coordinate information to data file.
%
%Input:
%    sevenDdata      : 7D data matrix [channels, time points, cells/trials, subjects, factors, freqs, relations]
%                 Trials are grouped by cell in single_trial format.
%    noise     : 4D matrix [channels, time points, cells/trials, subjects] mirroring .data.
%                 This is the +/- reference average (Schimmel, 1967) which provides an estimate of the noise level
%                 in averaged data by flipping every other trial.  Not applicable to spectral data.
%    stanDev       : 6D matrix [channels, time points, cells/trials, subjects, 1, freqs] mirroring .data.
%                 This is the standard deviation from averaging trials.
%    .stanDevCM     : 6D matrix [channels, time points, cells/trials, subjects, 1, freqs].
%                 This is the Cousineau-Morey std (std of data for a given time point mean-corrected across the cells of a given subject).
%                 This is empty for all but grand average data.
%    covMatrix : 3D matrix [subject, channels, channels] containing original covariance matrix from averaging step.
%                 This is the standard deviation from averaging trials.
%    fileFormat: The original file format.
%    dataType  : The type of the data: 'continuous', 'single_trial', or 'average' (default: average)
%    chanNames : The channel names.
%    timeNames : The msec of the sample onset with respect to the stimulus onset time.
%                 For time-frequency data, the msec of the middle of the .5 sec window.
%    subNames  : The subject names
%    cellNames : The cell names (once for each trial for single_trial files).
%    trialNames: The trial number ID per cell (starting from 1). (single_trial data only)
%    freqNames : The frequency at the middle of the frequency bin.
%    facNames  : The factor names (only for factor files)
%    chanTypes : The type of the channel: EEG, MEG, ECG, ANS (autonomic), REG (regional average)
%                 The CED file can also specify FID (fiducial) and BAD (delete) channel types,
%                 but they will not end up as a chanType.
%    ced       : The name of the .ced file for electrode coordinates.  Can also have the path.
%    eloc      : The electrode location information, one for each channel (see readlocs header)
%                 eloc is the same length as the channels, with REG channels having a blank entry.
%    facVecS   : For spatial PCA factor files, the factor scalp topography.  Used to compress the data.(rows=chans,cols=factors)
%    badChans : Array of corrected bad channels (subject,cell/trial,channel). -1 in a session or continuous file means still bad.
%                   Negative numbers in an average file means number of still bad channels that went into the average
%                   (or rather, were left out of the average).  NaN in an average file means still bad.
%    reference
%        .original    : Original recording reference(s): 1-2 numbers in array
%        .current     : Current reference channel(s): 1-2 numbers in array
%        .type        : Current reference type: REG (regular), AVG (average reference, even if no longer adds to zero), CSD (current source density)
%    hdr      : Header information from original data file, if available.
%    impedances
%        .channels   : impedance values of the channels (chan,subject) (can be empty)
%        .ground     : impedance values of the ground electrode (subject) (can be empty)

%History
%  by Joseph Dien (1/15/14)
%  jdien07@mac.com
%
% bugfix 2/17/14 JD
% Fixed crash when reading ced file with REF channel type indicated.
%
% bugfix 3/2/14 JD
% Fixed REF channel type not being changed to EEG for files with one reference channel.
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% bugfix 3/24/14 JD
% Fixed crash when loading file type that has internal channel names (like Neuroscan) and the only mismatch between it
% and the ced file is a single implicit REF channel or there is no mismatch and there is a single explicit REF channel.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% bugfix 5/28/14 JD
% Fixed crash when the ced file has no type field.
%
% modified 1/7/15 JD
% Handle alternate EGI names for the vertex channel when matching channel names with the ced channel names.
%
% bugfix 5/21/15 JD
% Fixed losing electrode coordinates of eeglab files when type field is empty.
%
% modified 10/9/15 JD
% When none of the CED channel names match but there are the same number of EEG
% channels, it will be assumed that they are the same channels and in the
% same order.  A warning message is provided.
%
% modified 11/22/17 JD
% If none of the channel names match but there are the same number of non-ref EEG channels, will
% not only use the CED channel names, will assume the ref channels are implicit and add them.
% Eliminated restrictions on location of CED files.
% Added support for impedances field.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% bufix 3/26/18 JD
% Fixed crash when adding a channel and there are impedance values.
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

covMatrix2=[];
impedances2.channels=[];
impedances2.ground=[];

if isempty(eloc)
    if strcmp(ced,'internal')
        disp('Electrode information is designated as being provided by the file but no electrode coordinates have been obtained.');
        return;
    end;
    if exist(ced,'file')
        whichCED=ced;
        [pathstr, name, fileSuffix] = fileparts(whichCED);
        ced=[name fileSuffix];
    else
        whichCED=which(ced);
    end;
    try
        disp(['Loading the ced file: ' whichCED]);
        evalc('eloc = readlocs([whichCED],''filetype'',''chanedit'');');
        origEloc=eloc;
        implicit=eloc(1);
        implicit(1)=[];
    catch
        disp(['The ced file ' ced ' did not work for some reason.  The error message was:']);
        disp(lasterr)
        ced = [];
        eloc = [];
        implicit=[];
    end;
else
    implicit=eloc(1);
    implicit(1)=[];
end;

if ~isempty(ced) && strcmp(ced(1:3),'GSN')
    %if is an EGI ced and there is a Cz channel in the eloc but not in the
    %set of channel names, then before adding Cz to the channels, check
    %first to make sure there isn't a VREF or an E129 (for 129 channels),
    %which are alternate EGI names for the vertex channel.  If there is,
    %then change the Cz eloc label to the alternate channel name.
    theCzChan=find(strcmp('Cz',{eloc.labels}));
    if ~isempty(theCzChan)
        if ~any(strcmp('Cz',chanNames))
            vertexChan=find(strcmp('VREF',chanNames));
            if ~isempty(vertexChan)
                eloc(theCzChan).labels=chanNames{vertexChan};
            else
                numElecs=length(find(strcmp('EEG',chanTypes)));
                if ismember(numElecs,[17,33,65,129,257])
                    altName=['E' num2str(numElecs)];
                    vertexChan=find(strcmp(altName,chanNames));
                    if ~isempty(vertexChan) && isempty(find(strcmp(altName,{eloc.labels})))
                        eloc(theCzChan).labels=chanNames{vertexChan};
                    end;
                end;
            end;
        end;
    end;
end;

if isfield(eloc,'type')
    for i=1:length(eloc)
        if ~isempty(eloc(i).labels) && ~isempty(eloc(i).theta) && isempty(eloc(i).type) && ~isnumeric(eloc(i).theta)
            eloc(i).type=eloc(i).theta; %if a CED file has just the label and the type filled out, the type info migrates over to the theta column for some reason.
            eloc(i).theta=[];
        end;
    end;
end;

if ~isempty(eloc)
    if length(eloc) ~= length(unique({eloc.labels}))
        disp(['The ced file ' ced ' had channel names that were not unique.']);
        ced = [];
        eloc = [];
    end
end;

if ~isempty(ced)
    if ~isfield(eloc,'type')
        eloc(1).type=[];
        implicit(1).type=[];
        disp('Type field missing from ced file.  Will assume all channels are EEG channels.');
    end;
    typeFlag=0;
    for i=1:length(eloc)
        if isempty(eloc(i).type)
            eloc(i).type='EEG'; %if type fields are empty, just assume they are EEG channels.
        end;
        if isnumeric(eloc(i).type)
            eloc(i).type='EEG';
            typeFlag=1;
        end;
    end;
    if typeFlag
        disp('Warning: the ''type'' field from the CED file contained numbers rather than labels like EEG or REF.');
    end;
    
    numEEG=length(find(ismember({eloc.type},{'EEG','MGM','MGA','MGP'})));
    numREF=length(find(strcmp('REF',{eloc.type})));
    numFID=length(find(strcmp('FID',{eloc.type})));
    numANS=length(find(ismember({eloc.type},{'ANS','ECG'})));
    numBAD=length(find(strcmp('BAD',{eloc.type})));
    numOther=length(eloc)-numEEG-numREF-numFID-numANS-numBAD;
    if numOther > 0
        disp('Warning: Unrecognized channel type in CED file, will be deleted:');
        disp(unique(setdiff({eloc.type},{'EEG','MGM','MGA','MGP','REF','FID','ANS','ECG','BAD'})))
        otherInd=find(ismember({eloc.type},setdiff({eloc.type},{'EEG','MGM','MGA','MGP','REF','FID','ANS','ECG','BAD'})));
        for i=1:length(otherInd)
            eloc(otherInd(i)).type='BAD';
        end;
        numBAD=numBAD+length(otherInd);
    end;
    
    if numFID > 0
        FIDs=find(strcmp('FID',{eloc.type}));
        temp=eloc(FIDs);
        [temp.type] = deal('FID');
        implicit(end+1:end+numFID) = temp;
        eloc=eloc(setdiff([1:length(eloc)],FIDs));
    end;
    
    %Formats that use a fixed channel order and have no channel labels in the header
    if any(strcmp(fileFormat,{'egi_egia','egi_egis','egi_sbin','text','ns_mat'}))
        nonFID=find(ismember({eloc.type},{'EEG','MGM','MGA','MGP','REF','ANS','ECG','BAD'}));
        temp={eloc.labels};
        if (length(chanNames) ~= length(nonFID)) && (length(chanNames)+1 ~= length(nonFID))
            disp(['Error: This is a file format where the number of channels should be fixed and yet the number of channels (' num2str(length(chanNames)) ') differs from the number of EEG channels in the CED file (' num2str(length(nonFID)) ').']);
            disp('Giving up on electrode coordinates.');
            eloc = [];
            return
        end;
        chanNames=temp(nonFID(1:length(chanNames)));
    end;
    
    %Formats that include channel labels in the header and where channels have flexible ordering or may even be absent
    
    %For such files channels that were left out due to being bad data or due to being an implicit reference
    %but are known to exist from the CED file will be added back to the data file.  Missing data channels
    %will be marked bad and be set to zero.  An implicit reference channel will be set to zero but not be marked
    %bad.  If there are two references (one implicit and one explicit) then it will be assumed to be a mean
    %mastoid reference or equivalent and the implicit one will be set to be the negative of the explicit one (as
    %is appropriate for such a reference scheme).  If there are two reference channels (according to the CED)
    %and both are missing (as is sometimes done) then both channels will be set to be bad as their contents
    %cannot be determined.  Fiducial channels will always be set to be implicit and there will be no
    %corresponding channel in the voltage data.
    
    %      else
    chanCount=0; %channels that are in the data and also in the CED
    nonBadChanCount=0; %channels that are in the data and also in the CED and are not BAD
    nonCED=[]; %channels that are in the data but not in the CED
    elocLabels={eloc.labels};
    elocTypes={eloc.type};
    extraDataNames=cell(0);
    for i=1:length(chanNames)
        if find(strcmpi(chanNames(i),elocLabels(find(ismember(elocTypes,{'EEG','MGM','MGA','MGP','REF','ANS','ECG'})))))
            chanCount=chanCount+1;
            nonBadChanCount=nonBadChanCount+1;
        elseif find(strcmpi(chanNames(i),elocLabels(find(strcmp('BAD',elocTypes)))))
            chanCount=chanCount+1;
        else
            nonCED(end+1)=i;
            extraDataNames{end+1}=chanNames{i};
        end;
    end;
    
    if (chanCount ==0) && (length(chanNames) == length(elocLabels(find(ismember(elocTypes,{'EEG','REF'})))))
        disp(['None of the channel names in the ced file ' ced ' match those in the file.']);
        disp('But the number of EEG channels are the same.  Will make the assumption that they are the same channels and in the same order.')
        disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
        disp('Will keep the channel names in the CED file.');
        chanCount=length(chanNames);
        chanNames=elocLabels(find(ismember(elocTypes,{'EEG','REF'})));
        nonCED=[];
        extraDataNames=cell(0);
    elseif (chanCount ==0) && (length(chanNames) == length(elocLabels(find(ismember(elocTypes,'EEG')))))
        disp(['None of the channel names in the ced file ' ced ' match those in the file.']);
        disp('But the number of EEG channels are the same as the non-reference channels.  Will make the assumption that they are the same channels and in the same order.')
        disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
        disp('Will keep the channel names in the CED file and also add the missing reference channels.');
        chanCount=length(chanNames);
        chanNames=elocLabels(find(ismember(elocTypes,'EEG')));
        nonCED=[];
        extraDataNames=cell(0);        
    elseif (chanCount ==0) && (length(chanNames) == length(elocLabels(find(ismember(elocTypes,{'EEG','BAD'})))))
        disp(['None of the channel names in the ced file ' ced ' match those in the file.']);
        disp('But the number of EEG and BAD channels are the same as the non-reference channels.  Will make the assumption that they are the same channels and in the same order.')
        disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
        disp('Will keep the channel names in the CED file and also add the missing reference channels.');
        chanCount=length(chanNames);
        chanNames=elocLabels(find(ismember(elocTypes,'EEG')));
        nonCED=[];
        extraDataNames=cell(0);        
    elseif (nonBadChanCount ==0) && (length(chanNames) == length(elocLabels(find(ismember(elocTypes,{'EEG','BAD'})))))
        disp(['None of the non-BAD channel names in the ced file ' ced ' match those in the file.']);
        disp('But the number of EEG channels are the same as the non-reference channels.  Will make the assumption that they are the same channels and in the same order.')
        disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
        disp('Will keep the channel names in the CED file and also add the missing reference channels.');
        chanCount=length(chanNames);
        chanNames=elocLabels(find(ismember(elocTypes,'EEG')));
        nonCED=[];
        extraDataNames=cell(0);        
    end;
    
    if chanCount > 0
        extraCED=0; %number of channels only in the CED file
        extraCEDNames=cell(0);
        chanIndex=find(ismember(elocTypes,{'EEG','MGM','MGA','MGP','REF','ANS','ECG','BAD'}));
        for newElocChan=1:length(chanIndex)
            oldElocChan=chanIndex(newElocChan);
            oldFileChan=find(strcmpi(eloc(oldElocChan).labels,chanNames));
            if isempty(oldFileChan)
                extraCED=extraCED+1;
                extraCEDNames{end+1}=eloc(oldElocChan).labels;
            end;
        end;
        %Set up new data arrays with the full set of channels
        numChan=length(chanNames)+extraCED;
        numPoints=length(timeNames);
        numCellTrials=length(cellNames);
        numFacs=size(sevenDdata,5);
        numFreqs=length(freqNames);
        numSubs=length(subNames);
        numRels=size(sevenDdata,7);
        if isempty(badChans)
            badChans=zeros(numSubs,numCellTrials,numChan);
        end;
        if ~isempty(facVecS)
            facVecS2=zeros(numChan,numFacs);
        else
            sevenDdata2=zeros(numChan,max(1,numPoints),numCellTrials,numSubs,numFacs,max(1,numFreqs),numRels);
        end;
        if ~isempty(noise)
            noise2=zeros(numChan,max(1,numPoints),numCellTrials,numSubs,numFacs);
        end;
        if ~isempty(stanDev)
            stanDev2=zeros(numChan,max(1,numPoints),numCellTrials,numSubs,numFacs,max(1,numFreqs));
        end;
        if ~isempty(stanDevCM)
            stanDevCM2=zeros(numChan,max(1,numPoints),numCellTrials,numSubs,numFacs,max(1,numFreqs));
        end;
        if ~isempty(impedances.channels)
            impedances2.channels=nan(numChan,numSubs);
        end;
        if ~isempty(covMatrix)
            covMatrix2=zeros(numSubs,numChan,numChan);
        end;
        chanNames2=cell(numChan,1);
        chanTypes2=cell(numChan,1);
        if ~isempty(badChans)
            badChans2=ones(numSubs,size(badChans,2),numChan); %left out channels will be marked bad
            if ~strcmp(dataType,{'single_trial','continuous'})
                badChans2=badChans2*(-1); %if an average file then change to -1 as the marker for a bad channel
            end;
        end;
        for newElocChan=1:length(chanIndex)
            oldElocChan=chanIndex(newElocChan);
            oldFileChan=find(strcmpi(eloc(oldElocChan).labels,chanNames));
            if ~isempty(oldFileChan) %channel is in both CED and data
                eloc2(newElocChan)=eloc(oldElocChan); %generate new full eloc array
                if ~isempty(facVecS)
                    facVecS2(newElocChan,:)=facVecS(oldFileChan,:);
                else
                    sevenDdata2(newElocChan,:,:,:,:,:,:)=sevenDdata(oldFileChan,:,:,:,:,:,:);
                end;
                if ~isempty(noise)
                    noise2(newElocChan,:,:,:,:)=noise(oldFileChan,:,:,:,:);
                end;
                if ~isempty(stanDev)
                    stanDev2(newElocChan,:,:,:,:,:)=stanDev(oldFileChan,:,:,:,:,:);
                end;
                if ~isempty(stanDevCM)
                    stanDevCM2(newElocChan,:,:,:,:,:)=stanDevCM(oldFileChan,:,:,:,:,:);
                end;
                if ~isempty(covMatrix)
                    covMatrix2(:,newElocChan,newElocChan)=covMatrix(:,oldFileChan,oldFileChan);
                end;
                if ~isempty(impedances.channels)
                    impedances2.channels(newElocChan,:)=impedances.channels(oldFileChan);
                end;
                chanNames2{newElocChan}=eloc(oldElocChan).labels;
                chanTypes2{newElocChan}=eloc(oldElocChan).type;
                if ~isempty(badChans)
                    badChans2(:,:,newElocChan)=badChans(:,:,oldFileChan);
                end;
                if strcmp(chanTypes2{newElocChan},'REF')
                    badChans2(:,:,newElocChan)=0; %implicit reference channels are not bad
                    %chanTypes2{newElocChan}='EEG';
                end;
            else %channel is only in the CED
                eloc2(newElocChan)=eloc(oldElocChan); %generate new full eloc array
                chanNames2{newElocChan}=eloc(oldElocChan).labels;
                chanTypes2{newElocChan}=eloc(oldElocChan).type;
            end;
        end;
        for i=1:length(nonCED) %channel is only in the data
            theNonCED=nonCED(i);
            eloc2(length(chanIndex)+i).labels=extraDataNames{i};
            eloc2(length(chanIndex)+i).theta=[]; %empty theta is marker for missing electrode coordinates
            if ~isempty(facVecS)
                facVecS2(length(chanIndex)+i,:)=facVecS(theNonCED,:);
            else
                sevenDdata2(length(chanIndex)+i,:,:,:,:,:,:)=sevenDdata(theNonCED,:,:,:,:,:,:);
            end;
            if ~isempty(noise)
                noise2(length(chanIndex)+i,:,:,:,:)=noise(theNonCED,:,:,:,:);
            end;
            if ~isempty(stanDev)
                stanDev2(length(chanIndex)+i,:,:,:,:,:)=stanDev(theNonCED,:,:,:,:,:);
            end;
            if ~isempty(stanDevCM)
                stanDev2(length(chanIndex)+i,:,:,:,:,:)=stanDevCM(theNonCED,:,:,:,:,:);
            end;
            if ~isempty(covMatrix)
                covMatrix2(:,length(chanIndex)+i,length(chanIndex)+i)=covMatrix(:,theNonCED,theNonCED);
            end;
            if ~isempty(impedances.channels)
                impedances2.channels(length(chanIndex)+i,:)=impedances.channels(theNonCED);
            end;
            chanNames2(length(chanIndex)+i)=chanNames(theNonCED);
            chanTypes2(length(chanIndex)+i)=chanTypes(theNonCED);
            if strcmp(fileFormat,'egi_mff_v1') && ~isempty(hdr)
                if strcmp(hdr.chantype(theNonCED),'ecg')
                    chanTypes2{length(chanIndex)+i}='ECG';
                end;
            end;
            if ~isempty(badChans)
                badChans2(:,:,length(chanIndex)+i)=badChans(:,:,theNonCED);
            end;
        end;
        
        %handle reference channels
        if numREF > 2
            disp('Warning: File has more than two reference channels marked in CED file.  Reference channel information ignored.');
        end;
        
        if numREF < 3
            impRefs=setdiff({eloc(find(strcmp('REF',{eloc.type}))).labels},chanNames); %implicit refs
            impRefNums=zeros(length(impRefs),1);
            for i=1:length(impRefs)
                impRefNums(i)=find(strcmpi(impRefs{i},{eloc.labels}));
            end;
            numImpRefs=length(impRefs);
            if (numREF==2) && (length(impRefs) == 1) %if mean mastoids with one explicit and one implicit
                impMastoid=find(strcmpi(eloc(impRefNums).labels,chanNames2));
                expMastoid=find(strcmpi(eloc(setdiff(find(strcmp('REF',{eloc.type})),impRefNums)).labels,chanNames2));
                if ~isempty(facVecS)
                    facVecS2(impMastoid,:)=-facVecS(expMastoid,:); %mean mastoid refs will be mirror images
                else
                    sevenDdata2(impMastoid,:,:,:,:,:,:)=-sevenDdata(expMastoid,:,:,:,:,:,:); %mean mastoid refs will be mirror images
                end;
                if ~isempty(noise)
                    noise2(impMastoid,:,:,:,:)=noise(expMastoid,:,:,:,:);
                end;
                if ~isempty(stanDev)
                    stanDev2(impMastoid,:,:,:,:,:)=stanDev(expMastoid,:,:,:,:,:);
                end;
                if ~isempty(stanDevCM)
                    stanDev2(impMastoid,:,:,:,:,:)=stanDevCM(expMastoid,:,:,:,:,:);
                end;
                if ~isempty(impedances.channels)
                    impedances2.channels(impMastoid,:)=impedances.reference(1);
                end;
                if ~isempty(covMatrix)
                    covMatrix2(:,impMastoid,impMastoid)=covMatrix(:,expMastoid,expMastoid);
                end;
                badChans2(:,:,impMastoid)=0; %implicit reference channels are not bad
                reference.original=[impMastoid expMastoid];
                reference.current=reference.original;
                chanTypes2{impMastoid}='EEG';
                chanTypes2{expMastoid}='EEG';
                disp('CED file indicates was originally mean mastoids reference and will assume is still so.');
            elseif (numREF==2) && (length(impRefs) == 2) %if mean mastoids and both were implicit
                % both ref channels will be left as being marked bad as their voltage values will be missing
                M1=find(strcmpi(eloc2(impRefNums).labels,chanNames2));
                M2=find(strcmpi(eloc2(setdiff(find(strcmp('REF',{eloc2.type})),M1)).labels,chanNames2));
                reference.original=[M1 M2];
                reference.current=reference.original;
                chanTypes2{M1}='EEG';
                chanTypes2{M2}='EEG';
                if ~isempty(impedances.channels)
                    impedances2.channels(M1,:)=impedances.reference(1);
                    impedances2.channels(M2,:)=impedances.reference(2);
                end;
                disp('CED file indicates was originally mean mastoids reference and will assume is still so.');
            elseif (numREF==2) && (length(impRefs) == 0) %if mean mastoids and both were explicit
                references=find(strcmp('REF',chanTypes2));
                M1=references(1);
                M2=references(2);
                reference.original=[M1 M2];
                reference.current=reference.original;
                chanTypes2{M1}='EEG';
                chanTypes2{M2}='EEG';
                disp('CED file indicates was originally physically linked mastoids reference and will assume is still so (but incorrect if was rereferenced to mean mastoids after data collection).');
            elseif (numREF==1) && (length(impRefs) == 1) %if a single implicit reference
                M1=find(strcmpi(eloc2(impRefNums).labels,chanNames2));
                reference.original=M1;
                reference.current=reference.original; %assume that if there are implicit references then they were the original reference and still are
                chanTypes2{M1}='EEG';
                badChans2(:,:,M1)=0; %implicit reference channels are not bad
                if ~isempty(impedances.channels)
                    impedances2.channels(M1,:)=impedances.reference(1);
                end;
            elseif (numREF==1) && isempty(impRefs) %if a single explicit reference
                theRef=find(strcmpi(eloc2(find(strcmp('REF',{eloc2.type}))).labels,chanNames2));
                reference.original=theRef;
                
                EEGchans=find(ismember(chanTypes2,{'EEG','REF'}));
                if ~squeeze(any(any(any(any(any(sevenDdata2(theRef,:,:,:,:,:)))))))
                    reference.current=reference.original;
                elseif ~squeeze(any(any(any(any(sum(sevenDdata2(EEGchans,:,:,:,:,:,:),1))))))
                    reference.current=[];
                    reference.type='AVG';
                else
                    reference.current=[];
                    reference.type='REG';
                end;
                chanTypes2{theRef}='EEG';
            elseif numREF==0
                disp('Note: No reference channel was found.');
                EEGchans=find(ismember(chanTypes2,{'EEG','REF'}));
                if ~squeeze(any(any(any(any(sum(sevenDdata2(EEGchans,:,:,:,:,:,:),1))))))
                    reference.current=[];
                    reference.type='AVG';
                end;
            else
                disp('Warning: EP Toolkit was unable to match up contents of the CED file to the data file.');
            end;
        else
            impRefs=[];
            impRefNums=[];
            numImpRefs=0;
        end;
        
        if ~isempty(facVecS)
            facVecS=facVecS2;
        else
            sevenDdata=sevenDdata2;
        end;
        if ~isempty(noise)
            noise=noise2;
        end;
        if ~isempty(stanDev)
            stanDev=stanDev2;
        end;
        if ~isempty(stanDevCM)
            stanDevCM=stanDevCM2;
        end;
        chanNames=chanNames2;
        chanTypes=chanTypes2;
        if ~isempty(badChans)
            badChans=badChans2;
        end;
        
        if ~isempty(nonCED)
            disp(['The following data file channels are not represented in the ced file:']);
            for i=1:length(extraDataNames)
                disp(extraDataNames{i});
            end;
        end;
        if extraCED > 0
            disp(['The following ced file channels are not represented in the data file:']);
            for i=1:length(extraCEDNames)
                disp(extraCEDNames{i});
            end;
        end;
    else
        disp(['None of the channel names in the ced file ' ced ' match those in the file']);
        disp('and the number of channels in the ced file and in the data are different.');
        disp('You need to change the channel names in the CED file or in the data file so they match up.');
        ced = [];
        eloc2 = [];
    end;
end;

if exist('eloc2','var')
    eloc=eloc2;
else
    eloc=[];
end;
