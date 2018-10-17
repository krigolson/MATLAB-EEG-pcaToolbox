function [data] = ep_PCAoutput(FactorResults, params, paramNames);
% ep_PCAoutput - [data] = ep_PCAoutput(FactorResults, params, paramNames); - output PCA results in EP format
%Inputs
%  FactorResults: Structured variable with results and information about analysis.
%    See doPCA and doPCAst documentation.
%
%  -----Optional------
%  params    : parameters to correlate against factor results to produce parametric PCA solutions (weights, parameter).
%              The rows should correspond to input cells and subjects (e.g., [C1S1; C1S2; C2S1; C2S2]).
%              Rows with a value of -999 will be dropped.
%  paramNames: Cell array of parameter names.
%
%Outputs
%  data  : Structured array with the data and accompanying information.  See readData.
%
%  Outputs EP data variable with the factor results.  The factor scores and the factor pattern are multiplied to produce
%  the portion of the grand average accounted for by the factoring.
%  For a two-step PCA, there are the number of the first-stage factors multiplied by the second stage factors.
%  The second stage factors vary most quickly (e.g., for temporospatial PCA, T1S1 T1S2 T1S3...T2S1...etc.).
%  The cells input array allows the output
%  cells to be any desired combination of input cells.  In addition, a grand average of all the input cells is
%  generated and added to the end of the file.  The observations are the factors.
%  For average data, an additional subject is added to the end that is the total amount of the grand average accounted
%  for by the factor analysis.  Also, a cell is added that is the average of all the cells and is labeled "all".
%  Optionally, may specify parameters that will be used for parametric PCA results,
%  which will be output as additional cells.  Parametric cells contain the microvolt change corresponding to a one standard
%  deviation change in the parameter.
%  When using parametric PCA, also saves a .mat file with "Scores" appended to the name containing the
%  scores matrix with factors*chans as cols and obs as rows to facilitate analysis.
%  For when outputting separate subject files, parametric correlations calculated separately for each subject.

%History
%  by Joseph Dien (10/1/00)
%  jdien07@mac.com
%
%  modified (1/9/01) JD
%  Added parametric PCA option.
%
%  bugfix(1/10/01) JD
%  Transpose spatial data blocks before writing out to disk.  Also, changed data input to varSD.
%
%  modified (2/27/01) JD
%  Added peakLatency output.
%
%  modified (5/29/01) JD
%  Added peakChan and peakSamp output.
%
%  corrected (8/22/01) JD
%  Gave added to factor reconstructions to assist identification.
%
%  modified (2/2/02) JD
%  Eliminated assumption that cell numbers of dataset are in ascending order.
%
%  modified (6/7/02) JD
%  Also, modified to eliminate warning message when no exclusion channels.
%
%  modified (6/26/02) JD
%  added option to output matlab matrices rather than EGIS file.
%
%  corrected (6/30/02) JD
%  Fixed subject specs to include subject ID number so NetStation will label observations properly and NAvg so BESA will have N= info.
%
%  bugfix (7/9/02) JD
%  Fixed bug in parametric PCA output.  Results were not being calculated properly.
%
%  upgraded (8/22/02) JD
%  Changed to work under Matlab 6.5
%
%  modified (11/3/02) JD
%  Missing number for parameter vectors is now 999.
%
%  bugfix (3/24/03) JD
%  Eliminated x10 scaling in parametric PCA output.
%
%  bugfix & modified (4/25/03) JD
%  Parametric output after first parametric cell was computed incorrectly.  Also,
%  will change bins to 500 for EGI format output since losing too much resolution
%  when bins equal 1 (as with NS output) since EGIS files are saved in integer format.
%  Also now saves scores matrix when using parametric option.
%
%  modified (5/11/03) JD
%  Added individual subject file output option.
%
% bugfix (5/1/04) JD
% Fixed naming of scores .mat file
%
% bugfix (1/6/05) JD
% Fixed fclose statement, individual file output option, and Mat output option.  Thanks to Robert Frank!
%
% modified (6/16/05) JD
% Added error check for cell numbers.
%
% modified (11/15/05) JD
% In parametric PCA procedure, drop cells corresponding to missing data
% parameters when generating the mean voltage map.
%
% bugfix (12/9/07) JD
% Fixed scaling error when input data is not 500 bins per microvolt
%
% modified (12/19/07) JD
% Changed input cell names and parameter names arrays to cell arrays.
% Changed cellcoll to cell array, thus eliminating need to pad it out with zeroes.
% Deleted old filetype function (only worked on OS 9 Macs).
% Eliminated indexdata.  Eliminated the +gave option.
%
% modified (1/27/08) JD
% EGIS files always written as big-endian since NetStation makes this assumption.
%
% modified (2/21/08) JD
% Copyright notice appears only once per session.  Input factor results are now packaged in a
% structured variable.  Accommodates two-step PCA results by itself.  Subjects mode fixed.
%
% bugfix (6/2/08) JD
% Fixed sign of parametric output.  For example, positive correlation
% should result in positive voltage for positive channels but negative for
% negative channels.  Instead it was producing positive voltages at both
% types of channels.  Sign is now determined solely by correlation.
%
% modified (9/9/08) JD
% Automatically sets File Type and Creator Type if on an OS X Mac.
%
% modified (11/5/08) JD
% Automatically sets montage type if on an OS X Mac.
%
% modified (2/24/09) JD
% Gets baseline, experiment name, and sampling rate from structured variable.  Deleted cal field.  Assumes bins equals 1.
%
%  modified (11/3/09) JD
%  Missing number for parameter vectors changed from 999 to -999 (less likely to have problems with freq parameter).
%
% modified (3/11/09) JD
% Changed name to PCAoutput.m and changed from saving EGIS file to outputing EP data structure.
% Made the input fields, other than the data structure, optional.
%
% bugfix (4/25/09) JD
% Subs option had been unintentionally eliminated by last modification.  Restored the option.
%
% modified 5/10/09 JD
% Added cellType and subType and facNames and facTypes fields.  Eliminated subs option in favor of a single EP file
% with both subject average and grand average information in a 5D data field.
%
% bugfix 7/17/09 JD
% cellTypes correctly set to be a vector.
%
% modified 6/26/09 JD
% Moved facName and facType to doPCA and doPCAst.  Addition of .pca field moved to this function.
% Data represented by factor loadings no longer expanded into full .data field in order to avoid memory limitations.
% Instead, factor loadings are only expanded as needed by ep_expFac function, using facVecS and facVecT.
% Added dataName and facData fields.
%
% bugfix 8/3/09 JD
% crash when doing PCA of single trial data fixed.
%
% bugfix 8/26/09 JD
% Made sumFacs more memory efficient to reduce number of out of memory errors.
%
%  modified 2/15/10 JD
%  Analysis fields no longer optional.
%
% bugfix 3/14/10 JD
% ensure that summary factor, subject, and cell names are not already taken.  If so, add a unique suffix to the name.
%
% bugfix 6/16/10 JD
% avoid adding grand average to subjects dimension if the data type is continuous.
%
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
%
% modified 1/24/12 JD
% Added reference field to keep track of original and current reference scheme.
%
% modified 3/21/12 JD
% Eliminated peak chan and latency and sample output as no longer needed.
%
% modified 4/11/12 JD
% Added support for 6th dimension of frequency for time-frequency analyses.  Unlogs frequency data for computations.
%
% bugfix 10/19/12 JD
% Fixed facnames and factypes not necessarily being column vectors.
%
% bugfix 1/21/13 JD
% Fixed error message when applying PCA to data with only one cell.
%
% bugfix 1/25/13 JD
% Fixed bug in reconstruction of PCA data from frequency-domain PCAs.
%
% bugfix 1/31/13 JD
% Fixed error message when generating PCA data from a single subject average.
%
% bugfix 9/26/13 JD
% Fixed PCA of continuous data generating error message.
%
% modified 10/9/13 JD
% Added recTime field.
%
% modified 2/27/14 JD
% Fields output in standard order.
%
% bufix 3/12/14 JD
% Handles decimal sampling rates gracefully.
%
% bufix 3/19/14 JD
% Fixed recTime field not including space for the 'all' cell in PCA output, resulting in crashes when edited.
%
% modified 3/23/14 JD
% Average numbers, trial specs, and events carried over to the PCA file.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% modified 5/1/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure, including complex numbers.
% Added ability to perform PCA on datasets with bad data.
% FacScr observations are now always arranged with permutations in order of the seven data dimensions.
% Dropped cellcoll and exclChan parameters.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% modified 10/16/16 JD
% Added .stims field.
%
% modified 11/14/16 JD
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

% peakLatency=[];
% peakSamp=[];
% peakChan=[];
data=[];

if ~isfield(FactorResults,'varSD')
    error('No varSD field.');
end
varSD=FactorResults.varSD;
if ~isfield(FactorResults,'numchan')
    error('No numchan field.');
end
numChans=FactorResults.numchan;
if ~isfield(FactorResults,'timepoints')
    error('No timepoints field.');
end
numPoints=FactorResults.timepoints;
if ~isfield(FactorResults,'numFreqs')
    error('No freqs field.');
end
numFreqs=FactorResults.numFreqs;
if ~isfield(FactorResults,'numSubs')
    error('No numSubs field.');
end
numSubs=FactorResults.numSubs;

if ~isfield(FactorResults,'numCells')
    error('No numCells field.');
end
numCells=FactorResults.numCells;

if ~isfield(FactorResults,'numRels')
    error('No numRels field.');
end
numRels=FactorResults.numRels;

if ~isfield(FactorResults,'montage')
    montage=ep_askForMontage;
else
    montage=FactorResults.montage;
end
if ~isfield(FactorResults,'Fs')
    error('No Fs field.');
end
SampleRate=FactorResults.Fs;
if ~isfield(FactorResults,'baseline')
    error('No baseline field.');
end
base=FactorResults.baseline*(1000/SampleRate);

if ~isfield(FactorResults,'ename')
    error('No ename field.');
end

if nargin < 2
    params = [];
    paramNames = [];
end

cnames = FactorResults.cellNames;

numOutputCells=length(cnames);

if length(paramNames) ~= size(params,2)
    error('Error - The number of Parameter Cell Names does not match the number of parameter output cell specifications');
end;

if isfield(FactorResults,'PCAmode3')  %if this is the result of a three-step PCA
    numfacs3 = FactorResults.numFacs3;	%number of factors retained in third PCA step
    numfacs2 = FactorResults.numFacs2;	%number of factors retained in second PCA step
    numfacs1 = FactorResults.numFacs;	%number of factors retained in first PCA step
    
    disp(['Generating PCA file with ' num2str(numfacs1) ' factors from the first PCA step and ' num2str(numfacs2) ' factors in the second and ' num2str(numfacs3) ' factors in the third.']);
    
    tempFacPatST=zeros(size(FactorResults.FacPatST,1),numfacs1*numfacs2*numfacs3);
    tempFacPat=zeros(size(FactorResults.FacPat,1),numfacs1*numfacs2*numfacs3);
    for i = 1:numfacs1
        for j = 1:numfacs2
            for k = 1:numfacs3
                tempFacPat3(:,(i-1)*numfacs2*numfacs3+(j-1)*numfacs3+k)=diag(FactorResults.varSD3((i-1)*numfacs2+j,:))*(FactorResults.FacPat3(:,(i-1)*numfacs2*numfacs3+(j-1)*numfacs3+k));
                tempFacPatST(:,(i-1)*numfacs2*numfacs3+(j-1)*numfacs3+k)=diag(FactorResults.varSDST(i,:))*(FactorResults.FacPatST(:,(i-1)*numfacs2+j));
                tempFacPat(:,(i-1)*numfacs2*numfacs3+(j-1)*numfacs3+k)=diag(FactorResults.varSD)*FactorResults.FacPat(:,i);
            end;
        end;
    end;
    facVecT=[];
    facVecS=[];
    facVecF=[];
    switch FactorResults.PCAmode3
        case 'temp'
            facVecT=tempFacPat3;
        case 'spat'
            facVecS=tempFacPat3;
        case 'freq'
            facVecF=tempFacPat3;
    end
    switch FactorResults.PCAmode2
        case 'temp'
            facVecT=tempFacPatST;
        case 'spat'
            facVecS=tempFacPatST;
        case 'freq'
            facVecF=tempFacPatST;
    end
    switch FactorResults.PCAmode
        case 'temp'
            facVecT=tempFacPat;
        case 'spat'
            facVecS=tempFacPat;
        case 'freq'
            facVecF=tempFacPat;
    end
    
    FacScr=FactorResults.FacScr3;
    if isfield(FactorResults,'badObs3')
        badObs=FactorResults.badObs3;
    else
        msg{1}='Error: This PCA data was generated prior to changes made in EP 2.44 so it cannot be used.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    facNames=FactorResults.facNames3;
    facTypes=FactorResults.facTypes3;
    facVar=FactorResults.facVar3;
    facVarQ=FactorResults.facVarQ3;
    
    ObsNum=1; %the number of channels or timepoints not accounted for by factor loadings and thus accounted for by factor scores.
    dataChanNum=1; %the number of channels represented by the factor scores
    dataTimeNum=1; %the number of time points represented by the factor scores
    dataFreqNum=1; %the number of frequency bins represented by the factor scores
elseif isfield(FactorResults,'PCAmode2')  %if this is the result of a two-step PCA
    varSDST=FactorResults.varSDST;
    PCAmode=FactorResults.PCAmode2;
    numfacs2 = FactorResults.numFacs2;	%number of factors retained in second PCA step
    numfacs1 = FactorResults.numFacs;	%number of factors retained in first PCA step
    
    disp(['Generating PCA file with ' num2str(numfacs1) ' factors from the first PCA step and ' num2str(numfacs2) ' factors in the second.']);
    
    tempFacPatST=zeros(size(FactorResults.FacPatST,1),numfacs1*numfacs2);
    tempFacPat=zeros(size(FactorResults.FacPat,1),numfacs1*numfacs2);
    for i = 1:numfacs1
        for k = 1:numfacs2
            tempFacPatST(:,(i-1)*numfacs2+k)=diag(FactorResults.varSDST(i,:))*(FactorResults.FacPatST(:,(i-1)*numfacs2+k));
            tempFacPat(:,(i-1)*numfacs2+k)=diag(varSD)*FactorResults.FacPat(:,i);
        end;
    end;
    facVecT=[];
    facVecS=[];
    facVecF=[];
    dataChanNum=numChans;
    dataTimeNum=numPoints;
    dataFreqNum=numFreqs;
    dataRelNum=numRels;
    switch FactorResults.PCAmode2
        case 'temp'
            facVecT=tempFacPatST;
            dataTimeNum=1;
        case 'spat'
            facVecS=tempFacPatST;
            dataChanNum=1;
        case 'freq'
            facVecF=tempFacPatST;
            dataFreqNum=1;
    end
    switch FactorResults.PCAmode
        case 'temp'
            facVecT=tempFacPat;
            dataTimeNum=1;
        case 'spat'
            facVecS=tempFacPat;
            dataChanNum=1;
        case 'freq'
            facVecF=tempFacPat;
            dataFreqNum=1;
    end
    ObsNum=dataTimeNum*dataChanNum*dataFreqNum*dataRelNum; %the number of channels or timepoints or frequency bins not accounted for by factor loadings and thus accounted for by factor scores.
    
    FacScr=FactorResults.FacScrST;
    if isfield(FactorResults,'badObsST')
        badObs=FactorResults.badObsST;
    else
        msg{1}='Error: This PCA data was generated prior to changes made in EP 2.44 so it cannot be used.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    facNames=FactorResults.facNames2;
    facTypes=FactorResults.facTypes2;
    facVar=FactorResults.facVarST;
    facVarQ=FactorResults.facVarQST;
    
elseif isfield(FactorResults,'PCAmode')
    PCAmode=FactorResults.PCAmode;
    if ~isfield(FactorResults,'FacScr')
        error('No FacScr field.');
    end
    FacScr=FactorResults.FacScr;
    if isfield(FactorResults,'badObs')
        badObs=FactorResults.badObs;
    else
        msg{1}='Error: This PCA data was generated prior to changes made in EP 2.44 so it cannot be used.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    if ~isfield(FactorResults,'FacPat')
        error('No FacPat field.');
    end
    
    facVecS=[];
    facVecT=[];
    facVecF=[];
    switch PCAmode
        case 'temp'
            facVecT=diag(varSD)*FactorResults.FacPat;
            dataChanNum=numChans;
            dataRelNum=numRels;
            dataTimeNum=1;
            dataFreqNum=numFreqs;
        case 'spat'
            facVecS=diag(varSD)*FactorResults.FacPat;
            dataChanNum=1;
            dataRelNum=numRels;
            dataTimeNum=numPoints;
            dataFreqNum=numFreqs;
        case 'freq'
            facVecF=diag(varSD)*FactorResults.FacPat;
            dataChanNum=numChans;
            dataRelNum=numRels;
            dataTimeNum=numPoints;
            dataFreqNum=1;
        otherwise
            error([PCAmode ' is not a recognized PCAmode.']);
    end;
    
    numfacs1 = FactorResults.numFacs;	%number of factors retained in first PCA step
    disp(['Generating PCA file with ' num2str(numfacs1) ' factors.']);
    
    facNames=FactorResults.facNames;
    facTypes=FactorResults.facTypes;
    facVar=FactorResults.facVar;
    facVarQ=FactorResults.facVarQ;
    
    ObsNum=dataTimeNum*dataChanNum*dataFreqNum*dataRelNum; %the number of channels or timepoints or frequency bins not accounted for by factor loadings and thus accounted for by factor scores.
    
else
    msg{1}='Not a PCA analysis file.';
    [msg]=ep_errorMsg(msg);
    return
end;

if mod(size(FacScr,1),(ObsNum * numSubs)) ~= 0
    error('Error - The number of factor scores is not evenly divided by the number of subjects and the channels/time points/frequency bins.');
end;

NumCell = size(FacScr,1)/(ObsNum * numSubs);	%Find number of cells

if ~isempty(params)
    if size(params,1) ~= (size(FacScr,1)/ObsNum)
        error('Error - number of weights for each parameter needs to be equal to the number of subjects times cells');
    end;
end;

numFacs = size(FacScr,2);
theFacScr=permute(reshape(FacScr,dataChanNum,dataTimeNum,NumCell,numSubs,dataFreqNum,dataRelNum,numFacs),[1 2 3 4 7 5 6]);
badData=reshape(full(badObs(:,1)),dataChanNum,dataTimeNum,NumCell,numSubs,1,dataFreqNum,dataRelNum); %each column of badObs should be identical

cellTypes=cell(numOutputCells,1);
for i = 1:numOutputCells
    cellTypes{i}='SGL';
end;

fprintf('%60s\n',' ' );

sevenDdata=zeros(dataChanNum,dataTimeNum,numOutputCells,numSubs,numFacs,dataFreqNum,dataRelNum);
%Set up the data output matrix populated by factor score information.  The factor loading information will remain in the facVec matrices.
% for iSub = 1:numSubs
%     for iCell = 1:length(cellcoll)
%         goodCells=find(~squeeze(any(any(any(any(any(badData(:,:,:,iSub,:,:,:),1),2),5),6),7)));
%         inputCells=intersect(cellcoll{iCell},goodCells);
%         if ~isempty(inputCells)
%             sevenDdata(:,:,iCell,iSub,:,:,:)=mean(theFacScr(:,:,inputCells,iSub,:,:,:),3);
%         end;
%     end;
% end;

for iCell = 1:numOutputCells
    for iChan=1:dataChanNum
        for iPoint=1:dataTimeNum
            for iSub=1:numSubs
                for iFac=1:numFacs
                    for iFreq=1:dataFreqNum
                        for iRel=1:dataRelNum
                            if badData(iChan,iPoint,iCell,iSub,1,iFreq,iRel)
                                sevenDdata(iChan,iPoint,iCell,iSub,iFac,iFreq,iRel)=NaN;
                            else
                                sevenDdata(iChan,iPoint,iCell,iSub,iFac,iFreq,iRel)=theFacScr(iChan,iPoint,iCell,iSub,iFac,iFreq,iRel);
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

%Construct parametric PCA output if so instructed
if ~isempty(params)
    
    %Compute correlations with factor scores
    scoresMatrix=reshape(permute(theFacScr,[4 3 1 2 6 7 5]),numSubs*numCells,[]);
    
    eval(['save ' fname 'S' num2str(sub) 'scores.mat scoresMatrix']);
    
    cellCounter=numCells;
    for i = 1:size(params,2)
        index = find(params(:,i)~=-999);	%observations with non-missing params and thus retained for correlation
        correl = corrcoef([scoresMatrix(index,:) params(index,i)]);	%column vector of correlations
        correl = correl(1:end-1,end);	%extract just the correlations between the parameter and the chans (e.g., chans,factors)
        corrMatrix = reshape(correl,ObsNum,numFacs);	%matrix of correlations (e.g., chans,factors)
        
        %merge param missing data info and bad data info
        paramBadData=badData;
        for iCellParam=1:numCells
            for iSubParam=1:numSubs
                if params((iCellParam-1)*numSubs+iSubParam,i)==-999
                    paramBadData(:,:,iCellParam,iSubParam,:,:,:)=1;
                end;
            end;
        end;
        
        cellCounter=cellCounter+1;
        for iChan=1:dataChanNum
            for iPoint=1:dataTimeNum
                for iSub=1:numSubs
                    for iFac=1:numFacs
                        for iFreq=1:dataFreqNum
                            for iRel=1:dataRelNum
                                corrNum=iChan+(iPoint-1)*dataChanNum+(iFreq-1)*dataChanNum*dataTimeNum+(iRel-1)*dataChanNum*dataTimeNum*dataFreqNum;
                                sevenDdata(iChan,iPoint,cellCounter,iSub,iFac,iFreq,iRel)=abs(mean(theFacScr(iChan,iPoint,find(~squeeze(paramBadData(iChan,iPoint,:,iSub,1,iFreq,iRel))),iSub,iFac,iFreq,iRel),3))*corrMatrix(corrNum,iFac);
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    numOutputCells=numOutputCells+size(params,2);
    cnames = [cnames; paramNames];
    [cellTypes{end+1:size(params,2),1}] = deal('CMB');
end;

data.data=sevenDdata;
data.noise=[];
data.std=[];
data.stdCM=[];
data.cov=[];
data.montage=FactorResults.montage;
data.fileFormat=FactorResults.fileFormat;
data.chanNames=FactorResults.chanNames;
data.timeNames=FactorResults.timeNames;
data.subNames=FactorResults.subNames;

data.cellNames=cnames;
data.trialNames=FactorResults.trialNames;
data.facNames=facNames;
data.freqNames=FactorResults.freqNames;
data.relNames=FactorResults.relNames;
data.chanTypes=FactorResults.chanTypes;
data.timeUnits=FactorResults.timeUnits;
data.subTypes=FactorResults.subTypes;

data.cellTypes=cellTypes;
data.facTypes=facTypes;

try
    EPver=ver('EP_Toolkit');
catch
    EPver='unavailable'; %workaround for bug in earlier version of Matlab
end;
data.EPver=EPver;
data.ver=ver;
data.date=date;
data.Fs=FactorResults.Fs;
data.baseline=FactorResults.baseline;

[pathstr, name, ext] = fileparts(FactorResults.fileName);

data.dataName=name;
data.ename=FactorResults.ename;
data.dataType=FactorResults.dataType;

if isfield(FactorResults,'trialSpecs')
    data.trialSpecs=FactorResults.trialSpecs;
else
    data.trialSpecs=cell(length(cnames),0);
end;
if isfield(FactorResults,'trialSpecNames')
    data.trialSpecNames=FactorResults.trialSpecNames;
else
    data.trialSpecNames=[];
end;

data.subjectSpecs=FactorResults.subjectSpecs;
data.subjectSpecNames=FactorResults.subjectSpecNames;

if isfield(FactorResults,'events')
    data.events=FactorResults.events;
else
    data.events=cell(numSubs,numOutputCells);
end;

if isfield(FactorResults,'avgNum')
    data.avgNum=FactorResults.avgNum;
else
    data.avgNum=zeros(numSubs,numOutputCells);
end;
if isfield(FactorResults,'subNum')
    data.subNum=FactorResults.subNum;
else
    data.subNum=ones(numSubs,numOutputCells)*numSubs;
end;
if isfield(FactorResults,'covNum')
    data.covNum=FactorResults.avgNum;
else
    data.covNum=zeros(numSubs,numOutputCells);
end;

data.fileName=FactorResults.fileName;
data.history={'PCAoutput',params, paramNames};
data.ced=FactorResults.ced;
data.eloc=FactorResults.eloc;
data.implicit=FactorResults.implicit;

data.facVecT=facVecT;
data.facVecS=facVecS;
data.facVecF=facVecF;
data.facData=[];
data.facVar=facVar;
data.facVarQ=facVarQ;
data.reference=FactorResults.reference;

if strcmp('continuous',FactorResults.dataType)
    numEpochs=floor(size(data.data,2)/ceil(data.Fs)); %excess time points are tacked onto final epoch
    data.analysis.blinkTrial=zeros(numSubs,numEpochs);
    data.analysis.saccadeTrial=zeros(numSubs,numEpochs);
    data.analysis.saccadeOnset=zeros(numSubs,numEpochs);
    data.analysis.moveTrial=zeros(numSubs,numEpochs);
    data.analysis.badTrials=zeros(numSubs,numEpochs);
    data.analysis.badChans=zeros(numSubs,numEpochs,numChans);
else
    data.analysis.blinkTrial=zeros(numSubs,numOutputCells);
    data.analysis.saccadeTrial=zeros(numSubs,numOutputCells);
    data.analysis.saccadeOnset=zeros(numSubs,numOutputCells);
    data.analysis.moveTrial=zeros(numSubs,numOutputCells);
    data.analysis.badTrials=zeros(numSubs,numOutputCells);
    data.analysis.badChans=zeros(numSubs,numOutputCells,numChans);
    for iCell=1:numCells
        for iSub=1:numSubs
            data.analysis.blinkTrial(iSub,iCell)=FactorResults.analysis.blinkTrial(iSub,iCell);
            data.analysis.saccadeTrial(iSub,iCell)=FactorResults.analysis.saccadeTrial(iSub,iCell);
            data.analysis.saccadeOnset(iSub,iCell)=FactorResults.analysis.saccadeOnset(iSub,iCell);
            data.analysis.moveTrial(iSub,iCell)=FactorResults.analysis.moveTrial(iSub,iCell);
            data.analysis.badTrials(iSub,iCell)=FactorResults.analysis.badTrials(iSub,iCell);
            data.analysis.badChans(iSub,iCell,:)=FactorResults.analysis.badChans(iSub,iCell,:);
            if all(reshape(badData(:,:,iCell,iSub,:,:,:),1,[])') && strcmp('average',FactorResults.dataType)
                data.avgNum(iSub,iCell)=-1;
                data.subNum(iSub,iCell)=-1;
            else
                for iChan=1:dataChanNum
                    if all(reshape(badData(iChan,:,iCell,iSub,:,:,:),1,[])')
                        if strcmp('average',FactorResults.dataType)
                            data.analysis.badChans(iSub,iCell,iChan)=NaN;
                        else
                            data.analysis.badChans(iSub,iCell,iChan)=-1;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

data.pca=FactorResults;
data.recTime=FactorResults.recTime;
data.stims=FactorResults.stims;
data.calibration=FactorResults.calibration;
data.impedances=FactorResults.impedances;

if numCells > 1
    %ensure that summary cell name is not already taken
    PCAname='all';
    sameName=1;
    suffix=0;
    PCAnameSuffix=PCAname;
    while sameName
        sameName=0;
        if any(strcmp(PCAnameSuffix,cnames))
            sameName=1;
        end;
        if sameName
            suffix=suffix+1;
            PCAnameSuffix=[PCAname '-' num2str(suffix)];
        end;
    end;
    [data]=ep_combineData(data,'cells',[1:numCells],ones(numCells,1),PCAnameSuffix,0); %add combined cells
end;

if numFacs > 1
    %ensure that summary factor name is not already taken
    PCAname='all';
    sameName=1;
    suffix=0;
    PCAnameSuffix=PCAname;
    while sameName
        sameName=0;
        if any(strcmp(PCAnameSuffix,facNames))
            sameName=1;
        end;
        if sameName
            suffix=suffix+1;
            PCAnameSuffix=[PCAname '-' num2str(suffix)];
        end;
    end;
    [data]=ep_combineData(data,'factors',[1:numFacs],ones(numFacs,1),PCAnameSuffix,0); %add combined factors
end;

if numSubs > 1
    if ~any(strcmp(FactorResults.dataType,{'continuous','single_trial'}))
        %ensure that summary subject name is not already taken
        PCAname='grand average';
        sameName=1;
        suffix=0;
        PCAnameSuffix=PCAname;
        while sameName
            sameName=0;
            if any(strcmp(PCAnameSuffix,data.subNames))
                sameName=1;
            end;
            if sameName
                suffix=suffix+1;
                PCAnameSuffix=[PCAname '-' num2str(suffix)];
            end;
        end;
        [data]=ep_combineData(data,'subjects',[1:numSubs],ones(numSubs,1),PCAnameSuffix,0); %add combined subjects
    end;
end;

fprintf('%s\r',' ' );
fprintf('%s\r','Done with reconstructing PCA waveforms...' );

[err]=ep_checkEPfile(data);

if err
    data=[];
    msg{1}='PCA file defective.  Contact programmer for support.';
    [msg]=ep_errorMsg(msg);
    return
end;
