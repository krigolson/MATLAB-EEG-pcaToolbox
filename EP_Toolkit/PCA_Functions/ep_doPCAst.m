function [FactorResultsST] = ep_doPCAst(FactorResults, ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, LOADING, PCAmode, crossVerify);
% ep_doPCAst - [FactorResultsST] = ep_doPCAst(FactorResults, ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, LOADING, PCAmode, crossVerify);
%Inputs
%  FactorResults: Structured variable with results and information about analysis.
%    See do PCA documentation.
%  ROTATION	: 
%       PCAs and has various options for the type of rotation, including
%       VMAX - Varimax
%       PMAX - Promax
%       IMAX - ICA Infomax (if EEGLab is also installed.  ICA is run with
%           the PCA subspace option on by default.  This program can be edited
%           to change its options.  Output is rescaled to be consistent with PCA conventions)
%       QMAX - Quartimax
%       QMIN - Quartimin
%       OMIN - Oblimin
%       CRFE - Crawford-Ferguson family
%       MINE - minimum entropy
%       IPSC - Bentler's invariant pattern simplicity criterion
%       TIIC - Comrey's tandem II criterion
%       GMIN - geomin
%       MMER - McCammon minimum entropy ratio
%  ROTOPT   : Rotation parameter for those having an optional parameter (Promax and Oblimin)
%  MAT_TYPE	: Matrix type (SCP: sums-of-squares-cross-product matrix, COV: variance-covariance matrix, COR: correlation matrix)
%  NUM_FAC	: Number of factors retained
%  LOADING	: Loading normalization ('K' Kaiser normalization, 'N' no normalization, 'C' covariance loadings, 'W' Cureton-Mulaik).
%  PCAmode  : Mode of new PCA ('temp': time points as the variables; 'spat': channels as the variables; 'freq': frequencies as variables)
%  crossVerify : factor pattern from previous PCA for cross-validation
%
%Outputs
% FactorResultsST: Structured variable with results and information about second PCA step.
%  .FacPatST	: Factor pattern matrix - produces standardized variables from scores, scaled by communality  (variables, factors)
%  .FacPat	: 1st step factor pattern matrix - produces standardized variables from scores, scaled by communality  (variables, factors)
%  .FacStrST	: Factor structure matrix - correlation between factors and variables (variables, factors)
%  .FacCofST : Factor scoring coefficients (rows=variables, cols=factors)
%  .FacScrST	: Factor scores, variance standardized, non-mean corrected (observations: cells*subjects, factors)
%  .FacCorST	: Factor correlations
%  .screeST    : The full set of unrotated eigenvalues, prior to truncation, for each input factor (output factors, input factors)
%  .facVarST   : %age of variance accounted for by each rotated factor, ignoring other factors.
%			  Calculated based on variance of relationship matrix (if COV or SCP, each variable weighted)
%  .facVarQST  : %age of variance uniquely accounted for by each rotated factor.
%			  Calculated based on variance of relationship matrix (if COV or SCP, each variable weighted)
%  .varSDST    : Standard deviations of the variables prior to the second PCA.
%  .varSD      : Standard deviations of the variables prior to the first PCA.
%  .PCAmode     : PCAmode from first step
%  .PCAmode2	: PCAmode from second step
%  .ROTATION    : Rotation from first step
%  .ROTATION2   : Rotation from second step
%  .MAT_TYPE    : Relationship matrix from first step
%  .MAT_TYPE2   : Relationship matrix from second step
%  .LOADING     : Loading option from first step
%  .LOADING2    : Loading option from second step
%  .RotOpt      : Rotation option from first step
%  .RotOpt2     : Rotation option from second step
%  .numFacs     : Number of factors from first step
%  .numFacs2    : Number of factors  from second step
%  .timepoints       : Number of time points in original data
%  .numchan          : Number of channels in original data
%  .numFreqs         : Number of frequency bins in the original data
%  .numRels          : Number of relations in the original data
%  .numSubs          : Number of subjects.
%  .numCells         : Number of cells
%  .montage          : The electrode montage.
%  .chanNames        : The channel names.
%  .timeNames        : The timepoint names.
%  .freqNames        : The frequency names.
%  .subNames         : The subject names.
%  .cellNames        : The cell names.
%  .relNames        : The relations names.
%  .chanTypes : The type of the channel: EEG, MGA (axial MEG), MGP (planar MEG), ECG, ANS (autonomic), REG (regional average)
%  .subTypes         : The type of the subject: RAW (single trial), AVG (subject average), GAV (grand average)
%  .cellTypes        : The type of the cell: SGL (one cell), CMB (combination of cells)
%  .facNames         : The factor names (only for factor files)
%  .facTypes         : The type of the factor: SGL (one factor), CMB (combination of factors)
%  .facNames2        : The factor names for the second step (only for factor files)
%  .facTypes2        : The type of the factor for the second step: SGL (one factor), CMB (combination of factors)
%  .cellNames        : The cell names.
%  .subjectSpecs     : Cell array (subject,spec) of specific information for subjects
%  .subjectSpecNames : Cell array of the name of each subject spec type.
%  .Fs               : The sampling frequency in Hz.
%  .baseline         : The number of samples in the baseline.
%  .ename            : The name of the experiment.
%  .dataType         : The type of the data: 'analysis'
%  .fileName         : Name of original file.
%  .history          : Command used to create current file.
%  .fileFormat       : The file format.
%  .ced       : The name of the .ced file for electrode coordinates.
%  .eloc      : The electrode location information, one for each channel (see readlocs header)
%  .implicit  : The electrode information for implicit references and fiducial locations (see readlocs header)
%  .implicit  : The electrode information for implicit references and fiducial locations (see readlocs header)
%    .reference
%        .original    : Original recording reference(s): 1-2 numbers in array
%        .current     : Current reference channel(s): 1-2 numbers in array
%        .type        : Current reference type: REG (regular), AVG (average reference, even if no longer adds to zero), CSD (current source density)
%
%  Based on doPCA.m.  Takes output of doPCA and performs a PCA on the factor scores in order to produce a spatiotemporal
%  or temporospatial PCA result.  A global number of factors is specified for the second stage PCA to keep things simple since
%  the statistics literature indicates the overretention of factors does not notably degrade solutions.
%  The resulting output has the same structure as the input data but with the number of factors multiplied by the
%  second stage factors.  The second stage factors vary most quickly (e.g.,
%  for temporospatial PCA, T1S1 T1S2 T1S3...T2S1...etc.).

%History
%  by Joseph Dien (1/10/01)
%  jdien07@mac.com
%
%  Modified (2/28/01) JD
%  Altered input names to facilitate use.
%
%  Modified (3/2/01)  JD
%  Dropped scoring coefficient output to avoid singularity problems when calculating it.
%
%  modified 7/27/02 JD
%  Implemented non mean-corrected factor scores.  Factor score output no longer has 1st stage factor patterns reintroduced.
%
%  modified 8/4/02 JD
%  Output generalized so that type of rotation can be specified.
%
%  modified 10/26/03 JD
%  Added option for unscaled loadings during rotation (not for Infomax).  Set Infomax
%  rotation to automatically use covariance matrix (i.e., mean-corrected
%  but not variance-corrected variables).
%
%  modified 3/23/04 JD
%  Folded LOADINGS option for covariance loadings into Kaiser parameter
%  since they're alternatives to each other.
%
%  bugfix 10/24/06 JD
%  If facVarQ is empty, as for current Infomax code, just don't calculate
%  it for the ST analysis.
%
%  modified (2/6/08) JD
%  Copyright notice appears only once per session.  Eliminated FacRef
%  output.  Output is now packaged in a structured variable.
%  Added RotOpt rotation parameter.  Changed KAISER to LOADING.
%
%  modified (11/7/08) JD
%  Input data comes in the form of a structured variable.  No longer
%  necessary to specify number of cells and subjects.  Includes information on montage and names for timepoints,
%  channels, subjects, and cells.
%
%  modified (1/31/09) JD
%  Added initEP.  Added version and date fields to output structure.
%
%  modified (3/12/09) JD
%  Added support for EP file format.
%
% modified 5/13/09 JD
% Added additional fields.
%
% modified 6/15/09 JD
% Added implicit, chanTypes, cellTypes, and subTypes fields.  FacVarST and FacVarQST
% now a vector rather than a 2D matrix.  Now accepts structured data variables without EP format info.
% Renamed fields to eliminate "first" prefix and replaced "second" prefix with a "2" suffix to make
% labels consistent with the doPCA fields.  Added matrices other than FacPat from the first step.
% Moved facName and facType to here from PCAoutput.
%
%  bugfix 8/3/09 JD
%  Fixed crash due to trialNames field when doing a PCA.
%
%  modified 2/21/10 JD
%  Analysis fields no longer optional.
%  Now zero-padding the factor names.
%
%  bugfix 5/22/09 JD
%  Fixed crash when a call to doPCA fails due to an error.
% 
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
% 
% modified 2/12/12 JD
% Added reference field to keep track of original and current reference scheme.
%
% modified 2/19/12 JD
% Added support for 6th dimension of frequency for time-frequency analyses
%
% modified 1/11/13 JD
% Added option to do internal calculations of frequency data in either amplitude or power form.
%
% bugfix 6/18/13 JD
% Fixed bug introduced in 2.30 wherein varSD for two-step PCAs applied incorrectly.
%
% modified 10/9/13 JD
% Added recTime field.
%
% bugfix 12/28/13 JD
% Fixed crash when performing two-step PCA.
%
% bugfix 2/27/14 JD
% facName and facType variables output as column vectors
%
% modified 3/23/14 JD
% Average numbers, trial specs, and events carried over to the PCA file.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% modified 4/28/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure, including complex numbers.
% Added ability to perform PCA on datasets with bad data.
%
% modified 11/25/14 JD
% Added cross-validation option.
%
% bugfix 12/17/15 JD
% Fixed crash when performing two-step PCA other than cross-validation.
% Fixed crash when performing second step of two-step PCA on older PCA datasets without a numRels field.
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
% modified 6/14/17 JD
% Added .timeUnits field.
%
% modified 11/22/17 JD
% Added support for impedances field.
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

FactorResultsST=[];

if nargin < 8
    crossVerify=[];
end

if ~isfield(FactorResults,'PCAmode')
    error('No PCAmode field.');
end
firstPCAmode=FactorResults.PCAmode;
if ~isfield(FactorResults,'FacScr')
    error('No FacScr field.');
end
FacScr=FactorResults.FacScr;
if ~isfield(FactorResults,'facVar')
    error('No facVar field.');
end
facVar=FactorResults.facVar;
if ~isfield(FactorResults,'facVarQ')
    error('No facVarQ field.');
end
facVarQ=FactorResults.facVarQ;
if ~isfield(FactorResults,'FacPat')
    error('No FacPat field.');
end
FacPat=FactorResults.FacPat;
if ~isfield(FactorResults,'varSD')
    error('No varSD field.');
end
varSD=FactorResults.varSD;
if ~isfield(FactorResults,'numchan')
    error('No numchan field.');
end
if ~isfield(FactorResults,'timepoints')
    error('No timepoints field.');
end
if ~isfield(FactorResults,'numFreqs')
    error('No numFreqs field.');
end
if ~isfield(FactorResults,'numSubs')
    error('No numSubs field.');
end
if ~isfield(FactorResults,'numCells')
    error('No numCells field.');
end
if ~isfield(FactorResults,'ROTATION')
    error('No ROTATION field.');
end
if ~isfield(FactorResults,'MAT_TYPE')
    error('No MAT_TYPE field.');
end
if ~isfield(FactorResults,'LOADING')
    error('No LOADING field.');
end
if ~isfield(FactorResults,'RotOpt')
    error('No RotOpt field.');
end
if ~isfield(FactorResults,'numFacs')
    error('No numFacs field.');
end
if ~isfield(FactorResults,'numRels')
    FactorResults.numRels=1;
end

if isfield(FactorResults,'dataType') %EP file format
    if ~isfield(FactorResults,'Fs')
        error('No Fs field.');
    end
    if ~isfield(FactorResults,'baseline')
        error('No baseline field.');
    end
    if ~isfield(FactorResults,'ename')
        error('No ename field.');
    end
end;

if strcmp(PCAmode,firstPCAmode)
    msg{1}=['A ' PCAmode ' PCA has already been applied to this data.'];
    [msg]=ep_errorMsg(msg);
    return
elseif isfield(FactorResults,'PCAmode2')
    if strcmp(PCAmode,FactorResults.PCAmode2)
        msg{1}=['A ' PCAmode ' PCA has already been applied to this data.'];
        [msg]=ep_errorMsg(msg);
        return
    end;
end;

dataChanNum=FactorResults.numchan;
dataTimeNum=FactorResults.timepoints;
dataFreqNum=FactorResults.numFreqs;
dataRelNum=FactorResults.numRels;

switch PCAmode
    case 'spat'
        NUM_VAR=FactorResults.numchan;
        varDim=1;
    case 'temp'
        NUM_VAR=FactorResults.timepoints;
        varDim=2;
    case 'freq'
        NUM_VAR=FactorResults.numFreqs;
        varDim=6;
    otherwise
        msg{1}='The PCAmode field must be set to either temp or spat or freq.';
        [msg]=ep_errorMsg(msg);
        return
end;

if isfield(FactorResults,'badObs')
    badObs=FactorResults.badObs;
else
    msg{1}='Error: This PCA data was generated prior changes made in EP 2.44 so it cannot be used for two-step PCA.  The first step will need to be rerun with the current EP Toolkit version.';
    [msg]=ep_errorMsg(msg);
    return
end;

FacPatST=[];
FacStrST=[];
FacCofST=[];
FacScrST=[];
FacCorST=[];
screeST=[];
facVarST=[];
facVarQST=[];
facVarTotST=[];
varSDST=[];
badObsST=[];

if isfield(FactorResults,'PCAmode')
    switch FactorResults.PCAmode
        case 'temp'
            dataTimeNum=1;
        case 'spat'
            dataChanNum=1;
        case 'freq'
            dataFreqNum=1;
    end
    if isfield(FactorResults,'PCAmode2')
        switch FactorResults.PCAmode2
            case 'temp'
                dataTimeNum=1;
            case 'spat'
                dataChanNum=1;
            case 'freq'
                dataFreqNum=1;
        end
    end;
end;

for factor = 1:size(FacScr,2)	%input factors
    %reshape factor vector so that the new mode is now the variables
    data=reshape(permute(reshape(FacScr(:,factor),dataChanNum,dataTimeNum,FactorResults.numCells,FactorResults.numSubs,1,dataFreqNum,dataRelNum),[varDim setdiff([1:7],varDim)]),NUM_VAR,[])';
    badData=reshape(permute(reshape(full(badObs),dataChanNum,dataTimeNum,FactorResults.numCells,FactorResults.numSubs,1,dataFreqNum,dataRelNum),[varDim setdiff([1:7],varDim)]),NUM_VAR,[])';
    theData=[];
    
    if isfield(FactorResults,'dataType') %EP file format
        theData.data=data;
        theData.montage=FactorResults.montage;
        theData.chanNames=FactorResults.chanNames;
        theData.timeNames=FactorResults.timeNames;
        theData.subNames=FactorResults.subNames;
        theData.cellNames=FactorResults.cellNames;
        theData.freqNames=FactorResults.freqNames;
        theData.relNames=FactorResults.relNames;
        theData.trialNames=FactorResults.trialNames;
        theData.dataType='PCAst';
        theData.fileName=FactorResults.fileName;
        theData.fileFormat=FactorResults.fileFormat;
        theData.eloc=FactorResults.eloc;
        theData.ced=FactorResults.ced;
        theData.implicit=FactorResults.implicit;
        theData.chanTypes=FactorResults.chanTypes;
        theData.timeUnits=FactorResults.timeUnits;
        theData.subTypes=FactorResults.subTypes;
        theData.cellTypes=FactorResults.cellTypes;
        theData.subjectSpecs=FactorResults.subjectSpecs;
        theData.subjectSpecNames=FactorResults.subjectSpecNames;
        theData.reference=FactorResults.reference;
        theData.recTime=FactorResults.recTime;
        theData.analysis=FactorResults.analysis;
        if isfield(FactorResults,'events')
            theData.events=FactorResults.events;
        end;
        if isfield(FactorResults,'avgNum')
            theData.avgNum=FactorResults.avgNum;
        end;
        if isfield(FactorResults,'covNum')
            theData.covNum=FactorResults.covNum;
        end;
        if isfield(FactorResults,'subNum')
            theData.subNum=FactorResults.subNum;
        end;
        if isfield(FactorResults,'trialSpecs')
            theData.trialSpecs=FactorResults.trialSpecs;
        end;
        if isfield(FactorResults,'trialSpecNames')
            theData.trialSpecNames=FactorResults.trialSpecNames;
        end;
    else
        theData=data;
    end;
    
    if ~isempty(crossVerify)
        theCrossVerify=crossVerify(:,(factor-1)*NUM_FAC+1:factor*NUM_FAC);
    else
        theCrossVerify=[];
    end;
    [stFactorResults] = ep_doPCA('asis', ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, theData, LOADING, 'N', badData, theCrossVerify);
    if isempty(stFactorResults)
        return
    end;
    FacPatST=[FacPatST stFactorResults.FacPat];
    FacStrST=[FacStrST stFactorResults.FacStr];
    FacCofST=[FacCofST stFactorResults.FacCof];    
    FacScrST=[FacScrST stFactorResults.FacScr];
    badObsST=[badObsST stFactorResults.badObs];
    FacCorST=[FacCorST stFactorResults.FacCor];
    screeST=[screeST stFactorResults.scree];
    facVarST = [facVarST (facVar(factor).* stFactorResults.facVar)];
    facVarQST = [facVarQST (facVarQ(factor).* stFactorResults.facVarQ)];
    facVarTotST=[facVarTotST stFactorResults.facVarTot];
    varSDST=[varSDST; stFactorResults.varSD'];
end

switch firstPCAmode
    case 'spat'
        label1='SF';
    case 'temp'
        label1='TF';
    case 'freq'
        label1='FF';
end

if isfield(FactorResults,'PCAmode2')
    switch PCAmode
        case 'spat'
            label3='SF';
        case 'temp'
            label3='TF';
        case 'freq'
            label3='FF';
    end
    switch FactorResults.PCAmode2
        case 'spat'
            label2='SF';
        case 'temp'
            label2='TF';
        case 'freq'
            label2='FF';
    end
else
    switch PCAmode
        case 'spat'
            label2='SF';
        case 'temp'
            label2='TF';
        case 'freq'
            label2='FF';
    end
    label3='';
end;

if isfield(FactorResults,'PCAmode2')
    factor1Digits=length(num2str(FactorResults.numFacs));
    factor2Digits=length(num2str(FactorResults.numFacs2));
    factor3Digits=length(num2str(NUM_FAC));
    
    for factor1=1:FactorResults.numFacs2
        for factor2=1:FactorResults.numFacs
            for factor3=1:NUM_FAC
                facNames{(factor3-1)*factor2*NUM_FAC+(factor2-1)*NUM_FAC+factor3,1}=sprintf('%s%s%s%s%s%s',label1,sprintf(['%0' num2str(factor1Digits) 'd'],factor1),label2,sprintf(['%0' num2str(factor2Digits) 'd'],factor2),label3,sprintf(['%0' num2str(factor3Digits) 'd'],factor3));
                facTypes{(factor3-1)*factor2*NUM_FAC+(factor2-1)*NUM_FAC+factor3,1}='SGL';
            end;
        end;
    end;
    
    FactorResultsST.FacPat3=FacPatST;
    FactorResultsST.FacPatST=FactorResults.FacPatST;
    FactorResultsST.FacPat=FactorResults.FacPat;
    
    FactorResultsST.FacStr3=FacStrST;
    FactorResultsST.FacStrST=FactorResults.FacStrST;
    FactorResultsST.FacStr=FactorResults.FacStr;
    
    FactorResultsST.FacScr3=FacScrST;
    FactorResultsST.FacScrST=FactorResults.FacScrST;
    FactorResultsST.FacScr=FactorResults.FacScr;
    
    FactorResultsST.badObs3=badObsST;
    FactorResultsST.badObsST=FactorResults.badObsST;
    FactorResultsST.badObs=FactorResults.badObs;

    FactorResultsST.FacCof3=FacCofST;
    FactorResultsST.FacCofST=FactorResults.FacCofST;
    FactorResultsST.FacCof=FactorResults.FacCof;
    
    FactorResultsST.FacCor3=FacCorST;
    FactorResultsST.FacCorST=FactorResults.FacCorST;
    FactorResultsST.FacCor=FactorResults.FacCor;
    
    FactorResultsST.screeST=FactorResults.screeST;
    FactorResultsST.scree3=screeST;
    
    FactorResultsST.facVar3=facVarST;
    FactorResultsST.facVarQ3=facVarQST;
    FactorResultsST.facVarTot3=facVarTotST;
    FactorResultsST.facVarST=FactorResults.facVarST;
    FactorResultsST.facVarQST=FactorResults.facVarQST;
    FactorResultsST.facVarTotST=FactorResults.facVarTotST;
    FactorResultsST.facVar=FactorResults.facVar;
    FactorResultsST.facVarQ=FactorResults.facVarQ;
    FactorResultsST.facVarTot=FactorResults.facVarTot;
    
    FactorResultsST.varSD3=varSDST;
    FactorResultsST.varSDST=FactorResults.varSDST;
    FactorResultsST.varSD=FactorResults.varSD;
    
    FactorResultsST.PCAmode=FactorResults.PCAmode;
    FactorResultsST.PCAmode2=FactorResults.PCAmode2;
    FactorResultsST.PCAmode3=PCAmode;
    
    FactorResultsST.ROTATION=FactorResults.ROTATION;
    FactorResultsST.ROTATION2=FactorResults.ROTATION2;
    FactorResultsST.ROTATION3=ROTATION;
    
    FactorResultsST.MAT_TYPE=FactorResults.MAT_TYPE;
    FactorResultsST.MAT_TYPE2=FactorResults.MAT_TYPE2;
    FactorResultsST.MAT_TYPE3=MAT_TYPE;
    
    FactorResultsST.LOADING=FactorResults.LOADING;
    FactorResultsST.LOADING2=FactorResults.LOADING2;
    FactorResultsST.LOADING3=LOADING;
    
    FactorResultsST.RotOpt=FactorResults.RotOpt;
    FactorResultsST.RotOpt2=FactorResults.RotOpt2;
    FactorResultsST.RotOpt3=ROTOPT;
    
    FactorResultsST.numFacs=FactorResults.numFacs;
    FactorResultsST.numFacs2=FactorResults.numFacs2;
    FactorResultsST.numFacs3=NUM_FAC;

    FactorResultsST.facNames=FactorResults.facNames;
    FactorResultsST.facTypes=FactorResults.facTypes;
    FactorResultsST.facNames2=FactorResults.facNames2;
    FactorResultsST.facTypes2=FactorResults.facTypes2;
    FactorResultsST.facNames3=facNames;
    FactorResultsST.facTypes3=facTypes;
else
    factor1Digits=length(num2str(FactorResults.numFacs));
    factor2Digits=length(num2str(NUM_FAC));
    
    for factor1=1:FactorResults.numFacs
        for factor2=1:NUM_FAC
            facNames{(factor1-1)*NUM_FAC+factor2,1}=sprintf('%s%s%s%s',label1,sprintf(['%0' num2str(factor1Digits) 'd'],factor1),label2,sprintf(['%0' num2str(factor2Digits) 'd'],factor2));
            facTypes{(factor1-1)*NUM_FAC+factor2,1}='SGL';
        end;
    end;
    
    FactorResultsST.FacPatST=FacPatST;
    FactorResultsST.FacPat=FactorResults.FacPat;
    FactorResultsST.FacStrST=FacStrST;
    FactorResultsST.FacStr=FactorResults.FacStr;
    FactorResultsST.FacScrST=FacScrST;
    FactorResultsST.FacScr=FactorResults.FacScr;
    FactorResultsST.badObsST=badObsST;
    FactorResultsST.badObs=FactorResults.badObs;
    FactorResultsST.FacCofST=FacCofST;
    FactorResultsST.FacCof=FactorResults.FacCof;
    FactorResultsST.FacCorST=FacCorST;
    FactorResultsST.FacCor=FactorResults.FacCor;
    FactorResultsST.screeST=screeST;
    FactorResultsST.facVarST=facVarST;
    FactorResultsST.facVarQST=facVarQST;
    FactorResultsST.facVarTotST=facVarTotST;
    FactorResultsST.facVar=FactorResults.facVar;
    FactorResultsST.facVarQ=FactorResults.facVarQ;
    FactorResultsST.facVarTot=FactorResults.facVarTot;
    FactorResultsST.varSDST=varSDST;
    FactorResultsST.varSD=FactorResults.varSD;
    FactorResultsST.PCAmode=firstPCAmode;
    FactorResultsST.PCAmode2=PCAmode;
    FactorResultsST.ROTATION=FactorResults.ROTATION;
    FactorResultsST.ROTATION2=ROTATION;
    FactorResultsST.MAT_TYPE=FactorResults.MAT_TYPE;
    FactorResultsST.MAT_TYPE2=MAT_TYPE;
    FactorResultsST.LOADING=FactorResults.LOADING;
    FactorResultsST.LOADING2=LOADING;
    FactorResultsST.RotOpt=FactorResults.RotOpt;
    FactorResultsST.RotOpt2=ROTOPT;
    FactorResultsST.numFacs=FactorResults.numFacs;
    FactorResultsST.numFacs2=NUM_FAC;
    FactorResultsST.facNames=FactorResults.facNames;
    FactorResultsST.facTypes=FactorResults.facTypes;
    FactorResultsST.facNames2=facNames;
    FactorResultsST.facTypes2=facTypes;
end

try
    EPver=ver('EP_Toolkit');
catch
    EPver='unavailable'; %workaround for bug in earlier version of Matlab
end;
FactorResultsST.EPver=EPver;
FactorResultsST.ver=ver;
FactorResultsST.date=date;
FactorResultsST.timepoints=FactorResults.timepoints;
FactorResultsST.numchan=FactorResults.numchan;
FactorResultsST.numFreqs=FactorResults.numFreqs;
FactorResultsST.numSubs=FactorResults.numSubs;
FactorResultsST.numCells=FactorResults.numCells;
FactorResultsST.numRels=FactorResults.numRels;
FactorResultsST.history={'doPCAST',ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, LOADING};

if isfield(FactorResults,'dataType') %EP file format
    FactorResultsST.montage=FactorResults.montage;
    FactorResultsST.chanNames=FactorResults.chanNames;
    FactorResultsST.timeNames=FactorResults.timeNames;
    FactorResultsST.subNames=FactorResults.subNames;
    FactorResultsST.cellNames=FactorResults.cellNames;
    FactorResultsST.trialNames=FactorResults.trialNames;
    FactorResultsST.chanTypes=FactorResults.chanTypes;
    FactorResultsST.timeUnits=FactorResults.timeUnits;    
    FactorResultsST.subTypes=FactorResults.subTypes;
    FactorResultsST.cellTypes=FactorResults.cellTypes;
    FactorResultsST.facNames=FactorResults.facNames;
    FactorResultsST.freqNames=FactorResults.freqNames;
    FactorResultsST.relNames=FactorResults.relNames;
    FactorResultsST.facTypes=FactorResults.facTypes;
    FactorResultsST.subjectSpecs=FactorResults.subjectSpecs;
    FactorResultsST.subjectSpecNames=FactorResults.subjectSpecNames;
    FactorResultsST.eloc=FactorResults.eloc;
    FactorResultsST.ced=FactorResults.ced;
    FactorResultsST.implicit=FactorResults.implicit;
    FactorResultsST.Fs=FactorResults.Fs;
    FactorResultsST.baseline=FactorResults.baseline;
    FactorResultsST.ename=FactorResults.ename;
    FactorResultsST.dataType=FactorResults.dataType;
    FactorResultsST.fileName=FactorResults.fileName;
    FactorResultsST.fileFormat=FactorResults.fileFormat;
    FactorResultsST.analysis.blinkTrial=FactorResults.analysis.blinkTrial;
    FactorResultsST.analysis.saccadeTrial=FactorResults.analysis.saccadeTrial;
    FactorResultsST.analysis.saccadeOnset=FactorResults.analysis.saccadeOnset;
    FactorResultsST.analysis.moveTrial=FactorResults.analysis.moveTrial;
    FactorResultsST.analysis.badTrials=FactorResults.analysis.badTrials;
    FactorResultsST.analysis.badChans=FactorResults.analysis.badChans;
    FactorResultsST.reference=FactorResults.reference;
    FactorResultsST.recTime=FactorResults.recTime;
    
    if isfield(FactorResults,'events')
        FactorResultsST.events=FactorResults.events;
    end;
    if isfield(FactorResults,'stims')
        FactorResultsST.stims=FactorResults.stims;
    end;
    if isfield(FactorResults,'calibration')
        FactorResultsST.calibration=FactorResults.calibration;
    end;
    if isfield(FactorResults,'impedances')
        FactorResultsST.impedances=FactorResults.impedances;
    end;
    if isfield(FactorResults,'avgNum')
        FactorResultsST.avgNum=FactorResults.avgNum;
    end;
    if isfield(FactorResults,'covNum')
        FactorResultsST.covNum=FactorResults.covNum;
    end;
    if isfield(FactorResults,'subNum')
        FactorResultsST.subNum=FactorResults.subNum;
    end;
    if isfield(FactorResults,'trialSpecs')
        FactorResultsST.trialSpecs=FactorResults.trialSpecs;
    end;
    if isfield(FactorResults,'trialSpecNames')
        FactorResultsST.trialSpecNames=FactorResults.trialSpecNames;
    end;
end;


