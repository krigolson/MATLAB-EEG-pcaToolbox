function [FactorResults] = ep_doPCA(PCAmode, ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, theData, LOADING, GAVE, badData, crossFacCof)

% ep_doPCA - [FactorResults] = ep_doPCA(PCAmode, ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, theData, LOADING, GAVE, badData, crossFacCof) -
%          do PCA and calculate factor scores
%Inputs
%  PCAmode	 : Primary mode for the PCA ('temp': time points as the variables; 'spat': channels as the variables; 'freq': frequencies as variables; 'asis': do not rearrange input data).
%  ROTATION	:
%       PCAs and has various options for the type of rotation, including
%       UNRT - Unrotated
%       VMAX - Varimax
%       PMAX - Promax (kappa=3 seems to give the best results for ERPs)
%       IMAX - ICA Infomax (if EEGLab is also installed.  ICA is run with
%           the PCA subspace option on by default.  This program can be edited
%           to change its options.  Output is rescaled to be consistent with PCA conventions)
%       QMAX - Quartimax
%       QMIN - Quartimin
%       OMIN - Oblimin (gamma=0 tends to be recommended)
%       CRFE - Crawford-Ferguson family
%       MINE - minimum entropy
%       IPSC - Bentler's invariant pattern simplicity criterion
%       TIIC - Comrey's tandem II criterion
%       GMIN - geomin
%       MMER - McCammon minimum entropy ratio
%       VOMN - Variable Oblimin (tries to choose optimal gamma value)
%  ROTOPT   : Rotation parameter for those having an optional parameter (Promax and Oblimin)
%  MAT_TYPE	: Matrix type (SCP: sums-of-squares-cross-product matrix, COV: variance-covariance matrix, COR: correlation matrix)
%  NUM_FAC	: Number of factors retained
%  theData  : Structured array with the data and accompanying information.  See readData.
%  LOADING	: Loading normalization ('K' Kaiser normalization, 'N' no
%  normalization, 'C' covariance loadings, 'W' Cureton-Mulaik weighting).
%
%  OPTIONAL
%  GAVE     : Convert the data to grand average form, perform the analysis,
%               then reconvert the factor scores back to average form, 'Y' or 'N'.
%  badData  : sparse matrix with same dimensions as theData where 1 denotes bad data.
%  crossFacCof : factor scoring coefficients from previous PCA for cross-verification
%
%Outputs
%  FactorResults: Structured variable with results and information about analysis.
%  .FacPat	: Factor pattern matrix - produces standardized variables from scores, scaled by communality  (rows=variables, cols=factors)
%  .FacStr	: Factor structure matrix - correlation between factors and variables (rows=variables, cols=factors)
%  .FacScr	: Factor scores, variance standardized, non-mean corrected  (rows=variables, cols=factors)
%             The ordering of the data dimensions in the rows (which ones vary the most quickly) is in the order of the
%             seven data dimensions (earlier varying faster), with the one serving as the column dimension left out.
%  .FacCof   : Factor scoring coefficients (rows=variables, cols=factors)
%  .FacCor	: Factor correlations
%  .scree    : The full set of unrotated eigenvalues, prior to truncation.
%  .facVar   : %age of variance accounted for by each rotated factor, ignoring other factors.
%			  Calculated based on variance of relationship matrix (if COV or SCP, each variable weighted).
%  .facVarQ  : %age of variance uniquely accounted for by each rotated factor.
%			  Calculated based on variance of relationship matrix (if COV or SCP, each variable weighted)
%  .facVarTot   : Total variance accounted for.
%  .varSD    : Standard deviations of the variables.
%  .PCAmode     : PCAmode used
%  .ROTATION    : Rotation used
%  .MAT_TYPE    : Relationship matrix used
%  .KAISER      : Loading option used
%  .timepoints  : Number of time points in original data
%  .numchan     : Number of channels in original data
%  .numFreqs    : Number of frequency bins in the original data
%  .numSubs     : Number of subjects.
%  .numCells    : Number of cells
%  .numFacs     : Number of factors
%  .montage     : The electrode montage.
%  .chanNames   : The channel names.
%  .timeNames   : The timepoint names.
%  .freqNames   : The frequency names.
%  .subNames    : The subject names.
%  .cellNames   : The cell names.
%  .facNames  : The factor names (only for factor files)
%  .chanTypes : The type of the channel: EEG, MGA (axial MEG), MGP (planar MEG), ECG, ANS (autonomic), REG (regional average)
%  .subTypes  : The type of the subject: RAW (single trial), AVG (subject average), GAV (grand average)
%  .cellTypes : The type of the cell: SGL (one cell), CMB (combination of cells)
%  .facTypes  : The type of the factor: SGL (one factor), CMB (combination of factors)
%  .subjectSpecs     : Cell array (subject,spec) of specific information for subjects
%  .subjectSpecNames : Cell array of the name of each subject spec type.
%  .Fs          : The sampling frequency in Hz.
%  .baseline    : The number of samples in the baseline.
%  .ename       : The name of the experiment.
%  .dataType    : The type of the data: 'analysis'
%  .fileName  : Name of original file.
%  .history   : Command used to create current file.
%  .fileFormat: The file format.
%  .ced       : The name of the .ced file for electrode coordinates.
%  .eloc      : The electrode location information, one for each channel (see readlocs header)
%  .implicit  : The electrode information for implicit references and fiducial locations (see readlocs header)
%    .reference
%        .original    : Original recording reference(s): 1-2 numbers in array
%        .current     : Current reference channel(s): 1-2 numbers in array
%        .type        : Current reference type: REG (regular), AVG (average reference, even if no longer adds to zero), CSD (current source density)
%  .badObs   : sparse column vector of same number of rows as the factor score matrix where 0=good data and 1=bad data.
%    .noise     : 4D matrix [channels, time points, cells/trials, subjects] mirroring .data.
%                 This is the +/- reference average (Schimmel, 1967) which provides an estimate of the noise level
%                 in averaged data by flipping every other trial.  Not applicable to spectral data.  Not from combining data other than across subjects.
%                 For grand average data, every other average is flipped.
%
%  doPCA runs the PCA on the data.  It can do either spatial or temporal,
%  using a number of rotations and rotation options.
%  See accompanying readme file for information about appropriate references for the rotations.

%History
%  by Joseph Dien (4/99)
%  jdien07@mac.com
%
%  with assistance by Bill Dunlap and dedicated to his memory.
%  also with thanks to Kris Preacher.
%
%  modified 6/00 JD
%  Scree output added.
%
%  modified 9/30/00 JD
%  Kaiser normalization can be turned off and %age variance accounted for calculated.
%
%  modified 1/10/00 JD
%  %age of variance accounted for by each Promax rotated factor returned.  Fixed error in Promax algorithm.
%  changed output variable names to more readable names.  Added varSD output.
%
%  corrected 2/27/01 JD
%  fixed error in factor score calculations for cov and sscp matrices.
%
%  modified 3/2/01 JD
%  Rotating factor scores directly to avoid singularity problems from generalized inverse of relation matrix.  Dropped scoring coefficient output.
%
%  modified 7/27/02 JD
%  Implemented non mean-corrected factor scores.
%
%  modified 8/4/02 JD
%  Output generalized so that type of rotation can be specified.  Variance correcting factor scores after varimax rotation.
%
%  modified 10/29/02 JD
%  Added option for unrotated output.
%
%  bugfix 12/5/02 JD
%  Prerotation factor score standardization was incorrect (recently
%  introduced error).  Not entirely sure that SSCP is implemented
%  correctly, although it matches up with JMP output.  No time for now to
%  work out.
%
%  modified 4/25/03 JD
%  Added ICA, also known as infomax rotation, to rotation options.  FacRef more properly set to empty
%  matrix for rotations other than Promax.
%
%  modified 10/26/03 JD
%  Added option for unscaled loadings during rotation (not for Infomax).  Set Infomax
%  rotation to automatically use covariance matrix (i.e., mean-corrected
%  but not variance-corrected variables).
%
%  bugfix 11/18/03 JD
%  Fixed bug to ICA factor score computation.
%
%  bugfix 2/15/04 JD
%  Added fix for when using both covariance loadings and varimax rotation
%  together.  Fixed unscaled loadings option.  Now specified as alternative to Kaiser normalization
%  and not operational during Promax step.
%
%  modified 12/29/04 JD
%  ICA option set to reinitialize random number generator to the starting
%  state with each use.  Also, set ICA routine to non-verbose output.
%
%  modified 2/14/05 JD
%  Added grand average analysis option.
%
%  modified 10/24/06 JD
%  Set FacVarQ to FacVar for Varimax.
%
%  bugfix 11/4/06 JD
%  scree for ICA should be based on the unrotated PCA if PCA is being used
%  as the pre-processing step.  Was previously doing it based on the rotated
%  ICA solution.  Furthermore, the scree vector for ICA was coming out as a row
%  vector rather than a column vector so doPCAst was turning the screeST
%  matrix into one long vector instead of separating out the screes for each
%  separate factor.  Now fixed.
%
%  modified (12/19/07) JD
%  Eliminated indexdata.  Replaced with numCells and numSubs inputs.
%  Changed input data to 3D array.  Shifted to GNU license.
%
%  modified (2/4/08) JD
%  Changed direct rotation of factor scores back to direct computation of
%  scores since it turned out to be more accurate.
%
%  bugfix (2/4/08) JD
%  Fixed bug in setting up data for spatial analyses.
%
%  modified (2/12/08) JD
%  Copyright notice appears only once per session.
%  Incorporated GPF rotations.  Output is now packaged in a
%  structured variable.  Added RotOpt rotation parameter.
%  Generalized variance accounted for calculation to all rotations.
%  Added fix for block size calculation of older versions of runICA.
%  Sorts factors in order of size (variance accounted for).
%  Moved loading weightings out of Varimax so they can apply to all rotations.
%  Changed KAISER to LOADING.  Looks for binary version of runica.
%
%  modified (3/5/08) JD
%  Tests for bad loadings and communalities.  Adding total variance accounted for to
%  the output variable.
%
%  bugfix (8/16/08) JD
%  Ensures that FacPat and FacStr are column vectors for one-factor
%  solutions to fix problem with Matlab version incompatibilities.
%
%  bugfix (8/22/08) JD
%  Fixed calculation of unique variance accounted for by factors.
%
%  modified (11/5/08) JD
%  Input data comes in the form of a structured variable.  No longer
%  necessary to specify number of cells and subjects.  Includes information on montage and names for timepoints,
%  channels, subjects, and cells.
%
%  bugfix (1/16/09) JD
%  Crashing for spatial PCA (unless done via doPCAst).
%
%  modified (2/4/09) JD
%  Added initEP.  Added version and date fields to output structure.  Accepts theData without fields other than .data.
%
%  modified (2/24/09) JD
%  Added support for sampling rate, baseline, and experiment name.
%
%  modified (3/12/09) JD
%  Added support for EP file format.
%
%  modified (5/6/09) JD
%  Added cellType and subType fields.
%
%  modified (6/15/09) JD
%  Use of EP file format input data is optional.  If not used, output will also not have EP file format elements.
%  Moved facName and facType to here from PCAoutput.
%
%  modified 8/28/09 JD
%  Modified to use new option for binary runICA to control the output file names so that they can be cleaned up
%  and so that unpredictable crash from bug causing strange file names can be circumvented.
%
%  bugfix (10/18/09) JD
%  Catches situation where ICA fails to converge on a solution.
%
%  bugfix 12/3/09 JD
%  Additional check for ICA failure (denoted by imaginary numbers for weights).
%
%  modified 2/15/10 JD
%  Analysis fields no longer optional.
%  Now zero-padding the factor names to max number of digits instead of always to three.
%
%  bugfix 3/5/10 JD
%  Fixed crash when using grand average option (only useable from command line).
%
%  bugfix 4/16/10 JD
%  Fixed crash when retaining only one factor for an Infomax rotation.
%
%  bugfix 7/22/10 JD
%  Fixed crash when there are variables which are flat (zero standard deviation).
%
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
%
% modified 1/25/12 JD
% Added reference field to keep track of original and current reference scheme.
%
% modified 2/17/12 JD
% Added support for 6th dimension of frequency for time-frequency analyses
%
% bugfix 4/29/12 JD
% Fixed PCA not restoring flat variables back to data after PCA, as in reference channels.
%
% modified 1/11/13 JD
% Added option to do internal calculations of frequency data in either amplitude or power form.
%
% modified 1/15/13 JD
% Detects when eigenvalue decomposition yields complex numbers due to rounding errors causing relationship matrix to not
% be symmetric and makes the relationship matrix symmetric.
%
% modified 1/16/13 JD
% Since eigenvalue decomposition is not always in ascending order, sort them first.
%
% modified 10/9/13 JD
% Added recTime field.
%
% bugfix 1/29/14 JD
% Fixed crash when no factors retained.
%
% bugfix 2/7/14 JD
% Fixed aborting PCA when factor loadings slightly over 1 due to rounding errors.
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
% FacScr observations are now always arranged with permutations in order of the seven data dimensions.
%
% modified 11/25/14 JD
% Added cross-validation option.
%
% bugfix 12/18/15 JD
% Preserves order of factors when cross-validating.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% modified 9/19/16 JD
% Changed runica call so that always initialized with same seed so that ICA results are fully replicable.
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
% modified 3/27/18 JD
% Added support for noise field.
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

FactorResults=[];

if nargin ==0
    error('You need to provide the input information. See the header of this file for more information.');
end;

if nargin < 7
    error('Too few input parameters in doPCA command.');
end

if nargin < 8
    GAVE = 'N';
end

if nargin < 9
    badData=[];
end

if nargin < 10
    crossFacCof=[];
end

if ~any(strcmp(ROTATION,{'QMIN', 'OMIN', 'CRFE', 'MINE', 'IPSC', 'TIIC', 'GMIN', 'IMAX', 'MMER', 'UNRT', 'VMAX','QMAX','PMAX', 'VOMN'}))
    error('Rotation not recognized.');
end;

switch ROTATION
    case 'QMIN'
        theRotationName = 'Quartimin';
    case 'OMIN'
        theRotationName = 'Oblimin';
    case 'CRFE'
        theRotationName = 'Crawford-Ferguson';
    case 'MINE'
        theRotationName = 'minimum entropy';
    case 'IPSC'
        theRotationName = 'Bentler''s invariant pattern simplicity criterion';
    case 'TIIC'
        theRotationName = 'Comrey''s tandem II criterion';
    case 'GMIN'
        theRotationName = 'geomin';
    case 'IMAX'
        theRotationName = 'Infomax';
    case 'MMER'
        theRotationName = 'McCammon minimum entropy ratio';
    case 'UNRT'
        theRotationName = 'none';
    case 'VMAX'
        theRotationName = 'Varimax';
    case 'QMAX'
        theRotationName = 'Quartimax';
    case 'PMAX'
        theRotationName = 'Promax';
    case 'VOMN'
        theRotationName = 'Variable Oblimin';
end;

if (strcmp(ROTATION,'PMAX') && ROTOPT <=0)
    disp('Promax kappa needs to be larger than zero.  Setting it to default value of 3.');
    ROTOPT=3;
end;

if isstruct(theData)
    data=theData.data;
    if isempty(badData)
        %mark the bad data
        badData=zeros(size(data));
        if strcmp(theData.dataType,'average')
            for iSub=1:size(data,4)
                for iCell=1:size(data,3)
                    badChans=find(isnan(theData.analysis.badChans(iSub,iCell,:)));
                    badData(badChans,:,iCell,iSub,:,:,:)=1;
                    if ~isempty(theData.relNames)
                        badData(:,:,iCell,iSub,:,:,badChans)=1;
                    end;
                    if theData.avgNum(iSub,iCell)==-1
                        badData(:,:,iCell,iSub,:,:,:)=1;
                    end;
                    if theData.subNum(iSub,iCell)==-1
                        badData(:,:,iCell,iSub,:,:,:)=1;
                    end;
                end;
            end;
        else
            badData(:,:,find(theData.analysis.badTrials(1,:)),:,:,:,:)=1;
            for iCell=1:size(data,3)
                badChans=find(theData.analysis.badChans(1,iCell,:)==-1);
                badData(badChans,:,iCell,:,:,:,:)=1;
                if ~isempty(theData.relNames)
                    badData(:,:,iCell,:,:,:,badChans)=1;
                end;
            end;
        end;
    end;
else
    data=theData;
    badData=zeros(size(data));
end;

badData(isnan(data))=1; %mark nan data bad, as in the reference channel in coherence data.

if strcmp(GAVE,'Y')
    rawData=data;
    data=mean(data,4);
elseif ~strcmp(GAVE,'N')
    error('GAVE field must equal either "N" or "Y".');
end;

numChans=size(data,1);
numPoints=size(data,2);
numCells=size(data,3);
numSubs=size(data,4);
numFreqs=size(data,6);
numRels=size(data,7);

%reshape the data matrix
if strcmp(PCAmode,'temp')
    data=reshape(permute(data,[2 1 3 4 5 6 7]),numPoints,[])';
    badData=reshape(permute(badData,[2 1 3 4 5 6 7]),numPoints,[])';
    if isempty(NUM_FAC)
        NUM_FAC=numPoints;
    end;
    if GAVE == 'Y'
        rawData=reshape(permute(rawData,[2 1 3 4 5 6 7]),numPoints,[])';
    end;
elseif strcmp(PCAmode,'spat')
    data=reshape(data,numChans,[])';
    badData=reshape(badData,numChans,[])';
    if isempty(NUM_FAC)
        NUM_FAC=numChans;
    end;
    if GAVE == 'Y'
        rawData=reshape(rawData,numChans,[])';
    end;
elseif strcmp(PCAmode,'freq')
    data=reshape(permute(data,[6 1 2 3 4 5 7]),numFreqs,[])';
    badData=reshape(permute(badData,[6 1 2 3 4 5 7]),numFreqs,[])';
    if isempty(NUM_FAC)
        NUM_FAC=numFreqs;
    end;
    if GAVE == 'Y'
        rawData=reshape(permute(rawData,[6 1 2 3 4 5 7]),numFreqs,[])';
    end;
elseif strcmp(PCAmode,'asis')
    if isempty(NUM_FAC)
        NUM_FAC=size(data,1);
    end;
else
    error('PCAmode must be set to either temp or spat or asis');
end;

dataAll=data;
goodVars=intersect(find(std(data) ~= 0),find(~all(badData)));
if ~any(goodVars)
    msg{1}='No good variables for PCA.';
    [msg]=ep_errorMsg(msg);
    return
end;
if ~all(goodVars)
    disp('At least one variable temporarily dropped from PCA to avoid crashing algorithm.  Will be restored to factor matrices with zero loadings.');
end;

if size(data,1) == 1
    msg{1}='It is not possible to do a PCA with only one observation.';
    [msg]=ep_errorMsg(msg);
    return
end;

goodObs=find(~any(badData(:,goodVars)'))';
badObs=sparse(any(badData(:,goodVars)'))';

if length(goodObs)<2
    msg{1}='Too much bad data to conduct PCA.';
    [msg]=ep_errorMsg(msg);
    return
end;

data=data(goodObs,goodVars);

if size(data,2) < NUM_FAC
    msg{1}=['Number of variables (' num2str(size(data,2)) ') is smaller than the number of factors to be retained (' num2str(NUM_FAC) '.'];
    [msg]=ep_errorMsg(msg);
    return
end;

if size(data,1) < NUM_FAC
    msg{1}=['Number of observations (' num2str(size(data,1)) ') is smaller than the number of factors to be retained (' num2str(NUM_FAC) ').'];
    [msg]=ep_errorMsg(msg);
    return
end;

if ~isreal(data) %if the data has an imaginary component, as in spectral data
    data=[real(data);imag(data)];
    isComplex=1;
else
    isComplex=0;
end;

varSD = std(data);
D = diag(varSD);  %diagonal matrix of standard deviations of variables

if strcmp(ROTATION,'IMAX')
    if (strcmp(MAT_TYPE,'COR') || strcmp(MAT_TYPE,'SCP'))
        MAT_TYPE = 'COV';
        disp('Infomax rotations automatically carried out with covariance rotation.');
    end;
    if exist('icadefs','file') ~= 2
        msg{1}='ICA not available.  You need to download EEGlab and then add it to Matlab''s path list.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    icadefs;
    if exist(ICABINARY) == 2 %is there a binary file of runica installed?
        if any(isspace(pwd))
            msg{1}=['Binary ICA does not work if any part of the pathname (' pwd ') has a space in it.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
    end;
end;

R=corrcoef(data);

if strcmp(MAT_TYPE,'SCP')
    S = ((data)' * (data))/(size(data,1) - 1); %SSCP
    Sd = diag(sqrt(diag(S)));  %diagonal matrix of standard deviations of variables as used to generate relationship matrix
    disp('Warning:algorithms involving SSCP relationship matrix have not been fully validated.  Use with caution.');
elseif strcmp(MAT_TYPE,'COV')
    S = cov(data); %covariance matrix
    Sd = diag(sqrt(diag(S)));  %diagonal matrix of standard deviations of variables as used to generate relationship matrix
elseif strcmp(MAT_TYPE,'COR')
    S = R; %correlation matrix
    Sd = diag(ones(size(data,2),1));  %diagonal of standard deviations of variables as used to generate relationship matrix
else
    error('Error - matrix type not specified');
end

% if strcmp(ROTATION,'IMAX') && ~strcmp(LOADING,'N')
%     disp(['Loading weights (the ' LOADING ' option) cannot be applied to this implementation of Infomax so this setting will be ignored.']);
% end;

if strcmp(ROTATION,'IMAX') && ~strcmp(MAT_TYPE,'COV')
    disp(['Matrix type (the ' MAT_TYPE ' option) will instead be treated as ''COV''.']);
end;

% if strcmp(ROTATION,'IMAX') && (ROTOPT > 0)
%     disp(['Rotation option (the ' num2str(ROTOPT) ' setting) is not relevant to Infomax rotation and will be ignored.']);
% end;

if ~isempty(crossFacCof)

    FacCof=crossFacCof;
    FacScr=(data)*FacCof;
    FacStr=corrcoef([FacScr data]); %compute factor structure matrix
    FacStr=FacStr(NUM_FAC+1:size(FacStr,1),1:NUM_FAC);
    FacCor = corrcoef(FacScr);
    FacPat = FacStr * inv(FacCor); %compute factor pattern matrix
    
    LargestCom=1;
    LargestLoading=1;
    [V,L] = eig(S);
    L = flipud(fliplr(L));
    scree = diag(L);

else
    if ~strcmp(ROTATION,'IMAX')
        [V,L] = eig(S);
        
        if ~all(isreal(V))
            %Rounding errors resulted in relationship matrix being non-symmetric, which in turn resulted in complex numbers.
            %To fix, lower half of relationship matrix will be flipped to upper half to form symmetric relationship matrix.
            disp('Making relationship matrix symmetric to address effects of rounding errors.');
            Sfix=tril(S)+tril(S)'-diag(diag(S));
            [V,L] = eig(S);
        end;
        
        %since the output of eig is not always sorted in order of ascending size, first sort them
        [B,IX] = sort(diag(L),'descend');
        L=L(IX,IX);
        V=V(:,IX);
        
        V = V(:,1:NUM_FAC);  %truncated eigenvector matrix
        scree = diag(L);
        L = L(1:NUM_FAC,1:NUM_FAC);  %truncated eigenvalue matrix
        
        
        %factor scores
        if strcmp(MAT_TYPE,'SCP')
            FacScr = (data) * V;	%V is factor scoring coefficient matrix
        elseif strcmp(MAT_TYPE,'COV')
            FacScr = (data) * V;    %factor scores, not mean corrected.
        elseif strcmp(MAT_TYPE,'COR')
            FacScr = (data) * inv(D) * V; %same as above but also corrected for variable variance
        else
            error('Error - matrix type not specified');
        end
        
        ScrDiag = diag(std(FacScr));
        A = inv(Sd) * (V * ScrDiag);  %unrotated factor loading matrix
        %Note, this equation is commonly cited in statistics books but is misleading in this form.  The full
        %form is A = inv(Sd) * (inv(V') * sqrt (L))
        %this means that A is not in fact a scaled form of V as is commonly implied.
        %This makes sense since A is a factor loading matrix (converts scores to raw data)
        %while V is a scoring coefficient matrix (converts raw data to scores).
        %inv(V') does reduce down to V, though, since X'=inv(X)
        %for an orthogonal matrix, X.
        
        C = sum((A.^2),2)';
        switch LOADING
            case 'K'
                A = (diag(sqrt(C).^-1)) * A;  %factor loadings Kaiser-normalized by communalities
            case 'C'
                A = D *A;
            case 'N'
            case 'W'
                A = (diag(sqrt(C).^-1)) * A;  %factor loadings Kaiser-normalized by communalities
                facReflect=sign(A(:,1)) + (sign(A(:,1)) == 0);
                A = diag(facReflect)*A;
                wpt1 = acos(sqrt(1/NUM_FAC));
                halfPi=pi/2;
                for f3 = 1:size(A,1)
                    A1 = A(f3,1);
                    if A1 >= sqrt(1/NUM_FAC)%if the first factor loading is greater than the cos of the angle between the test vector and its first axis...
                        w1 = cos(((wpt1 - acos(A1))/wpt1) * halfPi)^2 + .001;
                    else
                        w1 = cos(((acos(A1) - wpt1)/(halfPi-wpt1)) * halfPi)^2 + .001;
                    end;
                    w(f3,1) = w1;
                end;
                A = diag(w)*A;
            otherwise
                error('Loading normalization needs to be specified as K or N or C or W');
        end;
    end;
    
    switch ROTATION
        case 'UNRT'
            FacPat=A;
            FacCor = diag(ones(NUM_FAC,1));
            FacStr = FacPat;
        case {'VMAX','PMAX'}
            [FacPat]= ep_doVarimax(A);
            FacCor = diag(ones(NUM_FAC,1));
            FacStr = FacPat;
        case 'VOMN'
            [FacPat, FacStr, FacCor] = ep_doVarOblimin(A);
        case 'IMAX'
            if size(data,1) == NUM_FAC
                msg{1}='Error: Infomax rotation code will not work when the number of observations exactly equals the number of retained factors.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            theRank=rank(data);
            if theRank < NUM_FAC
                disp(['Warning: the rank of the data ( ' num2str(theRank) ' ) is less than the number of retained variables ( ' num2str(NUM_FAC) ' ).']);
                disp('Infomax will tend to generate noisy factors when this is the case.  It is recommended that you set the number of factors to equal the rank.');
            end;
            rng(0,'twister'); %sets random number generator to standardized start to ensure replicability of ICA results.
            blocksize=ceil(min(5*log(size(data,1)),0.3*size(data,1))); %earlier versions of runICA needed fixing
            if blocksize < 2
                msg{1}=['Error: Too few observations (' num2str(size(data,1)) ') to conduct ICA.'];
                [msg]=ep_errorMsg(msg);
                return
            end;
            icadefs;
            if any(isspace(pwd)) && exist(ICABINARY,'file') == 2
                disp(['Binary ICA does not work if any part of the pathname (' pwd ') has a space in it.']);
                disp('Will use regular ICA instead.');
            end;
            
            ICAsuccess=0;
            if exist(ICABINARY,'file') == 2 && ~any(isspace(pwd)) %is there a binary file of runica installed?
                try
                    [weights,sphere] = binica(data', 'pca', NUM_FAC, 'verbose', 'on','block', blocksize,'filenum',9501);
                    delete('binica9501.wts');
                    delete('binica9501.sph');
                    delete('bias_after_adjust');
                    delete('binica9501.sc');
                    ICAsuccess=1;
                catch
                    ICAsuccess=0;
                    disp('Binary ICA failed to run.');
                end;
            end;
            if ~ICAsuccess
                try
                    warning('off')
                    [weights,sphere] = el_runica(data', 'pca', NUM_FAC, 'verbose', 'off','block', blocksize);
                    warning('on')
                catch
                    msg{1}='Error computing ICA.  The error message was:';
                    msg{2}=lasterr;
                    [msg]=ep_errorMsg(msg);
                    return
                end;
            end;
            
            if NUM_FAC > 1
                if all(diag(weights) == 1) || ~all(isreal(weights))
                    msg{1}='Error: ICA was unable to converge on a solution.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
            end;
            
            FacScr = (weights * sphere * data')';
            FacStr=corrcoef([FacScr data]); %compute factor structure matrix
            FacStr=FacStr(NUM_FAC+1:size(FacStr,1),1:NUM_FAC);
            FacCor = corrcoef(FacScr);
            FacPat = FacStr * inv(FacCor); %compute factor pattern matrix
            varSD=std(data);
            [V,L] = eig(S);
            L = flipud(fliplr(L));
            scree = diag(L); %scree is based on the factors from the initial unrotated PCA.  Since the PCA Toolbox uses
            %the assumption that PCA is always used as a first step with ICA, this is the appropriate scree.
        otherwise
            [FacPat, FacStr, FacCor] = ep_doGPFrotation(A, ROTATION, ROTOPT);
    end;
    
    if ~strcmp(ROTATION,'IMAX')
        switch LOADING
            case 'K'
                FacPat = diag(sqrt(C)) * FacPat;  %renormalize factor loadings by original communalities
                FacStr = diag(sqrt(C)) * FacStr;
            case 'C'
                FacPat = inv(D) * FacPat;
                FacStr = inv(D) * FacStr;
            case 'N'
            case 'W'
                FacPat = inv(diag(w))*FacPat;
                FacPat = diag(facReflect)*FacPat;
                FacPat = diag(sqrt(C)) * FacPat;  %renormalize factor loadings by original communalities
                FacStr = inv(diag(w))*FacPat;
                FacStr = diag(facReflect)*FacPat;
                FacStr = diag(sqrt(C)) * FacPat;  %renormalize factor loadings by original communalities
            otherwise
                error('Loading normalization needs to be specified as K or N or C or W');
        end;
    end;
    
    if strcmp(ROTATION,'PMAX')
        [FacPat, FacCor] = ep_doPromax(FacPat, ROTOPT); %Only apply loading weighting to the Varimax step
        %to match SAS output and to avoid rounding errors.
        FacStr = FacPat * FacCor;	%factor structure matrix (Harman, eq. 12.19, p. 268)
    end;
    
    if NUM_FAC == 1 %ensure that matrices have the right orientation when there is only one factor.
        if size(FacPat,1) < size(FacPat,2)
            FacPat=FacPat';
        end;
        if size(FacStr,1) < size(FacStr,2)
            FacStr=FacStr';
        end;
    end;
    
    %Deal with loadings that are too large
    LargestLoading=max(max(abs(FacStr))); %factor pattern loadings can go over 1 for oblique rotations
    if round(LargestLoading*100) > 100 %allow very small violation of factor loading limit due to rounding errors
        msg{1}=['Loadings are over the permitted maximum of 1.  It appears this rotation has crashed.'];
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    LargestCom=max(max(abs(sum(FacPat.*FacStr,2))));
    if round(LargestCom*100) > 100 %allow very small violation of communality limit due to rounding errors
        msg{1}=['Communalities are over the permitted maximum of 1.  It appears this rotation has crashed.'];
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    invR=pinv(R);
    FacCof=invR*FacStr;
    FacScr=(data)*FacCof;
end;

facVar=[];
facVarQ=[];
Comm=[];
    
if LargestCom ~=0 && LargestLoading~=0

    FacScr=(FacScr)*inv(diag(std(FacScr))); %Standardize factor scores, not mean corrected.

    Var=Sd.^2;

    %calculate communalities (equivalent to calculating R-Squared in multiple regression).

    %Cohen, J. & Cohen, P. (1983). Applied multiple regression/correlation
    %analysis for the behavioral sciences. Hillsdale, NJ: Lawrence Erlbaum Associates.

    Comm=sum((Var*FacPat.*FacStr),2)/sum(diag(Var)); % p. 100

%     if strcmp(theRotationName,'none')
%         theRotationName2='unrotated';
%     else
%         theRotationName2=theRotationName;
%     end;
%     
%     disp(['Amount of variance accounted for by the ' theRotationName2 ' solution is ' sprintf('%3.2f%%.',100*sum(Comm,1))]);

    facVar = sum((Var*FacPat.*FacStr))/sum(diag(Var));

    %calculate unique variance of factors by calculating semi-partial correlation
    %between the factor and the variable.  p. 470 of C&C (1983)
    
    facVarQ=sum(Var*(FacPat*inv(diag(sqrt(diag(inv(FacCor)))))).^2)/sum(diag(Var));
    
    if GAVE == 'Y'
        FacCof = pinv(corrcoef(data))*FacStr;
        FacScr = rawData*FacCof;
    end;
    
    if isempty(crossFacCof)
        [dummy index]=sort(facVar); %sort factors in order of size
        index = fliplr(index);
        
        FacPat=FacPat(:,index);
        FacStr=FacStr(:,index);
        FacCof=FacCof(:,index);
        FacScr=FacScr(:,index);
        FacCor=FacCor(index,index);
        facVar=facVar(index);
        facVarQ=facVarQ(index);
        
        for f1 = 1:NUM_FAC
            if sum(FacPat(:,f1)) < 0	%flip factor loadings so mostly positive
                FacPat(:,f1) = FacPat(:,f1) .* (-1);
                FacStr(:,f1) = FacStr(:,f1) .* (-1);
                FacCof(:,f1) = FacCof(:,f1) .* (-1);
                FacScr(:,f1) = FacScr(:,f1) .* (-1); %if the loading is flipped, the scores must be flipped too.
                FacCor(:,f1) = FacCor(:,f1) .* (-1);
                FacCor(f1,:) = FacCor(f1,:) .* (-1);
            end;
        end;
    end;
end;

if isstruct(theData)
    if ~isfield(theData,'montage');
        theData.montage=[];
    end;
    if ~isfield(theData,'chanNames');
        theData.chanNames=[];
    end;
    if ~isfield(theData,'timeNames');
        theData.timeNames=[];
    end;
    if ~isfield(theData,'subNames');
        theData.subNames=[];
    end;
    if ~isfield(theData,'cellNames');
        theData.cellNames=[];
    end;
    if ~isfield(theData,'freqNames');
        theData.freqNames=[];
    end;
    if ~isfield(theData,'Fs');
        theData.Fs=[];
    end;
    if ~isfield(theData,'baseline');
        theData.baseline=[];
    end;
    if ~isfield(theData,'ename');
        theData.ename=[];
    end;
end;

switch PCAmode
    case 'spat'
        theLabel='SF';
    case 'temp'
        theLabel='TF';
    case 'freq'
        theLabel='FF';
end;

facNames=[];
facTypes=[];
if ~strcmp(PCAmode,'asis')
    factorDigits=length(num2str(NUM_FAC));
    for i=1:NUM_FAC
        facNames{i,1}=[theLabel sprintf(['%0' num2str(factorDigits) 'd'],i)];
        facTypes{i,1}='SGL';
    end;
end;

FactorResults.FacPat=zeros(size(badData,2),size(FacPat,2));
FactorResults.FacScr=zeros(size(badData,1),size(FacScr,2));
FactorResults.FacPat(goodVars,:)=FacPat;
FactorResults.FacStr=zeros(size(badData,2),size(FacStr,2));
FactorResults.FacStr(goodVars,:)=FacStr;
FactorResults.FacScr=zeros(size(dataAll,1),NUM_FAC);
if isComplex
    FactorResults.FacScr(goodObs,:)=complex(FacScr(1:size(FacScr,1)/2,:),FacScr(size(FacScr,1)/2+1:end,:));
else 
    FactorResults.FacScr(goodObs,:)=FacScr;
end;
FactorResults.FacCor=FacCor;
FactorResults.FacCof=zeros(size(badData,2),size(FacCof,2));
FactorResults.FacCof(goodVars,:)=FacCof;
FactorResults.scree=zeros(size(badData,2),1);
FactorResults.scree(1:length(scree))=scree;
FactorResults.facVar=facVar;
FactorResults.facVarQ=facVarQ;
FactorResults.facVarTot=sum(Comm,1);
FactorResults.varSD=zeros(size(badData,2),1);
FactorResults.varSD(goodVars)=varSD;
FactorResults.PCAmode=PCAmode;
FactorResults.ROTATION=ROTATION;
FactorResults.MAT_TYPE=MAT_TYPE;
FactorResults.LOADING=LOADING;
FactorResults.RotOpt=ROTOPT;
FactorResults.timepoints=numPoints;
FactorResults.numchan=numChans;
FactorResults.numFreqs=numFreqs;
FactorResults.numRels=numRels;
FactorResults.numSubs=numSubs;
FactorResults.numCells=numCells;
FactorResults.numFacs=NUM_FAC;
FactorResults.facNames=facNames;
FactorResults.facTypes=facTypes;
try
    EPver=ver('EP_Toolkit');
catch
    EPver='unavailable'; %workaround for bug in earlier version of Matlab
end;
FactorResults.EPver=EPver;
FactorResults.ver=ver;
FactorResults.date=date;
FactorResults.history={'doPCA',PCAmode, ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, LOADING, GAVE};
FactorResults.badObs=badObs;

if isstruct(theData)
    FactorResults.montage=theData.montage;
    FactorResults.eloc=theData.eloc;
    FactorResults.ced=theData.ced;
    FactorResults.implicit=theData.implicit;
    FactorResults.chanNames=theData.chanNames;
    FactorResults.timeNames=theData.timeNames;
    FactorResults.subNames=theData.subNames;
    FactorResults.cellNames=theData.cellNames;
    FactorResults.trialNames=theData.trialNames;
    FactorResults.freqNames=theData.freqNames;
    FactorResults.relNames=theData.relNames;
    FactorResults.chanTypes=theData.chanTypes;
    FactorResults.timeUnits=theData.timeUnits;
    FactorResults.subTypes=theData.subTypes;
    FactorResults.cellTypes=theData.cellTypes;
    FactorResults.subjectSpecs=theData.subjectSpecs;
    FactorResults.subjectSpecNames=theData.subjectSpecNames;
    FactorResults.Fs=theData.Fs;
    FactorResults.baseline=theData.baseline;
    FactorResults.ename=theData.ename;
    FactorResults.dataType=theData.dataType;
    FactorResults.fileName=theData.fileName;
    FactorResults.fileFormat=theData.fileFormat;
    FactorResults.analysis.blinkTrial=theData.analysis.blinkTrial;
    FactorResults.analysis.saccadeTrial=theData.analysis.saccadeTrial;
    FactorResults.analysis.saccadeOnset=theData.analysis.saccadeOnset;
    FactorResults.analysis.moveTrial=theData.analysis.moveTrial;
    FactorResults.analysis.badTrials=theData.analysis.badTrials;
    FactorResults.analysis.badChans=theData.analysis.badChans;
    FactorResults.reference=theData.reference;
    FactorResults.recTime=theData.recTime;
    
    if isfield(theData,'noise')
        FactorResults.noise=theData.noise;
    end;
    if isfield(theData,'events')
        FactorResults.events=theData.events;
    end;
    if isfield(theData,'stims')
        FactorResults.stims=theData.stims;
    end;
    if isfield(theData,'calibration')
        FactorResults.calibration=theData.calibration;
    end;
    if isfield(theData,'impedances')
        FactorResults.impedances=theData.impedances;
    end;
    if isfield(theData,'avgNum')
        FactorResults.avgNum=theData.avgNum;
    end;
    if isfield(theData,'covNum')
        FactorResults.covNum=theData.covNum;
    end;
    if isfield(theData,'subNum')
        FactorResults.subNum=theData.subNum;
    end;
    if isfield(theData,'trialSpecs')
        FactorResults.trialSpecs=theData.trialSpecs;
    end;
    if isfield(theData,'trialSpecNames')
        FactorResults.trialSpecNames=theData.trialSpecNames;
    end;

end;
