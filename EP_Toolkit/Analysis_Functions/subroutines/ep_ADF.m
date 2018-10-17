function ep_ADF(Y, NX, OPT1, PER, OPT2, NUMSIM, REPS, SEED, MISSING, factorNames, levelNames, factorGroupNames, levelGroupNames, elecFactors, alpha, followup, OPT3, SCALE, LOC1, LOC2, outfid, followupWithin, followupBetween);
%function ep_ADF(Y, NX, OPT1, PER, OPT2, NUMSIM, REPS, SEED, MISSING, factorNames, levelNames, factorGroupNames, levelGroupNames, elecFactors, alpha, followup, OPT3, SCALE, LOC1, LOC2, outfid, followupWithin, followupBetween);
%ADF algorithms made available by Harvey Keselman (kesel@Ms.UManitoba.CA)
%Department of Psychology, University of Manitoba, Winnipeg, Canada.
%http://www.umanitoba.ca/cgi-bin/psychology/hpg_main.cgi?data/Keselman.txt
% 
%The ADF function is a front end for the WJGLMml function.  It
%automatically generates the appropriate contrasts for up to six within and
%six between participant factors.  It also computes the averaged trimmed
%cell means for each of the contrasts.  Finally, it presents the full set
%of tests, along with the appropriate cell means.
%
%Front end made by Joseph Dien (jdien07@mac.com)
%
%Keselman, H. J., Wilcox, R. R., & Lix, L. M. (2003). A generally robust approach
%to hypothesis testing in independent and correlated groups designs
%Psychophysiology, 40, 586-596.
%
%McCarthy, G., & Wood, C. C. (1985). Scalp distribution of event-related potentials: An ambiguity associated with
%analysis of variance models. Electroencephalography and Clinical Neurophysiology, 62, 203-208.

%Inputs
%  Y	    : Input data (rows=subjects, columns = cells).  The first NX rows will be assigned to the first group and so forth.
%  NX   	: Number of subjects in each group.  Set to empty set if no between group factors.
%  OPT1     : Activate trimming option.  0 = no and 1 = yes.
%  PER  	: Percentage to trim the means.  .05 is the recommended number for ERP data.
%  OPT2 	: Activate Welch-James statistic.  0 = no and 1 = yes.  "No" option is not currently available.
%  NUMSIM	: Number of simulations used to generate bootstrapping statistic.  p-values will be unstable if too low. 4999 informally recommended.
%  REPS     : Number of repetitions to determine degree of p-value variability.  Eleven informally recommended.
%  SEED 	: Seed for random number generation.  0 specifies random seed. 1000 arbitrarily suggested as seed to ensure results are replicable.
%  MISSING  : Number to be treated as a missing value.  Observations with missing values are dropped from the analysis.
%  factorNames: The names of the within factors.  Each should be three letters long.
%               The order should correspond to the structure of the data with the first varying the slowest.
%               Set to empty set if no within factors.
%  levelNames : The names of the levels within each within factor.  One row for each factor.  Each level name should be a single letter.
%               The order should correspond to the order of the levels in the data.
%  factorGroupNames: The names of the between factors.  Each should be three letters long.
%                    The order should correspond to the structure of the data with the first varying the slowest.
%  levelGroupNames : The names of the levels within each between factor.  One row for each factor.  Each level name should be a single letter.
%                    The order should correspond to the order of the levels in the data.
%  elecFactors     : An array indicating which within group factors are electrode factors and hence need the vector
%                    scaling test.  0=not electrode factor, 1=electrode factor.  To disable, just set all to zeroes.
%  alpha :   The alpha thresholds for significance
%    .uncorrected:  Alpha for a priori interest (.05)
%    .corrected:  Alpha for posthoc multiple comparison corrected.
%  followup: How many levels deep a follow-up this is with zero being not a follow-up.
%  OPT3     : Provide effect sizes. 0 = no and 1 = yes.
%  SCALE    : "SCALE is a scalar indicator to control the use of a scaling factor for the effect size estimator 
%             (i.e., .642 for 20% symmetric trimming) when robust estimators are adopted. 
%             It takes a value of 0 or 1; a zero indicates that no scaling factor will be used, 
%             while a 1 indicates that a scaling factor will be adopted. The default is SCALE=1.
%  LOC1&LOC2: "If the user specifies LOC1 = 0 and LOC2 = 0, the square root of the average of the variances over the cells 
%             involved in the contrast is used as the standardizer. If the
%             user specifies LOC1 = 99 and LOC2 = 99, no standardizer is selected."
%  outfid:  The fid of the output file.  Will request output file if not provided.
%  followupWithin: The within group factors being followed up on.
%  followupBetween: The between group factors being followed up on.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bugfix 3/25/08 JD
% Fix for p-value bug where it would output just "0" if there were no zeroes in the decimal.
%
% modified 10/29/08 JD
% Less decimals for DF and test statistic outputs
%
% modified 11/8/08 JD
% When p-value less than .00000001, print out less than sign.  Improved statistics output.
%
% modified 3/25/09 JD
% changed level names to cell array.  Added alpha criteria variable so can highlight results.
% saves output file in html format.
%
% modified 7/24/09 JD
% Added testing of electrode factor effects using the McCarthy & Woods (1985) vector scaling test.
% Automatic follow-up tests to significant interactions
%
% bugfix 11/6/09 JD
% Fixes to spacing of cell mean output.
%
% bugfix 2/6/10 JD
% Fixed bug in code for follow-up ANOVAs that could cause crashes or the wrong factors to be used (but still correctly
% labeled).
%
% bugfix 10/23/14 JD
% Fixed crash when running ANOVA with only between-group factors.
%
% modified 9/20/15 JD
% Added effect sizes.
%
% modified 12/11/15 JD
% Checks if p-value variability exceeds two standard deviations over the threshold.
% Includes warning summary of %age of singular and nearly singular matrices in the output.
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

if nargin < 17
    OPT3 =1;
end;

if nargin < 18
    SCALE =0;
end;

if nargin < 19
    outfid =0;
end;

if nargin < 20
    followupWithin =[];
end;

if nargin < 21
    followupBetween =[];
end;

if ischar(outfid)
    msg{1}=[outfid ' is a string, not a file ID number.'];
    [msg]=ep_errorMsg(msg);
    return
end;

for i=1:length(factorNames)
    if length(factorNames{i}) ~= 3
        msg{1}=['The factor name ' factorNames{i} ' is not three letters long.'];
        [msg]=ep_errorMsg(msg);
        return
    end;
end;

for i=1:length(factorGroupNames)
    if length(factorGroupNames{i}) ~= 3
        msg{1}=['The factor name ' factorGroupNames{i} ' is not three letters long.'];
        [msg]=ep_errorMsg(msg);
        return
    end;
end;

if isempty(elecFactors)
    elecFactors=zeros(length(factorNames));
end;

if isempty(NX)
    NX = size(Y,1);
end;

if size(NX,1) > size(NX,2)
    NX=NX';
end;

varNum = size(Y,2);
obsNum = size(Y,1);
grpNum = length(NX);
maxDigits = length(num2str(floor(max(max(Y)))))+2;

if ~isempty(MISSING) %calculate how many observations with missing data points
    good = [];
    newNX= zeros(1,length(NX));
    count=0;
    for group = 1:length(NX)
        for obs = 1:NX(group)
            count=count+1;
            if isempty(find(Y(count,:)==MISSING))
                good=[good count];
                newNX(group)=newNX(group)+1;
            end;
        end;
    end;
    if ~isempty(find(newNX==0))
        msg{1}='A group has zero members after dropping missing data points';
        [msg]=ep_errorMsg(msg);
        return
    end;
end;

if outfid == 0
    fileName=uiputfile('results.html','File name for output');
    if fileName == 0
        msg{1}='Something went wrong.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    fileName=[pwd filesep fileName];
    outfid=fopen(fileName,'w');
    internalFID = 1;
else
    internalFID = 0;
end;

if outfid == -1
    msg{1}='Results output file not successfully created.';
    [msg]=ep_errorMsg(msg);
    return
end;

%determine which columns go with which withingroup levels

repFactor=1;
repCount=1;
lvlCount=1;
numFactors=length(levelNames);
factorTable=cell(numFactors,varNum);
for theFactor=numFactors:-1:1
    for col=1:varNum
        factorTable{theFactor,col}=levelNames{theFactor}(lvlCount);
        repCount=repCount+1;
        if repCount > repFactor
            repCount=1;
            lvlCount=lvlCount+1;
            if lvlCount > length(levelNames{theFactor})
                lvlCount=1;
            end;
        end;
    end;
    repFactor=repFactor*length(levelNames{theFactor});
end;

%determine which columns go with which betweengroup levels

repFactor=1;
repCount=1;
lvlCount=1;
numGroupFactors=length(levelGroupNames);
factorGroupTable=cell(grpNum,numGroupFactors);
for theFactor=numGroupFactors:-1:1
    for row=1:grpNum
        factorGroupTable{row,theFactor}=levelGroupNames{theFactor}(lvlCount);
        repCount=repCount+1;
        if repCount > repFactor
            repCount=1;
            lvlCount=lvlCount+1;
            if lvlCount > length(levelGroupNames{theFactor})
                lvlCount=1;
            end;
        end;
    end;
    repFactor=repFactor*length(levelGroupNames{theFactor});
end;

%%%%%%compute within group contrasts%%%%%
facNum = length(factorNames);
levelNum = zeros(facNum,1);

for i = 1:facNum
            levelNum(i)=length(levelNames{i});
end;

if varNum ~= prod(levelNum)
    msg{1}=['Number of variables (' num2str(varNum) ') must equal product of within-group factor levels (' num2str(prod(levelNum)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end;

if length(factorNames)~= length(levelNames) && ~isempty(factorNames)
    msg{1}=['Number of within factors given names (' num2str(length(factorNames)) ') must equal number of within factors given level names (' num2str(length(levelNames)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end;

if sum(NX) ~= size(Y,1)
    msg{1}=['Number of observations (' num2str(size(Y,1)) ') must equal the sum of the between group sizes (' num2str(sum(NX)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end;

if length(factorGroupNames)~= length(levelGroupNames) && ~isempty(factorGroupNames)
    msg{1}=['Number of between factors given names (' num2str(length(factorGroupNames)) ') must equal number of between factors given level names (' num2str(length(levelGroupNames)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end;

if facNum > 6
    msg{1}='Program only handles up to six factors.';
    [msg]=ep_errorMsg(msg);
    return
end;

contrast = [];
basicContrast = [];
facCombs = zeros(facNum,facNum);
contrast{1} = ones(1,varNum);
contrastName = 'NO WITHIN EFFECTS';
contrastLevelNames{1} = ' ';
contrastFactors=zeros(1,facNum); %list of which factors are participating in a contrast

for i = 1:facNum
    facCombs(i,i)=1;
    contrastName = char(contrastName, [factorNames{i} ' MAIN EFFECT']);
    contrastFactors(size(contrastName,1),i)=1;
    basicContrast{i} = [ones(levelNum(i)-1,1) -eye(levelNum(i)-1)];
    if isempty(basicContrast{i})
        basicContrast{i}=1;
    end;
    if i ~= 1
        contrast{i+1}  = ones(1,levelNum(1));
    else
        contrast{i+1}  = basicContrast{1};
    end;
    
    for factor = 2:facNum
        if i ~= factor
            contrast{i+1} = kron(contrast{i+1}, ones(1,levelNum(factor)));
        else
            contrast{i+1} = kron(contrast{i+1}, basicContrast{factor});
        end;
    end;
end;

[MUHAT, SIGMA, RESULTS]=ep_WJGLMml(Y, NX, ones(1,grpNum), contrast{1}', OPT1, PER, OPT2, 1, SEED, MISSING, OPT3, alpha, SCALE, LOC1, LOC2); %generate MUHAT
if isempty(RESULTS)
    return
end;

cellMeans=zeros(grpNum,1,varNum); %means of the cells with groups by contrast by levels in the contrast
cellMeans(:,1,1)=mean(reshape(MUHAT',varNum,grpNum),1);

for factor = 1:facNum
    theNames = [];
    for level = 1:levelNum(factor)
        theNames = [theNames; levelNames{factor}(level)];
    end
    contrastLevelNames{factor+1} = theNames;
    theContrast = contrast{factor+1};
    for group = 1:grpNum
        theMUHAT = MUHAT((group-1)*varNum+1:group*varNum,1);
        cellMeans1 = sum(((theContrast(1,:) == 1)*diag(theMUHAT))')/sum(theContrast(1,:) == 1);   %add together the cell means for the first condition
        cellMeans2 = sum(((theContrast == -1)*diag(theMUHAT))')/sum(theContrast(1,:) == -1); %generate separate cell mean sums for each of the other conditions
        cellMeans3 = [cellMeans1 cellMeans2];
        cellMeans(group,factor+1,1:length(cellMeans3)) = cellMeans3;
    end;
end;

if facNum > 1
    for i1 = 1:(facNum-1)
        for i2 = (i1 + 1):facNum
            facCombs = [facCombs; zeros(1,facNum)];
            facCombs(size(facCombs,1),i1) = 1;
            facCombs(size(facCombs,1),i2) = 1;
            contrastName = char(contrastName, [factorNames{i1} ' * ' factorNames{i2} ' INTERACTION EFFECT']);
            contrastFactors(size(contrastName,1),[i1 i2])=1;
            
            contrastNumber = size(contrast,2)+1;
            if facCombs(size(facCombs,1),1) == 0
                contrast{contrastNumber}  = ones(1,levelNum(1));
            else
                contrast{contrastNumber}  = basicContrast{1};
            end;
            
            for factor = 2:facNum
                if facCombs(size(facCombs,1),factor) == 0
                    contrast{contrastNumber} = kron(contrast{contrastNumber}, ones(1,levelNum(factor)));
                else
                    contrast{contrastNumber} = kron(contrast{contrastNumber}, basicContrast{factor});
                end;
            end;
            
            theContrast1 = contrast{i1+1};
            cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
            theContrast2 = contrast{i2+1};
            cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
            
            theNames = [];
            for level1 = 1:levelNum(i1)
                for level2 = 1:levelNum(i2)
                    theNames = [theNames; [levelNames{i1}(level1) levelNames{i2}(level2)]];
                    for group = 1:grpNum
                        theCellMeans{group}=[];
                        theMUHAT = MUHAT((group-1)*varNum+1:group*varNum,1);
                        theCell = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*theMUHAT';
                        cellMeans(group,contrastNumber,size(theNames,1)) = sum(theCell)/sum(theCell ~= 0);
                    end;
                end;
            end
            contrastLevelNames{contrastNumber} = theNames;
        end;
    end;
end;


if facNum > 2
    for i1 = 1:(facNum-2)
        for i2 = (i1 + 1):(facNum-1)
            for i3 = (i2 + 1):facNum
                facCombs = [facCombs; zeros(1,facNum)];
                facCombs(size(facCombs,1),i1) = 1;
                facCombs(size(facCombs,1),i2) = 1;
                facCombs(size(facCombs,1),i3) = 1;
                contrastName = char(contrastName, [factorNames{i1} ' * ' factorNames{i2} ' * ' factorNames{i3} ' INTERACTION EFFECT']);
                contrastFactors(size(contrastName,1),[i1 i2 i3])=1;
                
                contrastNumber = size(contrast,2)+1;
                if facCombs(size(facCombs,1),1) == 0
                    contrast{contrastNumber}  = ones(1,levelNum(1));
                else
                    contrast{contrastNumber}  = basicContrast{1};
                end;
                
                for factor = 2:facNum
                    if facCombs(size(facCombs,1),factor) == 0
                        contrast{contrastNumber} = kron(contrast{contrastNumber}, ones(1,levelNum(factor)));
                    else
                        contrast{contrastNumber} = kron(contrast{contrastNumber}, basicContrast{factor});
                    end;
                end;
                
                theContrast1 = contrast{i1+1};
                cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                theContrast2 = contrast{i2+1};
                cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                theContrast3 = contrast{i3+1};
                cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                
                theCellMeans = [];
                theNames = [];
                for level1 = 1:levelNum(i1)
                    for level2 = 1:levelNum(i2)
                        for level3 = 1:levelNum(i3)
                            theNames = [theNames; [levelNames{i1}(level1) levelNames{i2}(level2) levelNames{i3}(level3)]];
                            for group = 1:grpNum
                                theCellMeans{group}=[];
                                theMUHAT = MUHAT((group-1)*varNum+1:group*varNum,1);
                                theCell = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:).*theMUHAT';
                                cellMeans(group,contrastNumber,size(theNames,1)) = sum(theCell)/sum(theCell ~= 0);
                            end;
                        end;
                    end
                end;
                contrastLevelNames{contrastNumber} = theNames;
            end;
        end;
    end;
end;


if facNum > 3
    for i1 = 1:(facNum-3)
        for i2 = (i1 + 1):(facNum-2)
            for i3 = (i2 + 1):(facNum-1)
                for i4 = (i3 + 1):facNum
                    facCombs = [facCombs; zeros(1,facNum)];
                    facCombs(size(facCombs,1),i1) = 1;
                    facCombs(size(facCombs,1),i2) = 1;
                    facCombs(size(facCombs,1),i3) = 1;
                    facCombs(size(facCombs,1),i4) = 1;
                    contrastName = char(contrastName, [factorNames{i1} ' * ' factorNames{i2} ' * ' factorNames{i3} ' * ' factorNames{i4} ' INTERACTION EFFECT']);
                    contrastFactors(size(contrastName,1),[i1 i2 i3 i4])=1;
                    
                    contrastNumber = size(contrast,2)+1;
                    if facCombs(size(facCombs,1),1) == 0
                        contrast{contrastNumber}  = ones(1,levelNum(1));
                    else
                        contrast{contrastNumber}  = basicContrast{1};
                    end;
                    
                    for factor = 2:facNum
                        if facCombs(size(facCombs,1),factor) == 0
                            contrast{contrastNumber} = kron(contrast{contrastNumber}, ones(1,levelNum(factor)));
                        else
                            contrast{contrastNumber} = kron(contrast{contrastNumber}, basicContrast{factor});
                        end;
                    end;
                    
                    theContrast1 = contrast{i1+1};
                    cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                    theContrast2 = contrast{i2+1};
                    cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                    theContrast3 = contrast{i3+1};
                    cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                    theContrast4 = contrast{i4+1};
                    cellMeansKey4 = [(theContrast4(1,:) == 1); (theContrast4 == -1)];
                    
                    theCellMeans = [];
                    theNames = [];
                    for level1 = 1:levelNum(i1)
                        for level2 = 1:levelNum(i2)
                            for level3 = 1:levelNum(i3)
                                for level4 = 1:levelNum(i4)
                                    theNames = [theNames; [levelNames{i1}(level1) levelNames{i2}(level2) levelNames{i3}(level3) levelNames{i4}(level4)]];
                                    for group = 1:grpNum
                                        theCellMeans{group}=[];
                                        theMUHAT = MUHAT((group-1)*varNum+1:group*varNum,1);
                                        theCell = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:).*cellMeansKey4(level4,:).*theMUHAT';
                                        cellMeans(group,contrastNumber,size(theNames,1)) = sum(theCell)/sum(theCell ~= 0);
                                    end;
                                end;
                            end
                        end;
                    end;
                    contrastLevelNames{contrastNumber} = theNames;
                end;
            end;
        end;
    end;
end;


if facNum > 4
    for i1 = 1:(facNum-4)
        for i2 = (i1 + 1):(facNum-3)
            for i3 = (i2 + 1):(facNum-2)
                for i4 = (i3 + 1):(facNum-1)
                    for i5 = (i4 + 1):facNum
                        facCombs = [facCombs; zeros(1,facNum)];
                        facCombs(size(facCombs,1),i1) = 1;
                        facCombs(size(facCombs,1),i2) = 1;
                        facCombs(size(facCombs,1),i3) = 1;
                        facCombs(size(facCombs,1),i4) = 1;
                        facCombs(size(facCombs,1),i5) = 1;
                        contrastName = char(contrastName, [factorNames{i1} ' * ' factorNames{i2} ' * ' factorNames{i3} ' * ' factorNames{i4} ' * ' factorNames{i5} ' INTERACTION EFFECT']);
                        contrastFactors(size(contrastName,1),[i1 i2 i3 i4 i5])=1;
                        
                        contrastNumber = size(contrast,2)+1;
                        if facCombs(size(facCombs,1),1) == 0
                            contrast{contrastNumber}  = ones(1,levelNum(1));
                        else
                            contrast{contrastNumber}  = basicContrast{1};
                        end;
                        
                        for factor = 2:facNum
                            if facCombs(size(facCombs,1),factor) == 0
                                contrast{contrastNumber} = kron(contrast{contrastNumber}, ones(1,levelNum(factor)));
                            else
                                contrast{contrastNumber} = kron(contrast{contrastNumber}, basicContrast{factor});
                            end;
                        end;
                        
                        theContrast1 = contrast{i1+1};
                        cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                        theContrast2 = contrast{i2+1};
                        cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                        theContrast3 = contrast{i3+1};
                        cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                        theContrast4 = contrast{i4+1};
                        cellMeansKey4 = [(theContrast4(1,:) == 1); (theContrast4 == -1)];
                        theContrast5 = contrast{i5+1};
                        cellMeansKey5 = [(theContrast5(1,:) == 1); (theContrast5 == -1)];
                        
                        
                        theCellMeans = [];
                        theNames = [];
                        for level1 = 1:levelNum(i1)
                            for level2 = 1:levelNum(i2)
                                for level3 = 1:levelNum(i3)
                                    for level4 = 1:levelNum(i4)
                                        for level5 = 1:levelNum(i5)
                                            theNames = [theNames; [levelNames{i1}(level1) levelNames{i2}(level2) levelNames{i3}(level3) levelNames{i4}(level4) levelNames{i5}(level5)]];
                                            for group = 1:grpNum
                                                theCellMeans{group}=[];
                                                theMUHAT = MUHAT((group-1)*varNum+1:group*varNum,1);
                                                theCell = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:).*cellMeansKey4(level4,:).*cellMeansKey5(level5,:).*theMUHAT';
                                                cellMeans(group,contrastNumber,size(theNames,1)) = sum(theCell)/sum(theCell ~= 0);
                                            end;                                            
                                        end;
                                    end
                                end;
                            end;
                        end;
                        contrastLevelNames{contrastNumber} = theNames;
                    end;
                end;
            end;
        end;
    end;
end;


if facNum > 5
    for i1 = 1:(facNum-5)
        for i2 = (i1 + 1):(facNum-4)
            
            for i3 = (i2 + 1):(facNum-3)
                for i4 = (i3 + 1):(facNum-2)
                    for i5 = (i4 + 1):(facNum-1)
                        for i6 = (i5 + 1):facNum
                            facCombs = [facCombs; zeros(1,facNum)];
                            facCombs(size(facCombs,1),i1) = 1;
                            facCombs(size(facCombs,1),i2) = 1;
                            facCombs(size(facCombs,1),i3) = 1;
                            facCombs(size(facCombs,1),i4) = 1;
                            facCombs(size(facCombs,1),i5) = 1;
                            facCombs(size(facCombs,1),i6) = 1;
                            contrastName = char(contrastName, [factorNames{i1} ' * ' factorNames{i2} ' * ' factorNames{i3} ' * ' factorNames{i4} ' * ' factorNames{i5} ' * ' factorNames{i6} ' INTERACTION EFFECT']);
                            contrastFactors(size(contrastName,1),[i1 i2 i3 i4 i5 i6])=1;
                            
                            contrastNumber = size(contrast,2)+1;
                            if facCombs(size(facCombs,1),1) == 0
                                contrast{contrastNumber}  = ones(1,levelNum(1));
                            else
                                contrast{contrastNumber}  = basicContrast{1};
                            end;
                            
                            for factor = 2:facNum
                                if facCombs(size(facCombs,1),factor) == 0
                                    contrast{contrastNumber} = kron(contrast{contrastNumber}, ones(1,levelNum(factor)));
                                else
                                    contrast{contrastNumber} = kron(contrast{contrastNumber}, basicContrast{factor});
                                end;
                            end;
                            
                            theContrast1 = contrast{i1+1};
                            cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                            theContrast2 = contrast{i2+1};
                            cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                            theContrast3 = contrast{i3+1};
                            cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                            theContrast4 = contrast{i4+1};
                            cellMeansKey4 = [(theContrast4(1,:) == 1); (theContrast4 == -1)];
                            theContrast5 = contrast{i5+1};
                            cellMeansKey5 = [(theContrast5(1,:) == 1); (theContrast5 == -1)];
                            theContrast6 = contrast{i6+1};
                            cellMeansKey6 = [(theContrast6(1,:) == 1); (theContrast6 == -1)];
                            
                            theCellMeans = [];
                            theNames = [];
                            for level1 = 1:levelNum(i1)
                                for level2 = 1:levelNum(i2)
                                    for level3 = 1:levelNum(i3)
                                        for level4 = 1:levelNum(i4)
                                            for level5 = 1:levelNum(i5)
                                                for level6 = 1:levelNum(i6)
                                                    theNames = [theNames; [levelNames{i1}(level1) levelNames{i2}(level2) levelNames{i3}(level3) levelNames{i4}(level4) levelNames{i5}(level5) levelNames{i6}(level6)]];
                                                    for group = 1:grpNum
                                                        theCellMeans{group}=[];
                                                        theMUHAT = MUHAT((group-1)*varNum+1:group*varNum,1);
                                                        theCell = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:).*cellMeansKey4(level4,:).*cellMeansKey5(level5,:).*cellMeansKey6(level6,:).*theMUHAT';
                                                        cellMeans(group,contrastNumber,size(theNames,1)) = sum(theCell)/sum(theCell ~= 0);
                                                    end;                                            
                                                end;
                                            end
                                        end;
                                    end;
                                end;
                            end;
                            contrastLevelNames{contrastNumber} = theNames;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

%%%%%%compute between group contrasts%%%%%

groupContrast{1} = ones(1,grpNum);
contrastGroupName = 'NO BETWEEN EFFECTS';
betweenContrastFactors=zeros(1,grpNum); %list of which factors are participating in a contrast
contrastLevelGroupNames{1} = ' ';
cellGroupMeans(1,:,:,:)=mean(cellMeans,1);%means of the cells with groupscontrast by levels in it by withincontrast by levels in the contrast

if grpNum >1
    
    facGroupNum = length(factorGroupNames);
    levelGroupNum = zeros(facGroupNum,1);
    contrastNum = size(cellMeans,2);
    
    for i = 1:facGroupNum
                levelGroupNum(i)=length(levelGroupNames{i});
    end;
    
    if grpNum ~= prod(levelGroupNum)
        msg{1}='Number of groups must equal product of between factor levels.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if facGroupNum > 6
        msg{1}='Program only handles up to six factors.  You''d go insane trying to interpet that many levels of interactions anyway.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    basicContrast = [];
    facCombs = zeros(facGroupNum,facGroupNum);
    for i = 1:facGroupNum
        facCombs(i,i)=1;
        contrastGroupName = char(contrastGroupName, [factorGroupNames{i} ' MAIN EFFECT']);
        betweenContrastFactors(size(contrastGroupName,1),i)=1;
        basicContrast{i} = [ones(levelGroupNum(i)-1,1) -eye(levelGroupNum(i)-1)];
        if i ~= 1
            groupContrast{i+1}  = ones(1,levelGroupNum(1));
        else
            groupContrast{i+1}  = basicContrast{1};
        end;
        
        for factor = 2:facGroupNum
            if i ~= factor
                groupContrast{i+1} = kron(groupContrast{i+1}, ones(1,levelGroupNum(factor)));
            else
                groupContrast{i+1} = kron(groupContrast{i+1}, basicContrast{factor});
            end;
        end;
    end;
    
    for factor = 1:facGroupNum;
        theNames = [];
        for level = 1:levelGroupNum(factor)
            theNames = [theNames; levelGroupNames{factor}(level)];
        end
        contrastLevelGroupNames{factor+1} = theNames;
        theContrast = groupContrast{factor+1};
        theSum = zeros(1,size(cellMeans,2),size(cellMeans,3));
        theNum = 0;
        for group = 1:grpNum
            if theContrast(1,group) == 1
                theSum = theSum + cellMeans(group,:,:);
                theNum = theNum +1;
            end;
        end;
        theCellGroupMeans = [];
        theCellGroupMeans(1,:,:) = theSum/theNum;

        for forContrast = 1:size(theContrast,1)
            theSum = zeros(1,size(cellMeans,2),size(cellMeans,3));
            theNum = 0;
            for group = 1:grpNum
                if theContrast(forContrast,group) == -1
                    theSum = theSum + cellMeans(group,:,:);
                    theNum = theNum +1;
                end;
            end;
            theCellGroupMeans(forContrast+1,:,:) = theSum/theNum;
        end;
        cellGroupMeans(factor+1,1:size(theCellGroupMeans,1),1:size(theCellGroupMeans,2),:) = theCellGroupMeans;
    end;
    
    if facGroupNum > 1
        for i1 = 1:(facGroupNum-1)
            for i2 = (i1 + 1):facGroupNum
                facCombs = [facCombs; zeros(1,facGroupNum)];
                facCombs(size(facCombs,1),i1) = 1;
                facCombs(size(facCombs,1),i2) = 1;
                contrastGroupName = char(contrastGroupName, [factorGroupNames{i1} ' * ' factorGroupNames{i2} ' INTERACTION EFFECT']);
                betweenContrastFactors(size(contrastGroupName,1),[i1 i2])=1;
                
                contrastNumber = size(groupContrast,2)+1;
                if facCombs(size(facCombs,1),1) == 0
                    groupContrast{contrastNumber}  = ones(1,levelGroupNum(1));
                else
                    groupContrast{contrastNumber}  = basicContrast{1};
                end;
                
                for factor = 2:facGroupNum
                    if facCombs(size(facCombs,1),factor) == 0
                        groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, ones(1,levelGroupNum(factor)));
                    else
                        groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, basicContrast{factor});
                    end;
                end;
                
                theContrast1 = groupContrast{i1+1};
                cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)]; %which groups represent the levels of factor 1
                theContrast2 = groupContrast{i2+1};
                cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)]; %which groups represent the levels of factor 2
                
                theCellMeans = [];
                theNames = [];
                theLevel=0; %level of the between group contrast
                for level1 = 1:levelGroupNum(i1)
                    for level2 = 1:levelGroupNum(i2)
                        theLevel=theLevel+1;
                        theNames = [theNames; [levelGroupNames{i1}(level1) levelGroupNames{i2}(level2)]];
                        theCellKey = cellMeansKey1(level1,:).*cellMeansKey2(level2,:); %which groups are at the intesection of the two factor levels
                        theSum = zeros(1,size(cellMeans,2),size(cellMeans,3));
                        theNum = 0;
                        for group = 1:grpNum %add together the different groups that are part of this between cell
                            if theCellKey(group) == 1
                                theSum = theSum + cellMeans(group,:,:);
                                theNum = theNum +1;
                            end;
                        end;
                        theCellGroupMeans(theLevel,:,:) = theSum/theNum;
                    end;
                end
                contrastLevelGroupNames{contrastNumber} = theNames;
                cellGroupMeans(contrastNumber,1:size(theCellGroupMeans,1),:,:) = theCellGroupMeans;
            end;
        end;
    end;
    
    if facGroupNum > 2
        for i1 = 1:(facGroupNum-2)
            for i2 = (i1 + 1):(facGroupNum-1)
                for i3 = (i2 + 1):facGroupNum
                    facCombs = [facCombs; zeros(1,facGroupNum)];
                    facCombs(size(facCombs,1),i1) = 1;
                    facCombs(size(facCombs,1),i2) = 1;
                    facCombs(size(facCombs,1),i3) = 1;
                    contrastGroupName = char(contrastGroupName, [factorGroupNames{i1} ' * ' factorGroupNames{i2} ' * ' factorGroupNames{i3} ' INTERACTION EFFECT']);
                    betweenContrastFactors(size(contrastGroupName,1),[i1 i2 i3])=1;
                    
                    contrastNumber = size(groupContrast,2)+1;
                    if facCombs(size(facCombs,1),1) == 0
                        groupContrast{contrastNumber}  = ones(1,levelGroupNum(1));
                    else
                        groupContrast{contrastNumber}  = basicContrast{1};
                    end;
                    
                    for factor = 2:facGroupNum
                        if facCombs(size(facCombs,1),factor) == 0
                            groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, ones(1,levelGroupNum(factor)));
                        else
                            groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, basicContrast{factor});
                        end;
                    end;
                    
                    theContrast1 = groupContrast{i1+1};
                    cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                    theContrast2 = groupContrast{i2+1};
                    cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                    theContrast3 = groupContrast{i3+1};
                    cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                    
                    theCellMeans = [];
                    theNames = [];
                    theLevel=0; %level of the between group contrast
                    for level1 = 1:levelGroupNum(i1)
                        for level2 = 1:levelGroupNum(i2)
                            for level3 = 1:levelGroupNum(i3)
                                theLevel=theLevel+1;
                                theNames = [theNames; [levelGroupNames{i1}(level1) levelGroupNames{i2}(level2) levelGroupNames{i3}(level3)]];
                                theCellKey = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:); %which groups are at the intesection of the two factor levels
                                theSum = zeros(1,size(cellMeans,2),size(cellMeans,3));
                                theNum = 0;
                                for group = 1:grpNum %add together the different groups that are part of this between cell
                                    if theCellKey(group) == 1
                                        theSum = theSum + cellMeans(group,:,:);
                                        theNum = theNum +1;
                                    end;
                                end;
                                theCellGroupMeans(theLevel,:,:) = theSum/theNum;
                            end;
                        end
                    end;
                    contrastLevelGroupNames{contrastNumber} = theNames;
                    cellGroupMeans(contrastNumber,1:size(theCellGroupMeans,1),:,:) = theCellGroupMeans;
                end;
            end;
        end;
    end;
    
    
    if facGroupNum > 3
        for i1 = 1:(facGroupNum-3)
            for i2 = (i1 + 1):(facGroupNum-2)
                for i3 = (i2 + 1):(facGroupNum-1)
                    for i4 = (i3 + 1):facGroupNum
                        facCombs = [facCombs; zeros(1,facGroupNum)];
                        facCombs(size(facCombs,1),i1) = 1;
                        facCombs(size(facCombs,1),i2) = 1;
                        facCombs(size(facCombs,1),i3) = 1;
                        facCombs(size(facCombs,1),i4) = 1;
                        contrastGroupName = char(contrastGroupName, [factorGroupNames{i1} ' * ' factorGroupNames{i2} ' * ' factorGroupNames{i3} ' * ' factorGroupNames{i4} ' INTERACTION EFFECT']);
                        betweenContrastFactors(size(contrastGroupName,1),[i1 i2 i3 i4])=1;
                        
                        contrastNumber = size(groupContrast,2)+1;
                        if facCombs(size(facCombs,1),1) == 0
                            groupContrast{contrastNumber}  = ones(1,levelGroupNum(1));
                        else
                            groupContrast{contrastNumber}  = basicContrast{1};
                        end;
                        
                        for factor = 2:facGroupNum
                            if facCombs(size(facCombs,1),factor) == 0
                                groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, ones(1,levelGroupNum(factor)));
                            else
                                groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, basicContrast{factor});
                            end;
                        end;
                        
                        theContrast1 = groupContrast{i1+1};
                        cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                        theContrast2 = groupContrast{i2+1};
                        cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                        theContrast3 = groupContrast{i3+1};
                        cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                        theContrast4 = groupContrast{i4+1};
                        cellMeansKey4 = [(theContrast4(1,:) == 1); (theContrast4 == -1)];
                        
                        theCellMeans = [];
                        theNames = [];
                        theLevel=0; %level of the between group contrast
                        for level1 = 1:levelGroupNum(i1)
                            for level2 = 1:levelGroupNum(i2)
                                for level3 = 1:levelGroupNum(i3)
                                    for level4 = 1:levelGroupNum(i4)
                                        theLevel=theLevel+1;
                                        theNames = [theNames; [levelGroupNames{i1}(level1) levelGroupNames{i2}(level2) levelGroupNames{i3}(level3) levelGroupNames{i4}(level4)]];
                                        theCellKey = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:).*cellMeansKey4(level4,:); %which groups are at the intesection of the two factor levels
                                        theSum = zeros(1,size(cellMeans,2),size(cellMeans,3));
                                        theNum = 0;
                                        for group = 1:grpNum %add together the different groups that are part of this between cell
                                            if theCellKey(group) == 1
                                                theSum = theSum + cellMeans(group,:,:);
                                                theNum = theNum +1;
                                            end;
                                        end;
                                        theCellGroupMeans(theLevel,:,:) = theSum/theNum;
                                    end;
                                end
                            end;
                        end;
                        contrastLevelGroupNames{contrastNumber} = theNames;
                        cellGroupMeans(contrastNumber,1:size(theCellGroupMeans,1),:,:) = theCellGroupMeans;
                    end;
                end;
            end;
        end;
    end;
    
    
    if facGroupNum > 4
        for i1 = 1:(facGroupNum-4)
            for i2 = (i1 + 1):(facGroupNum-3)
                for i3 = (i2 + 1):(facGroupNum-2)
                    for i4 = (i3 + 1):(facGroupNum-1)
                        for i5 = (i4 + 1):facGroupNum
                            facCombs = [facCombs; zeros(1,facGroupNum)];
                            facCombs(size(facCombs,1),i1) = 1;
                            facCombs(size(facCombs,1),i2) = 1;
                            facCombs(size(facCombs,1),i3) = 1;
                            facCombs(size(facCombs,1),i4) = 1;
                            facCombs(size(facCombs,1),i5) = 1;
                            contrastGroupName = char(contrastGroupName, [factorGroupNames{i1} ' * ' factorGroupNames{i2} ' * ' factorGroupNames{i3} ' * ' factorGroupNames{i4} ' * ' factorGroupNames{i5} ' INTERACTION EFFECT']);
                            betweenContrastFactors(size(contrastGroupName,1),[i1 i2 i3 i4 i5])=1;
                            
                            contrastNumber = size(groupContrast,2)+1;
                            if facCombs(size(facCombs,1),1) == 0
                                groupContrast{contrastNumber}  = ones(1,levelGroupNum(1));
                            else
                                groupContrast{contrastNumber}  = basicContrast{1};
                            end;
                            
                            for factor = 2:facGroupNum
                                if facCombs(size(facCombs,1),factor) == 0
                                    groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, ones(1,levelGroupNum(factor)));
                                else
                                    groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, basicContrast{factor});
                                end;
                            end;
                            
                            theContrast1 = groupContrast{i1+1};
                            cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                            theContrast2 = groupContrast{i2+1};
                            cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                            theContrast3 = groupContrast{i3+1};
                            cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                            theContrast4 = groupContrast{i4+1};
                            cellMeansKey4 = [(theContrast4(1,:) == 1); (theContrast4 == -1)];
                            theContrast5 = groupContrast{i5+1};
                            cellMeansKey5 = [(theContrast5(1,:) == 1); (theContrast5 == -1)];
                            
                            
                            theCellMeans = [];
                            theNames = [];
                            theLevel=0; %level of the between group contrast
                            for level1 = 1:levelGroupNum(i1)
                                for level2 = 1:levelGroupNum(i2)
                                    for level3 = 1:levelGroupNum(i3)
                                        for level4 = 1:levelGroupNum(i4)
                                            for level5 = 1:levelGroupNum(i5)
                                                theLevel=theLevel+1;
                                                theNames = [theNames; [levelGroupNames{i1}(level1) levelGroupNames{i2}(level2) levelGroupNames{i3}(level3) levelGroupNames{i4}(level4) levelGroupNames{i5}(level5)]];
                                                theCellKey = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:).*cellMeansKey4(level4,:).*cellMeansKey5(level5,:); %which groups are at the intesection of the two factor levels
                                                theSum = zeros(1,size(cellMeans,2),size(cellMeans,3));
                                                theNum = 0;
                                                for group = 1:grpNum %add together the different groups that are part of this between cell
                                                    if theCellKey(group) == 1
                                                        theSum = theSum + cellMeans(group,:,:);
                                                        theNum = theNum +1;
                                                    end;
                                                end;
                                                theCellGroupMeans(theLevel,:,:) = theSum/theNum;
                                            end;
                                        end
                                    end;
                                end;
                            end;
                            contrastLevelGroupNames{contrastNumber} = theNames;
                            cellGroupMeans(contrastNumber,1:size(theCellGroupMeans,1),:,:) = theCellGroupMeans;
                        end;
                    end;
                end;
            end;
        end;
    end;
    
    
    if facGroupNum > 5
        for i1 = 1:(facGroupNum-5)
            for i2 = (i1 + 1):(facGroupNum-4)
                
                for i3 = (i2 + 1):(facGroupNum-3)
                    for i4 = (i3 + 1):(facGroupNum-2)
                        for i5 = (i4 + 1):(facGroupNum-1)
                            for i6 = (i5 + 1):facGroupNum
                                facCombs = [facCombs; zeros(1,facGroupNum)];
                                facCombs(size(facCombs,1),i1) = 1;
                                facCombs(size(facCombs,1),i2) = 1;
                                facCombs(size(facCombs,1),i3) = 1;
                                facCombs(size(facCombs,1),i4) = 1;
                                facCombs(size(facCombs,1),i5) = 1;
                                facCombs(size(facCombs,1),i6) = 1;
                                contrastGroupName = char(contrastGroupName, [factorGroupNames{i1} ' * ' factorGroupNames{i2} ' * ' factorGroupNames{i3} ' * ' factorGroupNames{i4} ' * ' factorGroupNames{i5} ' * ' factorGroupNames{i6} ' INTERACTION EFFECT']);
                                betweenContrastFactors(size(contrastGroupName,1),[i1 i2 i3 i4 i5 i6])=1;
                                
                                contrastNumber = size(groupContrast,2)+1;
                                if facCombs(size(facCombs,1),1) == 0
                                    groupContrast{contrastNumber}  = ones(1,levelGroupNum(1));
                                else
                                    groupContrast{contrastNumber}  = basicContrast{1};
                                end;
                                
                                for factor = 2:facGroupNum
                                    if facCombs(size(facCombs,1),factor) == 0
                                        groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, ones(1,levelGroupNum(factor)));
                                    else
                                        groupContrast{contrastNumber} = kron(groupContrast{contrastNumber}, basicContrast{factor});
                                    end;
                                end;
                                
                                theContrast1 = groupContrast{i1+1};
                                cellMeansKey1 = [(theContrast1(1,:) == 1); (theContrast1 == -1)];
                                theContrast2 = groupContrast{i2+1};
                                cellMeansKey2 = [(theContrast2(1,:) == 1); (theContrast2 == -1)];
                                theContrast3 = groupContrast{i3+1};
                                cellMeansKey3 = [(theContrast3(1,:) == 1); (theContrast3 == -1)];
                                theContrast4 = groupContrast{i4+1};
                                cellMeansKey4 = [(theContrast4(1,:) == 1); (theContrast4 == -1)];
                                theContrast5 = groupContrast{i5+1};
                                cellMeansKey5 = [(theContrast5(1,:) == 1); (theContrast5 == -1)];
                                theContrast6 = groupContrast{i6+1};
                                cellMeansKey6 = [(theContrast6(1,:) == 1); (theContrast6 == -1)];
                                
                                theCellMeans = [];
                                theNames = [];
                                theLevel=0; %level of the between group contrast
                                for level1 = 1:levelGroupNum(i1)
                                    for level2 = 1:levelGroupNum(i2)
                                        for level3 = 1:levelGroupNum(i3)
                                            for level4 = 1:levelGroupNum(i4)
                                                for level5 = 1:levelGroupNum(i5)
                                                    for level6 = 1:levelGroupNum(i6)
                                                        theLevel=theLevel+1;
                                                        theNames = [theNames; [levelGroupNames{i1}(level1) levelGroupNames{i2}(level2) levelGroupNames{i3}(level3) levelGroupNames{i4}(level4) levelGroupNames{i5}(level5) levelGroupNames{i6}(level6)]];
                                                        theCellKey = cellMeansKey1(level1,:).*cellMeansKey2(level2,:).*cellMeansKey3(level3,:).*cellMeansKey4(level4,:).*cellMeansKey5(level5,:).*cellMeansKey6(level6,:); %which groups are at the intesection of the two factor levels
                                                        theSum = zeros(1,size(cellMeans,2),size(cellMeans,3));
                                                        theNum = 0;
                                                        for group = 1:grpNum %add together the different groups that are part of this between cell
                                                            if theCellKey(group) == 1
                                                                theSum = theSum + cellMeans(group,:,:);
                                                                theNum = theNum +1;
                                                            end;
                                                        end;
                                                        theCellGroupMeans(theLevel,:,:) = theSum/theNum;
                                                    end;
                                                end
                                            end;
                                        end;
                                    end;
                                end;
                                contrastLevelGroupNames{contrastNumber} = theNames;
                                cellGroupMeans(contrastNumber,1:size(theCellGroupMeans,1),:,:) = theCellGroupMeans;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    
end;

if ~followup
    fprintf(outfid,'WELCH-JAMES APPROXIMATE DF SOLUTION</BR>');
    if OPT1==0
        fprintf(outfid,'LEAST SQUARES MEANS & VARIANCES</BR>');
    end;
    if OPT1==1
        fprintf(outfid,'TRIMMED MEANS & WINSORIZED VARIANCES</BR>');
    end;
    if OPT1==1
        fprintf(outfid,['PERCENTAGE OF TRIMMING: ' sprintf('%4.2f', PER) '</BR>']);
    end;
    if OPT2==0
        fprintf(outfid,'F DISTRIBUTION CRITICAL VALUE</BR>');
    end;
    if OPT2==1
        fprintf(outfid,'BOOTSTRAP CRITICAL VALUE FOR SINGLE TEST STATISTIC</BR>');
    end;
    if OPT2==1
        fprintf(outfid,['NUMBER OF BOOTSTRAP SAMPLES: ' sprintf('%4.0f',NUMSIM) '</BR>']);
        fprintf(outfid,['MEDIAN OF BOOTSTRAP RUNS: ' sprintf('%4.0f',REPS) '</BR>']);
        fprintf(outfid,['STARTING SEED (multiplied by the run number): ' sprintf('%15.0f',SEED) '</BR>']);
    end;
    
    fprintf(outfid,['Number of subjects in each group: ' num2str(NX) '</BR>']);
    if OPT1==1
        fprintf(outfid,['Number trimmed from each end of the groups: ' num2str(floor(PER*NX)) '</BR>']);
    end;
    if sum(NX) ~= sum(newNX)
        fprintf(outfid,['Number of subjects in each group after missing data dropped: ' num2str(newNX) '</BR>']);
    end;
    
    fprintf(outfid,['Uncorrected alpha criteria: ' num2str(alpha.uncorrected) '</BR>']);
    fprintf(outfid,['Corrected alpha criteria: ' num2str(alpha.corrected) '</BR>']);
end;

%generate corrected dataset for vector scaling test in case it is needed

%make the assumption that the electrode factors are the lowest factors
totElecLevels=1;
for i=1:length(levelNum)
    if elecFactors(i)
        totElecLevels=totElecLevels*levelNum(i);
    end;
end;

%scaling is done within each cell (defined by within factors other than electrode factors and between factors) with only electrode levels within each such
%cell.  The scaling is based on the mean value for each cell per McCarthy & Wood (1985) to get "a more stable
%estimate."  The mean values are the trimmed means in order to keep the scaling procedure consistent with the
%robust statistic.  The statement made in Dien & Santuzzi (2005) that Haig, Gordon, & Hook (1997) argued for calculating
%the vector scores based on the individual subjects was an error.

Yscaled=zeros(size(Y));
for col=1:totElecLevels:varNum
    for grp=1:grpNum
        Yscaled(1+sum(NX(1:grp-1)):sum(NX(1:grp)),col:col+totElecLevels-1)=Y(1+sum(NX(1:grp-1)):sum(NX(1:grp)),col:col+totElecLevels-1)/sqrt(sum(MUHAT((grp-1)*varNum+col:(grp-1)*varNum+col+totElecLevels-1).^2));
    end;
end;

%if this is a follow-up analysis rather than the primary analysis
if ~isempty(followupBetween) || ~isempty(followupWithin)
    theContrast=zeros(length(followupBetween),size(contrastGroupName,1));
    for theFactor = 1:length(followupBetween)
        theMatches=regexp(cellstr(contrastGroupName),['^' followupBetween{theFactor} '\s\w*']);
        theMatches2=regexp(cellstr(contrastGroupName),['\w*\s' followupBetween{theFactor} '\s\w*']);
        for i=1:size(contrastGroupName,1)
            theContrast(theFactor,i)= ~isempty(theMatches{i}) || ~isempty(theMatches2{i});
        end;
    end;
    betweenContrastList=min(find(sum(theContrast,1) == length(followupBetween))); %the lowest level contrast that contains all the factor names desired
else
    betweenContrastList=1:size(groupContrast,2);
end;

%scale contrast weights so that the positive and negative terms sum to +1
%and -1 respectively since this affects the ES and the CI statistics
%although not the overall test statistic.
%assume sum of positive and negative weights equal each other
%except for case where all weights are positive.

for iContrast=1:length(groupContrast)
    posWeights=find(sign(groupContrast{iContrast})==1);
    scaleWeights=sum(groupContrast{iContrast}(posWeights));
    groupContrast{iContrast}=groupContrast{iContrast}/scaleWeights;
end;
for iContrast=1:length(contrast)
    posWeights=find(sign(contrast{iContrast})==1);
    scaleWeights=sum(contrast{iContrast}(posWeights));
    contrast{iContrast}=contrast{iContrast}/scaleWeights;
end;

for betweenGroup = betweenContrastList
    fprintf(outfid,[repmat('&nbsp;',1,followup*2) '</BR>' repmat('#',1,60) '</BR>']);
    fprintf(outfid,[repmat('&nbsp;',1,followup*2) contrastGroupName(betweenGroup,:) '</BR>']);
    
    if ~isempty(followupBetween) || ~isempty(followupWithin)
        theContrast=zeros(length(followupWithin),size(contrastName,1));
        for theFactor = 1:length(followupWithin)
            theMatches=regexp(cellstr(contrastName),['^' followupWithin{theFactor} '\s\w*']);
            theMatches2=regexp(cellstr(contrastName),['\w*\s' followupWithin{theFactor} '\s\w*']);
            for i=1:size(contrastName,1)
                theContrast(theFactor,i)= ~isempty(theMatches{i}) || ~isempty(theMatches2{i});
            end;
        end;
        withinContrastList=min(find(sum(theContrast,1) == length(followupWithin))); %the lowest level contrast that contains all the factor names desired
    else
        if betweenGroup == 1
            firstWithin = 2; %analysis dropping both between and within factors would be meaningless
        else
            firstWithin = 1;
        end;
        withinContrastList=firstWithin:size(contrast,2);
    end;
    
    for withinGroup = withinContrastList
        if (length(contrast{withinGroup}) > 1) || ((betweenGroup > 1) && (withinGroup == 1))
            
            fprintf(outfid,[repmat('&nbsp;',1,followup*2) '</BR>' repmat('-',1,60) '</BR>']);
            fprintf(outfid,[repmat('&nbsp;',1,followup*2) contrastName(withinGroup,:) '</BR>']);
            resultsTot=cell(REPS,1);
            pList=zeros(REPS,1);
            for iRep=1:REPS
                [MUHAT, SIGMA, RESULTS]=ep_WJGLMml(Y, NX, groupContrast{betweenGroup}, contrast{withinGroup}', OPT1, PER, OPT2, NUMSIM, SEED*iRep, MISSING, OPT3, alpha, SCALE, LOC1, LOC2);
                resultsTot{iRep}=RESULTS;
                pList(iRep)=RESULTS(4);
            end;
            [A B]=sort(pList);
            RESULTS=resultsTot{find(B==ceil(length(pList)/2))};
            if isempty(RESULTS)
                return
            end;
            alpha.pvar=std(pList)*2;
            ep_printADFResults(RESULTS, alpha, followup, outfid, NUMSIM);
            
            %print out the cell means
            fprintf(outfid,[repmat('&nbsp;',1,followup*2) 'Averaged Trimmed Cell Means: </BR>']);
            theCellNames = contrastLevelNames{withinGroup};
            line = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
            for theCell = 1:size(theCellNames,1)
                theName = theCellNames(theCell,:);
                line = [line theName repmat('&nbsp;',1,maxDigits-length(theName))  '&nbsp;&nbsp;'];
            end;
            fprintf(outfid,[repmat('&nbsp;',1,followup*2) line '</BR>']); %the line with the level labels
            theWithinLevels=length(find(cellGroupMeans(betweenGroup,1,withinGroup,:)));
            theBetweenLevels=length(find(cellGroupMeans(betweenGroup,:,withinGroup,1)));
            
            theBetweenNames = contrastLevelGroupNames{betweenGroup};
            
            for forBetween = 1:theBetweenLevels
                theData = [theBetweenNames(forBetween,:) repmat('&nbsp;',1,5-length(theBetweenNames(forBetween,:)))]; %the between group level name
                theData = strrep(theData,' ','&nbsp;');
                for forWithin = 1:theWithinLevels
                    %cellGroupMeans is: (the between group contrast, the means for each cell in this contrast, the within
                    %group contrast,the means for each cell in this contrast)
                    eval(['theNumber = sprintf(''%+' num2str(maxDigits) '.2f'', cellGroupMeans(betweenGroup,forBetween,withinGroup,forWithin));']);
                    theData = [theData theNumber];
                    theData = [theData repmat('&nbsp;',1,maxDigits-length(theNumber)+2)];
                end;
                fprintf(outfid,[repmat('&nbsp;',1,followup*2) theData '</BR>']); %the cell means
            end;
            
            %perform vector scaling test.  The resulting cell means are not reported because, as demonstrated by Haig,
            %Gordon, & Hook (1997), they are not meaningful in of themselves.
            if RESULTS(4) <= alpha.uncorrected && ~isempty(elecFactors)
                if any(and(contrastFactors(withinGroup,:),elecFactors')) && any((contrastFactors(withinGroup,:)-elecFactors')>0) %both electrode & non-electrode factors?
                    resultsTot=cell(REPS,1);
                    pList=zeros(REPS,1);
                    for iRep=1:REPS
                        [MUHAT, SIGMA, RESULTSMW]=ep_WJGLMml(Yscaled, NX, groupContrast{betweenGroup}, contrast{withinGroup}', OPT1, PER, OPT2, NUMSIM, SEED*iRep, MISSING, OPT3, alpha, SCALE, LOC1, LOC2);
                        resultsTot{iRep}=RESULTSMW;
                        pList(iRep)=RESULTSMW(4);
                    end;
                    [A B]=sort(pList);
                    RESULTSMW=resultsTot{find(B==ceil(length(pList)/2))};
                    if isempty(RESULTSMW)
                        return
                    end;
                    alpha.pvar=std(pList)*2;
                    fprintf(outfid,[repmat('&nbsp;',1,followup*2) '</BR><i>Vector Scaling Test (is the scalp topography effect genuine?):</BR>']);
                    ep_printADFResults(RESULTSMW, alpha, followup, outfid, NUMSIM);
                    fprintf(outfid,'</i>');
                end;
            end;
            
            %perform follow-up ANOVAs of interactions by recursively calling upon this function
            %For the follow-up ANOVAs, there are three types of ANOVA factors to consider:
            % 1) the factor where one level is being held constant and the other levels are excluded.
            % 2) the factors which are entered in the normal way (retained factors)
            % 3) the factors which are not being considered and whose levels are therefore combined together.
            % The factor structure for the ANOVA is based on the full model but only the contrast corresponding to the
            % follow-up is actually computed.
            
            if RESULTS(4) <= alpha.uncorrected && ((sum(contrastFactors(withinGroup,:)) + sum(betweenContrastFactors(betweenGroup,:)))>1)
                betweenFactorList=find(betweenContrastFactors(betweenGroup,:));
                withinFactorList=find(contrastFactors(withinGroup,:));
                for betweenFact=1:length(betweenFactorList)
                    theBetweenFactor=betweenFactorList(betweenFact); %the one being kept constant
                    retainedFactors=setdiff(betweenFactorList,betweenFactorList(betweenFact)); %the factors being analyzed in the follow-up
                    otherFacs=(setdiff([1:length(factorGroupNames)],theBetweenFactor)); %the set of factors other than the one being held constant
                    for betweenLevel=1:length(levelGroupNames{theBetweenFactor})
                        theConstantGrps=find(strcmp(levelGroupNames{theBetweenFactor}(betweenLevel),factorGroupTable(:,theBetweenFactor))); %the between cells that are entering the follow-up ANOVAs
                        constantFactorGroupTable=factorGroupTable(theConstantGrps,:); %the part of the factor table that is being held constant and thus entered into the follow-up
                        constantNX=NX(theConstantGrps); %the part of NX that is being held constant and thus entered into the follow-up
                        constantY=[];
                        
                        %retain the appropriate between groups and collapse over the irrelevant between group factors
                        counter=1;
                        for i=1:length(NX)
                            if any(theConstantGrps == i)
                                constantY=[constantY; Y(counter:NX(i)+counter-1,:)]; %the part of Y that is being held constant and thus entered into the follow-up
                            end;
                            counter=counter+NX(i);
                        end;
                        
                        %do the follow-up ANOVA
                        fprintf(outfid,'<font size="1">');
                        fprintf(outfid,[repmat('&nbsp;',1,(followup+1)*2) 'Holding level ' levelGroupNames{theBetweenFactor}(betweenLevel) ' of factor ' factorGroupNames{theBetweenFactor} ' constant.']);
                        ep_ADF(constantY, constantNX, OPT1, PER, OPT2, NUMSIM, REPS, SEED, MISSING, factorNames, levelNames, factorGroupNames(otherFacs), levelGroupNames(otherFacs), elecFactors, alpha, followup+1, OPT3, SCALE, LOC1, LOC2, outfid, factorNames(withinFactorList), factorGroupNames(retainedFactors));
                        fprintf(outfid,'<font size="2">');
                    end;
                end;
                for withinFact=1:length(withinFactorList)
                    theWithinFactor=withinFactorList(withinFact);
                    retainedFactors=setdiff(withinFactorList,withinFactorList(withinFact)); %the factors being analyzed in the follow-up
                    otherFacs=setdiff([1:length(factorNames)],theWithinFactor); %the set of factors other than the one being held constant
                    for withinLevel=1:length(levelNames{theWithinFactor})
                        theConstantLvls=find(strcmp(levelNames{theWithinFactor}(withinLevel),factorTable(theWithinFactor,:))); %the within cells that are entering the follow-up ANOVAs
                        constantFactorTable=factorTable(:,theConstantLvls); %the part of the factor table that is being held constant and thus entered into the follow-up
                        constantY=[];
                        
                        %retain the appropriate within groups and collapse over the irrelevant within group factors
                        counter=1;
                        for i=1:varNum
                            if any(theConstantLvls == i)
                                constantY(:,end+1)=Y(:,i); %the part of Y that is being held constant and thus entered into the follow-up
                            end;
                        end;
                        
                        fprintf(outfid,'<font size="1">');
                        fprintf(outfid,[repmat('&nbsp;',1,(followup+1)*2) 'Holding level ' levelNames{theWithinFactor}(withinLevel) ' of factor ' factorNames{theWithinFactor} ' constant.']);
                        ep_ADF(constantY, NX, OPT1, PER, OPT2, NUMSIM, REPS, SEED, MISSING, factorNames(otherFacs), levelNames(otherFacs), factorGroupNames, levelGroupNames, elecFactors(otherFacs), alpha, followup+1, OPT3, SCALE, LOC1, LOC2, outfid, factorNames(retainedFactors), factorGroupNames(betweenFactorList));
                        fprintf(outfid,'<font size="2">');
                    end;
                end;
            end;
        end;
    end;
end;

fprintf(outfid,[repmat('&nbsp;',1,followup*2) '</BR>' repmat('#',1,60) '</BR>']);
if internalFID
    fclose(outfid);
end;