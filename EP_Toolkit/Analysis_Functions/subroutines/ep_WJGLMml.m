function [MUHAT, SIGMA, RESULTS]=ep_WJGLMml(Y, NX, C, U, OPT1, PER, OPT2, NUMSIM, SEED, MISSING, OPT3, ALPHA, SCALE, LOC1, LOC2);
%function [MUHAT, SIGMA, RESULTS]=ep_WJGLMml(Y, NX, C, U, OPT1, PER, OPT2, NUMSIM, SEED, MISSING, OPT3, ALPHA, SCALE, LOC1, LOC2);
%Function for calculating robust statistics using Welch-James ADF, boostrapping, and trimmed
%means plus winsorized variances and covariances.
%
%Based on SAS/IML code made available by Lisa Lix at:
%http://www.usaskhealthdatalab.ca/sas-programs/
%
%Converted to Matlab 9/14/15 by Joseph Dien (jdien07@mac.com)
%
%Wilcox, R. R. (2001). Fundamentals of modern statistical methods. New York: Springer-Verlag.
%
%Keselman, H. J., Wilcox, R. R., & Lix, L. M. (2003). A generally robust approach
%to hypothesis testing in independent and correlated groups designs
%Psychophysiology, 40, 586-596.
%
%Keselman, H. J., Algina, J., Lix, L. M., Wilcox, R. R., & Deering, K. N. (2008). 
%A generally robust approach for testing hypotheses and setting confidence intervals for effect sizes. 
%Psychological Methods, 13(2), 110.

%Inputs
%  Y	    : Input data (rows=subjects, columns = cells).  The first NX rows will be assigned to the first group and so forth.
%  NX   	: Number of subjects in each group.  If empty set then will assume a single group.
%  C    	: Contrast row vector for between factors.  Set to 1 in the case where there is only one group.
%  U    	: Contrast column vector for within factors.  Numbers should sum to zero.  Set to empty set [] for analyses with no within-factors.
%  OPT1     : Activate trimming option (rounded down).  0 = no and 1 = yes.
%  PER  	: Percentage to trim the means.  .05 is the number recommended for ERP data by Dien.
%  OPT2 	: Activate Welch-James statistic.  0 = no and 1 = yes.
%  NUMSIM	: Number of simulations used to generate bootstrapping statistic.  p-values will be unstable if too low.  50000 informally recommended.
%  SEED 	: Seed for random number generation.  0 specifies random SEED. 1000 arbitrarily suggested as SEED to ensure RESULTS are replicable.
%  MISSING  : Number to be treated as a missing value.  Observations with missing values are dropped from the analysis.
%  OPT3     : Provide effect sizes. 0 = no and 1 = yes.
%  ALPHA    : Alpha significance level.
%	.corrected : Alpha level corrected for multiple comparisons (not used)
%	.uncorrected : Alpha level not corrected for multiple comparisons
%  SCALE    : "SCALE is a scalar indicator to control the use of a scaling factor for the effect size estimator 
%             (i.e., .642 for 20% symmetric trimming) when robust estimators are adopted. 
%             It takes a value of 0 or 1; a zero indicates that no scaling factor will be used, 
%             while a 1 indicates that a scaling factor will be adopted. The default is SCALE=1.
%  LOC1&LOC2: "If the user specifies LOC1 = 0 and LOC2 = 0, the square root of the average of the variances over the cells 
%             involved in the contrast is used as the standardizer. If the
%             user specifies LOC1 = 99 and LOC2 = 99, no standardizer is selected."

%Outputs
%  MUHAT	: Vector of trimmed means used in calculations. 
%             For between group analyses, within group conditions vary fastest and between group conditions vary slowest.
%  SIGMA	: Winsorized covariance matrix
%  RESULTS	: (1) = test statistic
%  RESULTS	: (2) = Numerator DF
%  RESULTS	: (3) = Denominator DF
%  RESULTS	: (4) = Significance
%  RESULTS	: (5) = Effect size using delta-hat-star standardizer (only for contrasts with 1 df).
%  RESULTS	: (6) = Lower confidence limit
%  RESULTS	: (7) = Upper confidence limit
%  RESULTS	: (8) = Scaling factor for effect size estimator if SCALE option is chosen.
%  RESULTS	: (9) = Number of singular matrices during calculation of statistics.
%  RESULTS	: (10) = Number of nearly singular matrices during calculation of statistics.

% modified 10/7/07 JD
% Missing number feature added.
%
% modified 2/11/09 JD
% Warning messages only output once.
%
% modified 6/20/11 JD
% Handles new warning code for ill-conditioned matrices.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% modified 9/14/15 JD
% Revised based on updated SAS/IML code by Lisa Lix to provide effect sizes.
%
% bugfix 10/27/15 JD
% Restored warning messages that were erroneously deleted for problematic analyses and sets p-value output to NaN (not a number).
% Handles case where the entire bootstrap sample is drawn from the same observation and treats as an F of inf.
%
% bugfix 11/30/15 JD
% No longer produces effect sizes for contrasts of more than 1 df.
% Handles case where almost the entire bootstrap sample is drawn from the same observation and was generating an NaN value. 
%
% bugfix 12/6/15 JD
% Fixed NaN effect sizes when PER is set to zero.
%
% modified 12/10/15 JD
% Modified so that singular and near-singular runs are dropped from the
% bootstrapping run and a warning is issued.
% If effect size d* is made positive, confidence interval signs also flipped.
% Added information on number of singular and nearly singular matrices during calculation of statistics to RESULTS.
%
% modified 12/27/15 JD
% Per new information from Lisa Lix, effect size output disabled for within-group
% contrasts.
%
% bugfix 6/16/17 JD
% Fixed apparently providing effect size for within contasts of more than two levels when should have been disabled.
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

%****DEFINE MODULES TO PERFORM ALL CALCULATIONS****;
%****COMPUTE WELCH-JAMES STATISTIC****;

errorflag.TestSingularMatrix=0;
errorflag.TestIllConditionedMatrix=0;
errorflag.TooFewSubjects=0;

if SEED == 0
    SEED=sum(100*clock);
end;

sprev = rng(SEED,'twister');

numsim_b=NUMSIM;
numsim_es=NUMSIM;
numsim_bc=NUMSIM;

alphaThresh=ALPHA.uncorrected;

if isempty(NX)
    NX = size(Y,1);
end;
if isempty(C)
    C = [1];
end;

if size(C,2) ~= length(NX)
    msg{1}=['Number of between group cells (' num2str(length(NX)) ') must equal number of terms in contrast C (' num2str(size(C,2)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end;

if size(U,1) ~= size(Y,2)
    msg{1}=['Number of within group cells (' num2str(size(Y,2)) ') must equal number of terms in contrast U (' num2str(size(U,1)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end;

if OPT1 ~=0 && OPT1 ~= 1
    msg{1}='OPT1 must equal zero or one.';
    [msg]=ep_errorMsg(msg);
    return
end;

if OPT3 ~=0 && OPT3 ~= 1
    msg{1}='OPT3 must equal zero or one.';
    [msg]=ep_errorMsg(msg);
    return
end;

if (size(U,2) > 1) || (size(C,1) > 1)
    OPT3=0; %cannot calculate effect sizes for more than 1 degree of freedom contrasts
end;

if ~isempty(MISSING) %drop observations with missing data points
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
    NX=newNX;
    Y=Y(good,:);
end;

%**define module to check initial specifications****;
if isempty(U)
    U=eye(size(Y,2));
end;
if size(U,2)>size(U,1)
    msg{1}='Possible Error: Number Of Columns Of U Exceeds Number Of Rows';
    [msg]=ep_errorMsg(msg);
    return
end;

if OPT1==1
    if isempty(PER)
        PER=.20;
    end;
    if PER > .49
        msg{1}='Error: Percentage Of Trimming Exceeds Upper Limit';
        [msg]=ep_errorMsg(msg);
        return
    end;
end;
    
if OPT2==1
    if isempty(numsim_b)
        numsim_b=999;
    end;
end;

if OPT3==1
    if isempty(numsim_es)
        numsim_es=999;
    end;
    if isempty(LOC1)
        LOC1=1;
    end;
    if isempty(LOC2)
        LOC2=1;
    end;
end;
        
if isempty(alphaThresh)
    alphaThresh=.05;
end;
if isempty(SCALE)
    SCALE=1;
end;
if isempty(SEED)
    SEED=0;
end;
if isempty(numsim_bc)
    numsim_bc=699;
end;

for I=1:size(NX,2);
    X1=ones(NX(I),1)*I;
    if I==1
        X=X1;
    else
        X=[X; X1];
    end;
end;
X=design(X); %full design matrix
NTOT=size(Y,1); %number of subjects
WOBS=size(Y,2); %total number of within cells
BOBS=size(X,2); %number of between groups
WOBS1=WOBS-1; %one less than total number of within cells
R=kron(C,U'); %contrast vector combining within and between contrasts

%****compute Welch-James statistic****;
[MUHAT, BHAT, BHATW, YT, DF, errorflag] = mnmod(Y, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
if errorflag.TooFewSubjects
    return;
end;
[SIGMA,STDIZER] = sigmod(YT,X,BHATW,DF,BOBS,WOBS,WOBS1,NX);
[FSTAT,DF1,DF2, errorflag] = testmod(SIGMA,MUHAT,R,DF, BOBS, WOBS, WOBS1, errorflag);

if OPT3==1
    [MULTP,EFFSZ]=wjeffsz(BOBS,WOBS,WOBS1,NTOT,NX,LOC1,LOC2,SCALE,OPT1,PER,R,MUHAT,STDIZER);
end;
if OPT2==1
    for simloop=1:numsim_b;
        [YB1] = bootdat(Y, BHAT, BOBS, NX);
        if ~any(std(YB1) > 10^(-12)) %Matlab is showing errors at about the 10^(-14) level
            FSTATB=inf; %all of the bootstrap observations were identical
        else
            [YB]=bootcen(YB1,BHAT,BOBS,NX);
            [FSTATB, errorflag] = bootstat(YB, OPT1, R, BOBS, WOBS, WOBS1, NTOT, NX, PER, X, SEED, errorflag);
        end;
        if simloop==1
            FMAT=FSTATB;
        else
            FMAT=[FMAT; FSTATB];
        end;
    end;
    FMAT=sort(FMAT);
end;

%***calculate significance level for welch-james statistic****;
RESULTS=zeros(10,1);
RESULTS(1)=FSTAT;
RESULTS(2)=DF1;
RESULTS(3)=DF2;
if OPT2==0
    RESULTS(4)=1-probf(RESULTS(1),DF1,DF2);
end;
RESULTS(5) = NaN;
RESULTS(6) = NaN;
RESULTS(7) = NaN;
RESULTS(8) = NaN;RESULTS(9) = errorflag.TestSingularMatrix;
RESULTS(10) = errorflag.TestIllConditionedMatrix;

% if errorflag.TestSingularMatrix>0
%     msg=['Warning: The matrix was singular to working precision ' sprintf('%4.2f%%',100*(errorflag.TestSingularMatrix/numsim_b)) ' of the time when computing the test statistic.  These runs were ignored and may thus the estimate may have some inaccuracy.'];
%     disp(msg);
% end;
% if errorflag.TestIllConditionedMatrix>0
%     msg=['Warning: Matrix was singular, close to singular or badly scaled ' sprintf('%4.2f%%',100*(errorflag.TestIllConditionedMatrix/numsim_b)) ' of the time when computing the test statistic. These runs were ignored and may thus the estimate may have some inaccuracy.'];
%     disp(msg);
% end;

if OPT3==1
        
    errorflag.TestSingularMatrix=0;
    errorflag.TestIllConditionedMatrix=0;

    for simloop=1:numsim_es;
        [YB1] = bootdat(Y, BHAT, BOBS, NX);
        [MULTP,EFFSZB,errorflag]=bootes(YB1,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag);
        if simloop==1
            esmat=EFFSZB;
        else esmat=[esmat; EFFSZB];
        end;
    end;
    esmat=sort(esmat);
end;
    
if OPT2==1
    avec=(FSTAT<=FMAT);
    numNaN=sum(isnan(FMAT));
    pval=sum(avec)/(numsim_b-numNaN);
    RESULTS(4)=pval;
end;

if OPT3==1
    if (DF1>1) || (WOBS>1)
        %Effect sizes only available for one DF between-group contrasts, pending further research by Dr. Lix.
        RESULTS(5) = NaN;
        RESULTS(6) = NaN;
        RESULTS(7) = NaN;
        RESULTS(8) = NaN;
    else
        numNaN=sum(isnan(esmat));
        numsim_esGood=numsim_es-numNaN;
        if numsim_esGood>0
            ind1=intIML(numsim_esGood*(alphaThresh/2))+1;
            ind2=numsim_esGood-intIML((numsim_esGood*(alphaThresh/2)));
            lcl=esmat(ind1);
            ucl=esmat(ind2);
            [MULTP,EFFSZB,errorflag]=bootes(Y,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag);
            RESULTS(5) = abs(EFFSZB); %convention for Cohen's d is to provide absolute value JD
            RESULTS(6) = lcl;
            RESULTS(7) = ucl;
            if sign(EFFSZB) == -1
                RESULTS(6)=-RESULTS(6);
                RESULTS(7)=-RESULTS(7);
            end;
            RESULTS(8) = MULTP;
            
%             if errorflag.TestSingularMatrix>0
%                 msg=['Warning: The matrix was singular to working precision ' sprintf('%4.2f%%',100*(errorflag.TestSingularMatrix/numsim_es)) ' of the time when computing the effect size.  These runs were ignored and may thus the estimate may have some inaccuracy.'];
%                 disp(msg);
%             end;
%             if errorflag.TestIllConditionedMatrix>0
%                 msg=['Warning: Matrix was singular, close to singular or badly scaled ' sprintf('%4.2f%%',100*(errorflag.TestIllConditionedMatrix/numsim_es)) ' of the time when computing the effect size.  These runs were ignored and may thus the estimate may have some inaccuracy.'];
%                 disp(msg);
%             end;
        end;
    end;
end;

% eval(['save test' sprintf('%03d',i) '.mat']);
% end;

function a = design(X);
a = zeros(size(X,1),length(unique(X)));
for i=1:length(X);
    a(i,X(i)) = 1;
end;

function a = intIML(X);
%perform equivalent of SAS/IML int function
a=X-rem(X,1);

%****define module to compute least squares or trimmed means****;
function [MUHAT, BHAT, BHATW, YT, DF, errorflag] = mnmod(Y, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
    
if OPT1==0
    BHAT=inv(X'*X)*X'*Y;
    BHATW=BHAT;
    YT=Y;
    DF=NX-1;
end;

if OPT1==1
    BHAT=zeros(BOBS,WOBS);
    BHATW=BHAT;
    YT=zeros(NTOT,WOBS);
    DF=zeros(1,BOBS);
    F=1;
    M=0;
    for J=1:size(NX,2); %loop through between groups
        SAMP=NX(J);%size of between group
        L=M+SAMP;
        G=intIML(PER.*SAMP);
        DF(J)=SAMP-2*G-1;
        for K=1:size(Y,2) %loop through within groups
            tempVar=Y(F:L,K);
            NV=tempVar;
            [TEMP2 index] = sort(NV);
            TRIMY=TEMP2(G+1:SAMP-G,:); %trimmed cell
            TRIMMN=sum(sum(TRIMY))/(DF(J)+1); %trimmed mean
            BHAT(J,K)=TRIMMN; %matrix of trimmed means
            MINT=min(min(TRIMY));
            MAXT=max(max(TRIMY));
            for P=1:size(NV,1)
                if NV(P)<=MINT
                    NV(P)=MINT;
                end;
                if NV(P)>=MAXT
                    NV(P)=MAXT;
                end;
            end;
            YT(F:L,K)=NV; %winsorized sample
            WINMN=sum(sum(NV))/SAMP;
            BHATW(J,K)=WINMN;
        end;
        M=L;
        F=F+NX(J);
    end;
end;

if DF == 0
    errorflag.TooFewSubjects=1;
    MUHAT=0;
    BHAT=0;
    BHATW=0;
    YT=0;
    msg{1}='Error, too few subjects.  Degrees of freedom is zero.';
    [msg]=ep_errorMsg(msg);
    return
end;

MUHAT=reshape(BHAT',[],BOBS*WOBS)';

%***DEFINE MODULE TO COMPUTE SIGMA MATRIX****;
function [SIGMA,STDIZER] = sigmod(YT,X,BHATW,DF,BOBS,WOBS,WOBS1,NX);
SIGMA=zeros(WOBS.*BOBS,WOBS.*BOBS);
STDIZER=SIGMA;
for I=1:BOBS
    SIGB=(diag(X(:,I))*YT-X(:,I)*BHATW(I,:))'*(diag(X(:,I))*YT-X(:,I)*BHATW(I,:))/((DF(I)+1)*DF(I));
    F=I*WOBS-WOBS1;
    L=I*WOBS;
    SIGMA(F:L,F:L)=SIGB;
    STDIZER(F:L,F:L)=SIGB*((DF(I)+1)*DF(I))/(NX(I)-1);
end;

%****DEFINE MODULE TO COMPUTE TEST STATISTIC****;
function [FSTAT, DF1, DF2, errorflag] = testmod(SIGMA, MUHAT, R, DF, BOBS, WOBS, WOBS1, errorflag);
warning off MATLAB:singularMatrix %turn off display of warning otherwise one gets dozens of messages
warning off MATLAB:illConditionedMatrix
warning off MATLAB:nearlySingularMatrix
lastwarn('');
T=(R*MUHAT)'*inv(R*SIGMA*R')*(R*MUHAT); %T stat squared and without df correction
A=0;
IMAT=eye(WOBS);
for I=1:BOBS;
    QMAT=zeros(BOBS*WOBS,BOBS*WOBS);
    F=I*WOBS-WOBS1;
    L=I*WOBS;
    QMAT(F:L,F:L)=IMAT;
    PROD=(SIGMA*R')*inv(R*SIGMA*R')*R*QMAT;
    A=A+(trace(PROD*PROD)+trace(PROD)^2)/DF(I);
end;
A=A/2; %for adjusted degrees of freedom
DF1=size(R,1);
DF2=DF1*(DF1+2)/(3*A);
CVAL=DF1+2*A-6*A/(DF1+2);
FSTAT=T/CVAL;
[errmsg errID]=lastwarn; %detect if a warning occurred.  The warning will be displayed at the end of the WJGLM routine.
if strcmp(errID, 'MATLAB:singularMatrix')
    errorflag.TestSingularMatrix=errorflag.TestSingularMatrix+1;
    FSTAT=NaN;
end
if strcmp(errID, 'MATLAB:illConditionedMatrix') || strcmp(errID, 'MATLAB:nearlySingularMatrix')
    errorflag.TestIllConditionedMatrix=errorflag.TestIllConditionedMatrix+1;
    FSTAT=NaN;
end

warning on MATLAB:singularMatrix
warning on MATLAB:illConditionedMatrix
warning on MATLAB:nearlySingularMatrix

%****DEFINE MODULES TO PERFORM BOOTSTRAP****;
%***DEFINE MODULE TO GENERATE BOOTSTRAP DATA AND CENTRE DATA****;
function [YB] = bootdat(Y, BHAT, BOBS, NX); 
F=1;
M=0;
tempVar=[];
for J=1:BOBS
    L=M+NX(J);
    tempVar=Y(F:L,:);
    BVAL=tempVar;
    for P=1:size(tempVar,1);
        RVAL=0;
        while RVAL ==0
            RVAL=rand(1);
        end;
        BVAL(P,:)=tempVar(ceil(size(tempVar,1)*RVAL),:);
    end;
    if J==1
        YB=BVAL;
    else
        YB=[YB; BVAL];
    end;
    M=L;
    F=F+NX(J);
end;

%****CENTRE THE BOOTSTRAP DATA****;
%center the bootstrap sample based on the means from the full sample
function [YB]=bootcen(YB1,BHAT,BOBS,NX);
M=0;
F=1;
YB=YB1;
for I=1:BOBS
    L=M+NX(I);
    MVAL=BHAT(I,:);
    for Q=F:1:L
        YB(Q,:)=YB1(Q,:)-MVAL;
    end;
    M=L;
    F=F+NX(I);
end;

%****DEFINE MODULE TO COMPUTE BOOTSTRAP STATISTIC****;
function [FSTATB, errorflag] = bootstat(YB,OPT1,R, BOBS, WOBS, WOBS1, NTOT, NX, PER, X, SEED, errorflag);
[MUHATB, BHATB, BHATBW, YTB, DFB, errorflag] = mnmod(YB, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
[SIGMAB,STDIZERB] = sigmod(YTB,X,BHATBW,DFB,BOBS,WOBS,WOBS1,NX);
[FSTATB,DF1B,DF2B, errorflag] = testmod(SIGMAB,MUHATB,R,DFB, BOBS, WOBS, WOBS1, errorflag);

%****define module to compute bootstrap effect size****;
function [MULTP,EFFSZB,errorflag]=bootes(YB,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag);
[MUHATB, BHATB, BHATBW, YTB, DFB, errorflag] = mnmod(YB, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
[SIGMAB,STDIZERB] = sigmod(YTB,X,BHATBW,DFB,BOBS,WOBS,WOBS1,NX);
[MULTP,EFFSZB]=wjeffsz(BOBS,WOBS,WOBS1,NTOT,NX,LOC1,LOC2,SCALE,OPT1,PER,R,MUHATB,STDIZERB);

%****compute measure of effect size and bootstrap confidence interval****;
function [MULTP,EFFSZ]=wjeffsz(BOBS,WOBS,WOBS1,NTOT,NX,LOC1,LOC2,SCALE,OPT1,PER,R,MUHAT,STDIZER);
if OPT1==0 
    MULTP=1;
end;
if OPT1==1
    if SCALE==0
        MULTP=1;
    end;
    if SCALE==1
        if PER ~= 0
            cut=sqrt(2)*erfinv(2*PER-1); %probit per wikipedia https://en.wikipedia.org/wiki/Probit 12/6/2015
            a=zeros(1,3);
            a(1,1)=cut;
            a(1,3)=-1*cut;
            %**fun defines the probability density function of the standard normal distribution**;
            fun = @(Z) (Z.^2).*(1/(sqrt(2*3.141592653589793)).*(exp(-(1/2).*Z.^2)));
            int=integral(fun,cut,-cut);
            b=sum(int);
            winvar=b+PER*(2*(cut^2));
            MULTP=sqrt(winvar);
        else
            MULTP=1;
        end;
    end;
end;
num=R*MUHAT;
if LOC1==99
    stdz=1;
end;
if LOC1==0
    r2=R.^2;
    rvec=reshape(r2',[],BOBS*WOBS);
    stdz1=diag(sqrt(diag(STDIZER)));
    stdz=(rvec*stdz1*rvec')/sum(rvec); %average of square root of variances
end;
if LOC1>0
    if LOC1 <99
        loc=LOC1*WOBS-(WOBS-LOC2);
        stdz1 = STDIZER(loc,loc);
        if stdz1> 0
            stdz=sqrt(stdz1);
        end;
        if stdz1==0
            stdz=.00001;
        end;
    end;
end;

if (length(num)>1) || (WOBS > 1)
    %Effect sizes only available for between-group contrasts, pending further research by Dr. Lix.
    EFFSZ=NaN; 
else
    EFFSZ=MULTP*(num/stdz);
end;


% EFFSZ=zeros(length(num),size(stdz,2));
% for iCol=1:size(stdz,2)
%     EFFSZ(:,iCol)=num./stdz(:,iCol);
% end;

% ?****print RESULTS****;
% ?print 'Welch-James Approximate DF Solution';
% ?if OPT1=0 then print 'Least Squares Means & Variances';
% ?if OPT1=1 then print 'Trimmed Means & Winsorized Variances';
% ?if OPT1=1 then print 'Percentage of Trimming:'PER[format=4.2];
% ?if OPT2=0 then if OPT3=0 then print 'F Distribution Critical Value';
% ?if OPT2=1 then if OPT3=0 then print 'Bootstrap Critical Value for Single Test Statistic';
% ?if OPT2=1 then if OPT3=0 then do;
% ?print 'Number of Bootstrap Samples:' numsim_b[format=5.0],;
% ?print 'Starting Seed:' SEED[format=15.0],;
% ?end;
% ?print'Contrast Matrix:';
% ?print r[format=4.1],;
% ?MUHAT=MUHAT`;
% ?print 'Mean Vector:';
% ?print MUHAT[format=10.4],;
% ?print 'Sigma Matrix:';
% ?print SIGMA[format=10.4],;
% ?reslab={"Test Statistic" "Numerator DF" "Denominator DF" "p-value"};
% ?if OPT3=0 then do;
% ?print 'Significance Test RESULTS:';
% ?print RESULTS[rowname=reslab format=10.4]/;
% ?end;
% ?if OPT3=1 then do;
% ?cilev = (1 - ALPHA)*100;
% ?print 'Effect Size:' EFFSZ[format=10.4];
% ?if SCALE=0 then print 'No Scaling Factor';
% ?if SCALE=1 then print 'Scaling Factor is:' MULTP[format=6.3];
% ?if LOC1=0 then print 'Standardizer Is Square Root of Average Variance';
% ?if LOC1 > 0 then if LOC1 < 99 then do;
% ?print 'Standardizer is:';
% ?print 'LOC1:' LOC1[format=3.0];
% 
% ?print 'LOC2:' LOC2[format=3.0];
% ?end;
% ?if LOC1=99 then print 'No Standardizer';
% ?print 'Number of Bootstrap Samples:' numsim_es[format=5.0];
% ?print 'Starting Seed:' SEED[format=15.0];
% ?print 'Confidence Level (%):' cilev[format=5.1];
% ?print 'Lower Confidence Limit:' lcl[format=10.4];
% ?print 'Upper Confidence Limit:' ucl[format=10.4]/;
% ?end;

% %****define module to compute bootstrap critical value for FWR control of ADF statistics****;
% function bootcom;
% initial(C,U,Y,OPT1,PER,OPT2,numsim_b,numbsim_bc,OPT3,numsim_es,LOC1,LOC2,SCALE,ALPHA,SEED,R,X);
% [MUHAT, BHAT, BHATW, YT, DF, errorflag] = mnmod(Y, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
% [SIGMA,STDIZER] = sigmod(YT,X,BHATW,DF,BOBS,WOBS,WOBS1);
% 
% for i=1:size(R,1);
%     cm=R(i,);
%     [FSTAT,DF1,DF2, errorflag] = testmod(SIGMA,MUHAT,cm,DF, BOBS, WOBS, WOBS1, errorflag);
%     if i==1
%         cmmat=FSTAT;
%         DF1mat=DF1;
%         DF2mat=DF2;
%     end;
%     if i>1
%         cmmat=[cmmat, FSTAT];
%     end;
% end;
% 
% DF1mat=[DF1mat, DF1];
% DF2mat=[DF2mat, DF2];
% 
% for simloop=1:numsim_bc;
%     [YB1] = bootdat(Y, BHAT, BOBS, NX);
%     [YB]=bootcen(YB1,BHAT,BOBS,NX);
%     for q=1:size(R,1);
%         cm=R(q,:);
%         [FSTATB, errorflag] = bootstat(YB, OPT1, cm, BOBS, WOBS, WOBS1, NTOT, NX, PER, X, SEED, errorflag);
%         if q==1
%             frow=FSTATB;
%         else
%             frow=[frow, FSTATB];
%         end;
%     end;
%     if simloop==1
%         FMAT=frow;
%     else
%         FMAT=[FMAT; frow];
%     end;
% end;
% 
% fmax=zeros(numsim_bc,1);
% for k=1:numsim_bc;
%     fmax(k)=max(FMAT(k,:));
% end;
% 
% if isempty(ALPHA)
%     ALPHA = .05;
% end;
% fmax=sort(fmax);
% RESULTS=zeros(3,size(R,1));
% RESULTS(1,:)=cmmat;
% RESULTS(2,:)=DF1mat; 
% RESULTS(3,:)=DF2mat;
% qcrit=round((1-ALPHA)*numsim_bc);
% critv=fmax(qcrit);
% ****print RESULTS****;
% print 'Welch-James Approximate DF Solution';
% if OPT1=0 then print 'Least Squares Means & Variances';
% if OPT1=1 then print 'Trimmed Means & Winsorized Variances';
% if OPT1=1 then do;
% print 'Percentage of Trimming:';
% print PER[format=4.2];
% end;
% if OPT2=1 then print 'Bootstrap Critical Value for FWR Control';
% if OPT2=1 then print 'Number of Bootstrap Samples for Test Statistic Critical Value:';
% if OPT2=1 then print numsim_bc[format=4.0];
% if OPT2=1 then do;
% print 'Starting Seed:';
% print SEED[format=15.0],;
% end;
% print 'Contrast Matrix:';
% print r[format=4.1],;
% MUHAT=MUHAT`;
% print 'Mean Vector:';
% print MUHAT[format=10.4],;
% print 'Sigma Matrix:';
% print SIGMA[format=10.4],;
% reslab={"Test Statistic" "Numerator DF" "Denominator DF" "Significance"};
% print 'Significance Test RESULTS:';
% print RESULTS[rowname=reslab format=10.4],;
% if OPT2=1 then print 'Critical Value:';
% if OPT2=1 then print critv[format=5.2]/;
% finish;

function x = probf(f,d1,d2)
%Based on code from matrixlab-examples.com, implements SAS probf function.

x = 1;
% Computes using inverse for small F-values
if f < 1
    s = d2;
    t = d1;
    z = 1/f;
else
    s = d1;
    t = d2;
    z = f;
end
j = 2/(9*s);
k = 2/(9*t); 

% Uses approximation formulas
y = abs((1 - k)*z^(1/3) - 1 + j)/sqrt(k*z^(2/3) + j);
if t < 4
    y = y*(1 + 0.08*y^4/t^3);
end 

a1 = 0.196854;
a2 = 0.115194;
a3 = 0.000344;
a4 = 0.019527;
x = 0.5/(1 + y*(a1 + y*(a2 + y*(a3 + y*a4))))^4;
x = floor(x*10000 + 0.5)/10000; 

% Adjusts if inverse was computed
if f < 1
    x = 1 - x;
end 

x = 1 - x; %SAS probf provides percentile not tail end value.





