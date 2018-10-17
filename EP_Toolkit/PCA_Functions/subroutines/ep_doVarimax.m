function [Ak]= ep_doVarimax(A);

%[Ak]= ep_doVarimax(A);  - Compute varimax rotation for PCA
%  Gives essentially identical results to SAS command: PROC FACTOR rot = varimax
%
%Inputs
%  A		: Unrotated factor loadings matrix (variables, factors)
%
%Outputs
%  Ak		: Rotated factor pattern matrix
%
%History
%
%  Cureton, E. E. & Mulaik, S. A. (1975).  The weighted varimax rotation and the promax rotation.
%  Psychometrika, 40(2) 183-195.
%
%  Kaiser, H. F. (1959).  Computer program for varimax rotation in factor analysis.
%  Educational and Psychological Measurement, 19(3) 413-420.
%
%  Harman, H. H.  (1976).  Modern factor analysis, 3rd edition.  Chicago:University of Chicago Press.
%
%  by Joseph Dien (4/99)
%  jdien07@mac.com
%
%  modified (9/30/00) JD
%  Modified to allow Kaiser normalization to be turned off.
%
%  bugfix (2/27/01) JD
%  Fixed bug in Kaiser correction code (not undoing it when turned on)
%
%  modified (3/1/01) JD
%  Added direct rotation of factor scores.
%
%  bugfix (3/24/01) JD
%  Factor scores reordered and flipped when factor loadings are.
%
%  modified (5/27/01) JD
%  Factor scoring coefficients added.
%
%  modified (8/14/01) JD
%  Deleted rotation of V to minimize effect of cumulative rounding errors
%
%  modified (7/26/02) JD
%  Added non-convergence warning.
%
%  modified (10/22/02) JD
%  Lowered minimum criteria for rotation to .00001 since Kaiser suggested
%  .00116 was too loose and was sometimes not rotating when it should.
%  Added warning if no rotation occurred.
%
%  modified (2/15/04) JD
%  Added accomodation for covariance loading option
%
%  modified (01/09/05) JD and Dan Beal
%  Added Weighted Varimax option (Cureton and Mulaik, 1975)
%
%  modified (2/4/08) JD
%  Removed manual rotation of factor scores as it turned out to be less
%  accurate.
%
%  modified (2/12/08) JD
%  Moved factor loading normalization and factor loading ordering
%  to doPCA so that they can be applied to all rotations.
%
%  modified (2/27/08) JD
%  Random starting rotation procedure added to avoid local minima.
%
%  bugfix (3/24/08) JD
%  Fixed orientation problem with Ak when only one factor.
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

NumReps=10;
NumFacs=size(A,2);
NumVars=size(A,1);
FacPatReps=zeros(NumReps,NumVars,NumFacs);
ThReps=zeros(NumReps,NumFacs,NumFacs);
PhiReps=zeros(NumReps,NumFacs,NumFacs);
fReps=zeros(NumReps,1);
rotTot=zeros(NumReps,1);
didRot=zeros(NumReps,1);
for repetition =1:NumReps

    MaxRotations = 1000;	%Maximum number of rotations to do before calling it a day
    C = sum((A.^2)')';      %vector of communalities
    NUM_FAC = size(A,2);	%Number of factors retained
    NUM_VAR = size(A,1);	%Number of variables
    Ak=A;
    
    [T,R]=qr(randn(NUM_FAC,NUM_FAC)); %random initial rotation
    Ak=A*inv(T)';

    %  Start rotation

    counter = NUM_FAC * (NUM_FAC -1)/2;  %initialize counter
    nrot = 0;	%initialize rotation counter
    rotstats = [];
    rotated = 0;

    while (counter > 0) && (nrot < MaxRotations)  %do rotations until has gone through all w/o rotating or exceeded max rotations.
        for f1 = 1:(NUM_FAC-1)
            for f2 = (f1 + 1):NUM_FAC
                A1 = Ak(:,f1);
                A2 = Ak(:,f2);
                Ar = [A1 A2];
                ru = A1.^2 - A2.^2;
                rv = 2 * A1 .* A2;
                rA = sum(ru);
                rB = sum(rv);
                rC = sum(ru.^2 - rv.^2);
                rD = 2*sum(ru .* rv);

                NUM = rD - (2 * rA * rB)/NUM_VAR;
                DEN = rC - (rA.^2 - rB.^2)/NUM_VAR;
                if abs(NUM/DEN) > .00001 %(would rotation be enough to bother with?)
                    rotated = 1;

                    G = sqrt(NUM^2+DEN^2);	%Following computational variant due Wingersky (Harman, p. 294)
                    cos4phi = DEN/G;
                    cos2phi = sqrt((1+cos4phi)/2);
                    cosphi = sqrt((1+cos2phi)/2);
                    sinphi = sqrt((1-cos2phi)/2);

                    if NUM < 0
                        sinphi = -abs(sinphi);
                    else
                        sinphi = abs(sinphi);
                    end;

                    rot = [cosphi -sinphi; sinphi cosphi];
                    Ar = Ar * rot;
                    Ak(:,f1) = Ar(:,1);
                    Ak(:,f2) = Ar(:,2);
                    counter = NUM_FAC * (NUM_FAC -1)/2;

                else
                    counter = counter -1; %no rotation occurred so reduce counter.
                end;
            end;
        end;
        nrot = nrot +1;	%Increment number of rotations.  It is expected that after a certain point just not worth it.
    end;

    FacPatReps(repetition,:,:)=Ak;
    fReps(repetition,1)=sum(var(Ak.^2),2); %The Varimax criterion
    rotTot(repetition)=nrot;
    didRot(repetition)=rotated;
end;

[B I]=sort(fReps); %Find the global maximum of the computed rotations.

Ak=squeeze(FacPatReps(I(NumReps),:,:));

if size(Ak,1)==1 %FacPat oriented wrong way when just one factor
    Ak=Ak';
end;

if (sum(didRot) == 0)
    if NumFacs == 1
        disp('Since there was only one factor, no rotation occurred.');
    else
        disp(['Warning - no rotation occurred.']);
    end;
end;

if (min(rotTot) >= MaxRotations)
    disp(['Warning - solution did not converge, meaning that it was not able to reach a stable solution within a reasonable number of rotations.']);
end;
