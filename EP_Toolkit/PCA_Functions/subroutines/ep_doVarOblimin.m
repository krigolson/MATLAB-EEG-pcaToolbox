function [FacPat, FacStr, FacCor] = ep_doVarOblimin(A);

%[FacPat, FacStr, phi] = ep_doVarOblimin(A)  - Compute promax rotation for PCA
%  Tries to find the optimal gamma value for Oblimin using a measure of factor complexity.
%
%Inputs
%  A        : Unrotated factor loading matrix (variables, factors)
%
%Outputs
%  FacPat	: Factor pattern matrix
%  FacStr	: Factor structure matrix
%  FacCor   : Correlations between factors
%
%History
%
%  Direct Oblimin:
%  Jennrich, R. I. & Sampson, P. F. (1966). Rotation for simple loadings. Psychometrika, 31(3), 313-323.
%
%  Original suggestion to systematically vary gamma values due to:
%  Gorsuch, R. L.  (1983).  Factor analysis, 2nd edition.  Hillsdale, NJ:Lawrence Erlbaum Associates.
%
%  Index of factor complexity:
%  Hofmann, R. J. (1977). Indices descriptive of factor complexity. The Journal of General Psychology, 96(1), 103-110.
%
%  by Joseph Dien (2/08)
%  jdien07@mac.com
%
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

gamma=[-4 -3 -2 -1 0 .2 .4 .6 .8];

VarComp=zeros(length(gamma),1);
FacComp=zeros(length(gamma),1);
for theGamma = 1:length(gamma)
    [FacPat, FacStr, FacCor] = ep_doGPFrotation(A, 'OMIN', gamma(theGamma));
    
    D2=diag(diag(FacPat*FacPat'));
    X=pinv(sqrt(D2))*FacPat;
    VarComp(theGamma)=mean(1./sum((X.^4),2));
    
    C2=diag(diag(X'*X));
    Z=X*pinv(sqrt(C2));
    FacComp(theGamma)=mean(1./sum((Z.^4),1));  
end;

[M I]=min(FacComp);

[FacPat, FacStr, FacCor] = ep_doGPFrotation(A, 'OMIN', gamma(I));




