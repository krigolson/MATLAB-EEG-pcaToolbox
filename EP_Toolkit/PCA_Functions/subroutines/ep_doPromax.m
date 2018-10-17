function [PmxPat, phi] = ep_doPromax(VmxPat, k);

%[PmxPat, phi] = ep_doPromax(VmxPat, k)  - Compute promax rotation for PCA
%  Gives similar results to SAS command: PROC FACTOR rot = promax
%
%Inputs
%  VmxPat	: Varimax rotated factor loading matrix (variables, factors)
%  k		: Power to raise loadings to produce target matrix.  Higher power results in more oblique solutions.
%				4 better for simple factor structures, 2 for more complex, 3 is a good compromise (SAS default value)
%
%Outputs
%  PmxPat	: Factor pattern matrix
%  phi  : Correlations between factors
%
%History
%
%  With assistance from Lew Goldberg and Jack Digman.
%
%  First proposed by:
%  Hendrickson, A. E. & White, P. O. (1964).  Promax: A quick method for
%  rotation to oblique simple structure.  The British Journal of Statistical
%  Psychology, 17:65-70.
%
%  This algorithm uses the original procedure of rotating the reference
%  vectors and is therefore an "indirect Promax" rotation.
%
%  Normalization process attributed to:
%  Cureton, E. E., & D'Agostino, R. B. (1983). Factor Analysis: An applied approach. Hillsdale, NJ: Lawrence Erlbaum and Associates.
%
%  Harman, H. H.  (1976).  Modern factor analysis, 3rd edition.  Chicago:University of Chicago Press.
%
%  Gorsuch, R. L.  (1983).  Factor analysis, 2nd edition.  Hillsdale, NJ:Lawrence Erlbaum Associates.
%
%  Dillon, W. R. & Goldstein, M. (1984).  Multivariate analysis: Methods and applications.  New York:Wiley & Sons.
%
%  by Joseph Dien (4/99)
%  jdien07@mac.com
%
%  12/7/00 JD
%  Fixed error in promax algorithm.  Modified to output factor correlations and reference structure.
%  Given the same varimax solution, now produces identical results to SAS 6 promax output.
%
%  modified (4/3/01) JD
%  Added manual rotation of factor scores
%
%  bugfix 1/12/03 JD
%  Fixed column-normalization of H to be on absolute value.
%
%  modified (2/4/08) JD
%  Removed manual rotation of factor scores as it turned out to be less
%  accurate.

%     Copyright (C) 1999-2008  Joseph Dien
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

NUM_FAC = size(VmxPat,2);	%Number of factors retained
NUM_VAR = size(VmxPat,1);	%Number of variables

H = diag(1./sqrt(sum(VmxPat'.^2)))*VmxPat;	%Normalize rows - equalize sum of squares of each variable's loadings
H = H * diag(1./max(abs(H)));			%Column-normalize by highest absolute value in each column - tends to equalize factor sizes

H = H.^k;		%compute target matrix by taking higher power of factor loadings

H = abs(H) .* sign(VmxPat);		%add the signs back to the target matrix since even powers eliminate them

lambda = inv(VmxPat'*VmxPat) * VmxPat' * H;	%nonnormalized transformation matrix between starting factor matrix and reference vector (H & W, p. 66)

lambda = lambda*diag(1./sqrt(sum(lambda .^2)));		%Normalize the columns of lambda (H & W, p. 66)

psi = lambda' * lambda;		%Correlations of reference vectors (Harman, eq. 12.29, p. 273)                               

r = inv(psi);					%Inverse of psi

D = diag(1./sqrt(diag(r)));	%relationship between reference and primary axes (Gorsuch, p. 226)

T = lambda*inv(D);  %Procrustean Transformation Matrix (Gorsuch, eq. 10.2.9, p. 226)

phi = inv(T)*inv(T)';	%Correlations between oblique primary factors (Gorsuch, eq. 10.1.15, p. 215)

PmxPat = VmxPat * T; 		%factor pattern matrix (Gorsuch, eq. 10.2.2, p. 220)


