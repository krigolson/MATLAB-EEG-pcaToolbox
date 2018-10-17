function [FacPat, FacStr, Phi] = ep_doGPFrotation(A, Rotation, RotOpt);

%[FacPat, FacStr, phi] = ep_doGPFrotation(FacPat, Rotation, ROTOPT)  - 
%  Perform rotations using gradient projection algorithm.
%  A variety of rotations is available.
%
%Inputs
%  A        : Initial factor loading matrix (variables, factors)
%  Rotation	: The rotation
%  RotOpt   : The rotation option, if exists, as in delta for Oblimin
%
%Outputs
%  FacPat	: Factor pattern matrix
%  FacStr	: Factor structure matrix
%  Phi  : Correlations between factors
%
%History
%
%  This code has been adapted from Matlab code made publicly available by
%  Coen A. Bernaards and Robert I. Jennrich with their permission.
%
%  Website: http://www.stat.ucla.edu/research
%
%  The software was presented in the following journal article:
%  Bernaards, C. A. & Jennrich, R. I. (2005). Gradient Projection
%  Algorithms and Software for Arbitrary Rotation Criteria in Factor Analysis. 
%  Educational and Psychological Measurement, 65(5), 676-696.
%
%  According to the website (2/3/08), this code has been made available on the
%  condition that:
%
% "The software may be distributed free of charge and used by anyone if
% credit is given. It has been tested in several situations, but it comes
% with no guarantees and the authors assume no liability for its use or
% misuse. By downloading and/or using this software you agree to these statements."
%
% The gradient projection rotation procedure was presented in:
% Jennrich, R. I. (2002). A simple general method for oblique rotation. Psychometrika, 67(1), 7-20.
% Jennrich, R. I. (2001). A simple general procedure for orthogonal rotation. Psychometrika, 66, 289?306.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  by Joseph Dien (2/08)
%  jdien07@mac.com
%
%  modified (2/27/08) JD
%  Random starting rotation procedure added to avoid local minima.
%  inv changed to pinv to deal with ill-conditioned matrices.
%
%     Copyright (C) 1999-2018  Joseph Dien
%     (to the portions of this code that differ from that distributed by
%      Bernaards and Jennrich)
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
for repetition =1:NumReps
    switch Rotation
        case 'QMAX'  %Quartimax
            T=randorth(NumFacs);
            [FacPat,Th,table,f]=GPForth(A,T,Rotation);
            Phi = diag(ones(size(A,2),1));
        case {'QMIN', 'OMIN', 'CRFE', 'MINE', 'IPSC', 'TIIC', 'GMIN', 'MMER'}
            T=randorth(NumFacs);
            [FacPat,Phi,Th,table,f]=GPFoblq(A,T,Rotation, RotOpt);
            %     case 'PMAX'  %Promax commented out because the implementation is not as good
            %         H=A.^3;
            %         H = abs(H) .* sign(A);		%add the signs back to the target matrix since even powers eliminate them
            %         [FacPat,Th]=procrustes(H,A);
            %         Phi = inv(Th)*inv(Th)';	%Correlations between oblique primary factors (Gorsuch, eq. 10.1.15, p. 215)
            %         FacStr = FacPat * Phi;	%factor structure matrix (Harman, eq. 12.19, p. 268)
        otherwise
            error('Unknown rotation specified');
    end;
    FacPatReps(repetition,:,:)=FacPat;
    ThReps(repetition,:,:)=Th;    
    PhiReps(repetition,:,:)=Phi;  
    fReps(repetition,1)=f; 
end;

[B I]=sort(fReps); %Find the global maximum of the computed rotations.

FacPat=squeeze(FacPatReps(I(1),:,:));
Th=squeeze(ThReps(I(1),:,:));
Phi=squeeze(PhiReps(I(1),:,:));
FacStr = FacPat * Phi;	%factor structure matrix (Harman, eq. 12.19, p. 268)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate random theta rotation matrix
function T=randorth(k)

[T,R]=qr(randn(k,k));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lh,Th,table,f]=GPForth(A,T,Rotation)

al=1;
L=A*T;
[f,Gq]=QMAX(L);
G=A'*Gq;
table=[];
for iter=0:500
    M=T'*G;
    S=(M+M')/2;
    Gp=G-T*S;
    s=norm(Gp,'fro');
    table=[table;iter f log10(s) al];
    if s<10^(-5),break,end

    al=2*al;
    for i=0:10
       X=T-al*Gp;
       [U,D,V]=svd(X,0);
       Tt=U*V';
       L=A*Tt;
       [ft,Gq]=QMAX(L);
       if ft<f-.5*s^2*al,break,end
       al=al/2;
       end
    T=Tt;
    f=ft;
    G=A'*Gq;
    end
    Th=T;
    Lh=A*T;

%--------------------------------------------------------------
%What follows is  a vgQ subroutine defining the quartimax
%criterion and gradient.
%--------------------------------------------------------------

function [q,Gq]=QMAX(L)

q=-norm(L.^2,'fro')^2/4;
Gq=-L.^3;

%--------------------------------------------------------------
%Turning to the oblique case what follows is the general GPF
%subrutine for oblique rotation. This is used without change.
%--------------------------------------------------------------

function [Lh,Phi,Th,table,f]=GPFoblq(A,T,Rotation, RotOpt)

al=1;
table=[];
Ti=inv(T);
L=A*Ti';

switch Rotation
    case 'QMIN'
        [f,Gq]=QMIN(L);
    case 'OMIN'
        [f,Gq]=OMIN(L, RotOpt);
    case 'CRFE'
        [f,Gq]=CRFE(L);
    case 'MINE'
        [f,Gq]=MINE(L);
    case 'IPSC'
        [f,Gq]=IPSC(L);
    case 'TIIC'
        [f,Gq]=TIIC(L);
    case 'GMIN'
        [f,Gq]=GMIN(L);
    case 'SMAX'
        [f,Gq]=SMAX(L);
    case 'MMER'
        [f,Gq]=MMER(L);
end;

G=-(L'*Gq*Ti)';
for iter=0:500
    Gp=G-T*diag(sum(T.*G));
    s=norm(Gp,'fro');
    table=[table;iter f log10(s) al];
    if s<10^(-5),break,end

    al=2*al;
    for i=0:10
       X=T-al*Gp;
       v=1./sqrt(sum(X.^2));
       Tt=X*diag(v);
       Ti=pinv(Tt);
       L=A*Ti';
       
       switch Rotation
           case 'QMIN'
               [ft,Gq]=QMIN(L);
           case 'OMIN'
               [ft,Gq]=OMIN(L, RotOpt);
           case 'CRFE'
               [ft,Gq]=CRFE(L);
           case 'MINE'
               [ft,Gq]=MINE(L);
           case 'IPSC'
               [ft,Gq]=IPSC(L);
           case 'TIIC'
               [ft,Gq]=TIIC(L);
           case 'GMIN'
               [ft,Gq]=GMIN(L);
           case 'SMAX'
               [ft,Gq]=SMAX(L);
           case 'MMER'
               [ft,Gq]=MMER(L);
       end;
       
       if ft<f-.5*s^2*al,break,end
       al=al/2;
       end
    T=Tt;
    f=ft;
    G=-(L'*Gq*Ti)';
    end
Th=T;
Lh=L;
Phi=T'*T;

%--------------------------------------------------------------
%What follows is the vgQ subroutine defining the quartimin
%criterion and its gradient.
%------------------------------------------------------------

function [q,Gq]=QMIN(L)

L2=L.^2;
[p,k]=size(L);
N=ones(k,k)-eye(k);
X=L2*N;
q=sum(sum(L2.*X))/4;
Gq=L.*X;

%--------------------------------------------------------------
%oblimin family
%--------------------------------------------------------------
function [q,Gq]=OMIN(L,RotOpt)

gm=RotOpt;

L2=L.^2;
[p,k]=size(L);
N=ones(k,k)-eye(k);
C=ones(p,p)/p;
X=(eye(p)-gm*C)*L2*N;
q=sum(sum(L2.*X))/4;
Gq=L.*X;

%--------------------------------------------------------------
%Crawford-Ferguson family.
%--------------------------------------------------------------
function [q,Gq]=CRFE(L)

kp=0;

L2=L.^2;
[p,k]=size(L);
M=ones(p,p)-eye(p);
N=ones(k,k)-eye(k);

X=(1-kp)*L2*N+kp*M*L2;
q=sum(sum(L2.*X))/4;
Gq=L.*X;

%--------------------------------------------------------------
%minimum entropy
%--------------------------------------------------------------
function [q,Gq]=MINE(L)

L2=L.^2;
q=-sum(sum(L2.*log(L2)))/2;
Gq=-L.*log(L2)-L;

%--------------------------------------------------------------
%Bentler's invariant pattern simplicity criterion.
%--------------------------------------------------------------
function [q,Gq]=IPSC(L)

L2=L.^2;
X=L2'*L2;
D=diag(diag(X));
q=-(log(det(X))-log(det(D)));
Gq=-L.*(L2*(inv(X)-inv(D)));

%--------------------------------------------------------------
%Comrey's tandem II criterion
%--------------------------------------------------------------
function [q,Gq]=TIIC(L)

L2=L.^2;
[p,k]=size(L);
X0=ones(p,p);
X1=L*L';
X2=L2*L2';
q=sum(sum(L2.*((X0-X1.^2)*L2)));
Gq=4*L.*((X0-X1.^2)*L2)-4*(X1.*X2)*L;

%--------------------------------------------------------------
%geomin
%--------------------------------------------------------------
function [q,Gq]=GMIN(L)

ep=.01;

L2=L.^2;
[p,k]=size(L);
u=ones(k,1);
e=exp(log(L2+ep)*u/k);
q=sum(e);
Gq=(2/k)*L.*(e*u')./(L2+ep);

%--------------------------------------------------------------
%target rotation
%--------------------------------------------------------------
% function [q,Gq]=vgQ(L)
% 
% global H
% 
% q=norm(L-H,'fro')^2;
% Gq=2*(L-H);

%--------------------------------------------------------------
%partially specified target rotation
%--------------------------------------------------------------
% function [q,Gq]=vgQ(L)
% 
% global H W
% 
% q=norm(W.*(L-H),'fro')^2;
% Gq=2*W.*(L-H);

%--------------------------------------------------------------
%simplimax
%--------------------------------------------------------------
% function [q,Gq]=SMAX(L)
% 
% global m;
% 
% L2=L.^2;
% [p,k]=size(L);
% ls=reshape(L2,p*k,1);
% ls=sort(ls);
% lsm=ls(m);
% W=(L2<=lsm);
% q=sum(sum(L2.*W))/2;
% Gq=L.*W;

%--------------------------------------------------------------
%infomax
%--------------------------------------------------------------
% The infomax rotation is commented out because the implementation is not
% as effective as that in EEGlab.
%
% function [q,Gq]=IMAX(L)
% 
% S=L.^2;
% 
% [p,k]=size(L);
% u=ones(p,1);
% v=ones(k,1);
% s1=S*v;
% s2=u'*S;
% s=s2*v;
% 
% P=S/s;
% p1=s1/s;
% p2=s2/s;
% 
% Q0=-sum(sum(P.*log(P)));
% Q1=-sum(p1.*log(p1));
% Q2=-sum(p2.*log(p2));
% 
% q=log(k)+Q0-Q1-Q2;
% 
% H=-(log(P)+1);
% al=u'*(S.*H)*v/s^2;
% G0=H/s-al*u*v';
% 
% h1=-(log(p1)+1);
% al1=s1'*h1/s^2;
% G1=h1*v'/s-al1*u*v';
% 
% h2=-(log(p2)+1);
% al2=h2*s2'/s^2;
% G2=u*h2/s-al2*u*v';
% 
% Gq=(2*L).*(G0-G1-G2);

%--------------------------------------------------------------
%McCammon minimum entropy ratio
%--------------------------------------------------------------
function [q,Gq]=MMER(L)

[p,k]=size(L);
S=L.^2;

u=ones(p,1);
v=ones(k,1);
M=u*u';
R=M*S;
P=S./R;
H=-(log(P)+1);

Q1=-sum(sum(P.*log(P)));
G1=H./R-M*(S.*H./R.^2);

s2=sum(S);
s=sum(s2);
p2=s2/s;
h=-(log(p2)+1);
al=h*p2';

Q2=-sum(p2.*log(p2));
G2=u*h/s-al*u*v';

q=log(Q1)-log(Q2);
Gq=2*L.*(G1/Q1-G2/Q2);

% --------------------------------------------------------------
% Promax rotation. What follows is not a vgQ subroutine, but
% rather a subroutine for procrustes rotation. This form of
% rotation is not based on optimizing a criterion and hence
% requires a special subroutine. The procrustes subroutine may
% be used for promax rotation by choosing an appropriate target
% H, for example the element-wise cube of a varimax rotation of
% A.
% --------------------------------------------------------------
function [L,T]=procrustes(H,A)

S=inv(A'*A)*A'*H;  %compute lambda which rotates to the reference vector
%(columns of H don't need to be normalized, Cureton & D'Agostino, 1983, pp 263-264)
d=sqrt(diag(inv(S'*S))); 
D=diag(d); %compute relationship between reference and primary axes
L=A*S*D; %compute factor pattern matrix 
T=inv(S*D)'; %compute Procrustean Transformation Matrix