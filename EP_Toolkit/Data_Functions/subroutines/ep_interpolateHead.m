function [sphereValues]=ep_interpolateHead(dataIn, elecInSphere, sphereCoords);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [sphereValues]=ep_interpolateHead(dataIn, nSides, eloc);
% Interpolates the voltage values over a sphere representing the head.
%
%Inputs
%   dataIn  : The data matrix, where first dimension being channels.  The nature of the other dimensions don't matter.
%   elecInSphere   : The sphere vertices corresponding to the electrodes.
%   sphereCoords   : Coordinates for the all the vertices of the interpolated sphere.
%
%Outputs
%   sphereValues   : Voltage values for all the vertices of the interpolated sphere.
%
% History:
%
% by Joseph Dien (8/11/14)
% jdien07@mac.com
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

uniqueSphereVertices=unique(elecInSphere);
elecValues=zeros(length(uniqueSphereVertices),1);

for iVertex=1:length(uniqueSphereVertices)
    theVertex=uniqueSphereVertices(iVertex);
    elecValues(iVertex)=mean(dataIn(find(elecInSphere==theVertex)));
end;
nonElecSphereVertices=setdiff([1:size(sphereCoords,1)],uniqueSphereVertices);

W=ft_sphericalSplineInterpolate(sphereCoords(uniqueSphereVertices,:)',sphereCoords(nonElecSphereVertices,:)');

sphereValues=zeros(length(sphereCoords),1);
sphereValues(uniqueSphereVertices)=elecValues;
sphereValues(nonElecSphereVertices) = W * sphereValues(uniqueSphereVertices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [center,radius] = SPHERE_CENRAD(data)
% Fits a sphere to a set of 3D data and calculates the center.
% Vibor Paravic, IORE
%713-793-1725
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/17811

function [center,radius] = sphere_cenrad(data)

x = data(:,1);
y = data(:,2);
z = data(:,3);


data = [x.^2+y.^2+z.^2 x y z];
b = ones(length(x),1);
%coeffs = inv(data'*data)*data'*b;

coeffs = data\b; % per suggestion by Vibor Paravic

center = [-coeffs(2)/(2*coeffs(1)) -coeffs(3)/(2*coeffs(1)) -coeffs(4)/(2*coeffs(1))];

radius = sqrt(1/coeffs(1) + center(1)^2 + center(2)^2 + center(3)^2);