function [sphereCoords, elecInSphere]=ep_sphereHead(nSides, eloc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [sphereCoords, elecInSphere]=ep_sphereHead(nSides, eloc);
% Fits a sphere to the electrode coordinates.
%
%Inputs
%   nSides  : Number of sides to the sphere representing the head
%   eloc    : The The electrode location information, one for each channel (see readlocs header).  By apparent EEGlab convention, a row vector. 
%                 eloc is the same length as the channels, with REG channels having a blank entry.
%
%Outputs
%   elecInSphere   : The sphere vertices corresponding to the electrodes.
%   sphereCoords   : Coordinates for the all the vertices of the interpolated sphere.
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

[X,Y,Z] = sphere(nSides);
sphereCoords=[X(:),Y(:),Z(:)];
elecCoords=[[eloc.X]',[eloc.Y]',[eloc.Z]'];
[center,radius] = sphere_cenrad(elecCoords);
elecCoords=[elecCoords(:,1)-center(1),elecCoords(:,2)-center(2),elecCoords(:,3)-center(3)];
elecInSphere = dsearchn(sphereCoords,elecCoords);

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