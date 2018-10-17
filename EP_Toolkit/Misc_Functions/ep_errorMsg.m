function [msg]=ep_errorMsg(msg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  [msg]=ep_errorMsg(msg);
%       issues formatted error message to the command line window and presents an amodal error message box.
%
%Inputs:
%  msg       : The error message in a cell array, with each line in a different cell.

%
%  Example:
%  [msg]=ep_errorMsg('Hello world.');

%History:
%  by Joseph Dien (10/4/10)
%  jdien07@mac.com
%
% modified 11/7/16 JD
% Error message appears to the side of the main pane.
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

global EPmain
beep();

disp(' ');
disp('**************************************************************');
for i=1:length(msg)
    disp(msg{i});
end
disp('**************************************************************');
disp(' ');

errFig=errordlg(msg,'EP Toolkit Error');
set(errFig,'Position',[201 EPmain.scrsz(4)-70 397 70]);
msg=[];