function ep_viewANOVA
% ep_viewANOVA -
% Presents the loaded ANOVA data as organized by the ANOVA factors.
% It has no active controls at this point and so is just a figure.

%History
%  by Joseph Dien (6/28/09)
%  jdien07@mac.com
%
%  bugfix 9/5/09 JD
%  Location of windows appearing partly off screen on some systems due to offset from having Dock (on Mac) or Task Bar
%  (on PC) at the bottom of the screen.
%
%  bugfix 1/19/11 JD
%  Fixed View ANOVA function not matching up ANOVA levels with data columns correctly when a subset of the columns are
%  selected.
%
%  bugfix 11/1/13 JD
%  Fixes font sizes on Windows.
%
%  bugfix 1/12/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.

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

scrsz = EPmain.scrsz;
windowHeight=700;

figure('Name', 'ANOVA Window', 'NumberTitle', 'off', 'Position',[201 scrsz(4)-windowHeight 700 windowHeight]);
colormap jet;

tableData=num2cell(EPmain.anova.data.data);

for i=1:size(EPmain.anova.data.betweenLvl,2)
    for j=1:size(EPmain.anova.data.betweenLvl,1)
        tableData{j,i+size(EPmain.anova.data.data,2)}=EPmain.anova.data.betweenLvl{j,i};
    end;
end;

totalFactorLevels=EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1;

factorLevels=cell(0);
for theFactor=1:length(EPmain.anova.data.factor)
	if ~isempty(EPmain.anova.data.factor{theFactor}) && ~isempty(EPmain.anova.data.levels{theFactor})
        factorLevels{end+1}=EPmain.anova.data.levels{theFactor};
    end;
end;

tableNames=EPmain.anova.data.columnNames;
repFactor=1;
repCount=1;
lvlCount=1;
for theFactor=length(factorLevels):-1:1
    for theCol=EPmain.anova.data.leftColumn:EPmain.anova.data.leftColumn+totalFactorLevels-1
        tableNames{theCol}=[factorLevels{theFactor}(lvlCount) '|' char(tableNames{theCol})];
        repCount=repCount+1;
        if repCount > repFactor
            repCount=1;
            lvlCount=lvlCount+1;
            if lvlCount > length(factorLevels{theFactor})
                lvlCount=1;
            end;
        end;
    end;
    repFactor=repFactor*length(factorLevels{theFactor});
end;
 
try
    EPmain.handles.anova.dataTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
        'Position',[10 windowHeight-600 650 500]);
catch
    h = uicontrol('Style','text','HorizontalAlignment','left','String', 'This function does not work with this version of Matlab.','FontSize',EPmain.fontsize,...
        'Position',[10 windowHeight-600 650 500]);
end;



