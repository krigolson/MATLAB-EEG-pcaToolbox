function [outFiles, pathname]=ep_getFilesUI(inputFormat)
% [outFiles, pathname]=ep_getFilesUI(suffix)
% Provides UI for selecting multiple files
%
%Input:
% inputFormat    : format of desired files
%
%Outputs:
% outFiles  : selected file names
% pathname  : path of selected files
%
%History
%  by Joseph Dien (10/13/13)
%  jdien07@mac.com
%
% bugfix 1/29/14 JD
% Fixed not returning full pathname on non-Mac systems.
%
% modified 12/8/14 JD
% On a PC, select the mff file rather than the directory containing the mff file.

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

outFiles=[];
pathname=[];

if ~ismac && any(strcmp(inputFormat,{'egi_mff_v1','egi_mff_v2'}))
    inDir = uigetdir(pwd,'Select an mff file within the directory that contains the mff files you wish to select.');
    if isempty(inDir)
        return
    end
    inDir=fileparts(inDir);
    temp=dir(inDir);
    temp2={temp.name};
    temp3=regexp(temp2,'.*mff$');
    theFiles=cell(0);
    for i=1:length(temp3)
        if ~isempty(temp3{i})
            [activeDirectory, name, ext] = fileparts(temp2{i});
            theFiles{end+1,1}=[name ext];
        end;
    end;
    
    if isempty(theFiles)
        disp('No mff files present.  Make sure to select an mff file within the directory that contains the mff files you wish to select.  Also make sure the name of the mff file/folder ends in .mmf so it can be identified.');
        return
    end;
    
    if length(theFiles) > 1
        [selected,ok] = listdlg('PromptString',['Choose mff files.'],'ListString',theFiles);
        outFiles=theFiles(selected);
    else
        outFiles=theFiles;
    end;
    
    pathname=[inDir filesep];
else
    [outFiles, pathname] = uigetfile('*.*','Open:','MultiSelect','on');
    if ~iscell(outFiles)
        temp=outFiles;
        outFiles=[];
        outFiles{1}=temp;
    end;
end;