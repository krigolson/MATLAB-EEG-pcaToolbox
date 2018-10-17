function ep_writeSubjectText(EPdata, fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ep_writeSubjectText(EPdata, fileName) - writes out subject specs text file
%
%	Writes out subject specs text file for file formats that do not support such
%	information.
%
%Inputs
%   EPdata         : Structured array with the data and accompanying information.  See readData.
%	fileName: file name for subject specs text file, including path and suffix.
%
%Outputs
%	Saves subject specs text file.
%   subject specs text file has one subject per row starting with the second.
%   Each column is a column header with the name of the subject spec.
%
% History:
%
% by Joseph Dien (11/16)
% jdien07@mac.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

numSpecs=length(EPdata.subjectSpecNames);
numSubs=length(EPdata.subNames);

outFID=fopen(fileName,'w');
if (outFID == -1)
    msg{1}='Error creating output file!';
    [msg]=ep_errorMsg(msg);
    return
end;

fprintf(outFID,'subject');
for iSpec=1:numSpecs
    fprintf(outFID,['\t' EPdata.subjectSpecNames{iSpec}]);
end;
fprintf(outFID,'\r');

for iSubject=1:numSubs
    fprintf(outFID,EPdata.subNames{iSubject});
    for iSpec=1:numSpecs
        fprintf(outFID,'\t%s',EPdata.subjectSpecs{iSubject,iSpec});
    end;
    fprintf(outFID,'\r');
end;
fprintf(outFID,'\r');
fclose(outFID);

