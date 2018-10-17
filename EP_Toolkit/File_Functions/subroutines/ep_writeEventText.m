function ep_writeEventText(EPdata, fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ep_writeEventText(EPdata, fileName) - writes out events text file
%
%	Writes out events text file for file formats that do not support such
%	information.
%
%Inputs
%   EPdata         : Structured array with the data and accompanying information.  See readData.
%	fileName: file name for events text file, including path and suffix.
%
%Outputs
%	Saves events text file.
%   events text file has the following columns: type, sample, value, duration,
%   and for each key code data datatype description with a number after it
%   (e.g., key1).
%
% History:
%
% by Joseph Dien (9/14)
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
numKeys=0;

for iCell=1:numel(EPdata.events)
    for iEvent=1:length(EPdata.events{iCell})
        numKeys=max([length(EPdata.events{iCell}(iEvent).keys),numKeys]);
    end;
end;

outFID=fopen(fileName,'w');
if (outFID == -1)
    msg{1}='Error creating output file!';
    [msg]=ep_errorMsg(msg);
    return
end;

fprintf(outFID,'subject\tcell\ttype\tsample\tvalue\tduration');
for iKey=1:numKeys
    fprintf(outFID,'\tcode%d\tdata%d\tdatatype%d\tdescription%d',iKey,iKey,iKey,iKey);
end;
fprintf(outFID,'\r');

for iSubject=1:size(EPdata.events,1)
    for iCell=1:size(EPdata.events,2)
        for iEvent=1:length(EPdata.events{iSubject,iCell})
            fprintf(outFID,'%s\t%s\t%s\t%s\t%s\t%s',sprintf('%03.0f',iSubject),sprintf('%03.0f',iCell),num2str(EPdata.events{iSubject,iCell}(iEvent).type),num2str(EPdata.events{iSubject,iCell}(iEvent).sample),num2str(EPdata.events{iSubject,iCell}(iEvent).value),num2str(EPdata.events{iSubject,iCell}(iEvent).duration));
            numThisEventKeys=length(EPdata.events{iSubject,iCell}(iEvent).keys);
            for iKey=1:numThisEventKeys
                fprintf(outFID,'\t%s\t%s\t%s\t%s',num2str(EPdata.events{iSubject,iCell}(iEvent).keys(iKey).code),num2str(EPdata.events{iSubject,iCell}(iEvent).keys(iKey).data),num2str(EPdata.events{iSubject,iCell}(iEvent).keys(iKey).datatype),num2str(EPdata.events{iSubject,iCell}(iEvent).keys(iKey).description));
            end;
            for iKey=1:(numKeys-numThisEventKeys)
                fprintf(outFID,'\t\t\t\t');
            end;
            fprintf(outFID,'\r');
        end;
    end;
end;
fprintf(outFID,'\r');
fclose(outFID);

