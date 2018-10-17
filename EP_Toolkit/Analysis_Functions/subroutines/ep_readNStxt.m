function [data]=ep_readNStxt(fname);
%[data]=ep_readNStxt(fname);
%Reads in a NetStation statistical extraction text file.
%
%Inputs
%  fname  : The text file with the NetStation output in it.
%Outputs
%  data   : The extracted data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%History
%  by Joseph Dien (11/8/08)
%  jdien07@mac.com
%
%  modified (12/18/08) JD
%  added startLine parameter since it can differ for different versions of NetStation
%
%  modified (4/21/09) JD
%  automatically detects line where data starts.

if isempty(fname)
    [infid fname] = get_fid('r');
else
    infid = fopen(fname);
end

if infid == -1
    msg{1}='File error';
    [msg]=ep_errorMsg(msg);
    return
end;

C = textscan(infid, '%s','delimiter','\n');
fclose(infid);

if ~isempty(C{1}{size(C{1},1),:}) %make sure the final line is a blank.
    C{1}{size(C{1},1)+1,:}='';
end;

i=1;
while isempty(str2num(C{1}{i}(1:3))) || isempty(str2num(C{1}{i+1}(1:3)))
    i=i+1;
end;
    
for line=i:size(C{1},1)-1
    data(line-i+1,:)=strread(C{1}{line,:});
end;

data=data(:,2:end-1); %drop subject numbers and final blank column