function [status]=wt_PCAave_hdr_v(fid,infhdr,inchdr,ename,cnames,fcom,ftext);
%[status]=wt_PCAave_hdr_v(fid,fhdr,chdr,ename,cnames,fcom,ftext)
%
%This function creates a new EGIS average file header in the file with id fid.
%The structures which are passed to this routine are the same as
%are returned from function rd_ses_hdr_v (c.f.) and which are
%passed to wt_ses_hdr_v.  
%This routine reformats the file header and recalculates its size.
%The total number of elements written is returned.
%Note that this number is NOT the header size in bytes.
%
%NOTE: 	The NObs (NTrials) and NAvg (number of trials in average) fields in 
%		the Average Cell Specifics must be set prior to calling this routine
%		The same is true for the LastDone field of fhdr (Average File Description)
%
%See also:rd_ses_hdr_v.m

%  Modification history:
%	4/23/00 JD  Dropped mods to average file (like scalebins to 500) and automatic addition of channel
%	4/22/00	JD  Changed ~= to ~isempty to avoid warning messages.
%				Created the file based on wt_ave_hdr_v
%   1/11/01  JD  Delete ses hdr and minimum LSpec code
%
% modified (1/27/08) JD
% EGIS files always written as big-endian since NetStation makes this
% assumption.

ave_hdr_offsets_v;
fhdr = infhdr;
chdr = inchdr;

%
% Perform necessary transformations on the header
%

%% Put in the proper BytOrd value -- this is changed by wt_csdm_hdr_v, 
%% so it has to be changed back when writing files derived from csdm files

	fhdr(BytOrd) = 16909060;

%% Modify Average File Description
	fhdr(HdrVer) = -1;
	
%% Recompute header length

%disp('Computing length variables in header of average EGIS file...')
L = ((chdr(:,LSpec).*chdr(:,NObs))'+90);
fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))) = L';
nlcellhdr = fhdr(NCells)*2;
lhdr = 130+nlcellhdr+ fhdr(NChan) * 4 + sum(L)+fhdr(LComment)+fhdr(LText);
fhdr(LPad) = (512-rem(lhdr,512));
lhdr = lhdr+fhdr(LPad);
fhdr(LHeader) = lhdr;
fhdr(LData) = sum((chdr(:,NPoints).*chdr(:,NObs))*fhdr(NChan))*2;

status=fseek(fid,0,'bof');
status=fwrite(fid,fhdr(BytOrd),'int32','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(HdrVer),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message), error(message),end

status=fwrite(fid,fhdr(LHeader),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(LData),'int32','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,abs(ename),'char','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(6,1),'int16','ieee-be');  % empty fields -- run date and time in ses files
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(LastDone:BaseDur),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(3,1),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(NCells:BrdGain),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(fhdr(NChan)*2, 1),'int16','ieee-be');  % fill with zeros
message=ferror(fid);
if ~isempty(message),error(message),end

for loop=1:fhdr(NCells)
	status=fwrite(fid,chdr(loop,CellID),'int16','ieee-be');
	message=ferror(fid);
	if ~isempty(message),error(message),end
	status=fwrite(fid,abs(cnames(loop,:)),'char','ieee-be');
	message=ferror(fid);
	if ~isempty(message),error(message),end
	lastrow=5+((chdr(loop,LSpec)/2)*chdr(loop,NObs));
	status=fwrite(fid,chdr(loop,NObs:lastrow),'int16','ieee-be');
	message=ferror(fid);
	if ~isempty(message),error(message),end
end

status=fwrite(fid,fcom,'char','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,ftext,'char','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,zeros(fhdr(LPad),1),'char','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status = 0;
disp(['Successfully wrote the average file header. ' '(' num2str(fhdr(LHeader)) ' bytes long)']);
