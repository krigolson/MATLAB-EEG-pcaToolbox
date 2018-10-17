function [status]=wt_ses_hdr_v2(fid,fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext);
%[status]=wt_ses_hdr_v2(fid,fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext)
%
%This function creates a new EGIS header in the file with id fid.
%The structures which are passed to this routine are the same as
%are returned from function rd_ses_hdr_v (c.f.) prepended with the
%output file id. The total number of elements written is returned.
%Note that this number is NOT the header size in bytes.
%
%See also:rd_ses_hdr_v.m
%
%  Modification history:
%	6/6/95	PJ	Fixed writing of LPad
%
% modified (11/20/08) JD
% EGIS files always written as big-endian since NetStation makes this
% assumption.

ses_hdr_offsets_v;
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

status=fwrite(fid,fhdr(RunDate:BrdGain),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))),'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,czeros,'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

status=fwrite(fid,cgains,'int16','ieee-be');
message=ferror(fid);
if ~isempty(message),error(message),end

for loop=1:fhdr(NCells)
	status=fwrite(fid,chdr(loop,1),'int16','ieee-be');
	message=ferror(fid);
	if ~isempty(message),error(message),end
	status=fwrite(fid,abs(cnames(loop,:)),'char','ieee-be');
	message=ferror(fid);
	if ~isempty(message),error(message),end
	lastrow=5+((chdr(loop,5)/2)*chdr(loop,2));
	status=fwrite(fid,chdr(loop,2:lastrow),'int16','ieee-be');
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
