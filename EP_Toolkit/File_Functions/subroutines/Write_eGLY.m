function [SUCCESS] = Write_eGLY(FileName, header_array, CateNames, CatLengths, EventCodes, segHdr, eventData, trialData)%  [SUCCESS] = Write_eGLY(FileName, header_array, CateNames, CatLength, EventCodes, segHdr, eventData, trialData)% 		This function will write data into a segmented NetStation Raw File version 3.%%	input:		FileName		.. filename to be written%	    		header_array	.. differs between versions, read code for details%				CateNames		.. category names%               CatLengths      .. length of category names%				eventCodes		.. if NEvent (from header_array) != 0, then array of 4-char event names%               segHdr          .. condition codes and time stamps for each segment% 				eventData		.. if NEvent != 0 then event state for each sample, else 'none'% 				trialData		.. matrix of EEG data.  numChannels x (numSamples x numSegments)%%	output:		SUCCESS		    .. result of function call.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 		Copyright:		2000 Electrical Geodesics, Inc.  ALL RIGHTS RESERVED% 		Property of:	Electrical Geodesics, Inc.% 						Eugene, OR 97403%% 		Programmer:		Joseph Dien% 		Date:			March 2003%		History:		Altered from Read_eGLY and UGLYFileWriter.%%		Known Bugs:		none%%%       3/10/09 Changed cateNames to cell array to allow for names with different lengths.%       8/25/10 Now allows for unsegmented files.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fid = fopen(FileName, 'w','ieee-be');if (fid == -1)    error('ERROR:	bad fid:  fid = %g  -- You probably have an incorrect path or filename \n', fid);    returnendversion		= header_array(1);if bitand(version,1) == 0    segmented = 0;else    segmented = 1;end;precision = bitand(version,6);if precision == 0    error('File precision is not defined.');end;if length(header_array) ~= 17    error('header array is the wrong length');end;NumCategors=header_array(14);NSegments=header_array(15);NSamples=header_array(16);NEvent =header_array(17);%		write header...fwrite(fid, header_array(1),'int32'); %versionfwrite(fid, header_array(2),'int16'); %yearfwrite(fid, header_array(3),'int16'); %monthfwrite(fid, header_array(4),'int16'); %dayfwrite(fid, header_array(5),'int16'); %hourfwrite(fid, header_array(6),'int16'); %minutefwrite(fid, header_array(7),'int16'); %secondfwrite(fid, header_array(8),'int32'); %millisecondfwrite(fid, header_array(9),'int16'); %sampling ratefwrite(fid, header_array(10),'int16'); %Number of channelsfwrite(fid, header_array(11),'int16'); %Gainfwrite(fid, header_array(12),'int16'); %bitsfwrite(fid, header_array(13),'int16'); %rangeif segmented    fwrite(fid, header_array(14),'int16'); %number of categories        for j = 1:NumCategors        fwrite(fid, CatLengths(j),'int8'); %length of category name        for i = 1:CatLengths(j)            fwrite(fid, CateNames{j}(i),'char'); %category name        end    end    fwrite(fid, header_array(15),'int16'); %Number of segmentsend;fwrite(fid, header_array(16),'int32'); %Number of samples per segmentfwrite(fid, header_array(17),'int16'); %Number of events per segmentfor j = 1:NEvent    fwrite(fid, EventCodes(j,1:4),'char'); %event codesendfor j = 1:NSegments    if segmented        fwrite(fid,segHdr(j,1),'int16'); %cell        fwrite(fid,segHdr(j,2),'int32'); %time stamp    end;    if (NEvent == 0)        switch precision            case 2                fwrite(fid,trialData(:,((j-1)*NSamples+1):j*NSamples),'int16');            case 4                fwrite(fid,trialData(:,((j-1)*NSamples+1):j*NSamples),'single');            case 6                fwrite(fid,trialData(:,((j-1)*NSamples+1):j*NSamples),'double');        end    else        data = cat(1,trialData(:,((j-1)*NSamples+1):j*NSamples), eventData(:,((j-1)*NSamples+1):j*NSamples));        switch precision            case 2                fwrite(fid,data,'int16');            case 4                fwrite(fid,data,'single');            case 6                fwrite(fid,data,'double');        end    endend%%%%%%%%%%%%%%%%%%%%%status = fclose(fid);SUCCESS = 1;if (status == -1)    SUCCESS = 0;endreturn