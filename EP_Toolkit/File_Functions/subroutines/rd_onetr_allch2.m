function trialdata = rd_onetr_allch2(fid, endian, cell_data_offset, trial_num, nchan, npoints, seek_start, start_samp, stop_samp)
% trialdata = rd_onetr_allch2(fileid, endian, cell_data_offset, trial_num, nchan, npoints, seek_start, start_samp, stop_samp)
%
%   Reads the raw (untransformed) data from all channels for a single trial 
%   from a single cell in the EGIS session file referred to by the file_id.
%
%   Returns channels in columns

% Modification history:
%
%  7/20/95 PJ -- version 1.1  added command line parameter 'seek_start' which enables one
%				 to read directly from the last location accessed in the file
%				instead of seeking from the beginning of the file.
%				seek_start = 'cof' for reading from current position and seek_start = 'bof'
%				for searching from beginning and then reading
%
%				Also added the optional parameters 'start_samp' and 'stop_samp' 
%				which allow one to extract a subset of samples from the trial
%
%  2/15/95 PJ -- started work on module
%  
%  2/10/08 JD -- added endian support
%
%  8/20/08 JD -- modified so endian is obtained by input field

bytes_to_end_of_trial = 0;

if nargin <= 6
  trial_offset = cell_data_offset + (trial_num-1)*nchan*npoints*2;
  seek_start = 'bof';
  start_samp = 1; stop_samp = npoints;
  samp_offset = 0;
elseif nargin == 7
  if seek_start == 'bof'
    trial_offset = cell_data_offset + (trial_num-1)*nchan*npoints*2;
  else
    trial_offset = 0;
  end
  start_samp = 1; stop_samp = npoints;
  samp_offset = 0;
elseif nargin == 9
  if seek_start == 'bof'
    trial_offset = cell_data_offset + (trial_num-1)*nchan*npoints*2;
  else
    trial_offset = 0;
  end
  samp_offset = (start_samp - 1)*nchan*2;
  bytes_to_end_of_trial = (npoints - stop_samp)*nchan*2;
end

nsamps = stop_samp - start_samp + 1;
%trial_offset, samp_offset, nsamps
fseek(fid, trial_offset + samp_offset, seek_start);

trialdata = fread(fid, [nchan, nsamps],'int16',endian);

%seek to the end of trial in case only a subset of samples was read.
%This aligns the file pointer at the proper position for a subsequent
%read in which the 'seek_start' parameter is 'cof'.  Otherwise the
%subsequent read would read in the wrong chunk of data.

fseek(fid, bytes_to_end_of_trial, 'cof');

trialdata = trialdata';
