%% Remove gradient artifacts from EEG recorded with fMRI
%
% amri_eeg_gac()
%
% Usage
%   [OUTEEG]=amri_eeg_gac(EEG,'key1','value1' ...)
%
% Inputs
%   EEG                 =   EEG data stored in an EEGLAB structure
%
% Keywords
%   'method'            =   'AAS': Average Artifact Subtraction [Allen 2000]
%                           'AAR': Average Artifact Regression 
%                           'MAS': Median Artifact Subtraction 
%                           'MAR': Median Artifact Regression
%                           'PCA': Principal Component Analysis [Liu 2011]
%                           'TSE': Talyor Series Expansion [Wan 2006]
%                           'BLF': Band-limited filter [Hoffmann 2000]
%                           {default: 'PCA'}
%   'correctby'         =   ['volume' | 'slice']
%                           {default: 'slice'}
%   'winsize'           =   moving window size in number of volumes or slices
%                           If set to be [], use all volumes or slices to compute 
%                           the average artifact template
%                           {default: 101}
%   'trigger.name'      =   scanner trigger name 
%                           {default: 'R128'}
%                           Example: 'R128',{'R128'},{'R128','R129'}
%   'trigger.type'      =   ['volume' | 'slice'] scanner trigger type 
%                           {default: 'slice'}
%   'fmri.nslice'       =   number of slices within each volume
%                           {default: NaN}
%   'fmri.nvolume'      =   number of volumes throughout the whole scan
%                           {default: NaN}
%   'fmri.tr'           =   volume tr in seconds
%                           {default: NaN, computed based on marker timing}
%   'fmri.te'           =   echo time (te) in seconds
%                           {default: NaN}
%   'verbose'           =   1|0
%
% Advanced options
%
%   'realign.flag'      =   run trigger realignment [0:NO; 1:YES]
%                           {default:0}
%   'realign.maxtshift' =   maximum time-shift for volume/slice re-alignment
%                           {default: 1}
%   'realign.refscan'   =   number of volume/slice to which all other volumes/slices 
%                           are re-aligned.
%                           {default: 10}
%   'realign.refchan'   =   number or name of the channel from which the
%                           recorde signal is used for volume re-alignment.
%                           If 'All', then all the channels will be used.
%                           {default: 1}
%                           Example:
%                           1, [1 2 3], 'Fp1','All',{'Fp1','Fp2'}
%   'blf.method'        =   'ifft': frequency-domain filter using ifft
%                           'butter': butterworth filter (???)
%                           'firls': least square linear phase FIR filter (???)
%   'detrend'           =   detrend option
%                           nan: no detrending
%                           0: dc detrend
%                           1: linear detrend
%   'upsample'          =   upsampling frequency (Hz) 
%                           rounded to a multiple of original sampling frequency
%                           {default: 5000}
%   'downsample'        =   downsampling frequency (Hz)
%                           {default: 250}
%   'highpass.method'   =   'firls': least square linear phase FIR filter
%                           'ifft': frequency-domain filter using ifft
%                           {default: 'ifft'}
%   'highpass.cutoff'   =   low cutoff frequency (Hz) of the highpass filter applied
%                           before removing gradient artifacts
%                           {default: 1}
%   'lowpass.method'    =   'firls': least square linear phase FIR filter
%                           'ifft': frequency-domain filter using ifft
%                           {default: 'ifft'}
%   'lowpass.cutoff'    =   high cutoff frequency (Hz) of the lowpass filter
%                           applied after removing gradient artifacts
%                           {default: 125 Hz}
%   'check'             =   check intermediate results in run time
%                           {default: 0}
%
% See also:
%  amri_sig_findpeaks, 
%  amri_sig_filtfft, 
%  amri_sig_corr, 
%  amri_sig_xcorr
%  amri_stat_iqr,
%  amri_stat_ttest, 
%  amri_stat_tcdf
%
% Version:
%   0.23
%
% Examples:
%   N/A
%
% Reference:
%   Liu Z, de Zwart JA, van Gelderen P, Kuo L-W, Duyn JH, Statistical
%   feature extraction for artifact removal from concurrent fMRI-EEG
%   recordings, NeuroImage (2011), doi:10.1016/j.neuroimage.2011.10.042

%% DISCLAIMER AND CONDITIONS FOR USE:
%     This software is distributed under the terms of the GNU General Public
%     License v3, dated 2007/06/29 (see http://www.gnu.org/licenses/gpl.html).
%     Use of this software is at the user's OWN RISK. Functionality is not
%     guaranteed by creator nor modifier(s), if any. This software may be freely
%     copied and distributed. The original header MUST stay part of the file and
%     modifications MUST be reported in the 'MODIFICATION HISTORY'-section,
%     including the modification date and the name of the modifier.
%
% CREATED:
%     Oct. 7, 2009
%     Zhongming Liu, PhD
%     Advanced MRI, NINDS, NIH

%% MODIFICATION HISTORY
%
% 0.00 - 10/07/2009 - ZMLIU - Implement IAR (Allen et al. NeuroImage, 2000)
% 0.01 - 10/16/2009 - ZMLIU - Channel-by-channel correction to save memory
% 0.02 - 10/17/2009 - ZMLIU - Implement aFIR (Wan et al. CNP, 2006)
% 0.03 - 10/19/2009 - ZMLIU - Use a moving window to compute artifact template
% 0.04 - 11/19/2009 - ZMLIU - use detrend(..,'constant') instead of 'linear'
%                           - fix bugs on scanner_event.name, type, ignore,
%                           - fix bugs on scanner_event.latencies, intervals
%                           - add 'realignment' keyword
% 0.05 - 12/15/2009 - ZMLIU - implement PCA based approaches
%                           - re-organize the code
% 0.06 - 12/19/2009 - ZMLIU - re-organized the code and add comments
% 0.07 - 12/21/2009 - ZMLIU - implement PCA corrected by slice
% 0.08 - 12/22/2009 - ZMLIU - rename 'TS' to 'TSE'
%                           - ensure signal continuity
%                           - piece-wise linear detrend
% 0.09 - 12/23/2009 - ZMLIU - fix a small bug
%                           - adjust a parameter for linear detrend
% 0.10 - 12/24/2009 - ZMLIU - add figures to check intermediate results
% 0.11 - 01/06/2010 - ZMLIU - add detrend option, default nan (no detrend)
%                           - track figure size/position
%                           - FFT and iFFT (Hoffmann 1999)
% 0.12 - 01/07/2010 - ZMLIU - include here subroutines (fitlfft, findpeaks) 
%                           - lowpass filter before/after artifact removal
%                           - replace filtfft by firls and filtfilt
%                           - design filters before processing
% 0.13 - 01/08/2010 - ZMLIU - add zero padding to bandstop ifft filter
%                           - add transition zeros to bandstop ifft filter
%                           - smooth bandstop ifft filter using cos, sin
%                           - implement iFFT for lowpass filter
%                           - implement iFFT for highpass filter
%                           - modify the way of handling incomplete epochs
% 0.14 - 03/17/2010 - ZMLIU - display 'help' when no input
%                           - use dwtest and autocorrelation to select
%                             artifactual component in the PCA-based method
% 0.15 - 03/20/2010 - ZMLIU - use 'PCA' as the default method
%                           - deal with incomplete epochs
%                           - round scanner event latency
% 0.16 - 06/18/2010 - ZMLIU - use amri_sig_findpeaks instead of findpeaks
%                           - use amri_sig_filtfft instead of amri_sig_filtfft
% 0.17 - 07/06/2010 - ZMLIU - use amri_sig_corr instead of corr
%                           - use amri_stat_iqr instead of iqr
%                           - use amri_stat_ttest instead of ttest
%                           - use amri_sig_xcorr instead of xcorr
%                           - use varargin instead of p#,v# pairs
% 0.18 - 07/07/2010 - ZMLIU - fix a minor bug re eeg.xmax in downsampling
% 0.19 - 07/08/2010 - ZMLIU - simplify text message
% 0.20 - 07/26/2011 - ZMLIU - add back the autocorrelation criterion for
%                             artifactual component determination.
%                             this criterion works well when scanner
%                             trigger and eeg sampling is not synchronized
% 0.21 - 11/10/2011 - ZMLIU - clear up code and comments before release
%        16/11/2011 - JAdZ  - v0.21 included in amri_eegfmri_toolbox v0.1
% 0.22 - 11/18/2011 - ZMLIU - fix a bug in filter design when only volume
%                             triggers are available and tr>=2 sec
%        2011/11/18 - JAdZ  - v0.22 included in amri_eegfmri_toolbox v0.1.1
% 0.23 - 07/30/2013 - ZMLIU - exclude the outlier epochs in data matrix
%                             applied with PCA
%                           - corrected my mis-understanding that SVD is
%                             robust against the presence of outliers
%                           - remove artifacts for volumes/slices with
%                             slightly inaccurate TR (briefly tested on
%                             Jen's data).
% 0.23.1 12/25/2013 -JDIEN  - fix for crash when final gradient artifact
%                             to close to end of EEG recording by simply
%                             not trying to correct it.
%%

function [outeeg] = amri_eeg_gac(eeg,varargin)


if nargin<1
    eval('help amri_eeg_gac');
    return
end 

%% ************************************************************************
% Defaults
% *************************************************************************

gac_method              =       {'PCA'};    % algorithm name

% trigger markers received from scanner

scanner_trigger.type    =       'slice';    % scanner trigger type (per slice or volume)
scanner_trigger.name    =       {'R128'};   % scanner trigger name 
scanner_trigger.count   =       0;          % total number of scanner triggers
scanner_trigger.index   =       [];         % array of indices to eeg.event
scanner_trigger.latency =       [];         % array of latencies in time points

% scanner events                            
                                            % specify how the artifact would be corrected 
scanner_event.type      =       'slice';    % (by volume or slice)? 
scanner_event.latency   =       [];         % latencies (in time points)
scanner_event.count     =       0;          % total number of scanner events
scanner_event.ignore    =       [];         % list of scanner events to be ignored
scanner_event.tshift    =       [];         % time shift of event latency
scanner_event.epochonset=       -5;         % epoch onset w.r.t. event latency

% fmri acquisition information

fmri.nslice             =       NaN;        % number of slices
fmri.nvolume            =       NaN;        % number of volumes
fmri.trsec              =       NaN;        % volume tr in seconds
fmri.trpnt              =       NaN;        % volume tr in sampling points

% moving window

moving_window.size      =       101;        % size (in number of scanner events)
moving_window.start     =       [];         % index to starting event of each window
moving_window.end       =       [];         % index to ending event of each window
moving_window.index     =       [];         % indices to all events within each window
moving_window.nowin     =       0.3;        % sec, the half window length to exclude neighboring epochs
moving_window.minsize   =       10;         % min size(in number of scanner events)

% parameters for scanner event latency realignment

realign.flag            =       0;          % realign (1) or not (0)
realign.maxtshift       =       1;          % range of time-shift in time points
realign.refscan         =       10;         % scanner event to be realigned to
realign.refchannel      =       1;          % channels used for realignment


% signal processing parameters

upsample_rate           =       1;          % upsample rate for pre-processing
downsample_rate         =       round(eeg.srate/250);          % downsample rate for post-processing

detrend_option          =       NaN;        % no detrend

% taylor series order 

ts.order                =       25;         % order of taylor series expansion

% filters

lowpass.method          =       'ifft';     % filter design method ['ifft'|'firls']
lowpass.filter          =       NaN;        % lowpass filter coefficient
lowpass.cutoff          =       125;        % cutoff frequency (Hz)
lowpass.trans           =       0.04;       % relative transition
lowpass.order           =       NaN;        % filter order

highpass.method         =       'ifft';     % filter design method ['ifft'|'firls']
highpass.filter         =       NaN;        % highpass filter coefficient
highpass.cutoff         =       NaN;        % cutoff frequency (Hz)
highpass.trans          =       0.15;       % relative transition
highpass.order          =       NaN;        % filter order

bandstop.method         =       'ifft';     % filter design method ['ifft'|'butter'|'firls']
bandstop.order          =       NaN;        % order of firls filter
bandstop.trans          =       0.15;       % relative transition (<=0.5)
bandstop.filter         =       NaN;        % filter coefficients
bandstop.width          =       1;          % bandwidth (Hz)
bandstop.lcutoff        =       NaN;        % lowcutoff freq
bandstop.hcutoff        =       NaN;        % highcutoff freq
bandstop.tstart         =       NaN;        % starting time of the time series to be filtered
bandstop.tend           =       NaN;        % ending time of the time series to be filtered

% setting for pca component selection

pca.selection           =       'auto';     % 'auto' | 'fixed'
pca.alpha               =       0.05;       % significance level

% checking intermediate results

checking                =       0;          % checking data
checkfig.position       =       NaN;        % figure position

% flags
flag_verbose            =       0;         % whether print out information

%% ************************************************************************
% Collect keyword-value pairs
% *************************************************************************

for i = 1:2:size(varargin,2) 

    Keyword = varargin{i};
    Value   = varargin{i+1};

    if ~ischar(Keyword)
        printf('amri_eeg_gac(): keywords must be strings'); 
        continue;
    end
    
    if strcmpi(Keyword,'method')
        
        if ischar(Value)
            gac_method{1}=Value;
        elseif iscell(Value)
            gac_method = Value;
        else
            printf('amri_eeg_gac(): invalid value for ''method''');
        end
        
    elseif strcmpi(Keyword,'correctby')
        
        if strcmpi(Value,'slice')
            scanner_event.type='slice';
        else
            scanner_event.type='volume';
        end

    elseif strcmpi(Keyword,'winsize') || strcmpi(Keyword,'win_size')
    
        if isnumeric(Value)
            moving_window.size=round(Value(1));
        elseif isempty(Value)
            moving_window.size=NaN;
        elseif isnan(Value)
            moving_window.size=NaN;
        else
            printf('amri_eeg_gac(): invalid ''winsize'' value');
        end

    elseif strcmpi(Keyword,'trigger.name') || ...
           strcmpi(Keyword,'trigger_name') || ...
           strcmpi(Keyword,'triggername')

        if ischar(Value)
            scanner_trigger.name{1}=Value;
        elseif iscell(Value)
            scanner_trigger.name=Value;
        else
            printf('amri_eeg_gac(): invalid value for ''trigger.name''');
        end

    elseif strcmpi(Keyword,'trigger.type') || ...
           strcmpi(Keyword,'trigger_type') || ...
           strcmpi(Keyword,'triggertype')
       
        if ~ischar(Value)
            printf('amri_eeg_gac(): ''trigger.type'' must be a string');
            return;
        else
            scanner_trigger.type = Value;
        end
        
    elseif strcmpi(Keyword,'ignore')
        scanner_event.ignore = Value;
        
    elseif strcmpi(Keyword,'fmri.nslice') || ...
           strcmpi(Keyword,'nslice') || ...
           strcmpi(Keyword,'numslice') || ...
           strcmpi(Keyword,'nrslice') || ...
           strcmpi(Keyword,'nrofslice') || ...
           strcmpi(Keyword,'nrofslices')

       if ~isempty(Value)
           fmri.nslice=round(Value);
       end
       
    elseif strcmpi(Keyword,'fmri.nvolume') || ...
           strcmpi(Keyword,'nvolume') || ...
           strcmpi(Keyword,'nvol') || ...
           strcmpi(Keyword,'numvolume') || ...
           strcmpi(Keyword,'numvol') || ...
           strcmpi(Keyword,'nrvolume') || ...
           strcmpi(Keyword,'nrvol') || ...
           strcmpi(Keyword,'nrofvolume') || ...
           strcmpi(Keyword,'nrofvol')
       
        if ~isempty(Value)
            fmri.nvolume=round(Value);
        end
        
    elseif strcmpi(Keyword,'fmri.tr') || ...
           strcmpi(Keyword,'tr')
        
       if ~isempty(Value)
            fmri.trsec = Value;
       end

    elseif strcmpi(Keyword,'detrend') || ...
           strcmpi(Keyword,'detrend_opt') || ...
           strcmpi(Keyword,'detrend_option')
       
       detrend_option = Value;
       
    elseif strcmpi(Keyword,'upsample') || ...
           strcmpi(Keyword,'upsample')
        
        if Value<=eeg.srate
            upsample_rate=1;
        else
            upsample_rate = round(Value/eeg.srate);
        end
        
    elseif strcmpi(Keyword,'downsample') || ...
           strcmpi(Keyword,'down_sample')
        
        if Value>=eeg.srate;
            downsample_rate=1;
        else
            downsample_rate = eeg.srate/Value;
        end
        
    elseif strcmpi(Keyword,'lowpass.method')
        
        lowpass.method = Value;
        
    elseif strcmpi(Keyword,'lowpass.trans')
        
        lowpass.trans = Value;
        

    elseif strcmpi(Keyword,'highcutoff') || ...   % compatible
           strcmpi(Keyword,'high_cutoff') || ...  % compatible
           strcmpi(Keyword,'lowpass.cutoff')
        
        if isnumeriic(Value)
            lowpass.cutoff=Value(1);
        else
            lowpass.cutoff = NaN;
        end
        
    elseif strcmpi(Keyword,'highpass.method')
        
        highpass.method = Value;
        
    elseif strcmpi(Keyword,'highpass.trans')
        
        highpass.trans  = Value;
        
    elseif strcmpi(Keyword,'lowcutoff') || ...   % compatible
           strcmpi(Keyword,'low_cutoff') || ...  % compatible
           strcmpi(Keyword,'highpass.cutoff')
        
        if isnumeric(Value)
            highpass.cutoff=Value(1);
        else
            highpass.cutoff=NaN;
        end
        
    elseif strcmpi(Keyword,'realign.flag') || ...
           strcmpi(Keyword,'realign') || ...
           strcmpi(Keyword,'realignment')
        
        if isnumeric(Value)
            realign.flag = Value(1);
        else
            printf('amri_eeg_gac(): ''realign.flag'' must a numeric value');
        end
        
    elseif strcmpi(Keyword,'realign.maxtshift') || ...
           strcmpi(Keyword,'maxtshift') || ...
           strcmpi(Keyword,'max_tshift')
       
        if ~isnumeric(Value)
            printf('amri_eeg_gac(): ''realign.maxtshift'' must be a numeric value');
            return;
        else
            realign.maxtshift=Value;
        end
        
    elseif strcmpi(Keyword,'realign.refscan') || ...
           strcmpi(Keyword,'realign.ref_scan') || ...
           strcmpi(Keyword,'align2scan') || ...
           strcmpi(Keyword,'refscan') || ...
           strcmpi(Keyword,'ref_scan') % Kept for compatibility
       
        if ~isinteger(Value)
            printf('amri_eeg_gac(): ''refscan'' must be a integer value');
            return;
        else
            realign.refscan=Value;    
        end
        
    elseif strcmpi(Keyword,'realign.refchannel') || ...
           strcmpi(Keyword,'realign.refchan') || ...
           strcmpi(Keyword,'realign.ref_channel') || ...
           strcmpi(Keyword,'realign.ref_chan') || ...
           strcmpi(Keyword,'align2chan') || ...
           strcmpi(Keyword,'align2channel') || ...
           strcmpi(Keyword,'refchannel') || ...
           strcmpi(Keyword,'refch') || ...
           strcmpi(Keyword,'ref_channel') || ...
           strcmpi(Keyword,'ref_ch')
       
       if ischar(Value)
           if strcmpi(Value,'all')
               realign.refchannel=1:size(eeg.data,1);
           else
               temp_channel=[];
               for ichan=1:size(eeg.data,1)
                   if strcmpi(Value,eeg.chanlocs(ichan).labels)
                       temp_channel=ichan;
                       break;
                   end
               end
               if isempty(temp_channel)
                   printf(['amri_eeg_gac(): ''refchannel'' ''' Value ''' is not found']);
                   printf(['amri_eeg_gac(): the default refchannel # ' ...
                       int2str(realign.refchannel) ' is used']);
               else
                   realign.refchannel=temp_channel;
               end
           end
           
       elseif isnumeric(Value)
           temp_channel=Value;
           temp_channel((temp_channel>size(eeg.data,1)) | (temp_channel<1))=[];
           if isempty(temp_channel)
               printf('amri_eeg_gac(): the specified ''refchannel'' is not found');
               printf(['amri_eeg_gac(): the default refchannel # ' ...
                       int2str(realign.refchannel) ' is used']);
           else
               realign.refchannel=temp_channel;
           end         
       elseif iscell(Value)
           temp_channel=ones(1,length(Value))*nan;
           for iValue=1:length(Value)
               alabel=Value{iValue};
               if strcmpi(alabel,'all')
                   temp_channel=1:size(eeg.data,1);
                   break;
               end
               for ichan=1:size(eeg.data,1)
                   if strcmpi(alabel,eeg.chanlocs(ichan).labels)
                       temp_channel(iValue)=ichan;
                       break;
                   end
               end
           end
           temp_channel(isnan(temp_channel))=[];
           if isempty(temp_channel)
               printf('amri_eeg_gac(): the specified ''refchannel'' is not found');
               printf(['amri_eeg_gac(): the default refchannel # ' ...
                       int2str(realign.refchannel) ' is used']);
           else
               realign.refchannel=temp_channel;
           end
       else
           printf('amri_eeg_gac(): ''refchannel'' should be either a string or integer');
       end

    elseif strcmpi(Keyword,'check') || strcmpi(Keyword,'checking')
        
        checking = Value;
    elseif strcmpi(Keyword,'verbose')
        if Value>0
            flag_verbose=1;
        else
            flag_verbose=0;
        end
    else
        printf(['amri_eeg_gac(): ' Keyword ' is unknown']);
    end   
end


%% ************************************************************************
%                           scanner_trigger   
% *************************************************************************

% find all scanner trigger events by matching eeg.event.type with the 
% trigger name specified through keyword 'trigger.name'
is_scanner_trigger = zeros(length(eeg.event),1);
for ievent=1:length(eeg.event)
    if ismember(eeg.event(ievent).type,scanner_trigger.name)
        is_scanner_trigger(ievent)=1;
    end
end

% retrieve eeg event index of scanner triggers
scanner_trigger.index = find(is_scanner_trigger);
scanner_trigger.count = length(scanner_trigger.index);

% return, if <2 scanner triggers
if scanner_trigger.count<2
    printf('amri_eeg_gac(): less than 2 scanner triggers');
    printf('amri_eeg_gac(): data are either too short or too few slices');
    return;
end

% retrieve latencies of all scanner triggers
scanner_trigger.latency = zeros(length(scanner_trigger.index),1);
for itrigger=1:length(scanner_trigger.index)
    ievent=scanner_trigger.index(itrigger);
    scanner_trigger.latency(itrigger)=eeg.event(ievent).latency;
end

% in case scanner_trigger.latency is not integer, convert it to integer
scanner_trigger.latency=round(scanner_trigger.latency);

% print info about scanner triggers
% if flag_verbose==1
%     printf(['amri_eeg_gac(): ' ...
%         'scanner triggers: type = ' scanner_trigger.type '; ' ...
%         'count = ' int2str(scanner_trigger.count) '; ' ...
%         'trigger interval = ' ...
%         num2str(median(diff(scanner_trigger.latency))/eeg.srate) ' s']);
% end
                     
%% ************************************************************************
%         sort out some information about fmri acquisition
%
% ------------      fmri.nslice, fmri.nvolume     -------------------------
% ------------      fmri.trsec,  fmri.trpnt       -------------------------
% *************************************************************************


% set or compute fmri.nslice, fmri.nvolumem fmri.trsec, fmri.trpnt

if strcmpi(scanner_trigger.type,'slice')
    
    sliceinterval  = round(median(diff(scanner_trigger.latency)));
    
    if      isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  &&  isnan(fmri.trsec) 

        fmri.nslice  = 1;     
        fmri.nvolume = floor(scanner_trigger.count / fmri.nslice);
        fmri.trsec   = fmri.nslice * sliceinterval / eeg.srate;

    elseif  isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
        fmri.nslice  = floor(fmri.trsec * eeg.srate / sliceinterval);     
        fmri.nvolume = floor(scanner_trigger.count / fmri.nslice);

    elseif  isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  &&  isnan(fmri.trsec)
        
        fmri.nslice  = floor(scanner_trigger.count / fmri.nvolume);
        fmri.trsec   = fmri.nslice * sliceinterval / eeg.srate;
        
    elseif  isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
        fmri.nslice  = floor(scanner_trigger.count / fmri.nvolume);
        
    elseif ~isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  &&  isnan(fmri.trsec)
        
        fmri.nvolume = floor(scanner_trigger.count / fmri.nslice);
        fmri.trsec   = fmri.nslice * sliceinterval / eeg.srate;
        
    elseif ~isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
        fmri.nvolume = floor(scanner_trigger.count / fmri.nslice);
        
    elseif ~isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  &&  isnan(fmri.trsec)
        
        fmri.trsec   = fmri.nslice * sliceinterval / eeg.srate;
        
    elseif ~isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
    end
    
elseif strcmpi( scanner_trigger.type,  'volume')
    
    volumeinterval = round(median(diff(scanner_trigger.latency)));
    
    if      isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  &&  isnan(fmri.trsec) 

        fmri.nslice  = 1;     
        fmri.nvolume = scanner_trigger.count;
        fmri.trsec   = volumeinterval / eeg.srate;

    elseif  isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
        fmri.nslice  = 1;     
        fmri.nvolume = scanner_trigger.count;

    elseif  isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  &&  isnan(fmri.trsec)
        
        fmri.nslice  = 1;
        fmri.trsec   = volumeinterval / eeg.srate;
        
    elseif  isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
        fmri.nslice  = 1;
        
    elseif ~isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  &&  isnan(fmri.trsec)
        
        fmri.nvolume = scanner_trigger.count;
        fmri.trsec   = volumeinterval / eeg.srate;
        
    elseif ~isnan(fmri.nslice)  &&  isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
        fmri.nvolume = scanner_trigger.count;
        
    elseif ~isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  &&  isnan(fmri.trsec)
        
        fmri.trsec   = volumeinterval / eeg.srate;
        
    elseif ~isnan(fmri.nslice)  && ~isnan(fmri.nvolume)  && ~isnan(fmri.trsec)
        
    end    
end

fmri.trpnt = fmri.trsec * eeg.srate;

% if flag_verbose>0
%     printf(['amri_eeg_gac(): fMRI params: nslice  = ' ...
%         int2str(fmri.nslice)  '; ' ...
%         'nvolume = ' int2str(fmri.nvolume) '; ' ...
%         'tr = ' num2str(fmri.trsec) ' s']);
% end


%% ************************************************************************
%                           scanner_event     
% *************************************************************************

if strcmpi(scanner_event.type,          'volume') && ...
   strcmpi(scanner_trigger.type,        'volume')
    
    % one event per one volume (i.e. trigger)
    scanner_event.count   = min([scanner_trigger.count ...
                                 fmri.nvolume]);
    scanner_event.latency = scanner_trigger.latency;
    scanner_event.latency = scanner_event.latency(1:scanner_event.count);
    
elseif strcmpi(scanner_event.type,      'volume') && ...
       strcmpi(scanner_trigger.type,    'slice')

    % one event per one volume
    % i.e. N triggers (N: # of slices per each volume)
    scanner_event.count   = min([fmri.nvolume ...
                              floor(scanner_trigger.count/fmri.nslice)]);
    scanner_event.latency = scanner_trigger.latency(1:fmri.nslice:end);
    scanner_event.latency = scanner_event.latency(1:scanner_event.count);
    
elseif strcmpi(scanner_event.type,      'slice') && ...
       strcmpi(scanner_trigger.type,    'slice')
   
    % one event per slice (i.e. trigger)
    scanner_event.count   = min([fmri.nvolume*fmri.nslice ...
                                 scanner_trigger.count]);
    scanner_event.latency = scanner_trigger.latency;
    scanner_event.latency = scanner_event.latency(1:scanner_event.count);
    
elseif strcmpi(scanner_event.type,      'slice') && ...
       strcmpi(scanner_trigger.type,    'volume')
    % todo...
end

% by default, tshift = 0
scanner_event.tshift = zeros(size(scanner_event.latency));

% if flag_verbose>0
%     printf(['amri_eeg_gac(): ' ...
%         'scanner events: type = ' scanner_event.type '; ' ...
%         'count = ' int2str(scanner_event.count) '; ' ...
%         'event interval = ' 
%         num2str(median(diff(scanner_event.latency))/eeg.srate) ' s']);
% end

if checking>0
    figure(1); 
    set(gcf,'name','timing (trigger and event)','color','w');
    if ~isnan(checkfig.position)
        set(gcf,'position',checkfig.position);
    end
    subplot('position',[0.1 0.65 0.8 0.25]);
    plot(2:scanner_trigger.count,diff(scanner_trigger.latency)/eeg.srate);
    xlim([2 scanner_trigger.count]);
    legend({'scanner trigger interval'},'location','Best');
    xlabel('Trigger #'); ylabel('Trigger Interval (sec)');
    box off; grid on;
    subplot('position',[0.1 0.3 0.8 0.25]);
    plot(2:scanner_event.count,diff(scanner_event.latency)/eeg.srate);
    xlim([2 scanner_event.count]);
    legend({'scanner event interval'},'location','Best');
    xlabel('Event #'); ylabel('Event Interval (sec)');    
    box off; grid on;
    reply = input(['continue checking mode [' int2str(checking) ']'],'s');
    if ~isempty(reply)
        checking = round(str2double(reply)); 
    end
    
    checkfig.position=get(gcf,'position');

    close(gcf);
end
                     
                     
%% ************************************************************************
%                       moving window
% *************************************************************************

% determine moving_window.size
if isnan(moving_window.size)
    % if the size of moving window is not specified 
    % assume a moving window that covers all scanner events
    moving_window.size = scanner_event.count;
else
    % assure moving_window.size > min_window_size but < scanner_event.count
    min_window_size = moving_window.minsize;
    if moving_window.size < min_window_size
        moving_window.size = min_window_size;
    end
    if moving_window.size > scanner_event.count
        moving_window.size = scanner_event.count;
    end
    % make moving_window.size an odd number
    if mod(moving_window.size,2)==0 && ...
       moving_window.size<scanner_event.count
        moving_window.size=moving_window.size+1;
    end
end

% set indices to the starting and ending events of each moving window
moving_window.start = ones(scanner_event.count,1) * NaN;
moving_window.end   = ones(scanner_event.count,1) * NaN;
moving_window.index = ones(scanner_event.count,moving_window.size) * NaN;

if moving_window.size == scanner_event.count
    moving_window.start = ones(scanner_event.count,1);
    moving_window.end   = ones(scanner_event.count,1) * scanner_event.count;
    moving_window.index = repmat(1:moving_window.size, scanner_event.count, 1);
else
    % note: moving_window.size is an odd number
    half_window_size    = fix(moving_window.size/2);
    for iscan = 1 : scanner_event.count
        moving_window.start(iscan) = max([iscan-half_window_size 1]);
        moving_window.end  (iscan) = min([iscan+half_window_size scanner_event.count]);
        if moving_window.end(iscan)-moving_window.start(iscan)+1<moving_window.size
            % decrease start index or increase end index
            if     moving_window.start(iscan)==1
                % cannot decrease start index
                moving_window.end(iscan) = min([scanner_event.count ...
                                                moving_window.end(iscan)+...
                                                half_window_size-(iscan-1)]);
            elseif moving_window.end  (iscan)==scanner_event.count
                % cannot increase end index
                moving_window.start(iscan) = max([1 ...
                                                  moving_window.start(iscan)-...
                                                  (half_window_size-scanner_event.count+iscan)]);
            end
        end
        moving_window.index(iscan,:) = moving_window.start(iscan) : moving_window.end(iscan);
    end
end

% if flag_verbose>0
%     printf(['amri_eeg_gac(): moving window of ' ...
%         int2str(moving_window.size) ' ' scanner_event.type 's']);
%     printf(['amri_eeg_gac(): note that a moving window may not be used '...
%             'for some algorithms']);
% end

%% ************************************************************************
%              Filter Setting and design
%
% Using firls for filter design requires signal processing toolbox included
% in the matlab. 
% The default inverse fft method does not require this toolbox. 
% *************************************************************************

% ------------------------
%     lowpass filter
% ------------------------

if isnan(lowpass.cutoff)
    lowpass.cutoff = eeg.srate/downsample_rate/2;
end

if strcmpi(lowpass.method,'firls')
    % design lowpass filter using firls
    minfac = 3;
    min_filtorder = 24;
    fs = eeg.srate;
    F = [0 lowpass.cutoff lowpass.cutoff*(1+lowpass.trans) fs/2];
    F = F * 2 / fs;
    A = [1 1 0 0];
    lowpass.order = max([minfac*fix(fs/lowpass.cutoff) min_filtorder]);
    if mod(lowpass.order,2)==1
        lowpass.order=lowpass.order+1;
    end
    lowpass.filter=firls(lowpass.order,F,A);
elseif strcmpi(lowpass.method,'ifft')
    % design lowpass filter using frequency domain filter

    % sampling frequency
    fs=eeg.srate;
    % number of frequency points
    nfft=2^nextpow2(size(eeg.data,2));
    % frequency vector up to Nygt frequency
    fv = fs/2*linspace(0,1,nfft/2);
    % res in frequency domain
    fres = (fv(end)-fv(1))/(nfft/2-1);
    % frequency-domain filter
    lowpass.filter = ones(1,nfft);
    idxh = round(lowpass.cutoff/fres)+1;
    idxhpt = round(lowpass.cutoff*(1+lowpass.trans)/fres)+1;
    
    idxh = max(min(idxh,nfft),1);
    idxhpt = max(min(idxhpt,nfft),1);
    
    lowpass.filter(idxh:idxhpt)=0.5*(1+sin(pi/2+linspace(0,pi,idxhpt-idxh+1)));
    lowpass.filter(idxhpt:nfft/2)=0;
    lowpass.filter(nfft/2+1:nfft-idxh+1)=lowpass.filter(nfft/2:-1:idxh);
end

% ------------------------
%     highpass filter
% ------------------------

if isnan(highpass.cutoff)
    highpass.cutoff = fmri.nslice/fmri.trsec/2;
end
if strcmpi(highpass.method,'firls')
    % design highpass filter using firls
    minfac=3;
    min_filtorder = 24;
    fs = eeg.srate;
    F = [0 highpass.cutoff*(1-highpass.trans) highpass.cutoff fs/2];
    F = F * 2 / fs;
    A = [0 0 1 1];
    highpass.order = max([minfac*fix(fs/highpass.cutoff) min_filtorder]);
    if mod(highpass.order,2)==1
        highpass.order=highpass.order+1;
    end
    highpass.filter = firls(highpass.order,F,A);
elseif strcmpi(highpass.method,'ifft')
    % design highpass filter for ifft
    % sampling frequency
    fs=eeg.srate;
    % number of frequency points
    nfft=2^nextpow2(size(eeg.data,2));
    % frequency vector up to Nygt frequency
    fv = fs/2*linspace(0,1,nfft/2);
    % res in frequency domain
    fres = (fv(end)-fv(1))/(nfft/2-1);
    % frequency-domain filter
    highpass.filter = ones(1,nfft);
    idxl = round(highpass.cutoff/fres)+1;
    idxlmt = round(highpass.cutoff*(1-highpass.trans)/fres)+1;
    
    idxl = max(min(idxl,nfft),1);
    idxlmt = max(min(idxlmt,nfft),1);
    
    highpass.filter(idxlmt:idxl)=0.5*(1+sin(-pi/2+linspace(0,pi,idxl-idxlmt+1)));
    highpass.filter(1:idxlmt)=0;
    highpass.filter(nfft-idxl+1:nfft)=highpass.filter(idxl:-1:1);
end

% ------------------------
%      bandstop filter
% ------------------------

if isnan(bandstop.tstart)
    bandstop.tstart = scanner_trigger.latency(1);
end
if isnan(bandstop.tend)
    bandstop.tend   = scanner_trigger.latency(end);
end

if strcmpi(bandstop.method,'firls')
    % todo... (seems not working)
    % design a linear phase least-squares FIR filter
    fs=eeg.srate;
    minfac = 3;
    min_filtorder = 15;
    basefreq = fmri.nslice/fmri.trsec;
    bandstop.filter  = cell(round(lowpass.cutoff/basefreq),1);
    bandstop.order   = zeros(length(bandstop.filter),1); 
    bandstop.lcutoff = zeros(length(bandstop.filter),1);
    bandstop.hcutoff = zeros(length(bandstop.filter),1);
    for i = 1 : length(bandstop.filter)
        fcenter=basefreq*i;
        bandstop.lcutoff(i) = fcenter-bandstop.width/2;
        bandstop.hcutoff(i) = fcenter+bandstop.width/2;
        trans  = bandstop.trans * bandstop.width;
        F = [0 bandstop.lcutoff(i)-trans/2 bandstop.lcutoff(i)+trans/2 ...
               bandstop.hcutoff(i)-trans/2 bandstop.hcutoff(i)+trans/2 fs/2];
        F = F * 2 / fs;
        A = [1 1 0 0 1 1];
        bandstop.order(i)=max([minfac*fix(fs/fcenter) min_filtorder]);
        if mod(bandstop.order(i),2)==1
            bandstop.order(i)=bandstop.order(i)+1;
        end
        bandstop.filter{i}=firls(bandstop.order(i),F,A);
    end
elseif strcmpi(bandstop.method,'ifft')
    % slice repetition frequency
    basefreq=fmri.nslice/fmri.trsec;
    % design a frequency-domain filter 
    fs=eeg.srate;
    % number of frequency points
    nfft = 2^nextpow2(bandstop.tend-bandstop.tstart+1);
    % frequency vector up to Nygt frequency
    fv = fs/2*linspace(0,1,nfft/2);
    % res in frequency domain
    fres = (fv(end)-fv(1))/(nfft/2-1);
    % design filter
    bandstop.filter = ones(1,nfft);
    % number of stop bands
    nrofstopbands = fix(fs/2/basefreq);
    bandstop.lcutoff = zeros(nrofstopbands,1);
    bandstop.hcutoff = zeros(nrofstopbands,1);
    % set 0 to stop bands, 1 elsewhere
    for i=1:nrofstopbands
        fcenter=basefreq*i;
        bandstop.lcutoff(i)=fcenter-bandstop.width/2;
        bandstop.hcutoff(i)=fcenter+bandstop.width/2;
        trans=bandstop.width*bandstop.trans;
        idxlmt = round((bandstop.lcutoff(i)-trans/2)/fres)+1;
        idxlpt = round((bandstop.lcutoff(i)+trans/2)/fres)+1;
        idxhmt = round((bandstop.hcutoff(i)-trans/2)/fres)+1;
        idxhpt = round((bandstop.hcutoff(i)+trans/2)/fres)+1;
        % ---------------  bug fix (ZMLIU, 2011-11-18)  ------------------
        % make sure the fourier indices are positive integers 
        idxlmt = max(min(idxlmt,nfft),1);
        idxlpt = max(min(idxlpt,nfft),1);
        idxhmt = max(min(idxhmt,nfft),1);
        idxhpt = max(min(idxhpt,nfft),1);
        % ----------------------------------------------------------------

        if trans==0
            bandstop.filter(idxlmt:idxhpt)=0;
        else
%             bandstop.filter(idxlmt:idxlpt)=linspace(1,0,idxlpt-idxlmt+1);
            bandstop.filter(idxlmt:idxlpt)=0.5*(1+sin(pi/2+linspace(0,pi,idxlpt-idxlmt+1)));
            bandstop.filter(idxlpt:idxhmt)=0;
            bandstop.filter(idxhmt:idxhpt)=0.5*(1+sin(-pi/2+linspace(0,pi,idxhpt-idxhmt+1)));
%             bandstop.filter(idxhmt:idxhpt)=linspace(0,1,idxhpt-idxhmt+1);
        end
        % symmetric change to the "negative" half of the fourier space
        bandstop.filter(nfft-idxhpt+1:nfft-idxlmt+1)=...
            bandstop.filter(idxhpt:-1:idxlmt);
    end
elseif strcmpi(bandstop.method,'butter')
    % todo... (seems not working)
    % design a linear phase least-squares FIR filter
    fs=eeg.srate;
    norder = 3;
    basefreq = fmri.nslice/fmri.trsec;
    bandstop.filter  = cell(round(lowpass.cutoff/basefreq),1);
    bandstop.order   = zeros(length(bandstop.filter),1); 
    bandstop.lcutoff = zeros(length(bandstop.filter),1);
    bandstop.hcutoff = zeros(length(bandstop.filter),1);
    for i = 1 : length(bandstop.filter)
        fcenter=basefreq*i;
        bandstop.lcutoff(i) = fcenter-bandstop.width/2;
        bandstop.hcutoff(i) = fcenter+bandstop.width/2;
        bandstop.order=norder;
        [B,A]= butter(bandstop.order,...
                  [bandstop.lcutoff(i) bandstop.hcutoff(i)]/fs*2,...
                  'stop');
        bandstop.filter{i}=[A;B];
    end    
end

%% ************************************************************************
% remove DC for each channel
% *************************************************************************
for ichan=1:size(eeg.data,1)
     eeg.data(ichan,:)=eeg.data(ichan,:)-mean(eeg.data(ichan,:));
end 

%% ************************************************************************
%          scanner event realignment (scanner_event.tshift)
% 1. choose one or multiple channels for realignment
% 2. choose one scanner event as reference
% 3. for each channel of choice
%     3.1 retrieve the whole time course
%     3.2 highpass filtering
%     3.3 upsampling
%     3.4 segment the time course into many epochs, each of which
%         corresponds to a scanner event. Every epoch has a length of slice
%         tr, and starts from 2 original (before upsampling) data points
%         before the scanner event latency. 
%     3.5 shift the epoch by a small time jitter ranging from -maxtshift to
%         maxtshift
%     3.6 compute the correlation between each epoch with a time jitter and 
%         the epoch corresponding to the reference event
%     3.7 find the best time jitter that gives rise to the highest
%         correlation between the reference epoch and the epoch of interest
% 4. if more than one reference channels, choose the median of the
%    "optimal" jitter obtained from each reference channel.
% 5. save the jitter information into scanner_event.tshift
% *************************************************************************

if realign.flag>=1
%     if flag_verbose>0
%         printf(['amri_eeg_gac(): ' scanner_event.type ' time correction']);
%     end
    max_tshift   = round(realign.maxtshift*upsample_rate);
    epoch_onset  = -2 * max_tshift;
    % set epoch_length according to tr
    if strcmpi(scanner_event.type,'slice')
        epoch_length = fix(fmri.trpnt*upsample_rate/fmri.nslice);
    else
        epoch_length = fix(fmri.trpnt*upsample_rate);
    end
    % time-shifted cross-scan correlations
    realign.cc = nan * ones(length(realign.refchannel), ...
                            scanner_event.count, ...
                            2*max_tshift+1);

    % reform realignment for each specified ref_channel
    for irefchan = 1 : length(realign.refchannel)
        
        ref_chan = realign.refchannel(irefchan);   % reference channel
        ref_scan = realign.refscan;                % reference scan event
        
        ts_orig = eeg.data(ref_chan,:);                     % original data
        if strcmpi(highpass.method,'firls')
            ts_new = filtfilt(highpass.filter,1,ts_orig);
        elseif strcmpi(highpass.method,'ifft')
             nfft = 2^nextpow2(length(ts_orig));
             ts_start=1;
             ts_end=length(ts_orig);
             X = fft(ts_orig,nfft);
             X = X.*highpass.filter;
             ts_new = real(ifft(X,nfft));
             ts_new = ts_new(ts_start:ts_end); % truncate data    
        end
        ts_new = sig_upsample(ts_new,upsample_rate);        % upsample

        % reference scan time series
        rts_start = scanner_event.latency(ref_scan)*upsample_rate+epoch_onset;
        rts_end   = rts_start+epoch_length-1;
        rts       = ts_new(rts_start:rts_end);

        % correlation between reference scan and each scan with tshift
        for iscan = 1 : scanner_event.count
            for tshift=-max_tshift:max_tshift
                ats_start = round(scanner_event.latency(iscan)*upsample_rate+epoch_onset+tshift);
                ats_end   = ats_start+epoch_length-1;
                if ats_start<1 || ats_end>length(ts_new)
                    realign.cc(irefchan,iscan,tshift+max_tshift+1)=0;
                    continue;
                end
                % time shifted epoch
                ats = ts_new(ats_start:ats_end);
                realign.cc(irefchan,iscan,tshift+max_tshift+1)=corr(ats',rts');
            end
        end
    end

    realign.cc(isnan(realign.cc(:,1,1)),:,:)=[];
    [realign.ccmax,optimal_tshift]=max(realign.cc,[],3);
    optimal_tshift=optimal_tshift-(max_tshift+1);
    realign.tshift_x_upsample_rate=median(optimal_tshift,1)';   

    scanner_event.tshift=realign.tshift_x_upsample_rate / upsample_rate;
    clear ts_orig ts_new rts ats;
end

if ~isempty(scanner_event.ignore)
    % for the volumes to be ignored, set tshift to zero
    realign.tshift_x_upsample_rate(scanner_event.ignore)=0; %#ok<STRNU>
    scanner_event.tshift(scanner_event.ignore)=0;%                           
end

%% ************************************************************************
%                       Lowpass Filtering (NOT USED)
% *************************************************************************

% if lowpass.cutoff<eeg.srate/2
%     fprintf(['amri_eeg_gac(): lowpass filtering (cutoff=' ...
%              num2str(lowpass.cutoff) 'Hz; ' ...
%              'order=' int2str(lowpass.order) ') ']);
%     for ichan = 1 : length(eeg.chanlocs)
%         if strcmpi(lowpass.method,'ifft')
%              nfft = 2^nextpow2(size(eeg.data,2));
%              ts = eeg.data(ichan,:);
%              ts_start=1;
%              ts_end=length(ts);
%              X = fft(ts,nfft);
%              X = X.*lowpass.filter;
%              ts = real(ifft(X,nfft));
%              ts = ts(ts_start:ts_end); % truncate data
%              eeg.data(ichan,:)=ts;
%         elseif strcmpi(lowpass.method,'firls')
%             eeg.data(ichan,:)=filtfilt(lowpass.filter,1,eeg.data(ichan,:));
%         end
%         fprintf([int2str(ichan) ' ']);
%     end
%     fprintf('\n');
% end


%% ************************************************************************
%                   remove gradient artifact template
% *************************************************************************
msg = '';

for imethod = 1 : length(gac_method)
    
    if flag_verbose>0
        fprintf(['amri_eeg_gac(): run ' gac_method{imethod} ': ']);
    end
    
    if strcmpi(gac_method{imethod},'AAS') || ...
       strcmpi(gac_method{imethod},'AAR') || ...
       strcmpi(gac_method{imethod},'MAS') || ...
       strcmpi(gac_method{imethod},'MAR') || ...
       strcmpi(gac_method{imethod},'PCA') || ...
       strcmpi(gac_method{imethod},'TSE') || ...
       strcmpi(gac_method{imethod},'AFIR')
       
        % set the epoching onset wrt scanner event latency
        epochs.onset  = scanner_event.epochonset * upsample_rate;
        % set the epoch length equals the slice (or volume) tr
        if strcmpi(scanner_event.type,'slice')
            epochs.length = round(fmri.trpnt*upsample_rate/fmri.nslice);
        else
            epochs.length = round(fmri.trpnt*upsample_rate);
        end

        % run artifact correct channel by channel
        for ichan = 1 : length(eeg.chanlocs)
            
            % upsample the signal if upsample_rate~=1
            ts_orig = eeg.data(ichan,:);                   % original data
            ts_new  = sig_upsample(ts_orig,upsample_rate); % upsampled data
            
            % highpass filter data
            if strcmpi(highpass.method,'firls')
                ts_art  = filtfilt(highpass.filter,1,ts_new);
            elseif strcmpi(highpass.method,'ifft')
                nfft = 2^nextpow2(length(ts_new));
                ts_start=1;
                ts_end=length(ts_new);
                X = fft(ts_new,nfft);
                X = X.*highpass.filter;
                ts_art = real(ifft(X,nfft));
                ts_art = ts_art(ts_start:ts_end); % truncate data 
            end
            
            % extract epochs, each corresponds to one scanner event
            epochs.data   = ones(scanner_event.count,epochs.length)*nan;    
            epochs.tstart = ones(scanner_event.count,1)*nan;               
            epochs.tend   = ones(scanner_event.count,1)*nan;               
            
            for iscan = 1 : scanner_event.count
            
                % use the previous calculated jitter
                tshift=round(scanner_event.tshift(iscan)*upsample_rate);
                
                % set the starting time point of the epoch (latency+jitter)
                iepochstart=scanner_event.latency(iscan)*upsample_rate+...
                    tshift+epochs.onset;
                
                % find the ending time point of the epoch by the staring
                % time point and the epoch length
                iepochend  =iepochstart+epochs.length-1;
                
                % extract the epochs from ts_art
                if iepochstart>=1 && iepochend<=length(ts_art)
                    
                    epochs.data(iscan,:)=ts_art(iepochstart:iepochend);
                
                    % detrend if specified
                    if ~isnan(detrend_option)
                        if detrend_option==0
                            epochs.data(iscan,:)=detrend(epochs.data(iscan,:),'constant');
                        elseif detrend_option==1
                            epochs.data(iscan,:)=detrend(epochs.data(iscan,:),'linear');
                        end
                    end
                    % save the starting and ending time points of the epoch
                    epochs.tstart(iscan)=iepochstart;
                    epochs.tend(iscan)=iepochend;                    
                end
                
            end

            % for incomplete epochs, which likely exist in the beginning
            % and end of the data, set them identical to the closest
            % epochs, instead of remove them. by doing so, we ensure the
            % number of epochs equals the number of scanner events
            % nevetheless, the incomplete epochs hardly happen
            
%             if isnan(epochs.data(1,1))    % if the first data point in the first epoch is nan
%                 epochs.data(1,:)=epochs.data(2,:); 
%             elseif isnan(epochs.data(end,end)) % if the last data point in the last epoch is nan
%                 epochs.data(end,:)=epochs.data(end-1,:);
%             end

%begin fix for crash when final gradient artifact is too close to the end
%of the EEG recording.  JD
             if isnan(epochs.data(1,1))    % if the first data point in the first epoch is nan
                 epochs.data(1,:)=[]; 
                 epochs.tstart(1)=[];
                 epochs.tend(1)=[];
             elseif isnan(epochs.data(end,end)) % if the last data point in the last epoch is nan
                 epochs.data(end,:)=[]; 
                 epochs.tstart(end)=[];
                 epochs.tend(end)=[];
             end

%end fix for crash when final gradient artifact is too close to the end
%of the EEG recording.  JD

% count the total number of epochs
            epochs.count = size(epochs.data,1);
            
            % compute the median (across epochs), which is robust against outliers
            epochs.median = median(epochs.data,1);
            % compute the mean
            epochs.mean   = mean(epochs.data,1);
            
            % compute correlation with the median
            epochs.cc=zeros(epochs.count,1);
            for iepoch=1:epochs.count
                epochs.cc(iepoch)=amri_sig_corr(epochs.data(iepoch,:)',epochs.median');
            end
            
            % find outliers with small cross-epoch correlation
            tt=sort(epochs.cc,'ascend'); 
            tt_1q=tt(round(length(tt)/4));
            epochs.outlier = find(epochs.cc<tt_1q-4*amri_stat_iqr(epochs.cc));
            
            % find epochs.ignore from scanner_event.ignore
            epochs.ignore = scanner_event.ignore;
            epochs.ignore(epochs.ignore<1) = [];
            epochs.ignore(epochs.ignore>epochs.count) = [];
            
            % checking point (debug only)
            if checking>0

                figure(ichan);
                set(gcf,...
                    'name',eeg.chanlocs(ichan).labels,...
                    'color','w');
                if ~isnan(checkfig.position)
                    set(gcf,'position',checkfig.position);
                end
                % draw an image of epochs.data
                subplot('position',[0.1 0.5 0.4 0.4]);
                imagesc(epochs.data); axis xy;
                xlabel('Time'); ylabel('Epoch');
                % plot mean, median and stdev across epochs
                subplot('position',[0.1 0.1 0.4 0.3]);
                plot(epochs.mean,'c');
                hold on;
                plot(epochs.median,'r');
                plot(std(epochs.data),'g');
                tmp_index = 1 : epochs.count;
                tmp_index(ismember(tmp_index,epochs.outlier)|....
                          ismember(tmp_index,epochs.ignore))=[];
                plot(std(epochs.data(tmp_index,:)),'b');
                xlim([1 epochs.length]);
                box off; grid on;
                xlabel('Time');
                ylabel('Ampl (\muV)');
                % plot cc with the median
                subplot('position',[0.6 0.5 0.3 0.4]);
                plot(epochs.cc,1:epochs.count);
                hold on;
                if ~isempty(epochs.outlier)
                    plot(epochs.cc(epochs.outlier),epochs.outlier,'or');
                end
                if ~isempty(epochs.ignore)
                    plot(epochs.cc(epochs.ignore), epochs.ignore, '*r');
                end
                ylim([1 epochs.count]);
                xlimval = get(gca,'xlim');
                xlim([xlimval(1) 1]);
                grid on; box off;
                xlabel('Corr Coeff');
                ylabel('Epoch');
                reply = input(['continue checking mode [' int2str(checking) ']'],'s');
                if ~isempty(reply)
                    checking = round(str2double(reply)); 
                end
                checkfig.position=get(gcf,'position');
                close(gcf);

            end
            
            % artifact correction
       
            if strcmpi(gac_method{imethod},'AAS') || ...
               strcmpi(gac_method{imethod},'AAR') || ...
               strcmpi(gac_method{imethod},'MAS') || ...
               strcmpi(gac_method{imethod},'MAR') || ...
               strcmpi(gac_method{imethod},'TSE') || ...
               strcmpi(gac_method{imethod},'AFIR') 
           
                % moving windows
                for iepoch=1:epochs.count
                    
                    % retrieve pre-calculated moving window
                    iepoch_win_index = moving_window.index(iepoch,:);
                    
                    % exclude outliers
                    iepoch_win_index = iepoch_win_index(~ismember(iepoch_win_index,epochs.outlier));

                    % exclude ignored volumes
                    iepoch_win_index = iepoch_win_index(~ismember(iepoch_win_index,epochs.ignore));
                    
                    % exclude itself 
                    iepoch_win_index = iepoch_win_index(iepoch_win_index~=iepoch);
                    
                    % exclude epochs less than 0.3 sec apart
                    minisi = round(moving_window.nowin*eeg.srate*upsample_rate/epochs.length);
                    iepoch_win_index = iepoch_win_index(abs(iepoch_win_index-iepoch)>minisi);
                    
                    % don't exceed limit
                    iepoch_win_index = iepoch_win_index(iepoch_win_index<=epochs.count);
                    iepoch_win_index = iepoch_win_index(iepoch_win_index>=1);
                    
                    % do nothing if no epoch is available 
                    if isempty(iepoch_win_index)
                        continue;
                    end
                    
                    if strcmpi(gac_method{imethod},'AAS')
                        
                        % subtract local moving average
                        maga=mean(epochs.data(iepoch_win_index,:),1);
                        
                        % detrend if specified
                        if ~isnan(detrend_option)
                            if detrend_option==0
                                maga=detrend(maga,'constant');
                            elseif detrend_option==1
                                maga=detrend(maga,'linear');
                            end
                        end
                        
                        % subtract maga from ts_new 
                        if ~isnan(epochs.tstart(iepoch)) && ...
                           ~isnan(epochs.tend(iepoch))
                            ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))=...
                                ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))-...
                                maga;
                        end
                        
                    elseif strcmpi(gac_method{imethod},'AAR')
                        
                        % regress out local moving average
                        maga=mean(epochs.data(iepoch_win_index,:),1);
                        
                        % detrend if specified
                        if ~isnan(detrend_option)
                            if detrend_option==0
                                maga=detrend(maga,'constant');                        
                            elseif detrend_option==1
                                maga=detrend(maga,'constant');
                            end
                        end
                        
                        % regress out maga from ts_new
                        kk = epochs.data(iepoch,:)'\maga';    % regression coefficient
                        
                        if ~isnan(epochs.tstart(iepoch)) && ...
                           ~isnan(epochs.tend(iepoch))                        
                            ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))=...
                                ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))-...
                                kk*maga;
                        end
                        
                    elseif strcmpi(gac_method{imethod},'MAS')
                        
                        % subtract local moving median
                        mmga=median(epochs.data(iepoch_win_index,:),1);
                        % detrend if specified
                        if ~isnan(detrend_option)
                            if detrend_option==0
                                mmga=detrend(mmga,'constant');                        
                            elseif detrend_option==1
                                mmga=detrend(mmga,'linear');
                            end
                        end
                        % subtract mmga from ts_new
                        if ~isnan(epochs.tstart(iepoch)) && ...
                           ~isnan(epochs.tend(iepoch))                          
                            ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))=...
                                ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))-...
                                mmga;
                        end
                        
                    elseif strcmpi(gac_method{imethod},'MAR')
                        
                        % regress out local moving median
                        mmga=median(epochs.data(iepoch_win_index,:),1);
                        % detrend if specified
                        if ~isnan(detrend_option)
                            if detrend_option==0
                                mmga=detrend(mmga,'constant');                        
                            elseif detrend_option==1
                                mmga=detrend(mmga,'linear');
                            end
                        end
                        % regress out mmga from ts_new
                        kk = epochs.data(iepoch,:)'\mmga';
                        if ~isnan(epochs.tstart(iepoch)) && ...
                           ~isnan(epochs.tend(iepoch))                          
                            ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))=...
                                ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))-...
                                kk*mmga;  
                        end
                        
                    elseif strcmpi(gac_method{imethod},'TSE')  || ...
                           strcmpi(gac_method{imethod},'AFIR') 
                       
                        % ts0 is the mean across epochs
                        ts0 = mean(epochs.data(iepoch_win_index,:),1);
                    
                        % derive 1 -> n-order derivatives of 
                        tsn = diffn(ts0,1:ts.order);
                        
                        % normalize tsn
                        for iorder = 1 : ts.order
                            tsn(iorder,:) = tsn(iorder,:) / norm(tsn(iorder,:));
                        end
                        
                        % first subtract ts0
                        ts_tmp=epochs.data(iepoch,:)-ts0;

                        % then regress out tsn
                        coeffn = (pinv(tsn')*ts_tmp')';
                        if ~isnan(epochs.tstart(iepoch)) && ...
                           ~isnan(epochs.tend(iepoch))                           
                            ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))= ...
                                ts_new(epochs.tstart(iepoch):epochs.tend(iepoch)) - ...
                                coeffn * tsn - ts0;
                        end
                    else
                        

                    end  
                    
                end

            elseif strcmpi(gac_method{imethod},'PCA') 
                
                % use all epochs despite the moving window setting
                epochindex=(1:epochs.count)';
                % exclude the ignored epochs but not outliers
                % 1) assure continuity of epochs
                % 2) svd is robust against occasional outliers
                
                epochindex=epochindex(~ismember(epochindex,epochs.ignore));
                
                if strcmpi(pca.selection,'auto')
                    % 07/31/2013, ZMLIU 
                    % Actually, svd is not SO robust against outliers
                    % it would improve the accurary of artifact templates if
                    % excluding the outlier epochs from the matrix 
                    epochindex_no_outlier= ...
                        epochindex(~ismember(epochindex,epochs.outlier));
                    % retrieve a data matrix
                    Z = epochs.data(epochindex_no_outlier,:);
                    
                    % run compact SVD
                    % U(:,i) - time series within an epoch
                    % S(i,i) - scaling factor 
                    % V(:,i) - variation across epochs
                    [U,S,V]=svd(Z','econ'); %#ok<ASGLU>
                    numcomp=min(size(Z));
                    
                    % run a number of statistical test on each colume of V
                    isgac  =   false(numcomp,1);
                    
                    h_ttest     =   false(numcomp,1); 
                    p_ttest     =   ones(numcomp,1);
                    s_acorr     =   1; % skip in sec
                    max_acorr   =   zeros(numcomp,1);
                    
%                     alpha = pca.alpha/numcomp;
                    alpha = pca.alpha;
                    
                    for icomp = 1 : numcomp
                        av = V(:,icomp);
                        
                        % ttest
                        [h_ttest(icomp),p_ttest(icomp)] = amri_stat_ttest(av,0,alpha);
                       
                        % autocorrelation
                        acf = amri_sig_xcorr(av,av);
                        
                        acf_skip = false(length(acf),1);
                        acf_center = fix(length(acf)/2)+1;
                        pt_skip = round(s_acorr*eeg.srate*upsample_rate/epochs.length);
                        acf_skip(acf_center-pt_skip:acf_center+pt_skip)=1;
                        
                        max_acorr(icomp)=abs(max(acf(~acf_skip)));
                        
                        % combine test outcomes
                        isgac(icomp) = h_ttest(icomp) | max_acorr(icomp)>0.9;
                    end
                    
                    if any(isgac)

                        % construct optimal basis function (obsfunc)
                        obsindex = find(isgac==1);

                        obsfunc  = zeros(epochs.length,length(obsindex));
                        for ii =  1 : length(obsindex)
                            obsfunc(:,ii)=U(1:epochs.length,obsindex(ii));
                            obsfunc(:,ii)=obsfunc(:,ii) - mean(obsfunc(:,ii));
                            obsfunc(:,ii)=obsfunc(:,ii)/norm(obsfunc(:,ii));
                        end

                        % regress out the basis functions
                        for iepoch=1:epochs.count
                            % use pinv to compute the regression coefficient of every basis function
                            % 07/31/2013
                            if ismember(iepoch,epochs.outlier)
                                maxvar=0;
                                for ishf=-10:10
                                    obsfunc_shf=circshift(obsfunc,[ishf 0]);
                                    obswght=pinv(obsfunc_shf)*epochs.data(iepoch,:)';
                                    afit=(obsfunc_shf*obswght)';
                                    avar=var(afit);
                                    if avar>maxvar
                                        maxvar=avar;
                                        gg=ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))-afit;
                                    end
                                end
                                ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))=gg;
                            else
                                obswght=pinv(obsfunc)*epochs.data(iepoch,:)';
                                % regress out the signal explained by the basis functions
                                if ~isnan(epochs.tstart(iepoch)) && ...
                                   ~isnan(epochs.tend(iepoch))                         
                                    ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))=...
                                        ts_new(epochs.tstart(iepoch):epochs.tend(iepoch))-...
                                        (obsfunc*obswght)';
                                end
                            end
                        end
                    end
                end
            else                    
                    
            end
            
            % resume original sampling rate
            ts_new = ts_new(1*upsample_rate:upsample_rate:length(ts_new));
            
            % update eeg
            eeg.data(ichan,:) = ts_new;
            if flag_verbose>0
                if ~isempty(msg)
                    fprintf(repmat('\b',1,length(msg)));
                end
            end
            
            msg = [int2str(ichan) '/' int2str(size(eeg.data,1))];
            if flag_verbose>0
                fprintf(msg);
            end
        end
        if flag_verbose>0
            fprintf('\n');
        end
        
    elseif strcmpi(gac_method{imethod},'blf')
        
        for ichan=1:length(eeg.chanlocs)
            
            if strcmpi(bandstop.method,'ifft')
                
                % bandstop filter using fft and ifft
                nfft = 2^nextpow2(bandstop.tend-bandstop.tstart+1);
                ts = eeg.data(ichan,bandstop.tstart:bandstop.tend);
                ts_start=1;
                ts_end=length(ts);
                % zero padding
                ldiff = nfft-length(ts);
                if ldiff==0
                    ts_start=1;
                    ts_end=nfft;
                elseif mod(ldiff,2)==1 && ldiff>0
                    ts_start=1+fix(ldiff/2)+1;
                    ts_end  =nfft-fix(ldiff/2);
                elseif mod(ldiff,2)==0 && ldiff>0
                    ts_start=1+ldiff/2;
                    ts_end  =nfft-ldiff/2;
                end
                ts = [zeros(1,ts_start-1) ts zeros(1,nfft-ts_end)]; %#ok<AGROW>
                X = fft(ts,nfft);
                X = X.*bandstop.filter;
                ts = real(ifft(X,nfft));
                ts = ts(ts_start:ts_end); % truncate data
                eeg.data(ichan,bandstop.tstart:bandstop.tend)=ts;
                
            elseif strcmpi(bandstop.method,'butter')
                % bandstop filter using butterworth filter
                for i=1:length(bandstop.filter)
                    AB=bandstop.filter{i};
                    eeg.data(ichan,bandstop.tstart:bandstop.tend)=...
                        filtfilt(AB(1,:),AB(2,:),...
                        eeg.data(ichan,bandstop.tstart:bandstop.tend));
                end
            elseif strcmpi(bandstop.method,'firls')
                % bandstop filter using least squares linera phase filter
                for i=1:length(bandstop.filter)
                    eeg.data(ichan,bandstop.tstart:bandstop.tend)=...
                        filtfilt(bandstop.filter{i},1,...
                        eeg.data(ichan,bandstop.tstart:bandstop.tend));
                end
            end
            if flag_verbose>0
                if ~isempty(msg)
                    fprintf(repmat('\b',1,length(msg)));
                end
            end
            msg = [int2str(ichan) '/' int2str(size(eeg.data,1))];
            if flag_verbose>0
                fprintf(msg);
            end
        end
        if flag_verbose>0
            fprintf('\n');
        end
    end
end

%% ************************************************************************
%                       Lowpass Filtering
% *************************************************************************

if lowpass.cutoff<eeg.srate/2
    msg='';
%     if flag_verbose>0
%         fprintf(['amri_eeg_gac(): lowpass filtering (cutoff=' ...
%                 num2str(lowpass.cutoff) 'Hz; ' ...
%                 'order=' int2str(lowpass.order) ') ']);
%     end
    for ichan = 1 : length(eeg.chanlocs)
        if strcmpi(lowpass.method,'ifft') 
             nfft = 2^nextpow2(size(eeg.data,2));
             ts = eeg.data(ichan,:);
             ts_start=1;
             ts_end=length(ts);
             X = fft(ts,nfft);
             X = X.*lowpass.filter;
             ts = real(ifft(X,nfft));
             ts = ts(ts_start:ts_end); % truncate data
             eeg.data(ichan,:)=ts;
        elseif strcmpi(lowpass.method,'firls')
            eeg.data(ichan,:)=filtfilt(lowpass.filter,1,eeg.data(ichan,:));
        end
       
%         if flag_verbose>0
%             if ~isempty(msg)
%                 fprintf(repmat('\b',1,length(msg)));
%             end
%         end
        msg = [int2str(ichan) '/' int2str(size(eeg.data,1))];
%         if flag_verbose>0
%             fprintf(msg);
%         end
    end
%     if flag_verbose>0
%         fprintf('\n');
%     end
end


%% ************************************************************************
%                           Downsample
% *************************************************************************

if downsample_rate>1
    if flag_verbose>0
        printf(['amri_eeg_gac(): downsample by ' int2str(downsample_rate)]);
    end
    eeg.data  = eeg.data(:,1:downsample_rate:size(eeg.data,2));
    eeg.srate = eeg.srate/downsample_rate;
    eeg.pnts  = floor(eeg.pnts/downsample_rate);
    if ~isfield(eeg,'xmin'),eeg.xmin=0; end;
    eeg.xmax = eeg.xmin + (eeg.pnts-1)/eeg.srate;
%     anc.refsig = anc.refsig(1:downsample_rate:length(anc.refsig));
    for ievent=1:length(eeg.event)
        eeg.event(ievent).latency=round(eeg.event(ievent).latency/downsample_rate);
        eeg.event(ievent).latency=max([1 eeg.event(ievent).latency]);
    end
    for ievent=1:length(eeg.urevent)
        eeg.urevent(ievent).latency=round(eeg.urevent(ievent).latency/downsample_rate);
        eeg.urevent(ievent).latency=max([1 eeg.urevent(ievent).latency]);
    end
end


%% ************************************************************************
% Linear Detrend
% *************************************************************************
for ichan=1:size(eeg.data,1)
    eeg.data(ichan,:)=detrend(eeg.data(ichan,:),'linear');
end

%%
outeeg = eeg;
return;

%% ************************************************************************
%                       upsampling 
% *************************************************************************

function ts_new = sig_upsample(ts_orig, upsample_rate)

if upsample_rate>1
    t_orig = 1:length(ts_orig);                                 % orig (coarse) time points
    t_new  = 1/upsample_rate:1/upsample_rate:length(ts_orig);   % new (fine) time points
    ts_new = spline(t_orig,ts_orig,t_new);
else
    ts_new = ts_orig;
end

return;



%% ************************************************************************
%                       calculate derivative model 
% *************************************************************************

function diffmodel = diffn(ts, orders)

ts = ts(:);
ls = length(ts);
diffmodel = zeros(length(ts),length(orders));

for i = 1 : length(orders)
    order = orders(i);
    si    = diff(ts, order);
    ns    = length(si);
    ms    = fix((ls-ns)/2);
    if mod(order,2)==1
        diffmodel(:,i) = ([zeros(ms,1);si;zeros(ms+1,1)]+[zeros(ms+1,1);si;zeros(ms,1)])/2;
    else
        diffmodel(:,i) = [zeros(ms,1);si;zeros(ms,1)];
    end     
end

diffmodel = diffmodel';

return;






