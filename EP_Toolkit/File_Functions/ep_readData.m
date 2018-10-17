function [EPdata origEloc]=ep_readData(varargin);
%  [EPdata origEloc]=ep_readData(varargin);
%       imports data from a variety of file formats relying on FieldTrip I/O code and converts into EP file format.
%
%Inputs:
%  OPTIONAL KEYWORDS (accepts either a single cell string or a series of variables)
%  Information in keywords override contents of the file.
%  'file'       : Name of the file to be imported.  It needs to have an
%                 appropriate suffix in order for it to be recognized properly.
%                 See the tutorial file.
%                 If none is specified it will assume that it is an EGI EGIA average file.
%
%  'format'     : File format of the file.
%  'samples'    : Samples included in analysis.
%  'cells'      : Cells included in analysis (trials for single_trial data).  Refers to order of appearance in the data file.
%                 Should be specified in order desired for output.  See example.
%  'subjects'   : Subjects included in analysis.  Should be specified in
%                 order desired for output.
%  'channels'   : Channels included in analysis.
%  'montage'    : The electrode montage if an EGIS file.  If not specified and the SLAY resource is not in the file, it
%                 will be requested.  The list of montages can be found in askForMontage.m.
%  'name'       :  The name of the experiment.
%  'labels'      : Renames the condition names subject names (e.g., 'labels',{{'cell1','cell2'},{'sub1','sub2'},{8,10}}).
%                  If only one cell label is specified, it will be applied to all the data.
%                  Frequency labels should be a number that indicates the Hz or an empty set if not frequency data.
%  'type'       : The type of the data.  Should be followed by the data type.
%                 'continuous' raw data from a single subject, not segmented into epochs.
%                 'single_trial' segmented single subject not averaged.
%                 'average' file with trials averaged so there is only one waveform per condition.
%                 'grand_average' file with multiple subjects averaged together.
%                 'factors' output of PCA results, where "subjects" are factors.
%  'NumberSubjects'   : The number of subjects or factors.
%  'prestim'   : The number of msec prior to the stimulus onset (prior to using the samples keyword if used).
%  'ced'      : The name of the .ced file for electrode coordinates.  Use 'none' if none.  Can be invalid ced file.
%  'eloc'      : The eloc information (only if good info).
%  'silent'   : If 'on' rather than 'off', do not provide output to command line for standard messages. (default off)
%  'origReference'  : Specify the original reference channel if present explicitly in the dataset.
%                 Followed by a field for the original reference(s), no more than two, in an array.
%  'currReference'  : Specify the current reference channel if present explicitly in the dataset.
%                 Followed by a field for the current reference(s), no more than two, in an array.
%                 The field can also be the key word 'AVG' for average reference and 'CSD' for current source
%                 density.
%  'textPrefs'  : Specifies the parameters for text files.  Structured variable:
%                 .firstRow : the first row of the data.
%                 .lastRow  : the last row or the data.  Zero equals last row.
%                 .firstCol : the first column of the data.
%                 .lastCol  : the last column of the data.  Zero equals last column.
%                 .orientation : the orientation of the data.  1=channels are columns, 2=channels are rows.
%
%  'elecPrefs'  : Specifies whether electrode coordinates should be rotated when reading from inside files like MFF and FIFF
%                 1=no, 2=90 degree clockwise, 3=180 degrees, 4=270 degrees.
%  'SMIsuffix'  : Specifies add any SMI files with the same file name stem and this suffix to the data.
%  'specSuffix'  : Specifies add any spec text files with the same file name stem and this suffix to the data.
%                  The first line should contain the spec names.
%  'subjectSpecSuffix'  : Specifies add any subject spec text files with the same file name stem and this suffix to the data.
%                  The first line should contain the subject spec names.
%
%  Reads in FieldTrip I/O outputs: data, hdr, eventHdr.  Assumes 'average' data is organized
%  with the subjects grouped together by cell (e.g., C1S1 C1S2 C1S3 C2S1 C2S2 C2S3).
%  Except for simple binary average files where it is assumed to be the other way around.
%  Recommended but not necessarily required is a .ced file (EEGlab format) for electrode coordinates.
%
%Outputs:
%  EPdata         : Structured array with the data and accompanying information.
%    .data      : 7D data matrix [channels, time points, cells/trials, subjects, factors, freqs, relations].
%                 Trials are grouped by cell in single_trial format.
%                 NaN is allowed in non-EEG channels only.
%    .noise     : 4D matrix [channels, time points, cells/trials, subjects] mirroring .data.  Can be empty.
%                 This is the +/- reference average (Schimmel, 1967) which provides an estimate of the noise level
%                 in averaged data by flipping every other trial.  Not applicable to spectral data.  Not from combining data other than across subjects.
%                 For grand average data, every other average is flipped.
%    .std       : 6D matrix [channels, time points, cells/trials, subjects, 1, freqs].  Has full dimensionality even for factor data.  Can be empty.
%                 This is the standard deviation from averaging trials.  Not from combining data other than across subjects.
%    .stdCM     : 6D matrix [channels, time points, cells/trials, subjects, 1, freqs].  Has full dimensionality even for factor data.  Can be empty.
%                 This is the Cousineau-Morey std (std of data for a given time point mean-corrected across the cells of a given subject).  Not from combining data other than across subjects.
%                 This is empty for all but grand average data.
%    .cov        (empty with no subfields if no info)
%      .covMatrix   : Variance-covariance matrix for channels, generated during averaging as indicator of noise levels.
%             Provides original covariance matrix.  Average referenced and not updated when rereferencing.  (subject,chan,chan)
%             NaN for unknown or invalid entries.  This is provided for use by the MNE software, which uses it to whiten
%             the data during source analyses (p.123 of 2.7.3 MNE manual).  Per p.89, epochs are first baseline corrected and then
%             the matrix is actually calculated as an SSCP (sums-of-squares-cross-products) matrix.
%       .Nq     : Number of observations going into the covariance matrix computation, for weighting when combining them. (subject)
%
%    .montage   : String with the montage information, if available.
%    .fileFormat: The original file format.
%    .chanNames : The channel names.
%    .timeNames : The msec of the sample onset with respect to the stimulus onset time.
%                 For time-frequency data, the msec of the middle of the .5 sec window.
%                 For flexible segments, the percentage is still left-sided (so for 5% increments, from zero to 95%).
%    .timeUnits : 'ms' for milliseconds, 'per' for flexible segments.
%    .subNames  : The subject names
%    .cellNames : The cell names (once for each trial for single_trial files).
%    .trialNames: The trial number ID per cell (starting from 1). (single_trial data only)
%                 These should be unique numbers for each cell.
%    .facNames  : The factor names (only for factor files)
%    .freqNames : The frequency at the middle of the frequency bin (numeric array).
%    .relNames  : The channel names for relational data (e.g., coherence).  Mirrors the channels dimension exactly.
%    .chanTypes : The type of the channel: EEG, MGM (magnetometer MEG), MGA (axial MEG), MGP (planar MEG), ECG, ANS (autonomic), REG (regional average)
%                 PPL (pupil dilation), XEY (x-coordinate of eye-tracking), YEY (y-coordinate of eye-tracking)
%                 The CED file can also specify FID (fiducial) and BAD (delete) channel types,
%                 but they will not end up as a chanType.
%    .subTypes  : The type of the subject: RAW (single trial), AVG (subject average), GAV (grand average)
%    .cellTypes : The type of the cell: SGL (one cell), CMB (combination of cells), STS (sample test statistics output)
%    .facTypes  : The type of the factor: SGL (one factor), CMB (combination of factors)
%    .EPver     : The EP Toolkit version information.
%    .ver       : The Matlab version information.
%    .date      : The date the file was created.
%    .Fs        : The sampling frequency in Hz at the current data resolution (not necessarily the original sampling rate).  Equals number of %age increments times 10 for flexible segments (e.g., "200" if 20 bins).
%    .baseline  : The number of samples prior to the trigger event (positive number).
%                 Thus, the ms in the timeNames refers to the offset of each sample.  For 250 Hz, a baseline of 50 samples means
%                 that the epoch started 200 ms (50 samples) prior to the stimulus onset.  A baseline of 1 sample means
%                 that the epoch started 4 ms (1 sample) prior to the stimulus onset.
%                 When baseline falls within the sample, shift to left side of sample.
%    .dataName  : A descriptive name for the dataset, used to differentiate the active datasets during analysis.
%    .ename     : The name of the experiment.
%    .dataType  : The type of the data: 'continuous', 'single_trial', or 'average' (default: average)
%    .trialSpecs       : Cell array (trial,spec,subject) of specific information for trials or their averages.  Can be empty if no specs.
%    .trialSpecNames   : Cell array of the name of each trial spec type.  Can be empty if no specs.
%    .subjectSpecs     : Cell array (subject,spec) of specific information for subjects.  Can be empty if no specs.
%    .subjectSpecNames : Cell array of the name of each subject spec type.  Can be empty if no specs.
%    .events    : Cell array of event structured variables (subject,cell/trial).  Events cells may be blank and thus have no fields.
%      .type      = string (events of "trial" type dropped as they are redundant with cell name information.)
%                   (if there was a recording stop, start of recording is marked by an event with the value 'boundary').
%      .sample    = expressed in samples, the first sample of an epoch is 1
%      .value     = number or string (e.g., "stm+")
%                   (if there was a recording stop, start of recording is marked by an event with the value 'boundary').
%      .duration  = expressed in samples.  For boundary events, represents the length of the recording break.
%      .keys     = additional information (subfields are required).  All field contents are strings even if the datatype is numeric.
%           .code (name)
%           .data (data)
%           .datatype (type of variable, e.g. 'short')
%           .description (notes)
%    .avgNum    : Number of waveforms going into averages (subject,cell)
%                 0 means unknown and -1 means bad.
%    .subNum    : Number of subjects going into averages (subject,cell)
%                 0 means unknown and -1 means bad.
%    .covNum    : Adjusted number of waveforms going into averages (subject,cell)
%                 0 means unknown and -1 means bad.
%                 For use with MNE software.  When computing linear combinations of waves, the effective sample size
%                 is modified to reflect the increase in noise levels, per p. 128 of the 2.7.3 MNE manual.
%    .fileName  : Name of original file.
%    .history   : Command used to create current file.
%    .ced       : The name of the .ced file for electrode coordinates.
%    .eloc      : The electrode location information, one for each channel (see readlocs header).  By apparent EEGlab convention, a row vector.
%                 eloc is the same length as the channels, with non-EEG channels having a blank entry.
%                 For mff files, if internal eloc are present, will be
%                 used.
%    .implicit  : The electrode information for fiducial locations (see readlocs header)
%    .facVecT   : For temporal PCA factor files, the factor waveform.  Used to compress the data. (rows=points,cols=factors)
%    .facVecS   : For spatial PCA factor files, the factor scalp topography.  Used to compress the data.(rows=chans,cols=factors)
%    .facVecF   : For frequency PCA factor files, the factor frequency spectrum.  Used to compress the data.(rows=frequencies,cols=factors)
%    .facData   : 7D data matrix [channels, time points, cells/trials, subjects, factors, freqs, relations] for combined factors (only when facVecS & facVecT & facVecF are used).
%    .facVar    : The variance accounted for by factors.  Unlike the contents of .pca, changed to reflect editing.
%    .facVarQ   : The unique variance accounted for by factors.  Unlike the contents of .pca, changed to reflect editing.
%                 They are represented separately from .data because they can't be compacted using the facVec mechanism.
%    .reference
%        .original    : Original recording reference(s): 1-2 numbers in array
%        .current     : Current reference channel(s): 1-2 numbers in array
%        .type        : Current reference type: REG (regular), AVG (average reference, even if no longer adds to zero),
%                       CSD (current source density), PAR (PARE-corrected average reference)
%    .analysis (for continuous files, divided into one second epochs and "trials" refer to these epochs. Excess time
%                points are tacked onto final epoch)
%       .blinkTrial : Array of blink-corrected trials (subject,cell/trial)
%       .saccadeTrial : Array of saccade-corrected trials (subject,cell/trial)
%       .saccadeOnset : Array of onset in samples of detected saccades (subject,cell/trial)
%       .moveTrial : Array of movement-corrected trials (subject,cell/trial)
%       .badTrials : Array of bad trials (subject,cell/trial).  1=bad.  For averages, it is the number of bad trials that went
%                       in.  See .avgNum for which are still bad after averaging.
%       .badChans : Array of corrected bad channels (subject,cell/trial,channel). -1 in a session or continuous file means still bad.
%                   Negative numbers in an average file means number of still bad channels that went into the average
%                   (or rather, were left out of the average).  NaN in an average file means still bad.
%    .pca
%        fields from ep_doPCA and ep_doPCAst steps.  See them for documentation.
%        Not affected by editing.
%    .recTime   : Time in samples of the start of the epoch (1-based) from the start of the session (cell)
%                 For averaged data, the earliest sample of the trials going into the average.  When unknown, first will be
%                 set arbitrarily to 1 and the remaining will be in order of the data file (e.g., 1, then 251, for 250 sample epochs, etc.).
%    .stims     : cell string containing images of the screens of stimuli used in the experiment as loaded in by the imread command.  The events structure
%                 provides the timing information, via a key field where the .code field contains "stim" and the .data field contains the name of the file.
%                 When empty, it should still have the subfields.
%       .name   : name of the stimulus file including the suffix
%       .image  : the image matrix itself.
%    .calibration (can be empty with no subfields)
%       .ET          : Calibration values for eye-tracking channels
%           .Xzero   : zero correction for XEY channel
%           .Yzero   : zero correction for YEY channel
%           .Xscale  : scale correction for XEY channel
%           .Yscale  : scale correction for YEY channel
%       .SAC         : Calibration values for saccade channels
%           .Xzero   : zero correction for Hsaccade channel
%           .Yzero   : zero correction for Vsaccade channel
%           .Xscale  : scale correction for Hsaccade channel
%           .Yscale  : scale correction for Vsaccade channel
%     .impedances
%        .channels   : impedance values of the channels (chan,subject) (can be empty)
%        .ground     : impedance values of the ground electrode (subject) (can be empty)
%  origEloc:    The raw eloc prior to any editing.
%
% With respect to the FieldTrip data structure:
%  According to Vladimir Litvak,
%  "In case of events with duration that define trials event.sample is the first sample of a trial and event.offset is
%  the offset of the trigger with respect to the trial. An offset of 0 means that the first sample of the trial
%  corresponds to the trigger. A  positive offset indicates that the first sample is later than the trigger,
%  a negative offset indicates that the trial begins before the trigger."
%
%  According to the ft_read_event header documentation,the convention is to have a "trial" event that essentially keeps track of the recording time
%  for the first sample of the epoch and thus can be used to relate the event time to the proper sample in each epoch.
%  The .offset field is used for the "trigger" event only and indicates when the trigger really happened with respect to
%  the epoch (given that the trial .sample is being used to indicate the real time start time of the epoch rather than
%  the timing of the trigger event per se.
%
%  Since this approach does not lend itself well to the EP Toolkit environment, as of 2.40 the convention will be for the .sample field to be based on epoch time.
%  A recTime field will enable translation back to recording time and file time if needed.
%
%  EEGlab files don't really have set conventions but to the extent that there is one,
%  the events should have both a 'value' field to denote the generic type of event,
%  as in 'trigger', and a 'type' field to denote the nature of this generic event,
%  as in the condition of the experiment.  EEGlab essentially ignores
%  'value' fields and seems to have them mostly for purposes of using
%  FieldTrip's fileio for importing files.
%  The unofficial norm in FieldTrip analyses seems to be to use the generic
%  'type' field to subselect the event type of interest (like 'trigger')
%  and then use the 'value' field to define the conditions, so essentially the opposite to EEGlab.
%  The EP Toolkit hews to the FieldTrip convention with regard to the
%  'value' and the 'type' fields.
%
% With respect to spectral data:
%  There are a lot of things about spectral data that are merely conventions so there doesn't seem to be a definition or central authority to refer to.
%  After canvassing the literature, I've chosen to recognize the following as being accepted practice.
%  The internal spectral data are in complex form.  This is so that the imaginary component of the FFT can be represented for its phase information.
%  Prior to read-outs or analysis outputs, the data are always converted to spectral density form.
%  If the power (pw) or dB option is chosen, then the data are converted to power.
%  If the dB option is chosen, then the data are converted to dB form.  By definition, dB is always power. 
%     Smith, SW (2003) Digital Signal Processing.  Newnes: Amsterdam.  pp 263-264.
%  By convention, the spectral density conversion is applied to the power form.
%  So for example, power spectral density is computed by dividing by bin herz width and amplitude spectral density is computed by dividing by the square root of bin herz width.
%  Spectral density is computed prior to conversion to dB.
%  Internally frequency data are represented as complex numbers to preserve the phase information (both the .data and the .freqData fields).
%  However, when data are added together (as in averaging) they need to be converted to absolute amplitude form otherwise opposite phase data will cancel out.
%  So such data are represented internally as absolute amplitude.  There is no flag other than the nature of the numbers themselves (complex or real).
%  The exception is adding channels together as it is appropriate for opposite phase to cancel out to reflect reference channel effects.
%
%  Example:
%  [EPdata]=readData;
%  [EPdata]=readData('file','NT9coll.ave_af.egis','cells',[1:3]);
%  [EPdata]=readData('format','ns_avg');

%History:
%  by Joseph Dien (2/7/08)
%  jdien07@mac.com
%
%  bugfix 3/31/08 JD
%  Handles formats other than EGIS properly by adding information about number of cells, subjects, and ordering.
%
%  modified 4/1/08 JD
%  Added support for EGI's segmented simple binary format.
%
%  modified 4/30/08 JD
%  accepts EGIS session files.
%
%  modified 11/4/08 JD
%  outputs structured variable including montage and file format and name
%  information.  Also, data is now a 4D data matrix.  Also, now autodetects
%  the electrode montage in EGIS files.
%
%  modified 11/10/08 JD
%  inputs changed to keyword list approach.
%
%  modified and bugfix (1/31/09) JD
%  Added initEP.  Added version and date fields to output structure.  Fixed reversal of order for subjects and cells
%  specified with format keyword.
%
%  modified (2/3/09) JD
%  Added Channels keyword.  Samples changed to set of numbers rather than just first and last.  Format no longer
%  includes fields for numbers of cells and subjects.
%
%  modified and bugfix (3/25/09) JD
%  Added information on date, version, sampling rate, baseline, experiment name, data type, and specs.
%  Subject names for EGIS average files now correctly detected.
%  Generalized to handle session files too with addition of trial specs, trial names, and events.
%  Also handles factor output files.  Tested out with Neuroscan files.
%  Added .ced field for electrode coordinate information.
%
%  modified (4/17/09) JD
%  Groups epochs in single trial files by cell.
%  eloc and chanType and implicit fields added to the files.
%
%  modified (4/27/09) JD
%  Added subNum field.  Merged 'subject_average', 'combined_average' and 'grand_average' data types into 'average'.
%  Deblanks cell names for EGIS files.
%
%  modified and bugfix 7/22/09 JD
%  Eliminated the latency field from the event structure.  Dropping initial samples results in event offset being
%  adjusted.   Dropped events of "trial" type as they are redundant with cell name information.  Added cellTypes field.
%  Fixed some crashes when handling simple binary files.  Added factors as 5th dimension.  Added support for text files.
%  Added support for EP files.  Added analysis subfields blinkTrial, moveTrial, badTrials, and badChans.
%  Added support for factor file compression via facVecT and facVecS fields.  Added dataName and facData fields.
%  EGIS subject specs converted to strings.
%
%  bugfix 8/4/09 JD
%  Added workaround for defective EGIS average headers generated by NetStation (the LastDone field which indicates number of subjects is incorrect).
%  Also, fix for NetStation not outputing subject numbers in the EGIS average header.
%
%  bugfix 8/29/09 JD
%  For "other" format files, not distinguishing cell names that start with the same characters (e.g., "1" and "11" not
%  distinguished).  Thanks to Grega Repovs.
%  For single_trial data, was incorrectly reporting on the command line the wrong number of cells being read in.
%
%  modified 8/30/09 JD
%  Added settings for text files.
%
%  bugfix 9/5/09 JD
%  Crash when reading factor file other than EGIS format.
%  Header defective (led to subsequent crash) when reading simple binary average file.
%  Crash when setting reference channel type.
%  Unable to find ced file when already specified in the file.
%  Text files with non-number columns sometimes not being read properly.
%
%  bugfix 9/16/09 JD
%  Text file import was skipping the first timepoint.
%  Not reading correct number of text columns.
%
%  bugfix 10/12/09 JD
%  Crash when reading "other" format files in single file mode where event values are empty.
%
%  bugfix 10/18/09 JD
%  Now aborts when NaN or inf values are in the data.
%
%  modified 10/31/09 JD
%  Reads in eloc information if present in eeglab_set file.
%
%  bugfix 11/9/09 JD
%  Number of subjects for continuous simple binary files incorrectly set.
%  Crash when reading continuous simple binary files.
%
%  modified 11/15/09 JD
%  When reading simple binary files, if the cell names cannot be deduced, put all of the segments into the same single
%  cell rather than just aborting.
%
% bugfix and modified 11/20/09 JD
% Replaced "union" commands with "unique" commands because certain situations caused the "union" command to crash in
% Matlab 2007.
% When importing a text file, ignore tabs at the end of the line.
% Drops CELL and TRSP events for all simple binary files since NetStation loses the associated information when exporting simple binary files.
%
%  modified 1/15/10 JD
%  When reading in an EGIS format file, if none of the EGI montages are chosen, then it will allow for a .ced file to be chosen.
%
%  modified 1/26/10 JD
%  When reading in file formats that label the channels, they will be reordered to match the order of the ced file.  Also, if any channels
%  are missing from the data file then they will be added as bad channels.  Also, can now accommodate channels that are in the data file
%  but not in the ced file.  The type field in the .ced files is now used since the
%  EEGlab bug was apparently fixed and it is now functional.  The type field must now be present in the .ced file and assumptions
%  will no longer be made about which ones are REF or FID types.
%
%  bugfix 1/30/10
%  Fixed incorrect subject ID being extracted from EGIS session headers.
%  Fixed crash due to events with empty value fields (now uses type field instead if value field is empty).
%  Fixed crash due to empty type fields in eloc information by assuming they are EEG channels.
%
%  bugfix 2/4/10
%  Fixed crash when loading factor file with CED information.
%  Fixed loss of some factor information when loading EP format factors file.
%
%  modified & bugfix 2/27/10 JD
%  When importing data, epochs that are entirely flat across all channels are marked as bad.
%  Fixed crash when two REF channels (as in M1-M2 mean mastoid channels).
%  Analysis fields no longer optional.
%  Fixed crash when loading in an EP file with ced information.
%  Fixed crash when loading in a .set file that does not have ced information included in its header.
%  Initializes subject spec names as cell rather than empty array.
%  When importing EP file format data, correctly checks for data type even for "factors" and "grand_average".
%  Eliminated chantype field for implicit channels.
%  Now reads the nsweeps field of Neuroscan AVG files to fill in the avgNum and subNum fields.
%  Implicit reference channels will no longer be marked as bad.
%
%  bugfix 3/5/10 JD
%  One dimensional cell array fields are now standardized to be column vectors.  Some file formats were causing crashes
%  in the Edit function because some of these fields were coming out as row vectors.
%
%  modified 3/15/10 JD
%  if .type field is missing from eloc information, then add it (assume channels are "EEG").
%
%  bugfix 3/27/10 JD
%  Fixed crash when loading in .txt file format data.
%  Fixed crash when loading in data where event information is empty.
%  Fixed text files treating all channels as being implicit when used with ced files with channel names different than
%  the default channel names.
%  For file formats with fixed channel orders, use channel names from ced file if available.
%  Tried to reinstate support for .edf file format.
%
%  bugfix 4/18/10 JD
%  Fixed wrong number of channel names when reading fixed order file formats (like EGIS) and the reference channel is implicit.
%
%  bugfix 5/1/10 JD
%  Fixed crash when reading data using single file mode with .set files which contain the .eloc information.
%  Fixed crash when reading data file that was originally in .set format and then was saved in EP file format.
%
%  modified 5/12/10 JD
%  For files with ced set to "eeglab", change to either name in chaninfo.filename field if present, else "none".
%
%  bugfix 5/22/10 JD
%  Fixed crash when importing .set file or EP file with unavailable or invalid ced file named in ced field.
%  Fixed channel selection not operating on channel coordinates eloc field.
%
%  bugfix 6/17/10 JD
%  Fixed not keeping FID channels in implicit channel info for fixed channel order file formats (e.g., EGIS, text).
%  Fixed electrode information not matching the data for second file onwards when reading files in single cell file
%  mode.
%  Can now be passed eloc information so don't have to access ced file on every pass for single cell file mode.
%  Fixed cell labels and sub labels not being applied to EP file formats.
%  Fixed crash for file formats with fixed-order channels when ced file has wrong number of channels.
%  Fixed not ignoring extra tab at end of line of text files, resulting in "not-a-number" errors.
%  Fixed crash when loading EP or .set file with name of original ced file in addition to eloc information.
%  Fixed losing the electrode coordinate information (eloc and and ced) when loading in a .study file.
%
%  bugfix 8/2/10 JD
%  Fixed crash if using "trial" events to determine name of cells and the .value field is empty.
%
%  bugfix 8/25/10 JD
%  Fixed not obtaining channel names from eloc info when provided by function call, as when merging files during Single File mode.
%  Fixed errors when importing .study files where the group or the session fields were left blank.
%
%  bugfix 10/3/10 JD
%  Fixed events in the final sample of a segment being assigned to the succeeding segment and causing a crash if the
%  segment was already the last one.
%  Fixed crash when the montage keyword was followed by a blank, as when preprocessing a batch of non-EGIS files
%  containing more than one data file.
%
%  modified 10/12/10 JD
%  For continuous files, analysis fields refer to one second epochs.
%
%  modified 10/16/10 JD
%  Added support for HEOG saccade correction.
%
%  bugfix 11/3/10 JD
%  Now accepts upper case file suffixes (e.g., .EGIS).
%
%  bugfix 12/6/10 JD
%  When reading in EGIS session files with implicit ref, adds it in explicitly, thus avoiding crash in eyeblink
%  correction routine.
%
%  modified 2/15/11 JD
%  Added support for CSV text files in addition to tab-delimited text files.
%
% modified 2/26/11 JD
% Added support for EGI Matlab files
%
% modified 5/21/11 JD
% Added support for 6th dimension of frequency for time-frequency analyses
%
%  bugfix 6/26/11 JD
%  Eliminated crash when loading EP format file with both electrode coordinates and a regional channel.
%
% bugfix 12/19/11 JD
% Fixed not assigning cell names to .set files generated by ERPlab.
%
% modified 1/26/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
%
% modified 1/29/12 JD
% When loading older EP files, if data had implicit reference and was EGI data, then converted to explicit reference.
%
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% modified 3/25/12 JD
% Added support for freq label when reading in text files.
%
% bugfix 6/6/12 JD
% Fixed crash when loading in EP format file with no implicit channels.
%
% bugfix 9/4/12 JD
% Fixed crash when loading older EP format files with factor data.
%
% modified 9/8/12 JD
% Improved ability to figure out the cell names of EEGlab files.
%
% modified 10/16/12 JD
% Added option to set last row to be imported in text files.
%
% bugfix 10/18/12 JD
% Fixed subNames field not necessarily being a column vector.
%
% modified 1/10/13 JD
% Added power field to the EP format.
%
% bugfix 1/17/13 JD
% Events assigned to wrong subject (off by one).
%
% modified 1/17/13 JD
% Added support for .set files generated by Widmann's pop_grandaverage function.
% Added support for ERPlab .erp file format.
%
% bugfix 1/18/13 JD
% Fixed erroneous "labels" error message when trying to load .study file.
%
% bugfix 1/23/13 JD
% Fixed problem that loading EP files with frequency PCA data results in damaged data file and error messages.
%
% bugfix 2/7/13 JD
% Fixed crash when loading in .study average file where the .set files have no prestimulus period.
%
% bugfix 2/20/13 JD
% Fixed crash when data file has no event information.
%
% bugfix 4/24/13 JD
% Better handles ced files where there are channel types other than EEG, FID, and REF.
%
% bugfix 5/8/13 JD
% Fixed waveforms in simple binary average files getting scrambled!  It no longer makes assumptions about the ordering
% of the cells.
% Fixed bug where if baseline was zero ms then instead it was being recorded as being 4 ms.
% Fixed not handling Neuroscan files with two physically linked explicit reference sites.
%
% modified 5/9/13 JD
% Added support for EGI's epoch-marked simple binary format for session files.
%
% bugfix 5/20/13 JD
% Fixed crashing when merging fixed channel files (like Neuroscan) where there are channels in the data that are not in the CED file.
%
% modified 5/20/13 JD
% Added option to eliminate unwanted channels, as in a GFP channel, in the data by marking the channel as BAD in the CED file.
% Added original eloc output.
%
% modified 5/26/13 JD
% Text files with multiple delimiters between values (as in space-space) now treated as a single delimiter.
%
% bugfix 6/24/13 JD
% Fixed crash when segmented simple binary session file is incorrectly specified to be an average file by the user.
%
% bugfix 9/18/13 JD
% Fixed incorrect inference of reference scheme when a single ref channel is designated.
%
% bugfix 9/20/13 JD
% Fixed channel type not changed from REF to EEG for file formats with flexible order channels.
% Fixed REF channels assumed to be last for fixed order channel file formats.
% General overhaul of CED code to better support MEG and ANS chan types, BAD CED code, and added ECG chan type.
%
% modified 9/22/13 JD
% Added support for mff files with recording stops and for ep_recordingStart
% events.
%
% bugfix 9/24/13 JD
% Workaround for EEGlab issue where if a CED file has just the label and the type filled out, the type info migrates over to the theta column for some reason.
%
% modified 9/26/13 JD
% Improved detection of Simple Binary files with scrambled cells.
%
% bugfix 10/2/13 JD
% Fixed crash if cell names are a mix of numbers and strings.
% Fixed crash if type field from CED file contains numbers for some reason.
%
% modified 10/10/13 JD
% Added recTime field.
% No longer rearranging single trial data to group by cell.
% Eliminated offset field from events structure.
% Fixed hitory field getting overwritten.
%
% modified 10/16/13 JD
% Epochs with boundary events are marked as bad.
%
% modified 10/21/13 JD
% Added support for reading rejected fields from EEGlab files.
% Added support for reading nTrials and chanlocs fields from ERPlab files.
%
% bugfix 10/21/13 JD
% Ensures that power comes after analysis field so order of fields is always the same.
%
% bugfix 10/29/13 JD
% Fixed crash when reading file with TRSP events that was not in mff format.
%
% modified 10/30/13 JD
% "boundary" in .type as well as .value fields for boundary events.
% Added .keys field to events.
% Added subject specs support for mff files.
%
% bugfix 11/14/13 JD
% Shifted data selection code to ep_selectData to update implementation of data selection protocols.
%
% bugfix 11/22/13 JD
% No longer redo eloc when loading in EP format files and ced is specified to be the same as what it already is.
% No longer extract trial spec names from mff continuous files.
%
% bugfix 12/24/13 JD
% Fixed crash when keyword list includes empty cells.
%
% bugfix 1/9/14 JD
% Fixed crash when reading mff file with more than one subject field.
%
% modified 1/19/14 JD
% When reading in boundary events in mff files, duration field contains the length of the recording pause.
%
% modified 2/26/14 JD
% pca no longer optional field.  No fields are optional.
%
% bugfix 2/27/14 JD
% Fixed crash for EGIS average files with custom cell header lengths.
%
% bugfix 3/13/14 JD
% Fixed crash when reading mff file with subject field where the field was left blank.
%
% modified 3/13/14 JD
% Added support for reading FIFF files.
% ced label for electrode coordinates provided by file (e.g., eeglab, MFF, FIFF formats) is "internal".
%
% modified 3/14/14 JD
% Eliminated file type check for EGIS files since NetStation generates average EGIS headers that are incorrectly marked as being session files.
%
% modified 3/16/14 JD
% Uses internal electrode coordinates provided by MFF and FIFF files.  Added elecPrefs.
%
% bufix 3/19/14 JD
% Fixed recTime field not including space for the 'all' cell in PCA output, resulting in crashes when edited.
%
% modified 3/24/14 JD
% Added .cov field.  Added bad channel info for FIFF files.
%
% bufix 4/8/14 JD
% Fixed keys field of events not being added when missing, resulting in EP files created by older versions of Toolkit
% not being usable.
% Fixed not putting factor variance information in correct location when loading PCA .ept files,
% resulting in "There are 0 factors that meet the minimum variance criterion" error messages when trying to autoPCA them.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% bufix 4/25/14 JD
% Added conversion of REF channel type to EEG for older EP files.
%
% bufix 4/28/14 JD
% Fixed not able to load in EP files with frequency PCA data.
%
% bufix 5/21/14 JD
% Fixed crash when loading an ept file with no theta values for the electrode coordinates.
%
% bufix 6/12/14 JD
% Fixed blank keys field of events being produced without .key (e.g., .keys.keyCode instead of .keys.key.keyCode)
%
% bufix 6/29/14 JD
% Added fix for crash when trying to read mff files that erroneously label their COM channel as being a reference channel.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bufix 7/31/14 JD
% Fixed .cov.Nq not being updated when stripping off subject adds.
%
% modified 8/24/14 JD
% Workaround for NetStation bug where EGIS average files have number of trials in cell in wrong spot in the cell header.
%
% modified 10/2/14 JD
% Added support for reading in event text files when present.
% Added support for bdf recording stop events (by marking them as 'boundary' events).
%
% bufix 3/20/15 JD
% Fixed crash when FontSize not provided as in reading .study files.
%
% modified 5/29/15 JD
% Added support for reading edf files with channels of varying sampling rates.
%
% bufix 7/5/15 JD
% Fixed crash when reading data with an event prior to first sample of data.
%
% bufix 9/25/15 JD
% Fixed crash when preferences set to rotate mff/eeglab electrode coordinates 180
% or 270 degrees.
%
% bufix 10/9/15 JD
% Fixed crash when mff average file has only one subject.
%
% modified 10/13/15 JD
% Standardized electrode labeling (e.g., E10) when reading mff average files.
%
% bufix 10/19/15 JD
% Fixed crash when continuous set file has a single eventHdr event and it
% has an NaN or numeric .value.
% Handles the NaN values that EEGlab seems to set channels to when it
% automatically edits them to being bad data.
%
% bufix 12/10/15 JD
% Fixed crash when cell name was undefined and was therefore set to "cell001"
%
% modified 12/18/15 JD
% SMI word files now expected to end in .txt rather than _smi.txt
%
% bufix 2/3/16 JD
% Check to see if SMI info has already been added.
%
% bufix 2/17/16 JD
% Fixed couldn't read text files if lastrow equalled zero and firstrow is larger than 1.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% modified 9/24/16 JD
% NaN values no longer set to zero and bad data flags set for non-EEG channels.
%
% modified 10/15/16 JD
% Added stims field.
%
% modified 11/5/16 JD
% Added support for reading subject spec text files.
%
% modified 11/13/16 JD
% Added support for .calibration field
%
% modified 3/16/17 JD
% Added support for New Segment events in brainvision files.
%
% bufix 5/19/17 JD
% Adds .stims and .calibration fields to older .ept PCA data files to avoid crashes in two-step PCA.
% When reading NS5 continuous mff files, now recognizes reference channel type correctly.
% Added option to flip the electrode locations for .mff and .fiff files.
% Fixed crash when rotating electrode coordinates 180 or 270 degrees for .mff and .fiff files.
%
% modified 5/30/17 JD
% When importing an average file with multiple cells with the same name, modifies the names to be unique rather than assuming single-file type.
%
% modified & bugfix 6/19/17 JD
% Added support for flexible segments with addition of the .timeUnits field.
% Presence of NaN in EEG channels no longer zeroed as bad data if original file was an .ept file (as fix for zeroing out plv transform).
%
% modified & bugfix  11/28/17 JD
% Eliminated restrictions on location of CED files.
% Added support for impedances field.
% Fixed not recognizing ECG channels in mff files.
%
% bufix 12/7/17 JD
% Fixed boundary events not being handled correctly in BrainVision eeg files.
% Fixed odd behavior due to Matlab bug in which str2num executes word strings which are function names.
%
% bugfix 3/1/18 JD
% Fixed unable to use eloc information from a batch of .ept files when they originally used different ced files.
%
% bugfix 3/18/18 JD
% Now handles situation where one of the directories has a space after its name.
% Fixed crash when .ept file has event with empty key.
%
% modified 5/13/18 JD
% No longer adds every read event to the .history field.
%
% modified 5/25/18 JD
% Added preference option to turn off BV header encoding.
%
% bugfix 6/12/18 JD
% Fixed not ensuring subject spec names and trial spec names are column vectors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EPdata=[];
origEloc=[];

samples=[];
cells=[];
subjects=[];
channels=[];
EPdata=[];
montage=[];
fileFormat=[];
dataType=[];
chanNames=[];
timeNames=[];
timeUnits='';
cellNames=[];
trialNames=[];
subNames=[];
facNames=[];
freqNames=[];
relNames=[];
ename=[];
dataName=[];
trialSpecs=[];
trialSpecNames=[];
subjectSpecs={};
subjectSpecNames={};
eventHdr=[];
avgNum=[];
covNum=[];
subNum=[];
baseline=[];
sampleRate=[];
events=[];
ced=[];
eloc=[];
implicit=eloc;
chanTypes=[];
cellTypes=[];
subTypes=[];
facTypes=[];
facVecT=[];
facVecS=[];
facVecF=[];
facVar=[];
facVarQ=[];
origRefChan=[];
currRefChan=[];

fileName=[];
fileFormat=[];
numSubs=[];
numFacs=[];
numFreqs=[];
numCells=[];
numWaves=[];
cellLabels=[];
subLabels=[];
freqLabels=[];
silent='off';
refChan=[];
blinkTrial=[];
saccadeTrial=[];
saccadeOnset=[];
moveTrial=[];
badTrials=[];
badChans=[];
noise=[];
stanDev=[];
stanDevCM=[];
facData=[];
textPrefs.firstRow=1;
textPrefs.lastRow=0;
textPrefs.firstCol=1;
textPrefs.lastCol=0;
textPrefs.orientation=1;
reference.original=[];
reference.current=[];
reference.type='REG'; %default assumption is that data is regular reference
subjectsGrouped=0; %assume average data normally grouped by cells unless specified otherwise
recTime=[];
history=cell(0);
pca=[];
elecPrefs=[];
covMatrix=[];
covNq=[];
SMIsuffix='';
specSuffix='';
subjectSpecSuffix='';
scrsz=[];
FontSize=[];
stims=struct('name',{},'image',{});
calib=[];
impedances.channels=[];
impedances.ground=[];
BVheader=0;

if ~isempty(varargin)
    if isa(varargin{1},'cell') && nargin==1 %if keywords were input as a single cell string
        inputSet=varargin{1};
    else
        inputSet=varargin;
    end;
else
    inputSet=[];
end;

argNum=length(inputSet);
if mod(argNum,2) ~= 0
    msg{1}='The keywords need to all be in pairs, with a keyword followed by the keyword information.';
    [msg]=ep_errorMsg(msg);
    return
end;

argCount=1;
while argCount <= argNum
    if isempty(inputSet{argCount})
        argCount=argCount+1;
    else
        switch inputSet{argCount}
            case 'file'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''file'' keyword must be followed by a file name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''file'' keyword must be followed by a file name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                fileName=inputSet{argCount};
                argCount=argCount+1;
            case 'format'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''format'' keyword must be followed by a format name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''format'' keyword must be followed by a format name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                fileFormat=inputSet{argCount};
                argCount=argCount+1;
            case 'samples'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''samples'' keyword must be followed by a list of samples.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''samples'' keyword must be followed by a set of numbers (e.g., [1:250]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                samples=inputSet{argCount};
                argCount=argCount+1;
            case 'cells'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''cells'' keyword must be followed by a list of cells.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''cells'' keyword must be followed by a set of numbers (e.g., [1:3]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                cells=inputSet{argCount};
                argCount=argCount+1;
            case 'subjects'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''subjects'' keyword must be followed by a list of subjects.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''subjects'' keyword must be followed by a set of numbers (e.g., [1:3]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                subjects=inputSet{argCount};
                argCount=argCount+1;
            case 'channels'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''channels'' keyword must be followed by a list of channels.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''channels'' keyword must be followed by a set of numbers (e.g., [1:100 105:129]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                channels=inputSet{argCount};
                argCount=argCount+1;
            case 'montage'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''montage'' keyword must be followed by a montage name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount}) && ~isempty(inputSet{argCount})
                    msg{1}='The ''montage'' keyword must be followed by a montage name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                montage=inputSet{argCount};
                argCount=argCount+1;
            case 'name'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''name'' keyword must be followed by an experiment name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''name'' keyword must be followed by an experiment name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                ename=inputSet{argCount};
                argCount=argCount+1;
            case 'labels'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''labels'' keyword must be followed by a set of condition, subject, and freq labels, all within a set.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~iscell(inputSet{argCount})
                    msg{1}='The ''labels'' keyword must be followed by a set of condition, subject, and freq labels, all within a set.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if length(inputSet{argCount}) ~=3
                    msg{1}='The ''labels'' keyword must be followed by a set of condition, subject, and freq labels, all within a set.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                cellLabels=inputSet{argCount}{1};
                subLabels=inputSet{argCount}{2};
                freqLabels=inputSet{argCount}{3};
                if ~iscell(cellLabels)
                    cellLabels=cellstr(cellLabels);
                end;
                if ~iscell(subLabels)
                    subLabels=cellstr(subLabels);
                end;
                if ~iscell(freqLabels)
                    freqLabels=cellstr(freqLabels);
                end;
                argCount=argCount+1;
            case 'type'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''type'' keyword must be followed by a data type.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''type'' keyword must be followed by a data type.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~any(strcmp(inputSet{argCount},{'continuous','single_trial','average','grand_average', 'factors'}))
                    msg{1}='The data type must be one of the following: continuous, single_trial, average, grand_average, factors.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                dataType=inputSet{argCount};
                argCount=argCount+1;
            case 'NumberSubjects'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''NumberSubjects'' keyword must be followed by number of subjects or factors.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''NumberSubjects'' keyword must be followed by number of subjects or factors.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                numSubs=inputSet{argCount};
                argCount=argCount+1;
            case 'prestim'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''prestim'' keyword must be followed by the number of msec prior to the event onset.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''prestim'' keyword must be followed by the number of msec prior to the event onset.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                baseline=inputSet{argCount};
                argCount=argCount+1;
            case 'ced'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''ced'' keyword must be followed by the name of the .ced file.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''ced'' keyword must be followed by the name of the .ced file.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                ced=inputSet{argCount};
                argCount=argCount+1;
            case 'eloc'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''eloc'' keyword must be followed by the electrode information.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isempty(inputSet{argCount}) && ~isstruct(inputSet{argCount})
                    msg{1}='The ''eloc'' keyword must be followed by the electrode information.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                eloc=inputSet{argCount};
                if ~isempty(eloc)
                    implicit=eloc(1); %setup up implicit to have the same structure as eloc.
                    implicit(1)=[];
                end;
                argCount=argCount+1;
            case 'silent'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''silent'' keyword must be followed by on or off.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''silent'' keyword must be followed by on or off.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~any(strcmp(inputSet{argCount},{'on','off'}))
                    msg{1}='The ''silent'' keyword must be followed by on or off.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                silent=inputSet{argCount};
                argCount=argCount+1;
            case 'origReference'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''origReference'' keyword must be followed by the number of the original recording reference channel(s).';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                origRefChan=inputSet{argCount};
                argCount=argCount+1;
            case 'currReference'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''currReference'' keyword must be followed by the number of the current recording reference channel(s).';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                
                currRefChan=inputSet{argCount};
                argCount=argCount+1;
            case 'textPrefs'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''textPrefs'' keyword must be followed by a structured variable with the preferences.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isstruct(inputSet{argCount})
                    msg{1}='The ''textPrefs'' keyword must be followed by a structured variable with the preferences.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isfield(inputSet{argCount},'firstRow') || ~isfield(inputSet{argCount},'lastRow') || ~isfield(inputSet{argCount},'firstCol') || ~isfield(inputSet{argCount},'lastCol') || ~isfield(inputSet{argCount},'orientation')
                    msg{1}='The structured variable must have all five fields included.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                textPrefs=inputSet{argCount};
                argCount=argCount+1;
            case 'elecPrefs'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''elecPrefs'' keyword must be followed by a structured variable with the preferences.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''elecPrefs'' keyword must be followed by a number between 1 and 8.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if (inputSet{argCount} < 1) || (inputSet{argCount} > 8)
                    msg{1}='The ''elecPrefs'' keyword must be followed by a number between 1 and 8.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                elecPrefs=inputSet{argCount};
                argCount=argCount+1;
            case 'SMIsuffix'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''SMIsuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''SMIsuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                SMIsuffix=inputSet{argCount};
                argCount=argCount+1;
            case 'specSuffix'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''specSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''specSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                specSuffix=inputSet{argCount};
                argCount=argCount+1;
            case 'subjectSpecSuffix'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''subjectSpecSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''subjectSpecSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                subjectSpecSuffix=inputSet{argCount};
                argCount=argCount+1;
            case 'screenSize'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 4
                    msg{1}='The ''screenSize'' keyword must be followed by a vector with the four screen size values in pixels.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                scrsz=inputSet{argCount};
                argCount=argCount+1;
            case 'FontSize'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 1
                    msg{1}='The ''FontSize'' keyword must be followed by a keyword must be followed by a number.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                FontSize=inputSet{argCount};
                argCount=argCount+1;
            case 'BVheader'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 1 || ~ismember(inputSet{argCount},[0 1])
                    msg{1}='The ''BVheader'' keyword must be followed by a keyword must be followed by either zero or one.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                BVheader=inputSet{argCount};
                argCount=argCount+1;
            otherwise
                keyword=inputSet{argCount};
                if isnumeric(inputSet{argCount})
                    keyword=num2str(inputSet{argCount});
                end;
                msg{1}=[keyword ' is not a keyword.'];
                [msg]=ep_errorMsg(msg);
                return
        end;
    end;
end;

if isempty(scrsz)
    scrsz = get(0,'ScreenSize');
    if length(scrsz) < 2
        scrsz=[1 1 800 600];
    end;
end;

if isempty(FontSize)
    FontSize = 10;
end;

origEloc=eloc;

if isempty(fileName)
    [fileName, pathname] = uigetfile('*.*','Open:');
    if ~ischar(fileName)
        fileName=[];
    end;
    
    if (isempty(fileName))
        msg{1}='No filename selected. You have to click on a name';
        [msg]=ep_errorMsg(msg);
        return
    end
    fileName = [pathname fileName];
else
    if ~exist(fileName,'file')
        [pathstr, name, fileSuffix] = fileparts(fileName);
        if exist([fileName ' '],'file')
            msg{1}=['The file ' fileName ' has a space after its name which needs to be deleted.'];
        elseif ~exist([pathstr],'dir')
            msg{1}=['One of the directories containing ' fileName ' may have a space after its name which needs to be deleted.'];
        else
            msg{1}=['The file ' fileName ' does not exist.'];
        end;
        [msg]=ep_errorMsg(msg);
        return
    end;
end;

if exist(fileName,'file') && ~exist(deblank(fileName),'file')
    msg{1}=['The file ' fileName ' has a space after its name which needs to be deleted.'];
    [msg]=ep_errorMsg(msg);
    return
end;

%determine file type
[pathstr, name, fileSuffix] = fileparts(fileName);
[pathstr2, name2, fileSuffix2] = fileparts(name);

if any(strcmp(fileFormat,{'egi_egia','egi_egis'})) && ~any(strcmpi(fileSuffix,{'.egis','.ave','.gav','.raw','.ses'}))
    msg{1}='EGIS file names must end in .egis, .ave, .gav, .raw, or .ses to be recognized.';
    if any(strcmpi(fileSuffix2,{'.egis','.ave','.gav','.raw','.ses'}))
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if strcmp(fileFormat,'eeglab_erp') && ~strcmpi(fileSuffix,'.erp')
    msg{1}='ERPlab file names must end in .erp to be recognized.';
    if strcmpi(fileSuffix,'.erp')
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if strcmp(fileFormat,'ns_cnt') && ~strcmpi(fileSuffix,'.cnt')
    msg{1}='Neuroscan continuous file names must end in .cnt to be recognized.';
    if strcmpi(fileSuffix,'.cnt')
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if strcmp(fileFormat,'ns_eeg') && ~strcmpi(fileSuffix,'.eeg')
    msg{1}='Neuroscan single-trial file names must end in .cnt to be recognized.';
    if strcmpi(fileSuffix,'.eeg')
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if strcmp(fileFormat,'ns_avg') && ~strcmpi(fileSuffix,'.avg')
    msg{1}='Neuroscan average file names must end in .avg to be recognized.';
    if strcmpi(fileSuffix,'.avg')
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if strcmp(fileFormat,'egi_sbin') && ~any(strcmpi(fileSuffix,{'.sbin','.raw'}))
    msg{1}='SBIN file names must end in .raw or .sbin to be recognized.';
    if any(strcmp(fileSuffix2,{'.sbin','.raw'}))
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if strcmp(fileFormat,'ep_mat') && ~strcmpi(fileSuffix,'.ept')
    msg{1}='EP file names must end in .ept to be recognized.';
    if strcmp(fileSuffix2,'.ept')
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if any(strcmp(fileFormat,{'egi_mff_v1';'egi_mff_v2'})) && ~strcmpi(fileSuffix,'.mff')
    msg{1}='mff file names must end in .mff to be recognized.';
    if strcmp(fileSuffix2,'.mff')
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if any(strcmp(fileFormat,'edf')) && ~strcmpi(fileSuffix,'.edf')
    msg{1}='edf file names must end in .edf to be recognized.';
    if strcmp(fileSuffix2,'.edf')
        msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

if isempty(fileFormat)
    switch fileSuffix
        case '.sbin'
            fileFormat='egi_sbin';
        case '.egis'
            fileFormat='egi_egia';
        case '.mat'
            fileFormat='ep_mat';
        otherwise
            if strcmp(silent,'off')
                disp('Assuming file format is EGIS format.');
            end;
            fileFormat='egi_egia';
    end;
end;

if strcmp(fileFormat,'egi_egis') && ~strcmp(dataType,'single_trial')
    dataType='single_trial';
    if strcmp(silent,'off')
        disp('EGIS session files are always single trial format.');
    end;
end;

if strcmp(fileFormat,'egi_egia') && any(strcmp(dataType,{'single_trial','continuous'}))
    dataType='average';
    if strcmp(silent,'off')
        disp('EGIS average files are never single trial or continuous data types.  Defaulting to assume the file is an average file.');
    end;
end;

if strcmp(fileFormat,'eeglab_erp') && any(strcmp(dataType,{'single_trial','continuous'}))
    dataType='average';
    if strcmp(silent,'off')
        disp('ERPlab files are never single trial or continuous data types.  Defaulting to assume the file is an average file.');
    end;
end;

if strcmp(fileFormat,'biosemi_bdf') && any(strcmp(dataType,{'average','grand_average'}))
    dataType='continuous';
    if strcmp(silent,'off')
        disp('Only continuous and single trial bdf files are supported.  Defaulting to assume the file is a continuous file.');
    end;
end;

if strcmp(fileFormat,'edf') && ~strcmp(dataType,'continuous')
    dataType='continuous';
    if strcmp(silent,'off')
        disp('Only continuous edf files are supported.  Defaulting to assume the file is a continuous file.');
    end;
end;

if strcmp(fileFormat,'ep_mat')
    try
        tempVar=load('-mat', fileName);
        if isfield(tempVar,'EPdata')
            EPdata=tempVar.EPdata;
        else
            msg{1}=['The file ' fileName 'did not contain EP Toolkit data in it.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
        clear tempVar;
    catch
        msg{1}=['The attempt to load in the file ' fileName ' resulted in the error:' lasterr];
        [msg]=ep_errorMsg(msg);
        return
    end;
    try
        EPver=ver('EP_Toolkit');
        if str2num(EPdata.EPver.Version) > str2num(EPver.Version)
            msg{1}=['The file ' fileName 'was created by a more recent version of the EP Toolkit (' EPdata.EPver.Version ') than is running on this computer (' EPver.Version ') and therefore cannot be read.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
    catch
        EPver='unavailable'; %workaround for bug in earlier version of Matlab
    end;
    
    %backward compatibility conversions
    
    if ~isfield(EPdata,'freqNames')
        EPdata.freqNames=[];
    end;
    if ~isfield(EPdata,'facVecF')
        EPdata.facVecF=[];
    end;
    
    if ~isfield(EPdata,'pca')
        EPdata.pca=[];
    end;
    
    if ~isfield(EPdata.pca,'stims')
        EPdata.pca.stims=struct('name',{},'image',{});
    end;

    if ~isfield(EPdata.pca,'calibration')
        EPdata.pca.calibration=[];
    end;

    if ~isfield(EPdata,'facVar')
        if isfield(EPdata.pca,'facVar')
            EPdata.facVar=EPdata.pca.facVar;
        else
            EPdata.facVar=[];
        end;
    end;
    
    if ~isfield(EPdata,'facVarQ')
        if isfield(EPdata.pca,'facVarQ')
            EPdata.facVarQ=EPdata.pca.facVarQ;
        else
            EPdata.facVarQ=[];
        end;
    end;
    
    if ~isfield(EPdata,'relNames')
        EPdata.relNames=[];
    end;
    
    if ~isfield(EPdata,'stims')
        EPdata.stims=struct('name',{},'image',{});
    end;

    if ~isfield(EPdata,'timeUnits')
        EPdata.timeUnits=[];
    end;

    refChans=find(strcmp('REF',EPdata.chanTypes)); %assume that REF normally indicates original reference channels
    if ~isfield(EPdata,'reference')
        EPdata.reference.original=reference.original;
        EPdata.reference.current=reference.current;
        EPdata.reference.type=reference.type;
        
        if ~isempty(EPdata.freqNames) && ~isempty(EPdata.facNames) %if not spectral data or factor data
            if ~any(mean(EPdata.data(strcmp('EEG',EPdata.chanTypes),:,:,:,:,:),2))
                EPdata.reference.type='AVG'; %channels sum to zero so must be average reference
            end;
        end;
        
        if ~isempty(refChans)
            if length(refChans) > 2
                disp('More than two reference channels indicated.  Will ignore those past first two.');
                refChans=refChans(1:2);
            end;
            EPdata.reference.original=refChans;
            EPdata.chanTypes{refChans}='EEG'; %assume all REF channels are EEG channels
            if ~any(sum(ep_expandFacs(EPdata,refChans,[],[],[],[],[]),2))
                EPdata.reference.current=refChans; %reference channels still sum to zero so still being reference channels
            end;
        end;
    else
        if ~isempty(refChans)
            EPdata.chanTypes{refChans}='EEG'; %assume all REF channels are EEG channels
        end;
    end;
    
    if isfield(EPdata,'power')
        if strcmp(EPdata.power,'Power')
            EPdata=sqrt(EPdata.data);
            disp('Converting power data to amplitude data.');
        end;
        EPdata=rmfield(EPdata,'power');
    end;
    
    if ~isfield(EPdata,'recTime')
        EPdata.recTime=[1:EPdata.Fs:EPdata.Fs*(length(EPdata.cellNames)-1)+1]';
    elseif ~isempty(EPdata.facNames)
        if length(EPdata.recTime) == (length(EPdata.cellNames)-1)
            EPdata.recTime(end+1)=0; %repair effects of bug
        end;
    end;
    
    if ~isempty(EPdata.events)
        for i=1:size(EPdata.events,1)
            for k=1:size(EPdata.events,2)
                if ~isempty(EPdata.events{i,k})
                    if ~isfield(EPdata.events{i,k},'keys')
                        EPdata.events{i,k}(1).keys=struct('code','','data','','datatype','','description','');
                    elseif isfield(EPdata.events{i,k}(1).keys,'key')
                        if isfield(EPdata.events{i,k}(1).keys(1).key,'keyCode')
                            for iEvent=1:length(EPdata.events{i,k})
                                for iKey=1:length(EPdata.events{i,k}(iEvent).keys)
                                    if ~isempty(EPdata.events{i,k}(iEvent).keys(iKey).key)
                                        EPdata.events{i,k}(iEvent).keys(iKey).code=EPdata.events{i,k}(iEvent).keys(iKey).key.keyCode;
                                        EPdata.events{i,k}(iEvent).keys(iKey).data=EPdata.events{i,k}(iEvent).keys(iKey).key.data.data;
                                        EPdata.events{i,k}(iEvent).keys(iKey).datatype=EPdata.events{i,k}(iEvent).keys(iKey).key.data.dataType;
                                        EPdata.events{i,k}(iEvent).keys(iKey).description='';
                                    end;
                                end;
                            end;
                        end;
                    elseif (length(EPdata.events{i,k}(1).keys)>0) && isfield(EPdata.events{i,k}(1).keys(1),'keyCode')
                        for iEvent=1:length(EPdata.events{i,k})
                            newKeys=[];
                            for iKey=1:length(EPdata.events{i,k}(iEvent).keys)
                                newKeys(iKey).code=EPdata.events{i,k}(iEvent).keys(iKey).keyCode;
                                newKeys(iKey).data=EPdata.events{i,k}(iEvent).keys(iKey).data.data;
                                if isfield(EPdata.events{i,k}(iEvent).keys(iKey).data,'dataType')
                                    newKeys(iKey).datatype=EPdata.events{i,k}(iEvent).keys(iKey).data.dataType;
                                else
                                    newKeys(iKey).datatype='';
                                end;
                                if isfield(EPdata.events{i,k}(iEvent).keys(iKey).data,'description')
                                    newKeys(iKey).description=EPdata.events{i,k}(iEvent).keys(iKey).data.description;
                                else
                                    newKeys(iKey).description='';
                                end;
                            end;
                            EPdata.events{i,k}(iEvent).keys=newKeys;
                        end;
                    end;
                end;
            end;
        end;
    end;
    
    MEGchans=find(strcmp('MEG',EPdata.chanTypes));
    if ~isempty(MEGchans)
        chanTypes(MEGchans)='MGA';
        disp('Converting MEG channel types to MGA (axial gradiometer MEG).  If any of the MEG channels are actually planar or magnetometers, you will need to use the Edit function to correct them.');
    end;
    
    if ~isfield(EPdata,'covNum') || isempty(EPdata.covNum)
        EPdata.covNum=EPdata.avgNum;
    end;
    
    if xor(isempty(EPdata.trialSpecNames),isempty(EPdata.trialSpecs))
        EPdata.trialSpecNames=[];
        EPdata.trialSpecs=[];
    end;
    
    if ~isfield(EPdata,'calibration')
        EPdata.calibration=[];
    end;

    if ~isfield(EPdata,'impedances')
        EPdata.impedances=[];
    end;
    if ~isfield(EPdata.impedances,'channels')
        EPdata.impedances.channels=[];
    end;
    if ~isfield(EPdata.impedances,'ground')
        EPdata.impedances.ground=[];
    end;

    %extract analysis fields and backward conversions of them
    
    if isfield(EPdata,'analysis')
        if isfield(EPdata.analysis,'blinkTrial')
            blinkTrial= EPdata.analysis.blinkTrial;
            if strcmp(dataType,'continuous')
                if (size(EPdata.data,2) >= EPdata.Fs*2) && size(blinkTrial,2) ==1
                    blinkTrial = []; %backward compatibility conversion
                end;
            end;
        end;
        if isfield(EPdata.analysis,'saccadeTrial')
            saccadeTrial= EPdata.analysis.saccadeTrial;
            if strcmp(dataType,'continuous')
                if (size(EPdata.data,2) >= EPdata.Fs*2) && size(saccadeTrial,2) ==1
                    saccadeTrial = []; %backward compatibility conversion
                end;
            end;
        end;
        if isfield(EPdata.analysis,'saccadeOnset')
            saccadeOnset= EPdata.analysis.saccadeOnset;
            if strcmp(dataType,'continuous')
                if (size(EPdata.data,2) >= EPdata.Fs*2) && size(saccadeOnset,2) ==1
                    saccadeOnset = []; %backward compatibility conversion
                end;
            end;
        end;
        if isfield(EPdata.analysis,'moveTrial')
            moveTrial= EPdata.analysis.moveTrial;
            if strcmp(dataType,'continuous')
                if (size(EPdata.data,2) >= EPdata.Fs*2) && size(moveTrial,2) ==1
                    moveTrial = []; %backward compatibility conversion
                end;
            end;
        end;
        if isfield(EPdata.analysis,'badTrials')
            badTrials= EPdata.analysis.badTrials;
            if strcmp(dataType,'continuous')
                if (size(EPdata.data,2) >= EPdata.Fs*2) && size(badTrials,2) ==1
                    badTrials = []; %backward compatibility conversion
                end;
            end;
        end;
        if isfield(EPdata.analysis,'badChans')
            badChans= EPdata.analysis.badChans;
            if strcmp(dataType,'continuous')
                if (size(EPdata.data,2) >= EPdata.Fs*2) && size(badChans,2) ==1
                    badChans = []; %backward compatibility conversion
                end;
            end;
        end;
    end;
    
    %extract contents of the EP file.
    
    if isfield(EPdata,'cov')
        if ~isempty(EPdata.cov)
            covMatrix=EPdata.cov.covMatrix;
            covNq=EPdata.cov.Nq;
        end;
    end;
    
    sevenDdata=EPdata.data;
    timeNames=EPdata.timeNames;
    timeUnits=EPdata.timeUnits;
    subNames=EPdata.subNames;
    freqNames=EPdata.freqNames;
    numSubs=length(subNames);
    cellNames=EPdata.cellNames;
    numWaves=length(cellNames);
    numCells=length(unique(cellNames));
    numFreqs=length(freqNames);
    dataType=EPdata.dataType;
    trialNames=EPdata.trialNames;
    trialSpecs=EPdata.trialSpecs;
    avgNum=EPdata.avgNum;
    covNum=EPdata.covNum;
    subNum=EPdata.subNum;
    subjectSpecs=EPdata.subjectSpecs;
    chanNames=EPdata.chanNames;
    facNames=EPdata.facNames;
    relNames=EPdata.relNames;
    facTypes=EPdata.facTypes;
    chanTypes=EPdata.chanTypes;
    cellTypes=EPdata.cellTypes;
    subTypes=EPdata.subTypes;
    sampleRate=EPdata.Fs;
    trialSpecNames=EPdata.trialSpecNames;
    subjectSpecNames=EPdata.subjectSpecNames;
    if ~iscell(subjectSpecNames)
        subjectSpecNames={};
    end;
    facVecS=EPdata.facVecS;
    facVecT=EPdata.facVecT;
    facVecF=EPdata.facVecF;
    facData=EPdata.facData;
    if isempty(baseline)
        baseline=EPdata.baseline;
    end;
    reference=EPdata.reference;
    recTime=EPdata.recTime;
    history=EPdata.history;
    events=EPdata.events;
    pca=EPdata.pca;
    facVar=EPdata.facVar;
    facVarQ=EPdata.facVarQ;
    stims=EPdata.stims;
    calib=EPdata.calibration;
    impedances=EPdata.impedances;
    eloc=EPdata.eloc;
    ced=EPdata.ced;

%     if any(strcmp(ced,{'internal','none'})) || isempty(ced)
%         ced=EPdata.ced;
%         eloc=EPdata.eloc;
%     elseif strcmp(ced,EPdata.ced)
%         eloc=EPdata.eloc;
%     end;
    if isempty(montage)
        montage=EPdata.montage;
    end;
    
    if isfield(EPdata,'noise')
        noise= EPdata.noise;
    end;
    if isfield(EPdata,'std')
        stanDev= EPdata.std;
    end;
    if isfield(EPdata,'stdCM')
        stanDevCM= EPdata.stdCM;
    end;
    if isfield(EPdata,'implicit')
        if ~isempty(EPdata.implicit)
            implicit=EPdata.implicit;
            if isfield(implicit,'chantype')
                implicit=rmfield(implicit,'chantype'); %backward compatibility for older files with a chantype field.
            end;
            implicitTypes={implicit.type};
            implicitChan=find(strcmp('REF',implicitTypes));
            if ~isempty(EPdata.montage)
                if strcmp(EPdata.montage(1:3),'GSN') && (length(implicitChan) == 1) %if implicit reference and was EGI data, make explicit
                    
                    sevenDdata(end+1,:,:,:,:)=zeros(size(sevenDdata,2),size(sevenDdata,3),size(sevenDdata,4),size(sevenDdata,5),size(sevenDdata,6),size(sevenDdata,7));
                    if ~isempty(facVecS)
                        facVecS(end+1,:)=zeros(size(facVecS,2));
                    end;
                    
                    if ~isempty(facData)
                        facData(end+1,:,:,:,:)=zeros(size(facData,2),size(facData,3),size(facData,4),size(facData,5),size(facData,6),size(facData,7));
                    end;
                    
                    badChans(:,:,end+1)=zeros(size(badChans,1),size(badChans,2));
                    
                    if ~isempty(noise)
                        noise(end+1,:,:,:,:)=zeros(size(noise,2),size(noise,3),size(noise,4),size(noise,5));
                    end;
                    if ~isempty(stanDev)
                        stanDev(end+1,:,:,:,:)=zeros(size(stanDev,2),size(stanDev,3),size(stanDev,4),size(stanDev,5),size(stanDev,6));
                    end;
                    
                    chanTypes{end+1}='EEG';
                    chanNames{end+1}=implicit(implicitChan).labels;
                    eloc(end+1)=implicit(implicitChan);
                    implicit(implicitChan)=[];
                    reference.current(end+1)=length(eloc);
                    reference.original(end+1)=length(eloc);
                end;
            end;
        end;
    end;
    
    EPdata=[];
    
elseif (strcmp(fileFormat,'eeglab_set') && strcmp(fileSuffix,'.study'))
    try
        eval(['load(''-mat'',  ''' fileName ''');']);
    catch
        msg{1}=['The attempt to load in the file ' fileName 'resulted in the error:' lasterr];
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    fileList=cell(length(STUDY.datasetinfo),1);
    
    if (length(fileList) > 1) && strcmp(dataType,'continuous')
        msg{1}='Multiple files cannot be combined into a single continuous data file.  You will need to load them in separately.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if (length(unique({STUDY.datasetinfo.subject})) > 1) && strcmp(dataType,'single_trial')
        msg{1}='Multiple subjects cannot be combined into a single single-trial data file.  You will need to load them in separately.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    for setFile=1:length(STUDY.datasetinfo)
        theFile=[STUDY.datasetinfo(setFile).filepath filesep STUDY.datasetinfo(setFile).filename];
        if ~exist(theFile,'file')
            theFile=STUDY.datasetinfo(setFile).filename;
            if ~exist(theFile,'file')
                msg{1}=['The file ' theFile ' is missing.'];
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        fileList{setFile}=theFile;
    end;
    mergeName=STUDY.name;
    if isempty(mergeName)
        mergeName='mergedData';
    end;
    
    mergeArg=[];
    for file=1:length(fileList)
        mergeArg{file,1}=fileList{file};
        mergeArg{file,2}='format';
        mergeArg{file,3}='eeglab_set';
        mergeArg{file,4}='type';
        mergeArg{file,5}=dataType;
        mergeArg{file,6}='labels';
        
        theSession=STUDY.datasetinfo(file).session;
        if isempty(theSession)
            theSession='Cond';
        end;
        theCondition=STUDY.datasetinfo(file).condition;
        if isempty(theCondition)
            theSession='1';
        end;
        theCellName=[theSession '-' theCondition];
        
        theGroup=STUDY.datasetinfo(file).group;
        if isempty(theGroup)
            theGroup='Sub';
        end;
        theSubject=STUDY.datasetinfo(file).subject;
        if isempty(theSubject)
            theSubject='1';
        end;
        theSubjectName=[theGroup '-' theSubject];
        
        mergeArg{file,7}={theCellName theSubjectName cell(0)};
        mergeArg{file,8}='FontSize';
        mergeArg{file,9}=FontSize;
    end;
    
    [EPdata]=ep_mergeEPfiles(mergeArg,mergeName);
    if isempty(EPdata)
        return;
    end;
    hdr=[];
    sevenDdata=EPdata.data;
    ced=EPdata.ced;
    eloc=EPdata.eloc;
    timeUnits=EPdata.timeUnits;
    timeNames=EPdata.timeNames;
    subNames=EPdata.subNames;
    numSubs=length(subNames);
    cellNames=EPdata.cellNames;
    numWaves=length(cellNames);
    numCells=length(unique(cellNames));
    events=EPdata.events;
    trialNames=EPdata.trialNames;
    trialSpecs=EPdata.trialSpecs;
    avgNum=EPdata.avgNum;
    covNum=EPdata.covNum;
    subNum=EPdata.subNum;
    subjectSpecs=EPdata.subjectSpecs;
    chanNames=EPdata.chanNames;
    facNames=EPdata.facNames;
    freqNames=EPdata.freqNames;
    numFreqs=length(freqNames);
    facTypes=EPdata.facTypes;
    chanTypes=EPdata.chanTypes;
    cellTypes=EPdata.cellTypes;
    subTypes=EPdata.subTypes;
    sampleRate=EPdata.Fs;
    trialSpecNames=EPdata.trialSpecNames;
    subjectSpecNames=EPdata.subjectSpecNames;
    if ~iscell(subjectSpecNames)
        subjectSpecNames={};
    end;
    facVecS=EPdata.facVecS;
    facVecT=EPdata.facVecT;
    facVecF=EPdata.facVecF;
    facData=EPdata.facData;
    if isempty(baseline)
        baseline=EPdata.baseline;
    end;
    recTime=EPdata.recTime;
    pca=EPdata.pca;
else
    
    %determine data type
    if strcmp(fileFormat,'egi_egis') && isempty(dataType)
        dataType='single_trial';
    end;
    
    if strcmp(fileFormat,'egi_egia') && isempty(dataType)
        dataType='average';
    end;
    
    if strcmp(fileFormat,'ns_avg') && isempty(dataType)
        dataType='average';
    end;
    
    if strcmp(fileFormat,'ns_eeg') && isempty(dataType)
        dataType='single_trial';
    end;
    
    if strcmp(fileFormat,'ns_cnt') && isempty(dataType)
        dataType='continuous';
    end;
    
    if strcmp(fileFormat,'biosemi_bdf')
        dataType='continuous';
    end;
    
    if isempty(dataType)
        dataType='single_trial';
        if strcmp(silent,'off')
            disp('Defaulting to assuming the file is single trial segmented data from one subject.');
        end;
    end;
    
    %read in data
    if strcmp(fileFormat,'text')
        fid=fopen(fileName);
        for i=1:textPrefs.firstRow
            tempVar=fgetl(fid);
        end;
        delim='\t';
        if ~isempty(strfind(tempVar,',')) && isempty(strfind(tempVar,'\t'))
            delim=','; %if there are commas and no tabs, assume it is a comma-delimited file.
        elseif ~isempty(strfind(tempVar,' ')) && isempty(strfind(tempVar,'\t'))
            delim=' '; %if there are spaces and no tabs, assume it is a space-delimited file.
        end;
        
        numcols=length(regexp(tempVar,[delim '+']))+1; %determine number of columns based on number of delimiters (treating repeats as a single delimiter)
        if regexp(tempVar,[delim '$'])
            numcols=numcols-1; %if there is an extra tab at the end of the line, drop it.
        end;
        
        frewind(fid);
        for i=1:textPrefs.firstRow-1
            tempVar=fgetl(fid);
        end;
        
        lastCol=textPrefs.lastCol;
        if lastCol ==0
            lastCol=numcols;
        end;
        
        theRawData=textscan(fid, [repmat('%*s',1,textPrefs.firstCol-1) repmat('%f',1,lastCol-textPrefs.firstCol+1) repmat('%*s',1,numcols-lastCol)],'Delimiter',delim, 'MultipleDelimsAsOne', 1);
        fclose(fid);
        theData=cell2mat(theRawData);
        
        if isempty(theData)
            msg{1}='No data were read.  Did you check the EP Toolkit preference settings for reading text files to make sure they were set correctly?';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        lastRow=textPrefs.lastRow-textPrefs.firstRow+1;
        if lastRow < 1
            lastRow=size(theData,1);
        end;
        theData=theData(1:lastRow,:);
        
        if textPrefs.orientation ==1
            theData=theData';
        end;
        
        hdr.nChans=size(theData,1);
        hdr.nSamples=size(theData,2);
        hdr.nSamplesPre=0;
        hdr.nTrials=0;
        hdr.Fs=[];
        if ~isfield(hdr,'label')
            for i = 1:hdr.nChans
                hdr.label{i}  = ['e' num2str(i)];
            end;
        end;
        if ~isempty(freqLabels)
            if ~isempty(freqLabels{1})
                freqNames(1)=freqLabels{1};
            end;
        end;
    elseif strcmp(fileFormat,'ns_mat')
        [theNSFData,theEvents]=ep_readNSF(fileName);
        if isempty(theNSFData)
            msg{1}=['The attempt to load in the file ' fileName 'resulted in no data'];
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        theData=zeros(size(theNSFData{1,1},1),size(theNSFData{1,1},2),size(theNSFData,2));
        switch dataType
            case {'continuous','single_trial'}
                cellNames=theNSFData(2,:);
                for i=1:length(cellNames)
                    theIndex=strfind(cellNames{i},'_Segment');
                    cellNames{i}=cellNames{i}(1:theIndex(end)-1);
                    theData(:,:,i)=theNSFData{1,i};
                end;
            case {'average','grand_average','factors'}
                subNames=theNSFData(2,:);
                cellNames=theNSFData(2,:);
                for i=1:length(cellNames)
                    theIndex=strfind(cellNames{i},'_Average');
                    cellNames{i}=cellNames{i}(1:theIndex(end)-1);
                    subNames{i}=subNames{i}(theIndex(end)+8:end);
                    theData(:,:,i)=theNSFData{1,i};
                end;
        end;
        hdr.nChans=size(theData,1);
        hdr.nSamples=size(theData,2);
        hdr.nSamplesPre=0;
        hdr.nTrials=0;
        if ~isempty(sampleRate)
            hdr.Fs=sampleRate;
        else
            hdr.Fs=250;
            disp('Assuming default sample rate of 250Hz.');
        end;
        
        for i = 1:hdr.nChans
            hdr.label{i}  = ['e' num2str(i)];
        end;
        
        eventHdr=[];
        if (theEvents{2,end}/(1000/hdr.Fs)) < (size(theData,2)*size(theData,3))
            for i = 1:size(theEvents,2)
                eventHdr(i).value=theEvents{1,i};
                eventHdr(i).type=[];
                eventHdr(i).sample=theEvents{2,i}/(1000/hdr.Fs); %convert latency to samples from milliseconds.
            end;
        else
            disp('Some events fall outside the range of data points.  One possible explanation is that when NetStation is set to drop bad segments, it only drops');
            disp('them from the voltage data, not the event data.  If this was done, then it is no longer possible');
            disp('to determine which event data go with which voltage data.  In any case, the event data will therefore be ignored.');
        end;
        
    else
        try
            %            events = read_mff_event(fileName, hdr);
            %              if strcmp(fileFormat,'egi_mff_v2')
            %                  fileFormat='egi_mff_v1';
            %                 [tempHdr]=ft_read_header(fileName,'dataFormat',fileFormat, 'headerformat','egi_mff_v1'); %EGI's v2 routines don't provide all the header information yet.
            %                 hdr.orig=tempHdr.orig;
            %              end;
            
            if strcmp(fileFormat,'edf')
                [hdr]=ft_read_header(fileName,'dataFormat',fileFormat, 'headerformat',fileFormat,'chanindx',[1]);
                chanList=setdiff([1:hdr.orig.NS],hdr.orig.annotation);
                if isempty(chanList)
                    hdr=[];
                    eventHdr=[];
                else
                    hdr2=hdr;
                    maxFs=max(hdr.orig.SampleRate);
                    if ~all(hdr.orig.SampleRate==hdr.orig.SampleRate(1))
                        disp('EDF file has channels with different sampling rates.  Interpolating all to that of the highest rate.');
                        disp('May result in edge artifacts at end of epoch, which will be zeroed out.');
                    end;
                    maxSampLength=1000/maxFs;
                    theData=zeros(length(chanList),hdr.nSamples);
                    hdr.Fs=maxFs;
                    hdr.nChans=length(chanList);
                    for iChan=1:length(chanList)
                        theChan=chanList(iChan);
                        sampLength=1000/hdr.orig.SampleRate(theChan);
                        [theData(iChan,:)]=interp1([1:sampLength:hdr.orig.NRec*hdr.orig.Dur*1000],ft_read_data(fileName,'dataformat',fileFormat, 'headerformat',fileFormat,'chanindx',theChan),[1:maxSampLength:hdr.orig.NRec*hdr.orig.Dur*1000]);
                        hdrTmp=ft_read_header(fileName,'dataFormat',fileFormat, 'headerformat',fileFormat,'chanindx',theChan);
                        hdr.label{iChan}=hdrTmp.label{1};
                        hdr.chantype{iChan}=hdrTmp.chantype{1};
                        hdr.chanunit{iChan}=hdrTmp.chanunit{1};
                    end;
                    theData(isnan(theData))=0;
                    [eventHdr]=ft_read_event(fileName,'dataformat',fileFormat, 'headerformat',fileFormat,'eventformat',fileFormat,'chanindx',1,'header',hdr2);
                end;
            else
                if strcmp(fileFormat,'egi_mff_v2') && ~strcmp(dataType,'average')
                    [hdr]=ft_read_header(fileName,'dataFormat',fileFormat, 'headerformat',fileFormat);
                    eventHdr = read_mff_event(fileName, hdr);
                    [hdr2]=ft_read_header(fileName,'dataFormat','egi_mff_v1', 'headerformat','egi_mff_v1'); %v2 code doesn't provide subject specs and electrode coordinates yet
                    hdr.orig=hdr2.orig;
                else
                    [eventHdr]=ft_read_event(fileName,'dataformat',fileFormat, 'headerformat',fileFormat,'eventformat',fileFormat);
                    [hdr]=ft_read_header(fileName,'dataFormat',fileFormat, 'headerformat',fileFormat);
                end;
                if strcmp(fileFormat,'egi_mff_v1')
                    COMchan=find(strcmp('COM',hdr.label));
                    if ~isempty(COMchan) && (hdr.nChans ~= length(hdr.label))
                        %fix bug present in some mff files, resulting in the COM channel being identified in the header as an EEG channel.
                        hdr.label(COMchan)=[];
                        hdr.chantype(COMchan)=[];
                        hdr.chanunit(COMchan)=[];
                    end;
                end;
                
                if strcmp(fileFormat,'egi_mff_v2') || strcmp(fileFormat,'egi_mff_v1')
                    for iChan=1:length(hdr.label) %standardize mff channel labels
                        if ~isempty(str2num(hdr.label{iChan}))
                            newName=sprintf('E%s', hdr.label{iChan});
                            if isempty(find(strcmp(newName,hdr.label)))
                                hdr.label{iChan}=newName;
                            end;
                        elseif any(strcmp(hdr.label{iChan},{'Cz','REF','VREF'}))
                            newName=sprintf('E%s', num2str(iChan));
                            if isempty(find(strcmp(newName,hdr.label)))
                                hdr.label{iChan}=newName;
                            end;
                        end;
                    end;
                end;
                
                if strcmp(fileFormat,'egi_mff_v1') && strcmp(dataType,'continuous') %workaround for continuous mff files with recording stops not being read
                    if isfield(hdr.orig,'epochdef')
                        theData=zeros(hdr.nChans,hdr.nSamples);
                        if isempty(eventHdr)
                            eventHdr=struct('type',{},'sample',{},'value',{},'duration',{},'keys',struct('code','','data','','datatype','','description',''));
                        end;
                        if (size(hdr.orig.epochdef,1) > 2) && ~any(diff(hdr.orig.epochdef(:,2)-hdr.orig.epochdef(:,1)))
                            disp('Looks like this might actually be single trial data.');
                            [theData]=ft_read_data(fileName,'dataformat',fileFormat, 'headerformat',fileFormat);
                            dataType='single_trial';
                        else
                            for i=1:size(hdr.orig.epochdef,1)
                                [tempData]=ft_read_data(fileName,'dataformat',fileFormat, 'headerformat',fileFormat,'begsample',hdr.orig.epochdef(i,1),'endsample',hdr.orig.epochdef(i,2));
                                theData(:,hdr.orig.epochdef(i,1):hdr.orig.epochdef(i,2))=tempData;
                                if i > 1
                                    startEvents=find([eventHdr.sample] >= hdr.orig.epochdef(i,1));
                                    if ~isempty(startEvents)
                                        newEvent=startEvents(1);
                                        eventHdr=[eventHdr(1:newEvent-1) eventHdr(1) eventHdr(newEvent:end)];
                                    else
                                        newEvent=length(eventHdr)+1;
                                        eventHdr(end+1)=eventHdr(1);
                                    end;
                                    eventHdr(newEvent).type='boundary';
                                    eventHdr(newEvent).sample=hdr.orig.epochdef(i,1);
                                    eventHdr(newEvent).duration=hdr.orig.epochdef(i,3); %duration field holds length of recording pause in samples
                                    eventHdr(newEvent).value='boundary';
                                    eventHdr(newEvent).keys=struct('key',struct('code','','data','','datatype','','description',''));
                                end;
                            end;
                        end;
                    else
                        [theData]=ft_read_data(fileName,'dataformat',fileFormat, 'headerformat',fileFormat);
                    end;
                else
                    [theData]=ft_read_data(fileName,'dataformat',fileFormat, 'headerformat',fileFormat);
                end;
                
                if strcmp(fileFormat,'brainvision_eeg')
                    segEvents=find(strcmp('New Segment',{eventHdr.type}));
                    for iEvent=1:length(segEvents)
                        theEvent=segEvents(iEvent);
                        if eventHdr(theEvent).sample > 1
                            eventHdr(theEvent).type = 'boundary';
                            eventHdr(theEvent).duration=0;
                            eventHdr(theEvent).value='boundary';
                        end;
                    end;
                end;
                
            end;
        catch
            msg{1}='No data were read.  The error message was:';
            msg{2}=lasterr;
            if any(strcmp(fileFormat,{'egi_egis','egi_egia'}))
                msg{3}='Could you have specified the wrong file type?';
            end;
            if strcmp(fileFormat,'eeglab_set')
                msg{3}='EEGlab allows you to manually change the name of the .set file but not the .fdt file that goes with it.  Is that what happened here?';
            end;
            [msg]=ep_errorMsg(msg);
            return
        end;
    end;
    
    if ~any(theData)
        msg{1}='Error: The data contained nothing but zeroes.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if any(any(any(isinf(theData))))
        msg{1}='Error: The data contained infinite numbers.';
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if isfield(eventHdr,'sample')
        if any([eventHdr.sample] < 0)
            if ~strcmp(fileFormat,'egi_mff_v2') || any([eventHdr.sample] < -2)
                disp('Warning: The events data contained negative samples.  These events will be dropped.');
            end;
            if strcmp(fileFormat,'egi_mff_v2') && any([eventHdr.sample] == -1)
                disp(['There were ' num2str(length(find([eventHdr.sample] == -1))) ' events falling between epochs, which will be dropped.']);
            end;
            if strcmp(fileFormat,'egi_mff_v2') && any([eventHdr.sample] == -2)
                disp(['There were ' num2str(length(find([eventHdr.sample] == -2))) ' samples falling after the last epoch, which will be dropped.']);
            end;
            eventHdr=eventHdr(find([eventHdr.sample] > 0));
        end;
    else
        eventHdr.sample=[];
    end;
    
    if ~isfield(eventHdr,'value')
        eventHdr.value=[];
    end;
    
    if ~isfield(eventHdr,'type')
        eventHdr.type=[];
    end;

    for event=1:length(eventHdr)
        if ~isempty(eventHdr(event).value)
            eventValues{event}=eventHdr(event).value;
        else
            eventValues{event}=eventHdr(event).type;
        end;
    end;
    
    if isempty(theData) && strcmp(fileFormat,'egi_egia') %if it is actually  egi_egis, may return empty data
        fileFormat='egi_egis';
        dataType='single_trial';
        [theData]=ft_read_data(fileName,'dataformat',fileFormat, 'headerformat',fileFormat);
    end;
    
    if isempty(theData)
        msg{1}='No data were read.  Perhaps file format was misspecified?';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    % commented out.  Unfortunately, NetStation does not do this part of the EGIS header correctly.
    %     if (strcmp(fileFormat,'egi_egia') && hdr.orig.fhdr(2) == 3)
    %         disp('The file indicates that this is actually an EGIS session file.');
    %         fileFormat='egi_egis';
    %         dataType = 'single_trial';
    %         numSubs=1;
    %         [theData]=ft_read_data(fileName,'dataformat',fileFormat, 'headerformat',fileFormat);
    %     end;
    
    if (strcmp(fileFormat,'egi_egis') && hdr.orig.fhdr(2) == -1)
        disp('The file indicates that this is actually an EGIS average file.');
        fileFormat='egi_egia';
        if strcmp(dataType,'single_trial')
            dataType = 'average';
        end;
    end;
    
%     if strcmp(fileFormat,'egi_sbin') && any(strcmp(dataType,{'average','grand_average'})) && (hdr.orig.header_array(14) == 0)
%         disp('The file indicates that this is actually a continuous file.');
%         dataType = 'continuous';
%     end;
%     
    if strcmp(fileFormat,'biosemi_bdf') && ~isempty(strcmp('Status',hdr.label))
        statusChan=find(strcmp('Status',hdr.label)); %remove status chan from bdf files
        hdr.label(statusChan)=[];
        theData(statusChan,:)=[];
        hdr.chantype(statusChan)=[];
        hdr.chanunit(statusChan)=[];
        hdr.nChans=hdr.nChans-1;
    end;
    
    if strcmp(fileFormat,'neuromag_fif')
        if ~any(strcmp(dataType,{'average','grand_average'})) && hdr.orig.isaverage
            disp('The file indicates that this is actually an average file.');
            dataType = 'average';
        end
        if ~strcmp(dataType,'continuous') && hdr.orig.iscontinuous && ~hdr.orig.isepoched
            disp('The file indicates that this is actually a continuous file.');
            dataType = 'continuous';
        end
        if ~strcmp(dataType,'continuous') && hdr.orig.isepoched && ~hdr.orig.isaverage
            disp('The file indicates that this is actually a single trial file.');
            dataType = 'single_trial';
        end
        
        if isfield(hdr,'chantype')
            chanTypes=cellstr(repmat('EEG',size(theData,1),1));
            MGMchans=find(strcmp('megmag',hdr.chantype));
            if ~isempty(MGMchans)
                [chanTypes{MGMchans}]=deal('MGM');
            end;
            MGAchans=find(strcmp('megaxial',hdr.chantype));
            if ~isempty(MGAchans)
                [chanTypes{MGAchans}]=deal('MGA');
            end;
            MGPchans=find(strcmp('megplanar',hdr.chantype));
            if ~isempty(MGPchans)
                [chanTypes{MGPchans}]=deal('MGP');
            end;
            EEGchans=find(strcmp('eeg',hdr.chantype));
            if ~isempty(EEGchans)
                [chanTypes{EEGchans}]=deal('EEG');
            end;
            EOGchans=find(strcmp('eog',hdr.chantype));
            if ~isempty(EOGchans)
                [chanTypes{EOGchans}]=deal('REG'); %assumed to be bipolar array
            end;
            ATchans=find(strcmp('analog trigger',hdr.chantype));
            if ~isempty(ATchans)
                [chanTypes{ATchans}]=deal('BAD');
            end;
            DTchans=find(strcmp('digital trigger',hdr.chantype));
            if ~isempty(DTchans)
                [chanTypes{DTchans}]=deal('BAD');
            end;
        end;
    end;
    
    if strcmp(fileFormat,'brainvision_eeg') && strcmp(dataType,'continuous') && BVheader
%         %convert from S15-R15 convention to S255 convention if necessary.  %assuming event list has been sorted in order of sample.
%         eventHdr2=eventHdr(1);
%         eventHdr2(1)=[];
%         iEvent=1;
%         while iEvent <= length(eventHdr)
%             theSample=eventHdr(iEvent).sample;
%             eventList=find([eventHdr.sample]==theSample);
%             if length(eventList)==1
%                 if (length(eventHdr(iEvent).value) > 1) && strcmp(eventHdr(iEvent).value(1:2),'R ')
%                     eventHdr2(end+1)=eventHdr(iEvent);
%                     eventHdr2(end).value=['S ' num2str(bitshift(str2double(eventHdr(iEvent).value(2:end)),4))];
%                 else
%                     eventHdr2(end+1)=eventHdr(iEvent);
%                 end;
%                 iEvent=iEvent+1;
%             elseif length(eventList)==2
%                 if (length(eventHdr(eventList(1)).value) > 1) && (length(eventHdr(eventList(2)).value) > 1)
%                     if xor(strcmp('R ',eventHdr(eventList(1)).value(1:2)),strcmp('R ',eventHdr(eventList(2)).value(1:2))) && xor(strcmp('S ',eventHdr(eventList(1)).value(1:2)),strcmp('S ',eventHdr(eventList(2)).value(1:2)))
%                         if strcmp('R ',eventHdr(eventList(1)).value(1:2))
%                             eventR=1;
%                             eventS=2;
%                         else
%                             eventR=2;
%                             eventS=1;
%                         end;
%                         eventHdr2(end+1)=eventHdr(eventList(eventS));
%                         eventHdr2(end).value=['S ' num2str(bitshift(str2double(eventHdr(eventList(eventR)).value(2:end)),4)+str2double(eventHdr(eventList(eventS)).value(2:end)))];
%                     else
%                         eventHdr2(end+1)=eventHdr(iEvent);
%                         eventHdr2(end+1)=eventHdr(iEvent+1);
%                     end;
%                 else
%                     eventHdr2(end+1)=eventHdr(iEvent);
%                     eventHdr2(end+1)=eventHdr(iEvent+1);
%                 end;
%                 iEvent=iEvent+2;
%             else
%                 eventR=0;
%                 eventS=0;
%                 for iEventList=1:length(eventList)
%                     if (length(eventHdr(eventList(iEventList)).value) > 1) && strcmp('R ',eventHdr(eventList(iEventList)).value(1:2))
%                         if ~eventR
%                             eventR=eventList(iEventList);
%                         else
%                             disp('warning: problem with conversion of event codes');
%                         end;
%                     end;
%                     if (length(eventHdr(eventList(iEventList)).value) > 1) && strcmp('S ',eventHdr(eventList(iEventList)).value(1:2))
%                         if ~eventS
%                             eventS=eventList(iEventList);
%                         else
%                             disp('warning: problem with conversion of event codes');
%                         end;
%                     end;
%                 end;
%                 if ~eventS && ~eventR
%                     eventHdr2(end+1)=eventHdr(eventList(eventS));
%                     eventHdr2(end).value=['S ' num2str(bitshift(str2double(eventHdr(eventList(eventR)).value(2:end)),4)+str2double(eventHdr(eventList(eventS)).value(2:end)))];
%                     eventHdr2(end+1:end+length(eventList)-2)=eventHdr(eventList(~ismember(eventList,[eventR eventS])));
%                 else
%                     eventHdr2(end+1:end+length(eventList))=eventHdr(eventList(iEvent:iEvent+length(eventList)-1));
%                 end;
%                 iEvent=iEvent+length(eventList);
%             end;
%         end;
        
        %look for subject header using EP convention
        hdrStart=0;
        badFlag=0;
        trialSpecNames=cell(0);
        iEvent=1;
        while iEvent < length(eventHdr)-7 %minimum size of a valid subject header is eight
            if strcmp(eventHdr(iEvent).value,'S104') && strcmp(eventHdr(iEvent+1).value,'S100') && strcmp(eventHdr(iEvent+2).value,'S114') && (hdrStart==0)
                hdrStart=iEvent;
                iEvent=iEvent+3;
                numFields=str2double(eventHdr(iEvent).value(2:end));
                iEvent=iEvent+1;
                if isnan(numFields) || isempty(numFields) || (numFields < 1) || (round(numFields) ~= numFields)
                    disp('Header error.  Aborting effort to read it.');
                    badFlag=1;
                    continue
                end;
                for iField=1:numFields
                    if ~strcmp(eventHdr(iEvent).value(2:end),'  3')
                        disp('Header error.  Aborting effort to read it.');
                        continue
                    end;
                    [iEvent,charOut]=readBVheaderField(eventHdr,iEvent);
                    if ~ischar(charOut) || isempty(charOut)
                        disp('Header error.  Aborting effort to read it.');
                        badFlag=1;
                        continue
                    else
                        subjectSpecNames{end+1,1}=charOut;
                    end;
                    [iEvent,fieldOut]=readBVheaderField(eventHdr,iEvent);
                    if isnan(fieldOut)
                        disp('Header error.  Aborting effort to read it.');
                        badFlag=1;
                        continue
                    else
                        subjectSpecs{1,end+1}=fieldOut;
                    end;
                end;
                if badFlag
                    continue
                end;
                numSpecs=str2double(eventHdr(iEvent).value(2:end));
                iEvent=iEvent+1;
                if isnan(numSpecs) || isempty(numSpecs) || (numSpecs < 1) || (round(numSpecs) ~= numSpecs)
                    disp('Header error.  Aborting effort to read it.');
                    badFlag=1;
                    continue
                end;
                for iSpec=1:numSpecs
                    if ~strcmp(eventHdr(iEvent).value(2:end),'  3')
                        disp('Header error.  Aborting effort to read it.');
                        badFlag=1;
                        continue
                    end;
                    [iEvent,charOut]=readBVheaderField(eventHdr,iEvent);
                    if ~ischar(charOut) || isempty(charOut)
                        disp('Header error.  Aborting effort to read it.');
                        badFlag=1;
                        continue
                    else
                        trialSpecNames{end+1,1}=charOut;
                    end;
                end;
                if badFlag
                    continue
                end;
                if ~strcmp(eventHdr(iEvent).value,'S114') || ~strcmp(eventHdr(iEvent+1).value,'S100') || ~strcmp(eventHdr(iEvent+2).value,'S104')
                        disp('Header error.  Aborting effort to read it.');
                        badFlag=1;
                        continue
                end;
                eventHdr(hdrStart:iEvent+2)=[];
                continue
            else
                iEvent=iEvent+1;
            end;
        end;
        if badFlag
            trialSpecNames=[];
            subjectSpecs={};
            subjectSpecNames={};
            badFlag=0;
        else
            iEvent=1;
            newEventHdr=[];
            TRSPcount=0;
            numKeys=length(trialSpecNames);
            while iEvent <= length(eventHdr)
                if iEvent <= (length(eventHdr)+2)
                    if strcmp(eventHdr(iEvent).value,'S104') && strcmp(eventHdr(iEvent+1).value,'S100') && strcmp(eventHdr(iEvent+2).value,'S114')
                        disp(['***Warning: there is more than one header present.  You may want to trim this file before segmenting.***']);
                    end;
                end;
                if strcmp(eventHdr(iEvent).value,'S255')
                    TRSPcount=TRSPcount+1;
                    hdrStart=iEvent;
                    iEvent=iEvent+1;
                    if (iEvent > length(eventHdr)) || (length(trialSpecNames) ~= str2double(eventHdr(iEvent).value(2:end)))
                        disp(['TRSP error.  Aborting effort to read TRSP #' num2str(TRSPcount) '.']);
                        badFlag=1;
                        continue
                    end;
                    iEvent=iEvent+1; %skip over field indicating number of keys
                    theTRSP=[];
                    theTRSP.type='TRSP';
                    theTRSP.sample=eventHdr(iEvent).sample;
                    theTRSP.value='TRSP';
                    theTRSP.duration=0;
                    for iKey=1:numKeys
                        theTRSP.keys(iKey)=struct('code','','data','','datatype','','description','');
                        theTRSP.keys(iKey).code=trialSpecNames{iKey};
                        switch eventHdr(iEvent).value
                            case {'S  1','S  4'}
                                theTRSP.keys(iKey).datatype='short';
                            case {'S  2','S  5','S  6','S  7'}
                                theTRSP.keys(iKey).datatype='long';
                            case {'S  3','S  8'}
                                theTRSP.keys(iKey).datatype='text';
                            otherwise
                                disp(['TRSP anomaly in TRSP #' num2str(TRSPcount) '.']);
                                badFlag=1;
                        end;
                        [iEvent,fieldOut]=readBVheaderField(eventHdr,iEvent);
                        if isnan(fieldOut)
                            if strcmp(theTRSP.keys(1).data,'257') && (iKey==numKeys) %error due to several 1's in a row in the trigger line not being distinguished by PyCorder
                                disp(['Attempting fix in TRSP #' num2str(TRSPcount) '.']);
                                badFlag=0;
                                fieldOut=theTRSP.keys(numKeys-1).data;
                                for keyFix=numKeys:-1:3
                                    theTRSP.keys(keyFix).data=theTRSP.keys(keyFix-1).data;
                                    theTRSP.keys(keyFix).datatype=theTRSP.keys(keyFix-1).datatype;
                                end;
                                theTRSP.keys(2).data=1;
                                theTRSP.keys(2).datatype='short';
                            else
                                disp(['TRSP error.  Aborting effort to read TRSP #' num2str(TRSPcount) '.']);
                                badFlag=1;
                                continue
                            end;
                        end;
                        theTRSP.keys(iKey).data=fieldOut;
                    end;
                    if ~strcmp(eventHdr(iEvent).value,'S255')
                        disp(['TRSP error.  Aborting effort to read TRSP #' num2str(TRSPcount) '.']);
                        badFlag=1;
                    end;
                    if badFlag
                        %                         newEventHdr(end+1:end+(iEvent-hdrStart+1))=eventHdr(hdrStart:iEvent); %add back into new eventHdr
                        badFlag=0;
                    else
                        newEventHdr(end+1).type=theTRSP.type;
                        newEventHdr(end).sample=theTRSP.sample;
                        newEventHdr(end).value=theTRSP.value;
                        newEventHdr(end).duration=theTRSP.duration;
                        newEventHdr(end).keys=theTRSP.keys;
                    end;
                    iEvent=iEvent+1;
                else
                    %not part of TRSP so include in new eventHdr
                    newEventHdr(end+1).type=eventHdr(iEvent).type;
                    newEventHdr(end).sample=eventHdr(iEvent).sample;
                    newEventHdr(end).value=eventHdr(iEvent).value;
                    newEventHdr(end).duration=eventHdr(iEvent).duration;
                    if isfield(eventHdr,'keys')
                        newEventHdr(end).keys=eventHdr(iEvent).keys;
                    end;
                    iEvent=iEvent+1;
                end;
            end;
            eventHdr=newEventHdr;
        end;
        
        %add impedance info
        if isfield(hdr,'orig') && isfield(hdr.orig,'impedances')
            if isfield(hdr.orig.impedances,'channels')
                impedances.channels=hdr.orig.impedances.channels;
            end;
            if isfield(hdr.orig.impedances,'ground')
                impedances.ground=hdr.orig.impedances.ground;
            end;
        end;
    end;
    
    %Determine number of subjects and factors
    if strcmp(fileFormat,'egi_egia') && isempty(numSubs)
        numSubs = hdr.orig.chdr(1,2); %since NetStation does not properly set the fhdr(11) field, use the number of subjects from the chdr instead
    end;
    
    if strcmp(fileFormat,'egi_sbin') && isempty(numSubs)
        if ~any(strcmp(dataType,{'single_trial','continuous'}))
            numSubs = size(theData,3)/hdr.orig.header_array(14);
            if floor(numSubs) ~= numSubs
                disp(['The number of segments (' num2str(size(theData,3)) ') is not an even multiple of the number of categories (' num2str(hdr.orig.header_array(14)) '), which suggests this is really a single trial file.']);
                dataType='single_trial';
                numSubs=1;
            end;
        else
            numSubs=1;
        end;
    end;
    
    if strcmp(dataType,'factors') && ~isempty(numSubs)
        numFacs=numSubs; %non-EP files will normally have only subjects OR factors.
        numSubs=1;
        if strcmp(fileFormat,'egi_egia')
            if numFacs == hdr.orig.chdr(1,2)-1;
                numFacs=numFacs+1; %it is assumed that factor files always have an additional summed factors waveform
            end;
        end;
    end;
    
    if isempty(numSubs)
        numSubs=1;
    end;
    
    if strcmp(fileFormat,'egi_egia')
        if ~strcmp(dataType,'factors')
            if (numSubs ~= hdr.orig.chdr(1,2)) && (hdr.orig.chdr(1,2)*hdr.orig.fhdr(18)==size(theData,3))
                numSubs = hdr.orig.chdr(1,2);
                disp(['According to the file, the number of subjects is actually: ' num2str(hdr.orig.chdr(1,2))]);
            end;
        else
            
            if (numFacs ~= hdr.orig.chdr(1,2)) && (hdr.orig.chdr(1,2)*hdr.orig.fhdr(18)==size(theData,3))
                numFacs = hdr.orig.chdr(1,2);
                disp(['According to the file, the number of factors is actually: ' num2str(hdr.orig.chdr(1,2))]);
            end;
        end;
    end;
    
    if strcmp(dataType,'single_trial') && size(theData,3) == 1 && isempty(subLabels) %if subLabels is being used, likely part of file merge
        disp('A file that has only a single ''trial'' is actually a ''continuous'' file.');
        dataType='continuous';
    end;
    
    if strcmp(dataType,'continuous')
        if numSubs > 1
            disp(['For continuous data, the number of subjects must be one, rather than ' num2str(numSubs) '.']);
        end;
        numSubs=1;
        avgNum=ones(1,1);
        covNum=ones(1,1);
        subNum=ones(1,1);
    end;
    
    if strcmp(dataType,'single_trial')
        if numSubs > 1
            disp(['For single trial data, the number of subjects must be one, rather than ' num2str(numSubs) '.']);
        end;
        numSubs=1;
        avgNum=ones(1,size(theData,3));
        covNum=ones(1,size(theData,3));
        subNum=ones(1,size(theData,3));
    end;
    
    %Determine other general information.
    if ~isempty(hdr.Fs)
        sampleRate=hdr.Fs;
    else
        sampleRate=250;
        disp('Assuming default sample rate of 250Hz.');
    end;
    
    chanNames=hdr.label;
    
    if length(chanNames) ~= length(unique(chanNames))
        msg{1}='Error: Channel names must be unique.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if isempty(baseline)
        baseline=hdr.nSamplesPre;
        if isempty(baseline)
            baseline=0;
        end;
    else
        if mod(baseline,(1000/sampleRate)) >0
            baseline=baseline-mod(baseline,(1000/sampleRate));
            disp(['The specified baseline was not evenly divided by the sample lengths of ' num2str(1000/sampleRate) ' msec.']);
            disp(['The baseline has been changed to ' num2str(baseline) ' msec.']);
        end;
        baseline=baseline/(1000/sampleRate); %change to samples
    end;
    
    numSamples=size(theData,2);
    timeNames=(((0:numSamples-1)-baseline)*(1000/sampleRate));
    
    %EGIS session file
    if strcmp(fileFormat,'egi_egis')
        trialSpecs=cell(hdr.orig.fhdr(18),(max(hdr.orig.chdr(:,5))/2));
        trialCounter=0;
        for theCell = 1:hdr.orig.fhdr(18)
            numSpecs=(hdr.orig.chdr(theCell,5)/2);
            for trial = 1:hdr.orig.chdr(theCell,2)
                trialCounter=trialCounter+1;
                cellNames{trialCounter} = deblank(hdr.orig.cnames{theCell});
                trialNames(trialCounter)=trial;
                for spec=1:numSpecs
                    trialSpecs{trialCounter,spec} = hdr.orig.chdr(theCell,5+(trial-1)*numSpecs+spec);
                end;
            end;
        end;
        numWaves=trialCounter;
        numCells=hdr.orig.fhdr(18);
        
        subNames{1}=['sub' sprintf('%03.0f',hdr.orig.fhdr(11))];
        
        for i=1:size(trialSpecs,2)
            trialSpecNames{i}=['spec' sprintf('%03.0f',i)];
        end;
        trialSpecNames{1}='edit';
        
        subjectSpecs=cell(1,13);
        subjectSpecNames{1}='RunDateMo';
        subjectSpecNames{2}='RunDateDay';
        subjectSpecNames{3}='RunDateYr';
        subjectSpecNames{4}='RunTimeHr';
        subjectSpecNames{5}='RunTimeMin';
        subjectSpecNames{6}='RunTimeSec';
        subjectSpecNames{7}='SubjID';
        subjectSpecNames{8}='Handed';
        subjectSpecNames{9}='Sex';
        subjectSpecNames{10}='Age';
        subjectSpecNames{11}='ExperID';
        subjectSpecNames{12}='Comment';
        subjectSpecNames{13}='Text';
        
        for i=1:11
            subjectSpecs{i}=num2str(hdr.orig.fhdr(i+4));
        end;
        subjectSpecs{12}=char(hdr.orig.fcom);
        subjectSpecs{13}=char(hdr.orig.ftext);
        
        %EGIS average file
    elseif strcmp(fileFormat,'egi_egia')
        numCells = hdr.orig.fhdr(18);
        numWaves=numCells; %EGIS average files are never session files
        
        cellHeaderInfoSize=size(hdr.orig.chdr,2)-numSubs*hdr.orig.chdr(1,5)/2; %normally 11 fields but in some Dien datasets there were less to save on limited header space
        
        cellNames=deblank(hdr.orig.cnames);
        if strcmp(dataType,'factors')
            for i=1:numFacs
                facNames{i,1}=['fac' sprintf('%03.0f',i)];
                facTypes{i,1}='SGL';
            end;
            facNames{numFacs,1}='summed factors'; %assume the last one is a summed factor since that is what the EP Toolkit generates by default
            facTypes{numFacs,1}='CMB';
            subNames{1}='summed subjects';
        else
            if ~any(hdr.orig.chdr(1,cellHeaderInfoSize+(1:numSubs)*hdr.orig.chdr(1,5)/2))
                for i=1:numSubs
                    subNames{i}=['sub' sprintf('%03.0f',i)]; %correct for NetStation bug of not outputing subject numbers
                end;
            else
                for i=1:numSubs
                    subNames{i}=['sub' sprintf('%03.0f',hdr.orig.chdr(1,cellHeaderInfoSize+(i-1)*hdr.orig.chdr(1,5)/2))];
                end;
            end;
        end;
        
        subjectSpecs=cell(numSubs,cellHeaderInfoSize+2);
        numSpecs=(hdr.orig.chdr(1,5)/2); %assume subject information is the same for each cell
        for theSub = 1:numSubs
            for spec=1:cellHeaderInfoSize
                subjectSpecs{theSub,spec} = num2str(hdr.orig.chdr(1,(theSub-1)*numSpecs+spec));
            end;
            subjectSpecs{theSub,cellHeaderInfoSize+1}=char(hdr.orig.fcom);
            subjectSpecs{theSub,cellHeaderInfoSize+2}=char(hdr.orig.ftext);
        end;
        
        if cellHeaderInfoSize ==11
            subjectSpecNames{1}='RunDateMo';
            subjectSpecNames{2}='RunDateDay';
            subjectSpecNames{3}='RunDateYr';
            subjectSpecNames{4}='RunTimeHr';
            subjectSpecNames{5}='RunTimeMin';
            subjectSpecNames{6}='RunTimeSec';
            subjectSpecNames{7}='SubjID';
            subjectSpecNames{8}='Handed';
            subjectSpecNames{9}='Sex';
            subjectSpecNames{10}='Age';
            subjectSpecNames{11}='ExperID';
            subjectSpecNames{12}='Comment';
            subjectSpecNames{13}='Text';
        else
            for i=1:cellHeaderInfoSize
                subjectSpecNames{i}=sprintf('spec%02.0f',i);
            end;
            subjectSpecNames{cellHeaderInfoSize+1}='Comment';
            subjectSpecNames{cellHeaderInfoSize+2}='Text';
        end;
        
        for theCell=1:numCells
            for theSub=1:numSubs
                numSpecs=(hdr.orig.chdr(theCell,5)/2);
                if strcmp(dataType,'factors') || numSubs == 1
                    subNum(theSub,theCell)=hdr.orig.chdr(theCell,5+(theSub-1)*numSpecs+1);
                    avgNum(theSub,theCell)=0;
                    covNum(theSub,theCell)=0;
                else
                    theNum=hdr.orig.chdr(theCell,5+(theSub-1)*numSpecs+1);
                    if theNum==-1
                        theNum=hdr.orig.chdr(theCell,5+(theSub-1)*numSpecs+2); %workaround for NetStation bug
                    end;
                    avgNum(theSub,theCell)=theNum;
                    covNum(theSub,theCell)=theNum;
                    subNum(theSub,theCell)=1;
                end
            end;
        end;
        
        %EGI Simple Binary
    elseif strcmp(fileFormat,'egi_sbin')
        disp('Reminder: I do not recommend using simple binary format if it can be avoided.');
        numCells = hdr.orig.header_array(14);
        if strcmp(dataType,'factors')
            for i=1:numFacs
                facNames{i}=['fac' sprintf('%03.0f',i)];
                facTypes{i}='SGL';
            end;
            facNames{numFacs}='summed factors';
            facTypes{numFacs}='CMB';
            subNames{1}='summed subjects';
            numWaves=size(theData,3)/numFacs;
        else
            for i=1:numSubs
                subNames{i}=['sub' sprintf('%03.0f',i)];
            end;
            numWaves=size(theData,3)/numSubs;
        end;
        
        if any(strcmp(dataType,{'average','grand_average'}))
            subjectSpecs=cell(numSubs,7);
            for spec=1:7
                subjectSpecs{numSubs,spec} = hdr.orig.header_array(spec+1);
            end;
            subjectSpecNames{1}='RunDateYr';
            subjectSpecNames{2}='RunDateMo';
            subjectSpecNames{3}='RunDateDay';
            subjectSpecNames{4}='RunTimeHr';
            subjectSpecNames{5}='RunTimeMin';
            subjectSpecNames{6}='RunTimeSec';
            subjectSpecNames{7}='RunTimeMsec';
        end;
        
        droppedEvent=0;
        eventIndex=ones(length(eventValues),1);
        for i=1:length(eventValues)
            if isempty(eventValues{i})
                eventIndex(i)=0;
            elseif any(strcmp(eventValues{i},{'CELL','TRSP'})) %the information associated with these two events are lost when NS exports the file.
                eventIndex(i)=0;
                droppedEvent=1;
            end;
        end;
        
        eventValues=eventValues(find(eventIndex));
        eventHdr=eventHdr(find(eventIndex));
        
        if droppedEvent
            disp('Dropping CELL and TRSP events since NetStation loses the associated information when it exports a simple binary file.');
        end;
        
        subjectsGrouped=1;
        if strcmp(dataType,'continuous')
            disp('Note that NetStation, at least as of version 4.1.2, does not import the event information associated with');
            disp('continuous simple binary files, so if you have problems with this, this is not the Toolkit''s fault.');
            numCells=1; %continuous data only has one cell
            cellNames{1}='cell01';
        else
            if ~isempty(cellLabels)
                numCells=length(unique(cellLabels));
                if length(cellLabels)==1
                    cellNames=repmat(cellLabels,numWaves,1);
                else
                    cellNames=cellLabels;
                end;
            elseif length(eventHdr)==numWaves
                numCells=length(unique(eventValues));
                cellNames=eventValues;
            elseif sum(strcmp('trial',{eventHdr.type})) ==numWaves*numSubs
                if (hdr.orig.header_array(14))==0 && (hdr.orig.header_array(15) > 1) %epoch-marked simple binary file format
                    numCells=1; %epoch-marked doesn't allow for separate cells
                    for i=1:numWaves
                        cellNames{i}='cell01';
                    end;
                    disp('Epoch-marked simple binary is not able to separate the cells so putting all the segments into the same cell.');
                    disp('You will need to use the Edit function to provide them with their correct cell names.');
                    disp('It will not be possible to use this format for combined subject average files.');
                else
                    cellNames={eventHdr(strcmp('trial',{eventHdr.type})).value}; %the category names in the SegHdrs have been stored here.
                    numCells=length(unique(cellNames));
                end;
            else
                theTrial=ceil([eventHdr.sample]/numSamples);
                theOffset=mod([eventHdr.sample],numSamples);
                if (length(find(theOffset==1))==numWaves) && length(theTrial(unique(find(theOffset==1))))==numWaves
                    numCells=length(unique(eventValues(find(theOffset==1)))); %an event is at the offset every time
                    cellNames=eventValues(find(theOffset==1));
                else
                    numCells=1; %give up trying to separate the cells and put them all into the same one.
                    for i=1:numWaves
                        cellNames{i}='cell01';
                    end;
                    disp('Can''t figure out the cell names so putting all the segments into the same cell.');
                end;
            end;
            
            uniqueCellNames=unique(cellNames);
            if any(strcmp(dataType,{'average','grand_average'}))
                subjectsGrouped=-1; %whether a subject is grouped such that all its cells are consecutive in the file
                for i=1:length(uniqueCellNames)
                    whichCells=find(strcmp(uniqueCellNames{i},cellNames));
                    if length(whichCells) ~= numSubs
                        disp('The subjects do not all have one of each cell so it is not possible to decrypt which waveform goes with which subject.');
                        return
                    end;
                    if length(whichCells) > 1
                        if all(diff(whichCells) == 1) %the waveforms are grouped by cell
                            if subjectsGrouped==-1
                                subjectsGrouped=0;
                            elseif subjectsGrouped==1
                                disp('The waveforms are not consistently organized in the file so it is not possible to decrypt which waveform goes with which subject.');
                                return
                            end;
                        elseif ~any(diff(diff(whichCells)))
                            %the waveforms are grouped by subject
                            if subjectsGrouped==-1
                                subjectsGrouped=1;
                            elseif subjectsGrouped==0
                                disp('The waveforms are not consistently organized in the file so it is not possible to decrypt which waveform goes with which subject.');
                                return
                            end;
                        else
                            disp('The waveforms are not consistently organized in the file so it is not possible to decrypt which waveform goes with which subject.');
                        end;
                    end;
                end;
                if subjectsGrouped
                    cellNames=cellNames(1:length(uniqueCellNames));
                else
                    cellNames=cellNames(1:numSubs:length(cellNames));
                end;
            else
                trialCounter=zeros(length(uniqueCellNames),1);
                for theCell=1:length(cellNames)
                    newCell=strcmp(cellNames(theCell),uniqueCellNames);
                    trialCounter(newCell)=trialCounter(newCell)+1;
                    trialNames(theCell)=trialCounter(newCell);
                end;
                trialSpecs=cell(numWaves,0);
            end;
        end;
        
        %any other sort of file other than EGIS or Simple Binary
    else
        
        if length(eventHdr)==1
            if strcmp(eventHdr.type,'average') && ~any(strcmp(dataType,{'average','grand_average'}))
                dataType='average';
                disp('The file indicates that it is actually an average file.');
            end;
        end;
        
        if strcmp(dataType,'factors')
            for i=1:numFacs
                facNames{i}=['fac' sprintf('%03.0f',i)];
                facTypes{i}='SGL';
            end;
            facNames{numFacs}='summed factors';
            facTypes{numFacs}='CMB';
            subNames{1}='summed subjects';
            numWaves=size(theData,3)/numFacs;
        else
            for i=1:numSubs
                subNames{i}=['sub' sprintf('%03.0f',i)];
            end;
            numWaves=size(theData,3)/numSubs;
        end;
        
        if strcmp(fileFormat,'ns_avg')
            if strcmp(dataType,'grand_average')
                subNum=hdr.orig.nsweeps;
            else
                avgNum=hdr.orig.nsweeps;
                covNum=hdr.orig.nsweeps;
            end;
        end;
        
        if strcmp(fileFormat,'egi_mff_v2') && any(strcmp(dataType,{'average','grand_average'}))
            if isempty(hdr.orig.epochSubjects)
                subNames{1}='sub001';
                numSubs=1;
            else
                if length(unique(hdr.orig.multiSubj)) >1
                    numSubs=length(hdr.orig.epochSubjects);
                elseif length(unique(hdr.orig.epochFilenames)) >1
                    numSubs=length(unique(hdr.orig.epochFilenames));
                else
                    disp('Unable to determine number of subjects.');
                end;
                for i=1:numSubs
                    if ~isempty(hdr.orig.epochSubjects{i})
                        subNames{i}=hdr.orig.epochSubjects{i};
                    elseif ~isempty(hdr.orig.epochFilenames{i})
                        subNames{i}=hdr.orig.epochFilenames{i};
                    else
                        subNames{i}=['sub' sprintf('%03.0f',i)];
                    end;
                end;
            end;
            numWaves=size(theData,3)/numSubs;
            cellNames=hdr.orig.epochLabels;
            cellNames=cellNames(1:numSubs:length(cellNames)); %averages are grouped by cell
            numCells=length(cellNames);
        end;
        
        if any(strcmp(fileFormat,{'egi_mff_v1';'egi_mff_v2'})) && any(strcmp(dataType,{'continuous','single_trial'}))
            %not sure yet how this structure works so may need to be revised
            for iSpec=1:length(hdr.orig.xml.subject.fields)
                subjectSpecNames{end+1}=hdr.orig.xml.subject.fields(iSpec).field.name;
                if isfield(hdr.orig.xml.subject.fields(iSpec).field.data,'data')
                    subjectSpecs{1,end+1}=hdr.orig.xml.subject.fields(iSpec).field.data.data;
                else
                    subjectSpecs{1,end+1}=[];
                end;
            end;
            %             if ~isempty(find(strcmp('SESS',eventValues))) && strcmp(dataType,'continuous')
            %                 theSess=find(strcmp('SESS',eventValues));
            %                 theSess=theSess(end); %if there are multiple SESS events, assume that the earlier ones were false starts.
            %                 for iKey=1:length(eventHdr(theSess).orig.keys)
            %                     subjectSpecNames{end+1}=eventHdr(theSess).orig.keys(iKey).key.keyCode;
            %                     subjectSpecs{1,end+1}=eventHdr(theSess).orig.keys(iKey).key.data.data;
            %                 end;
            %             end;
        end;
        
        if strcmp(fileFormat,'eeglab_set')
            if isfield(hdr.orig.event,'setname')
                if ~strcmp(dataType,'average')
                    dataType='average';
                    disp('The header indicates this is actually an average file generated by Widmann''s pop_grandaverage function.')
                end;
                numSubs=numWaves; %third dimension of data is subjects.  one cell per file.
                numWaves=1;
                numCells=1;
                if isempty(cellNames)
                    cellNames{1}=hdr.orig.event(1).setname;
                end;
                avgNum=shiftdim(avgNum);
                covNum=shiftdim(covNum);
                subNum=ones(numSubs,numCells);
                for i=1:numSubs
                    subNames{i}=['sub' sprintf('%03.0f',i)];
                    avgNum(i,1)=hdr.orig.epoch(i).eventtrials;
                    covNum(i,1)=hdr.orig.epoch(i).eventtrials;
                end;
            elseif strcmp(dataType,'single_trial') && length(hdr.orig.event) == numWaves
                recTime=cell2mat({hdr.orig.event(:).latency});
            end;
            if strcmp(dataType,'single_trial')
                if isfield(hdr.orig,'reject')
                    if isfield(hdr.orig.reject,'manualE')
                        badChans(1,:,:)=hdr.orig.reject.manualE';
                    end;
                    if isfield(hdr.orig.reject,'manual')
                        badTrials(1,:)=hdr.orig.reject.manual';
                    end;
                end;
            end;
        end;
        
        if strcmp(fileFormat,'eeglab_erp')
            if isempty(subNames)
                subNames{1}=hdr.orig.erpname; %erplab files have one subject/multiple cells whereas Widmann eeglab variant has multiple subjects/one cell
            end;
            for i=1:numWaves
                avgNum(1,i)=hdr.orig.ntrials.accepted(i);
                covNum(1,i)=hdr.orig.ntrials.accepted(i);
            end;
            if length(hdr.orig.EVENTLIST)>1 %grand average file
                subNum=ones(1,numWaves)*length(hdr.orig.EVENTLIST);
                if ~strcmp(dataType,'grand_average')
                    disp('Header indicates this is actually a grand average file.');
                    dataType='grand_average';
                end;
            end;
            avgNum(1,:)=hdr.orig.ntrials.accepted;
            covNum(1,:)=hdr.orig.ntrials.accepted;
            badTrials(1,:)=hdr.orig.ntrials.rejected;
        end;
        
        if strcmp(fileFormat,'neuromag_fif')
            if any(strcmp(dataType,{'average','grand_average','single_trial'}))
                if isfield(hdr.orig,'evoked')
                    theNames=fieldnames(hdr.orig.evoked);
                    tempVar=struct2cell(hdr.orig.evoked);
                    theAspects=cell2mat(squeeze(tempVar(find(strcmp('aspect_kind',theNames)),1,:)));
                    theEpochs=find(ismember(theAspects,[100 102 103 104]));
                    if isempty(theEpochs)
                        disp('Error: No averages in this dataset.');
                        return
                    else
                        theData=theData(:,:,theEpochs);
                        numWaves=length(theEpochs);
                        numCells=numWaves;
                        cellTypes=cellTypes(theEpochs);
                        events=events(theEpochs);
                    end;
                    if any(cell2mat(squeeze(tempVar(find(strcmp('is_smsh',theNames)),1,:))))
                        disp('Warning: Data collected using MaxShield active shielding needs to be processed with MaxFilter(tm) to produce reliable results.');
                    end;
                    avgNum=squeeze(tempVar(find(strcmp('nave',theNames)),1,:));
                    covNum=squeeze(tempVar(find(strcmp('nave',theNames)),1,:));
                    if isempty(cellNames)
                        cellNames=squeeze(tempVar(find(strcmp('comment',theNames)),1,:));
                    end;
                end;
            end;
        end;
        
        if ~strcmp(dataType,'single_trial')
            if isempty(numCells)
                numCells = size(theData,3); %number of cells is not ambiguous if not a single_trial file.
            end;
            if isempty(cellNames)
                if ~isempty(eventHdr)
                    if sum(strcmp('trial',{eventHdr.type})) ==numCells
                        cellNames={eventHdr(strcmp('trial',{eventHdr.type})).value};
                    elseif sum(strcmp('trigger',{eventHdr.type})) ==numCells
                        cellNames={eventHdr(strcmp('trigger',{eventHdr.type})).value};
                    else
                        theTrial=ceil([eventHdr.sample]/numSamples);
                        theOffset=mod([eventHdr.sample],numSamples);
                        if (length(find(theOffset==1))==numWaves) && length(theTrial(unique(find(theOffset==1))))==numWaves && ~isempty(cell2mat(eventValues))
                            numCells=length(unique(eventValues(find(theOffset==1)))); %an event is at the offset every time
                            cellNames=eventValues(find(theOffset==1));
                        else
                            for i=1:numCells
                                cellNames{i}=['cell' sprintf('%03.0f',i)];
                            end;
                        end;
                    end;
                    for i=1:length(cellNames)
                        if isempty(cellNames{i})
                            cellNames{i}=['cell' sprintf('%03.0f',i)];
                        elseif isnumeric(cellNames{i})
                            if isnan(cellNames{i})
                                cellNames{i}=['cell' sprintf('%03.0f',i)];
                            else
                                cellNames{i}=num2str(cellNames{i});
                            end;
                        end;
                    end;
                else
                    for i=1:numCells
                        cellNames{i}=['cell' sprintf('%03.0f',i)];
                    end;
                end;
            end;
            
            if length(unique(cellNames)) ~= length(cellNames)
                cellList=unique(cellNames);
                for iCell=1:length(cellList)
                    cellIndex=find(strcmp(cellList{iCell},cellNames));
                    if length(cellIndex) > 1
                        for iList=1:length(cellIndex)
                            cellNames{cellIndex(iList)}=[cellNames{cellIndex(iList)} '-' sprintf('%03d',iList)];
                        end;
                    end;
                end;
                disp('EP Toolkit requires each average condition to have a unique name.  The names are being modified to make them unique.');
            elseif numSubs==1
                if all(strncmp(cellNames,'Sub',3)) %if all the cellNames start with 'Sub' then may be a combined subject average file
                    tempVar=char(cellNames);
                    subNums=str2double(tempVar(:,4:6));
                    tempCellNames=cellstr(deblank(tempVar(:,7:end)));
                    if all(subNums) %all of them must also have three digit numbers after the Sub
                        if all(hist(subNums,unique(subNums))==length(unique(tempCellNames))) %if the putative cells and subject names are fully crossed
                            cellNames=unique(tempCellNames);
                            tempSubNames=cellstr(num2str(subNums));
                            subNames=unique(tempSubNames);
                            numSubs=length(subNames);
                            numCells=length(cellNames);
                            numWaves=numCells;
                            theData2=zeros(length(chanNames),numSamples,numCells,numSubs);
                            for i=1:numCells
                                theData2(:,:,strcmp(tempCellNames{i},cellNames),strcmp(tempSubNames{i},subNames))=theData(:,:,i);
                            end;
                            for i=1:numSubs
                                subNames{i}=['sub' strtrim(subNames{i})];
                            end;
                        end;
                    end;
                end;
            end;
        end;
        
        if strcmp(dataType,'single_trial')
            if isempty(cellNames)
                if ~isempty(cellLabels)
                    if length(cellLabels)==1
                        cellNames=repmat(cellLabels,numWaves,1);
                    else
                        cellNames=cellLabels;
                    end;
                elseif strcmp(fileFormat,'ns_eeg')
                    for i=1:numWaves
                        cellNames{i}='cell001';
                    end;
                elseif length(eventHdr)==1 && strcmp(eventHdr.type,'average')
                    for i=1:numWaves
                        cellNames{i}='cell001';
                    end;
                elseif length(find(~strcmp(eventValues,'trial')))==numWaves
                    cellNames=eventValues(find(~strcmp(eventValues,'trial')))';
                elseif length(find(strcmp({eventHdr.type},'trial')))==numWaves
                    cellNames=eventValues(find(strcmp({eventHdr.type},'trial')))';
                    if strcmp(fileFormat,'egi_mff_v1')
                        recTime=[hdr.orig.epochdef(:,3)];
                    end;
                else
                    if ~isempty(eventHdr)
                        theTrial=ceil([eventHdr.sample]/numSamples);
                        theOffset=mod([eventHdr.sample],numSamples);
                        if (length(find(theOffset==1))==numWaves) && length(theTrial(unique(find(theOffset==1))))==numWaves && all(cellfun(@isempty,eventValues))
                            cellNames=eventValues(find(theOffset==1));
                        end;
                        n = hist(theOffset,length(unique(theOffset)));
                        if length(find(n==numWaves))==1 %if there is one and only one sample for which the number of events equals the number of trials
                            uniqueEventSamples=unique(theOffset);
                            for i=1:length(uniqueEventSamples)
                                whichEvents=find(theOffset==uniqueEventSamples(i));
                                if length(whichEvents) == numWaves
                                    if strcmp(eventHdr(whichEvents(1)).type,'trigger')
                                        if ~isempty(eventHdr(whichEvents(1)).value)
                                            cellNames=eventValues(whichEvents);
                                        end;
                                        
                                    end;
                                    continue %no need to go through entire loop since found the one
                                end;
                            end;
                        end;
                    end;
                end;
            end;
            
            if isempty(cellNames) %if was unable to determine cell names assume just one cell
                disp('Was unable to identify distinct cell names from the event header so assuming all just one condition.');
                for i=1:numWaves
                    cellNames{i}='cell01';
                end;
            end;
            
            for i=1:length(cellNames)
                if isnumeric(cellNames{i})
                    cellNames{i}=num2str(cellNames{i});
                end;
            end;
            
            uniqueCellNames=unique(cellNames);
            numCells=length(uniqueCellNames);
            trialCounter=zeros(length(uniqueCellNames),1);
            for theCell=1:length(cellNames)
                newCell=strmatch(cellNames(theCell),uniqueCellNames,'exact');
                trialCounter(newCell)=trialCounter(newCell)+1;
                trialNames(theCell)=trialCounter(newCell);
            end;
            trialSpecs=cell(numWaves,0);
        end
        
    end;
    
    if ~isempty(freqNames) && length(timeNames) == 1
        timeNames=[]; %may assume that data is spectral data so no time points
    end;
    
    if isempty(subNum)
        if strcmp(dataType,'factors')
            subNum=zeros(numSubs,numCells); %number of subjects is unknown
        else
            subNum=ones(numSubs,numCells);
        end;
    end;
    
    if isempty(avgNum)
        if any(strcmp(dataType,{'continuous','single_trial'}))
            avgNum=ones(numSubs,numCells);
        else
            avgNum=zeros(numSubs,numCells); %number in average is unknown
        end;
    end;
    
    if isempty(covNum)
        if any(strcmp(dataType,{'continuous','single_trial'}))
            covNum=ones(numSubs,numCells);
        else
            covNum=zeros(numSubs,numCells); %number in average is unknown
        end;
    end;
    
    if (mod(size(theData,3),numSubs) ~= 0) && strcmp(dataType,'average')
        msg{1}=['Number of subjects (' num2str(numSubs) ') does not divide evenly into the number of epochs (' num2str(size(theData,3)) ').'];
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if strcmp(dataType,'factors')
        if (mod(size(theData,3),numFacs) ~= 0)
            msg{1}=['Number of factors (' num2str(numFacs) '), including summed factor, does not divide evenly into the number of epochs (' num2str(size(theData,3)) ').'];
            [msg]=ep_errorMsg(msg);
            return
        end
    end;
    
    if isempty(ename)
        if any(strcmp(fileFormat,{'egi_egia','egi_egis','factors'}))
            ename = hdr.orig.ename;
        else
            ename = [];
        end;
    end;
    
    if size(ename,1) > size(ename,2)
        ename=ename';
        if size(ename,1) >1
            ename=[];
        end;
    end;
    
    if ~isempty(ename)
        ename=strtrim(ename);
    end;
    
    if sum(subjects)==0
        subjects=[];
    end
    
    if isempty(recTime)
        recTime=[1:sampleRate:sampleRate*(numWaves-1)+1]';
    end;
    
    if ~isempty(cells) && (max(cells) > numCells)
        msg{1}=['Largest cell to be kept (' num2str(max(cells)) ') is larger than number of cells (' num2str(numCells) ') in dataset.'];
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if ~strcmp(dataType,'factors')
        if (numWaves*numSubs ~= size(theData,3))
            msg{1}=['Number of epochs (' num2str(numWaves) ') times subjects (' num2str(numSubs) ') does not equal number of observations (' num2str(size(theData,3)) ').'];
            [msg]=ep_errorMsg(msg);
            return
        end;
    elseif (numWaves*numFacs ~= size(theData,3))
        msg{1}=['Number of epochs (' num2str(numWaves) ') times factors (' num2str(numFacs) ') does not equal number of observations (' num2str(size(theData,3)) ').'];
        [msg]=ep_errorMsg(msg);
        return
    end;
    
    if length(trialNames)==1
        trialStr='trial';
    else
        trialStr='trials';
    end
    if numCells==1
        cellStr='cell';
    else
        cellStr='cells';
    end
    if strcmp(dataType,'factors')
        if numFacs==1
            subStr='factor';
        else
            subStr='factors';
        end
    else
        if numSubs==1
            subStr='subject';
        else
            subStr='subjects';
        end
    end;
    
    if strcmp(dataType,'single_trial')
        disp(['Read in ' num2str(length(unique(cellNames))) ' ' cellStr ' with a total of ' num2str(length(trialNames)) ' ' trialStr '.']);
    elseif  strcmp(dataType,'factors')
        disp(['Read in ' num2str(numCells) ' ' cellStr ' and ' num2str(numFacs) ' ' subStr '.']);
        disp('Making the assumption that the last factor is the combination of the other factors.');
    else
        disp(['Read in ' num2str(numCells) ' ' cellStr ' and ' num2str(numSubs) ' ' subStr '.']);
    end;
    
    if strcmp(dataType,{'continuous'}) && ~isempty(cells)
        disp('Cannot select subset of cells in a continuous raw data file.');
        cells=[];
    end;
    
    if ~isempty(samples)
        if max(samples) > size(theData,2)
            msg{1}='Largest sample is larger than number of timepoints in data.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        if min(samples) < 0
            error('No negative samples.');
        end;
        
    end;
    
    if ~isempty(channels)
        if max(channels) > size(theData,1)
            msg{1}='Largest channel is larger than number of channels in data.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if min(channels) < 0
            error('No negative channels.');
        end;
    end;
    
    if ~isempty(subjects)
        if min(subjects) < 0
            error('No negative subjects.');
        end;
        if ~isempty(subjects) && (max(subjects) > numSubs)
            msg{1}=['Largest subject to be kept (' num2str(max(subjects)) ') is larger than number of subjects (' num2str(numSubs) ') indicated for dataset.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
    end;
    
    %Reorganize waveforms and events and extract events from EventHdr
    %mff is file time whereas everything else from FieldTrip is real time
    
    if ~isempty(eventHdr)
        if ~strcmp(dataType,'factors')
            events=cell(numSubs,numWaves);
            samples_trials = [eventHdr(find(strcmp('trial', {eventHdr.type}))).sample]; %find "trial" events
            if any(strcmp('TRSP',{eventHdr.value})) && isfield(eventHdr,'orig') && ~strcmp(dataType,'continuous')
                for i=1:length(eventHdr(min(find(strcmp('TRSP',{eventHdr.value})))).orig.keys) %assume TRSP keys same for all trials
                    trialSpecNames{end+1}=eventHdr(min(find(strcmp('TRSP',{eventHdr.value})))).orig.keys(i).code;
                end;
            end;
            for theEvent=1:length(eventHdr)
                if ~isempty(eventHdr(theEvent).sample) && ~isnan(eventHdr(theEvent).sample)
                    if isempty(samples_trials)
                        epoch=floor((eventHdr(theEvent).sample-1)/numSamples)+1;
                        sample=mod((eventHdr(theEvent).sample-1),numSamples)+1;
                    else
                        epoch=max(find(samples_trials<=eventHdr(theEvent).sample));
                        sample=eventHdr(theEvent).sample-samples_trials(epoch)+1;
                    end;
                    if strcmp(fileFormat,'egi_sbin') && subjectsGrouped
                        %cells are grouped together by subject.
                        subject=floor((epoch-1)/numWaves)+1;
                        cellTrial=mod((epoch-1),numWaves)+1;
                    else
                        %assume that subjects are grouped together by cell.
                        cellTrial=floor((epoch-1)/numSubs)+1;
                        subject=mod((epoch-1),numSubs)+1;
                    end;
                    if (epoch > numWaves*numSubs) || epoch < 1
                        disp('Warning: event falls outside recorded data.  Discarded.');
                    else
                        if strcmp(eventHdr(theEvent).value,'TRSP') && strcmp(dataType,'single_trial') && isfield(eventHdr,'orig')
                            for trsp=1:length(eventHdr(theEvent).orig.keys)
                                trialSpecs{cellTrial,find(strcmp(eventHdr(theEvent).orig.keys(trsp).code,trialSpecNames))}=eventHdr(theEvent).orig.keys(trsp).data;
                            end;
                        elseif ~isempty(trialSpecNames)
                            oneEvent=[];
                            oneEvent.type=eventHdr(theEvent).type;
                            oneEvent.sample=sample;
                            oneEvent.value=eventHdr(theEvent).value;
                            oneEvent.duration=eventHdr(theEvent).duration;
                            if isfield(eventHdr,'keys')
                                oneEvent.keys=eventHdr(theEvent).keys;
                            else
                                oneEvent.keys=struct('code','','data','','datatype','','description','');
                            end;
                            events{subject,cellTrial}=[events{subject,cellTrial} oneEvent];
                        elseif any(strcmp(fileFormat,{'eeglab_set','eeglab_erp'}))
                            oneEvent=[];
                            oneEvent.type=eventHdr(theEvent).type;
                            oneEvent.sample=sample;
                            oneEvent.value=eventHdr(theEvent).value;
                            oneEvent.duration=eventHdr(theEvent).duration;
                            oneEvent.keys=struct('code','','data','','datatype','','description','');
                            nameList=fieldnames(eventHdr(theEvent));
                            nameList=setdiff(nameList,{'type','value','sample','offset','duration'});
                            for iKey=1:length(nameList)
                                eval(['fieldEmpty=isempty(eventHdr(theEvent).' nameList{iKey} ');']);
                                if ~fieldEmpty
                                    fieldLabel=nameList{iKey};
                                    if ~isempty(strfind(fieldLabel,'hash_'))
                                        fieldLabel=strrep(fieldLabel,'hash_','#');
                                    end;
                                    if ~isempty(strfind(fieldLabel,'plus_'))
                                        fieldLabel=strrep(fieldLabel,'plus_','+');
                                    end;
                                    eval(['oneEvent.keys(end+1).code=''' fieldLabel ''';']);
                                    eval(['oneEvent.keys(end).data=eventHdr(theEvent).' nameList{iKey} ';']);
                                end;
                            end;
                            events{subject,cellTrial}=[events{subject,cellTrial} oneEvent];  
                        else
                            oneEvent=[];
                            oneEvent.type=eventHdr(theEvent).type;
                            oneEvent.sample=sample;
                            oneEvent.value=eventHdr(theEvent).value;
                            oneEvent.duration=eventHdr(theEvent).duration;
                            if isempty(oneEvent.value)
                                oneEvent.value=oneEvent.type;
                                oneEvent.type='trigger';
                            end;
                            if strcmp(fileFormat,'egi_mff_v2')
                                if strcmp(oneEvent.type,'break cnt') || strcmp(oneEvent.value,'break cnt')
                                    oneEvent.type='boundary';
                                    oneEvent.value='boundary';
                                end;
                            end;
                            if strcmp(fileFormat,'biosemi_bdf')
                                if strcmp(oneEvent.type,'trigger') && strcmp(oneEvent.value,'Epoch')
                                    oneEvent.type='boundary';
                                    oneEvent.value='boundary';
                                end;
                            end;
                            if isfield(eventHdr,'orig')
                                if isfield(eventHdr(theEvent).orig,'keys')
                                    if strcmp(fileFormat,'egi_mff_v1')
                                        if isfield(eventHdr(theEvent).orig.keys,'key')
                                            for iKey=1:length(eventHdr(theEvent).orig.keys)
                                                oneEvent.keys(iKey).code=eventHdr(theEvent).orig.keys(iKey).key.keyCode;
                                                oneEvent.keys(iKey).data=eventHdr(theEvent).orig.keys(iKey).key.data.data;
                                                oneEvent.keys(iKey).datatype=eventHdr(theEvent).orig.keys(iKey).key.data.dataType;
                                                oneEvent.keys(iKey).description='';
                                            end;
                                        else
                                            oneEvent.keys=struct('code','','data','','datatype','','description','');
                                        end;
                                    else
                                        oneEvent.keys=cell2mat(eventHdr(theEvent).orig.keys);
                                    end;
                                else
                                    oneEvent.keys=struct('code','','data','','datatype','','description','');
                                end;
                            else
                                oneEvent.keys=struct('code','','data','','datatype','','description','');
                            end;
                            events{subject,cellTrial}=[events{subject,cellTrial} oneEvent];
                        end;
                    end;
                end;
            end;
        else
            events=cell(numSubs,length(cellNames));
            for theEvent=1:length(eventHdr)/numFacs
                if ~isempty(eventHdr(theEvent).sample)
                    epoch=floor((eventHdr(theEvent).sample-1)/numSamples)+1;
                    %assume that the events for all the factors are duplicates of the first as all from same set of subject(s).
                    cellTrial=epoch;
                    subject=1;
                    oneEvent=eventHdr(theEvent);
                    events{subject,cellTrial}=[events{subject,cellTrial} oneEvent];
                end;
            end;
        end;
    else
        events=cell(numSubs,length(cellNames));
    end;
    
    %reorganize data into 6D array
    if ~strcmp(dataType,'factors')
        sevenDdata=zeros(size(theData,1),size(theData,2),numWaves,numSubs,1,1,1);
        if strcmp(fileFormat,'egi_sbin') && subjectsGrouped
            %cells are grouped together by subject.
            for sub=1:numSubs
                sevenDdata(:,:,:,sub,1,1,1)=theData(:,:,1+(sub-1)*numWaves:sub*numWaves);
            end;
        else
            %assuming that subjects are grouped together by cell.
            for sub=1:numSubs
                sevenDdata(:,:,:,sub,1,1,1)=theData(:,:,sub:numSubs:numSubs*numWaves);
            end;
        end;
    else
        sevenDdata=zeros(size(theData,1),size(theData,2),numWaves,numSubs,numFacs,1,1);
        for fac=1:numFacs
            sevenDdata(:,:,:,:,fac,1,1)=theData(:,:,fac:numFacs:numFacs*numWaves);
        end;
    end;
end;

numChan=length(chanNames);

% apply cell and sub names if separately specified
if ~isempty(cellLabels)
    if length(cellLabels)==1
        if ~strcmp(dataType,'average') || (length(cellNames)==1)
            cellNames=repmat(cellLabels,numWaves,1);
        end;
    elseif numWaves == length(cellLabels)
        cellNames=cellLabels;
    else
        disp(['The number of epochs (' num2str(numWaves) ') did not match the number of condition labels (' num2str(length(cellLabels)) ').']);
    end;
end;

if ~isempty(subLabels)
    if length(subLabels) ==1
        subNames{1}=subLabels{1};
    elseif numSubs == length(subLabels)
        subNames=subLabels;
    else
        disp(['The number of subjects (' num2str(numSubs) ') did not match the number of subject labels (' num2str(length(subLabels)) ').']);
    end;
end;

numEpochs=numWaves;
if strcmp(dataType,'continuous')
    numEpochs=floor(size(sevenDdata,2)/sampleRate); %excess time points are tacked onto final epoch
    if numEpochs == 0
        numEpochs =1;
    end;
end;

%add analysis fields
if isempty(blinkTrial)
    blinkTrial=zeros(numSubs,numEpochs);
end;
if isempty(saccadeTrial)
    saccadeTrial=zeros(numSubs,numEpochs);
end;
if isempty(saccadeOnset)
    saccadeOnset=zeros(numSubs,numEpochs);
end;
if isempty(moveTrial)
    moveTrial=zeros(numSubs,numEpochs);
end;
if isempty(badTrials)
    badTrials=zeros(numSubs,numEpochs);
end;
if isempty(badChans)
    badChans=zeros(numSubs,numEpochs,numChan);
end;

%add Type information

if isempty(chanTypes)
    chanTypes=cellstr(repmat('EEG',numChan,1));
end;

if isempty(cellTypes)
    cellTypes=cellstr(repmat('SGL',numWaves,1));
end;

if isempty(subTypes)
    switch dataType
        case {'single_trial', 'continuous'}
            subTypes=cellstr(repmat('RAW',numSubs,1));
        case 'average'
            subTypes=cellstr(repmat('AVG',numSubs,1));
        case 'grand_average'
            subTypes=cellstr(repmat('GAV',numSubs,1));
        case 'factors'
            subTypes=cellstr(repmat('GAV',numSubs,1));
    end;
end;

if isempty(timeUnits)
    timeUnits='ms';
end;

if ~strcmp(fileFormat,'ep_mat')
    %if an EEG channel has any NaN values, as from EEGlab automatic editing, it
    %will be set to bad data.
    if any(any(any(any(any(any(any(isnan(sevenDdata))))))))
        for iSub=1:size(sevenDdata,4)
            for iCell=1:size(sevenDdata,3)
                for iChan=1:size(sevenDdata,1)
                    if any(strcmp(chanTypes{iChan},{'EEG','REG'}))
                        if any(any(any(any(isnan(squeeze(sevenDdata(iChan,:,iCell,iSub)))))))
                            sevenDdata(iChan,:,iCell,iSub,:,:,:)=0;
                            if any(strcmp(dataType,{'continuous','single_trial'}))
                                badChans(iSub,iCell,iChan)=-1;
                            elseif strcmp(dataType,'average')
                                badChans(iSub,iCell,iChan)=NaN;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

%channel coordinate information
try
    if any(strcmp(fileFormat,{'egi_egia','egi_egis'})) && ismac && isempty(montage)
        [err, montage]=ep_identifySLAY(fileName);
        if ~isempty(montage)
            if err == 1
                disp('Unable to identify montage information in EGIS file because:');
                disp(montage);
                montage=[];
            end;
        end;
    end;
catch
end;

% if strcmp(fileFormat,'egi_mff') && isempty(montage)
%     montage=hdr.orig.xml.sensorLayout.name;
% end;

if isempty(montage) && any(strcmp(fileFormat,{'egi_egia','egi_egis'}))
    montage=ep_askForMontage;
end;

if any(strcmp(fileFormat,{'eeglab_set','eeglab_erp'}))
    if isfield(hdr,'orig')
        if isfield(hdr.orig,'chanlocs')
            if isfield(hdr.orig,'chaninfo')
                if isfield(hdr.orig.chaninfo,'filename')
                    [pathstr, name, ext] = fileparts(hdr.orig.chaninfo.filename);
                    ced=[name ext];
                else
                    ced='internal';
                end;
            else
                ced='internal';
            end;
            eloc=hdr.orig.chanlocs;
            if isempty(eloc)
                ced=[];
            end;
        end;
    end;
end;

if strcmp(fileFormat,'neuromag_fif')
    %the .chs coordinates will always be Neuromag, meaning in meters, with +Y through the nose and +X through the right pre-auricular fiducial.
    %but since the .chs coordinates don't include the fiducials, which are important, the .dig coordinates will be used
    %instead unless they are not available to keep all the coordinates consistent with each other.
    %since the .dig coordinates cannot be assumed to have any specific coordinate system or orientation, the head rotation preferences
    %will be applied.
    
    if isempty(eloc) && isempty(ced)
        implicitCoords=[];
        if isfield(hdr.orig,'dig')
            if isfield(hdr.orig.dig,'kind')
                FIDchans=find([hdr.orig.dig.kind]==1);
                if length(FIDchans)==3
                    if isempty(implicit)
                        implicitCoords(1).labels='LPA';
                        implicitCoords(1).X=hdr.orig.dig(FIDchans(1)).r(1);
                        implicitCoords(1).Y=hdr.orig.dig(FIDchans(1)).r(2);
                        implicitCoords(1).Z=hdr.orig.dig(FIDchans(1)).r(3);
                        implicitCoords(2).labels='nasion';
                        implicitCoords(2).X=hdr.orig.dig(FIDchans(2)).r(1);
                        implicitCoords(2).Y=hdr.orig.dig(FIDchans(2)).r(2);
                        implicitCoords(2).Z=hdr.orig.dig(FIDchans(2)).r(3);
                        implicitCoords(3).labels='RPA';
                        implicitCoords(3).X=hdr.orig.dig(FIDchans(3)).r(1);
                        implicitCoords(3).Y=hdr.orig.dig(FIDchans(3)).r(2);
                        implicitCoords(3).Z=hdr.orig.dig(FIDchans(3)).r(3);
                    end;
                end;
                %I don't understand MEG well enough yet to implement reading its coordinate information so will be set to [0 0 0].
                EEGchans=find(strcmp(chanTypes,'EEG'));
                badDigFlag=0;
                if (length(EEGchans)+1) == length(find([hdr.orig.dig.kind]==3))
                    %there is an implicit reference so add it to the data as an explicit channel.
                    theRef=intersect(find([hdr.orig.dig.ident]==0),find([hdr.orig.dig.kind]==3));
                    addRefFlag=0;
                    if length(theRef)==1
                        if ~isempty(facVecS)
                            facVecS(end+1,:)=zeros(size(facVecS,2),1);
                        end;
                        sevenDdata(end+1,:,:,:,:,:,:)=0;
                        if ~isempty(noise)
                            noise(end+1,:,:,:,:)=zeros(1,size(noise,2),size(noise,3),size(noise,4),size(noise,5));
                        end;
                        if ~isempty(stanDev)
                            stanDev(end+1,:,:,:,:,:)=zeros(1,size(stanDev,2),size(stanDev,3),size(stanDev,4),size(stanDev,5),size(stanDev,6));
                        end;
                        badChans(:,:,end+1)=zeros(size(badChans,1),size(badChans,2),1); %implicit reference channels are not bad
                        chanNames{end+1}='REF';
                        chanTypes{end+1}='EEG';
                        reference.original=length(chanNames);
                        reference.current=reference.original;
                        addRefFlag=1;
                    else
                        disp('Was unable to determine which channel was the implicit reference.');
                        badDigFlag=1;
                    end;
                elseif length(EEGchans) ~= length(find([hdr.orig.dig.kind]==3))
                    badDigFlag=1;
                end;
                if ~badDigFlag
                    digEEGchans=find([hdr.orig.dig.kind]==3);
                    tempEloc=struct('labels',chanNames);
                    for iChan=1:length(EEGchans)
                        theChan=EEGchans(iChan);
                        tempEloc(theChan).X=hdr.orig.dig(digEEGchans(iChan)).r(1);
                        tempEloc(theChan).Y=hdr.orig.dig(digEEGchans(iChan)).r(2);
                        tempEloc(theChan).Z=hdr.orig.dig(digEEGchans(iChan)).r(3);
                    end;
                    if addRefFlag
                        tempEloc(end).X=hdr.orig.dig(theRef).r(1);
                        tempEloc(end).Y=hdr.orig.dig(theRef).r(2);
                        tempEloc(end).Z=hdr.orig.dig(theRef).r(3);
                    end;
                    %rotate head so it is pointing upwards
                    eloc=tempEloc;
                    for iChan=1:length(eloc)
                        switch elecPrefs
                            case {1,5}
                                eloc(iChan).X=tempEloc(iChan).X;
                                eloc(iChan).Y=tempEloc(iChan).Y;
                            case {2,6}
                                eloc(iChan).X=tempEloc(iChan).Y;
                                eloc(iChan).Y=tempEloc(iChan).X;
                            case {3,7}
                                eloc(iChan).X=-tempEloc(iChan).X;
                                eloc(iChan).Y=-tempEloc(iChan).Y;
                            case {4,8}
                                eloc(iChan).X=-tempEloc(iChan).Y;
                                eloc(iChan).Y=-tempEloc(iChan).X;
                        end;
                        if elecPrefs > 4
                            eloc(iChan).Y=-eloc(iChan).Y;
                        end;
                    end;
                end;
            end;
        end;
        if isempty(eloc) && isfield(hdr.orig,'chs') %if it didn't work out with .dig, try .chs for electrode coordinates
            if isfield(hdr.orig.chs,'eeg_loc')
                tempEloc=struct('labels',hdr.label);
                refCoord=[];
                for iChan=1:length(hdr.label)
                    tempVar=hdr.orig.chs(iChan).eeg_loc;
                    if ~isempty(tempVar)
                        tempEloc(iChan).X=tempVar(1,1);
                        tempEloc(iChan).Y=tempVar(2,1);
                        tempEloc(iChan).Z=tempVar(3,1);
                        if isempty(refCoord) && (size(tempVar,2) == 2)
                            refCoord.X=tempVar(1,2);
                            refCoord.Y=tempVar(2,2);
                            refCoord.Z=tempVar(3,2);
                        end;
                    end;
                end;
                
                %add the reference channel, which is implicit if coordinates are given.  Assume a common reference site.
                if ~isempty(refCoord) && any([refCoord.X refCoord.Y refCoord.Z]) %if no reference, then the coords will be [0,0,0] in FIFF files
                    M1=length(hdr.label)+1;
                    reference.original=M1;
                    reference.current=reference.original; %assume that if there are implicit references then they were the original reference and still are
                    tempEloc(M1).X=refCoord.X;
                    tempEloc(M1).Y=refCoord.Y;
                    tempEloc(M1).Z=refCoord.Z;
                    tempEloc(M1).labels='REF';
                    chanTypes{M1}='EEG';
                    sevenDdata(end+1,:,:,:,:,:,:)=0;
                end;
            end;
        end;
        if ~isempty(eloc)
            ced='internal';
            try
                eloc = convertlocs(eloc, 'cart2all');
                [eloc.type] = chanTypes{:};
            catch
                disp('Warning: Reading of electrode coordinates from file has failed.');
                eloc=[];
            end;
            if isempty(implicit) && ~isempty(implicitCoords)
                try
                    implicit = convertlocs(implicitCoords, 'cart2all');
                    [implicit.type]=deal('FID');
                catch
                    disp('Warning: Reading of fiducial coordinates from file has failed.');
                end;
            end;
        end;
        %convert to microvolts
        for iChan=1:length(hdr.orig.chs)
            if hdr.orig.chs(iChan).kind==2 %EEG channels
                sevenDdata(iChan,:,:,:,:,:,:)=sevenDdata(iChan,:,:,:,:,:,:)*10^double(6+hdr.orig.chs(iChan).unit_mul);
            end;
        end;
    end;
    %not sure if cov matrices are normally included in the fif average file but just in case.
    try
        [tempVar] = mne_read_noise_cov(fileName);
        tempVar=tempVar*10^12; %convert to microvolts
        covMatrix(1,:,:)=tempVar;
        covNq=NaN;
    catch
        %no cov file
    end;
    for iBad=1:length(hdr.orig.bads)
        if any(strcmp(dataType,{'average','grand_average'}))
            badChans(:,:,find(strcmp(hdr.orig.bads{iBad},chanNames)))=NaN;
        else
            badChans(:,:,find(strcmp(hdr.orig.bads{iBad},chanNames)))=-1;
        end;
    end;
end;

if strcmp(fileFormat,'egi_mff_v1') || (strcmp(fileFormat,'egi_mff_v2') && ~any(strcmp(dataType,{'average','grand_average'})))
    %electrode coordinates included in mff assumed to be Polhemus Cartesian coordinates, as in .elp files.
    %Their absolute values should be equal to or less than one.
%    if isempty(eloc) && (isempty(ced) || strcmp(ced,'internal'))
        if isfield(hdr.orig.xml,'coordinates')
            if isfield(hdr.orig.xml.coordinates,'sensorLayout')
                tempEloc=struct('labels',hdr.label);
                tempImplicit=struct('labels',cell(0));
                for iChan=1:length(hdr.orig.xml.coordinates.sensorLayout.sensors)
                    switch hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.type
                        case '0' %EEG channel
                            tempEloc(iChan).X=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.x);
                            tempEloc(iChan).Y=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.y);
                            tempEloc(iChan).Z=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.z);
                        case '1' %reference channel
                            tempEloc(iChan).X=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.x);
                            tempEloc(iChan).Y=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.y);
                            tempEloc(iChan).Z=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.z);
                            if isempty(reference.original)
                                reference.original=iChan;
                            end;
                        case '2' %fiducial
                            if isempty(implicit)
                                tempImplicit(end+1).labels=hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.name;
                                tempImplicit(end).X=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.x);
                                tempImplicit(end).Y=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.y);
                                tempImplicit(end).Z=str2double(hdr.orig.xml.coordinates.sensorLayout.sensors(iChan).sensor.z);
                                tempImplicit(end).type='FID';
                            end;
                    end;
                end;
                %rotate head so it is pointing upwards
                eloc=tempEloc;
                implicitCoords=tempImplicit;
                if ~isempty(elecPrefs)
                    for iChan=1:length(eloc)
                        switch elecPrefs
                            case {1,5}
                                eloc(iChan).X=tempEloc(iChan).X;
                                eloc(iChan).Y=tempEloc(iChan).Y;
                            case {2,6}
                                eloc(iChan).X=tempEloc(iChan).Y;
                                eloc(iChan).Y=tempEloc(iChan).X;
                            case {3,7}
                                eloc(iChan).X=-tempEloc(iChan).X;
                                eloc(iChan).Y=-tempEloc(iChan).Y;
                            case {4,8}
                                eloc(iChan).X=-tempEloc(iChan).Y;
                                eloc(iChan).Y=-tempEloc(iChan).X;
                        end;
                        if elecPrefs > 4
                            eloc(iChan).Y=-eloc(iChan).Y;
                        end;
                    end;
                    for iChan=1:length(implicitCoords)
                        switch elecPrefs
                            case {1,5}
                                implicitCoords(iChan).X=tempImplicit(iChan).X;
                                implicitCoords(iChan).Y=tempImplicit(iChan).Y;
                            case {2,6}
                                implicitCoords(iChan).X=tempImplicit(iChan).Y;
                                implicitCoords(iChan).Y=tempImplicit(iChan).X;
                            case {3,7}
                                implicitCoords(iChan).X=-tempImplicit(iChan).X;
                                implicitCoords(iChan).Y=-tempImplicit(iChan).Y;
                            case {4,8}
                                implicitCoords(iChan).X=-tempImplicit(iChan).Y;
                                implicitCoords(iChan).Y=-tempImplicit(iChan).X;
                        end;
                        if elecPrefs > 4
                            implicitCoords(iChan).Y=-tempImplicit(iChan).Y;
                        end;
                    end;
                end;
                ced='internal';
                try
                    eloc = convertlocs(eloc, 'cart2all');
                    [eloc.type]=deal('EEG');
                catch
                    disp('Warning: Reading of electrode coordinates from file has failed.');
                    eloc=[];
                end;
                if isempty(implicit) && ~isempty(implicitCoords)
                    try
                        implicit = convertlocs(implicitCoords, 'cart2all');
                        [implicit.type]=deal('FID');
                    catch
                        disp('Warning: Reading of fiducial coordinates from file has failed.');
                    end;
                end;
            end;
        end;
%    end;
end;

if any(strcmp(fileFormat,{'egi_mff_v1','egi_mff_v2'}))
    for iChan=1:length(chanTypes)
        switch hdr.chantype{iChan}
            case 'eeg'
                chanTypes{iChan}='EEG';
            case {'ecg','ECG'}
                chanTypes{iChan}='ECG';
            otherwise
                chanTypes{iChan}='ANS'; %default to ANS if unknown or unrecognized
        end;
        if any(strcmp(hdr.label{iChan},{'REF','VREF','Cz','vertex reference'}))
            chanTypes{iChan}='EEG';
        end;
    end;
end;


%test out the ced file and set to null if it's no good.
if ~isempty(ced) && ~any(strcmp(ced,{'none','internal'})) && isempty(eloc)
    if exist(ced,'file')
        whichCED=ced;
        [pathstr, name, fileSuffix] = fileparts(whichCED);
        ced=[name fileSuffix];
    else
        whichCED=which(ced);
    end;
    try
        evalc('eloc = readlocs([whichCED],''filetype'',''chanedit'');');
        implicit=eloc(1); %setup up implicit to have the same structure as eloc.
        implicit(1)=[];
        eloc=[];
    catch
        disp(['The ced file ' ced ' did not work for some reason.  The error message was:']);
        disp(lasterr)
        ced = [];
        eloc = [];
        implicit=[];
    end;
end;

if any(strcmp(ced,{'none','internal'})) && isempty(eloc)
    ced=[];
end;

if isempty(ced)
    if any(strcmp(fileFormat,{'egi_egia','egi_egis'}))
        [ced]=ep_whichEGIced(montage);
    end;
    whichCED=ced;
    if isempty(ced)
        thisFile = mfilename('fullpath');
        [pathstr, name, ext] = fileparts(thisFile);
        theSeps=findstr(pathstr,filesep);
        pathstr=[pathstr(1:theSeps(end)-1) filesep 'electrodes'];
        [ced, pathstr] = uigetfile('*.ced',['Electrode Coordinate file (' num2str(numChan) ' channels):'],pathstr);
        whichCED=[pathstr ced];
    end;
    if isnumeric(ced)
        ced='none';
        whichCED=ced;
    end;
else
    whichCED=ced;
end;

if (~isempty(ced) && ~any(strcmp(ced,{'none','internal'})) && ~strcmp(fileFormat,'ep_mat')) || (~isempty(ced) && (isempty(eloc) || all(isempty([eloc.radius])) || all(isempty([eloc.theta]))) && ~strcmp('none',ced))
    if ~exist('hdr','var')
        hdr=[];
    end;
    if ~isempty(eloc)
        if  all(isempty([eloc.theta])) && all(~isempty([eloc.radius]))
            disp('Warning: There may be something wrong with your electrode coordinate information.');
            disp('All your theta information is missing.');
            disp('Check your ced file, if that is what you used, to verify that it is properly formed.');
        end;
        eloc(cellfun(@isempty,{eloc.labels}))=[]; %if there are blank extra elocs added by something like SMI remove them
    end;
    [eloc,chanNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM,badChans,reference,facVecS,covMatrix,impedances]=ep_addEloc(whichCED,eloc,fileFormat,dataType,chanNames,timeNames,cellNames,subNames,freqNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM,covMatrix,badChans,reference,facVecS,hdr,impedances);
end;

if strcmp(ced,'internal') && any(strcmp(fileFormat,{'egi_mff_v1';'egi_mff_v2'}))
    %since mff files have internal electrode coordinates, will not normally be using ced files to designate the channels types.
    %mff file format is not sufficiently well documented yet for me to implement a more generalized fix.
    if any(strcmp('ECG',chanNames))
        chanTypes{find(strcmp('ECG',chanNames))}='ECG';
    end;
    if any(strcmp('s2_unknown259',chanNames))
        chanTypes{find(strcmp('s2_unknown259',chanNames))}='BAD';
    end;
end;

%Ensure that one-dimensional string arrays are column vectors
chanNames=chanNames(:);
timeNames=timeNames(:);
subNames=subNames(:);
cellNames=cellNames(:);
trialNames=trialNames(:);
facNames=facNames(:);
freqNames=freqNames(:);
chanTypes=chanTypes(:);
subTypes=subTypes(:);
cellTypes=cellTypes(:);
facTypes=facTypes(:);
recTime=recTime(:);
trialSpecNames=trialSpecNames(:);
subjectSpecNames=subjectSpecNames(:);

%input settings for function override the data file
if ~isempty(origRefChan)
    reference.original=origRefChan;
end;
if ~isempty(currRefChan)
    if ischar(currRefChan)
        reference.type=currRefChan;
    else
        reference.current=currRefChan;
        reference.type='REG';
    end;
end;

%Detect bad epochs and mark them
for theSub=1:numSubs
    for theWave=1:numWaves
        if ~any(any(any(sevenDdata(:,:,theWave,theSub,:,:,:),1),2),5) %do not include epochs that are zeroed out
            if any(strcmp(dataType,{'single_trial','continuous'}))
                badTrials(theSub,theWave)=1;
            else
                avgNum(theSub,theWave)=-1;
                covNum(theSub,theWave)=-1;
                subNum(theSub,theWave)=-1;
            end;
        end;
    end;
end;

%Display summary of results

if strcmp(silent,'off')
    disp(['The name of the experiment is: ' ename]);
    if ~isempty(baseline) && ~isempty(sampleRate)
        disp(['The pre-stimulus period of the data is: ' num2str(baseline * (1000/sampleRate)) ' msec.']);
        if min(samples)>1
            disp('after dropping samples as instructed.');
        end;
        if any(strcmp(fileFormat,{'egi_egis','text'}))
            disp(['Note that the data file format (' fileFormat ') is unable to specify the baseline period so it may be in error, in which case you will need to manually fix it using the Edit function.']);
        end;
        if strcmp(fileFormat,'egi_sbin')
            if ~(hdr.orig.header_array(14))==0 && (hdr.orig.header_array(15) > 1)
                disp(['Note that the data file format (' fileFormat ') is unable to specify the baseline period so it may be in error, in which case you will need to manually fix it using the Edit function.']);
            end;
        end;
    end;
end;

%Final construction of output file
if any(strcmp(dataType,{'grand_average','factors'}))
    dataType = 'average';
end;

EPdata.data=sevenDdata;
EPdata.noise=noise;
EPdata.std=stanDev;
EPdata.stdCM=stanDevCM;
if ~isempty(covMatrix)
    EPdata.cov.covMatrix=covMatrix;
    EPdata.cov.Nq=covNq;
else
    EPdata.cov=[];
end;

if isempty(dataName)
    [pathstr, name, ext] = fileparts(fileName);
    dataName=name;
end;

EPdata.montage=montage;
EPdata.fileFormat=fileFormat;
EPdata.chanNames=chanNames;
EPdata.timeNames=timeNames;
EPdata.subNames=subNames;
EPdata.cellNames=cellNames;
EPdata.trialNames=trialNames;
EPdata.facNames=facNames;
EPdata.freqNames=freqNames;
EPdata.relNames=relNames;
EPdata.chanTypes=chanTypes;
EPdata.timeUnits=timeUnits;
EPdata.subTypes=subTypes;
EPdata.cellTypes=cellTypes;
EPdata.facTypes=facTypes;
try
    EPver=ver('EP_Toolkit');
catch
    EPver='unavailable'; %workaround for bug in earlier version of Matlab
end;
EPdata.EPver=EPver;
EPdata.ver=ver;
EPdata.date=date;
EPdata.Fs=sampleRate;
EPdata.baseline=baseline;
EPdata.dataName=dataName;
EPdata.ename=ename;
EPdata.dataType=dataType;
EPdata.trialSpecs=trialSpecs;
EPdata.trialSpecNames=trialSpecNames;
EPdata.subjectSpecs=subjectSpecs;
EPdata.subjectSpecNames=subjectSpecNames;
EPdata.events=events;
EPdata.avgNum=avgNum;
EPdata.subNum=subNum;
EPdata.covNum=covNum;
EPdata.fileName=fileName;
EPdata.history=history;
% EPdata.history{end+1}={'readData',varargin};
EPdata.ced=ced;
EPdata.eloc=eloc;
EPdata.implicit=implicit;
EPdata.facVecT=facVecT;
EPdata.facVecS=facVecS;
EPdata.facVecF=facVecF;
EPdata.facData=facData;
EPdata.facVar=facVar;
EPdata.facVarQ=facVarQ;
EPdata.reference=reference;
EPdata.analysis.blinkTrial=blinkTrial;
EPdata.analysis.saccadeTrial=saccadeTrial;
EPdata.analysis.saccadeOnset=saccadeOnset;
EPdata.analysis.moveTrial=moveTrial;
EPdata.analysis.badTrials=badTrials;
EPdata.analysis.badChans=badChans;
EPdata.pca=pca;
EPdata.recTime=recTime;
EPdata.stims=stims;
EPdata.calibration=calib;
EPdata.impedances=impedances;

%Add secondary files to the data file if present

if ~isempty(specSuffix) %add in trial events information if desired
    [pathstr, name, fileSuffix] = fileparts(fileName);
    specFileName=[pathstr filesep name specSuffix];
    EPdata = ep_readEventText(EPdata, specFileName);
end;

if ~isempty(subjectSpecSuffix) %add in subject specs information if desired
    [pathstr, name, fileSuffix] = fileparts(fileName);
    specFileName=[pathstr filesep name subjectSpecSuffix];
    EPdata = ep_readEventText(EPdata, specFileName);
end;

%drop events of "trial" type as they are redundant with cell name information
for subject=1:numSubs
    for wave=1:length(cellNames)
        eventList=[];
        for event=1:length(events{subject,wave})
            if ~strcmp(events{subject,wave}(event).type,'trial')
                eventList=[eventList event];
            end;
        end;
        events{subject,wave}=events{subject,wave}(eventList);
    end;
end;

%mark epochs with a boundary event as being bad
if ~isempty(events)
    if strcmp(dataType,'continuous')
        if ~isempty(events{1})
            epochSamples=[events{1}(:).sample];
            badTrials(1,ceil(epochSamples(find(strcmp('boundary',{events{1}(:).value})))/sampleRate))=1;
        end;
    else
        for subject=1:numSubs
            for wave=1:length(cellNames)
                if ~isempty(events{subject,wave})
                    if any(strcmp('boundary',{events{subject,wave}(:).value}))
                        badTrials(subject,wave)=1;
                    end;
                end;
            end;
        end;
    end;
end;

if ~isempty(SMIsuffix) && ~all(ismember({'pupil','x-eye','y-eye'},chanNames)) %add in SMI eye-tracker information if desired and not already added
    if strcmp(dataType,'continuous')
        SMItableData=[];
        [pathstr, name, fileSuffix] = fileparts(fileName);
        SMIfileName=[pathstr filesep name SMIsuffix];
        SMImatchFileName=[pathstr filesep name '_SMImatch.txt'];
        if exist(SMIfileName,'file')
            %if already previously compiled word files then don't need to do so again.
        elseif exist([pathstr filesep name '_word01.txt'],'file') || exist([pathstr filesep name '_word01' SMIsuffix],'file')
            fileCounter=1;
            wordFid=0;
            theWordData=[];
            while wordFid ~= -1
                wordFile=[pathstr filesep name '_word' sprintf('%02d',fileCounter) '.txt'];
                wordFile2=[pathstr filesep name '_word' sprintf('%02d',fileCounter) SMIsuffix];
                if exist(wordFile,'file')
                    wordFid=fopen(wordFile);
                elseif exist(wordFile2,'file')
                    wordFid=fopen(wordFile2);
                else
                    wordFid=-1;
                end;
                if wordFid ~= -1
                    if fileCounter == 1                        
                        theData=textscan(wordFid,'%s','Delimiter','\b'); %backspace as delimiter as don't want real delimiter
                        fclose(wordFid);
                    else
                        commentLine=1;
                        numComments=0;
                        while commentLine
                            tempLine=fgetl(wordFid);
                            if strcmp(tempLine(1:2),'##')
                                numComments=numComments+1;
                            else
                                commentLine=0;
                            end;
                        end;

                        frewind(wordFid);
                        for i=1:numComments+1 %comments plus the header line
                            tempLine=fgetl(wordFid);
                        end;
                        
                        theData=textscan(wordFid, '%s','Delimiter','\b');
                        fclose(wordFid);
                    end;
                    theWordData=[theWordData;theData{1}];
                end;
                fileCounter=fileCounter+1;
            end;
            writetable(cell2table(theWordData),SMIfileName,'Delimiter','\t','WriteVariableNames',false,'QuoteStrings',false);
        else
            if exist(SMImatchFileName,'file')
                disp(['Error: No SMI files available.']);
            end;
            SMIfileName=[];
        end
        if ~isempty(SMIfileName)
            disp(['Adding SMI eye-tracking data to: ' fileName]);
            SMIfid=fopen(SMIfileName);
            if SMIfid==-1
                disp(['Error: Unable to open the SMI file: ' SMIfileName]);
            else
                if exist(SMImatchFileName,'file')
                    fid=fopen(SMImatchFileName);
                    if fid ~=-1
                        theData=textscan(fid, '%s%s','Delimiter','\t');
                        SMItableData=cell(length(theData{1}),length(theData));
                        for iCol=1:length(theData)
                            for iRow=1:length(theData{iCol})
                                SMItableData{iRow,iCol}=theData{iCol}{iRow};
                            end
                        end
                        fclose(fid);
                    end
                end
                if isempty(SMItableData) %if there is no match file, then present GUI
                    commentLine=1;
                    numComments=0;
                    while commentLine
                        tempLine=fgetl(SMIfid);
                        if strcmp(tempLine(1:2),'##')
                            numComments=numComments+1;
                        else
                            commentLine=0;
                        end;
                    end;
                    
                    delim='\t';
                    numcols=length(regexp(tempLine,delim))+1; %determine number of columns based on number of delimiters
                    
                    frewind(SMIfid);
                    for i=1:numComments
                        tempLine=fgetl(SMIfid);
                    end;
                    
                    theData=textscan(SMIfid, repmat('%s',1,numcols),'Delimiter',delim, 'MultipleDelimsAsOne', 1);
                    fclose(SMIfid);
                    MSGlines=find(strcmp('MSG',theData{2}));
                    MSGtypes=unique(theData{4}(MSGlines));
                    EEGevents={events{1}.value};
                    EEGtypes=unique(cellfun(@num2str,EEGevents(find(~cellfun(@isempty,EEGevents))),'UniformOutput',false'));
                    
                    for i=1:length(EEGtypes)
                        SMItableData{i,1}=EEGtypes{i};
                        SMItableData{i,2}='none';
                    end;
                    
                    tableNames{1}='EEG';
                    tableNames{2}='SMI';
                    columnEditable =  [false true];
                    ColumnFormat{1}=[];
                    ColumnFormat{2}=[MSGtypes', 'none'];
                    
                    handles.read.SMI = figure('Name', 'Match SMI events', 'NumberTitle', 'off', 'Position',[1 scrsz(4)-1100 200 500], 'MenuBar', 'none');
                    colormap jet;
                    
                    handles.read.SMItable = uitable('Data',SMItableData,'ColumnName',tableNames,'FontSize',FontSize,...
                        'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,'Position',[0 50 200 450]);
                    
                    handles.read.SMItableDone = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',FontSize,...
                        'Position', [ 0 0 60 20], 'Callback', 'uiresume(gcbf)');
                    
                    uiwait(handles.read.SMI);
                    SMItableData=get(handles.read.SMItable,'Data');
                    close(handles.read.SMI);
                    drawnow;
                end
                
                if isempty(SMItableData)
                    disp(['Error: Addition of SMI eye-tracker information failed.']);
                else
                    [EPdataSMI logMsg]=ep_readSMI(EPdata,SMIfileName,SMItableData);
                    if isempty(EPdataSMI)
                        disp(['Error: Addition of SMI eye-tracker information failed.']);
                    else
                        EPdata=EPdataSMI;
                    end;
                end;
            end
        end;
    end;
end;

%Select subsets of data

badCEDchans=find(strcmp('BAD',chanTypes));
if ~isempty(badCEDchans)
    if isempty(channels)
        channels=[1:numChan];
    end;
    channels=setdiff(channels,badCEDchans); %drop "channels" that the CED file indicates should be deleted entirely.
    for i=1:length(badCEDchans)
        msgString=['Dropping channels: '];
        for i=1:length(badCEDchans)
            msgString=[msgString chanNames{badCEDchans(i)} ';'];
        end;
    end;
    disp(msgString);
end;

EPdata=ep_selectData(EPdata,{channels,samples,cells,subjects,[],[]});
if isempty(EPdata)
    return
end;

if ~isempty(cells) || ~isempty(subjects)
    if length(trialNames)==1
        trialStr='trial';
    else
        trialStr='trials';
    end
    if length(cellNames)==1
        cellStr='cell';
    else
        cellStr='cells';
    end
    if strcmp(dataType,'factors')
        if numFacs==1
            subStr='summed factor';
        else
            subStr='factors (including the summed factor)';
        end
    else
        if numSubs==1
            subStr='subject';
        else
            subStr='subjects';
        end
    end;
    if strcmp(silent,'off')
        if strcmp(dataType,'single_trial')
            disp(['After selecting, there are ' num2str(numCells) ' ' cellStr ' with a total of ' num2str(length(trialNames)) ' ' trialStr '.']);
        elseif strcmp(dataType,'factors')
            disp(['After selecting, there are ' num2str(numCells) ' ' cellStr ' and ' num2str(numFacs) ' ' subStr '.']);
        else
            disp(['After selecting, there are ' num2str(numCells) ' ' cellStr ' and ' num2str(numSubs) ' ' subStr '.']);
        end;
    end;
end;

[err]=ep_checkEPfile(EPdata);

if err
    EPdata=[];
    msg{1}='Defective file will not be loaded.';
    if strcmp(fileFormat,'ep_mat')
        msg{2}='File can be fixed by loading manually into Matlab and editing, as in load(''filename.ept'',''-mat''); and then save(''filename.ept'',''-mat'');';
    end;
    [msg]=ep_errorMsg(msg);
    return
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indexOut fieldOut]=readBVheaderField(eventHdr,indexIn) 
%translate a field of BV header info
fieldType=str2double(eventHdr(indexIn).value(2:end));
indexOut=indexIn;
if isempty(fieldType) || isnan(fieldType) || (fieldType < 1) || (round(fieldType) ~= fieldType)
    fieldOut=NaN;
    return
end;
switch fieldType
    case 1 %uint
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end;
        fieldOut=num2str(str2double(eventHdr(indexIn+1).value(2:end)));
        indexOut=indexOut+2;
    case 2 %slong
        if (indexIn+2) > length(eventHdr)
            fieldOut=NaN;
            return
        end;
        fieldOut=num2str(str2double(eventHdr(indexIn+1).value(2:end))+bitshift(str2double(eventHdr(indexIn+2).value(2:end)),8));
        indexOut=indexOut+3;
    case 3 %text
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end;
        textLength=str2double(eventHdr(indexIn+1).value(2:end));
        if isnan(textLength) || isempty(textLength) || (textLength < 1) || (round(textLength) ~= textLength)
            fieldOut=NaN;
            return
        end;
        theChar='';
        indexIn=indexIn+2;
        for iChar=1:textLength
            [indexIn,charOut]=readBVheaderField(eventHdr,indexIn);
            if ~isnan(charOut)
                theChar(iChar)=char(str2double(charOut));
            else
                theChar(iChar)=' ';
            end;
        end;
        indexOut=indexIn;
        fieldOut=theChar;
    case 4 %zero unint
        fieldOut='0';
        indexOut=indexOut+1;
    case 5 %low byte slong zero
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end;
        fieldOut=num2str(bitshift(str2double(eventHdr(indexIn+1).value(2:end)),8));
        indexOut=indexOut+2;
    case 6 %hi byte slong zero
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end;
        fieldOut=num2str(str2double(eventHdr(indexIn+1).value(2:end)));
        indexOut=indexOut+2;
    case 7 %slong zero
        fieldOut='0';
        indexOut=indexOut+1;
    case 8 %empty text
        fieldOut='';
        indexOut=indexOut+1;
    otherwise
        fieldOut=NaN;
        return
end;



