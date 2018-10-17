function ep(varargin);
% ep - function ep(varargin);
% GUI interface for the EP Toolkit.
%

%History
%  by Joseph Dien (7/30/09)
%  jdien07@mac.com
%
%  bugfix 8/2/09 JD
%  EPdataset.dataset initialized as an empty string rather than as a cell, which was causing a crash when adding the
%  first dataset.
%
%  modified 8/3/09 JD
%  Consolidated EGI EGIS average and session entries in the file format list.
%
%  bugfix 8/11/09 JD
%  For artifact correction, bug resulted in EGI EGIS session format always being used as input format.  Thanks
%  to Grega Repov.
%
%  bugfix 8/27/09 JD
%  Location of windows appearing partly off screen on some systems.
%  ANOVA preferences pane crash fixed and crash when doing ANOVA with only one between group and add option on.  Thanks to Tim Curran.
%  File title missing from first ANOVA in output files.  Thanks to Tim Curran.
%
%  modified 8/30/09 JD
%  Added preferences for importing text files.
%
%  bugfix and modified 9/5/09 JD
%  Crash when generating spatial PCA.  Thanks to Tim Curran.
%  Move factors preference were not used for preprocessing. 
%  Preprocess controls now remember settings.
%  Single cell files option for preprocessing were not functional.
%  Added reference channel field to preprocessing pane to deal with ICA problems when processing mean mastoid data.
%  Blink file was not passed on to blink correction function for fileTemplate option.
%  Crash when saving preferences file on a non-Mac.
%  Dataset information like the cells were not updating when changing datasets in the Window Data Pane.
%
%  bugfix 9/10/09 JD
%  Fixed not saving file when averaging session files and outputting as EGIS format.
%  Fixed crash when preprocessing data and using autoTemplate for blinks.
%
%  bugfix 9/23/09 JD
%  Not finding EPwork directory when in the search path but not in the active directory.
%
%  bugfix 10/12/09 JD
%  Fixed crash when using single file mode with Read function.
%
%  modified & bugfix 10/28/09 JD
%  Added option to disable preprocessing figure for low memory situations.
%  Movement correction preferences not registering.
%
%  modified 10/31/09 JD
%  When reading multiple files, assume all have the same electrode location and montage information rather than
%  repeatedly asking for them.
%  Cached more information in EPdataset to avoid slowdowns in main pane and preprocess pane.
%
%  modified 10/31/09 JD
%  checks to see if EPwork cache is current and regenerates it if not.
%
%  bugfix 11/2/09 JD
%  Fixed crash when using single file mode with Read function.
%
%  bugfix 11/4/09 JD
%  Fixed crash when running PCA scree and no noise PCA was conducted.  Thanks to Jia Wu.
%
%  bugfix 11/12/09 JD
%  Fixed bug where Toolkit was saving all variables into EPwork cache, not just EPdataset, resulting in erratic problems
%  when old values for variables were loaded back into the workspace.
%
%  bugfix 12/3/09 JD
%  Added workaround for sporadic Matlab menu bug that was causing menu to disappear or crash when menu item that is selected is
%  same as before.
%
%  bugfix 12/10/09 JD
%  Fixed crash when viewing Window pane for a dataset with less than four cells under Matlab versions prior to 2008.
%
%  bugfix 1/8/10 JD
%  Fixed crash when entering between group name that is not three letters long.
%  Fixed failure due to Matlab bug to redraw View pane when changing the dataset for one of the colors.
%
%  bugfix & modified 1/17/10 JD
%  Fixed crash when examining a file in the View pane with less than four cells.
%  Modified ANOVA module to accommodate non-ERP data by adding support for "behavioral" keyword.
%
%  modified 1/28/10 JD
%  Can leave out subject and cell fields when reading files in single file mode
%  (will use default values instead and will assume all the subjects/cells are the same).
%  Can accomodate file formats such as Neuroscan .avg where channels can be deleted, not just zeroed out, if bad and
%  channels can appear in an order different than that of the .ced file.  When importing a file with missing channels,
%  will now add them in as zeroed out channels and will mark them as bad.  Will also reorder channels to match that of
%  the ced file.
%  The type field in the .ced files is now used since the EEGlab bug was apparently fixed and it is now functional.  The
%  type field must now be present in the .ced file and assumptions will no longer be made about which ones are REF or FID types.
%  Added ability to plot data from continuous data files by showing only one second at a time.
%  Epoch ms now displayed as from beginning of first sample to end of last sample.
%  Fixed bug where checkboxes and subject and cell fields on Read and Preprocess panes impossible to deselect.
%  Refreshes Edit pane when Done button is pressed so when Data Name is changed, the name in the list of datasets is updated.
%
%  bugfix and modified 2/28/10 JD
%  Average function now correctly checks to see if output file name already exists for EP format files.
% analysis fields no longer optional.
% Now has option to use file's default reference settings in the preprocess pane.
% When running ANOVAs with between subject factors, fixed the addition of grand averages corresponding to the levels of the factors to the dataset.
% Made autoPCA more memory efficient so it wouldn't run out of memory.
% Fixed crash in ANOVA function when level names of between group factors were of different lengths.
% Close files after running batch of ANOVAs to avoid "too many files open" error.
% When running ANOVAs with between subject factors using the output of the autoPCA option, fixed the addition of grand averages corresponding to the levels of the factors to the dataset.
% Corrected the peak channel and peak time point identification of factors made by the autoPCA option of the window function.
% Fixed bug that prevented one from fully deleting the contents of edit fields, such as "baseline" on the preprocessing pane.
% Fixed bug that prevented one from changing the contents of the edit fields on the postprocessing pane.
% Fixed inability to set minimum variance setting for autoPCA preference to anything other than zero.
% Added option to turn off adding montage information to EGIS file format files due to incompatibility issues with some
% versions of NetStation.
%
% bugfix 3/6/10 JD
% When running ANOVAs with between subject factors, fixed having only the first between group waveform being added to
% the datafile.
%
% modified 3/27/10 JD
% Added ability to save data using .set file format.
%
% bugfix 4/26/10 JD
% Fixed wrong peak chans and time points for AutoPCA when not two-step PCA.
% Added topos function to the View Pane.
%
% bugfix 5/5/10 JD
% Fixed crash when in Window Data pane using noTable mode (older versions of Matlab) and there are less than five
% subject specs.
% Fixed not changing table of cells and table of specs when in Window Data pane using noTable mode (older versions of
% Matlab).
% Fixed crash when using Topo button of View EEG pane and not all four colors are being used.
% Added option to use unrotated solution for PCAs.
%
% modified 5/22/09 JD
% Number of factors and title of PCA no longer being changed to blank when changing other PCA settings.
% Added support for EEGlab .study files.
%
% bugfix 6/5/10 JD
% Fixed crash when postprocessing.
% View function no longer crashes when average file erroneously has non-empty trial names field.
% Fixed files sometimes not being recognized as being selected when names are in uppercase.
%
% modified & bugfix 6/17/10 JD
% Marks when a file isn't saved yet.
% Changed the suffix of EP files to ".ept".
% Fixed crash when windowing adds regional channel and there is no electrode coordinate information (eloc).
% Fixed crash when conducting robust ANOVA and there are no spec columns in the ANOVA data file.
%
% bugfix 6/28/10 JD
% Fixed when no within factors for ANOVAs, pane indicated needed one within cell rather than zero wthin cells.
%
% bugfix 7/2/10 JD
% In order fix crash due to Matlab bug on Windows, put in brief pause after spec name requestor diaglog box.
%
% bugfix 7/23/10 JD
% Fixed view pane listing all the trial names of the dataset rather than just those for a specific cell, resulting in
% crashes when they were selected.
% Fixed crash when viewing Waves and data are all zero.
% Fixed grand averages, which are added when computing ANOVAs with between factors, being computed incorrectly (only
% affected waveforms, not the ANOVA results).
%
% bugfix 8/27/10 JD
% Fixed unable to output averages in EGIS format.
% 
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
% 
%  modified 10/28/10 JD
%  Added support for saccade preprocessing preference settings.  Removed the timepoints option from the preprocessing
%  pane.
%
% bugfix 12/7/10 JD
% Fixed crash in Window function after switching to a dataset with fewer cells.
% Fixed crash in ANOVA function when performing ANOVA with no between group factors under Matlab 2008b.
%
% bugfix 1/5/11 JD
% Fixed crash in ANOVA function when used prior to data being loaded into the work set yet (due to badly initialized variable).
% Now putting preferences file and work directory at default user directory if old one cannot be found.
% 
% modified 1/18/11 JD
% Added support for selecting timepoints in preprocessing (previous implementation was non-functional).
% 
% modified 1/20/11 JD
% Added support for manually specifying EOG channels in the preferences.
%
% bugfix 2/3/11 JD
% Fixed crash when starting program on a computer with spaces in the default user path, as in "Documents and Settings"
%
% bugfix 2/4/11 JD
% Hopefully fixed crash when starting program with certain configurations of EPwork and EPprefs.
%
% bugfix 2/8/11 JD
% Hopefully fixed crashes from buggy reset of preferences.
%
% bugfix 2/10/11 JD
% Fixed crash when performing eyeblink correction due to change in preferences file.
% Fixed not using saccade preference settings.
% Fixed crashing when loading EPdataset files.
% 
% modified 1/20/11 JD
% Added commands to change or create EPwork directory from the EP menu.
% Only stores EPpref files in EPwork directories.
%
% modified 2/26/11 JD
% Added support for EGI Matlab files
%
% bugfix 4/18/11 JD
% Fixed output file from the Window function having the suffix ".txt.txt"
%
% bugfix 4/22/11 JD
% Fixed crash when seeking to form manual blink templates and the only available files have missing channel coordinates.
%
% bugfix 5/5/11 JD
%Fixed not saving reset preference values when preference file found to be out-of-date.
%
% modified 5/5/11 JD
% Added support for FFT analyses.
%
% bugfix 6/2/11 JD
% Fixed crash when userpath is blank and Matlab is version 2007 or presumably earlier.
%
% modified 6/7/11 JD
% Added support for specification of current reference in preprocessing pane.
%
% bugfix 6/20/11 JD
% Fixed not saving preference changes and now looks for the preferences file in the current working directory
%
% bugfix 7/8/11 JD
% Changed separation character for added grand average when computing robust ANOVAs from "|" to "_" as the former was
% causing crash on PCs when saving files to text format.
% 
% modified 1/24/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
%
% bugfix 1/27/12 JD
% Fixed Window pane not allowing channel groups to be chosen using PCA results when the data file has a regional
% channel.
% 
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% modified 4/11/12 JD
% Added support for 6th dimension of frequency.  Unlogs frequency data for computations.
%
% modified 4/23/12 JD
% Changed default for channel group windowing to collapsing channels first
% then taking measures.  Added preference option so can use original
% approach too of measuring first.
%
% modified 5/24/12 JD
% Added support for missing data when adding channels, cells, or subjects.
%
% modified 5/24/12 JD
% Eliminated "no ref" option for preprocessing as no longer needed and confusing.
%
% bugfix 6/4/12 JD
% Fixed regional channel waveform added to dataset when windowing not computed correctly 
% when there are more than one channel group and the last one to be edited by the channel group function is not the one being used. 
%
% modified 8/6/12 JD
% Added suport for coherence and phase locking measures.
%
% bugfix 9/3/12 JD
% Fixed crash when trying to add subject ANOVA waveform after computing robust ANOVA and subjects have been trimmed.
% Fixed crash when trying to add subject ANOVA waveform after computing robust ANOVA and data is not from PCA output.
%
% modified 10/16/12 JD
% Added option to set last row to be imported in text files.
%
% modified 1/11/13 JD
% Added option to do internal calculations of frequency data in either amplitude or power form.
%
% bugfix 1/11/13 JD
% Fixed erroneous error message and crash when trying to PCA average file where bad channels were dropped rather than replaced.
%
% modified 1/17/13 JD
% Added support for ERPlab .erp file format.
%
% bugfix 1/18/13 JD
% Fixed erroneous "labels" error message when trying to load .study file.
%
% bugfix 1/23/13 JD
% Fixed problem where information for expanding channels in waveform plot is lost under some circumstances.
%
% modified 1/28/13 JD
% Added manual editing option for artifact correction.
%
% modified 1/29/13 JD
% Added frequency-domain peak measures to windowing function.
%
% modified 1/30/13 JD
% When saving a dataset, name of dataset changes if new name is chosen for saved dataset.
%
% bugfix 2/1/13 JD
% Fixed crash in View pane under when changing a dataset under certain circumstances.
%
% modified and bugfix 2/1/13 JD
% Clearing volt, hz, and sample parameters in View pane no longer crashes.  If value is manually set, will not change
% until a new dataset is chosen or until the value is cleared (in which case it will be replaced with automatic value).
% Added support for reading/writing ERPlab datasets.
%
% bugfix 2/4/13 JD
% Fixed error when setting file format preferences to ERPlab files.
% Fixed failure to save preferences when save preferences button clicked.
% Fixed ms window information on Window pane so that it ranges from onset of sample to offset of sample rather than
% offset to offset (e.g., 4-4 ms rather than 0-4 ms for first sample) after changes to the windowing settings.
%
% modified 2/6/13 JD
% Allow View Scan function to operate on average files.
%
% modified 2/21/13 JD
% Added version number to preferences to allow for automatized updating processes.
%
% modified 4/1/13 JD
% Improved the controls for the Windowing pane.
%
% bugfix 4/2/13 JD
% Markers in waveform plots can be set at zero ms.
%
% modified 4/5/13 JD
% Addressed change to userpath function in Matlab 2013a.
%
% bugfix 4/12/13 JD
% Scaling of topos now obeys the values on the View pane.
%
% bugfix 4/24/13 JD
% Fixed choosing "auto" as Edit Mode setting in Preprocessing pane resulting in error message.
%
% modified 5/6/13 JD
% Implemented adaptive mean in Windows pane by adding width setting for peak measure.
%
% bugfix 5/8/13 JD
% Fixed detrend option in Preprocess Data pane not working.
%
% bugfix 5/9/13 JD
% Fixed Scan allowed to be activated if the scaling is from zero to zero, as with a bad trial, causing a crash.
%
% bugfix 5/18/13 JD
% Fixed Single File Mode in Preprocessing pane crashing when values entered.
%
% bugfix 5/20/13 JD
% Fixed Single File Mode in Preprocessing pane crashing.
%
% modified 5/20/13 JD
% Single-Trial Files from multiple subjects can now be selected using the Read pane's Single File Mode and read in as separate files.
%
% bugfix 6/18/13 JD
% Fixed windowing pane crashing when Herz bins changed.
%
% bugfix 8/20/13 JD
% Fixed Change Work Directory menu function not working.
%
%  modified 9/14/13 JD
%  Initializes mff at the outset so that global variables bug doesn't crash the toolkit down the line.
%
% modified 9/16/13 JD
% Added support for Biosemi bdf files
%
% bugfix 9/17/13 JD
% Fixed not making factor data available for setting channel groups when first in analysis set had different number of channels.
% 
% modified 10/10/13 JD
% Supports presence of non-EEG channels.
%
%  bugfix 11/1/13 JD
%  Fixes font sizes on Windows.
%  Fixed About EP Toolkit spawned new main window.
%
%  modified 11/6/13 JD
%  Scan and Waves functions can now present event markings.
%
%  modified 11/13/13 JD
%  Added Segment Data function.
%
%  modified 11/27/13 JD
%  Added fMRI artifact correction option to Preprocess data function using both EEGlab fMRIb plugin and AMRI EEG fMRI Toolbox.
%
%  modified 12/23/13 JD
%  Added windowing function option to window 'all' the channels.
%
%  bugfix 12/24/13 JD
%  Fixed between group ANOVAs not being calculated correctly when the
%  between group variable is not sorted alphabetically.
%
%  bugfix 1/12/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
%  bugfix 1/22/14 JD
%  Added Trim Data option to the Segment Data function.
%
% modified 2/26/14 JD
% Added View function option to plot or erpimage all trials and all subjects.
%
%  bugfix 2/27/14 JD
%  Fixed PCA not correctly recognizing that a file has already undergone two-step PCA.
%
% bugfix 3/5/14 JD
% Fixed when using single file mode to read in single-trial files, all the resulting files are identical to the very first subject.
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% modified 3/16/14 JD
% Uses internal electrode coordinates provided by MFF and FIFF files.  Added elecPrefs.
%
% modified 3/19/14 JD
% Added support for saving data in Neuromag FIFF file format.
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
% Added option to window single-trial data.
% Windowing outputs cells in actual order rather than alphabetical order.
% Eliminated noTable option for old versions of Matlab.
%
% bugfix 3/19/14 JD
% Fixed segmenting table + button mirroring the first line of the table rather than the current settings above the table.
% Fixed segmenting table - button deleting all but second to last line rather than just the last line.
%
% bugfix 3/23/14 JD
% Fixed waveforms added during ANOVAs to correspond to trimmed cell means being all flat.
%
% bugfix 3/26/14 JD
% Fixed crash when running an ANOVA on a PCA dataset results in trimmed cell means being added to it.
%
% modified 4/7/14 JD
% Added "all" and "erpimage" options to the Factors list in View.
%
% bugfix 4/8/14 JD
% Fixed AutoPCA generating nothing but missing data values when maxcentroid and mincentroid measures chosen.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% modified 4/24/14 JD
% Added coherence and phase-locking options and allows choosing either hanning or multi-taper methods for spectral measures.
%
% bugfix 4/24/14 JD
% Fixed 'all' option for Window Data leaving out last channel.
%
% bugfix 4/29/14 JD
% Fixed when spectral range changed in View function and dB or psd options are on, the  values are immediately further transformed.
%
% bugfix 5/13/14 JD
% Fixed crash when running an ANOVA on a windowed text file generated by autoPCA and the adds option is on. 
%
% bugfix 5/20/14 JD
% Fixed minimum and maximum voltages in View pane not reflecting correct values for single trial data.
%
% modified 5/20/14 JD
% Added support for BrainVision EEG files
%
% bugfix 5/29/14 JD
% Allow for manual windowing of PCA files rather than just autoPCA.
%
% bugfix 6/1/14 JD
% Fixed list of trials for single-trial data in View pane not correct.
%
% modified 6/3/14 JD
% Added support for SMI eye tracking files.
%
% modified 6/18/14 JD
% Added starts, ends, contains, and -follows- keywords to the Segment function.
%
% modified 6/19/14 JD
% Added support for sample-by-sample t-tests, including STS chanType.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% modified 8/12/14 JD
% Added support for PARE-corrected average reference.
%
% bugfix 8/27/14 JD
% Fixed crash for files where event values are all numbers rather than
% strings.
%
% modified 8/31/14 JD
% Added field to Transform panel indicating how much of a delay has been added to ERP latencies by one-pass filters.
% Added support for adding additional keys information to events in continuous data and trial specs in single-trial
% data.
%
% modified 9/4/14 JD
% Added delay field to segment function.
%
% modified 9/15/14 JD
% Added contents of type field to time-lock events list in segment function.
%
% modified 10/16/14 JD
% Passes screen size to ep_readData call.
%
% bugfix 10/21/14 JD
% Fixed cannot update the preferences for text event and SMI file suffixes.
%
% bugfix 10/24/14 JD
% Fixed crash when using Transform pane and rereference set to "none."
%
% bugfix 12/8/14 JD
% Fixed crash when selecting mff files on a PC.
%
% bugfix 12/14/14 JD
% Fixed crash in sampTest function when there are bad trials present in single-trial data.
%
% modified 12/23/14 JD
% Added convert option to Save function.
%
% bugfix 1/24/15 JD
% Fixed crash in Save function when name of dataset to be saved is changed.
%
% bugfix 3/25/15 JD
% Fixed crash when reading in multiple files.
%
% bugfix 5/29/15 JD
% Fixed crash in Segment pane when dataset has no events.
%
% modified 5/29/15 JD
% Changed min and max TopoPlot scale to be set by plotted data unless overriden.
%
% bugfix & modified 9/25/15 JD
% Fixed Window button on main pane not becoming active for single-trial
% data being present in the active set.
% Fixed crash when working set newly initialized and running an ANOVA on
% behavioral data with adds option activated.
% Fixed crash when ANOVA conducted on data with characters instead of
% numbers.  Now treated as missing data.
% Added capability to resegment single-trial data to reassign cell
% conditions.
% Added capability to specify OR criteria for a condition by just giving
% them the same name.
% Added trial specs for average files.
% Changed default head rotation for electrode coordinates of mff and fiff files to 90 degrees.
%
% modified 10/13/15 JD
% On the assumption that EGI users have largely migrated from GSN200 to
% Hydrocel nets, the default right mastoid channel has been changed from
% 101 to 100.
%
% bugfix 11/19/15 JD
% Corrected erroneous error message that degrees of freedom of robust ANOVA
% is calculated per multivariate rather than univariate approach.
%
% modified 11/25/14 JD
% Added cross-validation option to PCA.
%
% modified 12/10/15 JD
% Checks if p-value variability exceeds two standard deviations over the threshold.
%
% modified 12/18/15 JD
% On further thought, changed cross-validation option to apply FacPat rather than FacCof.
%
% bugfix 1/15/16 JD
% Fixed when clicking on Views channels, can only get expanded channel window if one clicks along the very top edge for FFT data.
% Fixed upper and lower amplitude/power being changed immediately after new values manually entered into the View panel for FFT data.
% Fixed min and max values displayed in Views pane divided by Hz bin rather than sqrt(Hz bin) when set to amplitudes.
% Now allows power scaled data to be displayed as amplitudes.
% Now handles complex FFT numbers.
%
% modified 1/22/16 JD
% Consolidated spectral unit controls so just four options (cmp, asd, psd, dB).
% If minimum value in View pane is zero, in dB scaling will set to -2 instead of ?inf to maintain useful range.
% Eliminated option to transform data in power rather than amplitude form to avoid potential confusion.
%
% modified 2/19/16 JD
% Fixed jack-knife not calculated for first step of PCAs.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% bugfix 3/21/16 JD
% Fixed PCA Woody option of SampleTest function not saving latencies and amplitudes to the correct trials.
%
% bugfix 5/22/16 JD
% Fixed crash when using average function and the active set is empty.
%
% bugfix 8/23/16 JD
% Fixed crash when changing the Trial Spec Model ofthe Latency-Lock procedure of the Average function.
%
% modified 9/16/16 JD
% Removed saccade and saccademin preferences.
%
% bugfix & modified 11/8/16 JD
% Added support for template Woody and for continuous data to Sample Test function.
% Added support for reading and writing subject spec text files.
% Fixed crash in segment function after changing line after -follows- to 'none'.
%
% modified 1/3/17 JD
% Added eyeTracker option to blink and saccade correction routines.
%
% bugfix 3/2/17 JD
% Fixed crash when switching to jitter-correct in Average Pane and there are no single-trial datasets in the working set.
%
% modified 4/19/17 JD
% When quiting, clear EPwork option leaves both the directory and the preferences file intact so there is no need to redo changes to the preference settings.
%
% bugfix 5/4/17 JD
% Fixed crash when setting the Mark fields in the View function.
%
% bugfix 5/19/17 JD
% Fixed crash when merging multiple subjects with single cell mode in Read function.
% Added option to flip the electrode locations for .mff and .fiff.
%
% bugfix & modified 6/19/17 JD
% Flexible segments implemented in Segmentation function.
% Added option to clear working set to File menu.
% Added support for averaging multi-file subjects.
% Now allows just one Hz bin of TFT data to be chosen for display using the View function.
% PCA now drops non-EEG channels and regional EEG channels.
% Added support for using sampleTest function with TFT data.
% Fixed error when applying sample method to average data in sampleTest function.
% Fixed Window pane crashing when first dataset in working set is not an average file.
% Fixed Window function saves both imaginary and real files for spectral data even when units are not set as complex.
% Scaling input boxes in View pane now allow up to four decimals to better accommodate power units for spectral data.
% Fixed spectral density conversion incorrectly applied to View pane controls for complex and amplitude scaling, resulting in the numbers being changed from what was input.
% Now allows just one sample of TFT data to be chosen for display using the View function.
%
% bugfix & modified 10/4/17 JD
% Added EMG correction option to Preprocessing function.
% Segment function user interface now providing correct default sample range for sampling rates other than 250 Hz.
% Segment function user interface now provides only unique set of values for integer trials specs.
%
% bugfix & modified 11/26/17 JD
% Added -precedes- crit option to the Segment function.
% Fixed crash when changing value of a -precedes- or -follows- crit to an event name that has no keys.
% Segment function user interface only includes events with identical Type and Value once in Event Lock list.
% Segment function user interface does not include TRSP in Event Lock list.
% Segment function user interface always provides "100" as an option in a list of TS-RT values.
% Fixed cross-validation option giving wrong results!  Also changed name to "cross-verification."
%
% modified 2/9/18 JD
% Added support for GFP plots and error bands.
% Added support for stdCM field.
%
% modified 2/23/18 JD
% When computing an add's std, no longer tries to combine std values.  Instead sets to zero.
% Added -all- and -erpimage- options for cells in View Waves function.
%
% bugfix 3/9/18 JD
% Fixed cells popupmenu listing all trials rather than just one of each name.
%
% bugfix 3/27/18 JD
% Fixed crash when performing scree on second step of a two-step PCA.
%
% modified 4/29/18 JD
% Added support for calculation of RT and accuracy summary scores during averaging.
% Added support for outputting behavioral data in the Window function.
% Added support for View Topos listing all trials or cells or subjects.
% Cell box for Window Data function now allows blanks so may edit freely to change order or to drop cells.  'none' no longer recognized as a keyword.
% Cell box for Window Data function no longer sorting cell names alphabetically for average files.
% Window Data pane settings no longer resetting every time one leaves the pane.
%
% bugfix 5/13/18 JD
% Fixed crash sometimes when starting Window Data pane.
%
% bugfix 5/18/18 JD
% Fixed failing to perform PCA scree test when noise field consists of zeros.
%
% bugfix & modified 6/5/18 JD
% Added preference option to turn off BV header encoding.
% Sorts batch files by name prior to running ANOVAs, averages, and preprocessing.
% Fixed not adding new cells as adds (if enabled) if new cell names specified during windowing.
% Fixed crash during ANOVA if a cell name is specified that is not present in the original data if adds are enabled.
% Fixed crash in Save function when clicking on table rather than using Convert mode.
% Don't remove star for unsaved data when save is cancelled for whatever reason.
% Added option to add marks for RT and selected events to View figures.
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

global EPdataset EPmain EPoverview EPwaves EPscree EPchanGrp EPtopos EPmanualEdit

refList={'ave ref','1-2 refs','CSD','default'};

if isempty(varargin)
    
    if ~exist('ep_initEP','file')
        disp(' ');
        disp('**************************************************************');
        disp('Error: The rest of the EP Toolkit is not on the path.  Did you use');
        disp('the ''Add with Subfolders'' button rather than the ''Add Folder'' button?');
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
        return
    end;
    
    if ~isempty(EPmain)
        if isfield(EPmain,'mode')
            if ~any(strcmp(EPmain.mode,{'agree','disagree'}))
                EPmain=[];
            end;
        else
            EPmain=[];
        end;
    end;
    if isempty(EPmain)
        splash
        return
    elseif strcmp(EPmain.mode,'disagree')
        close('About EP Toolkit')
        EPmain=[];
        return
    else
        close('About EP Toolkit')
        EPmain=[];
    end
    
    abortEP=ep_initEP;
    if abortEP
        return
    end;
    eval(['global EPdataset EPmain EPoverview EPwaves EPscree EPchanGrp EPtopos EPmanualEdit']); %restore the global variable linkages after the mff bug has eliminated them.
    
    EPmain=[];
    EPmain.fileFormatReadList={'EP (.ept)','EGI EGIS (.egis)','EGI Simple Binary (.sbin)','EEGlab (.set/.study)','Biosig EDF (.edf)','Neuroscan (.cnt/.eeg/.avg)','text (.txt)','EGI Matlab (.nsf)','EGI MFF (.mff)','ERPlab (.erp)','Biosemi (.bdf)','Neuromag (.fif)','BrainVision (.eeg/.dat/.seg)'};
    EPmain.fileFormatSaveList={'EP (.ept)','EGI EGIS (.egis)','EGI Simple Binary (.sbin)','EEGlab (.set/.study)','text (.txt)','ERPlab (.erp)','Neuromag (.fif)'};
    
    scrsz = get(0,'ScreenSize');
    if max(scrsz) <2
        msg{1}='Error: Matlab is currently unable to determine the size of the monitor.  Please restart Matlab.';
        [msg]=ep_errorMsg(msg);
        return
    else
        EPmain.scrsz=scrsz;
    end;
    
    PathName=userpath;
    if ~isempty(PathName)
        if any(strcmp(PathName(end),{';',':'}))
            PathName=PathName(1:end-1);
        end;
        if ismac
            colonPos=findstr(PathName,':');
            if ~isempty(colonPos)
                PathName=PathName(1:colonPos-1);
            end;
        end;
    end;
    %set up work directory
    disp('Loading summary information about datasets in working set (The more datasets in the working set, the longer this will take).');
    if exist([pwd filesep 'EPwork' filesep 'EPdataset.mat'],'file') %if there is already a work directory in the active directory, use it.
        eval(['tempVar=load(''EPwork' filesep 'EPdataset.mat'');']);
        if isfield(tempVar,'EPdataset')
            EPdataset=tempVar.EPdataset;
        end;
        clear tempVar;
        EPdataset.EPwork=pwd;
    elseif exist(['~/Documents/EPwork' filesep 'EPdataset.mat'],'file') %on a Mac, work directories kept in this location
        eval(['tempVar=load(''~/Documents/EPwork' filesep 'EPdataset.mat'');']); %if there is already a work directory, use it.
        if isfield(tempVar,'EPdataset')
            EPdataset=tempVar.EPdataset;
        end;
        clear tempVar;
        EPdataset.EPwork='~/Documents';
    elseif exist([PathName filesep 'EPwork' filesep 'EPdataset.mat'],'file')
        eval(['tempVar=load(''' PathName filesep 'EPwork' filesep 'EPdataset.mat'');']); %if there is already a work directory, use it.
        if isfield(tempVar,'EPdataset')
            EPdataset=tempVar.EPdataset;
        end;
        clear tempVar;
        EPdataset.EPwork=PathName;
    else
        if isempty(PathName)
            try
                userpath('reset')
                disp('Resetting userpath as it was set to root');
                PathName=userpath;
                if any(strcmp(PathName(end),{';',':'}))
                    PathName=PathName(1:end-1);
                end; 
            catch
                %userpath('reset') doesn't work for Matlab 2007 and presumably earlier
                disp('Putting EP Toolkit work directory in the active directory.');
                PathName=pwd;
            end;
        end;
        
        mkdir([PathName filesep 'EPwork']);
        EPdataset=[];
        EPdataset.EPwork=PathName;
        EPdataset.dataset=cell(0);
        eval(['save ''' PathName filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
    end;
    disp('Done loading summary information.');
    
    %set up preferences
    if exist([EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat'],'file')
        checkEPdataset=EPdataset;
        tempVar=load([EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat']);
        if isfield(tempVar,'prefs')
            prefs=tempVar.prefs;
        end;
        clear tempVar;
        if length(EPdataset.dataset) ~= length(checkEPdataset.dataset)
            resetPrefs;
            prefs=EPmain.preferences;
            eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs'' prefs']);
            disp('EPprefs file is corrupted.  Resetting preferences file.');
            EPdataset=checkEPdataset;
        end;
        EPmain.preferences=prefs;
    else
        resetPrefs;
        prefs=EPmain.preferences;
        eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs'' prefs']);
    end;
    
    err=checkPrefs;
    
    if err
        disp('Preferences are corrupt or out-of-date.  Resetting to default values.');
        resetPrefs;
        prefs=EPmain.preferences;
        eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs'' prefs']);
    end;
    
    EPdataset=ep_checkEPworkCache(EPdataset);
    
    varargin{1}='start';
    EPoverview=[];
    EPmain.handles.hMainWindow =[];
    EPmain.mode='main';
    
    if ispc
        EPmain.fontsize=8;
    else
        EPmain.fontsize=10;
    end;
    
    EPmain.segment=[];
    EPmain.segment.importFormat=EPmain.preferences.general.sessionImportFormat;
    EPmain.segment.outputFormat=EPmain.preferences.general.outputFormat;
    EPmain.segment.delay=0;
    
    EPmain.average.importFormat=EPmain.preferences.general.sessionImportFormat;
    EPmain.average.type=2;
    EPmain.average.outputFormat=EPmain.preferences.general.outputFormat;
    EPmain.average.method=EPmain.preferences.average.method;
    EPmain.average.freqMethod=1;
    EPmain.average.smoothing=1;
    EPmain.average.minLatency=[];
    EPmain.average.maxLatency=[];
    EPmain.average.RTmethod=1;
    EPmain.average.minRT=100;
    EPmain.average.maxRT=2;
    EPmain.average.dropBad=1;
    EPmain.average.dropError=1;
    EPmain.average.dropTimeout=1;
    
    EPmain.transform.importFormat=EPmain.preferences.general.importFormat;
    EPmain.transform.type=3;
    EPmain.transform.outputFormat=EPmain.preferences.general.outputFormat;
    EPmain.transform.reference=EPmain.preferences.transform.reference;
    EPmain.transform.refChan1=EPmain.preferences.transform.refChan1;
    EPmain.transform.refChan2=EPmain.preferences.transform.refChan2;
    EPmain.transform.baselineStart=EPmain.preferences.transform.baselineStart;
    EPmain.transform.baselineEnd=EPmain.preferences.transform.baselineEnd;
    EPmain.transform.preStim=-EPmain.preferences.transform.baselineStart;
    EPmain.transform.delay=0;
    EPmain.transform.domain=1;
    EPmain.transform.method=1;
    EPmain.transform.smoothing=1;
    EPmain.transform.waveletwidth=7;
    EPmain.transform.filterPass=1;
    EPmain.transform.filterType=2;
    EPmain.transform.filter1=[];
    EPmain.transform.filter2=[];
    EPmain.transform.filterOrder=6;
    EPmain.transform.detrend=0;
    EPmain.transform.dataMode=1;
    EPmain.transform.whiteNoise=rand(1,250);
    
    EPmain.read.check=0;
    EPmain.read.subject='';
    EPmain.read.cell='';
    EPmain.read.freq='';
    EPmain.read.format=EPmain.preferences.general.importFormat;
    EPmain.read.type=3;   
    
    EPmain.preprocess.importFormat=EPmain.preferences.general.sessionImportFormat;
    EPmain.preprocess.outputFormat=EPmain.preferences.general.sessionOutputFormat;
    EPmain.preprocess.timepoints='';
    EPmain.preprocess.baseline='';
    EPmain.preprocess.detrend=0;
    EPmain.preprocess.EMG=0;
    EPmain.preprocess.alpha=0;
    EPmain.preprocess.blinkTemplate=1;
    EPmain.preprocess.saccadeTemplate=5;
    EPmain.preprocess.channelMode=1;
    EPmain.preprocess.trialMode=1;
    EPmain.preprocess.check=0;
    EPmain.preprocess.subject='';
    EPmain.preprocess.cell='';
    EPmain.preprocess.type=2;
    EPmain.preprocess.origReference='';
    EPmain.preprocess.origRefType=4;
    EPmain.preprocess.currReference='';
    EPmain.preprocess.currRefType=4;
    EPmain.preprocess.editMode=3;
    EPmain.preprocess.fMRI=0;
    
    EPmain.view=[];
    
    EPmain.pca.facNum=0;
    EPmain.pca.name='';
    EPmain.pca.mode=EPmain.preferences.pca.mode;
    EPmain.pca.rotation=EPmain.preferences.pca.rotation;
    EPmain.pca.rotopt=EPmain.preferences.pca.rotopt;
    EPmain.pca.rel=EPmain.preferences.pca.rel;
    EPmain.pca.loadings=EPmain.preferences.pca.loadings;
    EPmain.pca.parametric=false;
    EPmain.pca.crossVerifyPCA=[];
    EPscree=[];
    
    EPmain.sampleTest=[];
    
    EPmain.window.minFacVar=[];
    EPmain.window.FFTunits=4;
    EPmain.window.sampAdapt=0;
    EPmain.window.datasetName='';
    
    EPmain.anova.numComps=0;
    EPmain.anova.data=[];
    
    EPmain.save.format=EPmain.preferences.general.outputFormat;
    EPmain.save.SGLchan=1;
    EPmain.save.REGchan=1;
    EPmain.save.SGLcell=1;
    EPmain.save.CMBcell=1;
    EPmain.save.RAW=1;
    EPmain.save.AVG=1;
    EPmain.save.GAV=1;
    EPmain.save.SGLfac=1;
    EPmain.save.CMBfac=1;
    EPmain.save.batch=1;
    
    EPmain.save.check=0;
    EPmain.save.subject='';
    EPmain.save.cell='';
    EPmain.save.freq='';
    EPmain.save.readFormat=EPmain.preferences.general.importFormat;
    EPmain.save.type=3;   
    
    EPchanGrp=[];
    
    if err
        ep('savePrefs'); %save reset preference settings
    end;
end;

scrsz = EPmain.scrsz;

switch varargin{1}
    case 'start'
        if ~isempty(EPmain.handles.hMainWindow)
            clf(EPmain.handles.hMainWindow)
            figure(EPmain.handles.hMainWindow)
        else
            EPmain.handles.hMainWindow = figure('Name', 'EP Toolkit', 'NumberTitle', 'off', 'Position',[1 scrsz(4)-550 200 500], 'MenuBar', 'none');
            colormap jet;
            drawnow
        end;
        
        switch EPmain.mode
            case 'main'
                ep('startMain');
                
            case 'segment'
                ep('startSegment');
                
            case 'preprocess'
                ep('startPreprocess');
                
            case 'average'
                ep('startAverage');
                
            case 'transform'
                ep('startTransform');
                
            case 'read'
                ep('startRead');
                
            case 'edit'
                ep('startEdit');
                
            case 'view'
                ep('startView');
                
            case 'sampleTest'
                ep('startSampleTest');
                
            case 'PCA'
                ep('startPCA');
                
            case 'window'
                ep('startWindow');
                
            case 'ANOVA'
                ep('startANOVA');
                
            case 'save'
                ep('startSave');
                
            case 'preferenceMain'
                ep('startPreferenceMain');
                
            case 'preferenceGeneral'
                ep('startPreferenceGeneral');
                
            case 'preferencePreprocess'
                ep('startPreferencePreprocess');
                
            case 'preferenceAverage'
                ep('startPreferenceAverage');
                
            case 'preferenceTransform'
                ep('startPreferenceTransform');
                
            case 'preferenceView'
                ep('startPreferenceView');
                
            case 'preferencePCA'
                ep('startPreferencePCA');
                
            case 'preferenceWindow'
                ep('startPreferenceWindow');
                
            case 'preferenceANOVA'
                ep('startPreferenceANOVA');
                
            otherwise
                error('Not a valid mode.');
        end;
        
    case 'startMain'
        
        screeFigure=findobj('Name', 'ScreeWindow');
        if ~isempty(screeFigure)
            close(screeFigure)
        end;
        
        chanGrpFigure=findobj('Name', 'Channel Group Window');
        if ~isempty(chanGrpFigure)
            close(chanGrpFigure)
        end;
        
        editDataFigure=findobj('Name', 'EditData');
        if ~isempty(editDataFigure)
            close(editDataFigure)
        end;
        
        set(EPmain.handles.hMainWindow,'Name', 'EP Toolkit');
        
        EPmain.handles.hMenu = uimenu('Label','File');
        uimenu(EPmain.handles.hMenu,'Label','About EP Toolkit','Callback',@splash);
        uimenu(EPmain.handles.hMenu,'Label','Preferences','Callback',['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        uimenu(EPmain.handles.hMenu,'Label','Change Work Directory','Callback',@changeWork);
        uimenu(EPmain.handles.hMenu,'Label','Create Work Directory','Callback',@createWork);
        uimenu(EPmain.handles.hMenu,'Label','Clear Work Directory','Callback',@clearWorkingSet);
        uimenu(EPmain.handles.hMenu,'Label','Quit','Callback',@quit);
        
        EPmain.handles.main.segment = uicontrol('Style', 'pushbutton', 'String', 'Segment','FontSize',EPmain.fontsize,...
            'Position', [20 450 100 30], 'Callback', ['global EPmain;','EPmain.mode=''segment'';','ep(''start'');']);
        
        EPmain.handles.main.preprocess = uicontrol('Style', 'pushbutton', 'String', 'Preprocess','FontSize',EPmain.fontsize,...
            'Position', [20 410 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preprocess'';','ep(''start'');']);
        
        EPmain.handles.main.transform = uicontrol('Style', 'pushbutton', 'String', 'Transform','FontSize',EPmain.fontsize,...
            'Position', [20 370 100 30], 'Callback', ['global EPmain;','EPmain.mode=''transform'';','ep(''start'');']);
        
        EPmain.handles.main.average = uicontrol('Style', 'pushbutton', 'String', 'Average','FontSize',EPmain.fontsize,...
            'Position', [20 330 100 30], 'Callback', ['global EPmain;','EPmain.mode=''average'';','ep(''start'');']);
        
        EPmain.handles.main.read = uicontrol('Style', 'pushbutton', 'String', 'Read','FontSize',EPmain.fontsize,...
            'Position', [20 290 100 30], 'Callback', ['global EPmain;','EPmain.mode=''read'';','ep(''start'');']);
        
        EPmain.handles.main.edit = uicontrol('Style', 'pushbutton', 'String', 'Edit','FontSize',EPmain.fontsize,...
            'Position', [20 250 100 30], 'Callback', ['global EPmain;','EPmain.mode=''edit'';','ep(''start'');']);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.main.edit,'enable','off');
        end;
        
        EPmain.handles.main.view = uicontrol('Style', 'pushbutton', 'String', 'View','FontSize',EPmain.fontsize,...
            'Position', [20 210 100 30], 'Callback', ['global EPmain;','EPmain.mode=''view'';','ep(''start'');']);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.main.view,'enable','off');
        end;
        
        EPmain.handles.main.PCA = uicontrol('Style', 'pushbutton', 'String', 'PCA','FontSize',EPmain.fontsize,...
            'Position', [20 170 100 30], 'Callback', ['global EPmain;','EPmain.mode=''PCA'';','ep(''start'');']);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.main.PCA,'enable','off');
        end;
        
        EPmain.handles.main.sampleTest = uicontrol('Style', 'pushbutton', 'String', 'Sample Test','FontSize',EPmain.fontsize,...
            'Position', [20 130 100 30], 'Callback', ['global EPmain;','EPmain.mode=''sampleTest'';','ep(''start'');']);
        
        if ~ft_hastoolbox('STATS', 0, 1)
            set(EPmain.handles.main.sampleTest,'enable','off');
        end;
        
        EPmain.handles.main.window = uicontrol('Style', 'pushbutton', 'String', 'Window','FontSize',EPmain.fontsize,...
            'Position', [20 90 100 30], 'Callback', ['global EPmain;','EPmain.mode=''window'';','ep(''start'');']);
        
        if ~isempty(EPdataset.dataset)
            isAve=false;
            for i=1:length(EPdataset.dataset)
                if any(strcmp(EPdataset.dataset(i).dataType,{'average','single_trial'})) && any(ismember(EPdataset.dataset(i).chanTypes,{'EEG','MGM','MGA','MGP'})) && any(strcmp(EPdataset.dataset(i).cellTypes,'SGL')) && (any(strcmp(EPdataset.dataset(i).facTypes,'SGL')) || isempty(EPdataset.dataset(i).facNames))
                    isAve=true;
                end;
            end;
            if ~isAve
                set(EPmain.handles.main.window,'enable','off');
            end;
        else
            set(EPmain.handles.main.window,'enable','off');
        end;
        
        EPmain.handles.main.ANOVA = uicontrol('Style', 'pushbutton', 'String', 'ANOVA','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 30], 'Callback', ['global EPmain;','EPmain.mode=''ANOVA'';','ep(''start'');']);
        
        EPmain.handles.main.save = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
            'Position', [20 10 100 30], 'Callback', ['global EPmain;','EPmain.mode=''save'';','ep(''start'');']);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.main.save,'enable','off');
        end;
        
        
    case 'startSegment'
        
        set(EPmain.handles.hMainWindow,'Name', 'Segment Data');
        refresh
        
        EPmain.segment.numSpecs=6;
        
        if ~isfield(EPmain.segment,'cellTable') %just entering segment function
            EPmain.segment.contData=[];
            for iFile=1:length(EPdataset.dataset)
                if any(strcmp(EPdataset.dataset(iFile).dataType,{'continuous','single_trial'})) && ~isempty(EPdataset.dataset(iFile).events{1})
                    EPmain.segment.contData(end+1)=iFile;
                end;
            end;
            
            EPmain.segment.relList={'=','~=','<','>','<=','>=','starts','ends','contains'};
            EPmain.segment.delay=0;
            EPmain.segment.cellTable=cell(1,6+EPmain.segment.numSpecs*3);
            EPmain.segment.importFormat=EPmain.preferences.general.sessionImportFormat;
            EPmain.segment.outputFormat=EPmain.preferences.general.outputFormat;
            EPmain.segment.flexible=0;
            EPmain.segment.flexEnd=1;
            EPmain.segment.flexLength=20;
            if ~isempty(EPmain.segment.contData)
                changeSegmentDataset;
            end
        end;
        
        if isempty(EPmain.segment.contData)
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','No continuous or single-trial dataset in working set to serve as template.','Position',[5 460 150 40]);
        else
            EPmain.handles.segment.dataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{EPdataset.dataset(EPmain.segment.contData).dataName},...
                'Value',find(EPmain.segment.dataset==EPmain.segment.contData),'Position',[5 480 200 20],...
                'Callback', @changeSegmentDataset);
            
            EPmain.handles.segment.flexible = uicontrol('Style','popupmenu',...
                'String',{'Time event','Flex event'},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Value',EPmain.segment.flexible+1,'Position',[5 460 110 20],...
                'Callback',['global EPmain;','EPmain.segment.flexible=get(EPmain.handles.segment.flexible,''Value'')-1;','ep(''start'')']);
            
            if ~isempty(EPmain.segment.eventLocks)
                EPmain.handles.segment.eventValues = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.eventLocks,...
                    'Value',EPmain.segment.event,'Position',[105 460 100 20],...
                    'Callback', @changeSegmentEvent);
            else
                EPmain.handles.segment.eventValues = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [5 460 80 20],'enable','off');
            end;
            
            if ~EPmain.segment.flexible
                
                uicontrol('Style','text',...
                    'String','samples','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 440 50 20]);
                
                uicontrol('Style','text',...
                    'String','ms','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[75 440 50 20]);
                
                uicontrol('Style','text',...
                    'String','delay ms','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[145 440 50 20]);
                
                EPmain.handles.segment.sampStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.sampStart,...
                    'TooltipString','First sample of epoch in relation to the event, where negative is before it.',...
                    'Position',[5 420 35 20],'Callback',@segmentSampStart);
                
                EPmain.handles.segment.sampEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.sampEnd,...
                    'TooltipString','Last sample of epoch in relation to the event, where negative is before it.',...
                    'Position',[40 420 35 20],'Callback',@segmentSampEnd);
                
                EPmain.handles.segment.msStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',round((EPmain.segment.sampStart)*(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)),...
                    'TooltipString','Last ms of epoch in relation to the event, where negative is before it.',...
                    'Position',[75 420 35 20],'Callback',@segmentSampStart);
                
                EPmain.handles.segment.msEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',round((EPmain.segment.sampEnd)*(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)),...
                    'TooltipString','First ms of epoch in relation to the event, where negative is before it.',...
                    'Position',[110 420 35 20],'Callback',@segmentSampEnd);
            else
                
                uicontrol('Style','text',...
                    'String','Flex end event','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 440 100 20]);
                
                if ~isempty(EPmain.segment.eventLocks)
                    EPmain.handles.segment.flexEnd = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.eventLocks,...
                        'Value',EPmain.segment.flexEnd,'Position',[105 440 100 20],...
                        'Callback',['global EPmain;','EPmain.segment.flexEnd=get(EPmain.handles.segment.flexEnd,''Value'');','refresh']);
                else
                    EPmain.handles.segment.flexEnd = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [5 440 80 20],'enable','off');
                end;
                
                uicontrol('Style','text',...
                    'String','Flex samp','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 420 50 20]);
                
                EPmain.handles.segment.flexLength= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',sprintf('%d',EPmain.segment.flexLength),...
                    'TooltipString','Number of samples in the flexible length segment.',...
                    'Position',[60 420 35 20],...
                    'Callback',['global EPmain;','EPmain.segment.flexLength=str2num(get(EPmain.handles.segment.flexLength,''String''));','refresh']);
                
                uicontrol('Style','text',...
                    'String','delay ms','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[95 420 50 20]);
                
            end;
            
            EPmain.handles.segment.delay= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',sprintf('%d',EPmain.segment.delay),...
                'TooltipString','Ms correction for nominal event timing (positive means delay).',...
                'Position',[145 420 35 20],...
                'Callback',['global EPmain;','EPmain.segment.delay=str2num(get(EPmain.handles.segment.delay,''String''));','refresh']);
            
            if ~isempty(EPmain.segment.critSpecNames)
                EPmain.handles.segment.trialSpecNames(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemNames{1},...
                    'Value',EPmain.segment.trialSpec(1),'Position',[5 400 80 20],...
                    'Callback',@changeSegmentCritName);
                EPmain.handles.segment.trialSpecRel(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.relList,...
                    'Value',EPmain.segment.trialSpecRel(1),'Position',[75 400 70 20],...
                    'Callback',['global EPmain;','EPmain.segment.trialSpecRel(1)=get(EPmain.handles.segment.trialSpecRel(1),''Value'');','refresh']);
                EPmain.handles.segment.trialSpecVal(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemValues{1},...
                    'Value',EPmain.segment.trialSpecVal(1),'Position',[135 400 70 20],...
                    'Callback',@changeSegmentCritValue);
            else
                EPmain.handles.segment.trialSpecNames(1) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [5 400 80 20],'enable','off');
                EPmain.handles.segment.trialSpecRel(1) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(1)},'FontSize',EPmain.fontsize,...
                    'Position', [75 400 70 20],'enable','off');
                EPmain.handles.segment.trialSpecVal(1) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [135 400 70 20],'enable','off');
            end;
            
            if ~isempty(EPmain.segment.critSpecNames)
                EPmain.handles.segment.trialSpecNames(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemNames{2},...
                    'Value',EPmain.segment.trialSpec(2),'Position',[5 380 80 20],...
                    'Callback',@changeSegmentCritName);
                EPmain.handles.segment.trialSpecRel(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.relList,...
                    'Value',EPmain.segment.trialSpecRel(2),'Position',[75 380 70 20],...
                    'Callback',['global EPmain;','EPmain.segment.trialSpecRel(2)=get(EPmain.handles.segment.trialSpecRel(2),''Value'');','refresh']);
                EPmain.handles.segment.trialSpecVal(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemValues{2},...
                    'Value',EPmain.segment.trialSpecVal(2),'Position',[135 380 70 20],...
                    'Callback',@changeSegmentCritValue);
            else
                EPmain.handles.segment.trialSpecNames(2) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [5 380 80 20],'enable','off');
                EPmain.handles.segment.trialSpecRel(2) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(2)},'FontSize',EPmain.fontsize,...
                    'Position', [75 380 70 20],'enable','off');
                EPmain.handles.segment.trialSpecVal(2) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [135 380 70 20],'enable','off');
            end;
            
            if ~isempty(EPmain.segment.critSpecNames)
                EPmain.handles.segment.trialSpecNames(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemNames{3},...
                    'Value',EPmain.segment.trialSpec(3),'Position',[5 360 80 20],...
                    'Callback',@changeSegmentCritName);
                EPmain.handles.segment.trialSpecRel(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.relList,...
                    'Value',EPmain.segment.trialSpecRel(3),'Position',[75 360 70 20],...
                    'Callback',['global EPmain;','EPmain.segment.trialSpecRel(3)=get(EPmain.handles.segment.trialSpecRel(3),''Value'');','refresh']);
                EPmain.handles.segment.trialSpecVal(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemValues{3},...
                    'Value',EPmain.segment.trialSpecVal(3),'Position',[135 360 70 20],...
                    'Callback',@changeSegmentCritValue);
                
            else
                EPmain.handles.segment.trialSpecNames(3) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [5 360 80 20],'enable','off');
                EPmain.handles.segment.trialSpecRel(3) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(3)},'FontSize',EPmain.fontsize,...
                    'Position', [75 360 70 20],'enable','off');
                EPmain.handles.segment.trialSpecVal(3) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [135 360 70 20],'enable','off');
            end;
            
            if ~isempty(EPmain.segment.critSpecNames)
                EPmain.handles.segment.trialSpecNames(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemNames{4},...
                    'Value',EPmain.segment.trialSpec(4),'Position',[5 340 80 20],...
                    'Callback',@changeSegmentCritName);
                EPmain.handles.segment.trialSpecRel(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.relList,...
                    'Value',EPmain.segment.trialSpecRel(4),'Position',[75 340 70 20],...
                    'Callback',['global EPmain;','EPmain.segment.trialSpecRel(4)=get(EPmain.handles.segment.trialSpecRel(4),''Value'');','refresh']);
                EPmain.handles.segment.trialSpecVal(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemValues{4},...
                    'Value',EPmain.segment.trialSpecVal(4),'Position',[135 340 70 20],...
                    'Callback',@changeSegmentCritValue);
                
            else
                EPmain.handles.segment.trialSpecNames(4) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [5 340 80 20],'enable','off');
                EPmain.handles.segment.trialSpecRel(4) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(4)},'FontSize',EPmain.fontsize,...
                    'Position', [75 340 70 20],'enable','off');
                EPmain.handles.segment.trialSpecVal(4) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [135 340 70 20],'enable','off');
            end;
            
            if ~isempty(EPmain.segment.critSpecNames)
                EPmain.handles.segment.trialSpecNames(5) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemNames{5},...
                    'Value',EPmain.segment.trialSpec(5),'Position',[5 320 80 20],...
                    'Callback',@changeSegmentCritName);
                EPmain.handles.segment.trialSpecRel(5) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.relList,...
                    'Value',EPmain.segment.trialSpecRel(5),'Position',[75 320 70 20],...
                    'Callback',['global EPmain;','EPmain.segment.trialSpecRel(5)=get(EPmain.handles.segment.trialSpecRel(5),''Value'');','refresh']);
                EPmain.handles.segment.trialSpecVal(5) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemValues{5},...
                    'Value',EPmain.segment.trialSpecVal(5),'Position',[135 320 70 20],...
                    'Callback',@changeSegmentCritValue);
                
            else
                EPmain.handles.segment.trialSpecNames(5) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [5 320 80 20],'enable','off');
                EPmain.handles.segment.trialSpecRel(5) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(5)},'FontSize',EPmain.fontsize,...
                    'Position', [75 320 70 20],'enable','off');
                EPmain.handles.segment.trialSpecVal(5) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [135 320 70 20],'enable','off');
            end;
            
            if ~isempty(EPmain.segment.critSpecNames)
                EPmain.handles.segment.trialSpecNames(6) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemNames{6},...
                    'Value',EPmain.segment.trialSpec(6),'Position',[5 300 80 20],...
                    'Callback',@changeSegmentCritName);
                EPmain.handles.segment.trialSpecRel(6) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.relList,...
                    'Value',EPmain.segment.trialSpecRel(6),'Position',[75 300 70 20],...
                    'Callback',['global EPmain;','EPmain.segment.trialSpecRel(6)=get(EPmain.handles.segment.trialSpecRel(6),''Value'');','refresh']);
                EPmain.handles.segment.trialSpecVal(6) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.critSpecItemValues{6},...
                    'Value',EPmain.segment.trialSpecVal(6),'Position',[135 300 70 20],...
                    'Callback',@changeSegmentCritValue);
                
            else
                EPmain.handles.segment.trialSpecNames(6) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [5 300 80 20],'enable','off');
                EPmain.handles.segment.trialSpecRel(6) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(6)},'FontSize',EPmain.fontsize,...
                    'Position', [75 300 70 20],'enable','off');
                EPmain.handles.segment.trialSpecVal(6) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [135 300 70 20],'enable','off');
            end;
            
            EPmain.handles.segment.addCell = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*2,...
                'Position', [10 270 40 25], 'Callback', @segmentAddCell);
            
            EPmain.handles.segment.delCell = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*2,...
                'Position', [55 270 40 25], 'Callback', @segmentDelCell);
            
        end;

        EPmain.handles.segment.preview = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
            'Position', [100 270 40 25], 'Callback', @segmentLoad);
        
        EPmain.handles.segment.preview = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
            'Position', [145 270 40 25], 'Callback', @segmentSave);
        
        tableData=EPmain.segment.cellTable;
        
        tableNames{1}='#';
        tableNames{2}='name';
        tableNames{3}='stim';
        if ~EPmain.segment.flexible
            tableNames{4}='prestim';
            tableNames{5}='poststim';
        else
            tableNames{4}='flexEnd';
            tableNames{5}='length';
        end;
        tableNames{6}='delay';
        
        for i=1:EPmain.segment.numSpecs
            tableNames{4+i*3}=['spec' num2str(i)];
            tableNames{5+i*3}=['rel' num2str(i)];
            tableNames{6+i*3}=['value' num2str(i)];
        end;
        
        columnEditable =  [false true(1,5+EPmain.segment.numSpecs*3)];
        ColumnFormat{1}='numeric';
        ColumnFormat{2}='char';
        ColumnFormat{3}='char';
        ColumnFormat{4}='char';
        ColumnFormat{5}='char';
        ColumnFormat{6}='numeric';
        for i=1:EPmain.segment.numSpecs
            ColumnFormat{4+i*3}='char';
            ColumnFormat{5+i*3}='char';
            ColumnFormat{6+i*3}='char';
        end;
        ColumnWidth=num2cell(repmat(50,1,6+EPmain.segment.numSpecs*3));
        
        EPmain.handles.segment.cellTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,'ColumnWidth',ColumnWidth,...
            'RearrangeableColumns','on',...
            'CellEditCallback',['global EPmain;','EPmain.segment.cellTable=get(EPmain.handles.segment.cellTable,''Data'');','ep(''start'');'],'Position',[5 120 180 150]);
        
        uicontrol('Style','text',...
            'String','Import File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 100 150 20]);
        
        EPmain.handles.segment.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.importFormat,''Value'');','if tempVar ~=0,EPmain.segment.importFormat=tempVar;end;','if isempty(tempVar),EPmain.segment.importFormat=tempVar;end;'],...
            'Value',EPmain.segment.importFormat,'Position',[20 80 160 20]);
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Output File Format','HorizontalAlignment','left',...
            'Position',[25 60 150 20]);
        
        EPmain.handles.segment.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.outputFormat,''Value'');','if tempVar ~=0,EPmain.segment.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.segment.outputFormat=tempVar;end;'],...
            'Value',EPmain.segment.outputFormat,'Position',[20 40 160 20]);
        
        EPmain.handles.segment.trim = uicontrol('Style', 'pushbutton', 'String', 'Trim','FontSize',EPmain.fontsize,...
            'Position', [2 0 50 35], 'Callback', ['global EPtrimData;','EPtrimData=[];','ep_trimData']);
        
        if isempty(EPmain.segment.contData) || ~isempty(EPdataset.dataset(EPmain.segment.dataset).freqNames) || ~isempty(EPdataset.dataset(EPmain.segment.dataset).facNames)
            set(EPmain.handles.segment.trim,'enable','off');
        end;

        EPmain.handles.segment.preview = uicontrol('Style', 'pushbutton', 'String', 'Preview','FontSize',EPmain.fontsize,...
            'Position', [52 0 50 35], 'Callback', ['EPmain.segment.preview=1;','ep(''segmentData'')']);
        
        if isempty(EPmain.segment.contData)
            set(EPmain.handles.segment.preview,'enable','off');
        end

        EPmain.handles.segment.segment = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
            'Position', [102 0 50 35], 'Callback', ['EPmain.segment.preview=0;','ep(''segmentData'')']);
        
        EPmain.handles.segment.main = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [152 0 50 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','EPmain.segment=[];','ep(''start'');']);
        
        if any(any(cellfun(@isempty,EPmain.segment.cellTable(:,2:4))))
            set(EPmain.handles.segment.preview,'enable','off');
            set(EPmain.handles.segment.segment,'enable','off');
        end;
        
    case 'segmentData'
        
        numSpecs=6;
        %check cell table
         EPmain.segment.cellTable(:,2)=cellfun(@deblank,EPmain.segment.cellTable(:,2),'UniformOutput',false);
        
        for iCell=1:size(EPmain.segment.cellTable,1)
            if iCell==1
                flexMode=strcmpi(EPmain.segment.cellTable{iCell,5}(1),'F');
            elseif flexMode && ~strcmpi(EPmain.segment.cellTable{iCell,5}(1),'F')
                msg{1}='Error: all cells need to be uniformly either flexmode or not.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            if flexMode
                flexLength=str2num(EPmain.segment.cellTable{iCell,5}(2:end));
                if isempty(flexLength) || isnan(flexLength)
                    msg{1}='Error: Flex mode requires a number in addition to the F prefix, as in F20.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if flexLength < 1
                    msg{1}='Error: Flex length needs to be at least one.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if floor(flexLength) ~= flexLength
                    msg{1}='Error: Flex length needs to be an integer.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
            else
                startTime=str2num(EPmain.segment.cellTable{iCell,4});
                endTime=str2num(EPmain.segment.cellTable{iCell,5});
                if isempty(startTime) || isnan(startTime)
                    msg{1}='Error: prestim value is not a number.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if isempty(endTime) || isnan(endTime)
                    msg{1}='Error: poststim value is not a number.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if startTime > endTime
                    msg{1}='Error: cell has starting sample after ending sample.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if (startTime == endTime) && ~EPmain.segment.flexible
                    msg{1}='Error: cell has a zero segment length.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if iCell==1
                    epochLength=endTime-startTime;
                elseif epochLength~=(endTime-startTime)
                    msg{1}='Error: cells do not have same epoch length.';
                    [msg]=ep_errorMsg(msg);
                    return
                end;
            end;
            
            delayTime=str2num(EPmain.segment.cellTable{iCell,6});
            if isempty(delayTime) || isnan(delayTime)
                msg{1}='Error: delay value is not a number.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            for iSpec=1:numSpecs
                if strcmp('none',EPmain.segment.cellTable{iCell,4+iSpec*3})
                    EPmain.segment.cellTable{iCell,4+iSpec*3}='';
                end;
                if isempty(strcmp(EPmain.segment.cellTable{iCell,5+iSpec*3},EPmain.segment.relList)) && ~isempty(EPmain.segment.cellTable{iCell,5+iSpec*3})
                    msg{1}=['Error: spec relationship (' EPmain.segment.cellTable{iCell,5+iSpec*3} ') does not match list of options.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                if strcmp('none',EPmain.segment.cellTable{iCell,6+iSpec*3})
                    EPmain.segment.cellTable{iCell,6+iSpec*3}='';
                end;
                if any(strcmp(EPmain.segment.cellTable{iCell,4+iSpec*3},{'-precedes-','-follows-'}))
                    if  iSpec ~= numSpecs
                        if any(strcmp(EPmain.segment.cellTable{iCell,4+(iSpec+1)*3},{'-precedes-','-follows-'}))
                            msg{1}='Error: Consecutive ''-precedes-'' and ''-follows-'' keywords are not allowed.  The criterion following a ''-precedes-'' or ''-follows-'' keyword must be used to finish its specification.';
                            [msg]=ep_errorMsg(msg);
                            return
                        end;
                    else
                        msg{1}='Error: The final criterion cannot use the ''-precedes-'' or ''-follows-'' keywords as a subsequent criterion is needed to finish their specifications.';
                        [msg]=ep_errorMsg(msg);
                        return
                    end;
                end;
            end;
        end;
        
        set(EPmain.handles.segment.trim,'enable','off');
        set(EPmain.handles.segment.preview,'enable','off');
        set(EPmain.handles.segment.segment,'enable','off');
        set(EPmain.handles.segment.main,'enable','off');
        drawnow

        if EPmain.segment.preview
            sessionFiles=EPdataset.dataset(EPmain.segment.dataset);
            importFormat='ep_mat';
            outputFormat=[];
        else
            [importSuffix,importFormatName,importFormat]=ep_fileFormats('continuous',EPmain.fileFormatReadList{EPmain.segment.importFormat});
            [outputSuffix,outputFormatName,outputFormat]=ep_fileFormats('single_trial',EPmain.fileFormatSaveList{EPmain.segment.outputFormat});
            [sessionFiles, activeDirectory]=ep_getFilesUI(importFormat);
            if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
                msg{1}='No filenames selected. You have to click on a name.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.segment.trim,'enable','on');
                set(EPmain.handles.segment.preview,'enable','on');
                set(EPmain.handles.segment.segment,'enable','on');
                set(EPmain.handles.segment.main,'enable','on');
                return
            end
            for iFile=1:size(sessionFiles,2)
                sessionFiles{iFile}=[activeDirectory sessionFiles{iFile}];
            end;
        end;
        
        cellNums=ep_segmentData(EPmain.segment.cellTable(:,2:end),sessionFiles,importFormat,outputFormat,EPmain.segment.preview);
        set(EPmain.handles.segment.trim,'enable','on');
        set(EPmain.handles.segment.preview,'enable','on');
        set(EPmain.handles.segment.segment,'enable','on');
        set(EPmain.handles.segment.main,'enable','on');        
        
        if ~isempty(cellNums)
            EPmain.segment.cellTable=get(EPmain.handles.segment.cellTable,'Data');
            if size(cellNums,1) ==1
                theSums=cellNums;
            else
                theSums=sum(cellNums);
            end;
            EPmain.segment.cellTable(:,1)=num2cell(theSums);
            set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
            ep('start');
        end;
        
    case 'startPreprocess'
        
        set(EPmain.handles.hMainWindow,'Name', 'Preprocess Data');
        
        uicontrol('Style','text',...
            'String','Import File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 470 150 20]);
        
        EPmain.handles.preprocess.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.importFormat,''Value'');','if tempVar ~=0,EPmain.preprocess.importFormat=tempVar;end;','if isempty(tempVar),EPmain.preprocess.importFormat=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.importFormat,'Position',[20 450 160 20]);
        
        uicontrol('Style','text',...
            'String','File Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 430 150 20]);
        
        EPmain.handles.preprocess.type= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.type,''Value'');','if tempVar ~=0,EPmain.preprocess.type=tempVar;end;','if isempty(tempVar),EPmain.preprocess.type=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.type,'Position',[20 410 160 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.preprocess.importFormat});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.preprocess.type,'enable','off');
        end;
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Output File Format','HorizontalAlignment','left',...
            'Position',[25 390 150 20]);
        
        EPmain.handles.preprocess.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.outputFormat,''Value'');','if tempVar ~=0,EPmain.preprocess.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.preprocess.outputFormat=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.outputFormat,'Position',[20 370 160 20]);
        
        EPmain.handles.preprocess.timepointsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Points','FontSize',EPmain.fontsize,...
            'Position',[25 350 50 20]);
        
        if isempty(EPmain.preprocess.timepoints)
            set(EPmain.handles.preprocess.timepointsLabel,'enable','off');
        elseif isempty(str2num(EPmain.preprocess.timepoints))
            set(EPmain.handles.preprocess.timepointsLabel,'enable','off');
        end;
        
        EPmain.handles.preprocess.timepoints = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.timepoints,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.timepoints,''String'');','if tempVar ~=0,EPmain.preprocess.timepoints=tempVar;end;','if isempty(tempVar),EPmain.preprocess.timepoints=tempVar;end;','ep(''start'');'],...
            'Position',[25 330 50 20],'TooltipString','For retaining a subset of timepoints - example 1:250');
        
        EPmain.handles.preprocess.fMRI= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','fMRI',...
            'CallBack',['global EPmain;','EPmain.preprocess.fMRI=get(EPmain.handles.preprocess.fMRI,''Value'');','ep(''start'');'],...
            'Value',EPmain.preprocess.fMRI,'Position',[130 350 65 20],'TooltipString','Gradient and BCG artifact correction.  Requires signal processing toolbox.');
        
        EPmain.handles.preprocess.detrend= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Detrend',...
            'CallBack',['global EPmain;','EPmain.preprocess.detrend=get(EPmain.handles.preprocess.detrend,''Value'');','ep(''start'');'],...
            'Value',EPmain.preprocess.detrend,'Position',[130 330 65 20],'TooltipString','Not recommended for segmented ERP trials as the ERP components would be attenuated.');
        
        EPmain.handles.preprocess.EMG= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','EMG',...
            'CallBack',['global EPmain;','EPmain.preprocess.EMG=get(EPmain.handles.preprocess.EMG,''Value'');','ep(''start'');'],...
            'Value',EPmain.preprocess.EMG,'Position',[130 310 65 20],'TooltipString','Removes EMG using BSS-CCA.');
        
        EPmain.handles.preprocess.alpha= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','alpha',...
            'CallBack',['global EPmain;','EPmain.preprocess.alpha=get(EPmain.handles.preprocess.alpha,''Value'');','ep(''start'');'],...
            'Value',EPmain.preprocess.alpha,'Position',[130 290 65 20],'TooltipString','Removes alpha using CWT-PCA.');
        
        EPmain.handles.preprocess.baselineLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline','FontSize',EPmain.fontsize,...
            'Position',[75 350 50 20]);
        
        if isempty(EPmain.preprocess.baseline)
            set(EPmain.handles.preprocess.baselineLabel,'enable','off');
        elseif isempty(str2num(EPmain.preprocess.baseline))
            set(EPmain.handles.preprocess.baselineLabel,'enable','off');
        end;
        
        EPmain.handles.preprocess.baseline = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.baseline,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.baseline,''String'');','if tempVar ~=0,EPmain.preprocess.baseline=tempVar;end;','if isempty(tempVar),EPmain.preprocess.baseline=tempVar;end;','ep(''start'');'],...
            'Position',[75 330 50 20],'TooltipString','Recommended - example 1:50 (if timepoints were dropped, 1st retained sample is sample #1)');
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Edit Mode','HorizontalAlignment','left',...
            'Position',[25 290 100 20]);
        
        EPmain.handles.preprocess.editMode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'auto','manual','both'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.editMode,''Value'');','if tempVar ~=0,EPmain.preprocess.editMode=tempVar;end;','if isempty(tempVar),EPmain.preprocess.editMode=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.editMode,'Position',[20 270 160 20],'TooltipString','Whether bad channel and trial correction is based on manual, automatic, or both types of criteria.');
        
        uicontrol('Style','text',...
            'String','Blinks','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 250 70 20]);
        
        EPmain.handles.preprocess.blinkTemplate = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'file','auto','both','eye-track','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.blinkTemplate,''Value'');','if tempVar ~=0,EPmain.preprocess.blinkTemplate=tempVar;end;','if isempty(tempVar),EPmain.preprocess.blinkTemplate=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.blinkTemplate,'Position',[20 230 90 20]);
        
        uicontrol('Style','text',...
            'String','Saccades','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[105 250 70 20]);
        
        EPmain.handles.preprocess.saccadeTemplate = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'file','auto','both','eye-track','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.saccadeTemplate,''Value'');','if tempVar ~=0,EPmain.preprocess.saccadeTemplate=tempVar;end;','if isempty(tempVar),EPmain.preprocess.saccadeTemplate=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.saccadeTemplate,'Position',[100 230 90 20]);
        
        uicontrol('Style','text',...
            'String','Bad Channel Correction','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 210 150 20]);
        
        EPmain.handles.preprocess.channelMode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'replace','mark','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.channelMode,''Value'');','if tempVar ~=0,EPmain.preprocess.channelMode=tempVar;end;','if isempty(tempVar),EPmain.preprocess.channelMode=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.channelMode,'Position',[20 190 160 20]);
        
        uicontrol('Style','text',...
            'String','Movement Correction','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 170 150 20]);
        
        EPmain.handles.preprocess.trialMode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'fix','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.trialMode,''Value'');','if tempVar ~=0,EPmain.preprocess.trialMode=tempVar;end;','if isempty(tempVar),EPmain.preprocess.trialMode=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.trialMode,'Position',[20 150 160 20]);
        
        uicontrol('Style','text',...
            'String','O','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 130 10 20]);
        
        EPmain.handles.preprocess.origRefType = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',refList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.origRefType,''Value'');','if tempVar ~=0,EPmain.preprocess.origRefType=tempVar;end;','if isempty(tempVar),EPmain.preprocess.origRefType=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.origRefType,'Position',[35 130 100 20],'TooltipString','Original reference channels explicitly represented as a waveform in the data file.');
        
        EPmain.handles.preprocess.origReference = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.origReference,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.origReference,''String'');','if tempVar ~=0,EPmain.preprocess.origReference=tempVar;end;','if isempty(tempVar),EPmain.preprocess.origReference=tempVar;end;','ep(''start'');'],...
            'Position',[130 130 45 20],'TooltipString','example: 20 21');
        
        if EPmain.preprocess.origRefType ~= 2
            set(EPmain.handles.preprocess.origReference,'enable','off');
        end;

        uicontrol('Style','text',...
            'String','C','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 110 10 20]);
        
        EPmain.handles.preprocess.currRefType = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',refList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.currRefType,''Value'');','if tempVar ~=0,EPmain.preprocess.currRefType=tempVar;end;','if isempty(tempVar),EPmain.preprocess.currRefType=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.preprocess.currRefType,'Position',[35 110 100 20],'TooltipString','Current reference channels explicitly represented as a waveform in the data file.');
        
        EPmain.handles.preprocess.currReference = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.currReference,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.currReference,''String'');','if tempVar ~=0,EPmain.preprocess.currReference=tempVar;end;','if isempty(tempVar),EPmain.preprocess.currReference=tempVar;end;','ep(''start'');'],...
            'Position',[130 110 45 20],'TooltipString','example: 20 21');
        
        if EPmain.preprocess.currRefType ~= 2
            set(EPmain.handles.preprocess.currReference,'enable','off');
        end;

        uicontrol('Style','frame',...
            'Position',[20 35 170 70]);
        
        EPmain.handles.preprocess.check= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Single Cell Files',...
            'CallBack',['global EPmain;','EPmain.preprocess.check=get(EPmain.handles.preprocess.check,''Value'');','ep(''start'');'],...
            'Value',EPmain.preprocess.check,'Position',[30 80 100 20]);
        
        EPmain.handles.preprocess.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preprocess.subject,...
            'CallBack',['global EPmain;','EPmain.preprocess.subject=get(EPmain.handles.preprocess.subject,''String'');','ep(''start'');'],...
            'Position',[30 60 50 20],'TooltipString','example 4:6');
        
        EPmain.handles.preprocess.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Subject','HorizontalAlignment','left',...
            'Position',[90 60 50 20]);

        
        if isempty(EPmain.preprocess.subject)
            set(EPmain.handles.preprocess.subjectLabel,'enable','off');
        elseif isempty(str2num(EPmain.preprocess.subject))
            set(EPmain.handles.preprocess.subjectLabel,'enable','off');
        end;
        
        EPmain.handles.preprocess.cell= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preprocess.cell,...
            'CallBack',['global EPmain;','EPmain.preprocess.cell=get(EPmain.handles.preprocess.cell,''String'');','ep(''start'');'],...
            'Position',[30 40 50 20],'TooltipString','example 7:9');
        
        EPmain.handles.preprocess.cellLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Cell','HorizontalAlignment','left',...
            'Position',[90 40 50 20]);
        
        if isempty(EPmain.preprocess.cell)
            set(EPmain.handles.preprocess.cellLabel,'enable','off');
        elseif isempty(str2num(EPmain.preprocess.cell))
            set(EPmain.handles.preprocess.cellLabel,'enable','off');
        end;
        
        EPmain.handles.preprocess.template = uicontrol('Style', 'pushbutton', 'String', 'Template','FontSize',EPmain.fontsize,...
            'Position', [10 0 60 35], 'Callback', ['global EPtemplate;','EPtemplate=[];','ep_template']);
        
        noData=true;
        for theData=1:length(EPdataset.dataset)
            EEGchans=find(strcmp('EEG',EPdataset.dataset(theData).chanTypes));
            if any(strcmp(EPdataset.dataset(theData).dataType,{'single_trial','continuous'}))
                if ~isempty(EPdataset.dataset(theData).eloc) && (length([EPdataset.dataset(theData).eloc(EEGchans).theta]) == length(EEGchans)) && (length([EPdataset.dataset(theData).eloc(EEGchans).radius]) == length(EEGchans))
                    noData=false;
                end;
            end;
        end;
        if noData
            set(EPmain.handles.preprocess.template,'enable','off');
            if length(EPdataset.dataset) > 1
                disp(['Could not use any of the datasets to form templates as they are either not raw data or they lack one or more electrode coordinates.']);
            end;
        end;
        
        EPmain.handles.preprocess.preprocess = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
            'Position', [70 0 60 35], 'Callback', 'ep(''preprocessData'')');
        
        EPmain.handles.preprocess.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'preprocessData' %start preprocessing the data
        
        set(EPmain.handles.preprocess.subject,'ForegroundColor','black');
        set(EPmain.handles.preprocess.cell,'ForegroundColor','black');
        
        set(EPmain.handles.preprocess.preprocess,'enable','off');
        set(EPmain.handles.preprocess.done,'enable','off');
        set(EPmain.handles.preprocess.template,'enable','off');
        drawnow
        
        textPrefs.firstRow=EPmain.preferences.general.firstRow;
        textPrefs.lastRow=EPmain.preferences.general.lastRow;
        textPrefs.firstCol=EPmain.preferences.general.firstCol;
        textPrefs.lastCol=EPmain.preferences.general.lastCol;
        textPrefs.orientation=EPmain.preferences.general.orientation;
        
        typeNum = get(EPmain.handles.preprocess.type,'value');
        switch typeNum
            case 1
                dataType='continuous';
            case 2
                dataType='single_trial';
            case 3
                dataType='average';
            case 4
                dataType='grand_average';
            case 5
                dataType='factors';
        end;
        
        importFormatNum = get(EPmain.handles.preprocess.importFormat,'value');        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{importFormatNum});

        
        outputFormatNum = get(EPmain.handles.preprocess.outputFormat,'value');
        [outputSuffix,outputFormatName,outputFormat]=ep_fileFormats(dataType,EPmain.fileFormatSaveList{outputFormatNum});
        
        EPmain.preprocess.format=importFormatNum;
        EPmain.preprocess.type=typeNum;
        
        timePoints=str2num(get(EPmain.handles.preprocess.timepoints,'string'));
        baseline=str2num(get(EPmain.handles.preprocess.baseline,'string'));
        detrend=get(EPmain.handles.preprocess.detrend,'value');
        fMRI=get(EPmain.handles.preprocess.fMRI,'value');
        if fMRI
            fMRI=EPmain.preferences.preprocess.fMRI;
        end;

        switch get(EPmain.handles.preprocess.editMode,'value')
            case 1
                editMode='automatic';
            case 2
                editMode='manual';
            case 3
                editMode='both';
        end;
                
        switch get(EPmain.handles.preprocess.blinkTemplate,'value')
            case 1
                blinkTemplate='fileTemplate';
            case 2
                blinkTemplate='autoTemplate';
            case 3
                blinkTemplate='bothTemplate';
            case 4
                blinkTemplate='eyeTracker';
            case 5
                blinkTemplate='none';
        end;
        switch get(EPmain.handles.preprocess.saccadeTemplate,'value')
            case 1
                saccadeTemplate='fileTemplate';
            case 2
                saccadeTemplate='autoTemplate';
            case 3
                saccadeTemplate='bothTemplate';
            case 4
                saccadeTemplate='eyeTracker';
            case 5
                saccadeTemplate='none';
        end;
        switch get(EPmain.handles.preprocess.channelMode,'value')
            case 1
                channelMode='replace';
            case 2
                channelMode='mark';
            case 3
                channelMode='none';
        end;
        switch get(EPmain.handles.preprocess.trialMode,'value')
            case 1
                trialMode='fix';
            case 2
                trialMode='none';
        end;
        
        [sessionFiles, activeDirectory]=ep_getFilesUI(importFormat);
        if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
            msg{1}='No filenames selected. You have to click on a name.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.preprocess.preprocess,'enable','on');
            set(EPmain.handles.preprocess.done,'enable','on');
            set(EPmain.handles.preprocess.template,'enable','on');
            return
        end
        
        for theFile=1:size(sessionFiles,2)
            sessionFiles{theFile}=[activeDirectory sessionFiles{theFile}];
        end;
        
        sessionFiles=sort(sessionFiles);
        
        if any(strcmp(blinkTemplate, {'fileTemplate','bothTemplate'}))
            if ismac
                h=msgbox('Please select blinks template file.'); %workaround for Matlab bug
            end;
            [FileName,PathName,FilterIndex] = uigetfile('*.mat','Blink Template','blinks.mat');
            if ismac
                close(h);
            end;
            if FileName == 0
                warndlg('No blink template file selected.');
                set(EPmain.handles.preprocess.preprocess,'enable','on');
                set(EPmain.handles.preprocess.done,'enable','on');
                set(EPmain.handles.preprocess.template,'enable','on');
                return;
            end;
            blinkFile=[PathName FileName];
        else
            blinkFile=[];
        end;
        
        if any(strcmp(saccadeTemplate, {'fileTemplate','bothTemplate'}))
            if ismac
                h=msgbox('Please select saccade template file.'); %workaround for Matlab bug
            end;
            [FileName,PathName,FilterIndex] = uigetfile('*.mat','Saccade Template','saccades.mat');
            if ismac
                close(h);
            end;
            if FileName == 0
                warndlg('No saccade template file selected.');
                set(EPmain.handles.preprocess.preprocess,'enable','on');
                set(EPmain.handles.preprocess.done,'enable','on');
                set(EPmain.handles.preprocess.template,'enable','on');
                return;
            end;
            saccadeFile=[PathName FileName];
        else
            saccadeFile=[];
        end;
        
        inArg=[];
        inArg{1}='files';
        inArg{2}=sessionFiles;
        inArg{3}='format';
        inArg{4}=importFormat;
        inArg{5}='type';
        inArg{6}=dataType;
        inArg{7}='outputFormat';
        inArg{8}=outputFormat;
        inArg{9}='template';
        inArg{10}=blinkTemplate;
        inArg{11}='channelMode';
        inArg{12}=channelMode;
        inArg{13}='saturation';
        inArg{14}=[-EPmain.preferences.preprocess.saturation EPmain.preferences.preprocess.saturation];
        inArg{15}='window';
        inArg{16}=EPmain.preferences.preprocess.window;
        inArg{17}='minmax';
        inArg{18}=EPmain.preferences.preprocess.minmax;
        inArg{19}='badnum';
        inArg{20}=EPmain.preferences.preprocess.badnum;
        inArg{21}='neighbors';
        inArg{22}=EPmain.preferences.preprocess.neighbors;
        inArg{23}='maxneighbor';
        inArg{24}=EPmain.preferences.preprocess.maxneighbor;
        inArg{25}='badchan';
        inArg{26}=EPmain.preferences.preprocess.badchan;
        inArg{27}='blink';
        inArg{28}=EPmain.preferences.preprocess.blink;
        inArg{29}='badtrials';
        inArg{30}=EPmain.preferences.preprocess.badtrials;
        inArg{31}='chunkSize';
        inArg{32}=EPmain.preferences.preprocess.chunkSize;
        inArg{33}='minTrialsPerCell';
        inArg{34}=EPmain.preferences.preprocess.minTrialsPerCell;
        inArg{35}='noadjacent';
        inArg{36}=EPmain.preferences.preprocess.noadjacent;
        inArg{37}='trialMode';
        inArg{38}=trialMode;
        inArg{39}='trialminmax';
        inArg{40}=EPmain.preferences.preprocess.trialminmax;
        inArg{41}='movefacs';
        inArg{42}=EPmain.preferences.preprocess.movefacs;
        inArg{43}='textPrefs';
        inArg{44}=textPrefs;
        inArg{45}='noFigure';
        inArg{46}=EPmain.preferences.preprocess.noFigure;
        inArg{47}='sacctemplate';
        inArg{48}=saccadeTemplate;
        inArg{49}='sacpot';
        inArg{50}=EPmain.preferences.preprocess.sacPot;
        inArg{51}='editMode';
        inArg{52}=editMode;
        inArg{53}='detrend';
        inArg{54}=detrend;
        inArg{55}='fMRI';
        inArg{56}=fMRI;
        inArg{57}='elecPrefs';
        inArg{58}=EPmain.preferences.general.rotateHead;
        inArg{59}='SMIsuffix';
        inArg{60}=EPmain.preferences.general.SMIsuffix;
        inArg{61}='screenSize';
        inArg{62}=EPmain.scrsz;
        inArg{63}='FontSize';
        inArg{64}=EPmain.fontsize;
        inArg{65}='EMG';
        inArg{66}=EPmain.preprocess.EMG;
        inArg{67}=EPmain.preferences.preprocess.EMGratio;
        inArg{68}=EPmain.preferences.preprocess.EMGthresh;
        inArg{69}='alpha';
        inArg{70}=EPmain.preprocess.alpha;
        if any(strcmp(blinkTemplate,{'bothTemplate','fileTemplate'}))
            inArg{end+1}='blinkFile';
            inArg{end+1}=blinkFile;
        end;
        if any(strcmp(saccadeTemplate,{'bothTemplate','fileTemplate'}))
            inArg{end+1}='saccadeFile';
            inArg{end+1}=saccadeFile;
        end;
        if ~isempty(EPmain.preferences.preprocess.EOGchans)
            inArg{end+1}='eog';
            inArg{end+1}=EPmain.preferences.preprocess.EOGchans;
        end;
        SMIsuffix=EPmain.preferences.general.SMIsuffix;
        if ~isempty(SMIsuffix)
            inArg{end+1}='SMIsuffix';
            inArg{end+1}=SMIsuffix;
        end;
        specSuffix=EPmain.preferences.general.specSuffix;
        if ~isempty(specSuffix)
            inArg{end+1}='specSuffix';
            inArg{end+1}=specSuffix;
        end;
        subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
        if ~isempty(subjectSpecSuffix)
            inArg{end+1}='subjectSpecSuffix';
            inArg{end+1}=subjectSpecSuffix;
        end;
        
        errorFlag=0;
        msg=cell(0);
        switch EPmain.preprocess.origRefType
            case 1
                msg{end+1}='Original recording reference could not have been average reference.';
                errorFlag=1;
            case 2
                refChan=str2num(EPmain.preprocess.origReference);
                if isempty(EPmain.preprocess.origReference)
                    msg{end+1}='Please specify explicit recording reference channel(s).';
                    errorFlag=1;
                elseif isempty(refChan)
                    msg{end+1}='Recording reference channel(s) need to be numbers.';
                    errorFlag=1;
                elseif length(refChan) > 2
                    msg{end+1}='No more than two recording reference channels can be indicated.';
                    errorFlag=1;
                end;
                inArg{end+1}='origReference';
                inArg{end+1}=refChan;
            case 3
                msg{end+1}='Cannot preprocess CSD data.';
                [msg]=ep_errorMsg(msg);
                return
            case 4
                if ~isempty(EPmain.preprocess.origReference)
                    msg{end+1}='Recording reference channel field should be empty.';
                    errorFlag=1;
                end;
        end;
        
        switch EPmain.preprocess.currRefType
            case 1
                if ~isempty(EPmain.preprocess.currReference)
                    msg{end+1}='Current reference channel field should be empty.';
                    errorFlag=1;
                else
                    inArg{end+1}='currReference';
                    inArg{end+1}='AVG';
                end;
%             case 2
%                 if ~isempty(EPmain.preprocess.currReference)
%                     msg{1}='Current reference channel field should be empty.';
%                     [msg]=ep_errorMsg(msg);
%                     return
%                 else
%                     inArg{end+1}='currReference';
%                     inArg{end+1}='none';
%                 end;
            case 2
                refChan=str2num(EPmain.preprocess.currReference);
                if isempty(EPmain.preprocess.currReference)
                    msg{end+1}='Please specify explicit current reference channel(s).';
                    errorFlag=1;
                elseif isempty(refChan)
                    msg{end+1}='Current reference channel(s) need to be numbers.';
                    errorFlag=1;
                elseif length(refChan) > 2
                    msg{end+1}='No more than two current reference channels can be indicated.';
                    errorFlag=1;
                end;
                inArg{end+1}='currReference';
                inArg{end+1}=refChan;
            case 3
                msg{end+1}='Cannot preprocess CSD data.';
                [msg]=ep_errorMsg(msg);
                return
            case 4
                if ~isempty(EPmain.preprocess.currReference)
                    msg{1}='Current reference channel field should be empty.';
                    errorFlag=1;
                end;
        end;
        
        if ~isempty(timePoints)
            if timePoints
                inArg{end+1}='timePoints';
                inArg{end+1}=timePoints;
            end;
        end;
        
        if ~isempty(baseline)
            if baseline
                inArg{end+1}='baseline';
                inArg{end+1}=baseline;
            end;
        end;
        
        if get(EPmain.handles.preprocess.check,'value') %single cell mode
            subPos=str2num(get(EPmain.handles.preprocess.subject,'string'));
            cellPos=str2num(get(EPmain.handles.preprocess.cell,'string'));
            
            if min(subPos) < 1
                beep
                set(EPmain.handles.preprocess.subject,'ForegroundColor','red');
                drawnow
                return;
            end;
            if min(cellPos) < 1
                beep
                set(EPmain.handles.preprocess.cell,'ForegroundColor','red');
                drawnow
                return;
            end;
            
            for theFile=1:size(sessionFiles,2)
                [pathstr, name, fileSuffix] = fileparts(sessionFiles{theFile});
                if ~isempty(subPos)
                    if max(subPos) > length(name)
                        msg{end+1}=['The file name ' [name fileSuffix] ' does not have ' num2str(max(subPos)) ' characters for the subject label.'];
                        errorFlag=1;
                    else
                        theSubs{theFile}=name(subPos);
                    end;
                else
                    theSubs{theFile}='S01';
                end;
                if ~isempty(cellPos)
                    if max(cellPos) > length(name)
                        msg{end+1}=['The file name ' [name fileSuffix] ' does not have ' num2str(max(cellPos)) ' characters for the cell label.'];
                        errorFlag=1;
                    else
                        theCells{theFile}=name(cellPos);
                    end;
                else
                    theCells{theFile}='C01';
                end;
            end;
            
            if errorFlag
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.preprocess.preprocess,'enable','on');
                set(EPmain.handles.preprocess.done,'enable','on');
                set(EPmain.handles.preprocess.template,'enable','on');
                return
            end;
            
            uniqueSubs=unique(theSubs);
            
            mergeName = char(inputdlg('Name of new merged dataset?','Dataset name'));
            pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
            
            inArg{2}=[];
            inArg{4}='ep_mat';
            eloc=[];
            ced=[];
            
            for sub=1:length(uniqueSubs)
                mergeArg=[];
                theFiles=find(strcmp(uniqueSubs(sub),theSubs));
                for file=1:length(theFiles)
                    mergeArg{file,1}=sessionFiles{theFiles(file)};
                    mergeArg{file,2}='format';
                    mergeArg{file,3}=importFormat;
                    mergeArg{file,4}='type';
                    mergeArg{file,5}=dataType;
                    mergeArg{file,6}='labels';
                    mergeArg{file,7}={theCells(theFiles(file)) theSubs(theFiles(file)) cell(0)};
                    mergeArg{file,8}='textPrefs';
                    mergeArg{file,9}=textPrefs;
                    mergeArg{file,10}='elecPrefs';
                    mergeArg{file,11}=EPmain.preferences.general.rotateHead;
                    mergeArg{file,12}='screenSize';
                    mergeArg{file,13}=EPmain.scrsz;
                    mergeArg{file,14}='FontSize';
                    mergeArg{file,15}=EPmain.fontsize;
                    if ~isempty(eloc)
                        mergeArg{file,16}='eloc';
                        mergeArg{file,17}=eloc;
                    end;
                    if ~isempty(ced)
                        mergeArg{file,18}='ced';
                        mergeArg{file,19}=ced;
                    end;
                    if ~isempty(EPmain.preferences.general.SMIsuffix)
                        mergeArg{file,20}='SMIsuffix';
                        mergeArg{file,21}=EPmain.preferences.general.SMIsuffix;
                    end;
                    if ~isempty(EPmain.preferences.general.specSuffix)
                        mergeArg{file,22}='specSuffix';
                        mergeArg{file,23}=EPmain.preferences.general.specSuffix;
                    end;
                    if ~isempty(EPmain.preferences.general.subjectSpecSuffix)
                        mergeArg{file,22}='subjectSpecSuffix';
                        mergeArg{file,23}=EPmain.preferences.general.subjectSpecSuffix;
                    end;
                end;
                [EPdata eloc]=ep_mergeEPfiles(mergeArg,mergeName);
                if isempty(EPdata)
                    beep();
                    set(EPmain.handles.preprocess.preprocess,'enable','on');
                    set(EPmain.handles.preprocess.done,'enable','on');
                    set(EPmain.handles.preprocess.template,'enable','on');
                    return
                end;
                ced=EPdata.ced;
                [err]=ep_writeData(EPdata,[mergeName num2str(sub)],EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,'ep_mat');
                if ~isempty(err)
                    inArg{2}{end+1}=[mergeName num2str(sub) '.ept'];
                    disp(['Generating temporary work file: ' inArg{2}{end} '.']);
                else
                    msg{1}='Writing temporary work file during file merging failed.  Aborting attempt.';
                    [msg]=ep_errorMsg(msg);
                    set(EPmain.handles.preprocess.preprocess,'enable','on');
                    set(EPmain.handles.preprocess.done,'enable','on');
                    set(EPmain.handles.preprocess.template,'enable','on');
                    return
                end;
            end;
        end;
        
        ep_artifactCorrection(inArg);
        
        if get(EPmain.handles.preprocess.check,'value') %single cell mode
            disp('Deleting temporary work files.');
            for i=1:length(inArg{2})
                delete(inArg{2}{i}); %delete temporary merged files
            end;
        end;
        
        try
            set(EPmain.handles.preprocess.preprocess,'enable','on');
            set(EPmain.handles.preprocess.done,'enable','on');
            set(EPmain.handles.preprocess.template,'enable','on');
        catch
            ep('start')
        end;
        drawnow
        
        ep('start');
        
    case 'startTransform'
        
        set(EPmain.handles.hMainWindow,'Name', 'Transform Data');
        
        uicontrol('Style','text',...
            'String','Input File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 470 150 20]);
        
        EPmain.handles.transform.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.importFormat,''Value'');','if tempVar ~=0,EPmain.transform.importFormat=tempVar;end;','if isempty(tempVar),EPmain.transform.importFormat=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.importFormat,'Position',[20 450 150 20]);
        
        uicontrol('Style','text',...
            'String','Input File Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 430 150 20]);
        
        EPmain.handles.transform.type= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.type,''Value'');','if tempVar ~=0,EPmain.transform.type=tempVar;end;','if isempty(tempVar),EPmain.transform.type=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.type,'Position',[20 410 160 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.transform.importFormat});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.transform.type,'enable','off');
        end;
        
        uicontrol('Style','text',...
            'String','Ouput File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 390 150 20]);
        
        EPmain.handles.transform.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.outputFormat,''Value'');','if tempVar ~=0,EPmain.transform.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.transform.outputFormat=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.outputFormat,'Position',[20 370 150 20]);
        
        uicontrol('Style','text',...
            'String','Data mode:','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 350 105 20]);
        
        EPmain.handles.transform.dataMode= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'EEG','MEG','pupil','ECG', 'ANS'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.dataMode,''Value'');','if tempVar ~=0,EPmain.transform.dataMode=tempVar;end;','if isempty(tempVar),EPmain.transform.dataMode=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.dataMode,'Position',[105 350 100 20]);
        
        EPmain.handles.transform.referenceLabel= uicontrol('Style','text',...
            'String','Rereference EEG','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 330 105 20]);
        
        EPmain.handles.transform.reference= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Average','Traditional','CSD','PARE','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.reference,''Value'');','if tempVar ~=0,EPmain.transform.reference=tempVar;end;','if isempty(tempVar),EPmain.transform.reference=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.reference,'Position',[105 330 100 20]);
        
        if EPmain.transform.dataMode ~= 1 %not EEG
            set(EPmain.handles.transform.referenceLabel,'enable','off');
            set(EPmain.handles.transform.reference,'enable','off');
        end;

        EPmain.handles.transform.referenceLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference Channel(s)','FontSize',EPmain.fontsize,...
            'Position',[25 310 150 20]);
        
        EPmain.handles.transform.refChan1 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.refChan1,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.refChan1,''String''),','if tempVar ~=0,EPmain.transform.refChan1=tempVar;end;','if isempty(tempVar),EPmain.transform.refChan1=tempVar;end;','ep(''start'');'],...
            'Position',[25 290 50 20]);
        
        EPmain.handles.transform.refChan2 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.refChan2,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.refChan2,''String'');','if tempVar ~=0,EPmain.transform.refChan2=tempVar;end;','if isempty(tempVar),EPmain.transform.refChan2=tempVar;end;','ep(''start'');'],...
            'Position',[85 290 50 20]);
        
        if (EPmain.transform.reference ~= 2) || (EPmain.transform.dataMode ~= 1) 
            set(EPmain.handles.transform.referenceLabel,'enable','off');
            set(EPmain.handles.transform.refChan1,'enable','off');
            set(EPmain.handles.transform.refChan2,'enable','off');
        end;
        
        EPmain.handles.transform.detrend= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Detrend',...
            'CallBack',['global EPmain;','EPmain.transform.detrend=get(EPmain.handles.transform.detrend,''Value'');','ep(''start'');'],...
            'Value',EPmain.transform.detrend,'Position',[20 270 160 20],'TooltipString','Detrend data (recommended only for continuous data as it will also detrend ERP components).');

        uicontrol('Style','text','HorizontalAlignment','left','String', 'Prestim','FontSize',EPmain.fontsize,...
            'Position',[5 250 40 20]);
        
        EPmain.handles.transform.baselineLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline','FontSize',EPmain.fontsize,...
            'Position',[55 250 50 20]);
        
        EPmain.handles.transform.delayLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Delay','FontSize',EPmain.fontsize,...
            'Position',[145 250 40 20]);
        
        EPmain.handles.transform.preStim = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.preStim,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.preStim,''String'');','if tempVar ~=0,EPmain.transform.preStim=tempVar;end;','if isempty(tempVar),EPmain.transform.preStim=tempVar;end;','ep(''start'');'],...
            'Position',[5 230 40 20],'TooltipString','Msec start of epoch with respect to stimulus with positive being prior to stimulus.');
        
        EPmain.handles.transform.baselineStart = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.baselineStart,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.baselineStart,''String'');','if tempVar ~=0,EPmain.transform.baselineStart=tempVar;end;','if isempty(tempVar),EPmain.transform.baselineStart=tempVar;end;','ep(''start'');'],...
            'Position',[55 230 40 20],'TooltipString','Msec start (left side of sample) of baseline period.');
        
        EPmain.handles.transform.baselineEnd = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.baselineEnd,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.baselineEnd,''String'');','if tempVar ~=0,EPmain.transform.baselineEnd=tempVar;end;','if isempty(tempVar),EPmain.transform.baselineEnd=tempVar;end;','ep(''start'');'],...
            'Position',[95 230 40 20],'TooltipString','Msec end (right side of sample) of baseline period.');
        
        EPmain.handles.transform.delay = uicontrol('Style','text','HorizontalAlignment','left','String', EPmain.transform.delay,'FontSize',EPmain.fontsize,...
            'Position',[145 230 20 20],'TooltipString','Msec latency delay from use of filter.  Correct during segmentation.');
        
        EPmain.handles.transform.filterPass= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Low Pass','High Pass','Band Pass','Band Stop','Notch'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.filterPass,''Value'');','if tempVar ~=0,EPmain.transform.filterPass=tempVar;end;','if isempty(tempVar),EPmain.transform.filterPass=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.filterPass,'Position',[5 210 110 20]);
        
        EPmain.handles.transform.filter1 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.filter1),'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.filter1,''String''));','if tempVar ~=0,EPmain.transform.filter1=tempVar;end;','if isempty(tempVar),EPmain.transform.filter1=tempVar;end;','ep(''start'');'],...
            'Position',[110 210 40 20],'TooltipString','Lower frequency limit.');
        
        EPmain.handles.transform.filter2 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.filter2),'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.filter2,''String''));','if tempVar ~=0,EPmain.transform.filter2=tempVar;end;','if isempty(tempVar),EPmain.transform.filter2=tempVar;end;','ep(''start'');'],...
            'Position',[155 210 40 20],'TooltipString','Upper frequency limit.');

        EPmain.handles.transform.filterType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'One-Pass Butterworth','Two-Pass Butterworth','One-Pass FIR','Two-Pass FIR','One-Pass FIRLS','Two-Pass FIRLS'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.filterType,''Value'');','if tempVar ~=0,EPmain.transform.filterType=tempVar;end;','if isempty(tempVar),EPmain.transform.filterType=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.filterType,'Position',[5 190 160 20]);
        
        EPmain.handles.transform.filterOrder = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.filterOrder),'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.filterOrder,''String''));','if tempVar ~=0,EPmain.transform.filterOrder=tempVar;end;','if isempty(tempVar),EPmain.transform.filterOrder=tempVar;end;','ep(''start'');'],...
            'Position',[160 190 30 20],'TooltipString','Order of the filter.');
        
        freqFilter=[];
        freqFilter.trial{1}=EPmain.transform.whiteNoise;
        freqFilter.label{1}='e1';
        freqFilter.time{1}=[-50:199]*.004;
        freqFilter.sampleinfo=[1 250];
        freqFilter.fsample=250;
        
        timeFilter=[];
        timeFilter.trial{1}=[repmat(0,1,100) 1 repmat(0,1,149)];
        timeFilter.label{1}='e1';
        timeFilter.time{1}=[-50:199]*.004;
        timeFilter.sampleinfo=[1 250];
        timeFilter.fsample=250;
        
        cfg=[];
        cfg.pad = 'nextpow2';
        switch EPmain.transform.filterType
            case 1 
                filterDirection='onepass';
                filterType='but';
            case 2 
                filterDirection='twopass';
                filterType='but';
            case 3 
                filterDirection='onepass';
                filterType='fir';
            case 4 
                filterDirection='twopass';
                filterType='fir';
            case 5 
                filterDirection='onepass';
                filterType='firls';
            case 6
                filterDirection='twopass';
                filterType='firls';
        end;
        switch EPmain.transform.filterPass
            case 1 %low pass
                cfg.lpfilter='yes';
                cfg.lpfreq=EPmain.transform.filter1;
                cfg.lpfiltdir=filterDirection;
                cfg.lpfilttype=filterType;
                cfg.lpfiltord=EPmain.transform.filterOrder;
                set(EPmain.handles.transform.filter2,'enable','off');
            case 2 %high pass
                cfg.hpfilter='yes';
                cfg.hpfreq=EPmain.transform.filter1;
                cfg.hpfiltdir=filterDirection;
                cfg.hpfilttype=filterType;
                cfg.hpfiltord=EPmain.transform.filterOrder;
                set(EPmain.handles.transform.filter2,'enable','off');
            case 3 %band pass
                cfg.bpfilter='yes';
                cfg.bpfreq=[EPmain.transform.filter1 EPmain.transform.filter2];
                cfg.bpfiltdir=filterDirection;
                cfg.bpfilttype=filterType;
                cfg.bpfiltord=EPmain.transform.filterOrder;
            case 4 %band stop
                cfg.bsfilter='yes';
                cfg.bsfreq=[EPmain.transform.filter1 EPmain.transform.filter2];
                cfg.bsfiltdir=filterDirection;
                cfg.bsfilttype=filterType;
                cfg.bsfiltord=EPmain.transform.filterOrder;
            case 5 %notch
                cfg.dftfilter='yes';
                cfg.dftfreq=EPmain.transform.filter1;
                cfg.dftfiltdir=filterDirection;
                cfg.dftfilttype=filterType;
                cfg.dftfiltord=EPmain.transform.filterOrder;
                set(EPmain.handles.transform.filter2,'enable','off');
        end;
        
        EPmain.transform.cfgFilter=cfg;
        cfg2=[];
        cfg2.taper='hanning';
        cfg2.method='mtmfft';
        cfg2.output='fourier';
        cfg2.pad = 'nextpow2';
        filterProb=0;
        
        if ~isempty(EPmain.transform.filter1) && ~(isempty(EPmain.transform.filter2) && ismember(EPmain.transform.filterPass,[3 4]))
            try
                evalc('[timeFiltered] = ft_preprocessing(cfg, timeFilter);');
                evalc('[freqFiltered] = ft_preprocessing(cfg, freqFilter);');
                evalc('[freqFilterFFT] = ft_freqanalysis(cfg2, freqFilter);');
                evalc('[freqFilteredFFT] = ft_freqanalysis(cfg2, freqFiltered);');
            catch
                disp('Filter settings did not work.  Try different settings.')
                filterProb=1;
                freqFilterFFT.freq=[0:125];
                freqFilterFFT.fourierspctrm=zeros(126,1);
                freqFilteredFFT.fourierspctrm=zeros(126,1);
                timeFiltered.trial{1}=timeFilter.trial{1};
            end
            [A originalLatency]=max(timeFilter.trial{1});
            [A newLatency]=max(timeFiltered.trial{1});
            EPmain.transform.delay=(newLatency-originalLatency)*(1000/timeFilter.fsample);
            set(EPmain.handles.transform.delay,'string',num2str(EPmain.transform.delay));
        else
            freqFilterFFT.freq=[0:125];
            freqFilterFFT.fourierspctrm=zeros(126,1);
            freqFilteredFFT.fourierspctrm=zeros(126,1);
            timeFiltered.trial{1}=timeFilter.trial{1};
        end;
        
        %Frequency domain plot of filter effects
        EPmain.transform.freqData=[abs(squeeze(freqFilteredFFT.fourierspctrm))';abs(squeeze(freqFilterFFT.fourierspctrm))'];
        EPmain.transform.freqScale=freqFilterFFT.freq;
        EPmain.transform.freqAxis=[freqFilterFFT.freq(1) freqFilterFFT.freq(end) -.1 max([.1; ([abs(squeeze(freqFilterFFT.fourierspctrm))])])];
        EPmain.handles.transform.frequencyFilter = axes('units','pixels','position',[15 130 75 50]);
        EPmain.handles.transform.freqWaves = plot(EPmain.transform.freqScale,EPmain.transform.freqData);
        axis(EPmain.transform.freqAxis);
        set(EPmain.handles.transform.frequencyFilter,'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.freqWaves(1),'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.freqWaves(2),'ButtonDownFcn',@expandTransformChan);
            
        %Time domain plot of filter effects
        EPmain.transform.timeData=[timeFiltered.trial{1};timeFilter.trial{1}];
        EPmain.transform.timeScale=timeFilter.time{1}*1000;
        EPmain.transform.timeAxis=[EPmain.transform.timeScale(1) EPmain.transform.timeScale(end) min([-1.2, timeFilter.time{1}]) max([1.2, timeFilter.time{1}])];
        EPmain.handles.transform.timeFilter = axes('units','pixels','position',[115 130 75 50]);
        EPmain.handles.transform.timeWaves = plot(EPmain.transform.timeScale,EPmain.transform.timeData);
        axis(EPmain.transform.timeAxis);
        set(EPmain.handles.transform.timeFilter,'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.timeWaves(1),'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.timeWaves(2),'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.hMainWindow,'DefaultAxesColorOrder',[[1 0 0;0 0 1]]);
        
        EPmain.handles.transform.domainLabel=uicontrol('Style','text',...
            'String','FFT','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 90 40 20]);
        
        EPmain.handles.transform.domain= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'no change','Frequency','Time-Frequency'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.domain,''Value'');','if tempVar ~=0,EPmain.transform.domain=tempVar;EPmain.transform.method=1;end;','if isempty(tempVar),EPmain.transform.domain=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.domain,'Position',[45 90 150 20]);
        
        EPmain.handles.transform.methodLabel=uicontrol('Style','text',...
            'String','Method','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 70 40 20]);
        
        switch EPmain.transform.domain
            case 1 %Time
                menuItems={'none'};
            case 2 %Frequency
                menuItems={'multi-taper','Hanning'};
            case 3 %Time-Frequency
                menuItems={'multi-taper','Hanning','wavelet multiplication','wavelet convolution'};
        end;
        
        EPmain.handles.transform.method= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',menuItems,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.method,''Value'');','if tempVar ~=0,EPmain.transform.method=tempVar;end;','if isempty(tempVar),EPmain.transform.method=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.transform.method,'Position',[45 70 150 20]);
        
        if EPmain.transform.dataMode >2
            set(EPmain.handles.transform.domainLabel,'enable','off');
            set(EPmain.handles.transform.methodLabel,'enable','off');
            set(EPmain.handles.transform.method,'enable','off');
            set(EPmain.handles.transform.domain,'enable','off');
        end;

        EPmain.handles.transform.smoothingLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Smoothing','FontSize',EPmain.fontsize,...
            'Position',[25 50 55 20]);
        EPmain.handles.transform.smoothing = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.smoothing,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.smoothing,''String'');','if tempVar ~=0,EPmain.transform.smoothing=tempVar;end;','if isempty(tempVar),EPmain.transform.smoothing=tempVar;end;','ep(''start'');'],...
            'Position',[25 30 50 20]);
        
        if ~ismember(EPmain.transform.domain,[2 3]) || (EPmain.transform.method ~= 1) || EPmain.transform.dataMode >2
            set(EPmain.handles.transform.smoothingLabel,'enable','off');
            set(EPmain.handles.transform.smoothing,'enable','off'); %smoothing only seems to be for multi-taper
        end;
        
        if (EPmain.transform.domain ~= 1) %Frequency domain analysis
            set(EPmain.handles.transform.baselineLabel,'enable','off');
            set(EPmain.handles.transform.baselineStart,'enable','off');
            set(EPmain.handles.transform.baselineEnd,'enable','off');
        end;
        
        EPmain.handles.transform.transform = uicontrol('Style', 'pushbutton', 'String', 'Transform','FontSize',EPmain.fontsize,...
            'Position', [20 0 80 30], 'Callback', 'ep(''transformData'')');
        
        if filterProb
            set(EPmain.handles.transform.transform,'enable','off');
        end;
        
        EPmain.handles.transform.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [100 0 80 30], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        

        
    case 'transformData'
        
        switch EPmain.transform.type
            case 1
                dataType='continuous';
            case 2
                dataType='single_trial';
            case 3
                dataType='average';
            case 4
                dataType='grand_average';
            case 5
                dataType='factors';
        end;
        
        importFormatNum = get(EPmain.handles.transform.importFormat,'value');
        [importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{importFormatNum});

        
        outputFormatNum = get(EPmain.handles.transform.outputFormat,'value');
        [outputSuffix,outputFormatName,outputFormat]=ep_fileFormats(dataType,EPmain.fileFormatSaveList{outputFormatNum});
        
        if EPmain.transform.dataMode ==1
            switch EPmain.transform.reference
                case 1
                    referenceMethod='Average';
                case 2
                    referenceMethod='Traditional';
                case 3
                    referenceMethod='CSD';
                case 4
                    referenceMethod='PARE';
                case 5
                    referenceMethod='none';
            end;
        else
            referenceMethod='none';
        end;
        
        
        switch EPmain.transform.dataMode
            case 1
                dataMode='EEG';
            case 2
                dataMode='MEG';
            case 3
                dataMode='pupil';
            case 4
                dataMode='ECG';
            case 5
                dataMode='ANS';
        end;
        
        EPmain.transform.detrend=get(EPmain.handles.transform.detrend,'value');
        
        if ~isempty(EPmain.transform.filter1) && ~(isempty(EPmain.transform.filter2) && ismember(EPmain.transform.filterPass,[3 4]))
            if (EPmain.transform.filterPass==3) && (EPmain.transform.filter1 >= EPmain.transform.filter2)
                msg{1}='The start of the bandpass window must come before the end of the bandpass window.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (EPmain.transform.filterPass==4) && (EPmain.transform.filter1 >= EPmain.transform.filter2)
                msg{1}='The start of the bandstop window must come before the end of the bandstop window.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (length(EPmain.transform.filter1) > 1) || (length(EPmain.transform.filter2) > 1) || (length(EPmain.transform.filterOrder) > 1)
                msg{1}='The filter setting in each field needs to be a single number.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (EPmain.transform.filter1 < 0) || (EPmain.transform.filterOrder < 0)
                msg{1}='Filter settings cannot be negative numbers.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if ~isempty(EPmain.transform.filter2) && (EPmain.transform.filter2 < 0)
                msg{1}='Filter settings cannot be negative numbers.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
        else
            EPmain.transform.cfgFilter=[];
        end;
        
        if EPmain.transform.dataMode <3
            switch EPmain.transform.domain
                case 1
                    domainName='Time';
                    methodName='None';
                case 2
                    domainName='Frequency';
                    switch EPmain.transform.method
                        case 1
                            methodName='multi-taper';
                        case 2
                            methodName='Hanning';
                    end;
                case 3
                    domainName='Time-Frequency';
                    switch EPmain.transform.method
                        case 1
                            methodName='multi-taper';
                        case 2
                            methodName='Hanning';
                        case 3
                            methodName='wavelet multiplication';
                        case 4
                            methodName='wavelet convolution';
                    end;
            end;
        else
            domainName='Time';
            methodName='None';
        end;
        
        EPmain.transform.importFormat=importFormatNum;
        EPmain.transform.outputFormat=outputFormatNum;
        EPmain.transform.refChan1=str2num(get(EPmain.handles.transform.refChan1,'String'));
        EPmain.transform.refChan2=str2num(get(EPmain.handles.transform.refChan2,'String'));
        EPmain.transform.baselineStart=str2num(get(EPmain.handles.transform.baselineStart,'String'));
        EPmain.transform.baselineEnd=str2num(get(EPmain.handles.transform.baselineEnd,'String'));
        EPmain.transform.preStim=str2num(get(EPmain.handles.transform.preStim,'String'));
        
        if EPmain.transform.baselineStart >= EPmain.transform.baselineEnd
            msg{1}='The start of the baseline period must come before the end of the baseline period.';
            if (EPmain.transform.baselineStart == 0) && (EPmain.transform.baselineEnd ==0)
                msg{2}='To disable baseline correction, set both fields to blank.';
            end;
            [msg]=ep_errorMsg(msg);
            return
        end
        
        [inputFiles, activeDirectory]=ep_getFilesUI(importFormat);
        if inputFiles{1}==0
            return %user hit cancel on file requestor
        end
        for theFile=1:size(inputFiles,2)
            inputFiles{theFile}=[activeDirectory inputFiles{theFile}];
        end;
        
        set(EPmain.handles.transform.transform,'enable','off');
        set(EPmain.handles.transform.done,'enable','off');
        drawnow
        ep_transformData(inputFiles,importFormat,dataType,outputFormat,referenceMethod,EPmain.transform,domainName,methodName,dataMode);
        try
            set(EPmain.handles.transform.transform,'enable','on');
            set(EPmain.handles.transform.done,'enable','on');
        catch
            ep('start')
        end;
        
        drawnow
        
        ep('start');
        
    case 'startAverage'
        
        set(EPmain.handles.hMainWindow,'Name', 'Average Data');
        
        if ~isfield(EPmain.average,'trialSpecDataset') %just entering average function
            EPmain.average.datasetList=[];
            EPmain.average.trialSpecDataset=1;
            EPmain.average.trialSpec=1;
            EPmain.average.channelDataset=1;
            EPmain.average.channel=1;
            EPmain.average.peakPolarity=1;
            EPmain.average.dropEvents=1;
            EPmain.average.subject='';
            for iFile=1:length(EPdataset.dataset)
                if ~isempty(EPdataset.dataset(iFile).trialSpecs) && strcmp(EPdataset.dataset(iFile).dataType,'single_trial')
                    EPmain.average.datasetList(end+1)=iFile;
                end;
            end;
        end;
        
        uicontrol('Style','text',...
            'String','Input File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 450 150 20]);
        
        EPmain.handles.average.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.importFormat,''Value'');','if tempVar ~=0,EPmain.average.importFormat=tempVar;end;','if isempty(tempVar),EPmain.average.importFormat=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.average.importFormat,'Position',[20 430 150 20]);
        
        uicontrol('Style','text',...
            'String','Input File Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 410 150 20]);
        
        EPmain.handles.average.type= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.type,''Value'');','if tempVar ~=0,EPmain.average.type=tempVar;end;','if isempty(tempVar),EPmain.average.type=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.average.type,'Position',[20 390 160 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.average.importFormat});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.average.type,'enable','off');
        end;
        
        uicontrol('Style','text',...
            'String','Ouput File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 370 150 20]);
        
        EPmain.handles.average.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.outputFormat,''Value'');','if tempVar ~=0,EPmain.average.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.average.outputFormat=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.average.outputFormat,'Position',[20 350 150 20]);

        EPmain.handles.average.dropEvents= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Drop events',...
            'CallBack',['global EPmain;','EPmain.average.dropEvents=get(EPmain.handles.average.dropEvents,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Value',EPmain.average.dropEvents,'Position',[30 330 150 20]);        
        
        EPmain.handles.average.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.average.subject,...
            'CallBack',['global EPmain;','EPmain.average.subject=get(EPmain.handles.average.subject,''String'');','ep(''start'');'],...
            'Position',[30 310 50 20],'TooltipString','For when a subject is spread across multiple files and its ID is encoded in the file name, example 4:6 for 4th-6th characters of ''sub001.ept''.');
        
        EPmain.handles.average.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Subject','HorizontalAlignment','left',...
            'Position',[90 310 50 20]);
        
        uicontrol('Style','text',...
            'String','Procedure','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 270 150 20]);
        
        EPmain.handles.average.method= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Average','Latency-Lock','Jitter-Correct','Frequency-Coherence','Time-Frequency-Coherence'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.method,''Value'');','if tempVar ~=0,EPmain.average.method=tempVar;end;','if isempty(tempVar),EPmain.average.method=tempVar;end;','EPmain.average.freqMethod=1;','ep(''start'');'],...
            'Value',EPmain.average.method,'Position',[20 240 150 20]);
        
        uicontrol('Style','text',...
            'String','Method','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 210 40 20]);
        
        switch EPmain.average.method
            case 1 %Average
                menuItems={'Mean','Median','Trimmed Mean'};
            case 2 %Latency Lock
                menuItems={'Mean','Median','Trimmed Mean'};
            case 3 %Jitter-Correct
                menuItems={'Mean','Median','Trimmed Mean'};
            case 4 %Frequency-Coherence
                menuItems={'multi-taper','Hanning'};
            case 5 %Time-Frequency-Coherence
                menuItems={'phase-lock wavelet'};
        end;
        
        EPmain.handles.average.freqMethod= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',menuItems,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.freqMethod,''Value'');','if tempVar ~=0,EPmain.average.freqMethod=tempVar;end;','if isempty(tempVar),EPmain.average.freqMethod=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.average.freqMethod,'Position',[45 210 150 20]);
        
        switch EPmain.average.method
            case 1
                
            case 2
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Model','FontSize',EPmain.fontsize,...
                    'Position',[5 180 40 20]);
                
                if ~isempty(EPmain.average.datasetList)
                    EPmain.handles.average.trialSpecDataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',{EPdataset.dataset(EPmain.average.datasetList).dataName},...
                        'Value',EPmain.average.trialSpecDataset,'Position',[60 180 120 20],'TooltipString','Model single-trial dataset for providing trial specs for jitter correction.',...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.trialSpecDataset,''Value'');','if tempVar ~=0,EPmain.average.trialSpecDataset=tempVar;end;','if isempty(tempVar),EPmain.average.trialSpecDataset=tempVar;end;','EPmain.average.trialSpec=1;','ep(''start'');']);
                else
                    EPmain.handles.average.trialSpecDataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'Value',1,'Position',[60 180 120 20],'TooltipString','Model single-trial dataset for providing trial specs for jitter correction.');
                end;
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Trial Specs','FontSize',EPmain.fontsize,...
                    'Position',[25 160 88 20]);
                
                if ~isempty(EPmain.average.datasetList)
                    EPmain.handles.average.trialSpec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames,...
                        'Value',EPmain.average.trialSpec,'Position',[20 140 160 20],'TooltipString','Trial spec to be used for jitter correction.',...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.trialSpec,''Value'');','if tempVar ~=0,EPmain.average.trialSpec=tempVar;end;','if isempty(tempVar),EPmain.average.trialSpec=tempVar;end;','ep(''start'');']);
                else
                    EPmain.handles.average.trialSpec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'Value',1,'Position',[20 140 160 20],'TooltipString','Trial spec to be used for jitter correction.');
                end;
                
                uicontrol('Style','text',...
                    'String','Trial Spec Range','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[25 120 88 20]);
                
                uicontrol('Style','text',...
                    'String','Prestim','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[125 120 88 20]);

                EPmain.handles.average.minLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.minLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.minLatency=str2num(get(EPmain.handles.average.minLatency,''String''));','ep(''start'');'],...
                    'Position',[25 100 50 20],'TooltipString','Lower limit of latencies to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.maxLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.maxLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.maxLatency=str2num(get(EPmain.handles.average.maxLatency,''String''));','ep(''start'');'],...
                    'Position',[75 100 50 20],'TooltipString','Upper limit of latencies to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                if ~isempty(EPmain.average.datasetList)
                    EPmain.handles.average.prestim = uicontrol('Style','text','HorizontalAlignment','left','String', EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).baseline,'FontSize',EPmain.fontsize,...
                    'Position',[125 100 50 20]);
                else
                    EPmain.handles.average.prestim = uicontrol('Style','text','HorizontalAlignment','left','String', '','FontSize',EPmain.fontsize,...
                    'Position',[125 100 50 20]);
                end;
            case 3
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Model','FontSize',EPmain.fontsize,...
                    'Position',[5 180 50 20]);
                
                EPmain.handles.average.channelDataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{EPdataset.dataset.dataName},...
                    'Value',EPmain.average.channelDataset,'Position',[45 180 150 20],'TooltipString','Model dataset for providing channel for jitter correction.',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.channelDataset,''Value'');','if tempVar ~=0,EPmain.average.channelDataset=tempVar;end;','if isempty(tempVar),EPmain.average.channelDataset=tempVar;end;','EPmain.average.channel=1;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Channel','FontSize',EPmain.fontsize,...
                    'Position',[5 160 50 20]);
                
                EPmain.handles.average.channel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPdataset.dataset(EPmain.average.channelDataset).chanNames,...
                    'Value',EPmain.average.channel,'Position',[45 160 150 20],'TooltipString','Channel used for jitter correction.',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.channel,''Value'');','if tempVar ~=0,EPmain.average.channel=tempVar;end;','if isempty(tempVar),EPmain.average.channel=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Polarity','FontSize',EPmain.fontsize,...
                    'Position',[5 140 50 20]);
                
                EPmain.handles.average.peakPolarity = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'Positive','Negative'},...
                    'Value',EPmain.average.peakPolarity,'Position',[45 140 150 20],'TooltipString','Use positive or negative peak for jitter correction.',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.peakPolarity,''Value'');','if tempVar ~=0,EPmain.average.peakPolarity=tempVar;end;','if isempty(tempVar),EPmain.average.peakPolarity=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text',...
                    'String','Range','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 120 50 20]);
                
                EPmain.handles.average.minLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.minLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.minLatency=str2num(get(EPmain.handles.average.minLatency,''String''));','ep(''start'');'],...
                    'Position',[45 120 50 20],'TooltipString','Lower limit of peaks to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.maxLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.maxLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.maxLatency=str2num(get(EPmain.handles.average.maxLatency,''String''));','ep(''start'');'],...
                    'Position',[100 120 50 20],'TooltipString','Upper limit of peaks to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');

            case {4,5}
                
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Smoothing','FontSize',EPmain.fontsize,...
                    'Position',[25 180 55 20]);
                EPmain.handles.average.smoothing = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.smoothing,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.smoothing,''String'');','if tempVar ~=0,EPmain.average.smoothing=tempVar;end;','if isempty(tempVar),EPmain.average.smoothing=tempVar;end;','ep(''start'');'],...
                    'Position',[25 160 50 20]);
                
                if ~ismember(EPmain.average.method,[4 5]) || (EPmain.average.freqMethod ~= 1)
                    set(EPmain.handles.average.smoothing,'enable','off'); %smoothing only seems to be for multi-taper
                end;
                
        end;
        
        uicontrol('Style','text',...
            'String','Drop:','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 80 30 20]);
        
        EPmain.handles.average.dropBad= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Bad',...
            'CallBack',['global EPmain;','EPmain.average.dropBad=get(EPmain.handles.average.dropBad,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Value',EPmain.average.dropBad,'Position',[35 80 50 20],'TooltipString','Drop bad trials.');           

        EPmain.handles.average.dropError= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Error',...
            'CallBack',['global EPmain;','EPmain.average.dropError=get(EPmain.handles.average.dropError,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Value',EPmain.average.dropError,'Position',[85 80 50 20],'TooltipString','Drop trials where ACC is recorded as an error trial with a 0 value.');                      
        
        EPmain.handles.average.dropTimeout= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Timeout',...
            'CallBack',['global EPmain;','EPmain.average.dropTimeout=get(EPmain.handles.average.dropTimeout,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Value',EPmain.average.dropTimeout,'Position',[135 80 70 20],'TooltipString','Drop trials where ACC is recorded as a timeout trial with a 2 value.');           
        
        uicontrol('Style','text',...
            'String','RT Method','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 60 60 20]);
        
        EPmain.handles.average.RTmethod= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Median','Mean','Trimmed Mean'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.RTmethod,''Value'');','if tempVar ~=0,EPmain.average.RTmethod=tempVar;end;','if isempty(tempVar),EPmain.average.RTmethod=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.average.RTmethod,'Position',[65 60 130 20],'TooltipString','Method for summarizing the reaction times.');
        
        uicontrol('Style','text',...
            'String','Min RT','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 40 40 20]);
        
        EPmain.handles.average.minRT= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.average.minRT,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.minRT,''String'');','if tempVar ~=0,EPmain.average.minRT=tempVar;end;','if isempty(tempVar),EPmain.average.minRT=tempVar;end;','ep(''start'');'],...
            'Position',[50 40 40 20],'TooltipString','RT less than this number is discarded from reaction time summary number.  Set to zero to disable this setting.');
        
        uicontrol('Style','text',...
            'String','Max RT SD','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[95 40 60 20]);
        
        EPmain.handles.average.maxRT= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.average.maxRT,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.maxRT,''String'');','if tempVar ~=0,EPmain.average.maxRT=tempVar;end;','if isempty(tempVar),EPmain.average.maxRT=tempVar;end;','ep(''start'');'],...
            'Position',[150 40 40 20],'TooltipString','RT more than this number in standard deviation units is discarded from reaction time summary number.  Set to zero to disable this setting.');
        
        EPmain.handles.average.average = uicontrol('Style', 'pushbutton', 'String', 'Average','FontSize',EPmain.fontsize,...
            'Position', [70 0 60 35], 'Callback', 'ep(''averageData'')');
        
        if (EPmain.average.method==2) && isempty(EPmain.average.datasetList)
            set(EPmain.handles.average.average,'enable','off');
        end;
        
        EPmain.handles.average.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.average=rmfield(EPmain.average,''trialSpecDataset'');','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'averageData'
        
        typeNum = get(EPmain.handles.average.type,'value');
        switch typeNum
            case 1
                dataType='continuous';
                dataOutType='average';
            case 2
                dataType='single_trial';
                dataOutType='average';
            case 3
                dataType='average';
                dataOutType='grand_average';
            case 4
                dataType='grand_average';
                dataOutType='grand_average';
            case 5
                dataType='factors';
                dataOutType='grand_average';
        end;
        
        importFormatNum = get(EPmain.handles.average.importFormat,'value');        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{importFormatNum});

        outputFormatNum = get(EPmain.handles.average.outputFormat,'value');
        [outputSuffix,outputFormatName,outputFormat]=ep_fileFormats(dataOutType,EPmain.fileFormatSaveList{outputFormatNum});
        
        if strcmp(importFormatName,'ep_mat')
            dataType=[];
        end;
        
        jitterChan=[];
        procedureNum = get(EPmain.handles.average.method,'value');
        methodNum = get(EPmain.handles.average.freqMethod,'value');
        switch procedureNum
            case 1
                averagingMethod='Average';
                switch methodNum
                    case 1
                        methodName='Mean';
                    case 2
                        methodName='Median';
                    case 3
                        methodName='Trimmed_Mean';
                end;
                suffix='_avg';
            case 2
                averagingMethod='Latency-Lock';
                switch methodNum
                    case 1
                        methodName='Mean';
                    case 2
                        methodName='Median';
                    case 3
                        methodName='Trimmed_Mean';
                end;
                suffix='_avg';
            case 3
                averagingMethod='Jitter-Correct';
                switch methodNum
                    case 1
                        methodName='Mean';
                    case 2
                        methodName='Median';
                    case 3
                        methodName='Trimmed_Mean';
                end;
                suffix='_avg';
                jitterChan=EPdataset.dataset(EPmain.average.channelDataset).chanNames{EPmain.average.channel};
            case 4
                averagingMethod='Frequency-Coherence';
                switch methodNum
                    case 1
                        methodName='multi-taper';
                    case 2
                        methodName='Hanning';
                end;
                suffix='_coh';
            case 5
                averagingMethod='Frequency-Phase Lock';
                methodName='phase-lock wavelet';
                suffix='_plv';
        end;
        
        smoothing=0;
        if methodNum > 3 %Frequency domain analysis
            smoothing=str2num(get(EPmain.handles.average.smoothing,'String'));
        end;

        EPmain.average.importFormat=importFormatNum;
        EPmain.average.type=typeNum;
        EPmain.average.outputFormat=outputFormatNum;
        
        minLatency=[];
        maxLatency=[];
        latencyName=[];
        if ismember(procedureNum,[2,3])
            minLatency=EPmain.average.minLatency;
            maxLatency=EPmain.average.maxLatency;
            if minLatency > maxLatency
                msg{1}='Latency minimum cannot be larger than the latency maximum.';
                [msg]=ep_errorMsg(msg);
                return
            end
            if procedureNum==2
                latencyName=EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames{EPmain.average.trialSpec};
            end;
        end;
        
        multiSessionSubject=[];
        if ~isempty(EPmain.average.subject)
            if ~isempty(findstr('-',EPmain.average.subject))
                msg{1}='The subject field cannot have a dash.  If the intent was to specify a range of numbers, a colon is needed.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            multiSessionSubject=str2num(EPmain.average.subject);
            if ~isnumeric(multiSessionSubject)
                msg{1}='The subject field needs to be a number or numbers or empty.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end
        
        [sessionFiles, activeDirectory]=ep_getFilesUI(importFormat);
        if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
            msg{1}='No filenames selected. You have to click on a name.';
            [msg]=ep_errorMsg(msg);
            return
        end
        for iFile=1:size(sessionFiles,2)
            if ~isempty(multiSessionSubject)
                [pathstr, name, ext] = fileparts(sessionFiles{iFile});
                if max(multiSessionSubject) > length(name)
                    msg{1}=['The file name ' name ' is shorter than specified by the Subject field.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end
            end;
            sessionFiles{iFile}=[activeDirectory sessionFiles{iFile}];
        end;
        
        sessionFiles=sort(sessionFiles);
        
        [pathstr, name, ext] = fileparts(sessionFiles{1});
        
        if strcmp(name(end-3:end),'_seg')
            name=name(1:end-4);
        end;
        outFileName=[pathstr filesep name suffix ext];

        [outFileName, pathname] = uiputfile('*.*','Save:',outFileName);
        if outFileName == 0
            msg{1}='No output name selected.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        outFileName=[pathname outFileName];
        if exist(outFileName,'file')
            delete(outFileName); %user must have clicked "yes" to whether to replace existing file
        end;
        
        set(EPmain.handles.average.average,'enable','off');
        set(EPmain.handles.average.done,'enable','off');
        drawnow
        averagedData=ep_averageData(sessionFiles,importFormat,dataType,averagingMethod,EPmain.preferences.average.trimLevel,methodName,smoothing,latencyName,minLatency,maxLatency,jitterChan,EPmain.average.peakPolarity,multiSessionSubject);
        if isfield(averagedData,'data')
            if EPmain.average.dropEvents
                averagedData.events=cell(size(averagedData.events));
            end;
            [err]=ep_writeData(averagedData,outFileName,EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,outputFormat);
        end
        try
            set(EPmain.handles.average.average,'enable','on');
            set(EPmain.handles.average.done,'enable','on');
        catch
            ep('start')
        end;
        
        drawnow
        
        ep('start');
        
    case 'startRead'
        
        set(EPmain.handles.hMainWindow,'Name', 'Read Data');
        
        uicontrol('Style','text',...
            'String','File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 450 100 20]);
        
        EPmain.handles.read.format = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.read.format,''Value'');','if tempVar ~=0,EPmain.read.format=tempVar;end;','if isempty(tempVar),EPmain.read.format=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.read.format,'Position',[20 420 150 20]);
        
        uicontrol('Style','text',...
            'String','File Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 390 100 20]);
        
        EPmain.handles.read.type= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.read.type,''Value'');','if tempVar ~=0,EPmain.read.type=tempVar;end;','if isempty(tempVar),EPmain.read.type=tempVar;end;','ep(''start'');'],...
            'Value',EPmain.read.type,'Position',[20 360 150 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.read.format});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.read.type,'enable','off');
        end;
        
        uicontrol('Style','frame',...
            'Position',[20 265 170 90]);
        
        EPmain.handles.read.check= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Single Cell Files',...
            'CallBack',['global EPmain;','EPmain.read.check=get(EPmain.handles.read.check,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Value',EPmain.read.check,'Position',[30 330 150 20]);
        
        EPmain.handles.read.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.read.subject,...
            'CallBack',['global EPmain;','EPmain.read.subject=get(EPmain.handles.read.subject,''String'');','ep(''start'');'],...
            'Position',[30 310 50 20],'TooltipString','example 4:6');
        
        EPmain.handles.read.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Subject','HorizontalAlignment','left',...
            'Position',[90 310 50 20]);
        
        if isempty(EPmain.read.subject)
            set(EPmain.handles.read.subjectLabel,'enable','off');
        elseif isempty(str2num(EPmain.read.subject))
            set(EPmain.handles.read.subjectLabel,'enable','off');
        end;
        
        EPmain.handles.read.cell= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.read.cell,...
            'CallBack',['global EPmain;','EPmain.read.cell=get(EPmain.handles.read.cell,''String'');','ep(''start'');'],...
            'Position',[30 290 50 20],'TooltipString','example 7:9');
        
        EPmain.handles.read.cellLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Cell','HorizontalAlignment','left',...
            'Position',[90 290 50 20]);

        if isempty(EPmain.read.cell)
            set(EPmain.handles.read.cellLabel,'enable','off');
        elseif isempty(str2num(EPmain.read.cell))
            set(EPmain.handles.read.cellLabel,'enable','off');
        end;
        
        EPmain.handles.read.freq= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.read.freq,...
            'CallBack',['global EPmain;','EPmain.read.freq=get(EPmain.handles.read.freq,''String'');','ep(''start'');'],...
            'Position',[30 270 50 20],'TooltipString','example 10:12');
        
        EPmain.handles.read.freqLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Freq','HorizontalAlignment','left',...
            'Position',[90 270 50 20]);
        
        if isempty(EPmain.read.freq)
            set(EPmain.handles.read.freqLabel,'enable','off');
        elseif isempty(str2num(EPmain.read.freq))
            set(EPmain.handles.read.freqLabel,'enable','off');
        end;
        
        if ~isempty(EPdataset.dataset)
            for i=1:length(EPdataset.dataset)
                fileName=EPdataset.dataset(i).dataName;
                if strcmp(EPdataset.dataset(i).saved,'no')
                    fileName=['*' fileName];
                end;
                tableData{i,1}=fileName;
            end;
        else
            tableData=[];
        end;
        
        tableNames{1}='data';
        columnEditable =  false;
        ColumnFormat{1}=[];
        
        EPmain.handles.read.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'CellSelectionCallback',@deleteReadData,'ForegroundColor','red',...
            'ColumnWidth',{300},'Position',[20 100 170 150]);
        
        EPmain.handles.read.hRead = uicontrol('Style', 'pushbutton', 'String', 'Read','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', 'ep(''readData'')');
        
        EPmain.handles.read.hQuitRead = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'readData'
        
        set(EPmain.handles.read.subject,'ForegroundColor','black');
        set(EPmain.handles.read.cell,'ForegroundColor','black');
        
        typeNum = get(EPmain.handles.read.type,'value');
        switch typeNum
            case 1
                dataType='continuous';
            case 2
                dataType='single_trial';
            case 3
                dataType='average';
            case 4
                dataType='grand_average';
            case 5
                dataType='factors';
        end;
        
        formatNum = get(EPmain.handles.read.format,'value');
        [importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{formatNum});
        
        if strcmp(importFormat,'ep_mat')
            dataType=''; %data type set by data file itself
        end;
        
        EPmain.read.format=formatNum;
        EPmain.read.type=typeNum;
        
        set(EPmain.handles.read.hRead,'enable','off');
        set(EPmain.handles.read.hQuitRead,'enable','off');
        drawnow
        
        EPmain.convertMode=0;
        theHandles=EPmain.handles.read;
        readFiles(theHandles,importFormat,dataType);

        try
            set(EPmain.handles.read.hRead,'enable','on');
            set(EPmain.handles.read.hQuitRead,'enable','on');
            drawnow
        catch
            ep('start')
        end;
        
        ep('start');
        
    case 'startEdit'
        
        if strcmp(EPoverview,'done')
            doneFlag=1;
        else
            doneFlag=0;
        end;
        EPoverview=[];
        tableData=[];
        
        if ~isempty(EPdataset.dataset)
            for i=1:length(EPdataset.dataset)
                fileName=EPdataset.dataset(i).dataName;
                if strcmp(EPdataset.dataset(i).saved,'no')
                    fileName=['*' fileName];
                end;
                tableData{i,1}=fileName;
            end;
        else
            tableData=[];
        end;
        
        tableNames{1}='data';
        columnEditable =  false;
        ColumnFormat{1}=[];
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Click on dataset to edit.','FontSize',EPmain.fontsize,...
            'Position',[20 250 160 20]);
        
        EPmain.handles.edit.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'CellSelectionCallback',@pickEditData,...
            'ColumnWidth',{300},'Position',[20 100 170 150]);
        
        EPmain.handles.edit.hQuitRead = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
        if length(EPdataset.dataset) == 1 && ~doneFlag %if just one dataset, don't need to wait for it to be selected.
            EPoverview.dataset=1;
            EPoverview.mode=[];
            ep_editData
        end;
        
    case 'startView'
        
        set(EPmain.handles.hMainWindow,'Name', 'View EEG');
        refresh
        
        if isempty(EPmain.view)
            
            %no option to switch to this panel unless there are data to view.
            for i=1:length(EPdataset.dataset)
                EPmain.view.theFileNames{i,1}=EPdataset.dataset(i).dataName;
                [uniqueCells waveOrder m]=unique(EPdataset.dataset(i).cellNames,'first'); %waveOrder is the true order of the unique cells (wave numbers of the first occurences)
                [a b cellOrder]=unique(waveOrder,'first'); %cellOrder is the true order of the unique cells (numbered by cells)
                EPmain.view.theCells{i}(cellOrder)=uniqueCells;
            end;
            
            EPmain.view.trialList=cell(4,1);
            EPmain.view.allTrials=zeros(4,1); %1=all trials/subs,2=trial/subs erpimage,3=all facs,4=fac erpimage,5=all cells,6=cells erpimage
            EPmain.view.correl=zeros(4,1);
            EPmain.view.STS=zeros(4,1);
            EPmain.view.rel=zeros(4,1);
            %initially set up the four colors with one for each cell
            for iColor=1:4
                if iColor <= length(unique(EPdataset.dataset(end).cellNames))                    
                    EPmain.view.dataset(iColor,1)=length(EPdataset.dataset);
                    EPmain.view.cell(iColor,1)=min(iColor,length(EPmain.view.theCells{EPmain.view.dataset(iColor)}));
                    EPmain.view.trial(iColor,1)=1;
                    EPmain.view.subject(iColor,1)=1;
                    if ~isempty(EPdataset.dataset(end).facNames)
                        EPmain.view.subject(iColor)=length(EPdataset.dataset(end).subNames);
                    end;
                    gavSubs=find(strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes,'GAV'));
                    if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average') && (length(gavSubs)==1)
                        EPmain.view.subject(iColor)=gavSubs;
                    end;
                    if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                        EPmain.view.trialList{iColor}=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
                    end;
                    EPmain.view.factor(iColor,1)=1;
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                        EPmain.view.factor(iColor,1)=length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)+1;
                    end;
                    
                    EPmain.view.rel(iColor,1)=~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).relNames);
                    EPmain.view.correl(iColor)=0;
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames)
                        if EPmain.view.rel(iColor)
                            EPmain.view.correl(iColor)=1; %coherence measure is correlation and phase-lock is unitless 0 to 1
                        end;
                    end;
                    if strcmp('STS',EPdataset.dataset(EPmain.view.dataset(iColor)).cellTypes{EPmain.view.cell(iColor)})
                        EPmain.view.STS(iColor)=1;
                    else
                        EPmain.view.STS(iColor)=0;
                    end;
                else
                    EPmain.view.dataset(iColor,1)=length(EPdataset.dataset)+1; %all the cells shown so just leave this one blank.
                end;
                EPmain.view.changeDatasetFlag(iColor)=0;
                EPmain.view.changeFlag(iColor)=1;
            end;
            EPmain.view.FFTunits=4;
            EPmain.view.marker1=[];
            EPmain.view.marker2=[];
            EPmain.view.plotMVmin=zeros(4,1);
            EPmain.view.plotMVmax=zeros(4,1);
            EPmain.view.plotMVminAbs=zeros(4,1);
            EPmain.view.plotMVmaxAbs=zeros(4,1);
            EPmain.view.startSamp=0;
            EPmain.view.endSamp=0;
            EPmain.view.startHz=0;
            EPmain.view.endHz=0;
            EPmain.view.binHz=0;
            EPmain.view.edited.bottomVolt=0;
            EPmain.view.edited.topVolt=0;
            EPmain.view.edited.startSamp=0;
            EPmain.view.edited.endSamp=0;
            EPmain.view.edited.startHz=0;
            EPmain.view.edited.endHz=0;
            EPmain.view.manual.bottomVolt=0;
            EPmain.view.manual.topVolt=0;
            EPmain.view.manual.startSamp=0;
            EPmain.view.manual.endSamp=0;
            EPmain.view.manual.startHz=0;
            EPmain.view.manual.endHz=0;
            EPmain.view.flexMode=strcmp(EPdataset.dataset(end).timeUnits,'per');
             
            if ~isempty(EPdataset.dataset(end).freqNames)
                if ~isempty(EPdataset.dataset(end).timeNames)
                    EPmain.view.dataTransform='TFT';
                else
                    EPmain.view.dataTransform='FFT';
                end
            else
                EPmain.view.dataTransform='VLT';
            end;
            
            EPmain.view.dataType=EPdataset.dataset(end).dataType;
                        
            EPmain.view.eventList=cell(0);
            EPmain.view.RT=0;
            for iColor=1:4
                if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                    if isempty(EPmain.view.eventList) || ~any(ismember(EPmain.view.dataset(iColor),EPmain.view.dataset(1:iColor-1)))
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                            theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{EPmain.view.subject,EPmain.view.trial};
                        else
                            theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{EPmain.view.subject,EPmain.view.cell};
                        end;
                        typeValues=cell(length(theEvents),1);
                        for iEvent=1:length(theEvents)
                            typeValues{iEvent}=[num2str(theEvents(iEvent).type) '-' num2str(theEvents(iEvent).value)];
                        end;
                        if any(strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).trialSpecNames,'RT'))
                            EPmain.view.RT=1;
                        end;
                    end;
                end;
                EPmain.view.eventList(end+1:end+length(theEvents),1)=typeValues;
            end;
            EPmain.view.eventList=unique(EPmain.view.eventList);
            EPmain.view.events=length(EPmain.view.eventList)+1;
            if EPmain.view.RT
                EPmain.view.events=EPmain.view.events+1; %add -RT-
            end;
            if EPmain.view.events > 2
                EPmain.view.events=EPmain.view.events+1;  %add -all-
            end;
        end;
        
        if any(EPmain.view.changeDatasetFlag) %if the dataset was just changed
            theColor=min(find(EPmain.view.changeDatasetFlag));
            if EPmain.view.dataset(theColor) <= length(EPdataset.dataset)
                if strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).dataType,'single_trial')
                    EPmain.view.cell(theColor)=1;
                    EPmain.view.trialList{theColor}=EPdataset.dataset(EPmain.view.dataset(theColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(theColor)}(EPmain.view.cell(theColor)),EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames));
                else
                    EPmain.view.cell(theColor)=min(theColor,length(EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames));
                end;
                EPmain.view.subject(theColor)=1;
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).facNames)
                    EPmain.view.subject(theColor)=length(EPdataset.dataset(EPmain.view.dataset(theColor)).subNames);
                end;
                gavSubs=find(strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).subTypes,'GAV'));
                if strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).dataType,'average') && (length(gavSubs)==1)
                    EPmain.view.subject(theColor)=gavSubs;
                end;
                EPmain.view.factor(theColor)=1;
                
                EPmain.view.trial(theColor)=1;
                
                EPmain.view.rel(theColor)=~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).relNames);
                
                %if changed transform type (voltage, FFT, or TFT), then change over all the other colors too.
                EPmain.view.correl(theColor)=0;
                if EPmain.view.rel(theColor)
                    EPmain.view.correl(theColor)=1; %coherence measure is correlation and phase-lock is unitless 0 to 1
                end;
                if strcmp('STS',EPdataset.dataset(EPmain.view.dataset(theColor)).cellTypes{EPmain.view.cell(theColor)})
                    EPmain.view.STS(theColor)=1;
                else
                    EPmain.view.STS(theColor)=0;
                end;
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).freqNames)
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).timeNames)
                        newType='TFT';
                    else
                        newType='FFT';
                    end
                else
                    newType='VLT';
                end;
                if ~strcmp(newType,EPmain.view.dataTransform) || (xor(strcmp(EPmain.view.dataType,'continuous'),strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).dataType,'continuous'))) || (xor(strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).timeUnits,'per'),EPmain.view.flexMode))
                    %if new dataset is different from existing ones then change them too to keep consistent
                    disp('The new dataset is different in type than current datasets so all four colors will all be replaced with the new dataset.');
                    EPmain.view.dataTransform=newType;
                    EPmain.view.dataType=EPdataset.dataset(EPmain.view.dataset(theColor)).dataType;
                    EPmain.view.flexMode=strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).timeUnits,'per');
                    for iColor=1:4
                        if iColor ~= theColor
                            if iColor <= length(EPmain.view.theCells{EPmain.view.dataset(theColor)})
                                EPmain.view.dataset(iColor)=EPmain.view.dataset(theColor);
                                EPmain.view.cell(iColor)=min(iColor,length(EPmain.view.theCells{EPmain.view.dataset(iColor)}));
                                if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                                    EPmain.view.trialList{iColor}=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
                                end;
                                EPmain.view.trial(iColor)=1;
                                EPmain.view.subject(iColor)=1;
                                if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).facNames)
                                    EPmain.view.subject(iColor)=length(EPdataset.dataset(theColor).subNames);
                                end;
                                EPmain.view.factor(iColor)=1;
                                EPmain.view.changeFlag(iColor)=1;
                                EPmain.view.correl(iColor)=EPmain.view.correl(theColor);
                                EPmain.view.STS(iColor)=EPmain.view.STS(theColor);
                                EPmain.view.rel(iColor)=EPmain.view.rel(theColor);
                            else
                                EPmain.view.dataset(iColor)=length(EPdataset.dataset)+1; %all the cells shown so just leave this one blank.
                            end;
                        end;
                    end;
%                     EPmain.view.eventList=cell(0);
%                     EPmain.view.RT=0;
%                     for iColor=1:4
%                         if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
%                             if isempty(EPmain.view.eventList) || ~any(ismember(EPmain.view.dataset(iColor),EPmain.view.dataset(1:iColor-1)))
%                                 if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
%                                     theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{EPmain.view.subject,EPmain.view.trial};
%                                 else
%                                     theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{EPmain.view.subject,EPmain.view.cell};
%                                 end;
%                                 typeValues=cell(length(theEvents),1);
%                                 for iEvent=1:length(theEvents)
%                                     typeValues{iEvent}=[num2str(theEvents(iEvent).type) '-' num2str(theEvents(iEvent).value)];
%                                 end;
%                                 if any(strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).trialSpecNames,'RT'))
%                                     EPmain.view.RT=1;
%                                 end;
%                             end;
%                         end;
%                         EPmain.view.eventList(end+1:end+length(theEvents),1)=typeValues;
%                     end;
%                     EPmain.view.eventList=unique(EPmain.view.eventList);
%                     EPmain.view.events=length(EPmain.view.eventList)+1;
%                     if EPmain.view.RT
%                         EPmain.view.events=EPmain.view.events+1; %add -RT-
%                     end;
%                     if EPmain.view.events > 2
%                         EPmain.view.events=EPmain.view.events+1;  %add -all-
%                     end;
                end;
            else
                EPmain.view.correl(theColor)=0;
                EPmain.view.rel(theColor)=0;
                EPmain.view.STS(theColor)=0;
            end;
            
            EPmain.view.eventList=cell(0);
            EPmain.view.RT=0;
            %set up the list of events
            for iColor=1:4
                if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                    if isempty(EPmain.view.eventList) || ~any(ismember(EPmain.view.dataset(iColor),EPmain.view.dataset(1:iColor-1)))
                        if ismember(EPmain.view.allTrials(1),[5 6])
                            cellList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames)];
                        else
                            cellList=EPmain.view.cell(iColor);
                        end;
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                            if ismember(EPmain.view.allTrials(1),[1 2])
                                cellList=EPmain.view.trialList{iColor};
                            else
                                cellList=EPmain.view.trial(iColor);
                            end;
                            subList=EPmain.view.subject(iColor);
                        else
                            if ismember(EPmain.view.allTrials(1),[1 2])
                                subList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames)];
                            else
                                subList=EPmain.view.subject(iColor);
                            end;
                        end;
                        typeValues=cell(0);
                        for iSub=1:length(subList)
                            theSub=subList(iSub);
                            for iCell=1:length(cellList)
                                theCell=cellList(iCell);
                                if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                                    theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{theSub,theCell};
                                else
                                    theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{theSub,theCell};
                                end;
                                for iEvent=1:length(theEvents)
                                    typeValues{end+1}=[num2str(theEvents(iEvent).type) '-' num2str(theEvents(iEvent).value)];
                                end;
                            end;
                        end;
                        if any(strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).trialSpecNames,'RT'))
                            EPmain.view.RT=1;
                        end;
                    end;
                end;
                EPmain.view.eventList(end+1:end+length(typeValues),1)=typeValues;
            end;
            EPmain.view.eventList=unique(EPmain.view.eventList);
            EPmain.view.events=length(EPmain.view.eventList)+1;
            if EPmain.view.RT
                EPmain.view.events=EPmain.view.events+1; %add -RT-
            end;
            if EPmain.view.events > 2
                EPmain.view.events=EPmain.view.events+1;  %add -all-
            end;
            
            EPmain.view.plotMVmin(theColor)=0;
            EPmain.view.plotMVmax(theColor)=0;
            EPmain.view.plotMVminAbs(theColor)=0;
            EPmain.view.plotMVmaxAbs(theColor)=0;
            EPmain.view.changeDatasetFlag(theColor)=0;
            EPmain.view.changeFlag(theColor)=1;
            EPmain.view.edited.bottomVolt=0;
            EPmain.view.edited.topVolt=0;
            EPmain.view.edited.startSamp=0;
            EPmain.view.edited.endSamp=0;
            EPmain.view.edited.startHz=0;
            EPmain.view.edited.endHz=0;
            EPmain.view.manual.bottomVolt=0;
            EPmain.view.manual.topVolt=0;
            EPmain.view.manual.startSamp=0;
            EPmain.view.manual.endSamp=0;
            EPmain.view.manual.startHz=0;
            EPmain.view.manual.endHz=0;
        end;
        
        marker1=num2str(EPmain.view.marker1);
        marker2=num2str(EPmain.view.marker2);
        
        %compute voltages and epoch for the changed colors
        if any(EPmain.view.changeFlag)
            for color=1:4
                if EPmain.view.changeFlag(color)
                    EPmain.view.changeFlag(color)=0;
                    EPmain.view.allTrials(color)=0;
                    
                    if EPmain.view.dataset(color) <= length(EPdataset.dataset)
                        
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(color)).dataType,'continuous') %continuous data
                            theCell=EPmain.view.trial(color);
                            EPmain.view.plotMVmin(color)=EPdataset.dataset(EPmain.view.dataset(color)).plotMVmin(theCell,1,EPmain.view.factor(color));
                            EPmain.view.plotMVmax(color)=EPdataset.dataset(EPmain.view.dataset(color)).plotMVmax(theCell,1,EPmain.view.factor(color));
                            EPmain.view.plotMVminAbs(color)=EPdataset.dataset(EPmain.view.dataset(color)).plotMVminAbs(theCell,1,EPmain.view.factor(color));
                            EPmain.view.plotMVmaxAbs(color)=EPdataset.dataset(EPmain.view.dataset(color)).plotMVmaxAbs(theCell,1,EPmain.view.factor(color));
                            
                        elseif strcmp(EPdataset.dataset(EPmain.view.dataset(color)).dataType,'average') %averaged data
                            numSubs=length(EPdataset.dataset(EPmain.view.dataset(color)).subNames);
                            if EPmain.view.subject(color) > numSubs
                                theSubject=[1:numSubs];
                                if EPmain.view.subject(color) == numSubs+1
                                    EPmain.view.allTrials(color)=1; %'all'
                                end;
                                if EPmain.view.subject(color) == numSubs+2
                                    EPmain.view.allTrials(color)=2; %'erpimage'
                                end;
                            else
                                theSubject=EPmain.view.subject(color);
                            end;
                            numFacs=length(EPdataset.dataset(EPmain.view.dataset(color)).facNames);
                            if (EPmain.view.factor(color) > numFacs) && (numFacs > 0)
                                theFactor=[1:numFacs];
                                if EPmain.view.factor(color) == numFacs+1
                                    EPmain.view.allTrials(color)=3; %'all'
                                end;
                                if EPmain.view.factor(color) == numFacs+2
                                    EPmain.view.allTrials(color)=4; %'erpimage'
                                end;
                            else
                                theFactor=EPmain.view.factor(color);
                            end;
                            numCells=length(EPdataset.dataset(EPmain.view.dataset(color)).cellNames);
                            theCell=EPmain.view.cell(color);
                            if theCell > numCells
                                theCell=[1:numCells];
                                if EPmain.view.cell(color) == numCells+1
                                    EPmain.view.allTrials(color)=5; %'all'
                                end;
                                if EPmain.view.cell(color) == numCells+2
                                    EPmain.view.allTrials(color)=6; %'erpimage'
                                end;
                            end;
                            EPmain.view.plotMVmin(color)=min(EPdataset.dataset(EPmain.view.dataset(color)).plotMVmin(theCell,theSubject,theFactor));
                            EPmain.view.plotMVmax(color)=max(EPdataset.dataset(EPmain.view.dataset(color)).plotMVmax(theCell,theSubject,theFactor));
                            EPmain.view.plotMVminAbs(color)=min(EPdataset.dataset(EPmain.view.dataset(color)).plotMVminAbs(theCell,theSubject,theFactor));
                            EPmain.view.plotMVmaxAbs(color)=max(EPdataset.dataset(EPmain.view.dataset(color)).plotMVmaxAbs(theCell,theSubject,theFactor));
                            
                        else  %if single_trial data
                            EPmain.view.trialList{color}=EPdataset.dataset(EPmain.view.dataset(color)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(color)}(EPmain.view.cell(color)),EPdataset.dataset(EPmain.view.dataset(color)).cellNames));
                            if EPmain.view.trial(color) > length(EPmain.view.trialList{color})
                                trialList=EPmain.view.trialList{color};
                                if EPmain.view.trial(color) == length(EPmain.view.trialList{color})+1
                                    EPmain.view.allTrials(color)=1; %'all'
                                end;
                                if EPmain.view.trial(color) == length(EPmain.view.trialList{color})+2
                                    EPmain.view.allTrials(color)=2; %'erpimage'
                                end;
                            else
                                trialList=EPmain.view.trial(color);
                            end;
                            numFacs=length(EPdataset.dataset(EPmain.view.dataset(color)).facNames);
                            if (EPmain.view.factor(color) > numFacs) && (numFacs > 0)
                                theFactor=[1:numFacs];
                                if EPmain.view.factor(color) == numFacs+1
                                    EPmain.view.allTrials(color)=3; %'all'
                                end;
                                if EPmain.view.factor(color) == numFacs+2
                                    EPmain.view.allTrials(color)=4; %'erpimage'
                                end;
                            else
                                theFactor=EPmain.view.factor(color);
                            end;
                            trialsInCell=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(color)}(EPmain.view.cell(color)),EPdataset.dataset(EPmain.view.dataset(color)).cellNames));
                            EPmain.view.plotMVmin(color)=min(EPdataset.dataset(EPmain.view.dataset(color)).plotMVmin(trialsInCell(trialList),1,theFactor));
                            EPmain.view.plotMVmax(color)=max(EPdataset.dataset(EPmain.view.dataset(color)).plotMVmax(trialsInCell(trialList),1,theFactor));
                            EPmain.view.plotMVminAbs(color)=min(EPdataset.dataset(EPmain.view.dataset(color)).plotMVminAbs(trialsInCell(trialList),1,theFactor));
                            EPmain.view.plotMVmaxAbs(color)=max(EPdataset.dataset(EPmain.view.dataset(color)).plotMVmaxAbs(trialsInCell(trialList),1,theFactor));
                        end;
                    end;
                end;
            end;
            
            
            if strcmp('FFT',EPmain.view.dataTransform)
                EPmain.view.startSamp=NaN;
                EPmain.view.endSamp=NaN;
            else
                EPmain.view.startSamp=0;
                EPmain.view.endSamp=0;
                % Time is counted for an epoch as, say, -200 to 800 ms.  In this case the samples have some length to them depending on
                % the digitization rate, such as 4ms for a 250Hz sampling rate.  The points prior to the baseline are counted from the
                % left side of the baseline sample and the points following the baseline are counted from the right side of the baseline
                % sample.  Thus, when the actual samples are calculated, one first adds the baseline to the numbers, yielding 0 to
                % 1000ms.  One then must add the length of the baseline to the first number to obtain 4 to 1000ms.  Then one would
                % divide the numbers by the sample size, resulting in 1 to 250 samples, which are the correct samples to use.
                % For continuous data or other such data where there is no baseline, one must by this convention start with 0 ms,
                % as in 0-1000ms and 1000-2000ms.
                for color=1:4
                    if EPmain.view.dataset(color) <= length(EPdataset.dataset)
                        if ~isnan(EPmain.view.startSamp) && ~isnan(EPmain.view.endSamp)
                            if ~EPmain.view.startSamp && ~EPmain.view.endSamp %if neither has been set yet
                                if strcmp(EPdataset.dataset(EPmain.view.dataset(color)).dataType,'continuous') %continuous data
                                    EPmain.view.startSamp=1000*(EPmain.view.trial(color)-1); %for continuous data, each "trial" is 1000 ms
                                    EPmain.view.endSamp=1000*(EPmain.view.trial(color));
                                else
                                    EPmain.view.startSamp=min(EPdataset.dataset(EPmain.view.dataset(color)).timeNames);
                                    EPmain.view.endSamp=max(EPdataset.dataset(EPmain.view.dataset(color)).timeNames)+(1000/EPdataset.dataset(EPmain.view.dataset(color)).Fs);
                                end;
                            else
                                if strcmp(EPdataset.dataset(EPmain.view.dataset(color)).dataType,'continuous') %continuous data
                                    if EPmain.view.startSamp ~= 1000*(EPmain.view.trial(color)-1)
                                        EPmain.view.startSamp=NaN;
                                    end;
                                    if EPmain.view.endSamp ~= 1000*(EPmain.view.trial(color))
                                        EPmain.view.endSamp=NaN;
                                    end;
                                else
                                    if EPmain.view.startSamp ~= min(EPdataset.dataset(EPmain.view.dataset(color)).timeNames)
                                        EPmain.view.startSamp=NaN;
                                    end;
                                    if (EPmain.view.endSamp ~= max(EPdataset.dataset(EPmain.view.dataset(color)).timeNames+(1000/EPdataset.dataset(EPmain.view.dataset(color)).Fs)))
                                        EPmain.view.endSamp=NaN;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
            
            if strcmp('VLT',EPmain.view.dataTransform)
                EPmain.view.startHz=NaN;
                EPmain.view.endHz=NaN;
            else
                EPmain.view.startHz=0;
                EPmain.view.endHz=0;
                for color=1:4
                    if EPmain.view.dataset(color) <= length(EPdataset.dataset)
                        if ~isnan(EPmain.view.startHz) && ~isnan(EPmain.view.endHz)
                            if ~EPmain.view.startHz && ~EPmain.view.endHz %if neither has been set yet
                                EPmain.view.startHz=min(EPdataset.dataset(EPmain.view.dataset(color)).freqNames);
                                EPmain.view.endHz=max(EPdataset.dataset(EPmain.view.dataset(color)).freqNames);
                                EPmain.view.binHz=EPdataset.dataset(EPmain.view.dataset(color)).freqNames(2)-EPdataset.dataset(EPmain.view.dataset(color)).freqNames(1); %making assumption all have the same bin resolution
                            else
                                if EPmain.view.startHz ~= min(EPdataset.dataset(EPmain.view.dataset(color)).freqNames)
                                    EPmain.view.startHz=NaN;
                                end;
                                if EPmain.view.endHz ~= max(EPdataset.dataset(EPmain.view.dataset(color)).freqNames)
                                    EPmain.view.endHz=NaN;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
        if all(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
            plotMVmin=min(EPmain.view.plotMVmin(EPmain.view.dataset <= length(EPdataset.dataset) & ~EPmain.view.STS));
            plotMVmax=max(EPmain.view.plotMVmax(EPmain.view.dataset <= length(EPdataset.dataset) & ~EPmain.view.STS));
        else
            %any any datasets are not correlations, then ignore correlations when calculating min and max values
            if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'})) && (EPmain.view.FFTunits > 1)
                plotMVmin=min(EPmain.view.plotMVminAbs(EPmain.view.dataset <= length(EPdataset.dataset) & ~EPmain.view.correl & ~EPmain.view.STS));
                plotMVmax=max(EPmain.view.plotMVmaxAbs(EPmain.view.dataset <= length(EPdataset.dataset) & ~EPmain.view.correl & ~EPmain.view.STS));
            else
                plotMVmin=min(EPmain.view.plotMVmin(EPmain.view.dataset <= length(EPdataset.dataset) & ~EPmain.view.correl & ~EPmain.view.STS));
                plotMVmax=max(EPmain.view.plotMVmax(EPmain.view.dataset <= length(EPdataset.dataset) & ~EPmain.view.correl & ~EPmain.view.STS));
            end;
        end
        
        menuStart=20;
        menuX=150;
        menuY=18;
        %blue        
        theDatasets=[EPmain.view.theFileNames; {'none'}];        
        EPmain.handles.view.dataset(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theDatasets,'ForegroundColor','blue',...
            'Value',EPmain.view.dataset(1),'Position',[menuStart 500-menuY menuX menuY],...
            'Callback', ['pause(.2);','global EPmain;','tempVar=get(EPmain.handles.view.dataset(1),''value'');','if tempVar ~=0,EPmain.view.dataset(1)=tempVar;end;','if isempty(tempVar),EPmain.view.dataset(1)=tempVar;end;','EPmain.view.changeDatasetFlag(1)=1;','ep(''start'')']);
        
        if EPmain.view.dataset(1) <= length(EPdataset.dataset)
            theCells=EPmain.view.theCells{EPmain.view.dataset(1)};
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theCells) > 1 && strcmp(EPdataset.dataset(EPmain.view.dataset(1)).dataType,'average') && ~ismember(EPmain.view.allTrials(1),[1 2 3 4])
                theCells{end+1}='-all-';
                theCells{end+1}='-erpimage-';
            end;
            EPmain.handles.view.cell(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theCells,'ForegroundColor','blue',...
                'Value',EPmain.view.cell(1),'Position',[menuStart 500-menuY*2 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.cell(1),''value'');','if tempVar ~=0,EPmain.view.cell(1)=tempVar;end;','if isempty(tempVar),EPmain.view.cell(1)=tempVar;end;','EPmain.view.changeFlag(1)=1;','ep(''start'')']);
            
            theSubs=EPdataset.dataset(EPmain.view.dataset(1)).subNames;
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theSubs) > 1 && ~ismember(EPmain.view.allTrials(1),[3 4 5 6])
                theSubs{end+1}='-all-';
                theSubs{end+1}='-erpimage-';
            end;
            EPmain.handles.view.subject(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theSubs,'ForegroundColor','blue',...
                'Value',EPmain.view.subject(1),'Position',[menuStart 500-menuY*3 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.subject(1),''value'');','if tempVar ~=0,EPmain.view.subject(1)=tempVar;end;','if isempty(tempVar),EPmain.view.subject(1)=tempVar;end;','EPmain.view.changeFlag(1)=1;','ep(''start'')']);
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(1)).trialNames) 
                theTrials=num2cell(EPmain.view.trialList{1});
                if ~strcmp('TFT',EPmain.view.dataTransform) && ~ismember(EPmain.view.allTrials(1),[3 4 5 6])
                    theTrials{end+1}='-all-';
                    theTrials{end+1}='-erpimage-';
                end;
                EPmain.handles.view.trial(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor','blue',...
                    'Value',EPmain.view.trial(1),'Position',[menuStart 500-menuY*4 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(1),''value'');','if tempVar ~=0,EPmain.view.trial(1)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(1)=tempVar;end;','EPmain.view.changeFlag(1)=1;','ep(''start'')']);
            elseif strcmp(EPdataset.dataset(EPmain.view.dataset(1)).dataType,'continuous')
                theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(1)).timeNames)/ceil(EPdataset.dataset(EPmain.view.dataset(1)).Fs))]'));
                for i=1:length(theEpochs)
                    theEpochs{i}=['Sec: ' theEpochs{i}];
                end;
                EPmain.handles.view.trial(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theEpochs,'ForegroundColor','blue',...
                    'Value',EPmain.view.trial(1),'Position',[menuStart 500-menuY*4 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(1),''value'');','if tempVar ~=0,EPmain.view.trial(1)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(1)=tempVar;end;','EPmain.view.changeFlag(1)=1;','ep(''start'')']);
            elseif strcmp('VLT',EPmain.view.dataTransform)
                theTrials=cell(0);
                theTrials{end+1}='-none-';
                theTrials{end+1}='-GFP-';
                theTrials{end+1}='-noise-';
                theTrials{end+1}='-StDev-';
                theTrials{end+1}='-CI-';
                EPmain.handles.view.trial(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor','blue',...
                    'Value',EPmain.view.trial(1),'Position',[menuStart 500-menuY*4 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(1),''value'');','if tempVar ~=0,EPmain.view.trial(1)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(1)=tempVar;end;','EPmain.view.changeFlag(1)=1;','ep(''start'')']);
            else
                EPmain.handles.view.trial(1) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Trials','FontSize',EPmain.fontsize,...
                    'ForegroundColor','blue','Position',[menuStart 500-menuY*4 menuX menuY]);
            end;
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(1)).facNames)
                theFacs=EPdataset.dataset(EPmain.view.dataset(1)).facNames;
                if ~strcmp('TFT',EPmain.view.dataTransform) && length(theFacs) > 1 && ~ismember(EPmain.view.allTrials(1),[1 2 5 6])
                    theFacs{end+1}='-all-';
                    theFacs{end+1}='-erpimage-';
                end;
                EPmain.handles.view.factor(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theFacs,'ForegroundColor','blue',...
                    'Value',EPmain.view.factor(1),'Position',[menuStart 500-menuY*5 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.factor(1),''value'');','if tempVar ~=0,EPmain.view.factor(1)=tempVar;end;','if isempty(tempVar),EPmain.view.factor(1)=tempVar;end;','EPmain.view.changeFlag(1)=1;','ep(''start'')']);
            else
                EPmain.handles.view.factor(1) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Factors','FontSize',EPmain.fontsize,...
                    'ForegroundColor','blue','Position',[menuStart 500-menuY*5 menuX menuY]);
            end;
        end;
        
        %red
        theDatasets=[EPmain.view.theFileNames; {'none'}];
        EPmain.handles.view.dataset(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theDatasets,'ForegroundColor','red',...
            'Value',EPmain.view.dataset(2),'Position',[menuStart 500-menuY*6 menuX menuY],...
            'Callback', ['pause(.2);','global EPmain;','tempVar=get(EPmain.handles.view.dataset(2),''value'');','if tempVar ~=0,EPmain.view.dataset(2)=tempVar;end;','if isempty(tempVar),EPmain.view.dataset(2)=tempVar;end;','EPmain.view.changeDatasetFlag(2)=1;','ep(''start'')']);
        
        if EPmain.view.dataset(2) <= length(EPdataset.dataset)
            theCells=EPmain.view.theCells{EPmain.view.dataset(2)};
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theCells) > 1 && strcmp(EPdataset.dataset(EPmain.view.dataset(2)).dataType,'average') && ~ismember(EPmain.view.allTrials(2),[1 2 3 4])
                theCells{end+1}='-all-';
                theCells{end+1}='-erpimage-';
            end;
            EPmain.handles.view.cell(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theCells,'ForegroundColor','red',...
                'Value',EPmain.view.cell(2),'Position',[menuStart 500-menuY*7 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.cell(2),''value'');','if tempVar ~=0,EPmain.view.cell(2)=tempVar;end;','if isempty(tempVar),EPmain.view.cell(2)=tempVar;end;','EPmain.view.changeFlag(2)=1;','ep(''start'')']);
            
            theSubs=EPdataset.dataset(EPmain.view.dataset(2)).subNames;
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theSubs) > 1 && ~ismember(EPmain.view.allTrials(2),[3 4 5 6])
                theSubs{end+1}='-all-';
                theSubs{end+1}='-erpimage-';
            end;
            EPmain.handles.view.subject(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theSubs,'ForegroundColor','red',...
                'Value',EPmain.view.subject(2),'Position',[menuStart 500-menuY*8 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.subject(2),''value'');','if tempVar ~=0,EPmain.view.subject(2)=tempVar;end;','if isempty(tempVar),EPmain.view.subject(2)=tempVar;end;','EPmain.view.changeFlag(2)=1;','ep(''start'')']);
            
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(2)).trialNames)
                theTrials=num2cell(EPmain.view.trialList{2});
                if ~strcmp('TFT',EPmain.view.dataTransform) && ~ismember(EPmain.view.allTrials(2),[3 4 5 6])
                    theTrials{end+1}='-all-';
                    theTrials{end+1}='-erpimage-';
                end;
                EPmain.handles.view.trial(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor','red',...
                    'Value',EPmain.view.trial(2),'Position',[menuStart 500-menuY*9 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(2),''value'');','if tempVar ~=0,EPmain.view.trial(2)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(2)=tempVar;end;','EPmain.view.changeFlag(2)=1;','ep(''start'')']);
            elseif strcmp(EPdataset.dataset(EPmain.view.dataset(2)).dataType,'continuous')
                theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(2)).timeNames)/ceil(EPdataset.dataset(EPmain.view.dataset(2)).Fs))]'));
                for i=1:length(theEpochs)
                    theEpochs{i}=['Sec: ' theEpochs{i}];
                end;
                EPmain.handles.view.trial(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theEpochs,'ForegroundColor','red',...
                    'Value',EPmain.view.trial(2),'Position',[menuStart 500-menuY*9 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(2),''value'');','if tempVar ~=0,EPmain.view.trial(2)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(2)=tempVar;end;','EPmain.view.changeFlag(2)=1;','ep(''start'')']);
            elseif strcmp('VLT',EPmain.view.dataTransform)
                theTrials=cell(0);
                theTrials{end+1}='-none-';
                theTrials{end+1}='-GFP-';
                theTrials{end+1}='-noise-';
                theTrials{end+1}='-StDev-';
                theTrials{end+1}='-CI-';
                EPmain.handles.view.trial(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor','red',...
                    'Value',EPmain.view.trial(2),'Position',[menuStart 500-menuY*9 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(2),''value'');','if tempVar ~=0,EPmain.view.trial(2)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(2)=tempVar;end;','EPmain.view.changeFlag(2)=1;','ep(''start'')']);
            else
                EPmain.handles.view.trial(2) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Trials','FontSize',EPmain.fontsize,...
                    'ForegroundColor','red','Position',[menuStart 500-menuY*9 menuX menuY]);
            end;
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(2)).facNames)
                theFacs=EPdataset.dataset(EPmain.view.dataset(2)).facNames;
                if ~strcmp('TFT',EPmain.view.dataTransform) && length(theFacs) > 1 && ~ismember(EPmain.view.allTrials(2),[1 2 5 6])
                    theFacs{end+1}='-all-';
                    theFacs{end+1}='-erpimage-';
                end;
                EPmain.handles.view.factor(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theFacs,'ForegroundColor','red',...
                    'Value',EPmain.view.factor(2),'Position',[menuStart 500-menuY*10 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.factor(2),''value'');','if tempVar ~=0,EPmain.view.factor(2)=tempVar;end;','if isempty(tempVar),EPmain.view.factor(2)=tempVar;end;','EPmain.view.changeFlag(2)=1;','ep(''start'')']);
            else
                EPmain.handles.view.factor(2) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Factors','FontSize',EPmain.fontsize,...
                    'ForegroundColor','red','Position',[menuStart 500-menuY*10 menuX menuY]);
            end;
        end;
        
        %green
        theDatasets=[EPmain.view.theFileNames; {'none'}];
        EPmain.handles.view.dataset(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theDatasets,'ForegroundColor',[0 .5 0],...
            'Value',EPmain.view.dataset(3),'Position',[menuStart 500-menuY*11 menuX menuY],...
            'Callback', ['pause(.2);','global EPmain;','tempVar=get(EPmain.handles.view.dataset(3),''value'');','if tempVar ~=0,EPmain.view.dataset(3)=tempVar;end;','if isempty(tempVar),EPmain.view.dataset(3)=tempVar;end;','EPmain.view.changeDatasetFlag(3)=1;','ep(''start'')']);
        
        if EPmain.view.dataset(3) <= length(EPdataset.dataset)
            theCells=EPmain.view.theCells{EPmain.view.dataset(3)};
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theCells) > 1 && strcmp(EPdataset.dataset(EPmain.view.dataset(3)).dataType,'average') && ~ismember(EPmain.view.allTrials(3),[1 2 3 4])
                theCells{end+1}='-all-';
                theCells{end+1}='-erpimage-';
            end;
            EPmain.handles.view.cell(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theCells,'ForegroundColor',[0 .5 0],...
                'Value',EPmain.view.cell(3),'Position',[menuStart 500-menuY*12 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.cell(3),''value'');','if tempVar ~=0,EPmain.view.cell(3)=tempVar;end;','if isempty(tempVar),EPmain.view.cell(3)=tempVar;end;','EPmain.view.changeFlag(3)=1;','ep(''start'')']);
            
            theSubs=EPdataset.dataset(EPmain.view.dataset(3)).subNames;
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theSubs) > 1 && ~ismember(EPmain.view.allTrials(3),[3 4 5 6])
                theSubs{end+1}='-all-';
                theSubs{end+1}='-erpimage-';
            end;
            EPmain.handles.view.subject(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theSubs,'ForegroundColor',[0 .5 0],...
                'Value',EPmain.view.subject(3),'Position',[menuStart 500-menuY*13 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.subject(3),''value'');','if tempVar ~=0,EPmain.view.subject(3)=tempVar;end;','if isempty(tempVar),EPmain.view.subject(3)=tempVar;end;','EPmain.view.changeFlag(3)=1;','ep(''start'')']);
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(3)).trialNames)
                theTrials=num2cell(EPmain.view.trialList{3});
                if ~strcmp('TFT',EPmain.view.dataTransform) && ~ismember(EPmain.view.allTrials(3),[3 4 5 6])
                    theTrials{end+1}='-all-';
                    theTrials{end+1}='-erpimage-';
                end;
                EPmain.handles.view.trial(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor',[0 .5 0],...
                    'Value',EPmain.view.trial(3),'Position',[menuStart 500-menuY*14 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(3),''value'');','if tempVar ~=0,EPmain.view.trial(3)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(3)=tempVar;end;','EPmain.view.changeFlag(3)=1;','ep(''start'')']);
            elseif strcmp(EPdataset.dataset(EPmain.view.dataset(3)).dataType,'continuous')
                theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(3)).timeNames)/ceil(EPdataset.dataset(EPmain.view.dataset(3)).Fs))]'));
                for i=1:length(theEpochs)
                    theEpochs{i}=['Sec: ' theEpochs{i}];
                end;
                EPmain.handles.view.trial(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theEpochs,'ForegroundColor',[0 .5 0],...
                    'Value',EPmain.view.trial(3),'Position',[menuStart 500-menuY*14 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(3),''value'');','if tempVar ~=0,EPmain.view.trial(3)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(3)=tempVar;end;','EPmain.view.changeFlag(3)=1;','ep(''start'')']);
            elseif strcmp('VLT',EPmain.view.dataTransform)
                theTrials=cell(0);
                theTrials{end+1}='-none-';
                theTrials{end+1}='-GFP-';
                theTrials{end+1}='-noise-';
                theTrials{end+1}='-StDev-';
                theTrials{end+1}='-CI-';
                EPmain.handles.view.trial(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor',[0 .5 0],...
                    'Value',EPmain.view.trial(3),'Position',[menuStart 500-menuY*14 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(3),''value'');','if tempVar ~=0,EPmain.view.trial(3)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(3)=tempVar;end;','EPmain.view.changeFlag(3)=1;','ep(''start'')']);
            else
                EPmain.handles.view.trial(3) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Trials','FontSize',EPmain.fontsize,...
                    'ForegroundColor',[0 .5 0],'Position',[menuStart 500-menuY*14 menuX menuY]);
            end;
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(3)).facNames)
                theFacs=EPdataset.dataset(EPmain.view.dataset(3)).facNames;
                if ~strcmp('TFT',EPmain.view.dataTransform) && length(theFacs) > 1 && ~ismember(EPmain.view.allTrials(3),[1 2 5 6])
                    theFacs{end+1}='-all-';
                    theFacs{end+1}='-erpimage-';
                end;
                EPmain.handles.view.factor(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theFacs,'ForegroundColor',[0 .5 0],...
                    'Value',EPmain.view.factor(3),'Position',[menuStart 500-menuY*15 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.factor(3),''value'');','if tempVar ~=0,EPmain.view.factor(3)=tempVar;end;','if isempty(tempVar),EPmain.view.factor(3)=tempVar;end;','EPmain.view.changeFlag(3)=1;','ep(''start'')']);
            else
                EPmain.handles.view.factor(3) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Factors','FontSize',EPmain.fontsize,...
                    'ForegroundColor',[0 .5 0],'Position',[menuStart 500-menuY*15 menuX menuY]);
            end;
        end;
        
        %black
        theDatasets=[EPmain.view.theFileNames; {'none'}];
        EPmain.handles.view.dataset(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theDatasets,'ForegroundColor','black',...
            'Value',EPmain.view.dataset(4),'Position',[menuStart 500-menuY*16 menuX menuY],...
            'Callback', ['pause(.2);','global EPmain;','tempVar=get(EPmain.handles.view.dataset(4),''value'');','if tempVar ~=0,EPmain.view.dataset(4)=tempVar;end;','if isempty(tempVar),EPmain.view.dataset(4)=tempVar;end;','EPmain.view.changeDatasetFlag(4)=1;','ep(''start'')']);
        
        if EPmain.view.dataset(4) <= length(EPdataset.dataset)
            theCells=EPmain.view.theCells{EPmain.view.dataset(4)};
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theCells) > 1 && strcmp(EPdataset.dataset(EPmain.view.dataset(4)).dataType,'average') && ~ismember(EPmain.view.allTrials(4),[1 2 3 4])
                theCells{end+1}='-all-';
                theCells{end+1}='-erpimage-';
            end;
            EPmain.handles.view.cell(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theCells,'ForegroundColor','black',...
                'Value',EPmain.view.cell(4),'Position',[menuStart 500-menuY*17 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.cell(4),''value'');','if tempVar ~=0,EPmain.view.cell(4)=tempVar;end;','if isempty(tempVar),EPmain.view.cell(4)=tempVar;end;','EPmain.view.changeFlag(4)=1;','ep(''start'')']);
            
            theSubs=EPdataset.dataset(EPmain.view.dataset(4)).subNames;
            if ~strcmp('TFT',EPmain.view.dataTransform) && length(theSubs) > 1 && ~ismember(EPmain.view.allTrials(4),[3 4 5 6])
                theSubs{end+1}='-all-';
                theSubs{end+1}='-erpimage-';
            end;
            EPmain.handles.view.subject(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theSubs,'ForegroundColor','black',...
                'Value',EPmain.view.subject(4),'Position',[menuStart 500-menuY*18 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.subject(4),''value'');','if tempVar ~=0,EPmain.view.subject(4)=tempVar;end;','if isempty(tempVar),EPmain.view.subject(4)=tempVar;end;','EPmain.view.changeFlag(4)=1;','ep(''start'')']);
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(4)).trialNames)
                theTrials=num2cell(EPmain.view.trialList{4});
                if ~strcmp('TFT',EPmain.view.dataTransform) && ~ismember(EPmain.view.allTrials(4),[3 4 5 6])
                    theTrials{end+1}='-all-';
                    theTrials{end+1}='-erpimage-';
                end;
                EPmain.handles.view.trial(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor','black',...
                    'Value',EPmain.view.trial(4),'Position',[menuStart 500-menuY*19 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(4),''value'');','if tempVar ~=0,EPmain.view.trial(4)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(4)=tempVar;end;','EPmain.view.changeFlag(4)=1;','ep(''start'')']);
            elseif strcmp(EPdataset.dataset(EPmain.view.dataset(4)).dataType,'continuous')
                theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(4)).timeNames)/ceil(EPdataset.dataset(EPmain.view.dataset(4)).Fs))]'));
                for i=1:length(theEpochs)
                    theEpochs{i}=['Sec: ' theEpochs{i}];
                end;
                EPmain.handles.view.trial(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theEpochs,'ForegroundColor','black',...
                    'Value',EPmain.view.trial(4),'Position',[menuStart 500-menuY*19 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(4),''value'');','if tempVar ~=0,EPmain.view.trial(4)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(4)=tempVar;end;','EPmain.view.changeFlag(4)=1;','ep(''start'')']);
            elseif strcmp('VLT',EPmain.view.dataTransform)
                theTrials=cell(0);
                theTrials{end+1}='-none-';
                theTrials{end+1}='-GFP-';
                theTrials{end+1}='-noise-';
                theTrials{end+1}='-StDev-';
                theTrials{end+1}='-CI-';
                EPmain.handles.view.trial(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theTrials,'ForegroundColor','black',...
                    'Value',EPmain.view.trial(4),'Position',[menuStart 500-menuY*19 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(4),''value'');','if tempVar ~=0,EPmain.view.trial(4)=tempVar;end;','if isempty(tempVar),EPmain.view.trial(4)=tempVar;end;','EPmain.view.changeFlag(4)=1;','ep(''start'')']);
            else
                
                EPmain.handles.view.trial(4) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Trials','FontSize',EPmain.fontsize,...
                    'ForegroundColor','black','Position',[menuStart 500-menuY*19 menuX menuY]);
            end;
            
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(4)).facNames)
                theFacs=EPdataset.dataset(EPmain.view.dataset(4)).facNames;
                if ~strcmp('TFT',EPmain.view.dataTransform) && length(theFacs) > 1 && ~ismember(EPmain.view.allTrials(4),[1 2 5 6])
                    theFacs{end+1}='-all-';
                    theFacs{end+1}='-erpimage-';
                end;
                EPmain.handles.view.factor(4) = uicontrol('Style','popupmenu',...
                    'String',theFacs,'ForegroundColor','black','FontSize',EPmain.fontsize,...
                    'Value',EPmain.view.factor(4),'Position',[menuStart 500-menuY*20 menuX menuY],...
                    'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.factor(4),''value'');','if tempVar ~=0,EPmain.view.factor(4)=tempVar;end;','if isempty(tempVar),EPmain.view.factor(4)=tempVar;end;','EPmain.view.changeFlag(4)=1;','ep(''start'')']);
            else
                EPmain.handles.view.factor(4) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Factors','FontSize',EPmain.fontsize,...
                    'ForegroundColor','black','Position',[menuStart 500-menuY*20 menuX menuY]);
            end;
        end;
        
        if ~isempty(EPmain.view.eventList) || EPmain.view.RT
            eventList=EPmain.view.eventList;
            if EPmain.view.RT
                eventList{end+1}='-RT-';
            end;
            if length(eventList) > 1
                eventList{end+1}='-all-';
            end;
            eventList{end+1}='-none-';
            EPmain.handles.view.events = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',eventList,'ForegroundColor','magenta',...
                'Value',EPmain.view.events,'Position',[menuStart 500-menuY*21 menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.events,''value'');','if tempVar ~=0,EPmain.view.events=tempVar;end;','if isempty(tempVar),EPmain.view.events=tempVar;end;','ep(''start'')']);
        else
            EPmain.handles.view.events = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Events','FontSize',EPmain.fontsize,...
                'ForegroundColor','magenta','Position',[menuStart 500-menuY*21 menuX menuY]);
        end;
        
        uicontrol('Style','frame',...
            'Position',[5 36 155 62]);
        
        theMax=plotMVmax;
        theMin=plotMVmin;
        
        theTextColor.bottomVolt='black';
        if EPmain.view.edited.bottomVolt
            theMin=EPmain.view.manual.bottomVolt;
            theTextColor.bottomVolt='blue';
        end;
        theTextColor.topVolt='black';
        if EPmain.view.edited.topVolt
            theMax=EPmain.view.manual.topVolt;
            theTextColor.topVolt='blue';
        end;
        
        if all(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
                EPmain.handles.view.FFTunits = uicontrol('Style','text','HorizontalAlignment','left','String', 'Corr','FontSize',EPmain.fontsize,...
                    'ForegroundColor','black','Position',[10 77 50 20]);
        else
            if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
                %the cached min and max numbers are based on the greater of the min and max of the real and imaginary parts
                theMax=theMax/sqrt(EPmain.view.binHz); %convert to spectral density
                theMin=theMin/sqrt(EPmain.view.binHz); %convert to spectral density
                if EPmain.view.FFTunits > 2
                    theMax=theMax.^2; %convert amplitude to power
                    theMin=theMin.^2; %convert amplitude to power
                end;
                if (EPmain.view.FFTunits == 4)
                    theMax=log10(abs(theMax))*10; %convert to dB log scaling
                    theMin=log10(abs(theMin))*10; %convert to dB log scaling
                    if isinf(theMin) && (theMax > -100)
                        theMin=-100; %log10 of zero is -inf.  Replace with -2 to maintain useful range.
                        theTextColor.bottomVolt='blue';
                    end;
                end;
                EPmain.handles.view.FFTunits = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'Value',EPmain.view.FFTunits,'Position',[1 77 65 20],...
                    'String',{'cm','am','pw','dB'},...
                    'Callback', ['global EPmain;','EPmain.view.FFTunits=get(EPmain.handles.view.FFTunits,''value'');','ep(''start'')'],...
                    'TooltipString','Units for spectral data.');
                
            else
                EPmain.handles.view.FFTunits = uicontrol('Style','text','HorizontalAlignment','left','String', 'Voltage','FontSize',EPmain.fontsize,...
                    'ForegroundColor','black','Position',[10 77 50 20]);
            end;
        end;
        
        
        if EPmain.view.edited.startSamp
            startSamp=EPmain.view.manual.startSamp;
            theTextColor.startSamp='blue';
        else
            startSamp=EPmain.view.startSamp;
            theTextColor.startSamp='black';
        end;
        if EPmain.view.edited.endSamp
            endSamp=EPmain.view.manual.endSamp;
            theTextColor.endSamp='blue';
        else
            endSamp=EPmain.view.endSamp;
            theTextColor.endSamp='black';
        end;
        
        if EPmain.view.edited.startHz
            startHz=EPmain.view.manual.startHz;
            theTextColor.startHz='blue';
        else
            startHz=EPmain.view.startHz;
            theTextColor.startHz='black';
        end;
        if EPmain.view.edited.endHz
            endHz=EPmain.view.manual.endHz;
            theTextColor.endHz='blue';
        else
            endHz=EPmain.view.endHz;
            theTextColor.endHz='black';
        end;

        EPmain.handles.view.topVolt = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.4f', theMax),'FontSize',EPmain.fontsize,...
            'Position',[10 58 40 20],'ForegroundColor',theTextColor.topVolt,'Callback',@checkViewSettings);
        EPmain.handles.view.bottomVolt = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.4f', theMin),'FontSize',EPmain.fontsize,...
            'Position',[10 38 40 20],'ForegroundColor',theTextColor.bottomVolt,'Callback',@checkViewSettings);
        
        if EPmain.view.flexMode
            theUnit='%';
        else
            theUnit='Ms';
        end;
        h = uicontrol('Style','text','HorizontalAlignment','left','String', theUnit,'FontSize',EPmain.fontsize,...
            'ForegroundColor','black','Position',[55 77 50 20]);
        EPmain.handles.view.startSamp = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(startSamp)),'FontSize',EPmain.fontsize,...
            'Position',[55 58 35 20],'ForegroundColor',theTextColor.startSamp,'Callback',@checkViewSettings);
        EPmain.handles.view.endSamp = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(endSamp)),'FontSize',EPmain.fontsize,...
            'Position',[55 38 35 20],'ForegroundColor',theTextColor.endSamp,'Callback',@checkViewSettings);
        if strcmp('FFT',EPmain.view.dataTransform)
            set(EPmain.handles.view.startSamp ,'enable','off');
            set(EPmain.handles.view.endSamp ,'enable','off');
        end;
        
        h = uicontrol('Style','text','HorizontalAlignment','left','String', 'Hz','FontSize',EPmain.fontsize,...
            'ForegroundColor','black','Position',[90 77 50 20]);
        EPmain.handles.view.startHz = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(startHz)),'FontSize',EPmain.fontsize,...
            'Position',[90 58 35 20],'ForegroundColor',theTextColor.startHz,'Callback',@checkViewSettings);
        EPmain.handles.view.endHz = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(endHz)),'FontSize',EPmain.fontsize,...
            'Position',[90 38 35 20],'ForegroundColor',theTextColor.endHz,'Callback',@checkViewSettings);
        if strcmp('VLT',EPmain.view.dataTransform)
            set(EPmain.handles.view.startHz ,'enable','off');
            set(EPmain.handles.view.endHz ,'enable','off');
        end;
        
        h = uicontrol('Style','text','HorizontalAlignment','left','String', 'Mark','FontSize',EPmain.fontsize,...
            'ForegroundColor','black','Position',[125 77 30 20]);
        
        EPmain.handles.view.marker1 = uicontrol('Style','edit','HorizontalAlignment','left','String', marker1,'FontSize',EPmain.fontsize,...
            'Position',[125 58 30 20],...
            'Callback',@checkViewSettings);
        
        EPmain.handles.view.marker2 = uicontrol('Style','edit','HorizontalAlignment','left','String', marker2,'FontSize',EPmain.fontsize,...
            'Position',[125 38 30 20],...
            'Callback',@checkViewSettings);
        
%          h = uicontrol('Style','text','HorizontalAlignment','left','String', 'evt','FontSize',EPmain.fontsize,...
%             'ForegroundColor','black','Position',[160 57 25 20]);
%         EPmain.handles.view.evt = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
%             'Value',EPmain.view.evt,'Position',[180 80 50 20],...
%             'Callback', ['global EPmain;','EPmain.view.evt=get(EPmain.handles.view.evt,''value'');','ep(''start'')'],...
%             'TooltipString','Include event lines in graphs.');
        
        if ~strcmp('VLT',EPmain.view.dataTransform)
            set(EPmain.handles.view.evt,'enable','off');
        end;

        EPmain.handles.view.waves= uicontrol('Style', 'pushbutton', 'String', 'Waves','FontSize',EPmain.fontsize,...
            'Position', [2 0 50 35], 'Callback', 'ep(''viewWaves'')');
        
        EPmain.handles.view.topos= uicontrol('Style', 'pushbutton', 'String', 'Topos','FontSize',EPmain.fontsize,...
            'Position', [52 0 50 35], 'Callback', 'ep(''viewTopos'')');
        
        EPmain.handles.view.scan = uicontrol('Style', 'pushbutton', 'String', 'Scan','FontSize',EPmain.fontsize,...
            'Position', [102 0 50 35], 'Callback', 'ep(''viewEdit'')');
        
        if EPmain.view.dataset(1) > length(EPdataset.dataset)
            set(EPmain.handles.view.scan ,'enable','off');
        end
        
        if strcmp('FFT',EPmain.view.dataTransform) || any(EPmain.view.allTrials == 2) || (EPmain.view.allTrials(1) == 1)
            set(EPmain.handles.view.scan ,'enable','off'); %can't scan with erpimages or FFT or if main is 'all'
        end;
        
        if any(ismember(EPmain.view.allTrials,[2 4 6]))
            set(EPmain.handles.view.topos ,'enable','off'); %can't view topos with erpimage options
        end;
        
        for iColor=1:4
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                if (strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average') && (EPmain.view.trial(iColor) > 1))
                    %can't view topos or scan with any band settings
                    set(EPmain.handles.view.topos ,'enable','off');
                    set(EPmain.handles.view.scan ,'enable','off');
                end;
            end;
        end;
        
        EPmain.handles.view.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [152 0 50 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','EPmain.view=[];','ep(''start'');']);
        
        if theMin >= theMax
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
        end;
        
        if startSamp > endSamp
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
        end;
        
        if (startHz > endHz) || (strcmp('FFT',EPmain.view.dataTransform) && (startHz == endHz))
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
        end;
        
        if (startSamp == endSamp) && (startHz == endHz)
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
        end;

        if ~isempty(marker1)
            if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
                if (str2num(marker1) < startSamp) || (str2num(marker1) > endSamp)
                    set(EPmain.handles.view.waves,'enable','off');
                    set(EPmain.handles.view.topos,'enable','off');
                end;
            else
                if (str2num(marker1) < startHz) || (str2num(marker1) > endHz)
                    set(EPmain.handles.view.waves,'enable','off');
                    set(EPmain.handles.view.topos,'enable','off');
                end;
            end;
        end;
        
        if ~isempty(marker2)
            if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
                if (str2num(marker2) < startSamp) || (str2num(marker2) > endSamp)
                    set(EPmain.handles.view.waves,'enable','off');
                    set(EPmain.handles.view.topos,'enable','off');
                end;
            else
                if (str2num(marker2) < startHz) || (str2num(marker2) > endHz)
                    set(EPmain.handles.view.waves,'enable','off');
                    set(EPmain.handles.view.topos,'enable','off');
                end;
            end;
        end;
        
        if EPmain.view.flexMode
            if (startSamp < 0) || (startSamp > 100)
                set(EPmain.handles.view.waves,'enable','off');
                set(EPmain.handles.view.topos,'enable','off');
            end;
        end;
        
    case 'viewWaves'
        set(EPmain.handles.view.waves,'enable','off');
        set(EPmain.handles.view.topos,'enable','off');
        set(EPmain.handles.view.scan,'enable','off');
        set(EPmain.handles.view.done,'enable','off');
        set(EPmain.handles.view.FFTunits,'enable','off');
        
        for color=1:4
            EPmain.view.dataset(color)=get(EPmain.handles.view.dataset(color),'value');
            if EPmain.view.dataset(color) <= length(EPdataset.dataset)
                EPmain.view.cell(color)=get(EPmain.handles.view.cell(color),'value');
                EPmain.view.subject(color)=get(EPmain.handles.view.subject(color),'value');
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(color)).trialNames)
                    EPmain.view.trial(color)=get(EPmain.handles.view.trial(color),'value');
                end;
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(color)).facNames)
                    EPmain.view.factor(color)=get(EPmain.handles.view.factor(color),'value');
                end;
            end;
        end;
        
        EPwaves=[];
        EPwaves.mode=[];
        EPwaves.direction=EPmain.preferences.view.positive;
        
        eventLines=ep_collateEventLines('EPwaves');

        if EPmain.view.edited.bottomVolt
            theMin=str2num(get(EPmain.handles.view.bottomVolt,'string'));
        else
            theMin=[];
        end;
        if EPmain.view.edited.topVolt
            theMax=str2num(get(EPmain.handles.view.topVolt,'string'));
        else
            theMax=[];
        end;
        ep_showWaves(theMin,theMax,...
            str2num(get(EPmain.handles.view.startSamp,'string')),str2num(get(EPmain.handles.view.endSamp,'string')),...
            str2num(get(EPmain.handles.view.startHz,'string')),str2num(get(EPmain.handles.view.endHz,'string')),...
            str2num(get(EPmain.handles.view.marker1,'string')),str2num(get(EPmain.handles.view.marker2,'string')),EPmain.view.FFTunits,eventLines,1);
        
        try
            set(EPmain.handles.view.waves,'enable','on');
            if ~any(EPmain.view.allTrials)
                set(EPmain.handles.view.topos,'enable','on');
                set(EPmain.handles.view.scan,'enable','on');
            end;
            set(EPmain.handles.view.done,'enable','on');
            if ~strcmp('VLT',EPmain.view.dataTransform)
                set(EPmain.handles.view.FFTunits,'enable','on');
            end;
            if strcmp('VLT',EPmain.view.dataTransform)
                set(EPmain.handles.view.evt,'enable','on');
            end;
        catch
            ep('start')
        end;
        
    case 'viewTopos'
        set(EPmain.handles.view.waves,'enable','off');
        set(EPmain.handles.view.topos,'enable','off');
        set(EPmain.handles.view.scan,'enable','off');
        set(EPmain.handles.view.done,'enable','off');
        set(EPmain.handles.view.marker1,'enable','off');
        set(EPmain.handles.view.marker2,'enable','off');
        set(EPmain.handles.view.events ,'enable','off');
        set(EPmain.handles.view.startSamp,'enable','off');
        set(EPmain.handles.view.endSamp,'enable','off');
        set(EPmain.handles.view.topVolt,'enable','off');
        set(EPmain.handles.view.bottomVolt,'enable','off');
        set(EPmain.handles.view.startHz,'enable','off');
        set(EPmain.handles.view.endHz,'enable','off');
        set(EPmain.handles.view.FFTunits,'enable','off');
        
        for iColor=1:4
            set(EPmain.handles.view.dataset(iColor),'enable','off');
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                set(EPmain.handles.view.cell(iColor),'enable','off');
                set(EPmain.handles.view.subject(iColor),'enable','off');
                set(EPmain.handles.view.trial(iColor),'enable','off');
                set(EPmain.handles.view.factor(iColor),'enable','off');
            end;
        end;
        
        for iColor=1:4
            EPmain.view.dataset(iColor)=get(EPmain.handles.view.dataset(iColor),'value');
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                EPmain.view.cell(iColor)=get(EPmain.handles.view.cell(iColor),'value');
                EPmain.view.subject(iColor)=get(EPmain.handles.view.subject(iColor),'value');
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                    EPmain.view.trial(iColor)=get(EPmain.handles.view.trial(iColor),'value');
                end;
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                    EPmain.view.factor(iColor)=get(EPmain.handles.view.factor(iColor),'value');
                end;
            end;
        end;
        
        eventLines=ep_collateEventLines('EPtopos');
        
        EPtopos=[];
        EPtopos.page=0;
        EPtopos.direction=EPmain.preferences.view.positive;
        if EPmain.view.edited.bottomVolt
            theMin=str2num(get(EPmain.handles.view.bottomVolt,'string'));
        else
            theMin=[];
        end;
        if EPmain.view.edited.topVolt
            theMax=str2num(get(EPmain.handles.view.topVolt,'string'));
        else
            theMax=[];
        end;
        
        ep_showTopos(str2num(get(EPmain.handles.view.startSamp,'string')),str2num(get(EPmain.handles.view.endSamp,'string')),...
            str2num(get(EPmain.handles.view.startHz,'string')),str2num(get(EPmain.handles.view.endHz,'string')),...
            str2num(get(EPmain.handles.view.marker1,'string')),str2num(get(EPmain.handles.view.marker2,'string')),EPmain.view.FFTunits,theMin,theMax,eventLines);

    case 'viewEdit'
        set(EPmain.handles.view.waves,'enable','off');
        set(EPmain.handles.view.topos,'enable','off');
        set(EPmain.handles.view.scan,'enable','off');
        set(EPmain.handles.view.done,'enable','off');
        set(EPmain.handles.view.FFTunits,'enable','off');
        
        set(EPmain.handles.view.startSamp ,'enable','off');
        set(EPmain.handles.view.endSamp ,'enable','off');
        set(EPmain.handles.view.startHz ,'enable','off');
        set(EPmain.handles.view.endHz ,'enable','off');
        set(EPmain.handles.view.topVolt ,'enable','off');
        set(EPmain.handles.view.bottomVolt ,'enable','off');
        set(EPmain.handles.view.marker1 ,'enable','off');
        set(EPmain.handles.view.marker2 ,'enable','off');
        set(EPmain.handles.view.events ,'enable','off');
        
        for iColor=1:4
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                set(EPmain.handles.view.dataset(iColor),'enable','off');
                set(EPmain.handles.view.subject(iColor),'enable','off');
                set(EPmain.handles.view.cell(iColor),'enable','off');
                set(EPmain.handles.view.trial(iColor),'enable','off');
                set(EPmain.handles.view.factor(iColor),'enable','off');
            end;
        end;
        
        for iColor=1:4
            EPmain.view.dataset(iColor)=get(EPmain.handles.view.dataset(iColor),'value');
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                EPmain.view.cell(iColor)=get(EPmain.handles.view.cell(iColor),'value');
                EPmain.view.subject(iColor)=get(EPmain.handles.view.subject(iColor),'value');
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                    EPmain.view.trial(iColor)=get(EPmain.handles.view.trial(iColor),'value');
                end;
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                    EPmain.view.factor(iColor)=get(EPmain.handles.view.factor(iColor),'value');
                end;
            end;
        end;
        
        eventLines=ep_collateEventLines('EPwaves');
        
        EPwaves.mode=[];
        EPwaves.direction=EPmain.preferences.view.positive;
        if strcmp(EPdataset.dataset(EPmain.view.dataset(1)).dataType,'continuous')
            ep_scanCont;
        else
            err=ep_showWaves(str2num(get(EPmain.handles.view.bottomVolt,'string')),str2num(get(EPmain.handles.view.topVolt,'string')),...
                str2num(get(EPmain.handles.view.startSamp,'string')),str2num(get(EPmain.handles.view.endSamp,'string')),...
                str2num(get(EPmain.handles.view.startHz,'string')),str2num(get(EPmain.handles.view.endHz,'string')),...
                str2num(get(EPmain.handles.view.marker1,'string')),str2num(get(EPmain.handles.view.marker2,'string')),EPmain.view.FFTunits,eventLines,0);
            
            if ~err %if the wave window did not error out
                EPmanualEdit=[];
                ep_manualEdit
            else
                try
                    set(EPmain.handles.view.waves,'enable','on');
                    set(EPmain.handles.view.topos,'enable','on');
                    set(EPmain.handles.view.scan,'enable','on');
                    set(EPmain.handles.view.done,'enable','on');
                    if ~strcmp('VLT',EPmain.view.dataTransform)
                        set(EPmain.handles.view.FFTunits,'enable','on');
                    end;
                    for iColor=1:4
                        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                            set(EPmain.handles.view.dataset(iColor),'enable','on');
                            set(EPmain.handles.view.subject(iColor),'enable','on');
                            set(EPmain.handles.view.cell(iColor),'enable','on');
                            set(EPmain.handles.view.trial(iColor),'enable','on');
                            set(EPmain.handles.view.factor(iColor),'enable','on');
                        end;
                    end;
                    set(EPmain.handles.view.startSamp ,'enable','on');
                    set(EPmain.handles.view.endSamp ,'enable','on');
                    set(EPmain.handles.view.startHz ,'enable','on');
                    set(EPmain.handles.view.endHz ,'enable','on');
                    set(EPmain.handles.view.topVolt ,'enable','on');
                    set(EPmain.handles.view.bottomVolt ,'enable','on');
                    set(EPmain.handles.view.marker1 ,'enable','on');
                    set(EPmain.handles.view.marker2 ,'enable','on');
                    set(EPmain.handles.view.events ,'enable','on');
                catch
                    ep('start')
                end;
            end;
        end;
        
    case 'startSampleTest'
        
        set(EPmain.handles.hMainWindow,'Name', 'sampleTest');
        noSampTest=0;
        
        if ~isfield(EPmain.sampleTest,'dataset') %just entering sampleTest function
            EPmain.sampleTest.datasetList=[];
            EPmain.sampleTest.PCAlist=[];
            EPmain.sampleTest.AVElist=[];
            for iFile=1:length(EPdataset.dataset)
                if ~isempty(EPdataset.dataset(iFile).timeNames) && isempty(EPdataset.dataset(iFile).facVecT) %cannot perform analysis on FFT data or temporal PCA data
                    switch EPdataset.dataset(iFile).dataType
                        case 'continuous'
                            EPmain.sampleTest.datasetList(end+1)=iFile;
                        case 'single_trial'
                            if length(unique(EPdataset.dataset(iFile).cellNames(find(strcmp(EPdataset.dataset(iFile).cellTypes,'SGL'))))) > 1
                                EPmain.sampleTest.datasetList(end+1)=iFile;
                            end;
                        case 'average'
                            if length(EPdataset.dataset(iFile).subNames(find(strcmp(EPdataset.dataset(iFile).subTypes,'AVG')))) > 1
                                EPmain.sampleTest.datasetList(end+1)=iFile;
                            end;
                    end;
                end;
            end;
            EPmain.sampleTest.dataset=length(EPmain.sampleTest.datasetList);
            EPmain.sampleTest.cell1=1;
            EPmain.sampleTest.cell2=2;
            if EPmain.sampleTest.dataset
                EPmain.sampleTest.cellNameList=unique(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).cellNames(find(strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).cellTypes,'SGL'))));
                EPmain.sampleTest.datasetNameList={EPdataset.dataset(EPmain.sampleTest.datasetList).dataName};
            else
                EPmain.sampleTest.cellNameList={'none'};
                EPmain.sampleTest.datasetNameList={'none'};
            end;
            EPmain.sampleTest.method=1;
            EPmain.sampleTest.contMethod=1;
            EPmain.sampleTest.test=1;
            EPmain.sampleTest.PCA=1;
            EPmain.sampleTest.AVE=1;
            EPmain.sampleTest.alpha=.05;
            EPmain.sampleTest.contiguous=4;
            EPmain.sampleTest.thresh=.5;
            EPmain.sampleTest.channel=1;
            EPmain.sampleTest.filterPass=1;
            EPmain.sampleTest.filterType=2;
            EPmain.sampleTest.filter1=[];
            EPmain.sampleTest.filter2=[];
            EPmain.sampleTest.filterOrder=6;
            EPmain.sampleTest.scale=300;
            EPmain.sampleTest.factor=1;
            EPmain.sampleTest.subject=1;
            EPmain.sampleTest.cell=1;
            EPmain.sampleTest.gridSize=67;
            EPmain.sampleTest.AVEdata=cell(0);
            EPmain.sampleTest.sampStart=[];
            EPmain.sampleTest.sampEnd=[];
            EPmain.sampleTest.changeFlag=1;
            EPmain.sampleTest.freqFlag=0;
            EPmain.sampleTest.freqBin=1;
            EPmain.sampleTest.freqScale=1;
            
            if EPmain.sampleTest.dataset
                EEGchans=find(strcmp('EEG',EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).chanTypes));
                EPmain.sampleTest.eloc=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).eloc(EEGchans);
                EPmain.sampleTest.dataType=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).dataType;
                EPmain.sampleTest.freqFlag=~isempty(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames);
                
                maxRad=0.5;
                [y,x] = pol2cart(([EPmain.sampleTest.eloc.theta]/360)*2*pi,[EPmain.sampleTest.eloc.radius]);  % transform electrode locations from polar to cartesian coordinates
                y=-y; %flip y-coordinate so that nose is upwards.
                plotrad = min(1.0,max([EPmain.sampleTest.eloc.radius])*1.02);            % default: just outside the outermost electrode location
                plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
                x = x*(maxRad/plotrad);
                y = y*(maxRad/plotrad);
                
                xmin = min(-maxRad,min(x));
                xmax = max(maxRad,max(x));
                ymin = min(-maxRad,min(y));
                ymax = max(maxRad,max(y));
                
                EPmain.sampleTest.x=round(((x/(xmax-xmin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));
                EPmain.sampleTest.y=round(((y/(ymax-ymin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));
                
                for iFile=1:length(EPdataset.dataset)
                    if ~isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc)
                        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
                        if length(EPmain.sampleTest.eloc) == length(newEloc)
                            if ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                                if (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames))
                                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames])
                                        EPmain.sampleTest.PCAlist(end+1)=iFile;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
                EPmain.sampleTest.PCAnameList={EPdataset.dataset(EPmain.sampleTest.PCAlist).dataName};
                
                for iFile=1:length(EPdataset.dataset)
                    if isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc) && strcmp(EPdataset.dataset(iFile).dataType,'average')
                        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
                        if length(EPmain.sampleTest.eloc) == length(newEloc)
                            if ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                                if (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames))
                                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames])
                                        EPmain.sampleTest.AVElist(end+1)=iFile;
                                        EPdata=ep_loadEPdataset(EPdataset,iFile);
                                        EPmain.sampleTest.AVEdata{end+1}=EPdata.data;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
                EPmain.sampleTest.AVEnameList={EPdataset.dataset(EPmain.sampleTest.AVElist).dataName};
            end;
        else
            EEGchans=[];
            EPmain.sampleTest.eloc=[];
        end;
        
        if EPmain.sampleTest.dataset
            EPmain.handles.sampleTest.dataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',EPmain.sampleTest.datasetNameList,...
                'Value',EPmain.sampleTest.dataset,'Position',[5 480 160 20],...
                'Callback', @changeSampleTestDataset);
            
            if ~strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).dataType,'continuous')
                EPmain.handles.sampleTest.cell1 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.sampleTest.cellNameList,...
                    'Value',EPmain.sampleTest.cell1,'Position',[5 460 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.cell1,''Value'');','if tempVar ~=0,EPmain.sampleTest.cell1=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.cell1=tempVar;end;','ep(''start'');']);
                
                EPmain.handles.sampleTest.cell2 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.sampleTest.cellNameList,...
                    'Value',EPmain.sampleTest.cell2,'Position',[5 440 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.cell2,''Value'');','if tempVar ~=0,EPmain.sampleTest.cell2=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.cell2=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Method','HorizontalAlignment','left',...
                    'Position',[10 420 150 20]);
                
                EPmain.handles.sampleTest.method = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'sample','CWT','Template Woody','PCA Woody'},...
                    'Value',EPmain.sampleTest.method,'Position',[5 400 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.method,''Value'');','if tempVar ~=0,EPmain.sampleTest.method=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.method=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Test','HorizontalAlignment','left',...
                    'Position',[10 380 150 20]);
                
                EPmain.handles.sampleTest.test = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'t-test','jack-knife'},...
                    'Value',EPmain.sampleTest.test,'Position',[5 360 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.test,''Value'');','if tempVar ~=0,EPmain.sampleTest.test=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.test=tempVar;end;','ep(''start'');']);
                
                switch EPmain.sampleTest.test
                    case 1 %t=test
                        
                        EPmain.handles.sampleTest.alphaLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Alpha','FontSize',EPmain.fontsize,...
                            'Position',[10 340 70 20]);
                        
                        EPmain.handles.sampleTest.alpha = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.alpha),'FontSize',EPmain.fontsize,...
                            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.alpha,''String''));','if tempVar ~=0,EPmain.sampleTest.alpha=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.alpha=tempVar;end;','ep(''start'');'],...
                            'Position',[10 320 70 20],'TooltipString','Threshold for sample by sample significance testing.');
                        
                        EPmain.handles.sampleTest.contiguousLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Contiguous','FontSize',EPmain.fontsize,...
                            'Position',[80 340 70 20]);
                        
                        EPmain.handles.sampleTest.contiguous = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.sampleTest.contiguous,'FontSize',EPmain.fontsize,...
                            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.contiguous,''String'');','if tempVar ~=0,EPmain.sampleTest.contiguous=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.contiguous=tempVar;end;','ep(''start'');'],...
                            'Position',[80 320 70 20],'TooltipString','How many contiguous significant samples are required to consider to be significant.');
                        
                    case 2 %jack-knife
                        
                        EPmain.handles.sampleTest.threshLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Threshold','FontSize',EPmain.fontsize,...
                            'Position',[10 340 70 20]);
                        
                        EPmain.handles.sampleTest.thresh = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.thresh),'FontSize',EPmain.fontsize,...
                            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.thresh,''String''));','if tempVar ~=0,EPmain.sampleTest.thresh=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.thresh=tempVar;end;','ep(''start'');'],...
                            'Position',[10 320 70 20],'TooltipString','Threshold for jack-knife estimation of onset latency.');
                        
                        EPmain.handles.sampleTest.chanLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Channel','FontSize',EPmain.fontsize,...
                            'Position',[80 340 70 20]);
                        
                        EPmain.handles.sampleTest.channel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                            'String',EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).chanNames(find(strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).chanTypes,'EEG'))),...
                            'Value',EPmain.sampleTest.channel,'Position',[80 320 100 20],...
                            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.channel,''Value'');','if tempVar ~=0,EPmain.sampleTest.channel=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.channel=tempVar;end;','ep(''start'');']);
                        
                end;
                
                EPmain.handles.sampleTest.filterPass= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'Low Pass','High Pass','Band Pass','Band Stop','Notch'},...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.filterPass,''Value'');','if tempVar ~=0,EPmain.sampleTest.filterPass=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filterPass=tempVar;end;','ep(''start'');'],...
                    'Value',EPmain.sampleTest.filterPass,'Position',[5 300 110 20]);
                
                EPmain.handles.sampleTest.filter1 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.filter1),'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.filter1,''String''));','if tempVar ~=0,EPmain.sampleTest.filter1=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filter1=tempVar;end;','ep(''start'');'],...
                    'Position',[110 300 40 20],'TooltipString','Lower frequency limit.');
                
                EPmain.handles.sampleTest.filter2 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.filter2),'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.filter2,''String''));','if tempVar ~=0,EPmain.sampleTest.filter2=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filter2=tempVar;end;','ep(''start'');'],...
                    'Position',[155 300 40 20],'TooltipString','Upper frequency limit.');
                
                EPmain.handles.sampleTest.filterType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'One-Pass Butterworth','Two-Pass Butterworth','One-Pass FIR','Two-Pass FIR','One-Pass FIRLS','Two-Pass FIRLS'},...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.filterType,''Value'');','if tempVar ~=0,EPmain.sampleTest.filterType=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filterType=tempVar;end;','ep(''start'');'],...
                    'Value',EPmain.sampleTest.filterType,'Position',[5 280 160 20]);
                
                EPmain.handles.sampleTest.filterOrder = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.filterOrder),'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.filterOrder,''String''));','if tempVar ~=0,EPmain.sampleTest.filterOrder=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filterOrder=tempVar;end;','ep(''start'');'],...
                    'Position',[160 280 30 20],'TooltipString','Order of the filter.');
                
                if ~isempty(EPdataset.dataset(EPmain.sampleTest.dataset).facVecT) %cannot filter temporal PCA data
                    set(EPmain.handles.sampleTest.filterPass,'enable','off');
                    set(EPmain.handles.sampleTest.filter1,'enable','off');
                    set(EPmain.handles.sampleTest.filter2,'enable','off');
                    set(EPmain.handles.sampleTest.filterType,'enable','off');
                    set(EPmain.handles.sampleTest.filterOrder,'enable','off');
                end;
                
                EPmain.handles.sampleTest.scaleLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Scale','FontSize',EPmain.fontsize,...
                    'Position',[10 260 50 20]);
                
                EPmain.handles.sampleTest.scale = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.sampleTest.scale,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.scale,''String'');','if tempVar ~=0,EPmain.sampleTest.scale=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.scale=tempVar;end;','ep(''start'');'],...
                    'Position',[70 260 50 20],'TooltipString','Width of the Mexican hat wavelet template.');
                
                if EPmain.sampleTest.method ~= 2
                    set(EPmain.handles.sampleTest.scaleLabel,'enable','off');
                    set(EPmain.handles.sampleTest.scale,'enable','off');
                end;
            else
                EPmain.handles.sampleTest.contMethod = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'Template Woody','PCA Woody'},...
                    'Value',EPmain.sampleTest.contMethod,'Position',[5 400 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.contMethod,''Value'');','if tempVar ~=0,EPmain.sampleTest.contMethod=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.contMethod=tempVar;end;','ep(''start'');']);              
            end;
            
            %AVE template
            if (strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).dataType,'continuous') && (EPmain.sampleTest.contMethod == 1)) || (~strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).dataType,'continuous') && (EPmain.sampleTest.method == 3))
                EPmain.handles.sampleTest.AVElabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Average Template','HorizontalAlignment','left',...
                    'Position',[10 240 150 20]);
                
                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.handles.sampleTest.AVE = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.AVEnameList,...
                        'Value',EPmain.sampleTest.AVE,'Position',[5 220 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.AVE,''Value'');','if tempVar ~=0,EPmain.sampleTest.AVE=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.AVE=tempVar;end;','EPmain.sampleTest.subject=1;','EPmain.sampleTest.cell=1;','ep(''start'');']);
                else
                    EPmain.handles.sampleTest.AVE = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 220 160 20]);
                end;
                
                EPmain.handles.sampleTest.subLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Subject','HorizontalAlignment','left',...
                    'Position',[10 200 150 20]);
                
                EPmain.handles.sampleTest.cellLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Cell','HorizontalAlignment','left',...
                    'Position',[10 180 150 20]);

                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.handles.sampleTest.subject = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).subNames,...
                        'Value',EPmain.sampleTest.subject,'Position',[50 200 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.subject,''Value'');','if tempVar ~=0,EPmain.sampleTest.subject=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.subject=tempVar;end;','ep(''start'');']);
                    
                    EPmain.handles.sampleTest.cell = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).cellNames,...
                        'Value',EPmain.sampleTest.cell,'Position',[50 180 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.cell,''Value'');','if tempVar ~=0,EPmain.sampleTest.cell=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.cell=tempVar;end;','ep(''start'');']);
                else
                    noSampTest=1;
                    EPmain.handles.sampleTest.subject = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 180 160 20]);
                    
                    EPmain.handles.sampleTest.cell = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 140 160 20]);
                end;
                %PCA Woody
            elseif (strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).dataType,'continuous') && (EPmain.sampleTest.contMethod == 2)) || (~strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).dataType,'continuous') && (EPmain.sampleTest.method == 4))
                EPmain.handles.sampleTest.PCAlabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','PCA','HorizontalAlignment','left',...
                    'Position',[10 240 150 20]);
                
                if ~isempty(EPmain.sampleTest.PCAlist)
                    EPmain.handles.sampleTest.PCA = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.PCAnameList,...
                        'Value',EPmain.sampleTest.PCA,'Position',[5 220 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.PCA,''Value'');','if tempVar ~=0,EPmain.sampleTest.PCA=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.PCA=tempVar;end;','EPmain.sampleTest.factor=1;','ep(''start'');']);
                else
                    EPmain.handles.sampleTest.PCA = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 220 160 20]);
                end;
                
                EPmain.handles.sampleTest.factorLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Factor','HorizontalAlignment','left',...
                    'Position',[10 200 150 20]);
                
                if ~isempty(EPmain.sampleTest.PCAlist)
                    EPmain.handles.sampleTest.factor = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facNames,...
                        'Value',EPmain.sampleTest.factor,'Position',[50 200 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.factor,''Value'');','if tempVar ~=0,EPmain.sampleTest.factor=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.factor=tempVar;end;','ep(''start'');']);
                else
                    noSampTest=1;
                    EPmain.handles.sampleTest.factor = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[50 200 160 20]);
                end;
            end;
            
            if EPmain.sampleTest.freqFlag
                EPmain.handles.sampleTest.freqLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Hz','HorizontalAlignment','left',...
                    'Position',[10 160 150 20]);
                EPmain.handles.sampleTest.freqScaleLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Hz Scale','HorizontalAlignment','left',...
                    'Position',[10 140 150 20]);
                theFreqs=num2cell(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames);
                theFreqs{end+1}='-all-';
                EPmain.handles.sampleTest.freqBin = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theFreqs,...
                    'Value',EPmain.sampleTest.freqBin,'Position',[50 160 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.freqBin,''Value'');','if tempVar ~=0,EPmain.sampleTest.freqBin=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.freqBin=tempVar;end;','ep(''start'');']);
                EPmain.handles.sampleTest.freqScale = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'cm-real','cm-imag','am','pw','dB'},...
                    'Value',EPmain.sampleTest.freqScale,'Position',[50 140 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.freqScale,''Value'');','if tempVar ~=0,EPmain.sampleTest.freqScale=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.freqScale=tempVar;end;','ep(''start'');']);
            end;
            
            EEGchans=find(strcmp('EEG',EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).chanTypes));
            EEGchans=EEGchans(~cellfun(@isempty,{EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).eloc(EEGchans).radius}));
            
            if strcmp(EPmain.sampleTest.dataType,'continuous')
                switch EPmain.sampleTest.contMethod
                    case 1 %AVE template
                        if ~isempty(EPmain.sampleTest.AVElist)
                            templateData=EPmain.sampleTest.AVEdata{EPmain.sampleTest.AVE}(EEGchans,:,EPmain.sampleTest.cell,EPmain.sampleTest.subject,:,EPmain.sampleTest.freqBin);
                            if EPmain.sampleTest.freqFlag
                                if (EPmain.sampleTest.freqScale == 1)
                                    templateData=real(templateData);
                                end;
                                if (EPmain.sampleTest.freqScale == 2)
                                    templateData=imag(templateData);
                                end;
                                if (EPmain.sampleTest.freqScale > 2)
                                    templateData(EEGchans,:,:,:)=abs(templateData(EEGchans,:,:,:)); %convert complex number to real number
                                end;
                                templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:)/mean(diff(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)); %convert to spectral density
                                if EPmain.sampleTest.freqScale > 3
                                    templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:).^2; %convert amplitude to power
                                end;
                                if (EPmain.sampleTest.freqScale == 5)
                                    if ~all(templateData(EEGchans,:,:,:) >=0)
                                        disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
                                    end;
                                    templateData(EEGchans,:,:,:)=log10(abs(templateData(EEGchans,:,:,:)))*10; %convert to dB log scaling
                                    tempVar=templateData(EEGchans,:,:,:);
                                    tempVar(isinf(tempVar))=-flintmax;
                                    templateData(EEGchans,:,:,:)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                                end;
                            end;
                            if ~isempty(EPmain.sampleTest.sampEnd)
                                templateData2=templateData(:,EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd);
                            else
                                templateData2=templateData;
                            end;
                            [A maxChan]=max(max(templateData2'));
                            [A maxPoint]=max(max(templateData2));
                            if ~isempty(EPmain.sampleTest.sampStart)
                                maxPoint=maxPoint+EPmain.sampleTest.sampStart-1;
                            end;
                            EPmain.sampleTest.templateWaveform=templateData(maxChan,:)';
                            EPmain.sampleTest.templateTopo=templateData(:,maxPoint);
                            tMs=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).timeNames;
                        else
                            templateData=zeros(length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames),length(EEGchans));
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames),1);
                            EPmain.sampleTest.templateTopo=zeros(1,length(EEGchans));
                            tMs=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames;
                        end;
                    case 2 %PCA
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT)
                            EPmain.sampleTest.templateWaveform=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT(:,EPmain.sampleTest.factor);
                            tMs=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).timeNames;
                        else
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames),1);
                            tMs=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames;
                        end;
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS)
                            templateTopo=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS(EEGchans,EPmain.sampleTest.factor);
                            if EPmain.sampleTest.test==2 %not just one channel
                                set(EPmain.handles.sampleTest.chanLabel,'enable','off');
                                set(EPmain.handles.sampleTest.channel,'enable','off');
                            end;
                        else
                            templateTopo=zeros(length(EEGchans));
                        end;
                    otherwise
                        disp('Oops: programmer error');
                        return;
                end;
                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.sampleTest.Fs=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).Fs;
                    EPmain.sampleTest.baseline=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).baseline;
                end;
            else
                switch EPmain.sampleTest.method
                    case 1 %sample
                        EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames),1);
                        EPmain.sampleTest.templateTopo=zeros(length(EEGchans));
                        tMs=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames;
                    case 2 %CWT
                        tMs=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames;
                        theTau=tMs(floor(length(tMs)/2));
                        EPmain.sampleTest.templateWaveform=(1-16.*(((tMs-theTau)/(EPmain.sampleTest.scale/1000))/1000).^2).*exp(-8.*(((tMs-theTau)/(EPmain.sampleTest.scale/1000))/1000).^2); %The Mexican Hat wavelet
                        EPmain.sampleTest.templateTopo=zeros(length(EEGchans));
                    case 3 %AVE template
                        if ~isempty(EPmain.sampleTest.AVElist)
                            templateData=EPmain.sampleTest.AVEdata{EPmain.sampleTest.AVE}(EEGchans,EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd,EPmain.sampleTest.cell,EPmain.sampleTest.subject,:,EPmain.sampleTest.freqBin);
                            if EPmain.sampleTest.freqFlag
                                if (EPmain.sampleTest.freqScale == 1)
                                    templateData=real(templateData);
                                end;
                                if (EPmain.sampleTest.freqScale == 2)
                                    templateData=imag(templateData);
                                end;
                                if (EPmain.sampleTest.freqScale > 2)
                                    templateData(EEGchans,:,:,:)=abs(templateData(EEGchans,:,:,:)); %convert complex number to real number
                                end;
                                templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:)/mean(diff(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)); %convert to spectral density
                                if EPmain.sampleTest.freqScale > 3
                                    templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:).^2; %convert amplitude to power
                                end;
                                if (EPmain.sampleTest.freqScale == 5)
                                    if ~all(templateData(EEGchans,:,:,:) >=0)
                                        disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
                                    end;
                                    templateData(EEGchans,:,:,:)=log10(abs(templateData(EEGchans,:,:,:)))*10; %convert to dB log scaling
                                    tempVar=templateData(EEGchans,:,:,:);
                                    tempVar(isinf(tempVar))=-flintmax;
                                    templateData(EEGchans,:,:,:)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                                end;
                            end;
                            if ~isempty(EPmain.sampleTest.sampEnd)
                                templateData2=templateData(:,EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd);
                            else
                                templateData2=templateData;
                            end;
                            [A maxChan]=max(max(templateData2'));
                            [A maxPoint]=max(max(templateData2));
                            if ~isempty(EPmain.sampleTest.sampStart)
                                maxPoint=maxPoint+EPmain.sampleTest.sampStart-1;
                            end;
                            EPmain.sampleTest.templateWaveform=templateData(maxChan,:)';
                            EPmain.sampleTest.templateTopo=templateData(:,maxPoint);
                            tMs=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).timeNames;
                        else
                            noSampTest=1;
                            templateData=zeros(length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames),length(EEGchans));
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames),1);
                            EPmain.sampleTest.templateTopo=zeros(1,length(EEGchans));
                            tMs=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames;
                        end;
                    case 4 %PCA
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT)
                            EPmain.sampleTest.templateWaveform=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT(:,EPmain.sampleTest.factor);
                            tMs=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).timeNames;
                        else
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames),1);
                            tMs=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames;
                        end;
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS)
                            EPmain.sampleTest.templateTopo=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS(EEGchans,EPmain.sampleTest.factor);
                            if EPmain.sampleTest.test==2 %not just one channel
                                set(EPmain.handles.sampleTest.chanLabel,'enable','off');
                                set(EPmain.handles.sampleTest.channel,'enable','off');
                            end;
                        else
                            noSampTest=1;
                            EPmain.sampleTest.templateTopo=zeros(length(EEGchans));
                        end;
                    otherwise
                        disp('Oops: programmer error');
                        return;
                end;
            end;
            
            if EPmain.sampleTest.changeFlag==1
                EPmain.sampleTest.changeFlag=0;
%                 if strcmp(EPmain.sampleTest.dataType,'continuous')
                    EPmain.sampleTest.sampStart=1;
                    EPmain.sampleTest.sampEnd=length(EPmain.sampleTest.templateWaveform);
%                 end;
            end;
            
            %Waveform plot of template
            EPmain.handles.sampleTest.templatePlot = axes('units','pixels','position',[5 88 100 50]);
            EPmain.handles.sampleTest.templateWaves = plot(tMs(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd),EPmain.sampleTest.templateWaveform(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd));
            axis([tMs(EPmain.sampleTest.sampStart) tMs(EPmain.sampleTest.sampEnd) min([-1; EPmain.sampleTest.templateWaveform]) max([1; EPmain.sampleTest.templateWaveform])]);
            
            %2D plot of template
            theEloc=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).eloc;
            warning off
            [Xi,Yi,Zi] = griddata(EPmain.sampleTest.x,EPmain.sampleTest.y,EPmain.sampleTest.templateTopo,[1:EPmain.sampleTest.gridSize]',[1:EPmain.sampleTest.gridSize],'linear');
            warning on
            EPmain.handles.sampleTest.templateTopo = axes('units','pixels','position',[110 88 50 50]);
            EPmain.handles.sampleTest.templateTopoImage = imagesc(Zi);
            set(gca,'XTickLabel','','YTickLabel','');
            
            if strcmp(EPmain.sampleTest.dataType,'continuous')
                uicontrol('Style','text',...
                    'String','samples','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 55 50 20]);
                
                uicontrol('Style','text',...
                    'String','ms','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[75 55 50 20]);
                
                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.handles.sampleTest.sampStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1,...
                        'TooltipString','First sample of template in relation to the event, where negative is before it.',...
                        'Position',[5 35 35 20],'Callback',@sampleTestSampStart);
                    
                    EPmain.handles.sampleTest.sampEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline,...
                        'TooltipString','Last sample of template in relation to the event, where negative is before it.',...
                        'Position',[40 35 35 20],'Callback',@sampleTestSampEnd);
                    
                    EPmain.handles.sampleTest.msStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',(EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1)*(1000/EPmain.sampleTest.Fs),...
                        'TooltipString','Last ms of template in relation to the event, where negative is before it.',...
                        'Position',[75 35 35 20],'Callback',@sampleTestSampStart);
                    
                    EPmain.handles.sampleTest.msEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',(EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline)*(1000/EPmain.sampleTest.Fs),...
                        'TooltipString','First ms of template in relation to the event, where negative is before it.',...
                        'Position',[110 35 35 20],'Callback',@sampleTestSampEnd);
                else
                    EPmain.handles.sampleTest.sampStart= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','First sample of template in relation to the event, where negative is before it.',...
                        'Position',[5 35 35 20]);
                    
                    EPmain.handles.sampleTest.sampEnd= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','Last sample of template in relation to the event, where negative is before it.',...
                        'Position',[40 35 35 20]);
                    
                    EPmain.handles.sampleTest.msStart= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','Last ms of template in relation to the event, where negative is before it.',...
                        'Position',[75 35 35 20]);
                    
                    EPmain.handles.sampleTest.msEnd= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','First ms of template in relation to the event, where negative is before it.',...
                        'Position',[110 35 35 20]);
                end;
            end;
            
            EPmain.handles.sampleTest.sampleTest = uicontrol('Style', 'pushbutton', 'String', 'sampleTest','FontSize',EPmain.fontsize,...
                'Position', [70 0 60 35], 'Callback', 'ep(''sampleTest'');');
            
            if noSampTest
                set(EPmain.handles.sampleTest.sampleTest,'enable','off');
            end;
            
        else
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','No suitable data',...
                'Position',[5 480 160 20]);
            
        end;
        
        EPmain.handles.sampleTest.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','EPmain.sampleTest=[];','ep(''start'');']);
        
    case 'sampleTest'
        
        set(EPmain.handles.sampleTest.sampleTest,'enable','off');
        set(EPmain.handles.sampleTest.done,'enable','off');
        drawnow
        
        EPdata=ep_loadEPdataset(EPdataset,EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset));
        
        if ~strcmp(EPmain.sampleTest.dataType,'continuous')
            templateTopo=[];
            templateWaveform=[];
            switch EPmain.sampleTest.method
                case 1
                    method='sample';
                case 2
                    method='CWT';
                case 3
                    method='Template Woody';
                    templateTopo=EPmain.sampleTest.templateTopo;
                    templateWaveform=EPmain.sampleTest.templateWaveform;
                case 4
                    method='PCA Woody';
                    facVecT=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT;
                    if ~isempty(facVecT)
                        templateWaveform=facVecT(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd,EPmain.sampleTest.factor);
                    else
                        templateWaveform=[];
                    end;
                    facVecS=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS;
                    if ~isempty(facVecS)
                        templateTopo=facVecS(:,EPmain.sampleTest.factor);
                    else
                        templateTopo=[];
                    end;
                otherwise
                    disp('Programmer error');
                    return;
            end;
            
            switch EPmain.sampleTest.test
                case 1
                    test='t-test';
                case 2
                    test='jack-knife';
                otherwise
                    disp('Programmer error');
                    return;
            end;
            
            cellList=find(ismember(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).cellNames,{EPmain.sampleTest.cellNameList{EPmain.sampleTest.cell1},EPmain.sampleTest.cellNameList{EPmain.sampleTest.cell2}}));
            if ~isempty(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)
                if EPmain.sampleTest.freqBin > length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)
                    freqList=[1:length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)];
                else
                    freqList=EPmain.sampleTest.freqBin;
                end;
            else
                freqList=[];
            end;
            EPdataST=ep_selectData(ep_stripAdds(EPdata),{[],[],cellList,[],[],[freqList]});
            if isempty(EPdataST)
                msg{1}='No regular data to analyze.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            EEGchans=find(strcmp('EEG',EPdataST.chanTypes));
            if isempty(EEGchans)
                msg{1}='No EEG channels to analyze.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            if isempty(EPdataset.dataset(EPmain.sampleTest.dataset).facVecT) && ~isempty(EPmain.sampleTest.filter1) %cannot filter temporal PCA data
                if (EPmain.sampleTest.filterPass==3) && (EPmain.sampleTest.filter1 >= EPmain.sampleTest.filter2)
                    msg{1}='The start of the bandpass window must come before the end of the bandpass window.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                
                if (EPmain.sampleTest.filterPass==4) && (EPmain.sampleTest.filter1 >= EPmain.sampleTest.filter2)
                    msg{1}='The start of the bandstop window must come before the end of the bandstop window.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                
                if (length(EPmain.sampleTest.filter1) > 1) || (length(EPmain.sampleTest.filter2) > 1) || (length(EPmain.sampleTest.filterOrder) > 1)
                    msg{1}='The filter setting in each field needs to be a single number.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                
                if (EPmain.sampleTest.filter1 < 0) || (EPmain.sampleTest.filter2 < 0) || (EPmain.sampleTest.filterOrder < 0)
                    msg{1}='Filter settings cannot be negative numbers.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                
                if strcmp(test,'jack-knife') && ~isempty(freqList)
                    msg{1}='Jack-knife can only be performed on a single frequency band.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                
                cfg=[];
                switch EPmain.sampleTest.filterType
                    case 1
                        filterDirection='onepass';
                        filterType='but';
                    case 2
                        filterDirection='twopass';
                        filterType='but';
                    case 3
                        filterDirection='onepass';
                        filterType='fir';
                    case 4
                        filterDirection='twopass';
                        filterType='fir';
                    case 5
                        filterDirection='onepass';
                        filterType='firls';
                    case 6
                        filterDirection='twopass';
                        filterType='firls';
                end;
                switch EPmain.sampleTest.filterPass
                    case 1 %low pass
                        cfg.lpfilter='yes';
                        cfg.lpfreq=EPmain.sampleTest.filter1;
                        cfg.lpfiltdir=filterDirection;
                        cfg.lpfilttype=filterType;
                        cfg.lpfiltord=EPmain.sampleTest.filterOrder;
                    case 2 %high pass
                        cfg.hpfilter='yes';
                        cfg.hpfreq=EPmain.sampleTest.filter1;
                        cfg.hpfiltdir=filterDirection;
                        cfg.hpfilttype=filterType;
                        cfg.hpfiltord=EPmain.sampleTest.filterOrder;
                    case 3 %band pass
                        cfg.bpfilter='yes';
                        cfg.bpfreq=[EPmain.sampleTest.filter1 EPmain.sampleTest.filter2];
                        cfg.bpfiltdir=filterDirection;
                        cfg.bpfilttype=filterType;
                        cfg.bpfiltord=EPmain.sampleTest.filterOrder;
                    case 4 %band stop
                        cfg.bsfilter='yes';
                        cfg.bsfreq=[EPmain.sampleTest.filter1 EPmain.sampleTest.filter2];
                        cfg.bsfiltdir=filterDirection;
                        cfg.bsfilttype=filterType;
                        cfg.bsfiltord=EPmain.sampleTest.filterOrder;
                    case 5 %notch
                        cfg.dftfilter='yes';
                        cfg.dftfreq=EPmain.sampleTest.filter1;
                        cfg.dftfiltdir=filterDirection;
                        cfg.dftfilttype=filterType;
                        cfg.dftfiltord=EPmain.sampleTest.filterOrder;
                end;
                EPdataST=ep_filterData(EPdataST,cfg,EEGchans);
            end;
            
            [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdataST,method,EPmain.sampleTest.alpha,EPmain.sampleTest.contiguous,EPmain.sampleTest.scale,templateWaveform,templateTopo,test,EPmain.sampleTest.thresh,EPmain.sampleTest.channel);
            if isempty(outputData)
                set(EPmain.handles.sampleTest.sampleTest,'enable','on');
                set(EPmain.handles.sampleTest.done,'enable','on');
                drawnow
                return
            end;
            
            %add the results to the dataset if t-test
            if strcmp(test,'t-test')
                sampleTestCellName=[method ':' EPmain.sampleTest.cellNameList{EPmain.sampleTest.cell1} '-' EPmain.sampleTest.cellNameList{EPmain.sampleTest.cell2}];
                
                if ~any(strcmp(sampleTestCellName,EPdata.cellNames))
                    EPadd=[];
                    EPadd.cellTypes{1}='STS';
                    if strcmp(EPdata.dataType,'single_trial')
                        EPadd.trialNames=1;
                    else
                        EPadd.trialNames=[];
                    end;
                    EPadd.cellNames{1}=sampleTestCellName;
                    [EPdata]=ep_addData(EPdata,EPadd,'cells');
                    if isempty(EPdata)
                        return
                    end;
                else
                    disp('Overwriting existing sample test results');
                end;
                if isempty(freqList)
                    theFreq=1;
                else
                    theFreq=freqList;
                end;
                if strcmp(EPdata.dataType,'average')
                    if ~any(strcmp('sampleTest',EPdata.subNames))
                        EPadd=[];
                        EPadd.subTypes{1}='GAV';
                        EPadd.subNames{1}='sampleTest';
                        [EPdata]=ep_addData(EPdata,EPadd,'subjects');
                        if isempty(EPdata)
                            return
                        end;
                    end
                    EPdata.data(EEGchans,:,strcmp(sampleTestCellName,EPdata.cellNames),strcmp('sampleTest',EPdata.subNames),:,theFreq,:)=outputData;
                else %single_trial
                    EPdata.data(EEGchans,:,strcmp(sampleTestCellName,EPdata.cellNames),:,:,theFreq,:)=outputData;
                end
                
                %add latency information if Woody PCA and if filter was spatial PCA and thus a single latency for all the channels
                if strcmp('PCA Woody',method) && ~isempty(templateTopo) && strcmp('single_trial',EPdata.dataType) && (length(freqList)==1)
                    cellNames=unique(EPdataST.cellNames);
                    goodCellTrials1=find(EPdata.analysis.badTrials(1,strcmp(cellNames{1},EPdata.cellNames))==0);
                    goodCellTrials2=find(EPdata.analysis.badTrials(1,strcmp(cellNames{2},EPdata.cellNames))==0);
                    ce11Trials1=find(strcmp(cellNames{1},EPdata.cellNames));
                    ce11Trials2=find(strcmp(cellNames{2},EPdata.cellNames));
                    goodTrials1=ce11Trials1(goodCellTrials1);
                    goodTrials2=ce11Trials2(goodCellTrials2);
                    
                    trialSpecName=['Woody PCA latency: ' EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facNames{EPmain.sampleTest.factor}];
                    
                    if ~any(strcmp(trialSpecName,EPdata.trialSpecNames))
                        EPdata.trialSpecNames{end+1}=trialSpecName;
                        EPdata.trialSpecs(:,end+1)=cell(length(EPdata.cellNames),1);
                    end;
                    
                    EPdata.trialSpecs(goodTrials1,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampLat1);
                    EPdata.trialSpecs(goodTrials2,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampLat2);
                    
                    trialSpecName=['Woody PCA amplitude: ' EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facNames{EPmain.sampleTest.factor}];
                    if ~any(strcmp(trialSpecName,EPdata.trialSpecNames))
                        EPdata.trialSpecNames{end+1}=trialSpecName;
                        EPdata.trialSpecs(:,end+1)=cell(length(EPdata.cellNames),1);
                    end;
                    EPdata.trialSpecs(goodTrials1,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampAmp1);
                    EPdata.trialSpecs(goodTrials2,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampAmp2);
                end;
            end;
        else %continuous
            
            switch EPmain.sampleTest.method
                case 1
                    method='Template Woody';
                    templateTopo=EPmain.sampleTest.templateTopo;
                    templateWaveform=EPmain.sampleTest.templateWaveform(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd);
                    sampleTestChanName=[method ':' EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).subNames{EPmain.sampleTest.subject} '-' EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).cellNames{EPmain.sampleTest.cell}];
                case 2
                    method='PCA Woody';
                    sampleTestChanName=[method ':' EPdata.facNames{EPmain.sampleTest.factor}];
                    facVecT=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT;
                    if ~isempty(facVecT)
                        templateWaveform=facVecT(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd,EPmain.sampleTest.factor);
                    else
                        templateWaveform=[];
                    end;
                    facVecS=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS;
                    if ~isempty(facVecS)
                        templateTopo=facVecS(:,EPmain.sampleTest.factor);
                    else
                        templateTopo=[];
                    end;
                otherwise
                    disp('Programmer error');
                    return;
            end;
            
            test='';
            
            [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdata,method,EPmain.sampleTest.alpha,EPmain.sampleTest.contiguous,EPmain.sampleTest.scale,templateWaveform,templateTopo,test,EPmain.sampleTest.thresh,EPmain.sampleTest.channel);
            if isempty(outputData)
                set(EPmain.handles.sampleTest.sampleTest,'enable','on');
                set(EPmain.handles.sampleTest.done,'enable','on');
                drawnow
                return
            end;
            
            if ~any(strcmp(sampleTestChanName,EPdata.chanNames))
                EPadd=[];
                EPadd.chanTypes{1}='REG';
                EPadd.chanNames{1}=sampleTestChanName;
                [EPdata]=ep_addData(EPdata,EPadd,'channels');
                if isempty(EPdata)
                    return
                end;
            else
                disp('Overwriting existing Woody waveform');
            end;
            
            EPdata.data(strcmp(sampleTestChanName,EPdata.chanNames),:,:,:,:,:,:)=outputData;
        end
        
        %update the copy of the dataset in the active set
        disp('Saving results to the working set.');
        whichData=EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset);
        delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(whichData).dataName '.mat']);
        EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],whichData));
        EPdataset=ep_saveEPdataset(EPdataset,EPdata,length(EPdataset.dataset)+1,'no');
        EPdataset.dataset=[EPdataset.dataset(1:whichData-1) EPdataset.dataset(end) EPdataset.dataset(whichData:end-1)];
        
        set(EPmain.handles.sampleTest.sampleTest,'enable','on');
        set(EPmain.handles.sampleTest.done,'enable','on');
        drawnow
        
    case 'startPCA'
        
        set(EPmain.handles.hMainWindow,'Name', 'PCA');
        
        h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Mode','HorizontalAlignment','left',...
            'Position',[25 470 150 20]);
        
        EPmain.handles.pca.mode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'spatial','temporal','frequency','cross-verification'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.mode,''Value'');','if tempVar ~=0,EPmain.pca.mode=tempVar;end;','if isempty(tempVar),EPmain.pca.mode=tempVar;end;EPmain.pca.crossVerifyPCA=[];','ep(''start'');'],...
            'Value',EPmain.pca.mode,'Position',[20 450 160 20]);
        
        if EPmain.pca.mode == 4
            %cross-verification mode
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','PCA Datasets','HorizontalAlignment','left',...
                'Position',[25 430 150 20]);
            
            crossVerifyTableData=cell(0);
            EPmain.pca.PCAdatasets=[];
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    if ~isempty(EPdataset.dataset(i).facNames)
                        fileName=EPdataset.dataset(i).dataName;
                        if strcmp(EPdataset.dataset(i).saved,'no')
                            fileName=['*' fileName];
                        end;
                        EPmain.pca.PCAdatasets(end+1)=i;
                        crossVerifyTableData{length(EPmain.pca.PCAdatasets),1}=fileName;
                    end;
                end;
            else
                crossVerifyTableData=cell(0);
            end;
            
            tableNames{1}='data';
            columnEditable =  false;
            ColumnFormat{1}=[];
            
            EPmain.handles.crossVerify.hTable = uitable('Data',crossVerifyTableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',@crossVerifyData,...
                'ColumnWidth',{300},'Position',[20 280 170 150]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Chosen:','HorizontalAlignment','left',...
                'Position',[25 250 45 20]);            

            if isempty(EPmain.pca.crossVerifyPCA)
                theText='none';
            else
                theText=EPdataset.dataset(EPmain.pca.crossVerifyPCA).dataName;
            end
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',theText,'HorizontalAlignment','left',...
                'Position',[75 250 150 20]);            

            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Title of PCA','HorizontalAlignment','left',...
                'Position',[25 210 150 20]);
            
            EPmain.handles.pca.name= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.name,...
                'Callback', ['global EPmain;','EPmain.pca.name=get(EPmain.handles.pca.name,''String'');'],...
                'Position',[25 190 150 20]);
            
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    fileName=EPdataset.dataset(i).dataName;
                    if strcmp(EPdataset.dataset(i).saved,'no')
                        fileName=['*' fileName];
                    end;
                    tableData{i,1}=fileName;
                end;
            else
                tableData=[];
            end;
            
            tableData=cell(0);
            EPmain.pca.targetDatasets=[];
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    dropFlag=0;
                    if ~isempty(EPmain.pca.crossVerifyPCA)
                        if ~isempty(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecT)
                            if ~isempty(EPdataset.dataset(i).facVecT)
                                dropFlag=1;
                            end;
                            if size(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecT,1) ~= length(EPdataset.dataset(i).timeNames)
                                dropFlag=1;
                            end;
                        end;
                        if ~isempty(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecS)
                            if ~isempty(EPdataset.dataset(i).facVecS)
                                dropFlag=1;
                            end;
                            if size(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecS,1) ~= length(EPdataset.dataset(i).chanNames(strcmp('EEG',EPdataset.dataset(i).chanTypes)))
                                dropFlag=1;
                            end;
                        end;
                        if ~isempty(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecF)
                            if ~isempty(EPdataset.dataset(i).facVecF)
                                dropFlag=1;
                            end;
                            if size(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecF,1) ~= length(EPdataset.dataset(i).freqNames)
                                dropFlag=1;
                            end;
                        end;
                    end;
                    if ~dropFlag
                        fileName=EPdataset.dataset(i).dataName;
                        if strcmp(EPdataset.dataset(i).saved,'no')
                            fileName=['*' fileName];
                        end;
                        EPmain.pca.targetDatasets(end+1)=i;
                        tableData{length(EPmain.pca.targetDatasets),1}=fileName;
                    end;
                end;
            else
                tableData=cell(0);
            end;

            tableNames{1}='data';
            columnEditable =  false;
            ColumnFormat{1}=[];
            
            EPmain.handles.pca.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',@pickPCAdata,...
                'ColumnWidth',{300},'Position',[20 40 170 150]);
        else
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Rotation','HorizontalAlignment','left',...
                'Position',[25 430 150 20]);
            
            EPmain.handles.pca.rotation= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'Varimax','Promax','Infomax','Quartimax','Quartimin','Oblimin','CRFE','MINE','IPSC','TIIC','Geomin','MMER','VOMN','unrotated'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.rotation,''Value'');','if tempVar ~=0,EPmain.pca.rotation=tempVar;end;','if isempty(tempVar),EPmain.pca.rotation=tempVar;end;','ep(''start'');'],...
                'Value',EPmain.pca.rotation,'Position',[20 410 160 20]);
            
            h = uicontrol('Style','text',...
                'String','Rotation Parameter','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 390 150 20]);
            
            EPmain.handles.pca.rotopt= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.rotopt,...
                'CallBack',['global EPmain;','EPmain.pca.rotopt=get(EPmain.handles.pca.rotopt,''String'');','ep(''start'');'],...
                'Position',[25 370 150 20]);
            
            h = uicontrol('Style','text',...
                'String','Relationship matrix','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 350 150 20]);
            
            EPmain.handles.pca.rel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'correlation','covariance'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.rel,''Value'');','if tempVar ~=0,EPmain.pca.rel=tempVar;end;','if isempty(tempVar),EPmain.pca.rel=tempVar;end;','ep(''start'');'],...
                'Value',EPmain.pca.rel,'Position',[20 330 160 20]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Loading Weighting','HorizontalAlignment','left',...
                'Position',[25 310 150 20]);
            
            EPmain.handles.pca.loadings = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'none','Kaiser','covariance','C-M'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.loadings,''Value'');','if tempVar ~=0,EPmain.pca.loadings=tempVar;end;','if isempty(tempVar),EPmain.pca.loadings=tempVar;end;','ep(''start'');'],...
                'Value',EPmain.pca.loadings,'Position',[20 290 160 20]);
            
            if get(EPmain.handles.pca.rotation,'Value') == 3 %infomax rotation
                set(EPmain.handles.pca.rotopt,'enable','off');
                set(EPmain.handles.pca.rel,'enable','off');
                set(EPmain.handles.pca.loadings,'enable','off');
            end;
            
            EPmain.handles.pca.parametric= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
                'String','Parametric Analysis',...
                'Value',EPmain.pca.parametric,'Position',[20 270 150 20]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','# Factors (0=scree)','HorizontalAlignment','left',...
                'Position',[25 250 150 20]);
            
            EPmain.handles.pca.facNum= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.facNum,...
                'Callback', ['global EPmain;','EPmain.pca.facNum=str2num(get(EPmain.handles.pca.facNum,''String''));'],...
                'Position',[25 230 150 20]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Title of PCA','HorizontalAlignment','left',...
                'Position',[25 210 150 20]);
            
            EPmain.handles.pca.name= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.name,...
                'Callback', ['global EPmain;','EPmain.pca.name=get(EPmain.handles.pca.name,''String'');'],...
                'Position',[25 190 150 20]);
            
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    fileName=EPdataset.dataset(i).dataName;
                    if strcmp(EPdataset.dataset(i).saved,'no')
                        fileName=['*' fileName];
                    end;
                    tableData{i,1}=fileName;
                end;
            else
                tableData=[];
            end;
            
            tableNames{1}='data';
            columnEditable =  false;
            ColumnFormat{1}=[];
            
            EPmain.handles.pca.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',@pickPCAdata,...
                'ColumnWidth',{300},'Position',[20 40 170 150]);
            
        end;
        
        EPmain.handles.pca.hQuitRead = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';EPmain.pca.crossVerifyPCA=[];','ep(''start'');']);
        
    case 'startWindow'
        
        set(EPmain.handles.hMainWindow,'Name', 'Window Data');
        
        %initialize Window Data parameters if entering in from the Main Menu.
        if ~isfield(EPmain.window,'measure') || (EPmain.window.dataset > length(EPdataset.dataset)) || ~strcmp(EPmain.window.datasetName,EPdataset.dataset(EPmain.window.dataset).dataName)
            
            if isfield(EPmain.window,'dataset')
                if (EPmain.window.dataset > length(EPdataset.dataset)) || (~isempty(EPmain.window.datasetName) && ~strcmp(EPmain.window.datasetName,EPdataset.dataset(EPmain.window.dataset).dataName))
                    disp('The working set has been changed so reinitializing the Window Data pane.')
                end;
            end;
            
            EPmain.window.minFacVar=EPmain.preferences.window.minFacVar;
            EPmain.window.measure=1;
            EPmain.window.chanGrp=1;
            EPmain.window.specSelect(1)=false;
            EPmain.window.factor=1;
            EPmain.window.FFTunits=4;
            EPmain.window.sampAdapt=0;
            if ~isfield(EPmain.window,'dataset')
                EPmain.window.dataset=1;
            end;
            
            EPmain.window.dataNames=[];
            EPmain.window.factorData=[];
            EPmain.window.aveData=[];
            factorDataNum=0;
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    if ~strcmp(EPdataset.dataset(i).dataType,'continuous') && ~isempty(find(strcmp('EEG',EPdataset.dataset(i).chanTypes))) &&...
                            ~isempty(find(ismember({'RAW','AVG'},EPdataset.dataset(i).subTypes))) &&...
                            ~isempty(find(strcmp('SGL',EPdataset.dataset(i).cellTypes)))
                        EPmain.window.aveData(end+1)=i; %keep track of which datasets are suitable for windowing
                        EPmain.window.dataNames{i,1}=EPdataset.dataset(i).dataName;
                    end;
                    
                    if isfield(EPdataset.dataset(i).pca,'PCAmode')
                        if isfield(EPdataset.dataset(i).pca,'PCAmode2')
                            if strcmp(EPdataset.dataset(i).pca.PCAmode2,'spat')
                                factorDataNum=factorDataNum+1;
                                EPmain.window.factorData(factorDataNum).name=EPdataset.dataset(i).dataName;
                                EPmain.window.factorData(factorDataNum).FacPat=EPdataset.dataset(i).pca.FacPatST;
                                facNames=EPdataset.dataset(i).facNames(find(strcmp('SGL',EPdataset.dataset(i).facTypes)));
                                if length(facNames)==size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                    EPmain.window.factorData(factorDataNum).facNames=facNames;
                                else
                                    for fac =1:size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                        EPmain.window.factorData(factorDataNum).facNames{fac}=['fac' num2str(fac)];
                                    end;
                                end;
                            end;
                        end;
                        if strcmp(EPdataset.dataset(i).pca.PCAmode,'spat')
                            factorDataNum=factorDataNum+1;
                            EPmain.window.factorData(factorDataNum).name=EPdataset.dataset(i).dataName;
                            EPmain.window.factorData(factorDataNum).FacPat=EPdataset.dataset(i).pca.FacPat;
                            facNames=EPdataset.dataset(i).facNames(find(strcmp('SGL',EPdataset.dataset(i).facTypes)));
                            if length(facNames)==size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                EPmain.window.factorData(factorDataNum).facNames=facNames;
                            else
                                for fac =1:size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                    EPmain.window.factorData(factorDataNum).facNames{fac}=['fac' num2str(fac)];
                                end;
                            end;
                        end;
                    end;
                end;
                EPmain.window.dataset=EPmain.window.aveData(1); %the active dataset
                EPmain.window.datasetName=EPdataset.dataset(EPmain.window.dataset).dataName;
            else
                EPmain.window.dataNames='none';
            end;
            EPmain.window.sampStart=1;
            EPmain.window.sampEnd=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
            if isempty(EPmain.window.sampEnd)
                EPmain.window.sampEnd=1;
            end;
            EPmain.window.HzStart=1;
            EPmain.window.HzEnd=length(EPdataset.dataset(EPmain.window.dataset).freqNames);
            if isempty(EPmain.window.HzEnd)
                EPmain.window.HzEnd=1;
            end;
            cellNames=EPdataset.dataset(EPmain.window.dataset).cellNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes)));
            [u i]=unique(cellNames,'first');
            EPmain.window.inCells=cellNames(sort(i));
            EPmain.window.outCells=EPmain.window.inCells;
        end;
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Dataset','HorizontalAlignment','left',...
            'Position',[25 470 140 20]);
        
        EPmain.handles.window.dataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.window.dataNames(EPmain.window.aveData),...
            'Value',find(EPmain.window.dataset==EPmain.window.aveData),'Position',[5 450 190 20],...
            'Callback', @changeWindowDataset);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.window.dataset,'enable','off');
        end;
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Measure','HorizontalAlignment','left',...
            'Position',[25 430 140 20]);
        
        if isempty(EPdataset.dataset(EPmain.window.dataset).timeNames) && ~isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
            measures={'mean','minHzPeak','maxHzPeak','minHzLatency','maxHzLatency','minHzCentroid','maxHzCentroid'};
        elseif ~isempty(EPdataset.dataset(EPmain.window.dataset).timeNames) && isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
            measures={'mean','minpeak','maxpeak','minlatency','maxlatency','mincentroid','maxcentroid'};
        else
            measures={'mean','minpeak','maxpeak','minlatency','maxlatency','mincentroid','maxcentroid','minHzPeak','maxHzPeak','minHzLatency','maxHzLatency','minHzCentroid','maxHzCentroid'};
        end;
        if ~isempty(EPdataset.dataset(EPmain.window.dataset).trialSpecNames) && any(ismember({'RT','ACC'},EPdataset.dataset(EPmain.window.dataset).trialSpecNames))
            measures{end+1}='behavioral';
        end;
        
        EPmain.handles.window.measure= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',measures,...
            'Value',EPmain.window.measure,'Position',[5 410 110 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.window.measure,''Value'');','if tempVar ~=0,EPmain.window.measure=tempVar;end;','if isempty(tempVar),EPmain.window.measure=tempVar;end;','ep(''start'');']);
        
        if strcmp(measures(EPmain.window.measure),'behavioral')
            
            
            
            
        else
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','samples','HorizontalAlignment','left',...
                'Position',[110 430 40 20]);
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','ms','HorizontalAlignment','left',...
                'Position',[155 430 40 20]);
            
            EPmain.handles.window.sampAdapt= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.sampAdapt,...
                'Position',[110 410 40 20],'Callback',@windowSampAdapt,...
                'TooltipString','Samples surrounding peak that are averaged together for peak measure (0 for peak only).');
            
            EPmain.handles.window.msAdapt= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',round(EPmain.window.sampAdapt*(1000/EPdataset.dataset(EPmain.window.dataset).Fs)),...
                'Position',[150 410 40 20],'Callback',@windowSampAdapt,...
                'TooltipString','Ms surrounding peak that are averaged together for peak measure (0 for peak only).');
            
            if ~any(strcmp(measures{EPmain.window.measure},{'minpeak','maxpeak'}))
                set(EPmain.handles.window.sampAdapt,'enable','off');
                set(EPmain.handles.window.msAdapt,'enable','off');
            end;
            
            uicontrol('Style','text',...
                'String','samples','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 390 140 20]);
            
            if strcmp(EPdataset.dataset(EPmain.window.dataset).timeUnits,'per')
                theUnit='%';
            else
                theUnit='ms';
            end;
            uicontrol('Style','text',...
                'String',theUnit,'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[120 390 140 20]);
            
            EPmain.handles.window.sampStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.sampStart,...
                'Position',[25 370 40 20],'Callback',@windowSampStart);
            
            EPmain.handles.window.sampEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.sampEnd,...
                'Position',[60 370 40 20],'Callback',@windowSampEnd);
            
            EPmain.handles.window.msStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',round((EPmain.window.sampStart-EPdataset.dataset(EPmain.window.dataset).baseline-1)*(1000/EPdataset.dataset(EPmain.window.dataset).Fs)),...
                'Position',[100 370 40 20],'Callback',@windowSampStart);
            
            EPmain.handles.window.msEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',round((EPmain.window.sampEnd-EPdataset.dataset(EPmain.window.dataset).baseline)*(1000/EPdataset.dataset(EPmain.window.dataset).Fs)),...
                'Position',[140 370 40 20],'Callback',@windowSampEnd);
            
            if isempty(EPdataset.dataset(EPmain.window.dataset).timeNames)
                set(EPmain.handles.window.sampStart,'enable','off');
                set(EPmain.handles.window.sampEnd,'enable','off');
                set(EPmain.handles.window.msStart,'enable','off');
                set(EPmain.handles.window.msEnd,'enable','off');
            end;
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Bins','HorizontalAlignment','left',...
                'Position',[25 350 140 20]);
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Hz','HorizontalAlignment','left',...
                'Position',[120 350 140 20]);
            
            EPmain.handles.window.binStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.HzStart,...
                'Position',[25 330 40 20],'Callback',@windowHzStart);
            
            EPmain.handles.window.binEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.HzEnd,...
                'Position',[60 330 40 20],'Callback',@windowHzEnd);
            
            HzStart=0;
            HzEnd=0;
            if ~isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
                HzStart=round(EPdataset.dataset(EPmain.window.dataset).freqNames(EPmain.window.HzStart)*10)/10;
                HzEnd=round(EPdataset.dataset(EPmain.window.dataset).freqNames(EPmain.window.HzEnd)*10)/10;
            end;
            EPmain.handles.window.HzStart= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',sprintf('%5.1f',HzStart),...
                'Position',[100 330 40 20]);
            
            EPmain.handles.window.HzEnd= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',sprintf('%5.1f',HzEnd),...
                'Position',[140 330 40 20]);
            
            h = uicontrol('Style','text','HorizontalAlignment','left','String', 'units','FontSize',EPmain.fontsize,...
                'ForegroundColor','black','Position',[160 310 25 20]);
            EPmain.handles.window.FFTunits = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'Value',EPmain.window.FFTunits,'Position',[155 290 65 20],...
                'String',{'cm','am','pw','dB'},...
                'Callback', ['global EPmain;','EPmain.window.FFTunits=get(EPmain.handles.window.FFTunits,''value'');','ep(''start'')'],...
                'TooltipString','Units for spectral data.');
            
            if isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
                set(EPmain.handles.window.HzStart,'enable','off');
                set(EPmain.handles.window.HzEnd,'enable','off');
                set(EPmain.handles.window.binStart,'enable','off');
                set(EPmain.handles.window.binEnd,'enable','off');
                set(EPmain.handles.window.FFTunits,'enable','off');
            elseif ~isempty(EPdataset.dataset(EPmain.window.dataset).relNames)
                set(EPmain.handles.window.FFTunits,'enable','off');
            end;
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Baseline:','HorizontalAlignment','left',...
                'Position',[10 280 70 20]);
            
            EPmain.handles.window.baseline= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',(EPdataset.dataset(EPmain.window.dataset).baseline)*(1000/EPdataset.dataset(EPmain.window.dataset).Fs),'HorizontalAlignment','left',...
                'Position',[60 280 30 20]);
            
            if isfield(EPchanGrp,'group')
                if isempty(EPchanGrp.group)
                    chanGrpNames='-all-';
                else
                    if isempty({EPchanGrp.group.name})
                        chanGrpNames='-all-';
                    else
                        chanGrpNames=[{EPchanGrp.group.name} '-all-'];
                    end;
                end;
            else
                chanGrpNames='-all-';
            end;
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Channels','HorizontalAlignment','left',...
                'Position',[10 310 50 20]);
            
            EPmain.handles.window.chanGrp = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',chanGrpNames,...
                'Value',EPmain.window.chanGrp,'Position',[60 310 100 20],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.window.chanGrp,''Value'');','if tempVar ~=0,EPmain.window.chanGrp=tempVar;end;','if isempty(tempVar),EPmain.window.chanGrp=tempVar;end;','ep(''start'');']);
            
            EPmain.handles.window.channels = uicontrol('Style', 'pushbutton', 'String', 'Channels','FontSize',EPmain.fontsize,...
                'Position', [90 275 60 35], 'Callback', ['global EPmain EPdataset;,ep_chanGrp(ep_loadEPdataset(EPdataset,EPmain.window.dataset),EPmain.window.factorData);']);
            
            if isempty(EPdataset.dataset)
                set(EPmain.handles.window.chanGrp,'enable','off');
            end;
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Factor','HorizontalAlignment','left',...
                'Position',[20 255 50 20]);
            
            if ~isempty(EPdataset.dataset(EPmain.window.dataset).facNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).facTypes))))
                EPmain.handles.window.factor = uicontrol('Style', 'popupmenu', 'String', EPdataset.dataset(EPmain.window.dataset).facNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).facTypes))),'FontSize',EPmain.fontsize,...
                    'Value',EPmain.window.factor,'Position', [70 240 100 35],...
                    'Callback', ['global EPmain;','EPmain.window.factor=get(EPmain.handles.window.factor,''Value'');','ep(''start'');']);
            else
                EPmain.handles.window.factors = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [70 240 100 35],'enable','off');
            end;
            
        end;
        
        %cell selection
        if ~isempty(find(strcmp('EEG',EPdataset.dataset(EPmain.window.dataset).chanTypes))) && ~isempty(find(ismember({'RAW','AVG'},EPdataset.dataset(EPmain.window.dataset).subTypes)))...
                && ~isempty(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes))) && ~(~isempty(EPdataset.dataset(EPmain.window.dataset).facNames) && isempty(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes))))
            cellNames=unique(EPdataset.dataset(EPmain.window.dataset).cellNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes))),'first');
            cellNames{end+1}='none';
            for i=1:length(cellNames)-1
                tableData(i,1)=EPmain.window.outCells(i);
                tableData(i,2)=EPmain.window.inCells(i);
            end;
            
            tableNames{1}='outCells';
            tableNames{2}='inCells';
            columnEditable=[true true];
            ColumnFormat{1}='char';
            ColumnFormat{2}='char';
            
            EPmain.handles.window.cellTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'ColumnWidth',{50 90},...
                'Position',[10 125 190 130],...
                'CellEditCallback', @windowCellTable);
            
        else
            uicontrol('Style','text','HorizontalAlignment','left','String', 'This average file cannot be windowed.',...
                'Position',[10 125 190 130]);
        end;
        
        %subject spec selection
        subjectSpecNames=EPdataset.dataset(EPmain.window.dataset).subjectSpecNames;
        if length(EPmain.window.specSelect) ~= length(subjectSpecNames)
            EPmain.window.specSelect=repmat(false,length(subjectSpecNames),1);
        end;
        if ~isempty(subjectSpecNames)
            tableData=[];
            for i=1:length(subjectSpecNames)
                tableData{i,1}=EPmain.window.specSelect(i);
                tableData{i,2}=subjectSpecNames{i};
            end;
        else
            tableData=[];
        end;
        
        tableNames{1}='Select';
        tableNames{2}='Specs';
        columnEditable=[true false];
        ColumnFormat{1}='logical';
        ColumnFormat{2}='char';
        
        EPmain.handles.window.specsTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'ColumnWidth',{40 100},...
            'Position',[10 40 190 80],...
            'CellEditCallback', @windowSpecsTable);
        
        EPmain.handles.window.window = uicontrol('Style', 'pushbutton', 'String', 'Window','FontSize',EPmain.fontsize,...
            'Position', [10 0 60 35], 'Callback', 'ep(''WindowData'')');
        
        EPmain.handles.window.automatic = uicontrol('Style', 'pushbutton', 'String', 'AutoPCA','FontSize',EPmain.fontsize,...
            'Position', [70 0 60 35], 'Callback', 'ep(''AutoPCA'')');
        
        if isempty(EPdataset.dataset(EPmain.window.dataset).facNames)
            set(EPmain.handles.window.automatic,'enable','off');
        end;
        
        EPmain.handles.window.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'WindowData'
        %output the windowed measure to a file for later ANOVA
        
        measureList=get(EPmain.handles.window.measure,'String');
        measureNum = EPmain.window.measure;
        measure=measureList{measureNum};
        
        set(EPmain.handles.window.window,'enable','off');
        set(EPmain.handles.window.done,'enable','off');
        if  ~strcmp(measure,'behavioral')
            set(EPmain.handles.window.channels,'enable','off');
            set(EPmain.handles.window.automatic,'enable','off');
            set(EPmain.handles.window.FFTunits,'enable','off');
        end;
        drawnow
        
        EPdata=ep_stripAdds(ep_loadEPdataset(EPdataset,EPmain.window.dataset));
        sampStart = EPmain.window.sampStart;
        sampEnd = EPmain.window.sampEnd;
        HzStart = EPmain.window.HzStart;
        HzEnd = EPmain.window.HzEnd;
        factor = EPmain.window.factor;
        theChanGrp = EPmain.window.chanGrp;
        
        if isempty(EPchanGrp) || theChanGrp > length(EPchanGrp.group)
            chanGrp.channel=[1:length(EPdata.chanNames)];
            chanGrp.name='-all-';
            chanGrp.areaName=[EPdata.chanNames; 'none'];
        else
            chanGrp=EPchanGrp.group(theChanGrp);
        end;
        
        areaList=setdiff(unique(chanGrp.channel),length(chanGrp.areaName));
        if length(areaList) > 1 && ~isempty(EPdata.facVecS)
            warndlg('Each spatial PCA factor is a virtual electrode covering the entire head.  The electrode regions would be 100% correlated so the ANOVA would fail.');
        end;
        
        if ~isempty(EPdataset.dataset(EPmain.window.dataset).relNames)
            disp('For coherence data, if there is one channel in an area then the correlations with all other channels will be used.');
            disp('The absolute value will be taken since otherwise they would tend to sum to zero (at least if average referenced).');
            disp('If a single common reference is used, it will be excluded since coherence with a flat channel produces a not-a-number result.');
            disp('If there is more than one channel in an area, then only the correlations between the channels in that area will be used.');
            disp('The absolute value will not be taken since the sign is meaningful and will not necessarily sum to zero.');
        end;
        
        inputcells=[];
        outCellNames=[];
        for iCell=1:length(EPmain.window.inCells)
            if ~isempty(EPmain.window.inCells{iCell}) && ~isempty(EPmain.window.outCells{iCell})
                if ~any(strcmp(EPmain.window.outCells(iCell),outCellNames))
                    outCellNames{end+1}=EPmain.window.outCells{iCell};
                    inputcells{end+1}=[];
                end;
                theCell=find(strcmp(EPmain.window.inCells(iCell),EPdata.cellNames));
                if isempty(theCell)
                    disp(['warning: The cell ' EPmain.window.inCells{iCell} ' is missing from the dataset and will be ignored.']);
                else
                    inputcells{find(strcmp(EPmain.window.outCells(iCell),outCellNames))}(end+1)=theCell;
                end;
            end;
        end;
        
        emptyCells=find(cellfun(@isempty,inputcells));
        if ~isempty(emptyCells)
            disp('warning: some outCells were not assigned inCells and are being dropped.')
            inputcells(emptyCells)=[];
            outCellNames(emptyCells)=[];
        end;
        
        subjectSpecs =find(EPmain.window.specSelect);
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Output Measures File',EPdataset.dataset(EPmain.window.dataset).dataName);
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.window.window,'enable','on');
                set(EPmain.handles.window.done,'enable','on');
                set(EPmain.handles.window.channels,'enable','on');
                set(EPmain.handles.window.automatic,'enable','on');
                drawnow
                return
            end;
        end;
        
        [pathstr, name, ext] = fileparts(FileName);
        
        if ~strcmp(ext,'.txt')
            ext=[ext '.txt'];
        end;
        
        if  strcmp(measure,'behavioral')
            if any(strcmp('RT',EPdataset.dataset(EPmain.window.dataset).trialSpecNames))
                ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'RT');
            end;
            if any(strcmp('ACC',EPdataset.dataset(EPmain.window.dataset).trialSpecNames))
                ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'ACC');
            end;
        else
            if ~isempty(EPdata.freqNames) && (EPmain.window.FFTunits ==1)
                ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'real');
                ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'imag');
            else
                ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'normal');
            end;
            
            if EPmain.preferences.window.adds
                if ~isempty(EPchanGrp) && (theChanGrp <= length(EPchanGrp.group))
                    windowAdds(inputcells, outCellNames, false);
                end;
            end;
        end;
        
        disp('done');
        ep('start')
        
    case 'AutoPCA'
        %output the windowed measure to a file for later ANOVA for an entire PCA dataset, setting channels and windows to
        %the peaks.
        
        set(EPmain.handles.window.window,'enable','off');
        set(EPmain.handles.window.done,'enable','off');
        set(EPmain.handles.window.channels,'enable','off');
        set(EPmain.handles.window.automatic,'enable','off');
        set(EPmain.handles.window.FFTunits,'enable','off');
        drawnow
        
        EPdata=ep_stripAdds(ep_loadEPdataset(EPdataset,EPmain.window.dataset));
        
        numChans=length(EPdata.chanNames);
        numPoints=length(EPdata.timeNames);
        numFacs=length(EPdata.facNames);
        numFreqs=length(EPdata.freqNames);
        
        inputcells=[];
        outCellNames=[];
        for iCell=1:length(EPmain.window.inCells)
            if ~isempty(EPmain.window.inCells{iCell}) && ~isempty(EPmain.window.outCells{iCell})
                if ~any(strcmp(EPmain.window.outCells(iCell),outCellNames))
                    outCellNames{end+1}=EPmain.window.outCells{iCell};
                    inputcells{end+1}=[];
                end;
                theCell=find(strcmp(EPmain.window.inCells(iCell),EPdata.cellNames));
                if isempty(theCell)
                    disp(['warning: The cell ' EPmain.window.inCells{iCell} ' is missing from the dataset and will be ignored.']);
                else
                    inputcells{find(strcmp(EPmain.window.outCells(iCell),outCellNames))}(end+1)=theCell;
                end;
            end;
        end;
        
        emptyCells=find(cellfun(@isempty,inputcells));
        if ~isempty(emptyCells)
            disp('warning: some outCells were not assigned inCells and are being dropped.')
            inputcells(emptyCells)=[];
            outCellNames(emptyCells)=[];
        end;
        
        subjectSpecs =find(EPmain.window.specSelect);
                
        disp('Starting AutoPCA.  Window size and adjoining samples option for peak measures ignored.');
        if isfield(EPdata,'facVar')
            factorList=find(EPdata.facVar >= EPmain.window.minFacVar);
            disp(['There are ' num2str(length(factorList)) ' factors that meet the minimum variance criterion of: ' num2str(EPmain.window.minFacVar) '.']);
        else
            factorList=[1:length(EPdata.facNames)];
        end;
        
        if isempty(factorList)
            disp('Procedure aborted as there are no factors to process.');
        else
            [FileName,PathName,FilterIndex] = uiputfile('*.*','PCA ANOVA file root name');
            
            if isnumeric(FileName)
                if FileName == 0
                    msg{1}='No file name specified.';
                    [msg]=ep_errorMsg(msg);
                    set(EPmain.handles.window.window,'enable','on');
                    set(EPmain.handles.window.done,'enable','on');
                    set(EPmain.handles.window.channels,'enable','on');
                    set(EPmain.handles.window.automatic,'enable','on');
                    drawnow
                    return
                end;
            end;
            
            [pathstr, name, ext] = fileparts(FileName);
            
            if ~strcmp(ext,'.txt')
                ext=[ext '.txt'];
            end;
            
            for theFactor=factorList
                gave=ep_expandFacs(EPdata,find(strcmp('EEG',EPdata.chanTypes)),[],[],[],theFactor,[]);
                if ~isempty(EPdata.facVecT)
                    [C peakLatency]=max(abs(EPdata.facVecT(:,theFactor)));
                else
                    [C peakLatency]=max(max(reshape(shiftdim(abs(mean(gave,4)),1),numPoints,[])'));
                end;
                
                startLatency=peakLatency;
                endLatency=peakLatency;
                
                if ~isempty(EPdata.facVecS)
                    [C peakChan]=max(abs(EPdata.facVecS(:,theFactor)));
                else
                    [C peakChan]=max(max(reshape(abs(mean(gave,4)),numChans,[])'));
                end;
                
                if ~isempty(EPdata.facVecF)
                    [C peakFreq]=max(abs(EPdata.facVecF(:,theFactor)));
                else
                    if numFreqs
                        [C peakFreq]=max(max(reshape(shiftdim(abs(mean(gave,4)),5),numFreqs,[])'));
                    else
                        peakFreq=1;
                    end;
                end;
                
                measureNum = get(EPmain.handles.window.measure,'Value');
                switch measureNum
                    case 1
                        measure='mean';
                    case 2
                        measure='minpeak';
                    case 3
                        measure='maxpeak';
                    case 4
                        measure='minlatency';
                    case 5
                        measure='maxlatency';
                    case 6
                        measure='mincentroid';
                        startLatency=1;
                        endLatency=numPoints;
                        disp('Note: AutoPCA will use the entire epoch for the centroid measure, which will not work well if there are multiple ERP components.  For better results, manually window using a smaller window customized for each factor.');
                    case 7
                        measure='maxcentroid';
                        startLatency=1;
                        endLatency=numPoints;
                        disp('Note: AutoPCA will use the entire epoch for the centroid measure, which will not work well if there are multiple ERP components.  For better results, manually window using a smaller window customized for each factor.');
                end;
                
                chanGrp=[];
                chanGrp.areaName(1)=EPdata.chanNames(peakChan);
                chanGrp.areaName{2}='none';
                chanGrp.name='autoPCA';
                chanGrp.channel=zeros(length(EPdata.chanNames),1);
                chanGrp.channel(peakChan)=1;
                
                if ~isreal(EPdata.data)
                    outputdata=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, startLatency, endLatency, measure, [PathName FileName '-' EPdata.facNames{theFactor} ext], theFactor, peakFreq, peakFreq, EPmain.window.FFTunits,EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, 0,'real');
                    outputdata=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, startLatency, endLatency, measure, [PathName FileName '-' EPdata.facNames{theFactor} ext], theFactor, peakFreq, peakFreq, EPmain.window.FFTunits,EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, 0,'imag');
                else
                    outputdata=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, startLatency, endLatency, measure, [PathName FileName '-' EPdata.facNames{theFactor} ext], theFactor, peakFreq, peakFreq, EPmain.window.FFTunits,EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, 0,'normal');
                end;
            end;
            
            if EPmain.preferences.window.adds
                windowAdds(inputcells, outCellNames, true);
            end;
            
            disp('done');
        end;
        
        ep('start')
        
    case 'startANOVA'
        %Set up ANOVA pane of main window.
        
        set(EPmain.handles.hMainWindow,'Name', 'Robust ANOVA');
        
        if isempty(EPmain.anova.data)
            EPmain.anova.data.leftColumn=1;
            EPmain.anova.data.rightColumn=1;
            EPmain.anova.data.columnNames{1}='none';
            EPmain.anova.data.totalLevels=0;
            EPmain.anova.numComps=0;
            
            for i=1:6
                EPmain.anova.data.between(i)=2;
                EPmain.anova.data.factor{i}='';
                EPmain.anova.data.levels{i}='';
                EPmain.anova.data.betweenName{i}='';
            end;
        end;
        
        betweenNames=EPmain.anova.data.columnNames;
        betweenNames{end+1}='none';
        
        uicontrol('Style','text',...
            'String','Data Columns','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 470 100 20]);
        
        EPmain.handles.anova.leftColumn = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.anova.data.columnNames,...
            'Value',EPmain.anova.data.leftColumn,'Position',[20 450 150 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.leftColumn,''Value'');','if tempVar ~=0,EPmain.anova.data.leftColumn=tempVar;end;','if isempty(tempVar),EPmain.anova.data.leftColumn=tempVar;end;','ep(''start'');']);
        
        if ~isfield(EPmain.anova.data,'name')
            set(EPmain.handles.anova.leftColumn,'enable','off');
        end;
        
        EPmain.handles.anova.rightColumn = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.anova.data.columnNames,...
            'Value',EPmain.anova.data.rightColumn,'Position',[20 430 150 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.rightColumn,''Value'');','if tempVar ~=0,EPmain.anova.data.rightColumn=tempVar;end;','if isempty(tempVar),EPmain.anova.data.rightColumn=tempVar;end;','ep(''start'');']);
        
        if ~isfield(EPmain.anova.data,'name')
            set(EPmain.handles.anova.rightColumn,'enable','off');
        end;
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String',sprintf('have: %d',EPmain.anova.data.totalLevels),'HorizontalAlignment','left',...
            'Position',[25 410 50 20]);
        
        if isfield(EPmain.anova.data,'name')
            needLevels=EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1;
            if needLevels == 1
                needLevels=0; %no within factors
            end;
        else
            needLevels=0;
        end;
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String',sprintf('need: %d',needLevels),'HorizontalAlignment','left',...
            'Position',[100 410 50 20]);
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Between Group','HorizontalAlignment','left',...
            'Position',[25 380 100 20]);
        
        for i=1:6
            
            betweenRange=[setdiff([1:length(betweenNames)],[1:EPmain.anova.data.rightColumn,EPmain.anova.data.between(setdiff(1:6,i))]) length(betweenNames)];
            %range of values betweenGroup can take.  Ones already taken are not an option.  "none" is always an option.
            
            if ~ismember(EPmain.anova.data.between(i),betweenRange)
                EPmain.anova.data.between(i)=length(betweenNames);
            end;
            
            EPmain.handles.anova.between(i) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',betweenNames(betweenRange),...
                'Value',find(EPmain.anova.data.between(i)==betweenRange),'Position',[10 380-i*20 130 20],...
                'Callback', @betweenGroup);
            
            if ~isfield(EPmain.anova.data,'name')
                set(EPmain.handles.anova.between(i),'enable','off');
            end;
            
            if (EPmain.anova.data.between(i) == length(betweenNames))
                EPmain.anova.data.betweenName{i}='';
            end;
            
            EPmain.handles.anova.betweenName(i)= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.anova.data.betweenName{i},...
                'Position',[130 380-i*20 50 20],...
                'Callback', @betweenName);
            
            if ~isfield(EPmain.anova.data,'name') || (EPmain.anova.data.between(i) == length(betweenNames))
                set(EPmain.handles.anova.betweenName(i),'enable','off');
            end;
        end;
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Factor','HorizontalAlignment','left',...
            'Position',[20 230 50 20]);
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Levels','HorizontalAlignment','left',...
            'Position',[80 230 50 20]);
        
        uicontrol('Style','frame',...
            'Position',[15 105 180 130]);
        
        for i=1:6
            
            EPmain.handles.anova.factor(i)= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.anova.data.factor{i},...
                'Position',[20 230-i*20 50 20],...
                'Callback', @ANOVATable);
            
            if ~isfield(EPmain.anova.data,'name')
                set(EPmain.handles.anova.factor(i),'enable','off');
            end;
            
            EPmain.handles.anova.levels(i)= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.anova.data.levels{i},...
                'Position',[80 230-i*20 100 20],...
                'Callback', @ANOVATable);
            
            if ~isfield(EPmain.anova.data,'name')
                set(EPmain.handles.anova.levels(i),'enable','off');
            end;
        end;
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Bonferroni','HorizontalAlignment','left',...
            'Position',[25 80 70 20]);
        
        EPmain.handles.anova.numComps= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.anova.numComps,...
            'Position',[25 60 50 20],...
            'Callback', ['global EPmain;','EPmain.anova.numComps=str2num(get(EPmain.handles.anova.numComps,''String''));']);
        
        EPmain.handles.anova.contrast = uicontrol('Style', 'pushbutton', 'String', 'Contrast','FontSize',EPmain.fontsize,...
            'Position', [100 60 80 30], 'Callback', 'ep_contrastANOVA');
        
        EPmain.handles.anova.load = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
            'Position', [20 30 80 30], 'Callback', 'ep(''loadANOVA'')');
        
        EPmain.handles.anova.view = uicontrol('Style', 'pushbutton', 'String', 'View','FontSize',EPmain.fontsize,...
            'Position', [100 30 80 30], 'Callback', 'ep_viewANOVA');
        
        EPmain.handles.anova.run = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
            'Position', [20 0 80 30], 'Callback', 'ep(''runANOVA'')');
        
        if (EPmain.anova.data.totalLevels ~= (EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1)) && ~(EPmain.anova.data.totalLevels ==0 && length(char([EPmain.anova.data.betweenName])) > 0)
            set(EPmain.handles.anova.run,'enable','off');
            set(EPmain.handles.anova.view,'enable','off');
        end;
        
        EPmain.handles.anova.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [100 0 80 30], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'loadANOVA'
        %Load in an initial ANOVA file in order to set up the structure of the analysis.
        
        [ANOVAfile, pathname] = uigetfile('*.txt','ANOVA File');
        
        if ANOVAfile == 0
            return;
        end;
        
        infid=fopen([pathname ANOVAfile],'r');
        
        if infid==-1
            msg{1}='There was a problem with loading the file.';
            [msg]=ep_errorMsg(msg);
            return
        end        
        
        ANOVA = textscan(infid, '%s','delimiter','\n');
        ANOVAdata.name=ANOVA{1}{1};
        ANOVAdata.window=ANOVA{1}{2};
        ANOVAdata.measure=ANOVA{1}{3};
        ANOVAdata.changrp=ANOVA{1}{4};
        ANOVAdata.factor=ANOVA{1}{5};
        
        remain = ANOVA{1}{6};
        numCols=length(strfind(remain,sprintf('\t')))+1;
        while strcmp(remain(end),sprintf('\t'))
            numCols=numCols-1;
            remain=remain(1:end-1);
        end;
        for i = 1:numCols
            [ANOVAdata.cellNames{i}, remain] = strtok(remain,sprintf('\t'));
        end
        remain = ANOVA{1}{7};
        for i = 1:numCols
            [ANOVAdata.areaNames{i}, remain] = strtok(remain,sprintf('\t'));
        end
        
        numSpecs=length(find(strcmp('spec',ANOVAdata.areaNames)));
        numDataCols=numCols-numSpecs;
        ANOVAdata.betweenLvl=[];
        
        for i=1:size(ANOVA{1},1)-7
            remain = ANOVA{1}{i+7};
            for k = 1:numDataCols
                [theNum, remain] = strtok(remain,sprintf('\t'));
                ANOVAdata.data(i,k)=str2num(theNum);
            end
            for k = 1:numSpecs
                [tempVar, remain] = strtok(remain,sprintf('\t'));
                ANOVAdata.betweenLvl{i,k}=tempVar;
            end
        end;
        
        EPmain.anova.data=[];
        
        if ~strcmp(ANOVAdata.name,'behavioral')
            for column=1:numCols
                EPmain.anova.data.columnNames{column}=[ANOVAdata.cellNames{column} '-' ANOVAdata.areaNames{column}];
            end;
        else
            for column=1:numCols
                EPmain.anova.data.columnNames{column}=[ANOVAdata.cellNames{column}];
            end;
        end;
        
        for i=1:6
            EPmain.anova.data.between(i)=numCols+1;
            EPmain.anova.data.factor{i}='';
            EPmain.anova.data.levels{i}='';
            EPmain.anova.data.betweenName{i}='';
        end;
        
        EPmain.anova.data.leftColumn=1;
        
        EPmain.anova.data.rightColumn=numDataCols;
        
        EPmain.anova.data.data=ANOVAdata.data;
        
        EPmain.anova.data.betweenLvl=ANOVAdata.betweenLvl;
        
        EPmain.anova.data.areaNames=ANOVAdata.areaNames;
        
        EPmain.anova.data.name=ANOVAdata.name;
        
        EPmain.anova.data.totalLevels=0;
        
        ep('start');
        
    case 'runANOVA'
        %run the ANOVAs on the data
        
        set(EPmain.handles.anova.contrast,'enable','off');
        set(EPmain.handles.anova.load,'enable','off');
        set(EPmain.handles.anova.view,'enable','off');
        set(EPmain.handles.anova.run,'enable','off');
        set(EPmain.handles.anova.done,'enable','off');
        drawnow
        
        %initial error checking of settings
        for i=1:6
            if xor(isempty(EPmain.anova.data.betweenName{i}),(EPmain.anova.data.between(i)>length(EPmain.anova.data.columnNames)))
                msg{1}='All between-group independent variables need to have a three letter label specified.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.anova.contrast,'enable','on');
                set(EPmain.handles.anova.load,'enable','on');
                set(EPmain.handles.anova.view,'enable','on');
                set(EPmain.handles.anova.run,'enable','on');
                set(EPmain.handles.anova.done,'enable','on');
                return
            end
            
            if ~isempty(EPmain.anova.data.factor{i}) && isempty(EPmain.anova.data.levels{i})
                msg{1}='There is a within-group factor with a name but no levels.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.anova.contrast,'enable','on');
                set(EPmain.handles.anova.load,'enable','on');
                set(EPmain.handles.anova.view,'enable','on');
                set(EPmain.handles.anova.run,'enable','on');
                set(EPmain.handles.anova.done,'enable','on');
                return
            end
        end;
        
        %determine the within group factors and which are electrode factors
        factorNames=cell(0);
        levelNames=cell(0);
        for i=1:6
            if ~isempty(EPmain.anova.data.factor{i})
                factorNames{end+1}=EPmain.anova.data.factor{i};
                levelNames{end+1}=EPmain.anova.data.levels{i};
            end;
        end;
        
        repFactor=1;
        repCount=1;
        lvlCount=1;
        
        if ~strcmp(EPmain.anova.data.name,'behavioral')
            elecFactors=ones(length(factorNames),1);
            for theFactor =length(factorNames):-1:1
                factorArea=cell(length(levelNames{theFactor}),1);
                for theCol=EPmain.anova.data.leftColumn:EPmain.anova.data.rightColumn
                    repCount=repCount+1;
                    if repCount > repFactor
                        repCount=1;
                        lvlCount=lvlCount+1;
                        if lvlCount > length(levelNames{theFactor})
                            lvlCount=1;
                        end;
                    end;
                    factorArea{lvlCount}{end+1}=EPmain.anova.data.areaNames{theCol};
                end;
                repFactor=repFactor*length(levelNames{theFactor});
                for theLevel=1:length(levelNames{theFactor});
                    if length(unique(factorArea{theLevel})) > 1
                        elecFactors(theFactor)=0; %not an electrode factor
                    end;
                end;
            end;
        else
            elecFactors=zeros(length(factorNames),1);
        end;
        
        if isempty(factorNames)
            factorNames{1}='   ';
        end
        if isempty(levelNames)
            levelNames{1}=' ';
        end
        
        %Select input and output files
        
        [outFileName, pathname] = uiputfile('*.*','ANOVA output:');
        if outFileName == 0
            ep('start');
            return %user hit cancel on file requestor
        end;
        outFileName=[pathname outFileName];
        if exist(outFileName,'file')
            delete(outFileName); %user must have clicked "yes" to whether to replace existing file
        end;
        [outPathstr, outFileName, outExt] = fileparts(outFileName);
        if ~strcmp(outExt,'.html')
            outExt='.html';
        end;
        
        outfid=fopen([outPathstr filesep outFileName outExt],'w');
        
        [ANOVAfiles, pathname] = uigetfile('*.txt','Open:','MultiSelect','on');
        activeDirectory=pathname;
        if ~iscell(ANOVAfiles)
            tempVar=ANOVAfiles;
            ANOVAfiles=[];
            ANOVAfiles{1}=tempVar;
        end;
        if ANOVAfiles{1}==0
            msg{1}='No filenames selected. You have to click on a name';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.anova.contrast,'enable','on');
            set(EPmain.handles.anova.load,'enable','on');
            set(EPmain.handles.anova.view,'enable','on');
            set(EPmain.handles.anova.run,'enable','on');
            set(EPmain.handles.anova.done,'enable','on');
            return
        end
        if ~iscell(ANOVAfiles)
            tempVar=ANOVAfiles;
            ANOVAfiles=[];
            ANOVAfiles{1}=tempVar;
        end;
        for theFile=1:size(ANOVAfiles,2)
            ANOVAfiles{theFile}=[activeDirectory ANOVAfiles{theFile}];
        end;
        
        ANOVAfiles=sort(ANOVAfiles);
        
        %start processing each ANOVA file
        
        fprintf(outfid,'<font face="Courier">');
        fprintf(outfid,'<font size="2">');
        disp('Commencing ANOVA run.  This may take some time.  You may monitor its progress by opening the html output file with a browser and periodically reloading it.'); 
        
        for theFile=1:length(ANOVAfiles)
            [pathstr, fileName, ext] = fileparts(ANOVAfiles{theFile});
            infid=fopen([pathname fileName ext],'r');
            if infid == -1
                msg{1}=['There was an error trying to open ' fileName '.'];
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.anova.contrast,'enable','on');
                set(EPmain.handles.anova.load,'enable','on');
                set(EPmain.handles.anova.view,'enable','on');
                set(EPmain.handles.anova.run,'enable','on');
                set(EPmain.handles.anova.done,'enable','on');
                return
            end
            
            fprintf(outfid,'%s </BR>',fileName);
            disp(fileName);
            
            ANOVA = textscan(infid, '%s','delimiter','\r');
            ANOVAdata.name=deblank(ANOVA{1}{1});
            ANOVAdata.window=deblank(ANOVA{1}{2});
            ANOVAdata.measure=deblank(ANOVA{1}{3});
            ANOVAdata.changrp=deblank(ANOVA{1}{4});
            ANOVAdata.factor=deblank(ANOVA{1}{5});
            ANOVAdata.betweenLvl=cell(0);
            remain = ANOVA{1}{6};
            numCols=length(strfind(remain,sprintf('\t')))+1;
            while strcmp(remain(end),sprintf('\t'))
                numCols=numCols-1;
                remain=remain(1:end-1);
            end;
            for i = 1:numCols
                [ANOVAdata.cellNames{i}, remain] = strtok(remain,sprintf('\t'));
            end
            remain = ANOVA{1}{7};
            for i = 1:numCols
                [ANOVAdata.areaNames{i}, remain] = strtok(remain,sprintf('\t'));
            end
            
            numSpecs=length(find(strcmp('spec',ANOVAdata.areaNames)));
            numDataCols=numCols-numSpecs;
            
            for i=1:size(ANOVA{1},1)-7
                remain = ANOVA{1}{i+7};
                for k = 1:numDataCols
                    [theStr, remain] = strtok(remain,sprintf('\t'));
                    theNum=str2num(theStr);
                    if isempty(theNum)
                        theNum=EPmain.preferences.anova.missing;
                    end;
                    ANOVAdata.data(i,k)=theNum;
                end
                for k = 1:numSpecs
                    [ANOVAdata.betweenLvl{i,k}, remain] = strtok(remain,sprintf('\t'));
                end
            end;
            
            if ~strcmp(EPmain.anova.data.name,'behavioral')
                for column=1:numCols
                    ANOVAcellNames{column}=[ANOVAdata.cellNames{column} '-' ANOVAdata.areaNames{column}];
                end;
            else
                for column=1:numCols
                    ANOVAcellNames{column}=[ANOVAdata.cellNames{column}];
                end;
            end;
            
            if (any(~ismember(EPmain.anova.data.columnNames,ANOVAcellNames)) && ~strcmp('autoPCA',ANOVAdata.changrp)) || (length(EPmain.anova.data.columnNames) ~= length(ANOVAcellNames))
                msg{1}=['The column names for ' fileName ' were different from the ANOVA structure that was set up.'];
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.anova.contrast,'enable','on');
                set(EPmain.handles.anova.load,'enable','on');
                set(EPmain.handles.anova.view,'enable','on');
                set(EPmain.handles.anova.run,'enable','on');
                set(EPmain.handles.anova.done,'enable','on');
                return;
            end;
            
            %check each file to make sure the between group levels are fully crossed
            betweenFactors=[];
            totBetweenLevels=1;
            for theFactor=1:6
                if ~isempty(EPmain.anova.data.betweenName{theFactor})
                    theLevels={ANOVAdata.betweenLvl{:,EPmain.anova.data.between(theFactor)-numDataCols}};
                    for i=1:length(theLevels)
                        theLevelsFirstLetter{i}=theLevels{i}(1);
                    end;
                    
                    if length(unique(theLevels)) ~= length(unique(theLevelsFirstLetter))
                        msg{1}=['The between group variable ' EPmain.anova.data.betweenName{theFactor} ' in ' fileName ' has different level labels with the same first letter.'];
                        [msg]=ep_errorMsg(msg);
                        set(EPmain.handles.anova.contrast,'enable','on');
                        set(EPmain.handles.anova.load,'enable','on');
                        set(EPmain.handles.anova.view,'enable','on');
                        set(EPmain.handles.anova.run,'enable','on');
                        set(EPmain.handles.anova.done,'enable','on');
                        return;
                    end;
                    
                    if length(unique({ANOVAdata.betweenLvl{:,EPmain.anova.data.between(theFactor)-numDataCols}})) ==1
                        msg{1}=['The between group variable ' EPmain.anova.data.betweenName{theFactor} ' in ' fileName ' has only a single level.'];
                        [msg]=ep_errorMsg(msg);
                        set(EPmain.handles.anova.contrast,'enable','on');
                        set(EPmain.handles.anova.load,'enable','on');
                        set(EPmain.handles.anova.view,'enable','on');
                        set(EPmain.handles.anova.run,'enable','on');
                        set(EPmain.handles.anova.done,'enable','on');
                        return;
                    end;
                    
                    betweenFactors(end+1)=EPmain.anova.data.between(theFactor);
                    totBetweenLevels=totBetweenLevels*length(unique(theLevels));
                end;
            end;
            
            ANOVAdata.betweenLvl=ANOVAdata.betweenLvl;
            
            if totBetweenLevels> 1
                %only use first letter of between levels
                for i=1:size(ANOVAdata.betweenLvl,1)
                    for j=1:size(ANOVAdata.betweenLvl,2)
                        ANOVAdata.betweenLvl{i,j}=ANOVAdata.betweenLvl{i,j}(1);
                    end;
                end;
                
                [sortedBetween,index] = sortrows(ANOVAdata.betweenLvl(:,betweenFactors-numDataCols));
                ANOVAdata.data=ANOVAdata.data(index,:);
                ANOVAdata.index=index;
                betweenCombos=cellstr(cell2mat(sortedBetween));
                if length(unique(betweenCombos)) ~= totBetweenLevels
                    msg{1}=['The between group factors are incompletely crossed.  Nested designs are not currently supported.'];
                    [msg]=ep_errorMsg(msg);
                    set(EPmain.handles.anova.contrast,'enable','on');
                    set(EPmain.handles.anova.load,'enable','on');
                    set(EPmain.handles.anova.view,'enable','on');
                    set(EPmain.handles.anova.run,'enable','on');
                    set(EPmain.handles.anova.done,'enable','on');
                    return;
                end;
                
                uniqueCombos=unique(betweenCombos);
                subjects=zeros(totBetweenLevels,1);
                for i=1:totBetweenLevels
                    subjects(i)=length(find(strcmp(uniqueCombos{i},betweenCombos)));
                end;
            else
                subjects=size(ANOVAdata.data,1);
                ANOVAdata.index=[1:subjects];
                betweenCombos=cellstr(repmat('gave',subjects,1));
            end;
            
            if any(subjects < 2)
                msg{1}=['Each between group cell needs at least two subjects to be able to estimate error variance.'];
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.anova.contrast,'enable','on');
                set(EPmain.handles.anova.load,'enable','on');
                set(EPmain.handles.anova.view,'enable','on');
                set(EPmain.handles.anova.run,'enable','on');
                set(EPmain.handles.anova.done,'enable','on');
                return;
            end;
            
            if any(subjects < 4)
                disp(['Warning: you have at least one between group cell where there is less than four subjects.  Such an analysis is likely to yield problematic results even if it is statistically significant.']);
            end;

            totDF=length(levelNames{1})-1;
            for i=2:length(levelNames)
                totDF=totDF*(length(levelNames{i})-1);
            end;
            if totDF > subjects-2*floor(EPmain.preferences.anova.trimming*subjects)
                msg{1}=['Product of (number of levels-1) of all within factors cannot exceed the number of participants minus trimming.'];
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.anova.contrast,'enable','on');
                set(EPmain.handles.anova.load,'enable','on');
                set(EPmain.handles.anova.view,'enable','on');
                set(EPmain.handles.anova.run,'enable','on');
                set(EPmain.handles.anova.done,'enable','on');
                return;
            end;
                
            factorGroupNames=cell(0);
            levelGroupNames=cell(0);
            for i=1:6
                if ~isempty(EPmain.anova.data.betweenName{i})
                    factorGroupNames{end+1}=EPmain.anova.data.betweenName{i};
                    levelGroupNames{end+1}=cell2mat(unique(sortedBetween(:,i)'));
                end;
            end;
            
            numComps=EPmain.anova.numComps;
            if numComps==0
                numComps=1;
            end;
            alpha.uncorrected=.05; %threshold for declaring statistical significance.
            alpha.corrected=alpha.uncorrected/numComps; %with Bonferroni correction
            
            Y=ANOVAdata.data(:,EPmain.anova.data.leftColumn:EPmain.anova.data.rightColumn);
            ep_ADF(Y, subjects, 1, EPmain.preferences.anova.trimming, 1, EPmain.preferences.anova.bootstrap, EPmain.preferences.anova.reps, EPmain.preferences.anova.seed, EPmain.preferences.anova.missing, factorNames, levelNames, factorGroupNames, levelGroupNames, elecFactors, alpha, 0, 1, 1, 0, 0, outfid,[],[]);
            
            if EPmain.preferences.anova.adds
                ANOVAAdds(ANOVAdata, Y, subjects, betweenCombos);
            end;
            fclose(infid);
        end;
        fclose(outfid);
        disp('Finished ANOVA run.');
        
        ep('start');
        
    case 'startSave'
        %Set up save pane of main window.
        
        set(EPmain.handles.hMainWindow,'Name', 'Save Data');
        
        uicontrol('Style','text',...
            'String','Save File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 450 100 20]);
        
        EPmain.handles.save.format = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'Value',EPmain.save.format,'Position',[20 420 150 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.save.format,''Value'');','if tempVar ~=0,EPmain.save.format=tempVar;end;','if isempty(tempVar),EPmain.save.format=tempVar;end;','ep(''start'');']);
        
        uicontrol('Style','text',...
            'String','','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 390 70 20]);
        
        uicontrol('Style','text',...
            'String','Single','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[80 390 50 20]);

        uicontrol('Style','text',...
            'String','Combined','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[120 390 50 20]);
        
        uicontrol('Style','text',...
            'String','Channels','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 370 70 20]);
        
        EPmain.handles.save.SGLchan= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.SGLchan,'Position',[80 370 50 20],...
            'Callback', ['global EPmain;','EPmain.save.SGLchan=get(EPmain.handles.save.SGLchan,''Value'');'],'TooltipString','Single channels.');
        
        EPmain.handles.save.REGchan= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.REGchan,'Position',[120 370 50 20],...
            'Callback', ['global EPmain;','EPmain.save.REGchan=get(EPmain.handles.save.REGchan,''Value'');'],'TooltipString','Regional channels.');
                
        uicontrol('Style','text',...
            'String','Cells','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 350 70 20]);
        
        EPmain.handles.save.SGLcell= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.SGLcell,'Position',[80 350 50 20],...
            'Callback', ['global EPmain;','EPmain.save.SGLcell=get(EPmain.handles.save.SGLcell,''Value'');'],'TooltipString','Single cells.');
        
        EPmain.handles.save.CMBcell= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.CMBcell,'Position',[120 350 50 20],...
            'Callback', ['global EPmain;','EPmain.save.CMBcell=get(EPmain.handles.save.CMBcell,''Value'');'],'TooltipString','Combination cells.');
        
        uicontrol('Style','text',...
            'String','Trials','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 330 70 20]);
        
        EPmain.handles.save.RAW= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.RAW,'Position',[80 330 90 20],'FontSize',EPmain.fontsize,...
            'Callback', ['global EPmain;','EPmain.save.RAW=get(EPmain.handles.save.RAW,''Value'');'],'TooltipString','Single trial data.');
        
        uicontrol('Style','text',...
            'String','Averages','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 310 70 20]);
        
        EPmain.handles.save.AVG= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.AVG,'Position',[80 310 50 20],...
            'Callback', ['global EPmain;','EPmain.save.AVG=get(EPmain.handles.save.AVG,''Value'');'],'TooltipString','Subject averages.');
        
        EPmain.handles.save.GAV= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.GAV,'Position',[120 310 50 20],...
            'Callback', ['global EPmain;','EPmain.save.GAV=get(EPmain.handles.save.GAV,''Value'');'],'TooltipString','Grand averages.');
        
        uicontrol('Style','text',...
            'String','Factors','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 290 70 20]);
        
        EPmain.handles.save.SGLfac= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.SGLfac,'Position',[80 290 50 20],...
            'Callback', ['global EPmain;','EPmain.save.SGLfac=get(EPmain.handles.save.SGLfac,''Value'');'],'TooltipString','Single Factors (including two-step factors).');
        
        EPmain.handles.save.CMBfac= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.CMBfac,'Position',[120 290 50 20],...
            'Callback', ['global EPmain;','EPmain.save.CMBfac=get(EPmain.handles.save.CMBfac,''Value'');'],'TooltipString','Combined Factors (grand factors, not two-step factors).');
        
         EPmain.handles.save.batch = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Single','Convert'},...
            'Value',EPmain.save.batch,'Position',[20 250 150 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.save.batch,''Value'');','if tempVar ~=0,EPmain.save.batch=tempVar;end;','if isempty(tempVar),EPmain.save.batch=tempVar;end;','ep(''start'');']);
       
        if EPmain.save.batch == 1
            
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    fileName=EPdataset.dataset(i).dataName;
                    if strcmp(EPdataset.dataset(i).saved,'no')
                        fileName=['*' fileName];
                    end;
                    tableData{i,1}=fileName;
                end;
            else
                tableData=[];
            end;
            
            tableNames{1}='data';
            columnEditable =  false;
            ColumnFormat{1}=[];
            
            EPmain.handles.save.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',@saveData,...
                'ColumnWidth',{300},'Position',[20 60 170 150]);
            
        else
            uicontrol('Style','text',...
                'String','Read File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 230 100 20]);
            
            EPmain.handles.save.readFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',EPmain.fileFormatReadList,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.save.readFormat,''Value'');','if tempVar ~=0,EPmain.save.readFormat=tempVar;end;','if isempty(tempVar),EPmain.save.readFormat=tempVar;end;','ep(''start'');'],...
                'Value',EPmain.save.readFormat,'Position',[20 210 150 20]);
            
            uicontrol('Style','text',...
                'String','File Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 190 100 20]);
            
            EPmain.handles.save.type= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'continuous','single_trial','average','grand_average','factors'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.save.type,''Value'');','if tempVar ~=0,EPmain.save.type=tempVar;end;','if isempty(tempVar),EPmain.save.type=tempVar;end;','ep(''start'');'],...
                'Value',EPmain.save.type,'Position',[20 170 150 20]);
            
            [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.save.readFormat});
            if strcmp(importFormat,'ep_mat')
                set(EPmain.handles.save.type,'enable','off');
            end;
            
            uicontrol('Style','frame',...
                'Position',[20 80 170 90]);
            
            EPmain.handles.save.check= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
                'String','Single Cell Files',...
                'CallBack',['global EPmain;','EPmain.save.check=get(EPmain.handles.save.check,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
                'Value',EPmain.save.check,'Position',[30 145 150 20]);
            
            EPmain.handles.save.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.save.subject,...
                'CallBack',['global EPmain;','EPmain.save.subject=get(EPmain.handles.save.subject,''String'');','ep(''start'');'],...
                'Position',[30 125 50 20],'TooltipString','example 4:6');
            
            EPmain.handles.save.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Subject','HorizontalAlignment','left',...
                'Position',[90 125 50 20]);
            
            if isempty(EPmain.save.subject)
                set(EPmain.handles.save.subjectLabel,'enable','off');
            elseif isempty(str2num(EPmain.save.subject))
                set(EPmain.handles.save.subjectLabel,'enable','off');
            end;
            
            EPmain.handles.save.cell= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.save.cell,...
                'CallBack',['global EPmain;','EPmain.save.cell=get(EPmain.handles.save.cell,''String'');','ep(''start'');'],...
                'Position',[30 105 50 20],'TooltipString','example 7:9');
            
            EPmain.handles.save.cellLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Cell','HorizontalAlignment','left',...
                'Position',[90 105 50 20]);
            
            if isempty(EPmain.save.cell)
                set(EPmain.handles.save.cellLabel,'enable','off');
            elseif isempty(str2num(EPmain.save.cell))
                set(EPmain.handles.save.cellLabel,'enable','off');
            end;
            
            EPmain.handles.save.freq= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.save.freq,...
                'CallBack',['global EPmain;','EPmain.save.freq=get(EPmain.handles.save.freq,''String'');','ep(''start'');'],...
                'Position',[30 85 50 20],'TooltipString','example 10:12');
            
            EPmain.handles.save.freqLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Freq','HorizontalAlignment','left',...
                'Position',[90 85 50 20]);
            
            if isempty(EPmain.save.freq)
                set(EPmain.handles.save.freqLabel,'enable','off');
            elseif isempty(str2num(EPmain.save.freq))
                set(EPmain.handles.save.freqLabel,'enable','off');
            end;
            
            EPmain.handles.save.convert = uicontrol('Style', 'pushbutton', 'String', 'Convert','FontSize',EPmain.fontsize,...
                'Position', [20 40 80 40], 'Callback', @convertFiles);
            
        end;
        
        EPmain.handles.save.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 80 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'startPreferenceMain'
        
        set(EPmain.handles.hMainWindow,'Name', 'Preferences');
        
        EPmain.handles.preferences.hGeneral = uicontrol('Style', 'pushbutton', 'String', 'Files',...
            'Position', [20 450 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceGeneral'';','ep(''start'');']);
        
        EPmain.handles.preferences.hPreprocess = uicontrol('Style', 'pushbutton', 'String', 'Preprocess',...
            'Position', [20 410 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferencePreprocess'';','ep(''start'');']);
        
        EPmain.handles.preferences.hAverage = uicontrol('Style', 'pushbutton', 'String', 'Average',...
            'Position', [20 370 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceAverage'';','ep(''start'');']);
        
        EPmain.handles.preferences.hTransform = uicontrol('Style', 'pushbutton', 'String', 'Transform',...
            'Position', [20 330 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceTransform'';','ep(''start'');']);
        
        EPmain.handles.preferences.hView = uicontrol('Style', 'pushbutton', 'String', 'View',...
            'Position', [20 290 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceView'';','ep(''start'');']);
        
        EPmain.handles.preferences.hPCA = uicontrol('Style', 'pushbutton', 'String', 'PCA',...
            'Position', [20 250 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferencePCA'';','ep(''start'');']);
        
        EPmain.handles.preferences.hWindow = uicontrol('Style', 'pushbutton', 'String', 'Window',...
            'Position', [20 210 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceWindow'';','ep(''start'');']);
        
        EPmain.handles.preferences.hANOVA = uicontrol('Style', 'pushbutton', 'String', 'ANOVA',...
            'Position', [20 170 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceANOVA'';','ep(''start'');']);
        
        EPmain.handles.preferences.done = uicontrol('Style', 'pushbutton', 'String', 'Main',...
            'Position', [20 100 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
        EPmain.handles.preferences.reset = uicontrol('Style', 'pushbutton', 'String', 'Reset',...
            'Position', [20 50 100 40], 'Callback', @resetPrefs);
        
        EPmain.handles.preferences.save = uicontrol('Style', 'pushbutton', 'String', 'Save',...
            'Position', [20 0 100 40], 'Callback', 'ep(''savePrefs'');');
        
    case 'startPreferenceGeneral'
        
        set(EPmain.handles.hMainWindow,'Name', 'Files Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Session Import Format','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.general.sessionImportFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'Value',EPmain.preferences.general.sessionImportFormat,'Position',[20 450 150 20],'TooltipString','Default file format for reading session files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Session Ouput Format','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.general.sessionOutputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'Value',EPmain.preferences.general.sessionOutputFormat,'Position',[20 410 150 20],'TooltipString','Default file format for saving session files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'File Import Format','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.general.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'Value',EPmain.preferences.general.importFormat,'Position',[20 370 150 20],'TooltipString','Default file format for reading average files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'File Ouput Format','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.general.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'Value',EPmain.preferences.general.outputFormat,'Position',[20 330 150 20],'TooltipString','Default file format for saving average files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'First Text Row','FontSize',EPmain.fontsize,...
            'Position',[20 310 90 20]);
        
        EPmain.handles.preferences.general.firstRow= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.firstRow),'FontSize',EPmain.fontsize,...
            'Position',[120 310 70 20],'TooltipString','First row of data (as opposed to header rows) when reading in text files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Last Text Row','FontSize',EPmain.fontsize,...
            'Position',[20 290 90 20]);
        
        EPmain.handles.preferences.general.lastRow= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.lastRow),'FontSize',EPmain.fontsize,...
            'Position',[120 290 70 20],'TooltipString','Last row of data when reading in text files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'First Text Column','FontSize',EPmain.fontsize,...
            'Position',[20 270 90 20]);
        
        EPmain.handles.preferences.general.firstCol= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.firstCol),'FontSize',EPmain.fontsize,...
            'Position',[120 270 70 20],'TooltipString','First column of data when reading in text files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Last Text Col (0=last)','FontSize',EPmain.fontsize,...
            'Position',[20 250 150 20]);
        
        EPmain.handles.preferences.general.lastCol= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.lastCol),'FontSize',EPmain.fontsize,...
            'Position',[20 230 70 20],'TooltipString','Last column of data when reading in text files (0=last).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'SMI suffix','FontSize',EPmain.fontsize,...
            'Position',[20 210 90 20]);
        
        EPmain.handles.preferences.general.SMIsuffix= uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preferences.general.SMIsuffix,'FontSize',EPmain.fontsize,...
            'Position',[120 210 70 20],'TooltipString','Suffix for SMI data file to merge into file being read.  Should have both an underscore and a dot suffix (e.g., _smi.txt)');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Event file suffix','FontSize',EPmain.fontsize,...
            'Position',[20 190 90 20]);
        
        EPmain.handles.preferences.general.specSuffix= uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preferences.general.specSuffix,'FontSize',EPmain.fontsize,...
            'Position',[120 190 70 20],'TooltipString','Suffix for event data file to merge into file being read.  Should have both an underscore and a dot suffix (e.g., _evt.txt)');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Subject spec suffix','FontSize',EPmain.fontsize,...
            'Position',[20 170 90 20]);
        
        EPmain.handles.preferences.general.subjectSpecSuffix= uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preferences.general.subjectSpecSuffix,'FontSize',EPmain.fontsize,...
            'Position',[120 170 70 20],'TooltipString','Suffix for subject spec text file to merge into file being read.  Should have both an underscore and a dot suffix (e.g., _sub.txt)');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Text File Orientation','FontSize',EPmain.fontsize,...
            'Position',[20 150 150 20]);
        
        EPmain.handles.preferences.general.orientation = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Chan Cols','Chan Rows'},...
            'Value',EPmain.preferences.general.orientation,'Position',[20 130 150 20],'TooltipString','When reading in text files, are channels the rows or the columns?');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'EGIS Montage','FontSize',EPmain.fontsize,...
            'Position',[20 110 150 20]);
        
        EPmain.handles.preferences.general.montage = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.preferences.general.montage,'Position',[180 110 150 20],'TooltipString','When saving files in EGIS file format, add montage information to it?');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Rotate electrodes','FontSize',EPmain.fontsize,...
            'Position',[20 90 90 20]);
        
        EPmain.handles.preferences.general.rotateHead = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'None';'90 clockwise';'180 clockwise';'270 clockwise';'flip';'90 clockwise+flip';'180 clockwise+flip';'270 clockwise+flip'},...
            'Value',EPmain.preferences.general.rotateHead,'Position',[110 90 90 20],'TooltipString','When reading in MFF or FIFF file coordinates, need to rotate and/or flip so head facing up?');

        uicontrol('Style','text','HorizontalAlignment','left','String', 'BV header','FontSize',EPmain.fontsize,...
            'Position',[20 70 150 20]);
        
        EPmain.handles.preferences.general.BVheader = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.preferences.general.BVheader,'Position',[180 70 150 20],'TooltipString','When loading BrainVision files, interpret event codes as being encoded EP Toolkit header and TRSP information?');
        
        EPmain.handles.preferences.general.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 30 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.general.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 30], 'Callback', 'ep(''getGeneralPrefs'');');
        
    case 'getGeneralPrefs' %retrieve preference settings from the general preference input fields
        
        tempVar=EPmain.preferences.general;
        
        EPmain.preferences.general.sessionImportFormat=get(EPmain.handles.preferences.general.sessionImportFormat,'Value');
        
        EPmain.preferences.general.sessionOutputFormat=get(EPmain.handles.preferences.general.sessionOutputFormat,'Value');
        
        EPmain.preferences.general.importFormat=get(EPmain.handles.preferences.general.importFormat,'Value');
        
        EPmain.preferences.general.outputFormat=get(EPmain.handles.preferences.general.outputFormat,'Value');
        
        EPmain.preferences.general.firstRow=str2num(get(EPmain.handles.preferences.general.firstRow,'String'));
        
        EPmain.preferences.general.lastRow=str2num(get(EPmain.handles.preferences.general.lastRow,'String'));

        EPmain.preferences.general.firstCol=str2num(get(EPmain.handles.preferences.general.firstCol,'String'));
        
        EPmain.preferences.general.lastCol=str2num(get(EPmain.handles.preferences.general.lastCol,'String'));
        
        EPmain.preferences.general.SMIsuffix=get(EPmain.handles.preferences.general.SMIsuffix,'String');
        
        EPmain.preferences.general.specSuffix=get(EPmain.handles.preferences.general.specSuffix,'String');
        
        EPmain.preferences.general.subjectSpecSuffix=get(EPmain.handles.preferences.general.subjectSpecSuffix,'String');

        EPmain.preferences.general.orientation=get(EPmain.handles.preferences.general.orientation,'Value');
        
        EPmain.preferences.general.montage=get(EPmain.handles.preferences.general.montage,'Value');
        
        EPmain.preferences.general.rotateHead=get(EPmain.handles.preferences.general.rotateHead,'Value');
        
        EPmain.preferences.general.BVheader=get(EPmain.handles.preferences.general.BVheader,'Value');

        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.general=tempVar;
        else
            EPmain.mode='preferenceMain';
            
            EPmain.average.importFormat=EPmain.preferences.general.sessionImportFormat;
            EPmain.average.outputFormat=EPmain.preferences.general.outputFormat;
            EPmain.transform.importFormat=EPmain.preferences.general.importFormat;
            EPmain.transform.outputFormat=EPmain.preferences.general.outputFormat;
            
            ep('start');
        end;
        
    case 'startPreferencePreprocess'
        
        set(EPmain.handles.hMainWindow,'Name', 'Preprocess Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Saturation','FontSize',EPmain.fontsize,...
            'Position',[20 470 100 20]);
        
        EPmain.handles.preferences.preprocess.saturation= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.saturation),'FontSize',EPmain.fontsize,...
            'Position',[130 470 70 20],'TooltipString','Saturation is the voltage level at which the amplifier reaches the maximum number that it is capable of reporting.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Moving Window','FontSize',EPmain.fontsize,...
            'Position',[20 450 100 20]);
        
        EPmain.handles.preferences.preprocess.window= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.window),'FontSize',EPmain.fontsize,...
            'Position',[130 450 70 20],'TooltipString','Moving Window is the number of milliseconds over which the artifact correction routines average the data in a form of low pass filtering.  The larger the number, the less sensitive it is to high frequency spikes.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Chan Min-Max','FontSize',EPmain.fontsize,...
            'Position',[20 430 100 20]);
        
        EPmain.handles.preferences.preprocess.minmax= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.minmax),'FontSize',EPmain.fontsize,...
            'Position',[130 430 70 20],'TooltipString','Chan Min-Max µv is the maximum allowed change in voltage levels for a channel during a trial before it is deemed to be a bad channel for that trial.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', '% Bad Channel','FontSize',EPmain.fontsize,...
            'Position',[20 410 100 20]);
        
        EPmain.handles.preferences.preprocess.badnum= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.badnum),'FontSize',EPmain.fontsize,...
            'Position',[130 410 70 20],'TooltipString','% Bad Channel is the maximum number of channels allowed to be bad in a trial before it is deemed to be a bad trial.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', '# Neigh Chans','FontSize',EPmain.fontsize,...
            'Position',[20 390 100 20]);
        
        EPmain.handles.preferences.preprocess.neighbors= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.neighbors),'FontSize',EPmain.fontsize,...
            'Position',[130 390 70 20],'TooltipString','# Neigh Chans is the number of channels considered to be a neighbor for purposes of the artifact correction algorithms.  The electrode coordinates are then used to figure out which channels are to be used.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Neigh Diff µv','FontSize',EPmain.fontsize,...
            'Position',[20 370 100 20]);
        
        EPmain.handles.preferences.preprocess.maxneighbor= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.maxneighbor),'FontSize',EPmain.fontsize,...
            'Position',[130 370 70 20],'TooltipString','Neigh Diff µv is the maximum voltage difference allowed between a channel and its neighbors before it is deemed to be a bad channel.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bad Chan Corr','FontSize',EPmain.fontsize,...
            'Position',[20 350 100 20]);
        
        EPmain.handles.preferences.preprocess.badchan= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%.3f', EPmain.preferences.preprocess.badchan),'FontSize',EPmain.fontsize,...
            'Position',[130 350 70 20],'TooltipString','Bad Chan Corr is the correlation criterion for determining whether a channel is a globally bad channel over the entire session.  If its best correlation with any neighbor is lower than this criterion then it is judged to be globally bad.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Blink Corr','FontSize',EPmain.fontsize,...
            'Position',[20 330 100 20]);
        
        EPmain.handles.preferences.preprocess.blink= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%.3f', EPmain.preferences.preprocess.blink),'FontSize',EPmain.fontsize,...
            'Position',[130 330 70 20],'TooltipString','Blink Corr is the correlation criterion for determining whether an ICA factor matches the blink template and should therefore be subtracted from the data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', '% Bad Trial','FontSize',EPmain.fontsize,...
            'Position',[20 310 100 20]);
        
        EPmain.handles.preferences.preprocess.badtrials= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.badtrials),'FontSize',EPmain.fontsize,...
            'Position',[130 310 70 20],'TooltipString','% Bad Trial Chan is the maximum number of trials a channel is allowed to be judged bad before it is deemed to be globally bad.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Size of Chunks','FontSize',EPmain.fontsize,...
            'Position',[20 290 100 20]);
        
        EPmain.handles.preferences.preprocess.chunkSize= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.chunkSize),'FontSize',EPmain.fontsize,...
            'Position',[130 290 70 20],'TooltipString','Size of Chunks is the number of time points that are read into each chunk (about 100,000 per GB of available RAM seems to generally work).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Warn Trials/Cell','FontSize',EPmain.fontsize,...
            'Position',[20 270 100 20]);
        
        EPmain.handles.preferences.preprocess.minTrialsPerCell= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.minTrialsPerCell),'FontSize',EPmain.fontsize,...
            'Position',[130 270 70 20],'TooltipString','Warn Trials/Cell is the minimum number of good trials that is considered to be sufficient for a cell.  Any cells dropping below this number will trigger a warning in the artifact correction log.  There is no other effect of this setting.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bad Neighbors','FontSize',EPmain.fontsize,...
            'Position',[20 250 100 20]);
        
        EPmain.handles.preferences.preprocess.noadjacent= uicontrol('Style','checkbox','HorizontalAlignment','left','value', EPmain.preferences.preprocess.noadjacent,'FontSize',EPmain.fontsize,...
            'Position',[130 250 70 20],'TooltipString','Bad Neighbors is an option where if two neighboring channels are marked as being locally bad then the trial is also marked bad (because this typically means that a movement artifact of some sort is present, as opposed to isolated bad channels).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'µv Move Fac','FontSize',EPmain.fontsize,...
            'Position',[20 230 100 20]);
        
        EPmain.handles.preferences.preprocess.trialminmax= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.trialminmax),'FontSize',EPmain.fontsize,...
            'Position',[130 230 70 20],'TooltipString','µv Move Fac is the maximum voltage difference (maximum-minimum) allowed by a factor by the movement artifact correction step.  Factors exceeding this limit are deeemed to reflect artifacts and are subtracted from the data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Move Corr Facs','FontSize',EPmain.fontsize,...
            'Position',[20 210 100 20]);
        
        EPmain.handles.preferences.preprocess.movefacs= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.movefacs),'FontSize',EPmain.fontsize,...
            'Position',[130 210 70 20],'TooltipString','Move Corr Facs is the number of factors to be retained by the movement artifact correction routine.  A larger number results in a more accurate but slower process.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'No Figure','FontSize',EPmain.fontsize,...
            'Position',[20 190 100 20]);
        
        EPmain.handles.preferences.preprocess.noFigure= uicontrol('Style','checkbox','HorizontalAlignment','left','value', EPmain.preferences.preprocess.noFigure,'FontSize',EPmain.fontsize,...
            'Position',[130 190 70 20],'TooltipString','No Figure is an option to not provide a summary figure for the artifact correction process.  While a very useful figure, it requires substantial memory and so dropping it can be helpful when encountering recalcitrant memory problems. ');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'EMG Ratio','FontSize',EPmain.fontsize,...
            'Position',[20 170 100 20]);
        
        EPmain.handles.preferences.preprocess.EMGratio= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.EMGratio),'FontSize',EPmain.fontsize,...
            'Position',[130 170 70 20],'TooltipString','The minimum ratio of signal power to EMG noise to retain during EMG correction.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'EMG Threshold','FontSize',EPmain.fontsize,...
            'Position',[20 150 100 20]);
        
        EPmain.handles.preferences.preprocess.EMGthresh= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.EMGthresh),'FontSize',EPmain.fontsize,...
            'Position',[130 150 70 20],'TooltipString','The Hz threshold considered to be the lower bound of possible EMG frquencies during EMG correction.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'EOG channels','FontSize',EPmain.fontsize,...
            'Position',[20 130 100 20]);
        
        EPmain.handles.preferences.preprocess.EOGchans= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%.3f', EPmain.preferences.preprocess.EOGchans),'FontSize',EPmain.fontsize,...
            'Position',[130 130 70 20],'TooltipString','EOG chans if automatic identification is not working [LUV RUV LLV RLV LH RH].  Enter in -1 for a missing electrode.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Saccade Potential','FontSize',EPmain.fontsize,...
            'Position',[20 110 100 20]);
        
        EPmain.handles.preferences.preprocess.sacPot= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%.3f', EPmain.preferences.preprocess.sacPot),'FontSize',EPmain.fontsize,...
            'Position',[130 110 70 20],'TooltipString','Threshold for saccade potential detection.');
        
        EPmain.handles.preferences.preprocess.fMRI = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'fMRIb OBS','AMRI ICA'},...
            'Value',EPmain.preferences.preprocess.fMRI,'Position',[20 90 150 20],'TooltipString','Algorithm for correcting for fMRI artifacts.');
        
        EPmain.handles.preferences.preprocess.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.preprocess.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getPreprocessPrefs'');');
        
    case 'getPreprocessPrefs' %retrieve preference settings from the preprocessing preference input fields
        
        tempVar=EPmain.preferences.preprocess;
        
        EPmain.preferences.preprocess.saturation=str2num(get(EPmain.handles.preferences.preprocess.saturation,'String'));
        
        EPmain.preferences.preprocess.window=str2num(get(EPmain.handles.preferences.preprocess.window,'String'));
        
        EPmain.preferences.preprocess.minmax=str2num(get(EPmain.handles.preferences.preprocess.minmax,'String'));
        
        EPmain.preferences.preprocess.badnum=str2num(get(EPmain.handles.preferences.preprocess.badnum,'String'));
        
        EPmain.preferences.preprocess.neighbors=str2num(get(EPmain.handles.preferences.preprocess.neighbors,'String'));
        
        EPmain.preferences.preprocess.maxneighbor=str2num(get(EPmain.handles.preferences.preprocess.maxneighbor,'String'));
        
        EPmain.preferences.preprocess.badchan=str2num(get(EPmain.handles.preferences.preprocess.badchan,'String'));
        
        EPmain.preferences.preprocess.blink=str2num(get(EPmain.handles.preferences.preprocess.blink,'String'));
        
        EPmain.preferences.preprocess.chunkSize=str2num(get(EPmain.handles.preferences.preprocess.chunkSize,'String'));
        
        EPmain.preferences.preprocess.minTrialsPerCell=str2num(get(EPmain.handles.preferences.preprocess.minTrialsPerCell,'String'));
        
        EPmain.preferences.preprocess.noadjacent=get(EPmain.handles.preferences.preprocess.noadjacent,'Value');
        
        EPmain.preferences.preprocess.trialminmax=str2num(get(EPmain.handles.preferences.preprocess.trialminmax,'String'));
        
        EPmain.preferences.preprocess.movefacs=str2num(get(EPmain.handles.preferences.preprocess.movefacs,'String'));
        
        EPmain.preferences.preprocess.noFigure=get(EPmain.handles.preferences.preprocess.noFigure,'Value');
        
        EPmain.preferences.preprocess.EMGratio=str2num(get(EPmain.handles.preferences.preprocess.EMGratio,'String'));
        
        EPmain.preferences.preprocess.EMGthresh=str2num(get(EPmain.handles.preferences.preprocess.EMGthresh,'String'));
        
        EPmain.preferences.preprocess.EOGchans=str2num(get(EPmain.handles.preferences.preprocess.EOGchans,'String'));
        
        EPmain.preferences.preprocess.sacPot=str2num(get(EPmain.handles.preferences.preprocess.sacPot,'String'));

        EPmain.preferences.preprocess.fMRI=get(EPmain.handles.preferences.preprocess.fMRI,'Value');

        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.preprocess=tempVar;
        else
            EPmain.mode='preferenceMain';
            ep('start');
        end;
        
    case 'startPreferenceAverage'
        
        set(EPmain.handles.hMainWindow,'Name', 'Average Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Averaging Method','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.average.method = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Mean','Median','Trimmed Mean'},...
            'Value',EPmain.preferences.average.method,'Position',[20 450 150 20],'TooltipString','Central tendency estimator used for averaging procedure.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Trimming Level','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.average.trimLevel = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.average.trimLevel,'Position',[20 410 150 20],'TooltipString','If trimmed mean chosen, proportion of trials trimmed from each tail of the distribution.');
        
        EPmain.handles.preferences.average.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.average.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getAveragePrefs'');');
        
    case 'getAveragePrefs' %retrieve preference settings from the general preference input fields
        
        tempVar=EPmain.preferences.average;
        
        EPmain.preferences.average.method=get(EPmain.handles.preferences.average.method,'Value');
        
        EPmain.preferences.average.trimLevel=str2num(get(EPmain.handles.preferences.average.trimLevel,'String'));
        
        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.average=tempVar;
        else
            EPmain.mode='preferenceMain';
            ep('start');
        end;
        
    case 'startPreferenceTransform'
        
        set(EPmain.handles.hMainWindow,'Name', 'Transform Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.transform.reference = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Average','Traditional','none'},...
            'Value',EPmain.preferences.transform.reference,'Position',[20 450 150 20],'TooltipString','Default type of referencing scheme to apply when transforming data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference Channel 1','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.transform.refChan1 = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.refChan1,'Position',[20 410 150 20],'TooltipString','If traditional reference chosen, default channel to use.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference Channel 2','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.transform.refChan2 = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.refChan2,'Position',[20 370 150 20],'TooltipString','If traditional reference chosen, default second channel to use, if any.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline Start','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.transform.baselineStart = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.baselineStart,'Position',[20 330 150 20],'TooltipString','Default msec start (left side of sample) of period to use for baseline correction when transforming data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline End','FontSize',EPmain.fontsize,...
            'Position',[20 310 150 20]);
        
        EPmain.handles.preferences.transform.baselineEnd = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.baselineEnd,'Position',[20 290 150 20],'TooltipString','Default msec end (right side of sample) of period to use for baseline correction when transforming data.');
        
        EPmain.handles.preferences.transform.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.transform.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getTransformPrefs'');');
        
    case 'getTransformPrefs' %retrieve preference settings from the transform preference input fields
        
        tempVar=EPmain.preferences.transform;
        
        EPmain.preferences.transform.reference=get(EPmain.handles.preferences.transform.reference,'Value');
        
        EPmain.preferences.transform.refChan1=str2num(get(EPmain.handles.preferences.transform.refChan1,'String'));
        
        EPmain.preferences.transform.refChan2=str2num(get(EPmain.handles.preferences.transform.refChan2,'String'));
        
        EPmain.preferences.transform.baselineStart=str2num(get(EPmain.handles.preferences.transform.baselineStart,'String'));
        
        EPmain.preferences.transform.baselineEnd=str2num(get(EPmain.handles.preferences.transform.baselineEnd,'String'));
        
        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.transform=tempVar;
        else
            EPmain.mode='preferenceMain';
            ep('start');
        end;
        
    case 'startPreferencePCA'
        
        set(EPmain.handles.hMainWindow,'Name', 'PCA Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Mode','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.pca.mode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'spatial','temporal'},...
            'Value',EPmain.preferences.pca.mode,'Position',[20 450 150 20],'TooltipString','When conducting PCA, whether default mode is temporal or spatial.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Rotation','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.pca.rotation = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Varimax','Promax','Infomax','Quartimax','Quartimin','Oblimin','CRFE','MINE','IPSC','TIIC','Geomin','MMER','VOMN'},...
            'Value',EPmain.preferences.pca.rotation,'Position',[20 410 150 20],'TooltipString','When conducting PCA, default rotation.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Parameter','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.pca.rotopt = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.pca.rotopt,'Position',[20 370 150 20],'TooltipString','When conducting PCA, default parameter for rotations that have a parameter.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Rel matrix','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.pca.rel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'correlation','covariance'},...
            'Value',EPmain.preferences.pca.rel,'Position',[20 330 150 20],'TooltipString','When conducting PCA, the default relationship matrix.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Loadings','FontSize',EPmain.fontsize,...
            'Position',[20 310 150 20]);
        
        EPmain.handles.preferences.pca.loadings = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'none','Kaiser','covariance','C-M'},...
            'Value',EPmain.preferences.pca.loadings,'Position',[20 290 150 20],'TooltipString','When conducting PCA, the default factor loading weighting scheme.');
        
        EPmain.handles.preferences.pca.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.pca.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getPCAPrefs'');');
        
    case 'getPCAPrefs' %retrieve preference settings from the PCA preference input fields
        
        tempVar=EPmain.preferences.pca;
        
        EPmain.preferences.pca.mode=get(EPmain.handles.preferences.pca.mode,'Value');
        
        EPmain.preferences.pca.rotation=get(EPmain.handles.preferences.pca.rotation,'Value');
        
        EPmain.preferences.pca.rotopt=str2num(get(EPmain.handles.preferences.pca.rotopt,'String'));
        
        EPmain.preferences.pca.rel=get(EPmain.handles.preferences.pca.rel,'Value');
        
        EPmain.preferences.pca.loadings=get(EPmain.handles.preferences.pca.loadings,'Value');
        
        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.pca=tempVar;
        else
            EPmain.mode='preferenceMain';
            ep('start');
        end;
        
    case 'startPreferenceView'
        
        set(EPmain.handles.hMainWindow,'Name', 'View Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Plot positive...','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.view.positive = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'up','down'},...
            'Value',EPmain.preferences.view.positive,'Position',[20 450 150 20],'TooltipString','When plotting waveforms, whether positive is up or down.');
        
        EPmain.handles.preferences.view.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.view.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getViewPrefs'');');
        
    case 'getViewPrefs' %retrieve preference settings from the view preference input fields
        
        tempVar=EPmain.preferences.view;
        
        EPmain.preferences.view.positive=get(EPmain.handles.preferences.view.positive,'Value');
        
        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.view=tempVar;
        else
            EPmain.mode='preferenceMain';
            ep('start');
        end;
        
    case 'startPreferenceWindow'
        
        set(EPmain.handles.hMainWindow,'Name', 'Window Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Min Factor Variance','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.window.minFacVar = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.window.minFacVar,'Position',[20 450 150 20],'TooltipString','When using autoPCA option, the minimum percent of variance for a factor to be included.');
        
        EPmain.handles.preferences.window.adds= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Adds',...
            'Value',EPmain.preferences.window.adds,'Position',[20 430 150 20],'TooltipString','When windowing data, whether to add summary cell and channel waveforms corresponding to the analysis to the dataset.');
        
        EPmain.handles.preferences.window.chanGrp = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Collapse First','Measure First'},...
            'Value',EPmain.preferences.window.chanGrp,'Position',[20 410 150 20],'TooltipString','Collapse channels then measure or vice versa.');
        
        EPmain.handles.preferences.window.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.window.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getWindowPrefs'');');
        
    case 'getWindowPrefs' %retrieve preference settings from the window preference input fields
        
        tempVar=EPmain.preferences.window;
        
        EPmain.preferences.window.minFacVar=str2num(get(EPmain.handles.preferences.window.minFacVar,'String'));
        
        EPmain.preferences.window.adds=get(EPmain.handles.preferences.window.adds,'Value');
        
        EPmain.preferences.window.chanGrp=get(EPmain.handles.preferences.window.chanGrp,'Value');
        
        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.window=tempVar;
        else
            EPmain.mode='preferenceMain';
            ep('start');
        end;
        
    case 'startPreferenceANOVA'
        
        set(EPmain.handles.hMainWindow,'Name', 'ANOVA Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Trimming','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.anova.trimming = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.trimming,'Position',[20 450 150 20],'TooltipString','Proportion of distribution of each cell''s tail to trim during ANOVA.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bootstrap Samples','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.anova.bootstrap = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.bootstrap,'Position',[20 410 150 20],'TooltipString','Number of times to run bootstrap to generate ANOVA distributions.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bootstrap Reps','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.anova.reps = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.reps,'Position',[20 370 150 20],'TooltipString','Number of times to run ANOVA reps to determine p-value variability.  Needs to be an odd number.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Seed','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.anova.seed = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.seed,'Position',[20 330 150 20],'TooltipString','The starting seed number to use for random number generator for bootstrapping.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Missing','FontSize',EPmain.fontsize,...
            'Position',[20 310 150 20]);
        
        EPmain.handles.preferences.anova.missing = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.missing,'Position',[20 290 150 20],'TooltipString','Value used to indicate a missing number.');
        
        EPmain.handles.preferences.anova.adds= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Adds',...
            'Value',EPmain.preferences.anova.adds,'Position',[20 270 150 20],'TooltipString','When conducting ANOVA, whether to add summary subject waveforms corresponding to the analysis to the dataset.');
        
        EPmain.handles.preferences.anova.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        
        EPmain.handles.preferences.anova.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getANOVAPrefs'');');
        
    case 'getANOVAPrefs' %retrieve preference settings from the view preference input fields
        
        tempVar=EPmain.preferences.anova;
        
        EPmain.preferences.anova.trimming=str2num(get(EPmain.handles.preferences.anova.trimming,'String'));
        
        EPmain.preferences.anova.bootstrap=str2num(get(EPmain.handles.preferences.anova.bootstrap,'String'));
        
        EPmain.preferences.anova.reps=str2num(get(EPmain.handles.preferences.anova.reps,'String'));
        
        EPmain.preferences.anova.seed=str2num(get(EPmain.handles.preferences.anova.seed,'String'));
        
        EPmain.preferences.anova.missing=str2num(get(EPmain.handles.preferences.anova.missing,'String'));
        
        EPmain.preferences.anova.adds=get(EPmain.handles.preferences.anova.adds,'Value');
        
        err=checkPrefs;
        
        if err
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.anova=tempVar;
        else
            EPmain.mode='preferenceMain';
            ep('start');
        end;
        
     case 'savePrefs' %save current preference settings
        
%         PathName=userpath;
%         if any(strcmp(PathName(end),{';',':'}))
%             PathName=PathName(1:end-1);
%         end;
        
        prefs=EPmain.preferences;
%         if exist('~/Library/Preferences/EPprefs.mat','file')
%             eval('save ~/Library/Preferences/EPprefs.mat prefs');
%         elseif exist('EPprefs.mat','file')
%             location = which('EPprefs.mat');
%             eval(['save ''' location ''' prefs']);
%         elseif exist('~/Library/Preferences','dir')
%             resetPrefs;
%             prefs=EPmain.preferences;
%             eval('save ~/Library/Preferences/EPprefs.mat prefs');
         if exist([EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat'],'file')
             eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat'' prefs']);
        else
%            [FileName,PathName,FilterIndex] = uiputfile('','Save Preferences File','EPprefs');
            eval(['save ''' EPdataset.EPwork 'EPprefs'' prefs']);
%             if ~strcmp(PathName,path)
%                 addpath(PathName,'-end'); %adds directory with preferences file to Matlab's path
%             end;
        end;
        
        EPmain.mode='main';
        ep('start');
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deleteReadData(src,eventdata)
global EPdataset EPmain

if isempty(eventdata.Indices) %if just deselecting
    return;
end;
theDataset=eventdata.Indices(1);

set(EPmain.handles.read.hTable,'enable','off');
set(EPmain.handles.read.hRead,'enable','off');
set(EPmain.handles.read.hQuitRead,'enable','off');
drawnow

delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(theDataset).dataName '.mat']);
EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],theDataset));
eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset'' EPdataset']);

set(EPmain.handles.read.hTable,'enable','on');
set(EPmain.handles.read.hRead,'enable','on');
set(EPmain.handles.read.hQuitRead,'enable','on');
drawnow

ep('start');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pickEditData(src,eventdata)
global EPoverview EPmain

if isempty(eventdata.Indices) %if just deselecting
    return;
end;

EPoverview.dataset=eventdata.Indices(1);

EPoverview.mode=[];
ep_editData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pickPCAdata(src,eventdata)
global EPmain EPdataset EPscree

if isempty(eventdata.Indices) %if just deselecting
    return;
end;

set(EPmain.handles.pca.hQuitRead,'enable','off');
drawnow

EPscree=[];
screeFigure=findobj('Name', 'ScreeWindow');
if ~isempty(screeFigure)
    close(screeFigure)
end;

theDataset=eventdata.Indices(1);

if EPmain.pca.mode == 4
    if isempty(EPmain.pca.crossVerifyPCA)
        return
    end;
    theDataset=EPmain.pca.targetDatasets(theDataset);
end;

theDataFull=ep_loadEPdataset(EPdataset,theDataset);

theData=ep_stripAdds(theDataFull);

if length(theData.subNames) ~= length(theDataFull.subNames)
    disp(['Stripping off ' num2str(length(theDataFull.subNames)-length(theData.subNames)) ' subject adds.']);
end;

if length(theData.cellNames) ~= length(theDataFull.cellNames)
    disp(['Stripping off ' num2str(length(theDataFull.cellNames)-length(theData.cellNames)) ' cell adds.']);
end;

EEGchans=find(strcmp('EEG',theData.chanTypes));
numChans=length(theData.chanNames);

if length(EEGchans) ~= numChans
    disp(['Stripping off ' num2str(numChans-length(EEGchans)) ' non-EEG channels and regional EEG channels.']);
end;

theData=ep_selectData(theData,{EEGchans,[],[],[],[],[]});

if isempty(theData.data)
    warndlg(['Error: The file had no data left after additions were removed.']);
    set(EPmain.handles.pca.hQuitRead,'enable','on');
    drawnow
    return;
end;

if strcmp(theData.dataType,'continuous') && (length(theData.timeNames) > 1000) && (EPmain.pca.mode == 2)
    warndlg(['Error: Attempting a temporal PCA on continuous data would likely crash the computer due to excessive memory usage.']);
    set(EPmain.handles.pca.hQuitRead,'enable','on');
    drawnow
    return;
end;

PCAname=get(EPmain.handles.pca.name,'String');
if isempty(PCAname)
    PCAname='PCA';
end;
sameName=1;
suffix=0;
PCAnameSuffix=PCAname;
while sameName
    sameName=0;
    for i=1:length(EPdataset.dataset)
        if strcmp(EPdataset.dataset(i).dataName,PCAnameSuffix)
            sameName=1;
        end;
    end;
    if sameName
        suffix=suffix+1;
        PCAnameSuffix=[PCAname '-' num2str(suffix)];
    end;
end;


if EPmain.pca.mode == 4
    crossVerifyPCA=ep_loadEPdataset(EPdataset,EPmain.pca.crossVerifyPCA);
    crossVerifyPCA=ep_stripAdds(crossVerifyPCA);
    
    %check to see if the data has been edited since the PCA
    if isfield(crossVerifyPCA.pca,'PCAmode')
        switch crossVerifyPCA.pca.PCAmode
            case 'spat'
                if size(crossVerifyPCA.pca.FacPat,1) ~= length(crossVerifyPCA.chanNames)
                    warndlg(['Error: The channels of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
            case 'temp'
                if size(crossVerifyPCA.pca.FacPat,1) ~= length(crossVerifyPCA.timeNames)
                    warndlg(['Error: The time points of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
            case 'freq'
                if size(crossVerifyPCA.pca.FacPat,1) ~= length(crossVerifyPCA.freqNames)
                    warndlg(['Error: The frequencies of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
        end;
    end;
    if isfield(crossVerifyPCA.pca,'PCAmode2')
        switch crossVerifyPCA.pca.PCAmode2
            case 'spat'
                if size(crossVerifyPCA.pca.FacPatST,1) ~= length(crossVerifyPCA.chanNames)
                    warndlg(['Error: The channels of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
            case 'temp'
                if size(crossVerifyPCA.pca.FacPatST,1) ~= length(crossVerifyPCA.timeNames)
                    warndlg(['Error: The time points of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
            case 'freq'
                if size(crossVerifyPCA.pca.FacPatST,1) ~= length(crossVerifyPCA.freqNames)
                    warndlg(['Error: The frequencies of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
        end;
    end;
    if isfield(crossVerifyPCA.pca,'PCAmode3')
        switch crossVerifyPCA.pca.PCAmode
            case 'spat'
                if size(crossVerifyPCA.pca.FacPat3,1) ~= length(crossVerifyPCA.chanNames)
                    warndlg(['Error: The channels of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
            case 'temp'
                if size(crossVerifyPCA.pca.FacPat3,1) ~= length(crossVerifyPCA.timeNames)
                    warndlg(['Error: The time points of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
            case 'freq'
                if size(crossVerifyPCA.pca.FacPat3,1) ~= length(crossVerifyPCA.freqNames)
                    warndlg(['Error: The frequencies of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return;
                end;
        end;
    end;
    
    if isfield(crossVerifyPCA.pca,'PCAmode')
        CVpca=crossVerifyPCA.pca.FacCof;
        FactorResults = ep_doPCA(crossVerifyPCA.pca.PCAmode, crossVerifyPCA.pca.ROTATION, crossVerifyPCA.pca.RotOpt, crossVerifyPCA.pca.MAT_TYPE, crossVerifyPCA.pca.numFacs, theData, crossVerifyPCA.pca.LOADING, 'N', [], CVpca);
        if isempty(FactorResults)
            msg{1}='Error: PCA was not successful.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return
        end;
        [theData] = ep_PCAoutput(FactorResults, [], []);
        if isempty(theData)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return
        end;
        theData.dataName=PCAnameSuffix;
    end;
    if isfield(crossVerifyPCA.pca,'PCAmode2')
        CVpca=crossVerifyPCA.pca.FacCofST;
        [FactorResults] = ep_doPCAst(theData.pca, crossVerifyPCA.pca.ROTATION2, crossVerifyPCA.pca.RotOpt2, crossVerifyPCA.pca.MAT_TYPE2, crossVerifyPCA.pca.numFacs2, crossVerifyPCA.pca.LOADING2,crossVerifyPCA.pca.PCAmode2,CVpca);
        if isempty(FactorResults)
            msg{1}='Error: PCA was not successful.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return
        end;
        [theData] = ep_PCAoutput(FactorResults, [], []);
        if isempty(theData)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return
        end;
        theData.dataName=PCAnameSuffix;
    end;
    if isfield(crossVerifyPCA.pca,'PCAmode3')
        CVpca=crossVerifyPCA.pca.FacCof3;
        [FactorResults] = ep_doPCAst(theData.pca, crossVerifyPCA.pca.ROTATION3, crossVerifyPCA.pca.RotOpt3, crossVerifyPCA.pca.MAT_TYPE3, crossVerifyPCA.pca.numFacs3, crossVerifyPCA.pca.LOADING3,crossVerifyPCA.pca.PCAmode3,CVpca);
        if isempty(FactorResults)
            msg{1}='Error: PCA was not successful.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return
        end;
        [theData] = ep_PCAoutput(FactorResults, [], []);
        if isempty(theData)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return
        end;
        theData.dataName=PCAnameSuffix;
    end;
    
    EPdataset=ep_saveEPdataset(EPdataset,theData,length(EPdataset.dataset)+1,'no');
    
else
    modeNum = get(EPmain.handles.pca.mode,'value');
    switch modeNum
        case 1
            PCAmode='spat';
        case 2
            PCAmode='temp';
        case 3
            PCAmode='freq';
    end;
    EPmain.pca.mode=modeNum;
    
    rotationNum = get(EPmain.handles.pca.rotation,'value');
    switch rotationNum
        case 1
            PCArotation='VMAX';
        case 2
            PCArotation='PMAX';
        case 3
            PCArotation='IMAX';
        case 4
            PCArotation='QMAX';
        case 5
            PCArotation='QMIN';
        case 6
            PCArotation='OMIN';
        case 7
            PCArotation='CRFE';
        case 8
            PCArotation='MINE';
        case 9
            PCArotation='IPSC';
        case 10
            PCArotation='TIIC';
        case 11
            PCArotation='GMIN';
        case 12
            PCArotation='MMER';
        case 13
            PCArotation='VOMN';
        case 14
            PCArotation='UNRT';
    end;
    
    if strcmp(PCArotation,'IMAX')
        disp('Using EEGlab function runica to perform Infomax rotation.');
    end;
    
    EPmain.pca.rotation=rotationNum;
    
    EPmain.pca.rotopt=str2num(get(EPmain.handles.pca.rotopt,'String'));
    
    relNum = get(EPmain.handles.pca.rel,'value');
    switch relNum
        case 1
            PCArel='COR';
        case 2
            PCArel='COV';
    end;
    EPmain.pca.rel=relNum;
    
    loadingNum = get(EPmain.handles.pca.loadings,'value');
    switch loadingNum
        case 1
            PCAloading='N';
        case 2
            PCAloading='K';
        case 3
            PCAloading='C';
        case 4
            PCAloading='W';
    end;
    EPmain.pca.loading=loadingNum;
    
    parametric = get(EPmain.handles.pca.parametric,'value');
    EPmain.pca.parametric=parametric;
    
    EPmain.pca.facNum=str2num(get(EPmain.handles.pca.facNum,'String'));
    
    if parametric
        [parametricFile, pathname] = uigetfile('*.*','Parametric File');
        if parametricFile ~= 0
            parametricData=importdata([pathname parametricFile],sprintf('\t'),1);
            if isempty(parametricData)
                msg{1}=['Was unable to retrieve the parametric data.  It needs to be a tab-delimited text file in which the first row are the labels.'];
                [msg]=ep_errorMsg(msg);
                return
            end;
        else
            parametricData.data=[];
            parametricData.colheaders=[];
        end;
    else
        parametricData.data=[];
        parametricData.colheaders=[];
    end;
    
    disp('Starting the PCA.');
    
    if ~isfield(theData.pca,'PCAmode')
        if strcmp(PCAmode,'temp')  && isempty(theData.timeNames)
            warndlg(['Error: The file has no time information with which to perform a temporal PCA.']);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return;
        end;
        
        if strcmp(PCAmode,'spat')  && isempty(theData.chanNames)
            warndlg(['Error: The file has no spatial information with which to perform a spatial PCA.']);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return;
        end;
        
        if strcmp(PCAmode,'freq')  && isempty(theData.freqNames)
            warndlg(['Error: The file has no frequency information with which to perform a frequency PCA.']);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return;
        end;
        
        if any(any(theData.analysis.badTrials)) && strcmp(theData.dataType,'single_trial')
            disp('Conducting PCA despite presence of bad data.  Observations with bad data will be left out from PCA calculations')
            disp('so they will not affect results but doing so will underweight their aspects of the dataset.');
        end;
        
        if any(any(any(theData.analysis.badChans < 0))) && strcmp(theData.dataType,'single_trial')
            disp('Conducting PCA despite presence of bad data.  Observations with bad data will be left out from PCA calculations')
            disp('so they will not affect results but doing so will underweight their aspects of the dataset.');
        end;
        
        if any(any(any(isnan(theData.analysis.badChans)))) && strcmp(theData.dataType,'average')
            disp('Conducting PCA despite presence of bad data.  Observations with bad data will be left out from PCA calculations')
            disp('so they will not affect results but doing so will underweight their aspects of the dataset.');
        end;
        
        if EPmain.pca.facNum==0 %scree test only needs unrotated solution
            disp('Scree: The PCA of the data.');
            FactorResults = ep_doPCA(PCAmode, 'UNRT', EPmain.pca.rotopt, PCArel, 1, theData, PCAloading);
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
            if ~isempty(theData.noise) && any(any(any(any(any(theData.noise)))))
                disp('Scree: The PCA of the noise data.');
                noiseResults = ep_doPCA(PCAmode, 'UNRT', EPmain.pca.rotopt, PCArel, 1, theData.noise, PCAloading);
                if isempty(noiseResults)
                    msg{1}='Error: PCA was not successful.';
                    [msg]=ep_errorMsg(msg);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return
                end;
            else
                noiseResults = [];
            end;
            disp('Scree: The PCA of the random data.');
            randResults = ep_doPCA(PCAmode, 'UNRT', EPmain.pca.rotopt, PCArel, 1, randn(size(theData.data)), PCAloading);
            if isempty(randResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
        else
            FactorResults = ep_doPCA(PCAmode, PCArotation, EPmain.pca.rotopt, PCArel, EPmain.pca.facNum, theData, PCAloading);
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
            
            %jack-knife analysis
            numSubs=length(theData.subNames);
            if numSubs > 1
                disp('Computing Jack-Knife PCA loadings.');
                fprintf('%60s\n',' ' );
                jackLoadings=zeros(size(FactorResults.FacPat,1),size(FactorResults.FacPat,2),numSubs);
                jackSD=zeros(length(FactorResults.varSD),numSubs);
                for theSubject=1:numSubs
                    fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d%s%4d.','Working on subject# ', theSubject, ' of ', numSubs))
                    theDataJK=theData;
                    subjectList=setdiff([1:numSubs],theSubject);
                    theDataJK.data=theDataJK.data(:,:,:,subjectList,:,:,:);
                    theDataJK.subNames=theDataJK.subNames(subjectList);
                    theDataJK.subTypes=theDataJK.subTypes(subjectList);
                    if ~isempty(theDataJK.subjectSpecs)
                        theDataJK.subjectSpecs=theDataJK.subjectSpecs(subjectList,:);
                    end;
                    theDataJK.avgNum=theDataJK.avgNum(subjectList,:);
                    theDataJK.covNum=theDataJK.covNum(subjectList,:);
                    theDataJK.subNum=theDataJK.subNum(subjectList,:);
                    theDataJK.events=theDataJK.events(subjectList,:);
                    theDataJK.analysis.badChans=theDataJK.analysis.badChans(subjectList,:,:);
                    theDataJK.analysis.moveTrial=theDataJK.analysis.moveTrial(subjectList,:);
                    theDataJK.analysis.blinkTrial=theDataJK.analysis.blinkTrial(subjectList,:);
                    theDataJK.analysis.saccadeTrial=theDataJK.analysis.saccadeTrial(subjectList,:);
                    theDataJK.analysis.saccadeOnset=theDataJK.analysis.saccadeOnset(subjectList,:);
                    theDataJK.analysis.badTrials=theDataJK.analysis.badTrials(subjectList,:);
                    FactorResultsJK = ep_doPCA(PCAmode, PCArotation, EPmain.pca.rotopt, PCArel, EPmain.pca.facNum, theDataJK, PCAloading);
                    if ~isempty(FactorResultsJK)
                        jackLoadings(:,:,theSubject)=FactorResultsJK.FacPat;
                        jackSD(:,theSubject)=FactorResultsJK.varSD;
                    end;
                end;
                FactorResults.jack.FacPat=jackLoadings;
                FactorResults.jack.varSD=jackSD;
                fprintf('%60s\n',' ' );
            end;
        end;
    elseif isfield(theData.pca,'PCAmode3')
        msg{1}='These data have already undergone a three-step PCA process.';
        [msg]=ep_errorMsg(msg);
        set(EPmain.handles.pca.hQuitRead,'enable','on');
        return
    elseif isfield(theData.pca,'PCAmode') || isfield(theData.pca,'PCAmode2')
        if isfield(theData.pca,'PCAmode')
            if strcmp(PCAmode,theData.pca.PCAmode)
                if strcmp(PCAmode,'temp')
                    msg{1}='These data have already undergone the temporal PCA step.';
                elseif strcmp(PCAmode,'spat')
                    msg{1}='These data have already undergone the spatial PCA step.';
                elseif strcmp(PCAmode,'freq')
                    msg{1}='These data have already undergone the frequency PCA step.';
                else
                    msg{1}='Error: PCA mode not recognized.';
                end;
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end
        end;
        if isfield(theData.pca,'PCAmode2')
            if strcmp(PCAmode,theData.pca.PCAmode2)
                if strcmp(PCAmode,'temp')
                    msg{1}='These data have already undergone the temporal PCA step.';
                elseif strcmp(PCAmode,'spat')
                    msg{1}='These data have already undergone the spatial PCA step.';
                elseif strcmp(PCAmode,'freq')
                    msg{1}='These data have already undergone the frequency PCA step.';
                else
                    msg{1}='Error: PCA mode not recognized.';
                end;
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end
        end;
        if EPmain.pca.facNum==0
            disp('Scree: The PCA of the data.');
            FactorResults = ep_doPCAst(theData.pca, 'UNRT', EPmain.pca.rotopt, PCArel, 1, PCAloading,PCAmode); %scree test only needs unrotated solution
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
            if isfield(theData.pca,'noise') && ~isempty(theData.pca.noise) && any(any(any(any(any(theData.pca.noise)))))
                %perform two-step PCA of noise data.  Since it is meant to show results when the signal is removed,
                %use the first step scoring coefficients so that they are within the same "window" as the regular data.
                
                disp('Scree: The PCA of the noise data.');
                noiseResults=theData.pca;
                
                numchan=size(noiseResults.noise,1);
                timePoints=size(noiseResults.noise,2);
                numCells=size(noiseResults.noise,3);
                numSubs=size(noiseResults.noise,4);
                
                if strcmp(noiseResults.PCAmode,'temp')
                    data2=zeros(numchan*numCells*numSubs,timePoints);
                    for cell = 1:numCells
                        for sub = 1:numSubs
                            data2((sub-1)*numchan+(cell-1)*numSubs*numchan+1:sub*numchan+(cell-1)*numSubs*numchan,:)=squeeze(noiseResults.noise(:,:,cell,sub));
                        end
                    end
                elseif strcmp(noiseResults.PCAmode,'spat')
                    data2=zeros(timePoints*numCells*numSubs,numchan);
                    for cell = 1:numCells
                        for sub = 1:numSubs
                            data2((sub-1)*timePoints+(cell-1)*numSubs*timePoints+1:sub*timePoints+(cell-1)*numSubs*timePoints,:)=squeeze(noiseResults.noise(:,:,cell,sub))';
                        end
                    end
                elseif strcmp(noiseResults.PCAmode,'asis')
                    data2=noiseResults.noise;
                else
                    error('PCAmode must be set to either temp or spat or asis');
                end;
                
                noiseResults.FacScr=data2*noiseResults.FacCof;
                noiseResults.FacScr=(noiseResults.FacScr)*inv(diag(std(noiseResults.FacScr))); %Standardize factor scores, not mean corrected.
                noiseResults = ep_doPCAst(noiseResults, 'UNRT', EPmain.pca.rotopt, PCArel, 1, PCAloading,PCAmode);
                if isempty(noiseResults)
                    msg{1}='Error: PCA was not successful.';
                    [msg]=ep_errorMsg(msg);
                    set(EPmain.handles.pca.hQuitRead,'enable','on');
                    return
                end;
            else
                noiseResults = [];
            end;
            disp('Scree: The PCA of the random data.');
            randResults = ep_doPCA(theData.pca.PCAmode, theData.pca.ROTATION, theData.pca.RotOpt, theData.pca.MAT_TYPE, theData.pca.numFacs, randn(theData.pca.numchan,theData.pca.timepoints,theData.pca.numCells,theData.pca.numSubs,1,theData.pca.numFreqs), theData.pca.LOADING);
            if isempty(randResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
            randResults = ep_doPCAst(randResults, 'UNRT', EPmain.pca.rotopt, PCArel, 1, PCAloading,PCAmode);
            if isempty(randResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
        else
            [FactorResults] = ep_doPCAst(theData.pca, PCArotation, EPmain.pca.rotopt, PCArel, EPmain.pca.facNum, PCAloading,PCAmode);
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
            
            %jack-knife analysis
            numSubs=length(theData.subNames);
            numCells=length(theData.cellNames);
            numVars=size(FactorResults.FacScr,1)/(numSubs*numCells);
            if numSubs > 1
                disp('Computing Jack-Knife PCA loadings.');
                fprintf('%60s\n',' ' );
                jackLoadingsST=zeros(size(FactorResults.FacPatST,1),size(FactorResults.FacPatST,2),numSubs);
                jackSDST=zeros(size(FactorResults.varSDST,1),size(FactorResults.varSDST,2),numSubs);
                for theSubject=1:numSubs
                    fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d%s%4d.','Working on subject# ', theSubject, ' of ', numSubs))
                    theDataJK=theData.pca;
                    subjectVec=ones(numSubs,1);
                    subjectVec(theSubject)=0;
                    subjectList=kron(ones(numCells,1),subjectVec);
                    subjectList=logical(kron(subjectList,ones(numVars,1)));
                    theDataJK.FacScr=theDataJK.FacScr(subjectList,:);
                    theDataJK.badObs=theDataJK.badObs(subjectList);
                    theDataJK.numSubs=numSubs-1;
                    [FactorResultsJK] = ep_doPCAst(theDataJK, PCArotation, EPmain.pca.rotopt, PCArel, EPmain.pca.facNum, PCAloading,PCAmode);
                    if ~isempty(FactorResultsJK)
                        jackLoadingsST(:,:,theSubject)=FactorResultsJK.FacPatST;
                        jackSDST(:,:,theSubject)=FactorResultsJK.varSDST;
                    end;
                end;
                FactorResults.jack.FacPatST=jackLoadingsST;
                FactorResults.jack.varSDST=jackSDST;
                fprintf('%60s\n',' ' );
            end;
        end;
    else
        msg{1}='These data have already undergone a PCA process but is corrupted in some manner.';
        [msg]=ep_errorMsg(msg);
        set(EPmain.handles.pca.hQuitRead,'enable','on');
        return
    end;
    
    if EPmain.pca.facNum==0
        ep_scree(FactorResults,noiseResults,randResults);
    else
        if parametric
            if size(parametricData.data,1) ~= size(FactorResults.FacScr,1)
                msg{1}='The number of observations in the parameters file does not match the data.  It needs to be a tab-delimited text file in which the first row are the labels.';
                [msg]=ep_errorMsg(msg);
                set(EPmain.handles.pca.hQuitRead,'enable','on');
                return
            end;
        end;
        
        [PCAoutput] = ep_PCAoutput(FactorResults, parametricData.data, parametricData.colheaders);
        if isempty(PCAoutput)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            set(EPmain.handles.pca.hQuitRead,'enable','on');
            return
        end;
        
        PCAoutput.dataName=PCAnameSuffix;
        EPdataset=ep_saveEPdataset(EPdataset,PCAoutput,length(EPdataset.dataset)+1,'no');
    end;
    
end;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resetPrefs(src,eventdata)
%reset the permanent preference settings

global EPmain EP_VER
EPmain.preferences.general.sessionImportFormat=1; %EP format
EPmain.preferences.general.sessionOutputFormat=1; %EP format
EPmain.preferences.general.importFormat=1; %EP format
EPmain.preferences.general.outputFormat=1; %EP format
EPmain.preferences.general.firstRow=1;
EPmain.preferences.general.lastRow=0;
EPmain.preferences.general.firstCol=1;
EPmain.preferences.general.lastCol=0;
EPmain.preferences.general.SMIsuffix='_smi.txt';
EPmain.preferences.general.specSuffix='_evt.txt';
EPmain.preferences.general.subjectSpecSuffix='_sub.txt';
EPmain.preferences.general.orientation=1;
EPmain.preferences.general.montage=1;
EPmain.preferences.general.rotateHead=1;
EPmain.preferences.general.BVheader=1;

EPmain.preferences.preprocess.saturation=1000;
EPmain.preferences.preprocess.window=80;
EPmain.preferences.preprocess.minmax=100;
EPmain.preferences.preprocess.badnum=10;
EPmain.preferences.preprocess.neighbors=6;
EPmain.preferences.preprocess.maxneighbor=30;
EPmain.preferences.preprocess.badchan=.4;
EPmain.preferences.preprocess.blink=.9;
EPmain.preferences.preprocess.badtrials=20;
EPmain.preferences.preprocess.chunkSize=200000;
EPmain.preferences.preprocess.minTrialsPerCell=15;
EPmain.preferences.preprocess.noadjacent=0;
EPmain.preferences.preprocess.trialminmax=200;
EPmain.preferences.preprocess.movefacs=20;
EPmain.preferences.preprocess.noFigure=0;
EPmain.preferences.preprocess.EOGchans=[];
EPmain.preferences.preprocess.sacPot=2;
EPmain.preferences.preprocess.fMRI=1;
EPmain.preferences.preprocess.EMGratio=9;
EPmain.preferences.preprocess.EMGthresh=15;    

EPmain.preferences.average.method=1;
EPmain.preferences.average.trimLevel=.25;

EPmain.preferences.transform.reference=1; %average reference
EPmain.preferences.transform.refChan1=57;
EPmain.preferences.transform.refChan2=100;
EPmain.preferences.transform.baselineStart=-200;
EPmain.preferences.transform.baselineEnd=0;
    
EPmain.preferences.view.positive=1;

EPmain.preferences.pca.mode=2;
EPmain.preferences.pca.rotation=2;
EPmain.preferences.pca.rotopt=3;
EPmain.preferences.pca.rel=2;
EPmain.preferences.pca.loadings=2;

EPmain.preferences.window.minFacVar=.005;
EPmain.preferences.window.adds=1;
EPmain.preferences.window.chanGrp=1;

EPmain.preferences.anova.trimming=.05;
EPmain.preferences.anova.missing=-999;
EPmain.preferences.anova.seed=1000;
EPmain.preferences.anova.bootstrap=4999;
EPmain.preferences.anova.reps=11;
EPmain.preferences.anova.adds=1;

EPmain.preferences.EPver=EP_VER;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err=checkPrefs(src,eventdata) %check the preference settings
global EPmain EP_VER

err=[];

if ~isfield(EPmain.preferences,'EPver')
    disp(['Updating preferences to version: ' EP_VER]);
    err = [err ;100];
end;    

if isfield(EPmain.preferences,'general')
    if ~isfield(EPmain.preferences.general,'sessionImportFormat')
        err = [err ;14];
        disp('The EPmain.preferences.general.sessionImportFormat field is missing.');
    elseif isempty(EPmain.preferences.general.sessionImportFormat)
        err = [err ;14];
        disp('The EPmain.preferences.general.sessionImportFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.sessionImportFormat)
        err = [err ;14];
        disp('The EPmain.preferences.general.sessionImportFormat field is not a number.');
    elseif EPmain.preferences.general.sessionImportFormat < 1 || EPmain.preferences.general.sessionImportFormat > length(EPmain.fileFormatReadList)
        err = [err ;14];
        disp(['The EPmain.preferences.general.sessionImportFormat field is not within the range of 1 to ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end;

    if ~isfield(EPmain.preferences.general,'sessionOutputFormat')
        err = [err ;15];
        disp('The EPmain.preferences.general.sessionOutputFormat field is missing.');
    elseif isempty(EPmain.preferences.general.sessionOutputFormat)
        err = [err ;15];
        disp('The EPmain.preferences.general.sessionOutputFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.sessionOutputFormat)
        err = [err ;15];
        disp('The EPmain.preferences.general.sessionOutputFormat field is not a number.');
    elseif EPmain.preferences.general.sessionOutputFormat < 1 || EPmain.preferences.general.sessionOutputFormat > length(EPmain.fileFormatSaveList)
        err = [err ;15];
        disp(['The EPmain.preferences.general.sessionOutputFormat field is not within the range of 1 to ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end;
    if ~isfield(EPmain.preferences.general,'importFormat')
        err = [err ;19];
        disp('The EPmain.preferences.general.importFormat field is missing.');
    elseif isempty(EPmain.preferences.general.importFormat)
        err = [err ;19];
        disp('The EPmain.preferences.general.importFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.importFormat)
        err = [err ;19];
        disp('The EPmain.preferences.general.importFormat field is not a number.');
    elseif EPmain.preferences.general.importFormat < 1 || EPmain.preferences.general.importFormat > length(EPmain.fileFormatReadList)
        err = [err ;19];
        disp(['The EPmain.preferences.general.importFormat field is not within the range of 1 to ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end;
    
    if ~isfield(EPmain.preferences.general,'outputFormat')
        err = [err ;20];
        disp('The EPmain.preferences.general.outputFormat field is missing.');
    elseif isempty(EPmain.preferences.general.outputFormat)
        err = [err ;20];
        disp('The EPmain.preferences.general.outputFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.outputFormat)
        err = [err ;20];
        disp('The EPmain.preferences.general.outputFormat field is not a number.');
    elseif EPmain.preferences.general.outputFormat < 1 || EPmain.preferences.general.outputFormat > length(EPmain.fileFormatSaveList)
        err = [err ;20];
        disp(['The EPmain.preferences.general.outputFormat field is not within the range of ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end;
    
    if ~isfield(EPmain.preferences.general,'firstRow')
        err = [err ;47];
        disp('The EPmain.preferences.general.firstRow field is missing.');
    elseif isempty(EPmain.preferences.general.firstRow)
        err = [err ;47];
        disp('The EPmain.preferences.general.firstRow field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.firstRow)
        err = [err ;47];
        disp('The EPmain.preferences.general.firstRow field is not a number.');
    elseif EPmain.preferences.general.firstRow < 1
        err = [err ;47];
        disp('The EPmain.preferences.general.firstRow field is less than 1.');
    end;
    
    if ~isfield(EPmain.preferences.general,'lastRow')
        err = [err ;58];
        disp('The EPmain.preferences.general.lastRow field is missing.');
    elseif isempty(EPmain.preferences.general.lastRow)
        err = [err ;58];
        disp('The EPmain.preferences.general.lastRow field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.lastRow)
        err = [err ;58];
        disp('The EPmain.preferences.general.lastRow field is not a number.');
    elseif EPmain.preferences.general.lastRow < 0
        err = [err ;58];
        disp('The EPmain.preferences.general.lastRow field is less than zero.');
    end;
    
    if ~isfield(EPmain.preferences.general,'firstCol')
        err = [err ;48];
        disp('The EPmain.preferences.general.firstCol field is missing.');
    elseif isempty(EPmain.preferences.general.firstCol)
        err = [err ;48];
        disp('The EPmain.preferences.general.firstCol field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.firstCol)
        err = [err ;48];
        disp('The EPmain.preferences.general.firstCol field is not a number.');
    elseif EPmain.preferences.general.firstCol < 1
        err = [err ;48];
        disp('The EPmain.preferences.general.firstCol field is less than 1.');
    end;
    
    if ~isfield(EPmain.preferences.general,'lastCol')
        err = [err ;49];
        disp('The EPmain.preferences.general.lastCol field is missing.');
    elseif isempty(EPmain.preferences.general.lastCol)
        err = [err ;49];
        disp('The EPmain.preferences.general.lastCol field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.lastCol)
        err = [err ;49];
        disp('The EPmain.preferences.general.lastCol field is not a number.');
    elseif EPmain.preferences.general.lastCol < 0
        err = [err ;49];
        disp('The EPmain.preferences.general.lastCol field is less than zero.');
    end;
    
    if ~isfield(EPmain.preferences.general,'SMIsuffix')
        err = [err ;60];
        disp('The EPmain.preferences.general.SMIsuffix field is missing.');
    elseif isempty(EPmain.preferences.general.SMIsuffix)
        err = [err ;60];
        disp('The EPmain.preferences.general.SMIsuffix field is empty.');
    elseif isnumeric(EPmain.preferences.general.SMIsuffix)
        err = [err ;60];
        disp('The EPmain.preferences.general.SMIsuffix field is a number.');
%     elseif ~strcmp(EPmain.preferences.general.SMIsuffix(1),'_')
%         err = [err ;60];
%         disp('The first character of the EPmain.preferences.general.SMIsuffix field is not an underscore.');
    end;
    
    if ~isfield(EPmain.preferences.general,'specSuffix')
        err = [err ;61];
        disp('The EPmain.preferences.general.specSuffix field is missing.');
    elseif isempty(EPmain.preferences.general.specSuffix)
        err = [err ;61];
        disp('The EPmain.preferences.general.specSuffix field is empty.');
    elseif isnumeric(EPmain.preferences.general.specSuffix)
        err = [err ;61];
        disp('The EPmain.preferences.general.specSuffix field is a number.');
    elseif ~strcmp(EPmain.preferences.general.specSuffix(1),'_')
        err = [err ;61];
        disp('The first character of the EPmain.preferences.general.specSuffix field is not an underscore.');
    end;
    
    if ~isfield(EPmain.preferences.general,'subjectSpecSuffix')
        err = [err ;62];
        disp('The EPmain.preferences.general.subjectSpecSuffix field is missing.');
    elseif isempty(EPmain.preferences.general.subjectSpecSuffix)
        err = [err ;62];
        disp('The EPmain.preferences.general.subjectSpecSuffix field is empty.');
    elseif isnumeric(EPmain.preferences.general.subjectSpecSuffix)
        err = [err ;62];
        disp('The EPmain.preferences.general.subjectSpecSuffix field is a number.');
    elseif ~strcmp(EPmain.preferences.general.subjectSpecSuffix(1),'_')
        err = [err ;62];
        disp('The first character of the EPmain.preferences.general.subjectSpecSuffix field is not an underscore.');
    end;

    if ~isfield(EPmain.preferences.general,'orientation')
        err = [err ;50];
        disp('The EPmain.preferences.general.orientation field is missing.');
    elseif isempty(EPmain.preferences.general.orientation)
        err = [err ;50];
        disp('The EPmain.preferences.general.orientation field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.orientation)
        err = [err ;50];
        disp('The EPmain.preferences.general.orientation field is not a number.');
    elseif EPmain.preferences.general.orientation < 1 || EPmain.preferences.general.orientation > 2
        err = [err ;50];
        disp('The EPmain.preferences.general.orientation field is not within the range of 1 to 2.');
    end;
    
    if ~isfield(EPmain.preferences.general,'montage')
        err = [err ;52];
        disp('The EPmain.preferences.general.montage field is missing.');
    elseif isempty(EPmain.preferences.general.montage)
        err = [err ;52];
        disp('The EPmain.preferences.general.montage field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.montage)
        err = [err ;52];
        disp('The EPmain.preferences.general.montage field is not a number.');
    elseif EPmain.preferences.general.montage ~= 0 && EPmain.preferences.general.montage ~= 1
        err = [err ;52];
        disp('The EPmain.preferences.general.montage field does not equal zero or one.');
    end;
    
    if ~isfield(EPmain.preferences.general,'rotateHead')
        err = [err ;53];
        disp('The EPmain.preferences.general.rotateHead field is missing.');
    elseif isempty(EPmain.preferences.general.rotateHead)
        err = [err ;53];
        disp('The EPmain.preferences.general.rotateHead field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.rotateHead)
        err = [err ;53];
        disp('The EPmain.preferences.general.rotateHead field is not a number.');
    elseif EPmain.preferences.general.rotateHead < 1 || EPmain.preferences.general.rotateHead > 8
        err = [err ;53];
        disp('The EPmain.preferences.general.rotateHead field is not within the range of 1 to 8.');
    end;
        
    if ~isfield(EPmain.preferences.general,'BVheader')
        err = [err ;54];
        disp('The EPmain.preferences.general.BVheader field is missing.');
    elseif isempty(EPmain.preferences.general.BVheader)
        err = [err ;54];
        disp('The EPmain.preferences.general.BVheader field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.BVheader)
        err = [err ;54];
        disp('The EPmain.preferences.general.BVheader field is not a number.');
    elseif EPmain.preferences.general.BVheader ~= 0 && EPmain.preferences.general.BVheader ~= 1
        err = [err ;54];
        disp('The EPmain.preferences.general.BVheader field does not equal zero or one.');
    end;
else
    err = [err ;17];
    disp('The EPmain.preferences.general field is missing.');
end

if isfield(EPmain.preferences,'preprocess')
    if ~isfield(EPmain.preferences.preprocess,'saturation')
        err = [err ;1];
        disp('The EPmain.preferences.preprocess.saturation field is missing.');
    elseif isempty(EPmain.preferences.preprocess.saturation)
        err = [err ;1];
        disp('The EPmain.preferences.preprocess.saturation field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.saturation)
        err = [err ;1];
        disp('The EPmain.preferences.preprocess.saturation field is not a number.');
    elseif EPmain.preferences.preprocess.saturation == 0
        err = [err ;1];
        disp('The EPmain.preferences.preprocess.saturation field is zero.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'window')
        err = [err ;2];
        disp('The EPmain.preferences.preprocess.window field is missing.');
    elseif isempty(EPmain.preferences.preprocess.window)
        err = [err ;2];
        disp('The EPmain.preferences.preprocess.window field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.window)
        err = [err ;2];
        disp('The EPmain.preferences.preprocess.window field is not a number.');
    elseif EPmain.preferences.preprocess.window < 1
        err = [err ;2];
        disp('The EPmain.preferences.preprocess.window field is less than one.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'minmax')
        err = [err ;3];
        disp('The EPmain.preferences.preprocess.minmax field is missing.');
    elseif isempty(EPmain.preferences.preprocess.minmax)
        err = [err ;3];
        disp('The EPmain.preferences.preprocess.minmax field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.minmax)
        err = [err ;3];
        disp('The EPmain.preferences.preprocess.minmax field is not a number.');
    elseif EPmain.preferences.preprocess.saturation < 1
        err = [err ;3];
        disp('The EPmain.preferences.preprocess.minmax field is less than one.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'badnum')
        err = [err ;4];
        disp('The EPmain.preferences.preprocess.badnum field is missing.');
    elseif isempty(EPmain.preferences.preprocess.badnum)
        err = [err ;4];
        disp('The EPmain.preferences.preprocess.badnum field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.badnum)
        err = [err ;4];
        disp('The EPmain.preferences.preprocess.badnum field is not a number.');
    elseif EPmain.preferences.preprocess.badnum < 0 || EPmain.preferences.preprocess.badnum > 100
        err = [err ;4];
        disp('The EPmain.preferences.preprocess.badnum field is not between 0 and 100.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'neighbors')
        err = [err ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is missing.');
    elseif isempty(EPmain.preferences.preprocess.neighbors)
        err = [err ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.neighbors)
        err = [err ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is not a number.');
    elseif EPmain.preferences.preprocess.neighbors < 0
        err = [err ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is less than zero.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'maxneighbor')
        err = [err ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is missing.');
    elseif isempty(EPmain.preferences.preprocess.maxneighbor)
        err = [err ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.maxneighbor)
        err = [err ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is not a number.');
    elseif EPmain.preferences.preprocess.maxneighbor < 0
        err = [err ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is less than zero.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'badchan')
        err = [err ;8];
        disp('The EPmain.preferences.preprocess.badchan field is missing.');
    elseif isempty(EPmain.preferences.preprocess.badchan)
        err = [err ;8];
        disp('The EPmain.preferences.preprocess.badchan field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.badchan)
        err = [err ;8];
        disp('The EPmain.preferences.preprocess.badchan field is not a number.');
    elseif EPmain.preferences.preprocess.badchan < 0 || EPmain.preferences.preprocess.badchan > 1
        err = [err ;8];
        disp('The EPmain.preferences.preprocess.badchan field is not between 0 and 1.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'blink')
        err = [err ;9];
        disp('The EPmain.preferences.preprocess.blink field is missing.');
    elseif isempty(EPmain.preferences.preprocess.blink)
        err = [err ;9];
        disp('The EPmain.preferences.preprocess.blink field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.blink)
        err = [err ;9];
        disp('The EPmain.preferences.preprocess.blink field is not a number.');
    elseif EPmain.preferences.preprocess.blink < 0 || EPmain.preferences.preprocess.blink > 1
        err = [err ;9];
        disp('The EPmain.preferences.preprocess.blink field is not between 0 and 1.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'badtrials')
        err = [err ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is missing.');
    elseif isempty(EPmain.preferences.preprocess.badtrials)
        err = [err ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.badtrials)
        err = [err ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is not a number.');
    elseif EPmain.preferences.preprocess.badtrials < 0 || EPmain.preferences.preprocess.badtrials > 100
        err = [err ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is not between 0 and 100.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'chunkSize')
        err = [err ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is missing.');
    elseif isempty(EPmain.preferences.preprocess.chunkSize)
        err = [err ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.chunkSize)
        err = [err ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is not a number.');
    elseif EPmain.preferences.preprocess.chunkSize < 1
        err = [err ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is less than one.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'minTrialsPerCell')
        err = [err ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is missing.');
    elseif isempty(EPmain.preferences.preprocess.minTrialsPerCell)
        err = [err ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.minTrialsPerCell)
        err = [err ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is not a number.');
    elseif EPmain.preferences.preprocess.minTrialsPerCell < 0
        err = [err ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is less than zero.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'noadjacent')
        err = [err ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field is missing.');
    elseif isempty(EPmain.preferences.preprocess.noadjacent)
        err = [err ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.noadjacent)
        err = [err ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field is not a number.');
    elseif EPmain.preferences.preprocess.noadjacent ~= 0 && EPmain.preferences.preprocess.noadjacent ~= 1
        err = [err ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field does not equal zero or one.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'trialminmax')
        err = [err ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is missing.');
    elseif isempty(EPmain.preferences.preprocess.trialminmax)
        err = [err ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.trialminmax)
        err = [err ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is not a number.');
    elseif EPmain.preferences.preprocess.trialminmax < 0
        err = [err ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is less than zero.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'movefacs')
        err = [err ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is missing.');
    elseif isempty(EPmain.preferences.preprocess.movefacs)
        err = [err ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.movefacs)
        err = [err ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is not a number.');
    elseif EPmain.preferences.preprocess.movefacs < 0
        err = [err ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is less than zero.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'noFigure')
        err = [err ;51];
        disp('The EPmain.preferences.preprocess.noFigure field is missing.');
    elseif isempty(EPmain.preferences.preprocess.noFigure)
        err = [err ;51];
        disp('The EPmain.preferences.preprocess.noFigure field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.noFigure)
        err = [err ;51];
        disp('The EPmain.preferences.preprocess.noFigure field is not a number.');
    elseif EPmain.preferences.preprocess.noFigure ~= 0 && EPmain.preferences.preprocess.noFigure ~= 1
        err = [err ;51];
        disp('The EPmain.preferences.preprocess.noFigure field does not equal zero or one.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'EMGratio')
        err = [err ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is missing.');
    elseif isempty(EPmain.preferences.preprocess.EMGratio)
        err = [err ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.EMGratio)
        err = [err ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is not a number.');
    elseif EPmain.preferences.preprocess.EMGratio < 0
        err = [err ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is negative.');
    end;
    
    if ~isfield(EPmain.preferences.preprocess,'EMGthresh')
        err = [err ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is missing.');
    elseif isempty(EPmain.preferences.preprocess.EMGthresh)
        err = [err ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.EMGthresh)
        err = [err ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is not a number.');
    elseif EPmain.preferences.preprocess.EMGthresh < 0
        err = [err ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is less than zero.');
    end;

    if ~isfield(EPmain.preferences.preprocess,'EOGchans')
        err = [err ;55];
        disp('The EPmain.preferences.preprocess.EOGchans field is missing.');
    elseif ~isnumeric(EPmain.preferences.preprocess.EOGchans)
        err = [err ;55];
        disp('The EPmain.preferences.preprocess.EOGchans field is not a number.');
    elseif length(EPmain.preferences.preprocess.EOGchans) ~= 6 && ~isempty(EPmain.preferences.preprocess.EOGchans)
        err = [err ;55];
        disp('The EPmain.preferences.preprocess.EOGchans field does not have six numbers in it.');
    end;

    if ~isfield(EPmain.preferences.preprocess,'sacPot')
        err = [err ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is missing.');
    elseif isempty(EPmain.preferences.preprocess.sacPot)
        err = [err ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.sacPot)
        err = [err ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is not a number.');
    elseif EPmain.preferences.preprocess.sacPot < 0
        err = [err ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is less than zero.');
    end;

    if ~isfield(EPmain.preferences.preprocess,'fMRI')
        err = [err ;58];
        disp('The EPmain.preferences.preprocess.fMRI field is missing.');
    elseif ~isnumeric(EPmain.preferences.preprocess.fMRI)
        err = [err ;58];
        disp('The EPmain.preferences.preprocess.fMRI field is not a number.');
    elseif EPmain.preferences.preprocess.fMRI < 1 || EPmain.preferences.preprocess.fMRI > 2
        err = [err ;58];
        disp('The EPmain.preferences.preprocess.fMRI field is not between 1 and 2.');
    end;

else
    err = [err ;16];
    disp('The EPmain.preferences.preprocess field is missing.');
end

if isfield(EPmain.preferences,'average')
    if ~isfield(EPmain.preferences.average,'method')
        err = [err ;23];
        disp('The EPmain.preferences.average.method field is missing.');
    elseif isempty(EPmain.preferences.average.method)
        err = [err ;23];
        disp('The EPmain.preferences.average.method field is empty.');
    elseif ~isnumeric(EPmain.preferences.average.method)
        err = [err ;23];
        disp('The EPmain.preferences.average.method field is not a number.');
    elseif EPmain.preferences.average.method < 1 && EPmain.preferences.average.method > 3
        err = [err ;23];
        disp('The EPmain.preferences.average.method field is not within the range of 1 to 3.');
    end;
    if ~isfield(EPmain.preferences.average,'trimLevel')
        err = [err ;24];
        disp('The EPmain.preferences.average.trimLevel field is missing.');
    elseif isempty(EPmain.preferences.average.trimLevel)
        err = [err ;24];
        disp('The EPmain.preferences.average.trimLevel field is empty.');
    elseif ~isnumeric(EPmain.preferences.average.trimLevel)
        err = [err ;24];
        disp('The EPmain.preferences.average.trimLevel field is not a number.');
    elseif ~(EPmain.preferences.average.trimLevel < .5 && EPmain.preferences.average.trimLevel > 0)
        err = [err ;24];
        disp('The EPmain.preferences.average.trimLevel field is not between 0 and .5.');
    end;
else
    err = [err ;22];
    disp('The EPmain.preferences.average field is missing.');
end

if isfield(EPmain.preferences,'transform')
    if ~isfield(EPmain.preferences.transform,'reference')
        err = [err ;26];
        disp('The EPmain.preferences.transform.reference field is missing.');
    elseif isempty(EPmain.preferences.transform.reference)
        err = [err ;26];
        disp('The EPmain.preferences.transform.reference field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.reference)
        err = [err ;26];
        disp('The EPmain.preferences.transform.reference field is not a number.');
    elseif EPmain.preferences.transform.reference < 1 && EPmain.transform.transform.reference > 3
        err = [err ;26];
        disp('The EPmain.preferences.transform.reference field is not within the range of 1 to 3.');
    end;
    if ~isfield(EPmain.preferences.transform,'refChan1')
        err = [err ;27];
        disp('The EPmain.preferences.transform.refChan1 field is missing.');
    elseif isempty(EPmain.preferences.transform.refChan1)
        err = [err ;27];
        disp('The EPmain.preferences.transform.refChan1 field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.refChan1)
        err = [err ;27];
        disp('The EPmain.preferences.transform.refChan1 field is not a number.');
    end;
    if ~isfield(EPmain.preferences.transform,'refChan2')
        err = [err ;28];
        disp('The EPmain.preferences.transform.refChan2 field is missing.');
    elseif isempty(EPmain.preferences.transform.refChan2)
        err = [err ;28];
        disp('The EPmain.preferences.transform.refChan2 field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.refChan2)
        err = [err ;28];
        disp('The EPmain.preferences.transform.refChan2 field is not a number.');
    end;
    if ~isfield(EPmain.preferences.transform,'baselineStart')
        err = [err ;29];
        disp('The EPmain.preferences.transform.baselineStart field is missing.');
    elseif isempty(EPmain.preferences.transform.baselineStart)
        err = [err ;29];
        disp('The EPmain.preferences.transform.baselineStart field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.baselineStart)
        err = [err ;29];
        disp('The EPmain.preferences.transform.baselineStart field is not a number.');
    end;
    if ~isfield(EPmain.preferences.transform,'baselineEnd')
        err = [err ;56];
        disp('The EPmain.preferences.transform.baselineEnd field is missing.');
    elseif isempty(EPmain.preferences.transform.baselineEnd)
        err = [err ;56];
        disp('The EPmain.preferences.transform.baselineEnd field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.baselineEnd)
        err = [err ;56];
        disp('The EPmain.preferences.transform.baselineEnd field is not a number.');
    end;
else
    err = [err ;25];
    disp('The EPmain.preferences.transform field is missing.');
end

if isfield(EPmain.preferences,'view')
    if ~isfield(EPmain.preferences.view,'positive')
        err = [err ;37];
        disp('The EPmain.preferences.view.positive field is missing.');
    elseif isempty(EPmain.preferences.view.positive)
        err = [err ;37];
        disp('The EPmain.preferences.view.positive field is empty.');
    elseif ~isnumeric(EPmain.preferences.view.positive)
        err = [err ;37];
        disp('The EPmain.preferences.view.positive field is not a number.');
    elseif EPmain.preferences.view.positive < 1 && EPmain.preferences.view.positive > 2
        err = [err ;37];
        disp('The EPmain.preferences.view.positive field is not within the range of 1 to 2.');
    end;
else
    err = [err ;36];
    disp('The EPmain.preferences.view field is missing.');
end

if isfield(EPmain.preferences,'pca')
    if ~isfield(EPmain.preferences.pca,'mode')
        err = [err ;31];
        disp('The EPmain.preferences.pca.mode field is missing.');
    elseif isempty(EPmain.preferences.pca.mode)
        err = [err ;31];
        disp('The EPmain.preferences.pca.mode field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.mode)
        err = [err ;31];
        disp('The EPmain.preferences.pca.mode field is not a number.');
    elseif EPmain.preferences.pca.mode < 1 && EPmain.preferences.pca.mode > 2
        err = [err ;31];
        disp('The EPmain.preferences.pca.mode field is not within the range of 1 to 2.');
    end;
    if ~isfield(EPmain.preferences.pca,'rotation')
        err = [err ;32];
        disp('The EPmain.preferences.pca.rotation field is missing.');
    elseif isempty(EPmain.preferences.pca.rotation)
        err = [err ;32];
        disp('The EPmain.preferences.pca.rotation field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.rotation)
        err = [err ;32];
        disp('The EPmain.preferences.pca.rotation field is not a number.');
    elseif EPmain.preferences.pca.rotation < 1 && EPmain.preferences.pca.rotation > 13
        err = [err ;32];
        disp('The EPmain.preferences.pca.rotation field is not within the range of 1 to 13.');
    end;
    if ~isfield(EPmain.preferences.pca,'rotopt')
        err = [err ;33];
        disp('The EPmain.preferences.pca.rotopt field is missing.');
    elseif isempty(EPmain.preferences.pca.rotopt)
        err = [err ;33];
        disp('The EPmain.preferences.pca.rotopt field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.rotopt)
        err = [err ;33];
        disp('The EPmain.preferences.pca.rotopt field is not a number.');
    end;
    if ~isfield(EPmain.preferences.pca,'rel')
        err = [err ;34];
        disp('The EPmain.preferences.pca.rel field is missing.');
    elseif isempty(EPmain.preferences.pca.rel)
        err = [err ;34];
        disp('The EPmain.preferences.pca.rel field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.rel)
        err = [err ;34];
        disp('The EPmain.preferences.pca.rel field is not a number.');
    elseif EPmain.preferences.pca.rel < 1 && EPmain.preferences.pca.rel > 2
        err = [err ;34];
        disp('The EPmain.preferences.pca.rel field is not within the range of 1 to 2.');
    end;
    if ~isfield(EPmain.preferences.pca,'loadings')
        err = [err ;35];
        disp('The EPmain.preferences.pca.loadings field is missing.');
    elseif isempty(EPmain.preferences.pca.loadings)
        err = [err ;35];
        disp('The EPmain.preferences.pca.loadings field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.loadings)
        err = [err ;35];
        disp('The EPmain.preferences.pca.loadings field is not a number.');
    elseif EPmain.preferences.pca.loadings < 1 && EPmain.preferences.pca.loadings > 4
        err = [err ;35];
        disp('The EPmain.preferences.pca.loadings field is not within the range of 1 to 4.');
    end;
else
    err = [err ;30];
    disp('The EPmain.preferences.pca field is missing.');
end

if isfield(EPmain.preferences,'window')
    if ~isfield(EPmain.preferences.window,'minFacVar')
        err = [err ;39];
        disp('The EPmain.preferences.window.minFacVar field is missing.');
    elseif isempty(EPmain.preferences.window.minFacVar)
        err = [err ;39];
        disp('The EPmain.preferences.window.minFacVar field is empty.');
    elseif ~isnumeric(EPmain.preferences.window.minFacVar)
        err = [err ;39];
        disp('The EPmain.preferences.window.minFacVar field is not a number.');
    elseif EPmain.preferences.window.minFacVar < 0 && EPmain.preferences.window.minFacVar > 1
        err = [err ;39];
        disp('The EPmain.preferences.window.minFacVar field is not within the range of 0 to 1.');
    end;
    if ~isfield(EPmain.preferences.window,'adds')
        err = [err ;46];
        disp('The EPmain.preferences.window.adds field is missing.');
    elseif isempty(EPmain.preferences.window.adds)
        err = [err ;46];
        disp('The EPmain.preferences.window.adds field is empty.');
    elseif ~isnumeric(EPmain.preferences.window.adds)
        err = [err ;46];
        disp('The EPmain.preferences.window.adds field is not a number.');
    elseif EPmain.preferences.window.adds ~= 0 && EPmain.preferences.window.adds ~=1
        err = [err ;46];
        disp('The EPmain.preferences.window.adds field must be either 0 or 1.');
    end;
    if ~isfield(EPmain.preferences.window,'chanGrp')
        err = [err ;57];
        disp('The EPmain.preferences.window.chanGrp field is missing.');
    elseif isempty(EPmain.preferences.window.chanGrp)
        err = [err ;57];
        disp('The EPmain.preferences.window.chanGrp field is empty.');
    elseif ~isnumeric(EPmain.preferences.window.chanGrp)
        err = [err ;57];
        disp('The EPmain.preferences.window.chanGrp field is not a number.');
    elseif EPmain.preferences.window.chanGrp ~= 1 && EPmain.preferences.window.chanGrp ~=2
        err = [err ;57];
        disp('The EPmain.preferences.window.adds field must be either 1 or 2.');
    end;
else
    err = [err ;38];
    disp('The EPmain.preferences.window field is missing.');
end

if isfield(EPmain.preferences,'anova')
    if ~isfield(EPmain.preferences.anova,'trimming')
        err = [err ;41];
        disp('The EPmain.preferences.anova.trimming field is missing.');
    elseif isempty(EPmain.preferences.anova.trimming)
        err = [err ;41];
        disp('The EPmain.preferences.anova.trimming field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.trimming)
        err = [err ;41];
        disp('The EPmain.preferences.anova.trimming field is not a number.');
    elseif EPmain.preferences.anova.trimming < 0 && EPmain.preferences.anova.trimming > .5
        err = [err ;41];
        disp('The EPmain.preferences.anova.trimming field is not within the range of 0 to .5.');
    end;
    if ~isfield(EPmain.preferences.anova,'bootstrap')
        err = [err ;42];
        disp('The EPmain.preferences.anova.bootstrap field is missing.');
    elseif isempty(EPmain.preferences.anova.bootstrap)
        err = [err ;42];
        disp('The EPmain.preferences.anova.bootstrap field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.bootstrap)
        err = [err ;42];
        disp('The EPmain.preferences.anova.bootstrap field is not a number.');
    end;
    if ~isfield(EPmain.preferences.anova,'reps')
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is missing.');
    elseif isempty(EPmain.preferences.anova.reps)
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.reps)
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is not a number.');
    elseif ceil(EPmain.preferences.anova.reps) ~= EPmain.preferences.anova.reps
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is not an integer.');
    elseif ceil(EPmain.preferences.anova.reps/2) == EPmain.preferences.anova.reps/2
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is not an odd number.');
    end;
    if ~isfield(EPmain.preferences.anova,'seed')
        err = [err ;43];
        disp('The EPmain.preferences.anova.seed field is missing.');
    elseif isempty(EPmain.preferences.anova.seed)
        err = [err ;43];
        disp('The EPmain.preferences.anova.seed field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.seed)
        err = [err ;43];
        disp('The EPmain.preferences.anova.seed field is not a number.');
    end;
    if ~isfield(EPmain.preferences.anova,'missing')
        err = [err ;44];
        disp('The EPmain.preferences.anova.missing field is missing.');
    elseif isempty(EPmain.preferences.anova.missing)
        err = [err ;44];
        disp('The EPmain.preferences.anova.missing field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.missing)
        err = [err ;44];
        disp('The EPmain.preferences.anova.missing field is not a number.');
    end;
    if ~isfield(EPmain.preferences.anova,'adds')
        err = [err ;45];
        disp('The EPmain.preferences.anova.adds field is missing.');
    elseif isempty(EPmain.preferences.anova.adds)
        err = [err ;45];
        disp('The EPmain.preferences.anova.adds field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.adds)
        err = [err ;45];
        disp('The EPmain.preferences.anova.adds field is not a number.');
    elseif EPmain.preferences.anova.adds ~= 0 && EPmain.preferences.anova.adds ~=1
        err = [err ;45];
        disp('The EPmain.preferences.anova.adds field must be either 0 or 1.');
    end;
else
    err = [err ;40];
    disp('The EPmain.preferences.anova field is missing.');
end

if err
    beep
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function splash(src,eventdata) %present splash page

global EPmain

hSPlash=figure('Name', 'About EP Toolkit', 'NumberTitle', 'off', 'MenuBar', 'none');

try
    EPver=ver('EP_Toolkit');
catch
    EPver.Version='unavailable'; %workaround for bug in earlier version of Matlab
end;

copyright=['ERP PCA Toolkit version: ' EPver.Version ' '...
    'Copyright (C) 1999-2018  Joseph Dien. ' ...
    sprintf('\n')...
    'This program is free software: you can redistribute it and/or modify '...
    'it under the terms of the GNU General Public License as published by '...
    'the Free Software Foundation, either version 3 of the License, or '...
    '(at your option) any later version. '...
    'This program is distributed in the hope that it will be useful, '...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of '...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE (especially medical).  See the '...
    'GNU General Public License for more details. '...
    'You should have received a copy of the GNU General Public License '...
    'along with this program.  If not, see <http://www.gnu.org/licenses/>. '...
    sprintf('\n')...
    'Note also that Dr. Dien does not have grant funding for this software '...
    'project and cannot therefore guarantee the same level of support '...
    'available for software projects that have full-time staff allocated '...
    'to them.  Nonetheless, Dr. Dien will endeavor to help users with '...
    'bug reports and questions to the best of his ability, given his time constraints. '...
    'He will be happy to help with implementing support for new file types. '...
    'Others are welcome to contribute code to this software project. '...
    'See files in the documentation directory for more information.'...
    sprintf('\n')...
    '***Again, this software is intended solely for scientific research and should never be used for medical purposes.***'];

if ispc
    fontsize=8;
else
    fontsize=10;
end;

uicontrol('Style','text','FontSize',fontsize,...
    'String',copyright,...
    'HorizontalAlignment','left','Position',[20 50 400 350]);

if isempty(EPmain)
EPmain.handles.splash.agree = uicontrol('Style', 'pushbutton', 'String', 'Agree','FontSize',fontsize,...
    'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''agree'';','ep']);

EPmain.handles.splash.disagree = uicontrol('Style', 'pushbutton', 'String', 'Disagree','FontSize',fontsize,...
    'Position', [120 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''disagree'';','ep']);
else
    uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',fontsize,...
    'Position', [120 0 100 40], 'Callback', ['close(''About EP Toolkit'');','ep(''start'');']);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSampAdapt(src,eventdata) %change period around peak that is averaged together for peak measures

global EPmain EPdataset

windowLength=EPmain.window.sampEnd-EPmain.window.sampStart+1;
Fs=EPdataset.dataset(EPmain.window.dataset).Fs;

if src==EPmain.handles.window.sampAdapt
    sampAdapt=str2num(get(EPmain.handles.window.sampAdapt,'String'));
else
    sampAdapt=round((str2num(get(EPmain.handles.window.msAdapt,'String'))/(1000/Fs)));
end;

if (sampAdapt*2+1) > windowLength
    sampAdapt=floor((windowLength-1)/2);
    disp(['For current window size, max peak width is: ' num2str(sampAdapt) '.']);
end;

if sampAdapt < 0
    sampAdapt=0;
    disp('Peak width cannot be negative.');
end;

set(EPmain.handles.window.sampAdapt,'String',sampAdapt);
set(EPmain.handles.window.msAdapt,'String',round(sampAdapt*(1000/Fs)));

EPmain.window.sampAdapt=sampAdapt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSampStart(src,eventdata) %change start of window

global EPmain EPdataset

numPoints=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
Fs=EPdataset.dataset(EPmain.window.dataset).Fs;
baseline=EPdataset.dataset(EPmain.window.dataset).baseline;

if src==EPmain.handles.window.sampStart
    sampStart=str2num(get(EPmain.handles.window.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.window.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.window.msStart,'String'))/(1000/Fs))+baseline+1);
    sampEnd=round((str2num(get(EPmain.handles.window.msEnd,'String'))/(1000/Fs))+baseline);
end;

if sampStart > numPoints
    sampStart=numPoints;
end;

if sampStart < 1
    sampStart=1;
end;

if sampStart > sampEnd
    sampEnd = sampStart;
end;

set(EPmain.handles.window.sampStart,'String',sampStart);
set(EPmain.handles.window.sampEnd,'String',sampEnd);
set(EPmain.handles.window.msStart,'String',(sampStart-baseline-1)*(1000/Fs));
set(EPmain.handles.window.msEnd,'String',(sampEnd-baseline)*(1000/Fs));

EPmain.window.sampStart=sampStart;
EPmain.window.sampEnd=sampEnd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSampEnd(src,eventdata) %change end of window

global EPmain EPdataset

numPoints=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
Fs=EPdataset.dataset(EPmain.window.dataset).Fs;
baseline=EPdataset.dataset(EPmain.window.dataset).baseline;

if src==EPmain.handles.window.sampEnd
    sampStart=str2num(get(EPmain.handles.window.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.window.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.window.msStart,'String'))/(1000/Fs))+baseline+1);
    sampEnd=round((str2num(get(EPmain.handles.window.msEnd,'String'))/(1000/Fs))+baseline);
end;

if sampEnd > numPoints
    sampEnd=numPoints;
end;

if sampEnd < 1
    sampEnd=1;
end;

if sampStart > sampEnd
    sampStart = sampEnd;
end;

set(EPmain.handles.window.sampStart,'String',sampStart);
set(EPmain.handles.window.sampEnd,'String',sampEnd);
set(EPmain.handles.window.msStart,'String',(sampStart-baseline-1)*(1000/Fs));
set(EPmain.handles.window.msEnd,'String',(sampEnd-baseline)*(1000/Fs));

EPmain.window.sampStart=sampStart;
EPmain.window.sampEnd=sampEnd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowHzStart(src,eventdata) %change start of Hz window

global EPmain EPdataset

HzStart=str2num(get(EPmain.handles.window.binStart,'String'));
HzEnd=str2num(get(EPmain.handles.window.binEnd,'String'));
numFreqs=length(EPdataset.dataset(EPmain.window.dataset).freqNames);

if HzStart > numFreqs
    HzStart=numFreqs;
end;

if HzStart < 1
    HzStart=1;
end;

if HzStart > HzEnd
    HzEnd = HzStart;
end;

set(EPmain.handles.window.binStart,'String',HzStart);
set(EPmain.handles.window.binEnd,'String',HzEnd);

EPmain.window.HzStart=HzStart;
EPmain.window.HzEnd=HzEnd;

set(EPmain.handles.window.HzStart,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzStart)*10)/10));
set(EPmain.handles.window.HzEnd,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzEnd)*10)/10));
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowHzEnd(src,eventdata) %change end of Hz window

global EPmain EPdataset

HzStart=str2num(get(EPmain.handles.window.binStart,'String'));
HzEnd=str2num(get(EPmain.handles.window.binEnd,'String'));
numFreqs=length(EPdataset.dataset(EPmain.window.dataset).freqNames);

if HzEnd > numFreqs
    HzEnd=numFreqs;
end;

if HzEnd < 1
    HzEnd=1;
end;

if HzStart > HzEnd
    HzStart = HzEnd;
end;

set(EPmain.handles.window.binStart,'String',HzStart);
set(EPmain.handles.window.binEnd,'String',HzEnd);

EPmain.window.HzStart=HzStart;
EPmain.window.HzEnd=HzEnd;

set(EPmain.handles.window.HzStart,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzStart)*10)/10));
set(EPmain.handles.window.HzEnd,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzEnd)*10)/10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSpecsTable(src,eventdata) %change specs table data in Window pane

global EPmain

specTable=get(EPmain.handles.window.specsTable,'Data');
EPmain.window.specSelect=cell2mat(specTable(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowCellTable(src,eventdata) %change cell table data in Window pane

global EPmain

cellTable=get(EPmain.handles.window.cellTable,'Data');
% if any(any(cellfun(@isempty,cellTable())))
%     disp('Cell names cannot be empty.');
%     ep('start');
%     return
% end;

EPmain.window.outCells=cellTable(:,1);
EPmain.window.inCells=cellTable(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveData(src,eventdata) 
%respond to click on save table

global EPmain EPdataset

if isempty(eventdata.Indices) %if just deselecting
    return;
end;
whichData=eventdata.Indices(1);

set(EPmain.handles.save.hTable,'enable','off');
set(EPmain.handles.save.done,'enable','off');
EPdata=ep_loadEPdataset(EPdataset,whichData);

[outFileName, pathname] = uiputfile('*.*','Save:',EPdataset.dataset(whichData).dataName);
if outFileName == 0
    set(EPmain.handles.save.done,'enable','on');
    set(EPmain.handles.save.hTable,'enable','on');
    drawnow
    return %user hit cancel on file requestor
end;

if exist([pathname outFileName],'file')
    delete([pathname outFileName]); %user must have clicked "yes" to whether to replace existing file
end;

if ~strcmp(EPdataset.dataset(whichData).dataName,outFileName) %if changed the dataset name
    button = questdlg('Did you want to change the name of the dataset in the working set as well?','Rename dataset?','Yes','No','Yes');
    if strcmp(button,'Yes')
        EPdata.dataName=outFileName;
        delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(whichData).dataName '.mat']);
        EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],whichData));
        EPdataset=ep_saveEPdataset(EPdataset,EPdata,whichData,'no');
    end;
end;

EPmain.convertMode=0;

saveFlag=saveFile(EPdata,[pathname outFileName]);

if saveFlag
    EPdataset.dataset(whichData).saved='yes';
end;

set(EPmain.handles.save.hTable,'enable','on');
set(EPmain.handles.save.done,'enable','on');
drawnow

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveFlag=saveFile(EPdata, fileName) 
%save a file

global EPmain EPdataset

saveFlag=0;

set(EPmain.handles.save.done,'enable','off');
drawnow

formatNum = get(EPmain.handles.save.format,'value');
switch formatNum
    case 1
        format='ep_mat';
    case 2
        format='egi_egis';
    case 3
        format='egi_sbin';
    case 4
        format='eeglab_set';      
    case 5
        format='text';
    case 6
        format='eeglab_erp';
    case 7
        format='neuromag_fif';
end;

if strcmp(format,'egi_egis') && ~strcmp(EPdata.dataType,'single_trial')
    format='egi_egia';
end;

if strcmp(format,'egi_egia') && EPdata.baseline ==0
    choice = questdlg('EGIS average file with a zero baseline will cause problems for NetStation.  Save anyway?', ...
        'Options', ...
        'Yes','Cancel','Yes');
    % Handle response
    switch choice
        case 'Yes'
        case 'Cancel'
            set(EPmain.handles.save.done,'enable','on');
            drawnow
            return
    end
end;

EPmain.save.format=formatNum;

if ~EPmain.save.SGLchan && ~EPmain.save.REGchan
    EPmain.save.SGLchan=1; %at least some type of channel needs to be kept
    disp('at least some type of channel needs to be kept so keeping single channels.');
end;

if ~EPmain.save.SGLcell && ~EPmain.save.CMBcell
    EPmain.save.SGLcell=1; %at least some type of cell needs to be kept
    disp('at least some type of cell needs to be kept so keeping single channels.');
end;

if ~EPmain.save.RAW && ~EPmain.save.AVG && ~EPmain.save.GAV
    EPmain.save.AVG=1; %at least some type of subject needs to be kept
    disp('at least some type of subject needs to be kept so keeping single channels.');
end;

if ~EPmain.save.SGLfac && ~EPmain.save.CMBfac
    EPmain.save.SGLfac=1; %at least some type of factor needs to be kept
    disp('at least some type of factor needs to be kept (if any) so keeping single channels.');
end;

adds=[];
if EPmain.save.SGLchan
    adds{end+1}='SGLchan';
end;
if EPmain.save.REGchan
    adds{end+1}='REGchan';
end;
if EPmain.save.SGLcell
    adds{end+1}='SGLcell';
end;
if EPmain.save.CMBcell
    adds{end+1}='CMBcell';
end;
if EPmain.save.RAW
    adds{end+1}='RAW';
end;
if EPmain.save.AVG
    if any(strcmp(format,{'egi_egis','egi_egia','egi_sbin','eeglab_set','text','neuromag_fif'})) && ~isempty(EPdata.facNames) &&  (length(EPdata.subNames)>1)
        choice = questdlg('Averaged Data option will result in an additional file for each subject.  Do anyway?', ...
            'Options', ...
            'Yes','No','Cancel','Yes');
        % Handle response
        switch choice
            case 'Yes'
                adds{end+1}='AVG';
            case 'No'
            case 'Cancel'
                set(EPmain.handles.save.done,'enable','on');
                drawnow
                return
        end
    else
        adds{end+1}='AVG';
    end;
end;
if EPmain.save.GAV
    adds{end+1}='GAV';
end;
if EPmain.save.SGLfac
    adds{end+1}='SGLfac';
end;
if EPmain.save.CMBfac
    adds{end+1}='CMBfac';
end;

prefs={};
if ~EPmain.preferences.general.montage
    prefs{end+1}='no montage';
end;

err=ep_writeData(EPdata,fileName,EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,format,adds, prefs);

if ~err
    saveFlag=1;
end;

if ~EPmain.convertMode
    set(EPmain.handles.save.done,'enable','on');
    drawnow
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowAdds(inputcells, outCellNames, autoPCA)
%add regional channels and combined cells to the active dataset that correspond to the windowing procedure.

global EPmain EPdataset EPchanGrp

EPdata=ep_loadEPdataset(EPdataset,EPmain.window.dataset);
numChans=length(EPdata.chanNames);
numPoints=length(EPdata.timeNames);
numCells=length(EPdata.cellNames);
numSubs=length(EPdata.subNames);
numFacs=length(EPdata.facNames);
numFreqs=length(EPdata.freqNames);
if numFacs==0
    numFacs=1;
end;
if ~isempty(EPdata.facData)
    numCMBfacs=size(EPdata.facData,5);
else
    numCMBfacs=0;
end;
numSGLfacs=numFacs-numCMBfacs;

newAreas=[];
tempData=[];
if ~autoPCA
    %add regional channels if any have new names
    areaList = setdiff(unique(EPchanGrp.group(EPmain.window.chanGrp).channel),EPchanGrp.numAreas+1); %identify area numbers that were used
    newAreas=find(~ismember(EPchanGrp.group(EPmain.window.chanGrp).areaName(areaList),EPdata.chanNames));
    if ~isempty(newAreas)
        disp(['Adding windowed regional channels to ' EPdataset.dataset(EPmain.window.dataset).dataName]);
        for theArea=1:length(newAreas)
            whichArea=areaList(newAreas(theArea));
            chanList=find(EPchanGrp.group(EPmain.window.chanGrp).channel == whichArea);
            tempData=ep_combineData(EPdata,'channels',chanList,ones(1,length(chanList)),EPchanGrp.group(EPmain.window.chanGrp).areaName{whichArea});
            if ~isempty(tempData)
                EPdata=tempData;
            end;
        end;
    end;
end;

%add combined output cells if any have new names

newCells=find(~ismember(outCellNames,EPdata.cellNames));
if ~isempty(newCells)
    disp(['Adding windowed combined cells to ' EPdataset.dataset(EPmain.window.dataset).dataName]);
    for theOutCell=1:length(newCells)
        cellList=inputcells{newCells(theOutCell)};
        tempData=ep_combineData(EPdata,'cells',cellList,ones(1,length(cellList)),outCellNames{newCells(theOutCell)});
        if ~isempty(tempData)
            EPdata=tempData;
        end;
    end;
end;

if ~isempty(newAreas) || ~isempty(newCells)
    [err]=ep_checkEPfile(EPdata);
    
    if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
        warndlg('The work directory cannot be found.')
        return
    end;
    
    delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(EPmain.window.dataset).dataName '.mat']);
    EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],EPmain.window.dataset));
    EPdataset=ep_saveEPdataset(EPdataset,EPdata,length(EPdataset.dataset)+1,'no');
    EPdataset.dataset=[EPdataset.dataset(1:EPmain.window.dataset-1) EPdataset.dataset(end) EPdataset.dataset(EPmain.window.dataset:end-1)];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ANOVAAdds(ANOVAdata, Y, subjects, betweenCombos);
%add new grand averages for each between group, such that the subjects that are trimmed correspond to the ones that will be dropped in the ANOVA.
%These will only be computed for channel areas and cells that are present in the ANOVA output.

global EPmain EPdataset

%search through active datasets to see if one has a name matching that in the ANOVA file
if isempty(EPdataset.dataset)
    return
end;
whichData=find(strcmp(ANOVAdata.name,{EPdataset.dataset.dataName}));
if isempty(whichData)
    return
end;

EPdata=ep_loadEPdataset(EPdataset,whichData);

if ~any(ismember(ANOVAdata.cellNames,EPdata.cellNames))
    disp('No cells in the ANOVA are present in the dataset so corresponding grand averages will not be added.');
    return
end;

if ~any(ismember(ANOVAdata.areaNames,EPdata.chanNames))
    disp('No channel regions in the ANOVA are present in the dataset so corresponding grand averages will not be added.');
    return
end;

numChans=length(EPdata.chanNames);
numPoints=length(EPdata.timeNames);
numCells=length(EPdata.cellNames);
numSubs=length(EPdata.subNames);
numFacs=length(EPdata.facNames);
numFreqs=length(EPdata.freqNames);
numRels=length(EPdata.relNames);
if numFacs==0
    numFacs=1;
end;
if ~isempty(EPdata.facData)
    numCMBfacs=size(EPdata.facData,5);
else
    numCMBfacs=0;
end;
numSGLfacs=numFacs-numCMBfacs;

theFactor=str2num(strtok(ANOVAdata.factor,'Factor:'));
if isempty(theFactor)
    theFactor =1;
end;

if strcmp(ANOVAdata.changrp,'autoPCA')
    windowName='autoPCA';
else
    windowName=ANOVAdata.window;
end;

subCounter=0;
for theBetween = 1:length(subjects)
    subList=[1+subCounter:subCounter+subjects(theBetween)];
    betweenName=[betweenCombos{1+subCounter} '_' windowName  '_' ANOVAdata.measure];
    
    if ~any(find(strcmp(betweenName,EPdata.subNames))) %add between group if not already present
        disp(['Adding trimmed grand average ' betweenName ' to ' EPdataset.dataset(whichData).dataName]);
        [EPdata]=ep_combineData(EPdata,'subjects',[1:numSubs],zeros(numSubs,1),betweenName,0);
        if isempty(EPdata)
            disp('Aborting addition to data.')
            return
        end;
        numSubs=numSubs+1;
    end;
    theSub=find(strcmp(betweenName,EPdata.subNames));
    
    for theColumn=EPmain.anova.data.leftColumn:EPmain.anova.data.rightColumn
        theCell=find(strcmp(ANOVAdata.cellNames{theColumn},EPdata.cellNames));
        chanList=find(strcmp(ANOVAdata.areaNames{theColumn},EPdata.chanNames));
        if ~isempty(EPdata.facVecS)
            chanList=1;
        else
            if strcmp(ANOVAdata.changrp,'autoPCA')
                chanList=[1:numChans];
            end;
        end;
        if ~isempty(theCell) && ~isempty(chanList) %if a corresponding cell and channel are present in the data then add the trimmed grand ave
            disp(['Adding ' EPdata.cellNames{theCell} ' waveform corresponding to trimmed ANOVA sample to ' EPdataset.dataset(whichData).dataName]);
            [B,index] = sort(ANOVAdata.data(ANOVAdata.index(subList),theColumn));
            numTrim=floor(EPmain.preferences.anova.trimming*length(subList));
            index=index(1+numTrim:length(subList)-numTrim); %determine which subjects should go into the grand average for this area/cell
            
            combineData=ANOVAdata.index(subList(index));
            combineWeights=ones(1,length(combineData));
            theOkaySubs=combineData(find(EPdata.avgNum(combineData,theCell) >= 0));
            
            for iChan=1:length(chanList)
                theChan=chanList(iChan);
                
                if strcmp(EPdata.dataType,'average')
                    goodSubList=theOkaySubs(find(~isnan(EPdata.analysis.badChans(theOkaySubs,theCell,theChan))));
                else
                    goodSubList=theOkaySubs(find(EPdata.analysis.badChans(theOkaySubs,theCell,theChan) >= 0));
                end;
                if ~isempty(goodSubList)
                    totalWeight=sum(abs(combineWeights(find(ismember(combineData,goodSubList)))));
                    newCellData=zeros(1,max(numPoints,1),1,1,1,max(numFreqs,1),max(numRels,1));
                    newCellNoise=zeros(1,max(numPoints,1),1,1,1);
                    newCellStd=zeros(1,max(numPoints,1),1,1,1,max(numFreqs,1));
                    newCellStdCM=zeros(1,max(numPoints,1),1,1,1,max(numFreqs,1));
                    newFacData=zeros(1,max(numPoints,1),1,1,numCMBfacs,max(numFreqs,1),max(numRels,1));
                    if ~isempty(EPdata.facVecT)
                        newCellData=newCellData(:,1,:,:,:,:);
                    end;
                    if ~isempty(EPdata.facVecF)
                        newCellData=newCellData(:,:,:,:,:,1);
                    end;
                    if ~isempty(EPdata.facVecS)
                        newCellData=newCellData(1,:,:,:,:,:);
                    end;
                    
                    for sub=1:length(goodSubList)
                        thisSub=goodSubList(sub);
                        theWeight=combineWeights(find(ismember(goodSubList(sub),combineData)));
                        newCellData(1,:,1,1,1,:,:)=newCellData(1,:,1,1,1,:)+(theWeight/totalWeight)*EPdata.data(theChan,:,theCell,thisSub,theFactor,:,:);
                        if ~isempty(EPdata.noise)
                            newCellNoise(1,:,1,1,1,:)=newCellNoise(1,:,1,1,1,:)+(theWeight/totalWeight)*EPdata.noise(theChan,:,theCell,thisSub,theFactor,:);
                        end;
                        EPdata.analysis.badChans(theSub,theCell,theChan)= EPdata.analysis.badChans(theSub,theCell,theChan)+EPdata.analysis.badChans(thisSub,theCell,theChan);
                    end;
                    EPdata.data(theChan,:,theCell,theSub,theFactor,:,:)=newCellData;
                    if ~isempty(EPdata.noise)
                        EPdata.noise(theChan,:,theCell,theSub,theFactor)=newCellNoise;
                    end;
                    if ~isempty(EPdata.std)
                        EPdata.std(theChan,:,theCell,theSub,theFactor,:)=newCellStd;
                    end;
                    if ~isempty(EPdata.stdCM)
                        EPdata.stdCM(theChan,:,theCell,theSub,theFactor,:)=newCellStdCM;
                    end;
                    EPdata.avgNum(theSub,theCell)=sum(EPdata.avgNum(theOkaySubs,theCell));
                    EPdata.covNum(theSub,theCell)=sum(EPdata.covNum(theOkaySubs,theCell));
                    EPdata.subNum(theSub,theCell)=sum(EPdata.subNum(theOkaySubs,theCell));
                    EPdata.analysis.blinkTrial(theSub,theCell)= sum(EPdata.analysis.blinkTrial(theOkaySubs,theCell));
                    EPdata.analysis.saccadeTrial(theSub,theCell)= sum(EPdata.analysis.saccadeTrial(theOkaySubs,theCell));
                    EPdata.analysis.saccadeOnset(theSub,theCell)= sum(EPdata.analysis.saccadeOnset(theOkaySubs,theCell));
                    EPdata.analysis.moveTrial(theSub,theCell)= sum(EPdata.analysis.moveTrial(theOkaySubs,theCell));
                    EPdata.analysis.badTrials(theSub,theCell)= sum(EPdata.analysis.badTrials(theOkaySubs,theCell));
                    
                else %none of the subjects have good data for this channel
                    if strcmp(EPdata.dataType,'average')
                        EPdata.analysis.badChans(theSub,theCell,theChan)=NaN;
                    else
                        EPdata.analysis.badChans(theSub,theCell,theChan)=-1;
                    end;
                    EPdata.avgNum(theSub,theCell)=-1;
                    EPdata.covNum(theSub,theCell)=-1;
                    EPdata.subNum(theSub,theCell)=-1;
                    EPdata.analysis.blinkTrial(theSub,theCell)= 0;
                    EPdata.analysis.saccadeTrial(theSub,theCell)= 0;
                    EPdata.analysis.saccadeOnset(theSub,theCell)= 0;
                    EPdata.analysis.moveTrial(theSub,theCell)= 0;
                    EPdata.analysis.badTrials(theSub,theCell)= 0;
                    EPdata.analysis.badChans(theSub,theCell,:)= 0;
                end;
            end;
        end;
    end;
    subCounter=subCounter+subjects(theBetween);
end;

[err]=ep_checkEPfile(EPdata);

if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
    warndlg('The work directory cannot be found.')
    return
end;

delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(whichData).dataName '.mat']);
EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],whichData));
EPdataset=ep_saveEPdataset(EPdataset,EPdata,length(EPdataset.dataset)+1,'no');
EPdataset.dataset=[EPdataset.dataset(1:whichData-1) EPdataset.dataset(end) EPdataset.dataset(whichData:end-1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ANOVATable(src,eventdata)
%process entry into ANOVA factor table.

global EPmain

theFactor=find(src == EPmain.handles.anova.factor);
theLevel=find(src == EPmain.handles.anova.levels);
if theFactor
    theEntry=get(EPmain.handles.anova.factor(theFactor),'string');
    if length(theEntry)==3
        EPmain.anova.data.factor{theFactor}=theEntry;
    else
        EPmain.anova.data.factor{theFactor}=[];
    end;
end;

if theLevel
    theEntry=get(EPmain.handles.anova.levels(theLevel),'string');
    EPmain.anova.data.levels{theLevel}=theEntry;
end;

ANOVAlevels=0;
for theFactor=1:6
    if ~isempty(EPmain.anova.data.levels{theFactor}) && ~isempty(EPmain.anova.data.factor{theFactor})
        if ANOVAlevels==0
            ANOVAlevels=1;
        end;
        ANOVAlevels=ANOVAlevels*length(EPmain.anova.data.levels{theFactor});
    end;
end;

EPmain.anova.data.totalLevels=ANOVAlevels;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function betweenName(src,eventdata)
%process entry of between factor name into ANOVA pane

global EPmain

theFactor=find(src == EPmain.handles.anova.betweenName);
if theFactor
    theEntry=get(EPmain.handles.anova.betweenName(theFactor),'string');
    if length(theEntry)==3
        EPmain.anova.data.betweenName{theFactor}=theEntry;
    else
        EPmain.anova.data.betweenName{theFactor}='';
    end;
end;

ep('start');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function betweenGroup(src,eventdata)
%process selection of between group level column in ANOVA pane

global EPmain

theFactor=find(src == EPmain.handles.anova.between);
if ~isempty(theFactor)
    tempVar=get(EPmain.handles.anova.between(theFactor),'value');
    if tempVar ~=0
        betweenNames=EPmain.anova.data.columnNames;
        betweenNames{end+1}='none';
        betweenRange=[setdiff([1:length(betweenNames)],[1:EPmain.anova.data.rightColumn,EPmain.anova.data.between(setdiff(1:6,theFactor))]) length(betweenNames)];
        EPmain.anova.data.between(theFactor)=betweenRange(tempVar);
    end;
    if isempty(tempVar)
        betweenNames=EPmain.anova.data.columnNames;
        betweenNames{end+1}='none';
        betweenRange=[setdiff([1:length(betweenNames)],[1:EPmain.anova.data.rightColumn,EPmain.anova.data.between(setdiff(1:6,theFactor))]) length(betweenNames)];
        EPmain.anova.data.between(theFactor)=betweenRange(tempVar);
    end;
end;
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function quit(src,eventdata)
%quit from EP Toolkit

global EPdataset

clearWorkingSet

close('EP Toolkit');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeWindowDataset(src,eventdata)
%change the dataset in the window function

global EPmain EPdata EPdataset

tempVar=get(EPmain.handles.window.dataset,'Value');
if tempVar ~=0
    EPmain.window.dataset=EPmain.window.aveData(tempVar);
end;
if isempty(tempVar)
    EPmain.window.dataset=EPmain.window.aveData(tempVar);
end;
EPmain.window.datasetName=EPdataset.dataset(EPmain.window.dataset).dataName;

if isempty(EPdataset.dataset(EPmain.window.dataset).timeNames)
    EPmain.window.measure=1;
end;

if strcmp(EPdataset.dataset(EPmain.window.dataset).dataType,'single_trial')
    [u i]=unique(EPdataset.dataset(EPmain.window.dataset).cellNames,'first');
    EPmain.window.inCells=EPdataset.dataset(EPmain.window.dataset).cellNames(sort(i));
else
    EPmain.window.inCells=EPdataset.dataset(EPmain.window.dataset).cellNames;
end;

EPmain.window.outCells=EPmain.window.inCells;

EPmain.window.specSelect=repmat(false,length(EPdataset.dataset(EPmain.window.dataset).subjectSpecNames),1);

EPmain.window.sampStart=1;
EPmain.window.sampEnd=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
EPmain.window.HzStart=1;
EPmain.window.HzEnd=length(EPdataset.dataset(EPmain.window.dataset).freqNames);

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkViewSettings(src,eventdata)
%check the settings in the view pane to see if the Waves control should be enabled.

global EPmain EPdataset

set(EPmain.handles.view.waves,'enable','on');

bottomVolt=str2num(get(EPmain.handles.view.bottomVolt,'string'));
topVolt=str2num(get(EPmain.handles.view.topVolt,'string'));
startSamp=str2num(get(EPmain.handles.view.startSamp,'string'));
endSamp=str2num(get(EPmain.handles.view.endSamp,'string'));
startHz=str2num(get(EPmain.handles.view.startHz,'string'));
endHz=str2num(get(EPmain.handles.view.endHz,'string'));
marker1=str2num(get(EPmain.handles.view.marker1,'string'));
marker2=str2num(get(EPmain.handles.view.marker2,'string'));

%convert Hz values to non-dB amplitude measures
if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'})) && ~all(EPmain.view.correl(EPmain.view.dataset <= length(EPdataset.dataset)))
    if (EPmain.view.FFTunits == 4)
        if bottomVolt == -flintmax
            bottomVolt=-inf; %log10 of zero is -inf.  Reverse replacing with maximum possible double-precision negative number.
        end;
        bottomVolt=10^(bottomVolt/10);
        if topVolt == -flintmax
            topVolt=-inf; %log10 of zero is -inf.  Reverse replacing with maximum possible double-precision negative number.
        end;
        topVolt=10^(topVolt/10);
    end;
    if (EPmain.view.FFTunits > 2)
        bottomVolt=sqrt(bottomVolt); %convert power to amplitude
        topVolt=sqrt(topVolt); %convert power to amplitude
    end;
    bottomVolt=bottomVolt*sqrt(EPmain.view.binHz); %convert from spectral density
    topVolt=topVolt*sqrt(EPmain.view.binHz); %convert from spectral density
end;

if isempty(startHz) || isempty(endHz) || isempty(startSamp) || isempty(endSamp) || isempty(bottomVolt) || isempty(topVolt)
    EPmain.view.changeFlag(1:4)=1;
end;

if src==EPmain.handles.view.bottomVolt
    if isempty(bottomVolt)
        EPmain.view.edited.bottomVolt=0;
    else
        EPmain.view.edited.bottomVolt=1;
        EPmain.view.manual.bottomVolt=bottomVolt;
    end;
end;

if src==EPmain.handles.view.topVolt
    if isempty(topVolt)
        EPmain.view.edited.topVolt=0;
    else
        EPmain.view.edited.topVolt=1;
        EPmain.view.manual.topVolt=topVolt;
    end;
end;

if src==EPmain.handles.view.startSamp
    if isempty(startSamp)
        EPmain.view.edited.startSamp=0;
    else
        EPmain.view.edited.startSamp=1;
        EPmain.view.manual.startSamp=startSamp;
    end;
end;

if src==EPmain.handles.view.endSamp
    if isempty(endSamp)
        EPmain.view.edited.endSamp=0;
    else
        EPmain.view.edited.endSamp=1;
        EPmain.view.manual.endSamp=endSamp;
    end;
end;

if src==EPmain.handles.view.startHz
    if isempty(startHz)
        EPmain.view.edited.startHz=0;
    else
        EPmain.view.edited.startHz=1;
        EPmain.view.manual.startHz=startHz;
    end;
end;

if src==EPmain.handles.view.endHz
    if isempty(endHz)
        EPmain.view.edited.endHz=0;
    else
        EPmain.view.edited.endHz=1;
        EPmain.view.manual.endHz=endHz;
    end;
end;

EPmain.view.startSamp=startSamp;
EPmain.view.endSamp=endSamp;
EPmain.view.startHz=startHz;
EPmain.view.endHz=endHz;
EPmain.view.marker1=marker1;
EPmain.view.marker2=marker2;

ep('startView');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeWork(src,eventdata)
%change work directory

global EPdataset EPmain

newDir = uigetdir(EPdataset.EPwork,'Select EPwork directory');
if newDir == 0
    msg{1}='No directory selected.';
    [msg]=ep_errorMsg(msg);
    return
end;

seps=findstr(newDir,filesep);
if ~isempty(seps)
    theDir=newDir(seps(end)+1:end);
    newDir=newDir(1:seps(end)-1);
else
    theDir=newDir;
    newDir='.';
end;

if ~strcmp(theDir,'EPwork')
    msg{1}='Directory not named EPwork.';
    [msg]=ep_errorMsg(msg);
    return
end;

if ~exist([newDir filesep 'EPwork' filesep 'EPdataset.mat'],'file')
    msg{1}='EPwork directory missing EPdataset file.';
    [msg]=ep_errorMsg(msg);
    return
end;

if ~exist([newDir filesep 'EPwork' filesep 'EPprefs.mat'],'file')
    msg{1}='EPwork directory missing EPprefs file.';
    [msg]=ep_errorMsg(msg);
    return
end;

tempVar=load([newDir filesep 'EPwork' filesep 'EPdataset.mat']);
if isfield(tempVar,'EPdataset')
    EPdataset=tempVar.EPdataset;
end;
clear tempVar;

tempVar=load([newDir filesep 'EPwork' filesep 'EPprefs.mat']);
if isfield(tempVar,'prefs')
    prefs=tempVar.prefs;
end;
clear tempVar;

EPmain.preferences=prefs; 
EPdataset.EPwork=newDir;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function createWork(src,eventdata)
%create work directory in current directory

global EPdataset EPmain

if exist([pwd filesep 'EPwork'],'dir')
    msg{1}='EPwork directory already present.';
    [msg]=ep_errorMsg(msg);
    return
end;

mkdir([pwd filesep 'EPwork']);
EPdataset=[];
EPdataset.EPwork=pwd;
EPdataset.dataset=cell(0);
eval(['save ''' pwd filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
prefs=EPmain.preferences;
eval(['save ''' pwd filesep 'EPwork' filesep 'EPprefs'' prefs']);

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentDataset(src,eventdata)
%change the dataset in the segment function

global EPmain EPdataset

if isfield(EPmain.segment,'dataset')
    theDataset=get(EPmain.handles.segment.dataset,'Value');
    theEvent=get(EPmain.handles.segment.eventValues,'Value');
else
    theDataset=length(EPmain.segment.contData);
    theEvent=1;
end;

EPmain.segment.dataset=EPmain.segment.contData(theDataset);
EPmain.segment.event=theEvent;

if strcmp(EPdataset.dataset(EPmain.segment.dataset).dataType,'continuous')
    allEvents=EPdataset.dataset(EPmain.segment.dataset).events;
else
    allEvents=cell(1);
    for iTrial=1:length(EPdataset.dataset(EPmain.segment.dataset).events)
        allEvents{1}=[allEvents{1} EPdataset.dataset(EPmain.segment.dataset).events{1,iTrial}];
    end;
end;

if ~isempty(EPdataset.dataset(EPmain.segment.dataset).events{1})
    eventTypes={allEvents{1}.type}';
else
    eventTypes=[];
end;
if ~isempty(EPdataset.dataset(EPmain.segment.dataset).events{1})
    eventValues={allEvents{1}.value}';
else
    eventValues=[];
end;

%[eventValues{end+1:end+length(tempValues)}]=tempValues;
EPmain.segment.eventTypes=unique(cellfun(@num2str,eventTypes(find(~cellfun(@isempty,eventTypes))),'UniformOutput',false'));
EPmain.segment.eventTypes=sort(EPmain.segment.eventTypes);
EPmain.segment.eventValues=unique(cellfun(@num2str,eventValues(find(~cellfun(@isempty,eventValues))),'UniformOutput',false'));
EPmain.segment.eventValues=sort(EPmain.segment.eventValues);

%remove from EventType if all of a Type also have the same name for the Value (e.g., 'boundary' type and 'boundary' value) and also TRSP events.
dropList=[];
for iEvent=1:length(EPmain.segment.eventTypes)
    theType=EPmain.segment.eventTypes{iEvent};
    theValues=unique(cellfun(@num2str,eventValues(find(~cellfun(@isempty,eventValues)&strcmp(theType,eventTypes))),'UniformOutput',false'));
    if (~isempty(theValues) && (length(unique(theValues))==1) && strcmp(theValues{1},theType)) || strcmp(theType,'TRSP')
        dropList(end+1)=iEvent;
    end;
end;
EPmain.segment.eventTypes(dropList)=[];
EPmain.segment.eventValues(strcmp('TRSP',EPmain.segment.eventValues))=[];

numEventTypes=length(EPmain.segment.eventTypes);
numEventValues=length(EPmain.segment.eventValues);
EPmain.segment.eventLocks=EPmain.segment.eventTypes;
EPmain.segment.eventLocks(end+1:end+numEventValues)=EPmain.segment.eventValues;

EPmain.segment.trialSpecNames=cell(0); %the names of the trial spec fields (cell array)
EPmain.segment.trialSpecValues=cell(0); %the values present for each of these trial spec fields (cell of cells)

EPmain.segment.stimSpecNames=cell(numEventTypes+numEventValues,1); %the list of key field names for each event value (cell of cell arrays)
EPmain.segment.stimSpecValues=cell(numEventTypes+numEventValues,1); %the values present for each of these  fields (cell of cell of cell arrys)

EPmain.segment.critSpecNames=cell(0); %the menu of names listed for each of the six crits
EPmain.segment.critSpecNames{1}='none';
EPmain.segment.critSpecNames{2}='-precedes-';
EPmain.segment.critSpecNames{3}='-follows-';
EPmain.segment.critSpecNames{4}='-secNum-';
EPmain.segment.critSpecNames{5}='-secTrial-';
EPmain.segment.critSpecNames{6}='-secTrialFromEnd-';

EPmain.segment.critSpecValues=cell(0); %the menu of values that would be listed for each of the names in the crit name menus.
EPmain.segment.critSpecValues{1}{1}='none';
EPmain.segment.critSpecValues{2}=EPmain.segment.eventLocks; %possible values for '-precedes-'
EPmain.segment.critSpecValues{3}=EPmain.segment.eventLocks; %possible values for '-follows-'
EPmain.segment.critSpecValues{4}={'0'};
EPmain.segment.critSpecValues{5}={'0'};
EPmain.segment.critSpecValues{6}={'0'};

if strcmp(EPdataset.dataset(EPmain.segment.dataset).dataType,'single_trial')
    for iTRSP=1:length(EPdataset.dataset(EPmain.segment.dataset).trialSpecNames)
        EPmain.segment.trialSpecNames{end+1}=EPdataset.dataset(EPmain.segment.dataset).trialSpecNames{iTRSP};
        EPmain.segment.trialSpecValues{end+1}=cell(0);
        for iTrial=1:length(EPdataset.dataset(EPmain.segment.dataset).trialSpecs)
            if ~any(find(strcmp(num2str(EPdataset.dataset(EPmain.segment.dataset).trialSpecs{iTrial,iTRSP}),EPmain.segment.trialSpecValues{end})))
                EPmain.segment.trialSpecValues{end}{end+1}=num2str(EPdataset.dataset(EPmain.segment.dataset).trialSpecs{iTrial,iTRSP});
            end;
        end;
    end;
else
    TRSPevents=find(strcmp('TRSP',eventValues));
    if ~isempty(TRSPevents)
        for iTRSP=1:length(TRSPevents)
            theTRSP=TRSPevents(iTRSP);
            for iKey=1:length(allEvents{1}(theTRSP).keys)
                if ~any(strcmp(['TS-' allEvents{1}(theTRSP).keys(iKey).code],EPmain.segment.trialSpecNames))
                    EPmain.segment.trialSpecNames{end+1}=['TS-' allEvents{1}(theTRSP).keys(iKey).code];
                    if strcmp(EPmain.segment.trialSpecNames{end},'TS-RT')
                        EPmain.segment.trialSpecValues{end+1}{1}='100'; %100ms always on the list of RT values.
                        if ~strcmp(num2str(allEvents{1}(theTRSP).keys(iKey).data),'100')
                            EPmain.segment.trialSpecValues{end}{2}=num2str(allEvents{1}(theTRSP).keys(iKey).data);
                        end;
                    else
                        EPmain.segment.trialSpecValues{end+1}{1}=num2str(allEvents{1}(theTRSP).keys(iKey).data);
                    end;
                else
                    theTrialSpec=find(strcmp(['TS-' allEvents{1}(theTRSP).keys(iKey).code],EPmain.segment.trialSpecNames));
                    if ~any(find(strcmp(num2str(allEvents{1}(theTRSP).keys(iKey).data),EPmain.segment.trialSpecValues{theTrialSpec})))
                        EPmain.segment.trialSpecValues{theTrialSpec}{end+1}=num2str(allEvents{1}(theTRSP).keys(iKey).data);
                    end;
                end;
            end;
        end;
    end;
end;

for iEventType=1:numEventTypes
    lockEvents=find(strcmp(EPmain.segment.eventTypes{iEventType},eventTypes));
    for iEvent=1:length(lockEvents)
        theEvent=lockEvents(iEvent);
        theValue=num2str(allEvents{1}(theEvent).value);
        if isempty(theValue)
            theValue='none';
        end;
        if iEvent==1
            EPmain.segment.stimSpecNames{iEventType}{end+1}='value';
            EPmain.segment.stimSpecValues{iEventType}{end+1}{1}=theValue;
        else
            if ~any(find(strcmp(theValue,EPmain.segment.stimSpecValues{iEventType}{1})))
                EPmain.segment.stimSpecValues{iEventType}{1}{end+1}=theValue;
            end;
        end;
        for iKey=1:length(allEvents{1}(theEvent).keys)
            if ~any(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{iEventType}))
                EPmain.segment.stimSpecNames{iEventType}{end+1}=allEvents{1}(theEvent).keys(iKey).code;
                EPmain.segment.stimSpecValues{iEventType}{end+1}{1}=allEvents{1}(theEvent).keys(iKey).data;
            else
                theTrialSpec=find(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{iEventType}));
                if ~any(find(strcmp(allEvents{1}(theEvent).keys(iKey).data,EPmain.segment.stimSpecValues{iEventType}{theTrialSpec})))
                    EPmain.segment.stimSpecValues{iEventType}{theTrialSpec}{end+1}=allEvents{1}(theEvent).keys(iKey).data;
                end;
            end;
        end;
    end;
end;

for iEventType=1:numEventValues
    lockEvents=find(strcmp(EPmain.segment.eventValues{iEventType},eventValues));
    for iEvent=1:length(lockEvents)
        theEvent=lockEvents(iEvent);
        for iKey=1:length(allEvents{1}(theEvent).keys)
            if ~any(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{numEventTypes+iEventType}))
                EPmain.segment.stimSpecNames{numEventTypes+iEventType}{end+1}=allEvents{1}(theEvent).keys(iKey).code;
                EPmain.segment.stimSpecValues{numEventTypes+iEventType}{end+1}{1}=allEvents{1}(theEvent).keys(iKey).data;
            else
                theTrialSpec=find(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{numEventTypes+iEventType}));
                if ~any(find(strcmp(allEvents{1}(theEvent).keys(iKey).data,EPmain.segment.stimSpecValues{numEventTypes+iEventType}{theTrialSpec})))
                    EPmain.segment.stimSpecValues{numEventTypes+iEventType}{theTrialSpec}{end+1}=allEvents{1}(theEvent).keys(iKey).data;
                end;
            end;
        end;
    end;
end;

EPmain.segment.trialSpec(1:EPmain.segment.numSpecs)=1;
EPmain.segment.trialSpecRel(1:EPmain.segment.numSpecs)=1;
EPmain.segment.trialSpecVal(1:EPmain.segment.numSpecs)=1;

EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.trialSpecNames))=EPmain.segment.trialSpecNames; %the names of the fields to be listed for a criterion that does not come after a "-follows-" (cell array)
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.trialSpecValues))=EPmain.segment.trialSpecValues; %the contents of the fields to be listed for a criterion that does not come after a "-follows-" (array of cell arrays)
EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.stimSpecNames{EPmain.segment.event}))=EPmain.segment.stimSpecNames{EPmain.segment.event};
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.stimSpecValues{EPmain.segment.event}))=EPmain.segment.stimSpecValues{EPmain.segment.event};

for iCrit=1:EPmain.segment.numSpecs
    EPmain.segment.critSpecItemNames{iCrit}=EPmain.segment.critSpecNames; %the actual list of field names to be listed for each of the five model criterion rows
    EPmain.segment.critSpecItemValues{iCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(iCrit)}; %the possible values for each of the actual field names to be listed for each of the five model criterion rows
end;
EPmain.segment.critSpec{EPmain.segment.numSpecs}=EPmain.segment.critSpecNames([1 4:end]); %drop the "-precedes-" and "-follows-" options

EPmain.segment.sampStart=-floor(200/(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)); %default 200 ms
EPmain.segment.sampEnd=floor(1000/(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)); %default 1000 ms


ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentEvent(src,eventdata)
%change the lock event in the segment function

global EPmain EPdataset

EPmain.segment.event=get(EPmain.handles.segment.eventValues,'Value');

EPmain.segment.critSpecNames=cell(0);
EPmain.segment.critSpecNames{1}='none';
EPmain.segment.critSpecNames{2}='-precedes-';
EPmain.segment.critSpecNames{3}='-follows-';
EPmain.segment.critSpecNames{4}='-secNum-';
EPmain.segment.critSpecNames{5}='-secTrial-';
EPmain.segment.critSpecNames{6}='-secTrialFromEnd-';
EPmain.segment.critSpecValues=cell(0);
EPmain.segment.critSpecValues{1}{1}='none';
EPmain.segment.critSpecValues{2}=EPmain.segment.eventLocks;
EPmain.segment.critSpecValues{3}=EPmain.segment.eventLocks;
EPmain.segment.critSpecValues{4}={'0'};
EPmain.segment.critSpecValues{5}={'0'};
EPmain.segment.critSpecValues{6}={'0'};

EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.trialSpecNames))=EPmain.segment.trialSpecNames;
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.trialSpecValues))=EPmain.segment.trialSpecValues;
EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.stimSpecNames{EPmain.segment.event}))=EPmain.segment.stimSpecNames{EPmain.segment.event};
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.stimSpecValues{EPmain.segment.event}))=EPmain.segment.stimSpecValues{EPmain.segment.event};

EPmain.segment.trialSpec(1:EPmain.segment.numSpecs)=1; %which field
EPmain.segment.trialSpecRel(1:EPmain.segment.numSpecs)=1; %which relation
EPmain.segment.trialSpecVal(1:EPmain.segment.numSpecs)=1; %which value

for iCrit=1:EPmain.segment.numSpecs
    EPmain.segment.critSpecItemNames{iCrit}=EPmain.segment.critSpecNames;
    EPmain.segment.critSpecItemValues{iCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(iCrit)};
end;
EPmain.segment.critSpec{EPmain.segment.numSpecs}=EPmain.segment.critSpecNames([1 4:end]);  %drop the "-precedes-" and "-follows-" options

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentCritName(src,eventdata)
%change a criterion name in the segment function

global EPmain

theCrit=find([EPmain.handles.segment.trialSpecNames] == src);

if EPmain.segment.trialSpec(theCrit) ~= get(EPmain.handles.segment.trialSpecNames(theCrit),'Value') %don't do anything unless the setting was actually changed
    oldName=EPmain.segment.trialSpec(theCrit);
    EPmain.segment.trialSpec(theCrit)=get(EPmain.handles.segment.trialSpecNames(theCrit),'Value');
    EPmain.segment.trialSpecRel(theCrit)=1;
    EPmain.segment.trialSpecVal(theCrit)=1;
    
    %if "-precedes-" or "-follows-" options, then reconfigure the options for the following criterion item
    if any(strcmp(EPmain.segment.critSpecItemNames{theCrit}{EPmain.segment.trialSpec(theCrit)},{'-precedes-','-follows-'}))
        EPmain.segment.trialSpec(theCrit+1)=1;
        EPmain.segment.trialSpecRel(theCrit+1)=1;
        EPmain.segment.trialSpecVal(theCrit+1)=1;
        EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.eventLocks;
        EPmain.segment.critSpecItemNames{theCrit+1}=EPmain.segment.stimSpecNames{EPmain.segment.trialSpecVal(theCrit)};
        EPmain.segment.critSpecItemNames{theCrit+1}{end+1}='none';
        if ~isempty(EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)})
            EPmain.segment.critSpecItemValues{theCrit+1}=EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)}{EPmain.segment.trialSpec(theCrit+1)};
        end;
        EPmain.segment.critSpecItemValues{theCrit+1}{end+1}='none';
    elseif any(strcmp(EPmain.segment.critSpecItemNames{theCrit}{oldName},{'-precedes-','-follows-'}))
        %if changing away from "-follows-" option, then reconfigure the options for the following criterion item
        EPmain.segment.trialSpec(theCrit+1)=1;
        EPmain.segment.trialSpecRel(theCrit+1)=1;
        EPmain.segment.trialSpecVal(theCrit+1)=1;
        EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(theCrit)};
        EPmain.segment.critSpecItemNames{theCrit+1}=EPmain.segment.critSpecNames;
        EPmain.segment.critSpecItemValues{theCrit+1}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(theCrit+1)};
    else
        if (theCrit>1) && any(strcmp(EPmain.segment.critSpecItemNames{theCrit-1}{EPmain.segment.trialSpec(theCrit-1)},{'-precedes-','-follows-'}))
            %if line after a -follows- criterion
            if strcmp('none',EPmain.segment.critSpecItemNames{theCrit}{EPmain.segment.trialSpec(theCrit)})
                EPmain.segment.critSpecItemValues{theCrit}=cell(0);
                EPmain.segment.critSpecItemValues{theCrit}{1}='none';
            else
                EPmain.segment.critSpecItemValues{theCrit}=cell(0);
                EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit-1)}{EPmain.segment.trialSpec(theCrit)};
            end;
        else
            EPmain.segment.critSpecItemValues{theCrit}=cell(0);
            EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(theCrit)};
        end;
    end;
end;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentCritValue(src,eventdata)
%change a criterion value in the segment function

global EPmain EPdataset

theCrit=find([EPmain.handles.segment.trialSpecVal] == src);

if EPmain.segment.trialSpecVal(theCrit) ~= get(EPmain.handles.segment.trialSpecVal(theCrit),'Value') %don't do anything unless the setting was actually changed
    EPmain.segment.trialSpecVal(theCrit)=get(EPmain.handles.segment.trialSpecVal(theCrit),'Value');
    
    %if "-precedes-" or "-follows-" options, then reconfigure the options for the following criterion item
    if any(strcmp(EPmain.segment.critSpecItemNames{theCrit}{EPmain.segment.trialSpec(theCrit)},{'-precedes-','-follows-'}))
        EPmain.segment.trialSpec(theCrit+1)=1;
        EPmain.segment.trialSpecRel(theCrit+1)=1;
        EPmain.segment.trialSpecVal(theCrit+1)=1;
        theName=EPmain.segment.stimSpecNames{EPmain.segment.trialSpecVal(theCrit)};
        if (length(theName)==1) && isempty(theName{1})
            EPmain.segment.critSpecItemNames{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemNames{theCrit+1}{1}='none';
        else
            EPmain.segment.critSpecItemNames{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemNames{theCrit+1}=theName;
            EPmain.segment.critSpecItemNames{theCrit+1}{end+1}='none';
        end;
        if ~isempty(EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)})
            theValue=EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)}{EPmain.segment.trialSpec(theCrit+1)};
        else
            theValue=cell(0);
        end;
        if (length(theValue)==1) && isempty(theValue{1})
            EPmain.segment.critSpecItemValues{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemValues{theCrit+1}{1}='none';
        else
            EPmain.segment.critSpecItemValues{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemValues{theCrit+1}=theValue;
            EPmain.segment.critSpecItemValues{theCrit+1}{end+1}='none';
        end;
    end;
end;

ep('start');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segmentSampStart(src,eventdata) %change start of segment

global EPmain EPdataset

Fs=EPdataset.dataset(EPmain.segment.dataset).Fs;

if src==EPmain.handles.segment.sampStart
    sampStart=str2num(get(EPmain.handles.segment.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.segment.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.segment.msStart,'String'))/(1000/Fs)));
    sampEnd=round((str2num(get(EPmain.handles.segment.msEnd,'String'))/(1000/Fs)));
end;

if sampStart > sampEnd
    sampEnd = sampStart;
end;

set(EPmain.handles.segment.sampStart,'String',sampStart);
set(EPmain.handles.segment.sampEnd,'String',sampEnd);
set(EPmain.handles.segment.msStart,'String',(sampStart)*(1000/Fs));
set(EPmain.handles.segment.msEnd,'String',(sampEnd)*(1000/Fs));

EPmain.segment.sampStart=sampStart;
EPmain.segment.sampEnd=sampEnd;
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentSampEnd(src,eventdata) %change end of segment

global EPmain EPdataset

Fs=EPdataset.dataset(EPmain.segment.dataset).Fs;

if src==EPmain.handles.segment.sampEnd
    sampStart=str2num(get(EPmain.handles.segment.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.segment.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.segment.msStart,'String'))/(1000/Fs)));
    sampEnd=round((str2num(get(EPmain.handles.segment.msEnd,'String'))/(1000/Fs)));
end;

if sampStart > sampEnd
    sampStart = sampEnd;
end;

set(EPmain.handles.segment.sampStart,'String',sampStart);
set(EPmain.handles.segment.sampEnd,'String',sampEnd);
set(EPmain.handles.segment.msStart,'String',(sampStart)*(1000/Fs));
set(EPmain.handles.segment.msEnd,'String',(sampEnd)*(1000/Fs));

EPmain.segment.sampStart=sampStart;
EPmain.segment.sampEnd=sampEnd;
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentAddCell(src,eventdata) %add cell to segmenting table

global EPmain

if all(cellfun(@isempty,EPmain.segment.cellTable(end,:)))
    theCell=size(EPmain.segment.cellTable,1);
else
    theCell=size(EPmain.segment.cellTable,1)+1;
end;

EPmain.segment.cellTable{theCell,1}=[];
EPmain.segment.cellTable{theCell,2}='';
EPmain.segment.cellTable{theCell,3}=EPmain.segment.eventLocks{EPmain.segment.event};
if ~EPmain.segment.flexible
    EPmain.segment.cellTable{theCell,4}=get(EPmain.handles.segment.msStart,'String');
    EPmain.segment.cellTable{theCell,5}=get(EPmain.handles.segment.msEnd,'String');
else
    EPmain.segment.cellTable{theCell,4}=EPmain.segment.eventLocks{EPmain.segment.flexEnd};
    EPmain.segment.cellTable{theCell,5}=['F' get(EPmain.handles.segment.flexLength,'String')];
end;
EPmain.segment.cellTable{theCell,6}=get(EPmain.handles.segment.delay,'String');

for iSpec=1:EPmain.segment.numSpecs
    EPmain.segment.cellTable{theCell,7+(iSpec-1)*3}=EPmain.segment.critSpecItemNames{iSpec}{EPmain.segment.trialSpec(iSpec)};
    EPmain.segment.cellTable{theCell,8+(iSpec-1)*3}=EPmain.segment.relList{EPmain.segment.trialSpecRel(iSpec)};
    EPmain.segment.cellTable{theCell,9+(iSpec-1)*3}=EPmain.segment.critSpecItemValues{iSpec}{EPmain.segment.trialSpecVal(iSpec)};
end;

set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentDelCell(src,eventdata) %delete cell of segmenting table

global EPmain

if size(EPmain.segment.cellTable,1) ==1
    for i=1:size(EPmain.segment.cellTable,2)
        EPmain.segment.cellTable{i}='';
    end;
else
    EPmain.segment.cellTable=EPmain.segment.cellTable(1:end-1,:);
end;

 set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
 ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentLoad(src,eventdata) %load segmenting table

global EPmain

[FileName,PathName,FilterIndex] = uigetfile('','Load Segmenting Table');
if FileName ~= 0
    eval(['tempVar=load( ''' PathName FileName ''');']);
    if isfield(tempVar,'EPsegmentTable')
        if size(tempVar.EPsegmentTable,2)==20 %backward compatibility
            tempVar.EPsegmentTable(:,7:EPmain.segment.numSpecs*3+1)=tempVar.EPsegmentTable(:,6:EPmain.segment.numSpecs*3);
            [tempVar.EPsegmentTable{:,6}]=deal(0);
        end        
        if size(tempVar.EPsegmentTable,2)==21 %backward compatibility
            tempVar.EPsegmentTable(:,22:22+(EPmain.segment.numSpecs-5)*3-1)=deal(0);
        end
        tempVar.EPsegmentTable=cellfun(@num2str,tempVar.EPsegmentTable,'UniformOutput',false);
        EPmain.segment.cellTable=tempVar.EPsegmentTable;
        set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
    else
        msg{1}='Error: not a segmenting table file.';
        [msg]=ep_errorMsg(msg);
    end;
    ep('start');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentSave(src,eventdata) %save segmenting table

global EPmain EPdataset

[fileName,PathName,FilterIndex] = uiputfile('*.mat','Save Segmenting Table',[EPdataset.dataset(EPmain.segment.dataset).dataName '_segmentTable']);

if fileName ~= 0
    EPsegmentTable=EPmain.segment.cellTable;
    eval(['save ''' PathName fileName ''' EPsegmentTable']);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convertFiles(src,eventdata) %batch convert files

global EPmain

set(EPmain.handles.save.convert,'enable','off');
set(EPmain.handles.save.done,'enable','off');

drawnow

typeNum = get(EPmain.handles.save.type,'value');
switch typeNum
    case 1
        dataType='continuous';
    case 2
        dataType='single_trial';
    case 3
        dataType='average';
    case 4
        dataType='grand_average';
    case 5
        dataType='factors';
end;

formatNum = get(EPmain.handles.save.readFormat,'value');
[importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{formatNum});

if strcmp(importFormat,'ep_mat')
    dataType=''; %data type set by data file itself
end;

EPmain.save.readFormat=formatNum;
EPmain.save.type=typeNum;

EPmain.convertMode=1;
theHandles=EPmain.handles.save;
readFiles(theHandles,importFormat,dataType);
set(EPmain.handles.save.convert,'enable','on');
set(EPmain.handles.save.done,'enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSampleTestDataset(src,eventdata)
%change the dataset in the sampleTest function

global EPmain EPdataset

EPmain.sampleTest.dataset=get(EPmain.handles.sampleTest.dataset,'Value');
EPmain.sampleTest.dataType=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).dataType;
EPmain.sampleTest.cell1=1;
EPmain.sampleTest.cell2=2;
EPmain.sampleTest.channel=1;
EPmain.sampleTest.cellNameList=unique(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).cellNames(find(strcmp(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).cellTypes,'SGL'))));
EEGchans=find(strcmp('EEG',EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).chanTypes));
EPmain.sampleTest.eloc=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).eloc(EEGchans);
EPmain.sampleTest.changeFlag=1;
EPmain.sampleTest.freqFlag=~isempty(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames);

maxRad=0.5;
[y,x] = pol2cart(([EPmain.sampleTest.eloc.theta]/360)*2*pi,[EPmain.sampleTest.eloc.radius]);  % transform electrode locations from polar to cartesian coordinates
y=-y; %flip y-coordinate so that nose is upwards.
plotrad = min(1.0,max([EPmain.sampleTest.eloc.radius])*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
x = x*(maxRad/plotrad);
y = y*(maxRad/plotrad);

xmin = min(-maxRad,min(x));
xmax = max(maxRad,max(x));
ymin = min(-maxRad,min(y));
ymax = max(maxRad,max(y));

EPmain.sampleTest.x=round(((x/(xmax-xmin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));
EPmain.sampleTest.y=round(((y/(ymax-ymin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));

EPmain.sampleTest.PCAlist=[];
for iFile=1:length(EPdataset.dataset)
    if ~isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc)
        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
        if length(EPmain.sampleTest.eloc) == length(newEloc)
            if ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                if (~strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) == length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames))) ||...
                        (strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames)))
                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames])
                        EPmain.sampleTest.PCAlist(end+1)=iFile;
                    end;
                end;
            end;
        end;
    end;
end;
EPmain.sampleTest.PCAnameList={EPdataset.dataset(EPmain.sampleTest.PCAlist).dataName};

EPmain.sampleTest.AVElist=[];
for iFile=1:length(EPdataset.dataset)
    if isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc) && strcmp(EPdataset.dataset(iFile).dataType,'average')
        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
        if length(EPmain.sampleTest.eloc) == length(newEloc)
            if ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                if (~strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) == length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames))) ||...
                        (strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).timeNames)))
                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).freqNames])
                        EPmain.sampleTest.AVElist(end+1)=iFile;
                        EPdata=ep_loadEPdataset(EPdataset,iFile);
                        EPmain.sampleTest.AVEdata{end+1}=EPdata.data;
                    end;
                end;
            end;
        end;
    end;
end;
EPmain.sampleTest.AVEnameList={EPdataset.dataset(EPmain.sampleTest.AVElist).dataName};

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function expandTransformChan(src,eventdata)
global EPmain

scrsz = EPmain.scrsz;

if isfield(EPmain.handles.transform,'expandedFigure')
    if ~isempty(EPmain.handles.transform.expandedFigure)
        if ishandle(EPmain.handles.transform.expandedFigure)
            close(EPmain.handles.transform.expandedFigure)
            EPmain.handles.transform.expandedFigure=[];
        end;
    end;
end;

if (src == EPmain.handles.transform.timeFilter) || (src == EPmain.handles.transform.timeWaves(1)) || (src == EPmain.handles.transform.timeWaves(2))
    name='Time Domain';
    theData=EPmain.transform.timeData;
    theScale=EPmain.transform.timeScale;
    theAxis=EPmain.transform.timeAxis;
    colorOrder=get(EPmain.handles.transform.timeFilter,'ColorOrder');
else
    name='Frequency Domain';
    theData=EPmain.transform.freqData;
    theScale=EPmain.transform.freqScale;
    theAxis=EPmain.transform.freqAxis;
    colorOrder=get(EPmain.handles.transform.frequencyFilter,'ColorOrder');
end

EPmain.handles.transform.expandedFigure = figure('Name', name, 'NumberTitle', 'off', 'Position',[scrsz(3)/2 scrsz(4)/2 600 400]);

set(EPmain.handles.transform.expandedFigure,'DefaultAxesColorOrder',colorOrder);
EPmain.handles.transform.expandedAxes = axes('position',[.05 .05 .9 .9]);

plot(theScale,theData);

axis(theAxis);

line([theScale(1) theScale(end)],[0 0],'Color','black','LineWidth',1) % zero line
if src == EPmain.handles.transform.frequencyFilter
    if ~isempty(EPmain.transform.filter1)
        axis(theAxis);
        line([EPmain.transform.filter1 EPmain.transform.filter1],[theAxis(3) theAxis(4)],'Color','black','LineWidth',1) % filter line
    end
    if ~isempty(EPmain.transform.filter2)
        axis(theAxis);
        line([EPmain.transform.filter2 EPmain.transform.filter2],[theAxis(3) theAxis(4)],'Color','black','LineWidth',1) % filter line
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readFiles(theHandles,importFormat,dataType)
%choose a set of files to read in
global EPmain EPdataset

textPrefs.firstRow=EPmain.preferences.general.firstRow;
textPrefs.lastRow=EPmain.preferences.general.lastRow;
textPrefs.firstCol=EPmain.preferences.general.firstCol;
textPrefs.lastCol=EPmain.preferences.general.lastCol;
textPrefs.orientation=EPmain.preferences.general.orientation;

singleCellMode=get(theHandles.check,'value');

if ~singleCellMode
    
    [fileNames, activeDirectory]=ep_getFilesUI(importFormat);
    if isempty(fileNames) || ((length(fileNames{1})==1) && (fileNames{1}==0))
        msg{1}='No filenames selected. You have to click on a name.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    ced=[];
    montage=[];
    for theFile=1:size(fileNames,2)
        inArg=[];
        inArg{1}='file';
        inArg{2}=[activeDirectory fileNames{theFile}];
        inArg{3}='format';
        inArg{4}=importFormat;
        inArg{5}='screenSize';
        inArg{6}=EPmain.scrsz;
        inArg{7}='FontSize';
        inArg{8}=EPmain.fontsize;
        if ~strcmp(importFormat,'ep_mat')
            inArg{9}='type';
            inArg{10}=dataType;
        end;
        inArg{end+1}='textPrefs';
        inArg{end+1}=textPrefs;
        inArg{end+1}='elecPrefs';
        inArg{end+1}=EPmain.preferences.general.rotateHead;
        if ~isempty(ced)
            inArg{end+1}='ced';
            inArg{end+1}=ced;
            disp('Assuming the electrode information is the same as for the last file.');
        end;
        if ~isempty(montage)
            inArg{end+1}='montage';
            inArg{end+1}=montage;
            disp('Assuming the montage information is the same as for the last file.');
        end;
        SMIsuffix=EPmain.preferences.general.SMIsuffix;
        if ~isempty(SMIsuffix)
            inArg{end+1}='SMIsuffix';
            inArg{end+1}=SMIsuffix;
        end;
        specSuffix=EPmain.preferences.general.specSuffix;
        if ~isempty(specSuffix)
            inArg{end+1}='specSuffix';
            inArg{end+1}=specSuffix;
        end;
        subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
        if ~isempty(subjectSpecSuffix)
            inArg{end+1}='subjectSpecSuffix';
            inArg{end+1}=subjectSpecSuffix;
        end;
        BVheader=EPmain.preferences.general.BVheader;
        if ~isempty(BVheader)
            inArg{end+1}='BVheader';
            inArg{end+1}=BVheader;
        end;
        
        disp(['Reading: ' fileNames{theFile}])
        EPdata=ep_readData(inArg);
        if isempty(EPdata)
            beep();
            disp(['Error: Was unable to read ' inArg{2}]);
        else            
            if ~isempty(EPdata.ced)
                ced=EPdata.ced;
            end;
            if ~isempty(EPdata.montage)
                montage=EPdata.montage;
            end;
            
            if EPmain.convertMode
                saveFile(EPdata,inArg{2});
            else
                EPdataset=ep_saveEPdataset(EPdataset,EPdata);
            end;
        end;
    end;
    
else %single cell mode
    
    subPosStr=get(theHandles.subject,'string');
    cellPosStr=get(theHandles.cell,'string');
    freqPosStr=get(theHandles.freq,'string');
    subPos=str2num(subPosStr);
    cellPos=str2num(cellPosStr);
    freqPos=str2num(freqPosStr);
    
    if min(subPos) < 1
        beep
        set(theHandles.subjectLabel,'ForegroundColor','red');
        drawnow
        return;
    end;
    if min(cellPos) < 1
        beep
        set(theHandles.cellLabel,'ForegroundColor','red');
        drawnow
        return;
    end;
    if min(freqPos) < 1
        beep
        set(theHandles.freqLabel,'ForegroundColor','red');
        drawnow
        return;
    end;
    
    if EPmain.convertMode
        EPmain.save.check=1;
        EPmain.save.subject=subPosStr;
        EPmain.save.cell=cellPosStr;
        EPmain.save.freq=freqPosStr;
    else
        EPmain.read.check=1;
        EPmain.read.subject=subPosStr;
        EPmain.read.cell=cellPosStr;
        EPmain.read.freq=freqPosStr;
    end;
    
    [sessionFiles, pathname]=ep_getFilesUI(importFormat);
    if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
        msg{1}='No filenames selected. You have to click on a name.';
        [msg]=ep_errorMsg(msg);
        return
    end
    for theFile=1:size(sessionFiles,2)
        [pathstr, name, fileSuffix] = fileparts(sessionFiles{theFile});
        if ~isempty(subPos)
            if max(subPos) > length(name)
                msg{1}=['The file name ' sessionFiles{theFile} ' does not have ' num2str(max(subPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end;
            theSubs{theFile}=sessionFiles{theFile}(subPos);
        else
            if any(strcmp(dataType,{'grand_average','factors'}))
                theSubs{theFile}='GAV';
            else
                theSubs{theFile}='S01';
            end;
        end;
        if ~isempty(cellPos)
            if max(cellPos) > length(name)
                msg{1}=['The file name ' sessionFiles{theFile} ' does not have ' num2str(max(cellPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end;
            theCells{theFile}=sessionFiles{theFile}(cellPos);
        else
            theCells{theFile}='C01';
        end;
        if ~isempty(freqPos)
            if max(freqPos) > length(name)
                msg{1}=['The file name ' sessionFiles{theFile} ' does not have ' num2str(max(freqPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end;
            theFreqs{theFile}=str2num(sessionFiles{theFile}(freqPos));
        else
            theFreqs{theFile}=[];
        end;
        sessionFiles{theFile}=[pathname sessionFiles{theFile}];
    end;
    
    mergeName = char(inputdlg('Name of new merged dataset?','Dataset name'));
    pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
    
    if isempty(mergeName)
        msg{1}='No filename given for the new merged file.  Read operation aborted.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    eloc=[];
    ced=[];
    montage=[];
    if strcmp(dataType,'average')
        uniqueSubs=[];
        numSubFiles=1;
    else
        uniqueSubs=unique(theSubs);
        numSubFiles=length(uniqueSubs);
    end;
    
    for sub=1:numSubFiles
        inArg=[];
        if strcmp(dataType,'average')
            theFiles=[1:length(theSubs)];
        else
            theFiles=find(strcmp(uniqueSubs(sub),theSubs));
        end;
        for theFile=1:length(theFiles)
            inArg{theFile,1}=sessionFiles{theFiles(theFile)};
            inArg{theFile,2}='format';
            inArg{theFile,3}=importFormat;
            inArg{theFile,4}='type';
            inArg{theFile,5}=dataType;
            inArg{theFile,6}='labels';
            inArg{theFile,7}={theCells(theFiles(theFile)) theSubs(theFiles(theFile)) theFreqs(theFiles(theFile))};
            inArg{theFile,8}='textPrefs';
            inArg{theFile,9}=textPrefs;
            inArg{theFile,10}='elecPrefs';
            inArg{theFile,11}=EPmain.preferences.general.rotateHead;
            inArg{theFile,12}='screenSize';
            inArg{theFile,13}=EPmain.scrsz;
            inArg{theFile,14}='FontSize';
            inArg{theFile,15}=EPmain.fontsize;
            if ~isempty(eloc)
                inArg{theFile,end+1}='eloc';
                inArg{theFile,end+1}=eloc;
            end;
            if ~isempty(ced)
                inArg{theFile,end+1}='ced';
                inArg{theFile,end+1}=ced;
                disp('Assuming the electrode information is the same as for the last file.');
            end;
            if ~isempty(montage)
                inArg{theFile,end+1}='montage';
                inArg{theFile,end+1}=montage;
                disp('Assuming the montage information is the same as for the last file.');
            end;
            SMIsuffix=EPmain.preferences.general.SMIsuffix;
            if ~isempty(SMIsuffix)
                inArg{theFile,end+1}='SMIsuffix';
                inArg{theFile,end+1}=SMIsuffix;
            end;
            specSuffix=EPmain.preferences.general.specSuffix;
            if ~isempty(specSuffix)
                inArg{theFile,end+1}='specSuffix';
                inArg{theFile,end+1}=specSuffix;
            end;
            subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
            if ~isempty(subjectSpecSuffix)
                inArg{theFile,end+1}='subjectSpecSuffix';
                inArg{theFile,end+1}=subjectSpecSuffix;
            end;
            BVheader=EPmain.preferences.general.BVheader;
            if ~isempty(BVheader)
                inArg{theFile,end+1}='BVheader';
                inArg{theFile,end+1}=BVheader;
            end;
        end;
        
        if length(uniqueSubs) > 1 && ~strcmp(dataType,'average')
            fullMergeName=[mergeName num2str(sub)];
        else
            fullMergeName=mergeName;
        end;
        [EPdata eloc]=ep_mergeEPfiles(inArg,fullMergeName);
        if isempty(EPdata)
            beep();
            disp(['Error: Was unable to read ' fullMergeName]);
        else
            if ~isempty(EPdata.ced)
                ced=EPdata.ced;
            end;
            if ~isempty(EPdata.montage)
                montage=EPdata.montage;
            end;
            if EPmain.convertMode
                saveFile(EPdata,[pathname fullMergeName]);
            else
                EPdataset=ep_saveEPdataset(EPdataset,EPdata);
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function crossVerifyData(src, eventdata)
%choose PCA dataset for cross-verification
global EPmain

if isempty(eventdata.Indices) %if just deselecting
    return;
end;

EPmain.pca.crossVerifyPCA=EPmain.pca.PCAdatasets(eventdata.Indices(1));
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sampleTestSampStart(src,eventdata) %change start of sample test template

global EPmain

if src==EPmain.handles.sampleTest.sampStart
    EPmain.sampleTest.sampStart=str2num(get(EPmain.handles.sampleTest.sampStart,'String'))+1+EPmain.sampleTest.baseline;
else
    EPmain.sampleTest.sampStart=round(str2num(get(EPmain.handles.sampleTest.msStart,'String'))/(1000/EPmain.sampleTest.Fs))+EPmain.sampleTest.baseline+1;
end;

if EPmain.sampleTest.sampStart < -EPmain.sampleTest.baseline+1
    EPmain.sampleTest.sampStart = -EPmain.sampleTest.baseline+length(EPmain.sampleTest.templateWaveform);
end;

if EPmain.sampleTest.sampStart > EPmain.sampleTest.sampEnd
    EPmain.sampleTest.sampEnd = sampStart;
end;

set(EPmain.handles.sampleTest.sampStart,'String',EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline+1);
set(EPmain.handles.sampleTest.sampEnd,'String',EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline);
set(EPmain.handles.sampleTest.msStart,'String',(EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline+1)*(1000/EPmain.sampleTest.Fs));
set(EPmain.handles.sampleTest.msEnd,'String',(EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline)*(1000/EPmain.sampleTest.Fs));

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sampleTestSampEnd(src,eventdata) 
%change end of sample test template

global EPmain

if src==EPmain.handles.sampleTest.sampEnd
    EPmain.sampleTest.sampEnd=str2num(get(EPmain.handles.sampleTest.sampEnd,'String'))+EPmain.sampleTest.baseline;
else
    EPmain.sampleTest.sampEnd=round((str2num(get(EPmain.handles.sampleTest.msEnd,'String')))/(1000/EPmain.sampleTest.Fs))+EPmain.sampleTest.baseline;
end;

if EPmain.sampleTest.sampEnd > -EPmain.sampleTest.baseline+length(EPmain.sampleTest.templateWaveform)
    EPmain.sampleTest.sampEnd = -EPmain.sampleTest.baseline+length(EPmain.sampleTest.templateWaveform);
end;

if EPmain.sampleTest.sampStart > EPmain.sampleTest.sampEnd
    EPmain.sampleTest.sampStart = EPmain.sampleTest.sampEnd;
end;

set(EPmain.handles.sampleTest.sampStart,'String',EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1);
set(EPmain.handles.sampleTest.sampEnd,'String',EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline);
set(EPmain.handles.sampleTest.msStart,'String',(EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1)*(1000/EPmain.sampleTest.Fs));
set(EPmain.handles.sampleTest.msEnd,'String',(EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline)*(1000/EPmain.sampleTest.Fs));

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clearWorkingSet(src,eventdata)
%clear working set?
global EPdataset

button = questdlg('Clear Working Set?');

switch button
    case 'Yes'
        fileList=dir([EPdataset.EPwork filesep 'EPwork']);
        for iFile=1:length(fileList)
            if ~any(strcmp(fileList(iFile).name,{'.','..','EPprefs.mat'}))
                delete([EPdataset.EPwork filesep 'EPwork' filesep fileList(iFile).name]);
            end;
        end;
        EPdataset.dataset=cell(0);
    case 'No'
        
    case 'Cancel'
end;

