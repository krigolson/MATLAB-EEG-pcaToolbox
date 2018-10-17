function ep_editData(varargin)
% ep_editData - ep_editData(varargin) -
% Allows one to examine and edit the structure of the data.
%
%Input:
%  EPdata         : Structured array with the data and accompanying information.  See readData.
%
%   EPoverview     : Structured array with the information for the current overview.
%     .lastData: data prior to last change.  To permit undo.  Structure same as for data.
%     .workData: current working copy of data.  Only made permanent if "save" is chosen.
%     .handles: handles for GUI objects
%        cells.hCellName    :Names of cells (cells)
%        cells.hDeleteCell  :delete button for each cell (cells)
%        cells.hPlusCell    :plus button for each cell (cells)
%        cells.hSubtractCell  :subtract button for each cell (cells)
%        cells.hCombinedCells :The text window for listing the cells to be combined.
%     .combineWeight   :The +/- for each unique cell.
%     .newCellName     :The name for the new cell.
%     .subjects.select  :The subjects that are currently selected.
%
%Output:
%  EPdataset      : Structured array with the data and accompanying information.
%                   See Input.
%  EPoverview     : Structured array with the information for the current analysis.
%                   See Input.

%History
%  by Joseph Dien (3/23/09)
%  jdien07@mac.com
%
%  modified 7/20/09 JD
%  Added ability to reorder cells and to rename cells.
%  Added pages for file overview, subjects, trials, channels, samples, event, and factors.  Most of the data can be edited.
%  Changed global data name from data to EPdataset
%
%  bugfix 9/5/09 JD
%  Location of windows appearing partly off screen on some systems due to offset from having Dock (on Mac) or Task Bar
%  (on PC) at the bottom of the screen.
%
%  bugfix 10/29/09 JD
%  Crash when editing single_trial data with the cells subpane.
%
%  modified 11/6/09 JD
%  Shortened the edit window so that it'll it into smaller monitors.
%
%  bugfix 11/9/09 JD
%  Crash when exporting data from QC and PCA subpanes when using a Matlab version predating 2008.
%
%  bugfix 1/15/10 JD
%  Error when changing ced field in overview subpane of edit pane.
%
%  bugfix 1/29/10 JD
%  Fixed crash when examining QC data on bad channels and trials.
%  Refreshes Edit pane when Done button is pressed so when Data Name is changed, the name in the list of datasets is updated.
%  Fixed crash when adding a cell in the edit pane whose type is 'SGL'.
%
% modified & bugfix 2/27/10 JD
% Changed bad channel field to negative numbers for still bad channels.  Divided QC subpane badChan into badChan and repChan
% subpanes.  Bad channel and rep channel QC figures no longer include bad trials.
% BadChans numbers of QC subpane will now use number of subjects to calculate proportion for grand averages.
% analysis fields no longer optional.
% Fixed crash when adding subject spec.
% Eliminated chantype field for implicit channels.
% Export file dialogs now correctly indicate that the file with be of type .txt.
% Fixed crash when examining factors that include a combined factor add.
%
% modified 5/16/10 JD
% Added summary page to the PCA subpane of the Edit function, including the total variance accounted for information.
%
% bugfix 5/17/10 JD
% Fixed location of the variance accounted for table for old versions of Matlab for two-step PCAs in the Summary
% subpane of the PCA pane.
% Fixed crash when reordering cells.
%
% bugfix 5/23/10 JD
% Fixed assumption that appended files will be in EGIS format.
%
% bugfix 5/27/10 JD
% Fixed baseline control on Samples subpane failing to change baseline value.
%
% modified 6/15/10 JD
% Marks when a file isn't saved yet.
%
% bugfix 7/2/10 JD
% In order fix crash due to Matlab bug on Windows, put in brief pause after spec name requestor diaglog box.
%
% bugfix 8/27/10 JD
% Fixed append cells and append subjects functions not working due to fallacious error messages.
% Fixed append cells and append subjects and append chans crashing when averaged data (typically ept files) contain noise or std information from having run
% the averaging with the Toolkit and the new file is not a factor file.
% Fixed append cells and append subjects and append chans not adding std and noise information from the new file.
% Fixed append cells failing if the files contained trial specs, as in EGIS files.
% Fixed not allowing subject rows to be changeable even when exceeding limits of subject subpane when using a Matlab
% version predating 2008.
%
% bugfix 10/5/10 JD
% Fixed not able to append subjects or cells to data with same names by generating new name.
% Fixed crash when appending cells to non-single trial data.
%
%  modified 10/12/10 JD
%  For continuous files, data now divided into one second epochs and can be artifact rejected in an epochwise fashion
%  in same fashion as segmented data.  Analysis fields reset to zero if sampling rate is changed.
% 
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
% 
%  modified 1/17/11 JD
%  When adding combinations of cells or subjects, the name of the new addition describes what went into it.
% 
% modified 1/24/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
% 
% bugfix 1/30/12 JD
% Changing Fs also changes time names.
% 
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% modified 4/11/12 JD
% Added support for 6th dimension of frequency. Unlogs frequency data for mean computation.
%
% modified 5/24/12 JD
% Added support for missing data when adding channels, cells, or subjects.
% 
% bugfix 5/30/12 JD
% Fixed crash when appending files to data each with just one entry (e.g., appending cells, resulting in two cells) 
% and then trying to append additional files.
% 
% bugfix 6/1/12 JD
% Fixed QC numbers not being computed correctly.
% 
% bugfix 6/6/12 JD
% Fixed changing prestimulus period modified timepoint names relative to prior values rather than absolutely.
%
% modified 7/17/12 JD
% Added option to weight cell combinations by number of trials in averages.
% 
% bugfix 8/4/12 JD
% Fixed crash when examining std table of the QC subpane for frequency data.
%
% modified 9/15/12 JD
% Changed Factors tab so not stripping out adds prior to displaying summary table to allow for datasets consisting only of adds.
% 
% bugfix 10/16/12 JD
% Fixed crash when changing the sampling rate in the samples subpane.
%
% bugfix 10/18/12 JD
% Fixed timeNames field not being a column vector after changing baseline.
%
% bugfix 10/18/12 JD
% Fixed crash when trimming low end Hz of spectral data.
%
% modified 1/29/13 JD
% Added controls to narrow range in time or frequency domain via ms and freq respectively.
%
% bugfix 2/5/13 JD
% Fixed time names not being calculated correctly when sampling rate changed using Edit function.
%
% bugfix 2/6/13 JD
% Fixed %age of blink, saccade, move, and bad trials in QC subpane being computed incorrectly for average data.
%
% modified 5/9/13 JD
% Added option to the Trials subpane to load a text file to rename the cell names of all the trials.
%
% bugfix 5/23/13 JD
% Fixed Export button on Factors tab not working.
%
% modified 9/16/13 JD
% Modified display to work with 640 pixels high screen.
%
% bugfix 9/17/13 JD
% Fixed overview pane graying out number of factors.
%
% modified 10/9/13 JD
% Added recTime field.  Eliminated offset field from events structure.
%
% bugfix 11/1/13 JD
% Fixes font sizes on Windows.
%
% bugfix 11/7/13 JD
% When channel type is changed to REG, ANS, or ECG, electrode coordinates are set to missing.
%
% bugfix 11/21/13 JD
% Fixed latency of events in Edit function off by 4ms.
%
% bugfix 11/28/13 JD
% Fixed added new subject or cell duplicating existing name resulted in corrupted file.
%
% bugfix 12/24/13 JD
% Fixed crash when adding factors together.
%
% bugfix 1/12/14 JD
% Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% modified 1/15/14 JD
% Added option to load in new electrode coordinates from Overview subpane.
%
% bugfix 2/27/14 JD
% Fixed subjects subpane of Edit function specifying average subject type as being AVE rather than AVG.
%
% bufix 3/11/14 JD
% Handles decimal sampling rates gracefully.
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% modified 3/19/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
% Eliminated noTable option for old versions of Matlab.
%
% bufix 3/19/14 JD
% Fixed crash when editing cells.
% Combining of factors now done as simple addition rather than as averaging and no error if add a factor combination to a PCA file with no combined factors.
% Fixed crash when appending channels.
%
% modified 3/24/14 JD
% Added .cov field.
%
% modified 3/28/14 JD
% When the dataset name has been changed on the Overview page, Edit will ask if the unedited dataset should be kept in addition to the edited version when Done is pressed.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% modified 4/16/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% bufix 5/4/14 JD
% Fixed peak latency of factors expressed as one sample too late.
%
% modified 6/1/14 JD
% Added STS cell type.
%
% bufix 6/22/14 JD
% Fixed Edit function reordering cells in single-trial data when anything clicked or changed in the Cells table and they were not already in alphabetical order.
% Fixed deleting wrong cell of single-trial data when cells were not in alphabetical order.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bufix 7/30/14 JD
% Fixed crash when reordering subjects.
%
% modified 9/18/14 JD
% Improved exported events text file.
%
% bufix 9/29/14 JD
% Fixed appending subjects not working.
%
% modified 10/16/14 JD
% Passes screen size to ep_readData call.
%
% bufix 10/21/14 JD
% Fixed appending cells really slow.
%
% bufix 10/24/14 JD
% Fixed crash when adding ced electrode coordinates info.
%
% modified 12/4/14 JD
% Added Import Events button to the Events subpane.
%
% bufix 12/9/14 JD
% Fixed crash when changing channel type and not staying changed when
% changing cell type.
%
% bufix 8/17/15 JD
% Fixed choosing Cancel after pressing Overview>Channel Coordinates
% automatically erased channel coordinates.
%
% modified 8/25/15 JD
% Added feature that when electrode coordinates are replaced with new ones,
% the EEG channels will be remapped via interpolation.
%
% bufix & modified 9/4/15 JD
% Fixed channel weight sum not updating when weights changed.
% Tables rearrangeable.
%
% modified 9/4/15 JD
% Added trial specs for average files and averages their contents.
%
% bufix 10/9/15 JD
% Fixed crash in cells tab for single subject average files with no trial
% specs.
%
% modified 11/3/15 JD
% Added Factors subpane now provides both negative and positive peak chans.
%
% bufix 11/3/15 JD
% Fixed crash when changing order of cells.
%
% bufix 1/22/16 JD
% Now doubles amplitude of spectral data when halving the bin size to correct for what it would have been
% and to ensure that subsequent spectral density scaling will be correct.
%
% bufix 2/19/16 JD
% Fixed computation of peak pos and neg channels and peak channel polarity
% in Factors subpane.
%
% bufix 5/5/16 JD
% New button inactive.
%
% modified 8/18/16 JD
% Append subjects can now select multiple files for appending.
%
% modified 11/5/16 JD
% Added support for reading and writing subject spec text files.
%
% bufix 5/19/17 JD
% Fixed cannot delete FID channels using the Edit function.
%
% modified 11/5/16 JD
% Added support for flexible segments.
%
% modified 11/22/17 JD
% Eliminated restrictions on location of CED files.
% Added support for impedances field.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% modified 3/13/18 JD
% Added option to double sampling rate in the Samples subpane.
%
% bufix 3/26/18 JD
% Fixed adding erroneous sph field to the eloc structure when making certain types of edits to the channels using the Edit function.
%
% bufix 4/26/18 JD
% Fixed crash when reordering cells of PCA data.
%
% bufix & modified 6/13/18 JD
% Added sort popupmenu for subjects, cells, and trials subpanes.
% Fixed crash when editing contents of trials subpane.
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

global EPdataset EPoverview EPmain

if isempty(EPoverview.mode)
    EPoverview.lastData=[];
    EPoverview.workData=ep_loadEPdataset(EPdataset,EPoverview.dataset);
    EPoverview.handles = [];
    EPoverview.handles.hEditWindow =[];
    EPoverview.mode='overview';
    EPoverview.QC.mode='avgNum';
    EPoverview.PCA.mode='Summary';
    EPoverview.cells.TrlWgts=0;
end;

numPoints=length(EPoverview.workData.timeNames);
numChans=length(EPoverview.workData.chanNames);
numCells=length(unique(EPoverview.workData.cellNames));
numWaves=length(EPoverview.workData.cellNames);
numSubs=length(EPoverview.workData.subNames);
numSpecs=length(EPoverview.workData.subjectSpecNames);
numImplicit=length(EPoverview.workData.implicit);
numFacs=length(EPoverview.workData.facNames);
numFreqs=length(EPoverview.workData.freqNames);
numRels=length(EPoverview.workData.relNames);
if ~isempty(EPoverview.workData.facData)
    numCMBfacs=size(EPoverview.workData.facData,5);
else
    numCMBfacs=0;
end;
numSGLfacs=numFacs-numCMBfacs;
refTypes={'REG','AVG','CSD'};
refLabels={'regular reference','average reference','current source density'};

numEvents=0;
for subject=1:numSubs
    for wave=1:numWaves
        numEvents=numEvents+length(EPoverview.workData.events{subject,wave});
    end;
end;

if nargin < 1
    EPoverview.subjects.select=repmat(false,numSubs,1);
    EPoverview.subjects.weights=repmat(0,numSubs,1);
    EPoverview.channels.select=repmat(false,numChans+numImplicit,1);
    EPoverview.channels.weights=repmat(0,numChans+numImplicit,1);
    EPoverview.cells.select=repmat(false,numCells,1);
    EPoverview.cells.weights=repmat(0,numCells,1);
    EPoverview.events.select=repmat(false,numEvents,1);
    EPoverview.events.weights=repmat(0,numEvents,1);
    EPoverview.factors.select=repmat(false,numFacs,1);
    EPoverview.factors.weights=repmat(0,numFacs,1);
end;

if isempty(varargin)
    varargin{1}='start';
end;

scrsz = EPmain.scrsz;

windowHeight=min(scrsz(4)-40,700);


switch varargin{1}
    case 'start'
        if ~isempty(EPoverview.handles.hEditWindow)
            clf(EPoverview.handles.hEditWindow)
        else
            [err]=ep_checkEPfile(EPoverview.workData);
            if err
                disp('File corrupted or out-of-date');
                return
            end;
            if isempty(findobj('name','editData'))
                h=findobj('name','EditData');
                for i=1:length(h)
                    close(h(i))
                end;
            end;
            EPoverview.handles.hEditWindow = figure('Name', 'EditData', 'NumberTitle', 'off', 'Position',[201 scrsz(4)-windowHeight 700 windowHeight], 'MenuBar', 'none');
            colormap jet;
        end;
        
        
        EPoverview.handles.hUndo = uicontrol('Style', 'pushbutton', 'String', 'Undo','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-100 60 30], 'Callback', 'ep_editData(''undo'');');
        
        EPoverview.handles.hDone = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-130 60 30], 'Callback', 'ep_editData(''done'');');
        
        EPoverview.handles.hCancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-160 60 30], 'Callback', 'ep_editData(''cancel'');');
        
        EPoverview.handles.hNew = uicontrol('Style', 'pushbutton', 'String', 'New','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-190 60 30], 'Callback', 'ep_editData(''new'');');
        
        
        
        EPoverview.handles.hOverview = uicontrol('Style', 'pushbutton', 'String', 'Overview','FontSize',EPmain.fontsize,...
            'Position', [100 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''overview'';','ep_editData(''start'');']);
        
        EPoverview.handles.hSubjects = uicontrol('Style', 'pushbutton', 'String', 'Subjects','FontSize',EPmain.fontsize,...
            'Position', [160 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''subjects'';','ep_editData(''start'');']);
        
        EPoverview.handles.hCells = uicontrol('Style', 'pushbutton', 'String', 'Cells','FontSize',EPmain.fontsize,...
            'Position', [220 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''cells'';','ep_editData(''start'');']);
        
        EPoverview.handles.hTrials = uicontrol('Style', 'pushbutton', 'String', 'Trials','FontSize',EPmain.fontsize,...
            'Position', [280 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''trials'';','ep_editData(''start'');']);
        
        if ~strcmp(EPoverview.workData.dataType,'single_trial')
            set(EPoverview.handles.hTrials,'enable','off');
        end;
        
        EPoverview.handles.hChannels = uicontrol('Style', 'pushbutton', 'String', 'Channels','FontSize',EPmain.fontsize,...
            'Position', [340 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''channels'';','ep_editData(''start'');']);
        
        EPoverview.handles.hSamples = uicontrol('Style', 'pushbutton', 'String', 'Samples','FontSize',EPmain.fontsize,...
            'Position', [400 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''samples'';','ep_editData(''start'');']);
        
        if isempty(EPoverview.workData.timeNames)
            set(EPoverview.handles.hSamples,'enable','off');
        end;
        
        EPoverview.handles.hfreqs = uicontrol('Style', 'pushbutton', 'String', 'Hz','FontSize',EPmain.fontsize,...
            'Position', [460 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''freqs'';','ep_editData(''start'');']);
                
        if isempty(EPoverview.workData.freqNames)
            set(EPoverview.handles.hfreqs,'enable','off');
        end;

        EPoverview.handles.hEvents = uicontrol('Style', 'pushbutton', 'String', 'Events','FontSize',EPmain.fontsize,...
            'Position', [520 windowHeight-40 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''events'';','ep_editData(''start'');']);
        
        EPoverview.handles.hFactors = uicontrol('Style', 'pushbutton', 'String', 'Factors','FontSize',EPmain.fontsize,...
            'Position', [100 windowHeight-80 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''factors'';','ep_editData(''start'');']);
        
        if isempty(EPoverview.workData.facNames)
            set(EPoverview.handles.hFactors,'enable','off');
        end;
        
        EPoverview.handles.hQC = uicontrol('Style', 'pushbutton', 'String', 'QC','FontSize',EPmain.fontsize,...
            'Position', [160 windowHeight-80 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''QC'';','ep_editData(''start'');']);
        
        EPoverview.handles.hPCA = uicontrol('Style', 'pushbutton', 'String', 'PCA','FontSize',EPmain.fontsize,...
            'Position', [220 windowHeight-80 60 30], 'Callback', ['global EPoverview;','EPoverview.mode=''PCA'';','ep_editData(''start'');']);
        
        if ~isfield(EPoverview.workData.pca,'PCAmode')
            set(EPoverview.handles.hPCA,'enable','off');
        end;
        
        
        switch EPoverview.mode
            case 'overview'
                set(EPoverview.handles.hOverview,'ForegroundColor','blue');
                ep_editData('startOverview');
                
            case 'subjects'
                set(EPoverview.handles.hSubjects,'ForegroundColor','blue');
                ep_editData('startSubjects');
                
            case 'cells'
                set(EPoverview.handles.hCells,'ForegroundColor','blue');
                ep_editData('startCells');
                
            case 'trials'
                set(EPoverview.handles.hTrials,'ForegroundColor','blue');
                ep_editData('startTrials');
                
            case 'channels'
                set(EPoverview.handles.hChannels,'ForegroundColor','blue');
                ep_editData('startChannels');
                
            case 'samples'
                set(EPoverview.handles.hSamples,'ForegroundColor','blue');
                ep_editData('startSamples');
                
            case 'freqs'
                set(EPoverview.handles.hfreqs,'ForegroundColor','blue');
                ep_editData('startFreqs');
                
            case 'events'
                set(EPoverview.handles.hEvents,'ForegroundColor','blue');
                ep_editData('startEvents');
                
            case 'factors'
                set(EPoverview.handles.hFactors,'ForegroundColor','blue');
                ep_editData('startFactors');
                
            case 'QC'
                set(EPoverview.handles.hQC,'ForegroundColor','blue');
                ep_editData('startQC');
                
            case 'PCA'
                set(EPoverview.handles.hPCA,'ForegroundColor','blue');
                ep_editData('startPCA');
                
            otherwise
                error('Not a valid mode.');
        end;
        
    case 'startOverview'
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Data Name'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-120 150 20]);
        EPoverview.handles.overview.hDataName = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%s', EPoverview.workData.dataName),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-140 150 20],'Callback', ['ep_editData(''renameDataName'');']);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Experiment Name'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-170 150 20]);
        EPoverview.handles.overview.hExpName = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%s', EPoverview.workData.ename),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-190 150 20],'Callback', ['ep_editData(''renameExp'');']);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Montage'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-220 150 20]);
        EPoverview.handles.overview.hMontage = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%s', EPoverview.workData.montage),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-240 150 20],'Callback', ['ep_editData(''montage'');']);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'File Format'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-270 150 20]);
        EPoverview.handles.overview.hFileFormat = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', EPoverview.workData.fileFormat),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-290 150 20]);
        
        if strcmp(EPoverview.workData.timeUnits,'per')
            uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of flexible samples'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-320 150 20]);
            EPoverview.handles.overview.hFs = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(EPoverview.workData.Fs/10)),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-340 150 20],'Callback', ['ep_editData(''Fs'');'],'TooltipString','Changes recorded sampling rate without changing the data by just changing the labels.');
            
            uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Prestimulus Period (percentage)'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-370 150 20]);
            EPoverview.handles.overview.hBaseline = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', (100/EPoverview.workData.Fs)*EPoverview.workData.baseline),'FontSize',EPmain.fontsize,...
                'TooltipString','Length of prestimulus recording.  Negative denotes stimulus was prior to start of the epoch.','Position',[100 windowHeight-390 150 20],'Callback', ['ep_editData(''baseline'');']);
        else
            uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Sampling Rate (Hz)'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-320 150 20]);
            EPoverview.handles.overview.hFs = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPoverview.workData.Fs),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-340 150 20],'Callback', ['ep_editData(''Fs'');'],'TooltipString','Changes recorded sampling rate without changing the data by just changing the labels.');
            
            uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Prestimulus Period (msec)'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-370 150 20]);
            EPoverview.handles.overview.hBaseline = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', (1000/EPoverview.workData.Fs)*EPoverview.workData.baseline),'FontSize',EPmain.fontsize,...
                'TooltipString','Length of prestimulus recording.  Negative denotes stimulus was prior to start of the epoch.','Position',[100 windowHeight-390 150 20],'Callback', ['ep_editData(''baseline'');']);
        end;
        
        [pathstr, name, ext] = fileparts(EPoverview.workData.fileName);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Original File Name'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-420 150 20]);
        EPoverview.handles.overview.hFileName = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', name),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-440 150 20]);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of Channels'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-470 150 20]);
        EPoverview.handles.overview.hChannels= uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.chanNames)),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-490 150 20]);
        
        EPoverview.handles.overview.hSamplesLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of Samples'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-520 150 20]);
        EPoverview.handles.overview.hSamples = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.timeNames)),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-540 150 20]);
        if isempty(EPoverview.workData.timeNames)
            set(EPoverview.handles.overview.hSamples,'enable','off');
            set(EPoverview.handles.overview.hSamplesLabel,'enable','off');
        end;
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of Cells'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-570 150 20]);
        EPoverview.handles.overview.hCells = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%d', length(unique(EPoverview.workData.cellNames))),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-590 150 20]);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of Subjects'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-520 150 20]);
        EPoverview.handles.overview.hSubjects = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.subNames)),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-540 150 20]);
        
        EPoverview.handles.overview.hFreqsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of Hz bands'),'FontSize',EPmain.fontsize,...
                'Position',[300 windowHeight-470 150 20]);
        EPoverview.handles.overview.hFreqs = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.freqNames)),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-490 150 20]);
        if isempty(EPoverview.workData.freqNames)
            set(EPoverview.handles.overview.hFreqs,'enable','off');
            set(EPoverview.handles.overview.hFreqsLabel,'enable','off');
        end;

        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Channel Coordinates'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-120 150 20]);
        EPoverview.handles.overview.hCed = uicontrol('Style','pushbutton','HorizontalAlignment','left','String', sprintf('%s', EPoverview.workData.ced),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-140 150 20],'Callback', ['ep_editData(''ced'');']);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Data Type'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-170 150 20]);
        EPoverview.handles.overview.hDataType = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', EPoverview.workData.dataType),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-190 150 20]);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Reference Type'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-220 150 20]);
        EPoverview.handles.overview.hRefType = uicontrol('Style','popupmenu','HorizontalAlignment','left','String', refLabels,'Value',find(strcmp(EPoverview.workData.reference.type,refTypes)),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-240 150 20],'Callback', ['ep_editData(''refType'');']);
        
        refchan1=[];
        refchan2=[];
        if length(EPoverview.workData.reference.original)>0
            refchan1=EPoverview.workData.reference.original(1);
        end;
        if length(EPoverview.workData.reference.original)>1
            refchan2=EPoverview.workData.reference.original(2);
        end;
        uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Original Reference Chans'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-270 150 20]);
        EPoverview.handles.overview.hOrigRefChans1 = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', refchan1),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-290 50 20],'Callback', ['ep_editData(''origRefChans1'');']);
        EPoverview.handles.overview.hOrigRefChans2 = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', refchan2),'FontSize',EPmain.fontsize,...
            'Position',[370 windowHeight-290 50 20],'Callback', ['ep_editData(''origRefChans2'');']);
        
        if strcmp(EPoverview.workData.reference.type,'REG')
            refchan1=[];
            refchan2=[];
            if length(EPoverview.workData.reference.current)>0
                refchan1=EPoverview.workData.reference.current(1);
            end;
            if length(EPoverview.workData.reference.current)>1
                refchan2=EPoverview.workData.reference.current(2);
            end;
            uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Current Reference Chans'),'FontSize',EPmain.fontsize,...
                'Position',[300 windowHeight-320 150 20]);
            EPoverview.handles.overview.hCurrRefChans1 = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', refchan1),'FontSize',EPmain.fontsize,...
                'Position',[300 windowHeight-340 50 20],'Callback', ['ep_editData(''currRefChans1'');']);
            EPoverview.handles.overview.hCurrRefChans2 = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', refchan2),'FontSize',EPmain.fontsize,...
                'Position',[370 windowHeight-340 50 20],'Callback', ['ep_editData(''currRefChans2'');']);
        end;
        
        EPoverview.handles.overview.hNumFacsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of Factors'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-370 150 20]);
        EPoverview.handles.overview.hNumFacs = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.facNames)),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-390 150 20]);
        if isempty(EPoverview.workData.facNames)
            set(EPoverview.handles.overview.hNumFacs,'enable','off');
            set(EPoverview.handles.overview.hNumFacsLabel,'enable','off');
        end;
        
        EPoverview.handles.overview.hNumTrialsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of Trials'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-420 150 20]);
        EPoverview.handles.overview.hNumTrials = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.cellNames)),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-440 150 20]);
        if ~strcmp(EPoverview.workData.dataType,'single_trial')
            set(EPoverview.handles.overview.hNumTrials,'enable','off');
            set(EPoverview.handles.overview.hNumTrialsLabel,'enable','off');
        end;
        
        EPoverview.handles.overview.hNumTrialsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Type of Data'),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-570 150 20]);
        
        if numFreqs > 1
            if numRels > 1
                if numPoints > 1
                    dataType='Phase Locking';
                else
                    dataType='Coherence';
                end;
            else
                if numPoints > 1
                    dataType='Wavelet';
                else
                    dataType='Spectral';
                end;
            end;
        else
            dataType='voltage';
        end;
                
        EPoverview.handles.overview.hNumTrials = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', dataType),'FontSize',EPmain.fontsize,...
            'Position',[300 windowHeight-590 150 20]);
        
        
    case 'startSubjects'
        
        for i=1:numSubs
            tableData{i,1}=i;
            tableData{i,2}=EPoverview.subjects.select(i);
            if ~strcmp(EPoverview.workData.dataType,'single_trial')
                tableData{i,3}=EPoverview.subjects.weights(i);
            end;
            tableData{i,4}=char(EPoverview.workData.subNames(i));
            tableData{i,5}=char(EPoverview.workData.subTypes(i));
            if numSpecs
                tableData(i,6:numSpecs+5)=EPoverview.workData.subjectSpecs(i,:);
            end;
        end;
        
        tableNames{1}='order';
        tableNames{2}='select';
        if ~strcmp(EPoverview.workData.dataType,'single_trial')
            tableNames{3}='weights';
        end;
        tableNames{4}='names';
        tableNames{5}='type';
        tableNames(6:length(EPoverview.workData.subjectSpecNames)+5)=EPoverview.workData.subjectSpecNames;
        numNames=size(tableNames,2);
        
        columnEditable =  repmat(true,1,length(EPoverview.workData.subjectSpecNames)+5);
        if strcmp(EPoverview.workData.dataType,'single_trial')
            columnEditable(3)=false;
        end;
        ColumnFormat{1}='numeric';
        ColumnFormat{2}='logical';
        ColumnFormat{3}='numeric';
        ColumnFormat{4}=[];
        ColumnFormat{5}={'RAW','AVG','GAV'};
        
        EPoverview.handles.subjects.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'RearrangeableColumns','on',...
            'CellEditCallback','ep_editData(''editSubs'');','Position',[100 windowHeight-600 500 min(500,windowHeight-160)]);
        
        
        EPoverview.handles.subjects.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Delete','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-250 60 30], 'Callback', 'ep_editData(''deleteSubs'');');
        if ~strcmp(EPoverview.workData.dataType,'single_trial')
            EPoverview.handles.subjects.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Add','FontSize',EPmain.fontsize,...
                'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''addSubs'');');
        end;
        
        tableData=get(EPoverview.handles.subjects.hTable,'Data');
        weights=cell2mat(tableData(:,3));
        
        EPoverview.handles.subjects.hWeightsSum = uicontrol('Style', 'text', 'String', sprintf('%f',sum(weights)),'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-310 60 30]);
        EPoverview.handles.subjects.hClear = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-370 60 30], 'Callback', 'ep_editData(''clearSubWeights'');');
        EPoverview.handles.subjects.hAll = uicontrol('Style', 'pushbutton', 'String', 'All','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-400 60 30], 'Callback', 'ep_editData(''allSubs'');');
        EPoverview.handles.subjects.hAppend = uicontrol('Style', 'pushbutton', 'String', 'Append','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-430 60 30], 'Callback', 'ep_editData(''appendSubs'');');
        EPoverview.handles.subjects.hImport = uicontrol('Style', 'pushbutton', 'String', 'Import','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''importSubs'');');
        EPoverview.handles.subjects.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-490 60 30], 'Callback', 'ep_editData(''exportSubs'');');
        EPoverview.handles.subjects.addSpec = uicontrol('Style', 'pushbutton', 'String', '+ Spec','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-520 60 30], 'Callback', 'ep_editData(''addSubSpec'');');
        EPoverview.handles.subjects.minusSpec = uicontrol('Style', 'pushbutton', 'String', '- Spec','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-550 60 30], 'Callback', 'ep_editData(''minusSubSpec'');');
        if isempty(EPoverview.workData.subjectSpecNames)
            set(EPoverview.handles.subjects.minusSpec,'enable','off');
        end;
        
        sortString=[EPoverview.workData.subjectSpecNames;'Names';'Sort'];
        EPoverview.handles.subjects.sort = uicontrol('Style', 'popupmenu', 'Value', length(sortString), 'String', sortString, 'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-580 80 30], 'Callback', 'ep_editData(''sortSubjects'');');
        
    case 'startCells'
        
        [uniqueCells waveOrder m]=unique(EPoverview.workData.cellNames,'stable');
        [n xout] = hist(m,length(uniqueCells));
        
        for iCell=1:numCells %where iCell is the display order
            tableData{iCell,1}=iCell;
            tableData{iCell,2}=EPoverview.cells.select(iCell);
            if ~strcmp(EPoverview.workData.dataType,'single_trial')
                tableData{iCell,3}=EPoverview.cells.weights(iCell);
            end;
            tableData{iCell,4}=char(uniqueCells{iCell});
            if strcmp(EPoverview.workData.dataType,'single_trial')
                tableData{iCell,5}=n(iCell);
            else
                tableData{iCell,5}=char(EPoverview.workData.cellTypes(iCell));
                %if length(EPoverview.workData.subNames)==1
                    %if ~isempty(EPoverview.workData.trialSpecs)
                        for iSpec=1:length(EPoverview.workData.trialSpecNames)
                            if any(cellfun(@isnumeric,EPoverview.workData.trialSpecs(iCell,iSpec,:)))
                                theNums=[];
                                for iSub=1:numSubs
                                    if isnumeric(EPoverview.workData.trialSpecs{iCell,iSpec,iSub}) && ~isempty(EPoverview.workData.trialSpecs{iCell,iSpec,iSub})&& ~isnan(EPoverview.workData.trialSpecs{iCell,iSpec,iSub})
                                        theNums(end+1,1)=EPoverview.workData.trialSpecs{iCell,iSpec,iSub};
                                    end;
                                end;
                                tableData{iCell,iSpec+5}=mean(theNums);
                            elseif all(strcmp(EPoverview.workData.trialSpecs(iCell,iSpec,1),EPoverview.workData.trialSpecs(iCell,iSpec,:)))
                                tableData(iCell,iSpec+5)=EPoverview.workData.trialSpecs(iCell,iSpec,1);
                            else
                                tableData{iCell,iSpec+5}='';
                            end;
                        end;
                    %end;
                %end;
            end;
        end;
        
        tableNames{1}='order';
        tableNames{2}='select';
        if ~strcmp(EPoverview.workData.dataType,'single_trial')
            tableNames{3}='weights';
        end;
        tableNames{4}='names';
        if strcmp(EPoverview.workData.dataType,'single_trial')
            tableNames{5}='trials';
        else
            tableNames{5}='types';
            %if length(EPoverview.workData.subNames)==1
                tableNames(6:length(EPoverview.workData.trialSpecNames)+5)=EPoverview.workData.trialSpecNames;
            %end;
        end
        
        columnEditable =  repmat(true,1,5);
        ColumnFormat{1}='numeric';
        ColumnFormat{2}='logical';
        ColumnFormat{3}='numeric';
        ColumnFormat{4}=[];
        if strcmp(EPoverview.workData.dataType,'single_trial')
            ColumnFormat{5}=[];
            columnEditable(1)=false;
            columnEditable(5)=false;
        else
            ColumnFormat{5}={'SGL','CMB','STS'};
            columnEditable(5)=true;
            ColumnFormat(6:length(EPoverview.workData.trialSpecNames)+5)=cell(1,length(EPoverview.workData.trialSpecNames));
            if length(EPoverview.workData.subNames)==1
                columnEditable(6:length(EPoverview.workData.trialSpecNames)+5)=true;
            else
                columnEditable(6:length(EPoverview.workData.trialSpecNames)+5)=false;
            end;
        end
        
        EPoverview.handles.cells.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'RearrangeableColumns','on',...
            'CellEditCallback',@editCells,'Position',[100 windowHeight-600 500 min(500,windowHeight-160)]);
        
        EPoverview.handles.cells.hDelete = uicontrol('Style', 'pushbutton', 'String', 'Delete','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-250 60 30], 'Callback', 'ep_editData(''deleteCells'');');
        if ~strcmp(EPoverview.workData.dataType,'single_trial')
            EPoverview.handles.cells.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Add','FontSize',EPmain.fontsize,...
                'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''addCells'');');
            
            tableData=get(EPoverview.handles.cells.hTable,'Data');
            weights=cell2mat(tableData(:,3));
            
            EPoverview.handles.cells.hWeights = uicontrol('Style', 'text', 'String', sprintf('%f',sum(weights)),'FontSize',EPmain.fontsize,...
                'Position', [20 windowHeight-310 60 30]);
            EPoverview.handles.cells.hTrialWeights = uicontrol('Style', 'togglebutton', 'String', 'TrlWgts', 'Value',EPoverview.cells.TrlWgts,'FontSize',EPmain.fontsize,...
                'Position', [20 windowHeight-340 60 30], 'Callback', 'EPoverview.cells.TrlWgts=get(EPoverview.handles.cells.hTrialWeights,''Value'');','FontSize',EPmain.fontsize,...
                'TooltipString','Weight by number of trials in average if available.');
        else
            EPoverview.handles.cells.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Average','FontSize',EPmain.fontsize,...
                'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''addTrials'');');
        end;
        EPoverview.handles.cells.hDelete = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-370 60 30], 'Callback', 'ep_editData(''clearCellWeights'');');
        EPoverview.handles.cells.hAll = uicontrol('Style', 'pushbutton', 'String', 'All','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-400 60 30], 'Callback', 'ep_editData(''allCells'');');
        EPoverview.handles.cells.hAppend = uicontrol('Style', 'pushbutton', 'String', 'Append','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-430 60 30], 'Callback', 'ep_editData(''appendCells'');');
        EPoverview.handles.cells.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''exportCells'');');
        
        if strcmp(EPoverview.workData.dataType,'average') && (length(EPoverview.workData.subNames)==1) && ~isempty(EPoverview.workData.trialSpecNames)
            sortString=[EPoverview.workData.trialSpecNames;'Names';'Sort'];
        else
            sortString={'Names';'Sort'};
        end;
        EPoverview.handles.cells.sort = uicontrol('Style', 'popupmenu', 'Value', length(sortString), 'String', sortString, 'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-580 80 30], 'Callback', 'ep_editData(''sortCells'');');
        
    case 'startTrials'
        
        tableData(:,1)=EPoverview.workData.cellNames;
        tableData(:,2)=num2cell(EPoverview.workData.trialNames);
        tableData(:,3)=num2cell(EPoverview.workData.recTime);
        tableData(:,4:length(EPoverview.workData.trialSpecNames)+3)=EPoverview.workData.trialSpecs;
        
        
        tableNames{1}='cell';
        tableNames{2}='trial';
        tableNames{3}='time';
        tableNames(4:length(EPoverview.workData.trialSpecNames)+3)=EPoverview.workData.trialSpecNames;
        
        
        columneditable =  repmat(true,1,length(EPoverview.workData.trialSpecNames));
        
        EPoverview.handles.trials.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columneditable,...
            'RearrangeableColumns','on',...
            'CellEditCallback','ep_editData(''trialSpecs'');','Position',[100 windowHeight-600 500 min(500,windowHeight-160)]);
                
        EPoverview.handles.trials.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''exportTrials'');');
        
        EPoverview.handles.trials.hExport = uicontrol('Style', 'pushbutton', 'String', 'Import','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-490 60 30], 'Callback', 'ep_editData(''importTrials'');',...
            'TooltipString','Load in text file with new  cell names for trials');
        
        sortString=[EPoverview.workData.trialSpecNames;'Cells';'Trials';'Sort'];
        EPoverview.handles.trials.sort = uicontrol('Style', 'popupmenu', 'Value', length(sortString), 'String', sortString, 'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-580 80 30], 'Callback', 'ep_editData(''sortTrials'');');
        
    case 'startChannels'
        
        for i=1:numChans
            tableData{i,1}=i;
            tableData{i,2}=EPoverview.channels.select(i);
            tableData{i,3}=EPoverview.channels.weights(i);
            tableData{i,4}=char(EPoverview.workData.chanNames{i});
            tableData{i,5}=char(EPoverview.workData.chanTypes{i});
        end;
        
        for i=1:numImplicit
            tableData{i+numChans,1}=i+numChans;
            tableData{i+numChans,2}=EPoverview.channels.select(i+numChans);
            tableData{i+numChans,3}=EPoverview.channels.weights(i+numChans);
            tableData{i+numChans,4}=char(EPoverview.workData.implicit(i).labels);
            tableData{i+numChans,5}=char(EPoverview.workData.implicit(i).type);
        end;
        
        tableNames{1}='order';
        tableNames{2}='select';
        tableNames{3}='weights';
        tableNames{4}='Name';
        tableNames{5}='Type';
        
        columnEditable =  [repmat(true,1,5)];
        ColumnFormat{1}='numeric';
        ColumnFormat{2}='logical';
        ColumnFormat{3}='numeric';
        ColumnFormat{4}=[];
        ColumnFormat{5}={'EEG','ANS','MGM','MGA','MGP','ECG','FID','REG'};
        
        if ~isempty(EPoverview.workData.impedances.channels)
            tableNames{6}='impedance';
            ColumnFormat{6}='numeric';
            columnEditable(6)=true;
            for i=1:numChans
                tableData{i,6}=mean(EPoverview.workData.impedances.channels(i,:),2);
            end;
        end;
        
        EPoverview.handles.channels.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'RearrangeableColumns','on',...
            'CellEditCallback','ep_editData(''editChannels'');','Position',[100 windowHeight-600 500 min(500,windowHeight-160)]);
        
        EPoverview.handles.channels.hDelete = uicontrol('Style', 'pushbutton', 'String', 'Delete','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-250 60 30], 'Callback', 'ep_editData(''deleteChans'');');
        EPoverview.handles.channels.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Add','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''addChans'');');
        
        tableData=get(EPoverview.handles.channels.hTable,'Data');
        weights=cell2mat(tableData(:,3));
        
        EPoverview.handles.channels.hWeights = uicontrol('Style', 'text', 'String', sprintf('%f',sum(weights)),'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-310 60 30]);
        EPoverview.handles.channels.hClear = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-370 60 30], 'Callback', 'ep_editData(''clearChanWeights'');');
        EPoverview.handles.channels.hAll = uicontrol('Style', 'pushbutton', 'String', 'All','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-400 60 30], 'Callback', 'ep_editData(''allChans'');');
        EPoverview.handles.channels.hAppend = uicontrol('Style', 'pushbutton', 'String', 'Append','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-430 60 30], 'Callback', 'ep_editData(''appendChans'');');
        EPoverview.handles.channels.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''exportChans'');');
        
    case 'startSamples'
        
        if strcmp(EPoverview.workData.timeUnits,'per')
            uicontrol('Style', 'text', 'String', '%','FontSize',EPmain.fontsize,...
                'Position', [100 windowHeight-120 60 15],'HorizontalAlignment','left');
        else
            uicontrol('Style', 'text', 'String', 'Ms','FontSize',EPmain.fontsize,...
                'Position', [100 windowHeight-120 60 15],'HorizontalAlignment','left');
        end;
        uicontrol('Style', 'text', 'String', 'Samples','FontSize',EPmain.fontsize,...
            'Position', [100 windowHeight-145 60 15],'HorizontalAlignment','left');
        
        theNumber=EPoverview.workData.timeNames(1);
        if ceil(theNumber)==theNumber
            theString=sprintf('%d', theNumber);
        else
            theString=sprintf('%f', theNumber);
        end;        
        EPoverview.handles.samples.startMs = uicontrol('Style','edit','HorizontalAlignment','left','String', theString,'FontSize',EPmain.fontsize,...
            'Position',[200 windowHeight-120 100 25],'Callback', ['ep_editData(''changeStartMs'');'],'TooltipString','msec of first sample onset');
        
        theNumber=EPoverview.workData.timeNames(end)+(1000/EPoverview.workData.Fs);
        if ceil(theNumber)==theNumber
            theString=sprintf('%d', theNumber);
        else
            theString=sprintf('%f', theNumber);
        end;        
        EPoverview.handles.samples.endMs = uicontrol('Style','edit','HorizontalAlignment','left','String', theString,'FontSize',EPmain.fontsize,...
            'Position',[350 windowHeight-120 100 25],'Callback', ['ep_editData(''changeEndMs'');'],'TooltipString','msec of last sample offset');
        
        EPoverview.handles.samples.startSamp = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', 1),'FontSize',EPmain.fontsize,...
            'Position',[200 windowHeight-145 100 25],'Callback', ['ep_editData(''changeStartSamp'');']);
        EPoverview.handles.samples.endSamp = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.timeNames)),'FontSize',EPmain.fontsize,...
            'Position',[350 windowHeight-145 100 25],'Callback', ['ep_editData(''changeEndSamp'');']);
        
        if strcmp(EPoverview.workData.timeUnits,'per')
            EPoverview.handles.samples.hFsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Number of flexible samples'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-270 150 20]);
            EPoverview.handles.samples.hFs = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(EPoverview.workData.Fs/10)),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-290 150 20],'Callback', ['ep_editData(''Fs'');'],'TooltipString','Changes recorded sampling rate without changing the data by just changing the labels.');
            
            EPoverview.handles.samples.hBaselineLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Prestimulus period (percentage)'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-320 150 20]);
            
            EPoverview.handles.samples.hBaseline = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', (100/EPoverview.workData.Fs)*EPoverview.workData.baseline),'FontSize',EPmain.fontsize,...
                'TooltipString','Length of prestimulus recording.  Negative denotes stimulus was prior to start of the epoch.','Position',[100 windowHeight-340 150 20],'Callback', ['ep_editData(''baseline'');']);
        else
            EPoverview.handles.samples.hFsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Sampling Rate (Hz)'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-270 150 20]);
            EPoverview.handles.samples.hFs = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPoverview.workData.Fs),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-290 150 20],'Callback', ['ep_editData(''Fs'');'],'TooltipString','Changes recorded sampling rate without changing the data by just changing the labels.');
            
            EPoverview.handles.samples.hBaselineLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Prestimulus period (msec)'),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-320 150 20]);
            
            EPoverview.handles.samples.hBaseline = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', (1000/EPoverview.workData.Fs)*EPoverview.workData.baseline),'FontSize',EPmain.fontsize,...
                'TooltipString','Length of prestimulus recording.  Negative denotes stimulus was prior to start of the epoch.','Position',[100 windowHeight-340 150 20],'Callback', ['ep_editData(''baseline'');']);
        end;
        
        EPoverview.handles.samples.halfRate = uicontrol('Style', 'pushbutton', 'String', '1/2 Rate','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''halfRate'');','TooltipString','Halves the effective sampling rate by merging adjoining time points.');
        EPoverview.handles.samples.doubleRate = uicontrol('Style', 'pushbutton', 'String', '2x Rate','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-310 60 30], 'Callback', 'ep_editData(''twiceRate'');','TooltipString','Doubles the effective sampling rate by interpolating between time points.');        
        
    case 'startFreqs'
        
        sampleTime=1000/EPoverview.workData.Fs;
        
        uicontrol('Style', 'text', 'String', 'Hz','FontSize',EPmain.fontsize,...
            'Position', [100 windowHeight-120 60 15],'HorizontalAlignment','left');
        uicontrol('Style', 'text', 'String', 'Bins','FontSize',EPmain.fontsize,...
            'Position', [100 windowHeight-145 60 15],'HorizontalAlignment','left');
        
        
        EPoverview.handles.freqs.startFreq = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPoverview.workData.freqNames(1)),'FontSize',EPmain.fontsize,...
            'Position',[200 windowHeight-120 100 25],'Callback', ['ep_editData(''changeStartFreq'');']);
        EPoverview.handles.freqs.endFreq = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPoverview.workData.freqNames(end)),'FontSize',EPmain.fontsize,...
            'Position',[350 windowHeight-120 100 25],'Callback', ['ep_editData(''changeEndFreq'');']);        
        
        EPoverview.handles.freqs.startBin = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', 1),'FontSize',EPmain.fontsize,...
            'Position',[200 windowHeight-145 100 25],'Callback', ['ep_editData(''changeStartBin'');']);
        EPoverview.handles.freqs.endBin = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', length(EPoverview.workData.freqNames)),'FontSize',EPmain.fontsize,...
            'Position',[350 windowHeight-145 100 25],'Callback', ['ep_editData(''changeEndBin'');']);        
        
        EPoverview.handles.freqs.hBinningLabel = uicontrol('Style','text','HorizontalAlignment','left','String', sprintf('%s', 'Binning (Hz)'),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-270 150 20]);
        EPoverview.handles.freqs.hBinning = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPoverview.workData.freqNames(2)-EPoverview.workData.freqNames(1)),'FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-290 150 20],'Callback', ['ep_editData(''Binning'');'],'TooltipString','Changes recorded frequency bin size without changing the data by just changing the labels.');
        
        EPoverview.handles.freqs.halfBinning = uicontrol('Style', 'pushbutton', 'String', '1/2 Binning','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''halfBinning'');','TooltipString','Halves the effective frequency binning by merging adjoining frequency bins.');
        
    case 'startEvents'
        
        keyNames=cell(0);
        tableData=[];
        eventCounter=0;
        if strcmp(EPoverview.workData.dataType,'single_trial')
            trialCol=1;
        else
            trialCol=0;
        end
        for subject=1:numSubs
            for wave=1:numWaves
                for event=1:length(EPoverview.workData.events{subject,wave})
                    eventCounter=eventCounter+1;
                    tableData{eventCounter,1}=EPoverview.events.select(eventCounter);
                    tableData{eventCounter,2}=EPoverview.workData.cellNames{wave};
                    if strcmp(EPoverview.workData.dataType,'single_trial')
                        tableData{eventCounter,3}=EPoverview.workData.trialNames(wave);
                    end;
                    tableData{eventCounter,3+trialCol}=EPoverview.workData.subNames{subject};
                    tableData{eventCounter,4+trialCol}=event;
                    tableData{eventCounter,5+trialCol}=EPoverview.workData.events{subject,wave}(event).type;
                    tableData{eventCounter,6+trialCol}=EPoverview.workData.events{subject,wave}(event).sample;
                    tableData{eventCounter,7+trialCol}=EPoverview.workData.events{subject,wave}(event).value;
                    tableData{eventCounter,8+trialCol}=(EPoverview.workData.events{subject,wave}(event).sample-EPoverview.workData.baseline-1)*(1000/EPoverview.workData.Fs);
                    tableData{eventCounter,9+trialCol}=EPoverview.workData.events{subject,wave}(event).duration;
                    for iKey=1:length(EPoverview.workData.events{subject,wave}(event).keys)
                        theKeyName=EPoverview.workData.events{subject,wave}(event).keys(iKey).code;
                        if ~isempty(theKeyName)
                            if ~any(strcmp(theKeyName,keyNames))
                                keyNames{end+1}=theKeyName;
                            end;
                            tableData{eventCounter,9+trialCol+find(strcmp(theKeyName,keyNames))}=EPoverview.workData.events{subject,wave}(event).keys(iKey).data;
                        end;
                    end;
                end;
            end;
        end;
        
        tableNames{1}='select';
        tableNames{2}='cells';
        nonEdit=3;
        if strcmp(EPoverview.workData.dataType,'single_trial')
            tableNames{3}='trials';
            nonEdit=4;
        end;
        tableNames{3+trialCol}='subjects';
        tableNames{4+trialCol}='events';
        tableNames{5+trialCol}='type';
        tableNames{6+trialCol}='sample';
        tableNames{7+trialCol}='value';
        tableNames{8+trialCol}='latency';
        tableNames{9+trialCol}='duration';
        numNonKey=length(tableNames);
        for i=1:length(keyNames)
            tableNames(9+trialCol+i)=keyNames(i);
        end;
        
        columnEditable = [true repmat(false,1,nonEdit) repmat(true,1,5) repmat(true,1,length(keyNames))];
%         if ~strcmp(EPoverview.workData.dataType,'single_trial')
%             columnEditable(3)=false;
%         end;
        
        ColumnFormat = cell(1,numNonKey);
        ColumnFormat{1}='logical';
        for i=2:numNonKey+length(keyNames)
            ColumnFormat{i}=[];
        end;
        
        EPoverview.handles.trials.hImport = uicontrol('Style', 'pushbutton', 'String', 'Import','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-490 60 30], 'Callback', 'ep_editData(''importEvents'');');

        if ~isempty(tableData)
            try
                EPoverview.handles.events.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                    'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                    'RearrangeableColumns','on',...
                    'CellEditCallback','ep_editData(''editEvents'');','Position',[100 windowHeight-600 500 min(500,windowHeight-160)]);
                EPoverview.handles.events.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Delete','FontSize',EPmain.fontsize,...
                    'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''deleteEvents'');');
                EPoverview.handles.events.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
                    'Position', [20 windowHeight-370 60 30], 'Callback', 'ep_editData(''clearEventWeights'');');
                EPoverview.handles.events.hAdd = uicontrol('Style', 'pushbutton', 'String', 'All','FontSize',EPmain.fontsize,...
                    'Position', [20 windowHeight-400 60 30], 'Callback', 'ep_editData(''allEvents'');');
                EPoverview.handles.trials.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
                    'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''exportEvents'');');
            catch
                EPoverview.handles.events.hTable = uicontrol('Style','text','HorizontalAlignment','left','String', 'This function does not work with this version of Matlab','FontSize',EPmain.fontsize,...
                    'Position',[130 windowHeight-200 200 30]);
            end;
        else
            EPoverview.handles.events.hTable = uicontrol('Style','text','HorizontalAlignment','left','String', 'No events in this data.','FontSize',EPmain.fontsize,...
                'Position',[130 windowHeight-200 200 30]);
        end;
        
    case 'startFactors'
        
        %[EPdata]=ep_stripAdds(EPoverview.workData,{'SGLfac','CMBfac'});
        EPdata=EPoverview.workData;
        
        numPoints=max(length(EPdata.timeNames),1);
        numChans=length(EPdata.chanNames);
        numFreqs=length(EPdata.freqNames);
        
        if ~isempty(EPdata.facVecT)
            [C peakLatency]=max(abs(EPdata.facVecT));
        else
            peakLatency=zeros(1,numSGLfacs);
            for theFactor=1:numSGLfacs
                [C peakLatency(theFactor)]=max(max(reshape(shiftdim(abs(mean(EPdata.data(:,:,:,:,theFactor,:,:),4)),1),numPoints,[])'));
            end;
        end;
        peakLatency=[peakLatency repmat(NaN,1,numCMBfacs)];
        
        peakNegChan=zeros(1,numSGLfacs);
        peakPosChan=zeros(1,numSGLfacs);
        peakPolarity=zeros(1,numSGLfacs);
        for theFactor=1:numSGLfacs
            if ~isempty(EPdata.facVecS)
                maxVal=max(max(max(max(max(mean(EPdata.data(:,:,:,:,theFactor,:,:),4))))));
                minVal=min(min(min(min(min(mean(EPdata.data(:,:,:,:,theFactor,:,:),4))))));
                if abs(maxVal) > abs(minVal)
                    theVal=maxVal;
                else
                    theVal=minVal;
                end;
                [minChanVal minChan]=min(EPdata.facVecS(:,theFactor));
                [maxChanVal maxChan]=max(EPdata.facVecS(:,theFactor));
                minChanVal=minChanVal*theVal;
                maxChanVal=maxChanVal*theVal;
                if sign(maxChanVal) == 1
                    peakPosChan(theFactor) = maxChan;
                    if sign(minChanVal) == -1
                        peakNegChan(theFactor) = minChan;
                    end;
                    if abs(maxChanVal) > abs(minChanVal)
                        peakPolarity(theFactor)=1;
                    else
                        peakPolarity(theFactor)=-1;
                    end;
                else
                    peakNegChan(theFactor) = maxChan;
                    if sign(minChanVal) == 1
                        peakPosChan(theFactor) = minChan;
                    end;
                    if abs(maxChanVal) > abs(minChanVal)
                        peakPolarity(theFactor)=-1;
                    else
                        peakPolarity(theFactor)=1;
                    end;
                end;
            else
                [C peakNegChan(theFactor)]=min(min(reshape(mean(EPdata.data(:,:,:,:,theFactor,:,:),4),numChans,[])'));
                [C peakPosChan(theFactor)]=max(max(reshape(mean(EPdata.data(:,:,:,:,theFactor,:,:),4),numChans,[])'));
            end;
            if ~isempty(EPdata.facVecT)
                [A maxPnt]=max(abs(EPdata.facVecT(:,theFactor)));
                if sign(maxPnt) == -1
                    t1=peakNegChan(theFactor);
                    t2=peakPosChan(theFactor);
                    peakNegChan(theFactor)=t2;
                    peakPosChan(theFactor)=t1;
                    peakPolarity(theFactor)=-peakPolarity(theFactor);
                end;
            end;
            if ~isempty(EPdata.facVecF)
                [A maxPnt]=max(abs(EPdata.facVecF(:,theFactor)));
                if sign(maxPnt) == -1
                    t1=peakNegChan(theFactor);
                    t2=peakPosChan(theFactor);
                    peakNegChan(theFactor)=t2;
                    peakPosChan(theFactor)=t1;
                    peakPolarity(theFactor)=-peakPolarity(theFactor);
                end;
            end;
        end;
            
        peakNegChan=[peakNegChan repmat(NaN,1,numCMBfacs)];
        peakPosChan=[peakPosChan repmat(NaN,1,numCMBfacs)];
        
        if ~isempty(EPdata.facVecF)
            [C peakHz]=max(abs(EPdata.facVecF));
        else
            peakHz=zeros(1,numSGLfacs);
            if numFreqs
                for theFactor=1:numSGLfacs
                    [C peakHz(theFactor)]=max(max(reshape(shiftdim(abs(mean(EPdata.data(:,:,:,:,theFactor,:,:),4)),5),numFreqs,[])'));
                end;
            end;
        end;
        peakHz=[peakHz repmat(NaN,1,numCMBfacs)];
        
%         if ~isempty(EPdata.facVecT)
%             timeSign=sign(EPdata.facVecT(peakLatency(1:numSGLfacs)));
%         else
%             timeSign=repmat(1,1,numSGLfacs);
%         end;
%         
%         
%         for theFactor=1:numSGLfacs
%             theMax=max(max(reshape(shiftdim(mean(EPdata.data(:,:,:,:,theFactor,:,:),4),1),size(EPdata.data,2),[])')); %to get signs for factor data, multiply signs of loadings by scores.
%             theMin=min(min(reshape(shiftdim(mean(EPdata.data(:,:,:,:,theFactor,:,:),4),1),size(EPdata.data,2),[])')); %to get signs for factor data, multiply signs of loadings by scores.
%             if abs(theMax) >= abs(theMin)
%                 peakPolarity(theFactor)=1;
%             else
%                 peakPolarity(theFactor)=-1;
%             end;
%             peakPolarity(theFactor)=peakPolarity(theFactor)*timeSign(theFactor);
%         end;
        peakPolarity=[peakPolarity repmat(NaN,1,numCMBfacs)];
        
        peakLatency=(peakLatency-EPdata.baseline-1)*(1000/EPdata.Fs);
        if isempty(EPdata.timeNames)  
            peakLatency=NaN(size(peakLatency));
            peakPolarity=NaN(size(peakPolarity));
        end;
        if isempty(EPdata.freqNames)
            peakHz=NaN(size(peakHz));
        else
            peakHz=EPdata.freqNames(peakHz(1:numSGLfacs))';
            peakHz=peakHz(:)';
            peakHz=[peakHz NaN(1,numCMBfacs)];
        end;
        
        for iFac=1:numSGLfacs
            if peakNegChan(iFac)
                peakNegChanNames{iFac,1}=EPdata.chanNames{peakNegChan(iFac)};
            else
                peakNegChanNames{iFac,1}=' ';
            end;
            if peakPosChan(iFac)
                peakPosChanNames{iFac,1}=EPdata.chanNames{peakPosChan(iFac)};
            else
                peakPosChanNames{iFac,1}=' ';
            end;            
        end;
        
        peakNegChanNames=[peakNegChanNames; cellstr(repmat(' ',numCMBfacs,1))];
        peakPosChanNames=[peakPosChanNames; cellstr(repmat(' ',numCMBfacs,1))];
        
        facVar=zeros(numFacs,1);
        facVarQ=zeros(numFacs,1);
        singleFacs=find(ismember(EPdata.facTypes,'SGL'));
        facVar(singleFacs)=EPdata.facVar;
        facVarQ(singleFacs)=EPdata.facVarQ;
        
        for i=1:numFacs
            tableData{i,1}=i;
            tableData{i,2}=EPoverview.factors.select(i);
            tableData{i,3}=EPoverview.factors.weights(i);
            tableData{i,4}=char(EPdata.facNames(i));
            tableData{i,5}=char(EPdata.facTypes(i));
            tableData{i,6}=num2str(peakLatency(i));
            tableData{i,7}=char(peakNegChanNames(i));
            tableData{i,8}=char(peakPosChanNames(i));
            tableData{i,9}=peakPolarity(i);
            tableData{i,10}=peakHz(i);
            tableData{i,11}=facVar(i);
            tableData{i,12}=facVarQ(i);
        end;
        
        tableNames{1}='order';
        tableNames{2}='select';
        tableNames{3}='weights';
        tableNames{4}='names';
        tableNames{5}='type';
        tableNames{6}='peakLatency';
        tableNames{7}='peakNegChan';
        tableNames{8}='peakPosChan';
        tableNames{9}='peakPolarity';
        tableNames{10}='peakHz';
        tableNames{11}='variance';
        tableNames{12}='unique variance';
        
        columnEditable =  [repmat(true,1,5) repmat(false,1,4)];
        ColumnFormat{1}='numeric';
        ColumnFormat{2}='logical';
        ColumnFormat{3}='numeric';
        ColumnFormat{4}=[];
        ColumnFormat{5}={'SGL','CMB','STS'};
        ColumnFormat{6}=[];
        ColumnFormat{7}=[];
        ColumnFormat{8}=[];
        ColumnFormat{9}='+';
        ColumnFormat{10}=[];
        ColumnFormat{11}='numeric';
        ColumnFormat{12}='numeric';
        
        EPoverview.handles.factors.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'RearrangeableColumns','on',...
            'CellEditCallback','ep_editData(''editFacs'');','Position',[100 windowHeight-600 500 min(500,windowHeight-160)]);
        
        EPoverview.handles.factors.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Delete','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-280 60 30], 'Callback', 'ep_editData(''deleteFacs'');');
        EPoverview.handles.factors.hAdd = uicontrol('Style', 'pushbutton', 'String', 'Add','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-310 60 30], 'Callback', 'ep_editData(''addFacs'');');
        
        tableData=get(EPoverview.handles.factors.hTable,'Data');
        weights=cell2mat(tableData(:,3));
        
        EPoverview.handles.factors.hWeightSum = uicontrol('Style', 'text', 'String', sprintf('%f',sum(weights)),'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-340 60 30]);
        EPoverview.handles.factors.hClear = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-370 60 30], 'Callback', 'ep_editData(''clearFacWeights'');');
        EPoverview.handles.factors.hAll = uicontrol('Style', 'pushbutton', 'String', 'All','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-400 60 30], 'Callback', 'ep_editData(''allFacs'');');
        EPoverview.handles.factors.hAppend = uicontrol('Style', 'pushbutton', 'String', 'Append','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-430 60 30], 'Callback', 'ep_editData(''appendFacs'');');
        EPoverview.handles.factors.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''exportFacs'');');
        
    case 'startQC'
        
        badTrials=EPoverview.workData.analysis.badTrials;
        
        %prepare contents of the table based on the QC setting
        switch EPoverview.QC.mode
            case 'avgNum'
                QCdata=EPoverview.workData.avgNum;
                if all(QCdata == 0)
                    QCdata=nan(size(QCdata));
                end;
            case 'subNum'
                QCdata=EPoverview.workData.subNum;
            case 'blinkTrial'
                QCdata=EPoverview.workData.analysis.blinkTrial;
            case 'saccadeTrial'
                QCdata=EPoverview.workData.analysis.saccadeTrial;
            case 'saccadeOnset'
                QCdata=EPoverview.workData.analysis.saccadeOnset;
            case 'moveTrial'
                QCdata=EPoverview.workData.analysis.moveTrial;
            case 'badTrials'
                QCdata=EPoverview.workData.analysis.badTrials;
            case 'badChans'
                switch EPoverview.workData.dataType
                    case {'single_trial','continuous'}
                        theNumbers=EPoverview.workData.analysis.badChans.*repmat(~badTrials,[1,1,numChans]); %count only good trials
                        if any(EPoverview.workData.analysis.badChans < 0)
                            QCdata=-squeeze(sum((theNumbers < 0).*theNumbers,3));
                            theLabel='bad';
                        else
                            QCdata=squeeze(sum(theNumbers,3));
                            theLabel='replaced';
                        end;
                    case 'average'
                        if any(isnan(EPoverview.workData.analysis.badChans))
                            theNumbers=isnan(EPoverview.workData.analysis.badChans);
                            QCdata=squeeze(sum(theNumbers,3));
                            theLabel='bad';
                        elseif any(EPoverview.workData.analysis.badChans < 0)
                            theNumbers=EPoverview.workData.analysis.badChans;
                            QCdata=-squeeze(sum((theNumbers < 0).*theNumbers,3));
                            theLabel='dropped';
                        else
                            theNumbers=EPoverview.workData.analysis.badChans;
                            QCdata=squeeze(sum(theNumbers,3));
                            theLabel='replaced';
                        end;
                    otherwise
                        disp('Oops - data type not recognized.')
                        return
                end;
                EPoverview.handles.QC.badChans = uicontrol('Style', 'text', 'String', theLabel,'FontSize',EPmain.fontsize,...
                'TooltipString','Type of bad channel info.',...
                'Position', [520 windowHeight-680 60 30]);
            case 'noise'
                QCdata=sqrt(squeeze(mean(mean(mean(EPoverview.workData.noise.^2,1),2),5)))'; %root mean square
            case 'std'
                QCdata=squeeze(mean(mean(mean(mean(abs(EPoverview.workData.std),1),2),5),6))';
        end;
        
        avgNum=EPoverview.workData.avgNum;
        subNum=EPoverview.workData.subNum;
        if all(avgNum(:) == 0) && all(subNum(:) == 0)
            avgNum=nan(size(avgNum)); %if number of trials going into each waveform is unavailable, default to nan.
        elseif all(avgNum(:) == 0) && ~all(subNum(:) == 0)
            avgNum=subNum; %if number of trials going into each waveform is unavailable but subject numbers are available (as in grand averages) then default to subject numbers.
        end;
        if strcmp(EPoverview.workData.dataType,'average')
            avgNum=avgNum+badTrials;
        end;
        
        %convert QC statistics to %if appropriate
        switch EPoverview.QC.mode
            case 'blinkTrial'
                QCdata=QCdata./(avgNum);
            case 'saccadeTrial'
                QCdata=QCdata./(avgNum);
            case 'moveTrial'
                QCdata=QCdata./(avgNum);
            case 'badTrials'
                QCdata=QCdata./(avgNum);
            case 'badChans'
                QCdata=QCdata./(numChans*(avgNum));
        end;
        
        [uniqueCells waveOrder m]=unique(EPoverview.workData.cellNames,'stable');
        
        %aggregate over trials for single trial data
        if any(strcmp(EPoverview.workData.dataType,{'single_trial','continuous'}))
            for i=1:numCells
                theTrials=find(strcmp(uniqueCells{i},EPoverview.workData.cellNames));
                theTrials=find(~badTrials(theTrials));
                
                switch EPoverview.QC.mode
                    case 'avgNum'
                        tempQC(i) = sum(QCdata(theTrials));
                    case 'subNum'
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'blinkTrial'
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'saccadeTrial'
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'saccadeOnset'
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'moveTrial'
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'badTrials'
                        theTrials=find(strcmp(uniqueCells{i},EPoverview.workData.cellNames));
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'badChans'
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'noise'
                        tempQC(i) = mean(QCdata(theTrials));
                    case 'std'
                        tempQC(i) = mean(QCdata(theTrials));
                end;
                
            end;
            QCdata=tempQC;
        end;
        
        for i=1:numSubs
            tableData{i,1}=char(EPoverview.workData.subNames(i));
            tableData{i,2}=char(EPoverview.workData.subTypes(i));
            tableData(i,3:2+numCells)=num2cell(QCdata(i,:));
        end;
        
        tableNames{1}='names';
        tableNames{2}='type';
        for i=1:numCells
            tableNames(2+i)=uniqueCells(i);
        end;
        
        columnEditable =  repmat(false,1,2+numCells);
        
        ColumnFormat=cell(1,numCells+2);
        for i=1:numCells
            ColumnFormat{i+2}='bank';
        end;
        
        EPoverview.handles.QC.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'RearrangeableColumns','on',...
            'Position',[100 100 500 min(500,windowHeight-200)]);
        
        EPoverview.handles.QC.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''exportQC'');');
        
        callBackEnd=['ep_editData(''start'');'];
        
        if ~isempty(EPoverview.workData.avgNum)
            EPoverview.handles.QC.avgNum = uicontrol('Style', 'pushbutton', 'String', 'Trials','FontSize',EPmain.fontsize,...
                'TooltipString','Number of trials going into each waveform.',...
                'Position', [100 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''avgNum'';',callBackEnd]);
        end;
        
        if ~isempty(EPoverview.workData.subNum)
            EPoverview.handles.QC.subNum = uicontrol('Style', 'pushbutton', 'String', 'Subs','FontSize',EPmain.fontsize,...
                'TooltipString','Number of subjects going into each waveform.',...
                'Position', [160 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''subNum'';',callBackEnd]);
        end;
        
        EPoverview.handles.QC.blinkTrial = uicontrol('Style', 'pushbutton', 'String', 'Blinks','FontSize',EPmain.fontsize,...
            'TooltipString','Number of blink-corrected trials going into each waveform.',...
            'Position', [220 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''blinkTrial'';',callBackEnd]);
        EPoverview.handles.QC.saccadeTrial = uicontrol('Style', 'pushbutton', 'String', 'Saccades','FontSize',EPmain.fontsize,...
            'TooltipString','Number of saccade-corrected trials going into each waveform.',...
            'Position', [280 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''saccadeTrial'';',callBackEnd]);
        EPoverview.handles.QC.saccadeOnset = uicontrol('Style', 'pushbutton', 'String', 'Sacc Ms','FontSize',EPmain.fontsize,...
            'TooltipString','Onset of saccades.',...
            'Position', [340 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''saccadeOnset'';',callBackEnd]);
        EPoverview.handles.QC.moveTrial = uicontrol('Style', 'pushbutton', 'String', 'Move','FontSize',EPmain.fontsize,...
            'TooltipString','Number of movement-corrected trials going into each waveform.',...
            'Position', [400 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''moveTrial'';',callBackEnd]);
        EPoverview.handles.QC.badTrials = uicontrol('Style', 'pushbutton', 'String', 'BadTrials','FontSize',EPmain.fontsize,...
            'TooltipString','Proportion of trials excluded from each averaged waveform due to being bad trials.',...
            'Position', [460 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''badTrials'';',callBackEnd]);
        %if any(any(any(EPoverview.workData.analysis.badChans(repmat(~badTrials,[1,1,numChans])) < 0)))
            EPoverview.handles.QC.badChans = uicontrol('Style', 'pushbutton', 'String', 'BadChan','FontSize',EPmain.fontsize,...
                'TooltipString','Proportion of bad channels in each set of waveforms, not include bad trials.',...
                'Position', [520 50 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''badChans'';',callBackEnd]);
        %end;

        if ~isempty(EPoverview.workData.noise)
            EPoverview.handles.QC.noise = uicontrol('Style', 'pushbutton', 'String', 'Noise','FontSize',EPmain.fontsize,...
                'TooltipString','Noise level in each averaged waveform, as computed by the +/- reference.',...
                'Position', [100 20 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''noise'';',callBackEnd]);
        end;
        if ~isempty(EPoverview.workData.std)
            EPoverview.handles.QC.std = uicontrol('Style', 'pushbutton', 'String', 'Std','FontSize',EPmain.fontsize,...
                'TooltipString','Mean standard deviation in each averaged waveform.',...
                'Position', [160 20 60 30], 'Callback', ['global EPoverview;','EPoverview.QC.mode=''std'';',callBackEnd]);
        end;
        
        switch EPoverview.QC.mode
            case 'avgNum'
                set(EPoverview.handles.QC.avgNum,'ForegroundColor','blue');
            case 'subNum'
                set(EPoverview.handles.QC.subNum,'ForegroundColor','blue');
            case 'blinkTrial'
                set(EPoverview.handles.QC.blinkTrial,'ForegroundColor','blue');
            case 'saccadeTrial'
                set(EPoverview.handles.QC.saccadeTrial,'ForegroundColor','blue');
            case 'saccadeOnset'
                set(EPoverview.handles.QC.saccadeOnset,'ForegroundColor','blue');
            case 'moveTrial'
                set(EPoverview.handles.QC.moveTrial,'ForegroundColor','blue');
            case 'badTrials'
                set(EPoverview.handles.QC.badTrials,'ForegroundColor','blue');
            case 'badChans'
                set(EPoverview.handles.QC.badChans,'ForegroundColor','blue');
            case 'noise'
                set(EPoverview.handles.QC.noise,'ForegroundColor','blue');
            case 'std'
                set(EPoverview.handles.QC.std,'ForegroundColor','blue');
        end;
        
    case 'startPCA'
        
        if strcmp(EPoverview.PCA.mode,'Summary')
            if isfield(EPoverview.workData.pca,'PCAmode2')
                title='First PCA';
            else
                title='The PCA';
            end;
            uicontrol('Style','text','HorizontalAlignment','left','String', title,'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-120 150 20]);
            switch EPoverview.workData.pca.PCAmode
                case 'spat'
                    theMode='spatial';
                case 'temp'
                    theMode='temporal';
                case 'freq'
                    theMode='frequency';
                otherwise
                    beep()
                    disp('Oops')
                    theMode='Oops';
            end;
            
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Mode','FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-140 150 20]);
            uicontrol('Style','text','HorizontalAlignment','left','String', theMode,'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-160 150 20]);
            
            switch EPoverview.workData.pca.ROTATION
                case 'QMIN'
                    theRotationName = 'Quartimin';
                case 'OMIN'
                    theRotationName = 'Oblimin';
                case 'CRFE'
                    theRotationName = 'Crawford-Ferguson';
                case 'MINE'
                    theRotationName = 'minimum entropy';
                case 'IPSC'
                    theRotationName = 'Bentler''s invariant pattern simplicity criterion';
                case 'TIIC'
                    theRotationName = 'Comrey''s tandem II criterion';
                case 'GMIN'
                    theRotationName = 'geomin';
                case 'IMAX'
                    theRotationName = 'Infomax';
                case 'MMER'
                    theRotationName = 'McCammon minimum entropy ratio';
                case 'UNRT'
                    theRotationName = 'none';
                case 'VMAX'
                    theRotationName = 'Varimax';
                case 'QMAX'
                    theRotationName = 'Quartimax';
                case 'PMAX'
                    theRotationName = 'Promax';
                case 'VOMN'
                    theRotationName = 'Variable Oblimin';
                otherwise
                    beep()
                    disp('Oops')
                    theRotationName = 'oops';
            end;
            
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Rotation','FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-180 150 20]);
            uicontrol('Style','text','HorizontalAlignment','left','String', theRotationName,'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-200 150 20]);
            
            switch EPoverview.workData.pca.MAT_TYPE
                case 'SCP'
                    theMat='Cross-Product';
                case 'COV'
                    theMat='Covariance';
                case 'COR'
                    theMat='Correlation';
                otherwise
                    beep()
                    disp('Oops')
                    theMat='Oops';
            end;
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Matrix Type','FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-220 150 20]);
            uicontrol('Style','text','HorizontalAlignment','left','String', theMat,'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-240 150 20]);
            
            switch EPoverview.workData.pca.LOADING
                case 'K'
                    theLoadings='Kaiser';
                case 'C'
                    theLoadings='Covariance';
                case 'N'
                    theLoadings='None';
                case 'W'
                    theLoadings='Cureton-Mulaik';
                otherwise
                    beep()
                    disp('Oops')
                    theLoadings='Oops';
            end;
            if strcmp(theRotationName,'Infomax')
                theLoadings='None';
            end;
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Loading Weighting','FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-260 150 20]);
            uicontrol('Style','text','HorizontalAlignment','left','String', theLoadings,'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-280 150 20]);
            
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Rotation Option','FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-300 150 20]);
            uicontrol('Style','text','HorizontalAlignment','left','String', num2str(EPoverview.workData.pca.RotOpt),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-320 150 20]);
            
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Factors Retained','FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-340 150 20]);
            uicontrol('Style','text','HorizontalAlignment','left','String', num2str(EPoverview.workData.pca.numFacs),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-360 150 20]);
            
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Total Variance','FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-380 150 20]);
            uicontrol('Style','text','HorizontalAlignment','left','String', num2str(EPoverview.workData.pca.facVarTot),'FontSize',EPmain.fontsize,...
                'Position',[100 windowHeight-400 150 20]);
            
            if isfield(EPoverview.workData.pca,'PCAmode2')
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Second Rotation','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-120 150 20]);
                switch EPoverview.workData.pca.PCAmode2
                    case 'spat'
                        theMode='spatial';
                    case 'temp'
                        theMode='temporal';
                    otherwise
                        beep()
                        disp('Oops')
                        theMode='Oops';
                end;
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Mode','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-140 150 20]);
                uicontrol('Style','text','HorizontalAlignment','left','String', theMode,'FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-160 150 20]);
                
                switch EPoverview.workData.pca.ROTATION2
                    case 'QMIN'
                        theRotationName = 'Quartimin';
                    case 'OMIN'
                        theRotationName = 'Oblimin';
                    case 'CRFE'
                        theRotationName = 'Crawford-Ferguson';
                    case 'MINE'
                        theRotationName = 'minimum entropy';
                    case 'IPSC'
                        theRotationName = 'Bentler''s invariant pattern simplicity criterion';
                    case 'TIIC'
                        theRotationName = 'Comrey''s tandem II criterion';
                    case 'GMIN'
                        theRotationName = 'geomin';
                    case 'IMAX'
                        theRotationName = 'Infomax';
                    case 'MMER'
                        theRotationName = 'McCammon minimum entropy ratio';
                    case 'UNRT'
                        theRotationName = 'none';
                    case 'VMAX'
                        theRotationName = 'Varimax';
                    case 'QMAX'
                        theRotationName = 'Quartimax';
                    case 'PMAX'
                        theRotationName = 'Promax';
                    case 'VOMN'
                        theRotationName = 'Variable Oblimin';
                    otherwise
                        beep()
                        disp('Oops')
                        theRotationName = 'oops';
                end;
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Rotation','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-180 150 20]);
                uicontrol('Style','text','HorizontalAlignment','left','String', theRotationName,'FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-200 150 20]);
                
                switch EPoverview.workData.pca.MAT_TYPE2
                    case 'SCP'
                        theMat='Cross-Product';
                    case 'COV'
                        theMat='Covariance';
                    case 'COR'
                        theMat='Correlation';
                    otherwise
                        beep()
                        disp('Oops')
                        theMat='Oops';
                end;
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Matrix Type','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-220 150 20]);
                uicontrol('Style','text','HorizontalAlignment','left','String', theMat,'FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-240 150 20]);
                
                switch EPoverview.workData.pca.LOADING2
                    case 'K'
                        theLoadings='Kaiser';
                    case 'C'
                        theLoadings='Covariance';
                    case 'N'
                        theLoadings='None';
                    case 'W'
                        theLoadings='Cureton-Mulaik';
                    otherwise
                        beep()
                        disp('Oops')
                        theLoadings='Oops';
                end;
                if strcmp(theRotationName,'Infomax')
                    theLoadings='None';
                end;
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Loading Weighting','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-260 150 20]);
                uicontrol('Style','text','HorizontalAlignment','left','String', theLoadings,'FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-280 150 20]);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Rotation Option','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-300 150 20]);
                uicontrol('Style','text','HorizontalAlignment','left','String', num2str(EPoverview.workData.pca.RotOpt2),'FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-320 150 20]);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Factors Retained','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-340 150 20]);
                uicontrol('Style','text','HorizontalAlignment','left','String', num2str(EPoverview.workData.pca.numFacs2),'FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-360 150 20]);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Total Variance','FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-380 150 20]);
                uicontrol('Style','text','HorizontalAlignment','left','String', num2str(sum(EPoverview.workData.pca.facVar.*EPoverview.workData.pca.facVarTotST)),'FontSize',EPmain.fontsize,...
                    'Position',[300 windowHeight-400 150 20]);
                
                PCAdata=[EPoverview.workData.pca.facVarTotST; EPoverview.workData.pca.facVar.*EPoverview.workData.pca.facVarTotST]';
                tableNames={'Factor Var','Total Var'};
                
                tableData=num2cell(PCAdata);
                
                columnEditable =  repmat(false,1,size(PCAdata,2));
                
                ColumnFormat = repmat([],1,size(PCAdata,2));
                
                EPoverview.handles.PCA.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                    'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                    'RearrangeableColumns','on',...
                    'Position',[500 windowHeight-400 200 200]);
                
            end;
            
        else
            
            %prepare contents of the table based on the PCA setting
            
            switch EPoverview.PCA.mode
                case 'FacPat'
                    PCAdata=EPoverview.workData.pca.FacPat;
                    tableNames=EPoverview.workData.pca.facNames;
                case 'FacStr'
                    PCAdata=EPoverview.workData.pca.FacStr;
                    tableNames=EPoverview.workData.pca.facNames;
                case 'FacScr'
                    PCAdata=EPoverview.workData.pca.FacScr;
                    tableNames=EPoverview.workData.pca.facNames;
                case 'FacCof'
                    PCAdata=EPoverview.workData.pca.FacCof;
                    tableNames=EPoverview.workData.pca.facNames;
                case 'FacCor'
                    PCAdata=EPoverview.workData.pca.FacCor;
                    tableNames=EPoverview.workData.pca.facNames;
                case 'FacPatST'
                    PCAdata=EPoverview.workData.pca.FacPatST;
                    tableNames=EPoverview.workData.pca.facNames2;
                case 'FacStrST'
                    PCAdata=EPoverview.workData.pca.FacStrST;
                    tableNames=EPoverview.workData.pca.facNames2;
                case 'FacScrST'
                    PCAdata=EPoverview.workData.pca.FacScrST;
                    tableNames=EPoverview.workData.pca.facNames2;
                case 'FacCofST'
                    PCAdata=EPoverview.workData.pca.FacCofST;
                    tableNames=EPoverview.workData.pca.facNames2;
                case 'FacCorST'
                    PCAdata=EPoverview.workData.pca.FacCorST;
                    tableNames=EPoverview.workData.pca.facNames2;
            end;
            
            
            tableData=num2cell(PCAdata);
            
            columnEditable =  repmat(false,1,size(PCAdata,2));
            
            ColumnFormat = repmat([],1,size(PCAdata,2));
            
            EPoverview.handles.PCA.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'RearrangeableColumns','on',...
                'Position',[100 100 500 min(500,windowHeight-200)]);
            
            EPoverview.handles.PCA.hExport = uicontrol('Style', 'pushbutton', 'String', 'Export','FontSize',EPmain.fontsize,...
                'Position', [20 windowHeight-460 60 30], 'Callback', 'ep_editData(''exportPCA'');');
            
        end;
        
        
        callBackEnd=['ep_editData(''start'');'];
        
        
        if isfield(EPoverview.workData.pca,'PCAmode')
            EPoverview.handles.PCA.Summary = uicontrol('Style', 'pushbutton', 'String', 'Summary','FontSize',EPmain.fontsize,...
                'Position', [100 50 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''Summary'';',callBackEnd]);
            if isfield(EPoverview.workData.pca,'FacPat')
                EPoverview.handles.PCA.FacPat = uicontrol('Style', 'pushbutton', 'String', 'FacPat','FontSize',EPmain.fontsize,...
                    'Position', [165 50 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacPat'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacStr')
                EPoverview.handles.PCA.FacStr = uicontrol('Style', 'pushbutton', 'String', 'FacStr','FontSize',EPmain.fontsize,...
                    'Position', [230 50 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacStr'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacScr')
                EPoverview.handles.PCA.FacScr = uicontrol('Style', 'pushbutton', 'String', 'FacScr','FontSize',EPmain.fontsize,...
                    'Position', [295 50 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacScr'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacCof')
                EPoverview.handles.PCA.FacCof = uicontrol('Style', 'pushbutton', 'String', 'FacCof','FontSize',EPmain.fontsize,...
                    'Position', [360 50 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacCof'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacCor')
                EPoverview.handles.PCA.FacCor = uicontrol('Style', 'pushbutton', 'String', 'FacCor','FontSize',EPmain.fontsize,...
                    'Position', [425 50 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacCor'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacPatST')
                EPoverview.handles.PCA.FacPatST = uicontrol('Style', 'pushbutton', 'String', 'FacPatST','FontSize',EPmain.fontsize,...
                    'Position', [165 20 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacPatST'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacStrST')
                EPoverview.handles.PCA.FacStrST = uicontrol('Style', 'pushbutton', 'String', 'FacStrST','FontSize',EPmain.fontsize,...
                    'Position', [230 20 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacStrST'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacScrST')
                EPoverview.handles.PCA.FacScrST = uicontrol('Style', 'pushbutton', 'String', 'FacScrST','FontSize',EPmain.fontsize,...
                    'Position', [295 20 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacScrST'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacCofST')
                EPoverview.handles.PCA.FacCofST = uicontrol('Style', 'pushbutton', 'String', 'FacCofST','FontSize',EPmain.fontsize,...
                    'Position', [360 20 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacCofST'';',callBackEnd]);
            end;
            if isfield(EPoverview.workData.pca,'FacCorST')
                EPoverview.handles.PCA.FacCorST = uicontrol('Style', 'pushbutton', 'String', 'FacCorST','FontSize',EPmain.fontsize,...
                    'Position', [425 20 65 30], 'Callback', ['global EPoverview;','EPoverview.PCA.mode=''FacCorST'';',callBackEnd]);
            end;
        end;
        
        switch EPoverview.PCA.mode
            case 'Summary'
                set(EPoverview.handles.PCA.Summary,'ForegroundColor','blue');
            case 'FacPat'
                set(EPoverview.handles.PCA.FacPat,'ForegroundColor','blue');
            case 'FacStr'
                set(EPoverview.handles.PCA.FacStr,'ForegroundColor','blue');
            case 'FacScr'
                set(EPoverview.handles.PCA.FacScr,'ForegroundColor','blue');
            case 'FacCof'
                set(EPoverview.handles.PCA.FacCof,'ForegroundColor','blue');
            case 'FacCor'
                set(EPoverview.handles.PCA.FacCor,'ForegroundColor','blue');
            case 'FacPatST'
                set(EPoverview.handles.PCA.FacPatST,'ForegroundColor','blue');
            case 'FacStrST'
                set(EPoverview.handles.PCA.FacStrST,'ForegroundColor','blue');
            case 'FacScrST'
                set(EPoverview.handles.PCA.FacScrST,'ForegroundColor','blue');
            case 'FacCofST'
                set(EPoverview.handles.PCA.FacCofST,'ForegroundColor','blue');
            case 'FacCorST'
                set(EPoverview.handles.PCA.FacCorST,'ForegroundColor','blue');
        end;
        
%     case 'rename'
%         if length(varargin) < 1
%             msg{1}='No cell specified for renaming.';
%             [msg]=ep_errorMsg(msg);
%             return
%         end;
%         if ~isa(varargin{2},'char')
%             msg{1}='Old name needs to be specified.';
%             [msg]=ep_errorMsg(msg);
%             return
%         end;
%         
%         EPoverview.lastData=EPoverview.workData;
%         
%         oldName=varargin{2};
%         
%         [uniqueCells m order]=unique(EPoverview.workData.cellNames,'first');
%         
%         theCell=order(strcmp(oldName,uniqueCells)); %the cell number as ordered alphabetically
%         
%         newName=get(EPoverview.handles.cells.hCellName(theCell),'String');
%         
%         EPoverview.workData.cellNames(strcmp(oldName,EPoverview.workData.cellNames))={newName};
%         
%         [err]=ep_checkEPfile(EPoverview.workData);
%         
%         ep_editData('start');
        
    case 'renameDataName'
        
        EPoverview.lastData=EPoverview.workData;
        
        newName=get(EPoverview.handles.overview.hDataName,'String');
        
        EPoverview.workData.dataName=newName;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'renameExp'
        
        EPoverview.lastData=EPoverview.workData;
        
        newName=get(EPoverview.handles.overview.hExpName,'String');
        
        EPoverview.workData.ename=newName;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'montage'
        
        EPoverview.lastData=EPoverview.workData;
        
        newMontage=get(EPoverview.handles.overview.hMontage,'String');
        
        EPoverview.workData.montage=newMontage;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'Fs'
        
        EPoverview.lastData=EPoverview.workData;
        
        if strcmp(EPoverview.mode,'overview')
            newFs=str2num(get(EPoverview.handles.overview.hFs,'String'));
        elseif strcmp(EPoverview.mode,'samples')
            newFs=str2num(get(EPoverview.handles.samples.hFs,'String'));
        else
            msg{1}='Sampling rate change not supported from this pane.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.workData.timeNames=(EPoverview.workData.timeNames/(1000/EPoverview.workData.Fs))*(1000/newFs);
        EPoverview.workData.Fs=newFs;
        
        if strcmp(EPoverview.workData.dataType,'continuous')
            numEpochs=floor(size(EPoverview.workData.data,2)/ceil(EPoverview.workData.Fs)); %excess time points are tacked onto final epoch
            if numEpochs == 0
                numEpochs =1;
            end;
            EPoverview.workData.analysis.blinkTrial=zeros(numSubs,numEpochs);
            EPoverview.workData.analysis.saccadeTrial=zeros(numSubs,numEpochs);
            EPoverview.workData.analysis.saccadeOnset=zeros(numSubs,numEpochs);
            EPoverview.workData.analysis.moveTrial=zeros(numSubs,numEpochs);
            EPoverview.workData.analysis.badTrials=zeros(numSubs,numEpochs);
            EPoverview.workData.analysis.badChans=zeros(numSubs,numEpochs,numChan);
            disp('Resetting artifact correction fields to zero since length of one-second epochs has been changed.');
        end;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'baseline'
        
        EPoverview.lastData=EPoverview.workData;
        
        if strcmp(EPoverview.mode,'overview')
            newBaseline=(str2num(get(EPoverview.handles.overview.hBaseline,'String')))/(1000/EPoverview.workData.Fs);
        elseif strcmp(EPoverview.mode,'samples')
            newBaseline=(str2num(get(EPoverview.handles.samples.hBaseline,'String')))/(1000/EPoverview.workData.Fs);
        else
            msg{1}='Baseline change not supported from this pane.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.workData.timeNames=(((0:length(EPoverview.workData.timeNames)-1)+EPoverview.workData.baseline)*(1000/EPoverview.workData.Fs));
        
        EPoverview.workData.baseline=newBaseline;
        
        EPoverview.workData.timeNames=(((0:length(EPoverview.workData.timeNames)-1)-newBaseline)*(1000/EPoverview.workData.Fs));
        
        EPoverview.workData.timeNames=EPoverview.workData.timeNames(:); %make sure is a column vector
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');

    case 'refType'
        
        EPoverview.lastData=EPoverview.workData;
        
        newRefType=get(EPoverview.handles.overview.hRefType,'Value');
        
        EPoverview.workData.reference.type=refTypes{newRefType};
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');

    case 'origRefChans1'
        
        EPoverview.lastData=EPoverview.workData;
        
        newChan=str2num(get(EPoverview.handles.overview.hOrigRefChans1,'String'));
        
        EPoverview.workData.reference.original(1)=newChan;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');

    case 'origRefChans2'
        
        EPoverview.lastData=EPoverview.workData;
        
        newChan=str2num(get(EPoverview.handles.overview.hOrigRefChans2,'String'));
        
        EPoverview.workData.reference.original(2)=newChan;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');

    case 'currRefChans1'
        
        EPoverview.lastData=EPoverview.workData;
        
        newChan=str2num(get(EPoverview.handles.overview.hCurrRefChans1,'String'));
        
        EPoverview.workData.reference.current(1)=newChan;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');

    case 'currRefChans2'
        
        EPoverview.lastData=EPoverview.workData;
        
        newChan=str2num(get(EPoverview.handles.overview.hCurrRefChans2,'String'));
        
        EPoverview.workData.reference.current(2)=newChan;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');

    case 'ced'
        
        thisFile = mfilename('fullpath');
        [pathstr, name, ext] = fileparts(thisFile);
        theSeps=findstr(pathstr,filesep);
        pathstr=[pathstr(1:theSeps(end)-1) filesep 'electrodes'];
        [ced, pathstr] = uigetfile('*.ced',['Electrode Coordinate file (' num2str(length(EPoverview.workData.chanNames)) ' channels):'],pathstr);
        
        EPoverview.lastData=EPoverview.workData;
        
        if isnumeric(ced)
            button = questdlg('Delete electrode coordinates entirely?','','Yes','No','No');
            if strcmp(button,'Yes')
                EPoverview.workData.ced='none';
            end;
        else
            interpChans = '';
            if ~isempty(EPoverview.workData.eloc)
                interpChans = questdlg('Remap EEG data to new coordinates as well?','','Yes','No','No');
            end;
            whichCED=[pathstr ced];
            if ~isempty(EPoverview.workData.eloc) && strcmp(interpChans,'Yes')
                disp('Remapping EEG data to new coordinates as well.');
                sevenDdataIn=EPoverview.workData.data;
                oldElocs=EPoverview.workData.eloc;
                try
                    evalc('newElocs = readlocs([whichCED],''filetype'',''chanedit'');');
                catch
                    msg{1}=['The ced file ' ced ' did not work for some reason.  The error message was:'];
                    msg{2}=lasterr;
                    [msg]=ep_errorMsg(msg);
                    return
                end;
                nonEEGworkData=ep_selectData(EPoverview.workData,{find(~strcmp('EEG',EPoverview.workData.chanTypes)),[],[],[],[],[]});
                
                fileFormat=nonEEGworkData.fileFormat;
                dataType=nonEEGworkData.dataType;
                chanNames=nonEEGworkData.chanNames;
                timeNames=nonEEGworkData.timeNames;
                cellNames=nonEEGworkData.cellNames;
                subNames=nonEEGworkData.subNames;
                freqNames=nonEEGworkData.freqNames;
                chanTypes=nonEEGworkData.chanTypes;
                sevenDdata=nonEEGworkData.data;
                noise=nonEEGworkData.noise;
                stanDev=nonEEGworkData.std;
                stanDevCM=nonEEGworkData.stdCM;
                if isempty(nonEEGworkData.cov)
                    covMatrix=[];
                else
                    covMatrix=nonEEGworkData.cov.covMatrix;
                end;
                badChans=nonEEGworkData.analysis.badChans;
                reference=nonEEGworkData.reference;
                facVecS=nonEEGworkData.facVecS;
                impedances=nonEEGworkData.impedances;
                
                [eloc,chanNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM, badChans,reference,facVecS,covMatrix2,impedances]=ep_addEloc(whichCED,newElocs,fileFormat,dataType,chanNames,timeNames,cellNames,subNames,freqNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM, covMatrix,badChans,reference,facVecS,[],impedances);
                
                [sevenDdataOut]=ep_interpChans(sevenDdata, oldElocs, eloc);
                EPoverview.workData.data=sevenDdataOut;
                [err]=ep_checkEPfile(EPoverview.workData);
                ep_editData('start');
                
            else
                
                eloc=[];
                fileFormat=EPoverview.workData.fileFormat;
                dataType=EPoverview.workData.dataType;
                chanNames=EPoverview.workData.chanNames;
                timeNames=EPoverview.workData.timeNames;
                cellNames=EPoverview.workData.cellNames;
                subNames=EPoverview.workData.subNames;
                freqNames=EPoverview.workData.freqNames;
                chanTypes=EPoverview.workData.chanTypes;
                sevenDdata=EPoverview.workData.data;
                noise=EPoverview.workData.noise;
                stanDev=EPoverview.workData.std;
                stanDevCM=EPoverview.workData.stdCM;
                if isempty(EPoverview.workData.cov)
                    covMatrix=[];
                else
                    covMatrix=EPoverview.workData.cov.covMatrix;
                end;
                badChans=EPoverview.workData.analysis.badChans;
                reference=EPoverview.workData.reference;
                facVecS=EPoverview.workData.facVecS;
                impedances=EPoverview.workData.impedances;
                
                [eloc,chanNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM, badChans,reference,facVecS,covMatrix2,impedances]=ep_addEloc(whichCED,eloc,fileFormat,dataType,chanNames,timeNames,cellNames,subNames,freqNames,chanTypes,sevenDdata,noise,stanDev,stanDevCM, covMatrix,badChans,reference,facVecS,[],impedances);
            end;
            EPoverview.workData.eloc=eloc;
            EPoverview.workData.chanNames=chanNames;
            EPoverview.workData.timeNames=timeNames;
            EPoverview.workData.cellNames=cellNames;
            EPoverview.workData.freqNames=freqNames;
            EPoverview.workData.chanTypes=chanTypes;
            EPoverview.workData.data=sevenDdata;
            EPoverview.workData.noise=noise;
            EPoverview.workData.std=stanDev;
            EPoverview.workData.stdCM=stanDevCM;
            EPoverview.workData.analysis.badChans=badChans;
            EPoverview.workData.reference=reference;
            EPoverview.workData.facVecS=facVecS;
            if ~isempty(EPoverview.workData.cov)
                EPoverview.workData.cov.covMatrix=covMatrix2;
            end;
            EPoverview.workData.impedances=impedances;
            EPoverview.channels.select=repmat(false,length(EPoverview.workData.chanNames)+length(EPoverview.workData.implicit),1);
            EPoverview.channels.weights=repmat(0,length(EPoverview.workData.chanNames)+length(EPoverview.workData.implicit),1);
            EPoverview.workData.ced=ced;
            [err]=ep_checkEPfile(EPoverview.workData);
        end;
        
        ep_editData('start');
        
        
    case 'editSubs'
        
        EPoverview.lastData=EPoverview.workData;
        
        tableData=get(EPoverview.handles.subjects.hTable,'Data');
        
        if any(any(cellfun(@isempty,tableData(:,4))))
            disp('Subject names cannot be empty.');
            ep_editData('start');
            return
        end;
        
        %subjects reordered?
        if any(cell2mat(tableData(:,1))-[1:size(tableData,1)]')
            [tableData,sortOrder] = sortrows(tableData,1);
            
            [tempData]=ep_reorderData(EPoverview.workData,'subjects',sortOrder);
            if ~isempty(tempData)
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=tempData;
            end;         
            ep_editData('start');
        else
            EPoverview.subjects.select=cell2mat(tableData(:,2));
            EPoverview.subjects.weights=cell2mat(tableData(:,3));
            EPoverview.workData.subNames =tableData(:,4);
            EPoverview.workData.subTypes =tableData(:,5);
            EPoverview.workData.subjectSpecs =tableData(:,6:5+numSpecs);
            set(EPoverview.handles.subjects.hWeightsSum,'String',sprintf('%f',sum(EPoverview.subjects.weights)));
            [err]=ep_checkEPfile(EPoverview.workData);
        end;        
        
    case 'deleteSubs'
        if sum(EPoverview.subjects.select) ==0
            msg{1}='No subjects specified for deletion.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        if sum(EPoverview.subjects.select) ==length(EPoverview.workData.subNames)
            msg{1}='There would be nothing left.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        keepSubjects=not(EPoverview.subjects.select);
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[],[],keepSubjects,[],[]});
        
        EPoverview.subjects.select=repmat(false,length(EPoverview.workData.subNames),1);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'addSubs'
        
        tableData=get(EPoverview.handles.subjects.hTable,'Data');
        weights=cell2mat(tableData(:,3));
        
        if ~any(weights)
            msg{1}='No weights specified for forming new subject average.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        subList=find(weights);
        subWeights=weights(subList);
        
        %make sure new subject doesn't duplicate existing subject names
        newName=[];
        posWeights=find(weights>0);
        negWeights=find(weights<0);
        for iWeight=1:length(posWeights)
            newName=[newName '+' EPoverview.workData.subNames{posWeights(iWeight)}];
        end;
        for iWeight=1:length(negWeights)
            newName=[newName '-' EPoverview.workData.subNames{negWeights(iWeight)}];
        end;
        
        comb=1;
        candidateName=newName;
        while any(strcmp(candidateName, EPoverview.workData.subNames))
            comb=comb+1;
            candidateName=[newName num2str(comb)];
        end;
        if comb > 1
            newName=[newName '_' num2str(comb)];
        end;
 
        tempData=ep_combineData(EPoverview.workData,'subjects',subList,subWeights,newName);
        if ~isempty(tempData)
            EPoverview.lastData=EPoverview.workData;
            EPoverview.workData=tempData;
        end;
        
        EPoverview.subjects.weights=repmat(0,length(EPoverview.workData.subNames),1);
        EPoverview.subjects.select=repmat(false,length(EPoverview.workData.subNames),1);
        
        ep_editData('start');
        
    case 'clearSubWeights'
        
        EPoverview.subjects.weights=repmat(0,length(EPoverview.workData.subNames),1);
        EPoverview.subjects.select=repmat(false,length(EPoverview.workData.subNames),1);
        
        ep_editData('start');
        
    case 'allSubs'
        
        EPoverview.subjects.weights=repmat(1,length(EPoverview.workData.subNames),1);
        EPoverview.subjects.select=repmat(true,length(EPoverview.workData.subNames),1);
        
        ep_editData('start');
        
    case 'appendSubs'
        
        disp('Making assumption that file format is same as that of the original file.');
        importFormat=EPoverview.workData.fileFormat;
        dataType=EPoverview.workData.dataType;
        if any(strcmp(dataType,{'continuous','single_trial'}))
            msg{1}='It is not possible to append additional subjects to this type of file.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        textPrefs.firstRow=EPmain.preferences.general.firstRow;
        textPrefs.lastRow=EPmain.preferences.general.lastRow;
        textPrefs.firstCol=EPmain.preferences.general.firstCol;
        textPrefs.lastCol=EPmain.preferences.general.lastCol;
        textPrefs.orientation=EPmain.preferences.general.orientation;
        
        [fileNames, activeDirectory]=ep_getFilesUI(importFormat);
        if isempty(fileNames) || ((length(fileNames{1})==1) && (fileNames{1}==0))
            msg{1}='No filenames selected. You have to click on a name.';
            [msg]=ep_errorMsg(msg);
            return
        end
        
        EPoverview.lastData=EPoverview.workData;
        
        ced=EPoverview.workData.ced;
        montage=EPoverview.workData.montage;
        newNames=EPoverview.workData.subNames;
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
            inArg{end+1}='ced';
            inArg{end+1}=ced;
            inArg{end+1}='montage';
            inArg{end+1}=montage;
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
            
            newData=ep_readData(inArg);
            errorFlag=0;
            if isempty(newData)
                errorFlag=1;
                disp(['Error: Was unable to read ' inArg{2}]);
            else                
                if ~strcmp(dataType,newData.dataType)
                    errorFlag=1;
                    disp(['Error: Mismatched data type for ' inArg{2}]);
                end;
                
                if length(EPoverview.workData.cellNames) ~= length(newData.cellNames)
                    errorFlag=1;
                    disp(['Error: Mismatched number of cells for ' inArg{2}]);
                end;
                
                if length(EPoverview.workData.chanNames) ~= length(newData.chanNames)
                    errorFlag=1;
                    disp(['Error: Mismatched number of channels for ' inArg{2}]);
                end;
                
                if length(EPoverview.workData.timeNames) ~= length(newData.timeNames)
                    errorFlag=1;
                    disp(['Error: Mismatched number of time points for ' inArg{2}]);
                end;
                
                if length(EPoverview.workData.facNames) ~= length(newData.facNames)
                    errorFlag=1;
                    disp(['Error: Mismatched number of factors for ' inArg{2}]);
                end;
                
                if size(EPoverview.workData.subjectSpecs,2) ~= size(newData.subjectSpecs,2)
                    errorFlag=1;
                    disp(['Error: Mismatched number of subject specs for ' inArg{2}]);
                end;
                
                for i=1:length(newData.subNames)
                    if any(strcmp(newData.subNames{i},newNames))
                        subNameSuffix='sub001';
                        sameName=1;
                        suffix=0;
                        while sameName
                            sameName=0;
                            for j=1:length(newNames)
                                if strcmp(newNames{j},subNameSuffix)
                                    sameName=1;
                                end;
                            end;
                            if sameName
                                suffix=suffix+1;
                                subNameSuffix=['sub' sprintf('%03d',suffix)];
                            end;
                        end;
                        newData.subNames{i,1}=subNameSuffix;
                    end;
                    newNames{end+1,1}=newData.subNames{i};
                end;
                
                if errorFlag
                    beep();
                else
                    workData=ep_addData(EPoverview.workData,newData,'subjects');
                end;
            end;
        end;
        
        if isempty(workData)
            disp('Error: no data.');
        else
            [err]=ep_checkEPfile(EPoverview.workData);
            if ~err
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=workData;
                EPoverview.subjects.select=repmat(false,length(EPoverview.workData.subNames),1);
                EPoverview.subjects.weights=repmat(0,length(EPoverview.workData.subNames),1);
            end;
        end;
        
        ep_editData('start');
        
    case 'exportSubs'
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Subject Specs',[EPoverview.workData.dataName EPmain.preferences.general.subjectSpecSuffix]);
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        ep_writeSubjectText(EPoverview.workData, [PathName filesep FileName]);
        
        ep_editData('start');
        
    case 'importSubs'
        
        [FileName,PathName,FilterIndex] = uigetfile('*.txt','Import Subject Specs',[EPoverview.workData.dataName EPmain.preferences.general.subjectSpecSuffix]);
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        EPoverview.workData = ep_readSubjectSpecText(EPoverview.workData, [PathName filesep FileName]);
        
        ep_editData('start');
        
    case 'addSubSpec'
        
        specName = inputdlg('Name of new spec?');
        pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
        if ~isempty(specName)
            EPoverview.lastData=EPoverview.workData;
            EPoverview.workData.subjectSpecNames(end+1)=specName;
            
            for i=1:numSubs
                EPoverview.workData.subjectSpecs{i,length(EPoverview.workData.subjectSpecNames)}='';
            end;
        end;
        
        ep_editData('start');
        
    case 'minusSubSpec'
        
        [Selection,ok] = listdlg('PromptString','Choose specs to delete','ListString',EPoverview.workData.subjectSpecNames);
        if ok
            EPoverview.lastData=EPoverview.workData;
            deleteList=setdiff([1:length(EPoverview.workData.subjectSpecNames)],Selection);
            EPoverview.workData.subjectSpecNames=EPoverview.workData.subjectSpecNames(deleteList);
            EPoverview.workData.subjectSpecs=EPoverview.workData.subjectSpecs(:,deleteList);
        end;
        
        ep_editData('start');
        
    case 'deleteCells'
        if sum(EPoverview.cells.select) ==0
            msg{1}='No cells specified for deletion.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        if sum(EPoverview.cells.select) ==length(unique(EPoverview.workData.cellNames))
            msg{1}='There would be nothing left.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        tableData=get(EPoverview.handles.cells.hTable,'Data');
        keepCells=ismember(EPoverview.workData.cellNames,tableData(not(EPoverview.cells.select),4));        
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[],keepCells,[],[],[]});
        
        EPoverview.cells.weights=repmat(0,length(EPoverview.workData.cellNames),1);
        EPoverview.cells.select=repmat(false,length(unique(EPoverview.workData.cellNames)),1);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'addCells'
        
        tableData=get(EPoverview.handles.cells.hTable,'Data');
        weights=cell2mat(tableData(:,3));
        
        if ~any(weights)
            msg{1}='No weights specified for forming new cell average.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        cellList=find(weights);
        cellWeights=weights(cellList);
        
        %make sure new cell doesn't duplicate existing cell names
        newName=[];
        posWeights=find(weights>0);
        negWeights=find(weights<0);
        for iWeight=1:length(posWeights)
            newName=[newName '+' EPoverview.workData.cellNames{posWeights(iWeight)}];
        end;
        for iWeight=1:length(negWeights)
            newName=[newName '-' EPoverview.workData.cellNames{negWeights(iWeight)}];
        end;
        
        comb=1;
        candidateName=newName;
        while any(strcmp(candidateName, EPoverview.workData.cellNames))
            comb=comb+1;
            candidateName=[newName num2str(comb)];
        end;
        if comb > 1
            newName=[newName '_' num2str(comb)];
        end;
 
        tempData=ep_combineData(EPoverview.workData,'cells',cellList,cellWeights,newName,EPoverview.cells.TrlWgts);
        if ~isempty(tempData)
            EPoverview.lastData=EPoverview.workData;
            EPoverview.workData=tempData;
        end;
        
        EPoverview.cells.weights=repmat(0,length(EPoverview.workData.cellNames),1);
        EPoverview.cells.select=repmat(false,length(EPoverview.workData.cellNames),1);
        
        ep_editData('start');
        
    case 'addTrials'
        
        tableData=get(EPoverview.handles.cells.hTable,'Data');
        selectedCell=find(cell2mat(tableData(:,2)));
        tempData=EPoverview.workData;
        
        for iCell=1:length(selectedCell)
            theCell=selectedCell(iCell);
            cellList=find(strcmp(tableData{theCell,4},EPoverview.workData.cellNames));
            cellWeights=ones(length(cellList),1);
            
            if length(cellList) ==1
                msg{1}='Only one waveform for averaging.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            %make sure new cell doesn't duplicate existing cell names
            newName=['AVG-' tableData{theCell,4}];
            
            comb=1;
            candidateName=newName;
            while any(strcmp(candidateName, EPoverview.workData.cellNames))
                comb=comb+1;
                candidateName=[newName num2str(comb)];
            end;
            if comb > 1
                newName=[newName '_' num2str(comb)];
            end;
            
            tempData2=ep_combineData(tempData,'cells',cellList,cellWeights,newName,0);
            if ~isempty(tempData2)
                tempData=tempData2;
            end;
        end;
        
        if ~isempty(tempData)
            EPoverview.lastData=EPoverview.workData;
            EPoverview.workData=tempData;
        end;
        
        EPoverview.cells.weights=repmat(0,length(EPoverview.workData.cellNames),1);
        EPoverview.cells.select=repmat(false,length(EPoverview.workData.cellNames),1);
        
        ep_editData('start');
    
    case 'clearCellWeights'
        
        EPoverview.cells.weights=repmat(0,length(EPoverview.workData.cellNames),1);
        EPoverview.cells.select=repmat(false,length(EPoverview.workData.cellNames),1);
        
        ep_editData('start');
        
    case 'allCells'
        
        EPoverview.cells.weights=repmat(1,length(EPoverview.workData.cellNames),1);
        EPoverview.cells.select=repmat(true,length(EPoverview.workData.cellNames),1);
        
        ep_editData('start');
        
    case 'appendCells'
        
        disp('Making assumption that file format is same as that of the original file.');
        importFormat=EPoverview.workData.fileFormat;
        fileType=EPoverview.workData.dataType;
        if strcmp(fileType,'continuous')
            msg{1}='It is not possible to append additional cells to this type of file.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        textPrefs.firstRow=EPmain.preferences.general.firstRow;
        textPrefs.lastRow=EPmain.preferences.general.lastRow;
        textPrefs.firstCol=EPmain.preferences.general.firstCol;
        textPrefs.lastCol=EPmain.preferences.general.lastCol;
        textPrefs.orientation=EPmain.preferences.general.orientation;
        
        inArg{1}='format';
        inArg{2}=importFormat;
        inArg{3}='type';
        inArg{4}=fileType;
        inArg{5}='elecPrefs';
        inArg{6}=EPmain.preferences.general.rotateHead;
        inArg{7}='textPrefs';
        inArg{8}=textPrefs;
        inArg{9}='screenSize';
        inArg{10}=EPmain.scrsz;
        inArg{11}='FontSize';
        inArg{12}=EPmain.fontsize;
        
        [newData]=ep_readData(inArg);
        
        if isempty(newData)
            return
        end;
        
        tableData=get(EPoverview.handles.cells.hTable,'Data');
        
        if ~strcmp(fileType,newData.dataType)
            msg{1}='The file types do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        if length(EPoverview.workData.subNames) ~= length(newData.subNames)
            msg{1}='The numbers of subjects do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(EPoverview.workData.chanNames) ~= length(newData.chanNames)
            msg{1}='The numbers of channels do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(EPoverview.workData.timeNames) ~= length(newData.timeNames)
            msg{1}='The numbers of time points do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(EPoverview.workData.facNames) ~= length(newData.facNames)
            msg{1}='The numbers of factors do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
       
        if ~strcmp(EPoverview.workData.dataType,'single_trial')
            newCells=EPoverview.workData.cellNames;
            for i=1:length(newData.cellNames)
                if any(strcmp(newData.cellNames{i},newCells))
                    cellNameSuffix='cell001';
                    sameCell=1;
                    suffix=0;
                    while sameCell
                        sameCell=0;
                        for j=1:length(newCells)
                            if strcmp(newCells{j},cellNameSuffix)
                                sameCell=1;
                            end;
                        end;
                        if sameCell
                            suffix=suffix+1;
                            cellNameSuffix=['cell' sprintf('%03d',suffix)];
                        end;
                    end;
                    newData.cellNames{i,1}=cellNameSuffix;
                end;
                newCells{end+1,1}=newData.cellNames{i};
            end;
            
            combinedCells = [EPoverview.workData.cellNames; newData.cellNames];
            if (length(unique(combinedCells)) ~= length(combinedCells))
                msg{1}='Cell names cannot be duplicated except in single trial data.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        workData=ep_addData(EPoverview.workData,newData,'cells');
        if isempty(workData)
            disp('Error: no data.');
        else
            [err]=ep_checkEPfile(EPoverview.workData);
            if ~err
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=workData;
                EPoverview.cells.select=[EPoverview.cells.select; repmat(false,length(unique(newData.cellNames)),1)];
                EPoverview.cells.weights=[EPoverview.cells.weights; repmat(0,length(unique(newData.cellNames)),1)];
            end;
        end;
        
        ep_editData('start');
        
    case 'exportCells'
        
        tableData=get(EPoverview.handles.cells.hTable,'Data');
        tableData=tableData(:,4:end);
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Table');
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        [pathstr, name, ext] = fileparts(FileName);
        if ~strcmp(ext,'.txt')
            FileName=[FileName '.txt'];
        end;
        
        fid=fopen([PathName FileName],'w');
        
        for i=1:size(tableData,1)
            for j=1:size(tableData,2)
                theCell=tableData{i,j};
                if isnumeric(theCell)
                    fprintf(fid,'%d', theCell);
                else
                    fprintf(fid,'%s', theCell);
                end;
                if j < size(tableData,2)
                    fprintf(fid,'\t');
                end;
            end;
            if i < size(tableData,1)
                fprintf(fid,'\r');
            end;
        end;
        
        ep_editData('start');
        
    case 'trialSpecs'
        
        EPoverview.lastData=EPoverview.workData;
        
        tableData =get(EPoverview.handles.trials.hTable,'Data');
        
        if any(any(cellfun(@isempty,tableData(:,1:2))))
            disp('Cell and trial names cannot be empty.');
            ep_editData('start');
            return
        end;

        EPoverview.workData.cellNames=tableData(:,1);
        EPoverview.workData.trialNames=cell2mat(tableData(:,2));
        EPoverview.workData.recTime=cell2mat(tableData(:,3));
        EPoverview.workData.trialSpecs=tableData(:,4:end);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
    case 'exportTrials'
        
        tableData =get(EPoverview.handles.trials.hTable,'Data');
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Table');
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        [pathstr, name, ext] = fileparts(FileName);
        if ~strcmp(ext,'.txt')
            FileName=[FileName '.txt'];
        end;
        
        fid=fopen([PathName FileName],'w');
        
        for i=1:size(tableData,1)
            for j=1:size(tableData,2)
                theCell=tableData{i,j};
                if isnumeric(theCell)
                    fprintf(fid,'%d', theCell);
                else
                    fprintf(fid,'%s', theCell);
                end;
                if j < size(tableData,2)
                    fprintf(fid,'\t');
                end;
            end;
            if i < size(tableData,1)
                fprintf(fid,'\r');
            end;
        end;
        
        ep_editData('start');
        
    case 'editChannels'
        
        EPoverview.lastData=EPoverview.workData;
        
        numChans=length(EPoverview.workData.chanNames);
        
        tableData =get(EPoverview.handles.channels.hTable,'Data');
        
        if any(any(cellfun(@isempty,tableData(:,4))))
            disp('Channel names cannot be empty.');
            ep_editData('start');
            return
        end;
        
        %channels reordered?  Can't edit implicit sites, which must also stay at the end.
        if any(cell2mat(tableData(1:numChans,1))-[1:numChans]')
            [tableData(1:numChans,:),sortOrder] = sortrows(tableData(1:numChans,:),1);
            [tempData]=ep_reorderData(EPoverview.workData,'channels',sortOrder);
            if ~isempty(tempData)
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=tempData;
            end;
            EPoverview.channels.select(1:numChans)=EPoverview.channels.select(sortOrder);
            EPoverview.channels.weights(1:numChans)=EPoverview.channels.weights(sortOrder);
            ep_editData('start');
        else
            EPoverview.channels.select=cell2mat(tableData(:,2));
            EPoverview.channels.weights(1:numChans)=cell2mat(tableData(1:numChans,3));
            set(EPoverview.handles.channels.hWeights,'String',sprintf('%f',sum(EPoverview.channels.weights)));
            EPoverview.workData.chanNames=tableData(1:numChans,4);
            EPoverview.workData.chanTypes=tableData(1:numChans,4);
            for iChan=1:length(EPoverview.workData.chanTypes)
                if ~strcmp(EPoverview.workData.chanTypes{iChan},tableData(iChan,5))
                	EPoverview.workData.chanTypes(iChan)=tableData(iChan,5);
                    if any(strcmp(tableData(iChan,5),{'REG','ANS','ECG'}))
                        EPoverview.workData.eloc(iChan).theta=[];
                        EPoverview.workData.eloc(iChan).radius=[];
                        EPoverview.workData.eloc(iChan).X=[];
                        EPoverview.workData.eloc(iChan).Y=[];
                        EPoverview.workData.eloc(iChan).Z=[];
                        EPoverview.workData.eloc(iChan).sph_theta=[];
                        EPoverview.workData.eloc(iChan).sph_phi=[];
                        EPoverview.workData.eloc(iChan).sph_radius=[];
                    end;
                end;
            end;
            [err]=ep_checkEPfile(EPoverview.workData);
        end;
        
    case 'deleteChans'
        
        if sum(EPoverview.channels.select) ==0
            msg{1}='No channels specified for deletion.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        if sum(EPoverview.channels.select(1:numChans)) == numChans
            msg{1}='There would be no data left.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        keepChans=not(EPoverview.channels.select(1:numChans));
        keepImplicit=not(EPoverview.channels.select(numChans+1:end));
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{find(keepChans),[],[],[],[],[]});
        
        EPoverview.workData.implicit=EPoverview.workData.implicit(keepImplicit);
        
        EPoverview.channels.select=repmat(false,length(keepChans)+length(keepImplicit),1);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'addChans'  %adds a regional channel
        
        tableData=get(EPoverview.handles.channels.hTable,'Data');
        
        weights=cell2mat(tableData(:,3));
        
        if ~any(weights)
            msg{1}='No weights specified for forming new channel average.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if ~isempty(find(ismember({'FID','MGM','MGA','MGP','ANS','ECG'},tableData(find(weights),5))))
            msg{1}='At present only EEG data, including existing regional channels, can be combined into regional channels.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
%         if sum([any(find(ismember(tableData(find(weights),5),{'EEG','REG'}))) any(find(ismember(tableData(find(weights),5),{'ANS','ECG'}))) any(find(strcmp('MEG',tableData(find(weights),5))))]) > 1
%             msg{1}='Unlike types of data cannot be combined into a regional channel.';
%             [msg]=ep_errorMsg(msg);
%             return
%         end;
        
        chanList=find(weights);
        chanWeights=weights(chanList);
        
        %make sure new channel doesn't duplicate existing channel names
        comb=1;
        while any(strcmp(['combined' num2str(comb)], EPoverview.workData.chanNames))
            comb=comb+1;
        end;
 
        tempData=ep_combineData(EPoverview.workData,'channels',chanList,chanWeights,['combined' num2str(comb)]);
        if ~isempty(tempData)
            EPoverview.lastData=EPoverview.workData;
            EPoverview.workData=tempData;
        end;
        
        EPoverview.channels.weights=repmat(0,numChans+numImplicit+1,1);
        EPoverview.channels.select=repmat(false,numChans+numImplicit+1,1);
        
        ep_editData('start');
        
    case 'clearChanWeights'
        
        EPoverview.channels.weights=repmat(0,numChans+numImplicit,1);
        EPoverview.channels.select=repmat(false,numChans+numImplicit,1);
        
        ep_editData('start');
        
    case 'allChans'
        
        EPoverview.channels.weights=repmat(1,numChans+numImplicit,1);
        EPoverview.channels.select=repmat(true,numChans+numImplicit,1);
        
        ep_editData('start');
        
    case 'appendChans'
        
        disp('Making assumption that file format is same as that of the original file.');
        importFormat=EPoverview.workData.fileFormat;
        fileType=EPoverview.workData.dataType;
        
        textPrefs.firstRow=EPmain.preferences.general.firstRow;
        textPrefs.lastRow=EPmain.preferences.general.lastRow;
        textPrefs.firstCol=EPmain.preferences.general.firstCol;
        textPrefs.lastCol=EPmain.preferences.general.lastCol;
        textPrefs.orientation=EPmain.preferences.general.orientation;
        
        inArg{1}='format';
        inArg{2}=importFormat;
        inArg{3}='type';
        inArg{4}=fileType;
        inArg{5}='elecPrefs';
        inArg{6}=EPmain.preferences.general.rotateHead;
        inArg{7}='textPrefs';
        inArg{8}=textPrefs;
        inArg{9}='screenSize';
        inArg{10}=EPmain.scrsz;
        inArg{11}='FontSize';
        inArg{12}=EPmain.fontsize;
        
        [newData]=ep_readData(inArg);
        
        if isempty(newData)
            return
        end;
        
        tableData=get(EPoverview.handles.channels.hTable,'Data');
        
        if ~strcmp(fileType,newData.dataType)
            msg{1}='The file types do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if numCells ~= length(newData.cellNames)
            msg{1}='The numbers of cells do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if numSubs ~= length(newData.subNames)
            msg{1}='The numbers of subjects do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if numPoints ~= length(newData.timeNames)
            msg{1}='The numbers of time points do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(EPoverview.workData.facNames) ~= length(newData.facNames)
            msg{1}='The numbers of factors do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if size(EPoverview.workData.facVecS,2) ~= size(newData.facVecS,2)
            msg{1}='The two datasets are not equally factor compressed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        if ~isempty(EPoverview.workData.facVecS)
            EPoverview.workData.facVecS(end+1:length(newData.chanNames),:)=newData.facVecS;
        else
            EPoverview.workData.data(end+1:end+length(newData.chanNames),:,:,:,:,:,:)=newData.data;
        end;
        
        if ~isempty(EPoverview.workData.facData)
            if ~isempty(newData.facData)
                EPoverview.workData.facData(end+1:end+length(newData.chanNames),:,:,:,:,:,:)=newData.facData;
            else
                EPoverview.workData.facData(end+1:end+length(newData.chanNames),:,:,:,:,:,:)= zeros(size(newData.facData,1),size(newData.facData,2),length(newData.cellNames),length(newData.subNames),max(length(newData.facNames),1),max(length(newData.freqNames),1),max(length(newData.relNames),1));
            end;
        end;
        
        if ~isempty(EPoverview.workData.noise)
            if ~isempty(newData.noise)
                EPoverview.workData.noise(end+1:end+length(newData.chanNames),:,:,:,:)=newData.noise;
            else
                EPoverview.workData.noise(end+1:end+length(newData.chanNames),:,:,:,:)= zeros(length(newData.chanNames),length(newData.timeNames),length(newData.cellNames),length(newData.subNames),max(length(newData.facNames),1),max(length(newData.freqNames),1));
            end;
        end;
        if ~isempty(EPoverview.workData.std)
            if ~isempty(newData.std)
                EPoverview.workData.std(end+1:end+length(newData.chanNames),:,:,:,:,:)=newData.std;
            else
                EPoverview.workData.std(end+1:end+length(newData.chanNames),:,:,:,:,:)= zeros(length(newData.chanNames),length(newData.timeNames),length(newData.cellNames),length(newData.subNames),max(length(newData.facNames),1),max(length(newData.freqNames),1));
            end;
        end;
        if ~isempty(EPoverview.workData.stdCM)
            if ~isempty(newData.stdCM)
                EPoverview.workData.stdCM(end+1:end+length(newData.chanNames),:,:,:,:,:)=newData.stdCM;
            else
                EPoverview.workData.stdCM(end+1:end+length(newData.chanNames),:,:,:,:,:)= zeros(length(newData.chanNames),length(newData.timeNames),length(newData.cellNames),length(newData.subNames),max(length(newData.facNames),1),max(length(newData.freqNames),1));
            end;
        end;
        if ~isempty(EPoverview.workData.cov)
            if ~isempty(newData.cov)
                EPoverview.workData.cov.covMatrix(:,end+1:end+length(newData.chanNames),end+1:end+length(newData.chanNames))=newData.cov.covMatrix;
                EPoverview.workData.cov.covMatrix(:,1:length(EPoverview.workData.chanNames),length(EPoverview.workData.chanNames)+1:end)=NaN;
                EPoverview.workData.cov.covMatrix(:,length(EPoverview.workData.chanNames)+1:end,1:length(EPoverview.workData.chanNames))=NaN;
            else
                EPoverview.workData.cov.covMatrix(:,end+1:end+length(newData.chanNames),end+1:end+length(newData.chanNames))= NaN(length(newData.subNames),length(newData.chanNames),length(newData.chanNames));
            end;
        end;
        EPoverview.workData.analysis.badChans(:,:,end+1:end+length(newData.chanNames))=newData.analysis.badChans;
        
        EPoverview.workData.chanNames(end+1:end+length(newData.chanNames),1)=newData.chanNames;
        EPoverview.workData.chanTypes(end+1:end+length(newData.chanNames),1)=newData.chanTypes;
        EPoverview.workData.eloc(end+1:end+length(newData.chanNames))=newData.eloc;
        
        EPoverview.channels.select=repmat(false,length(EPoverview.workData.chanNames)+length(EPoverview.workData.implicit),1);
        EPoverview.channels.weights=repmat(0,length(EPoverview.workData.chanNames)+length(EPoverview.workData.implicit),1);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'exportChans'
        
        tableData=get(EPoverview.handles.chans.hTable,'Data');
        tableData=tableData(:,4:end);
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Table');
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        [pathstr, name, ext] = fileparts(FileName);
        if ~strcmp(ext,'.txt')
            FileName=[FileName '.txt'];
        end;
        
        fid=fopen([PathName FileName],'w');
        
        for i=1:size(tableData,1)
            for j=1:size(tableData,2)
                theCell=tableData{i,j};
                if isnumeric(theCell)
                    fprintf(fid,'%d', theCell);
                else
                    fprintf(fid,'%s', theCell);
                end;
                if j < size(tableData,2)
                    fprintf(fid,'\t');
                end;
            end;
            if i < size(tableData,1)
                fprintf(fid,'\r');
            end;
        end;
        
        ep_editData('start');
        
    case 'sortSubjects'
        
        sortSpec=get(EPoverview.handles.subjects.sort,'Value');
        specList=get(EPoverview.handles.subjects.sort,'String');
        
        if sortSpec == length(specList)-1
            [B, sortOrder]=sort(EPoverview.workData.subNames);
            
            [tempData]=ep_reorderData(EPoverview.workData,'subjects',sortOrder);
            if ~isempty(tempData)
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=tempData;
            end;
        elseif sortSpec < length(specList)-1
            theSpec=find(strcmp(specList{sortSpec},EPoverview.workData.subjectSpecNames));
            if isempty(theSpec)
                disp('Sort spec not present.')
            else
                sortArray=cellfun(@num2str,EPoverview.workData.subjectSpecs(:,theSpec),'UniformOutput',false);
                numArray=cellfun(@str2double,sortArray,'UniformOutput',false);
                numArray=cell2mat(numArray);
                if ~any(isnan(numArray))
                    sortArray=numArray;
                end;
                [B, sortOrder]=sort(sortArray);
                
                [tempData]=ep_reorderData(EPoverview.workData,'subjects',sortOrder);
                if ~isempty(tempData)
                    EPoverview.lastData=EPoverview.workData;
                    EPoverview.workData=tempData;
                end;
            end;
        end;
        
        ep_editData('start');
        
    case 'sortCells'
        
        sortSpec=get(EPoverview.handles.cells.sort,'Value');
        specList=get(EPoverview.handles.cells.sort,'String');
        
        if sortSpec == length(specList)-1
            [B, sortOrder]=sort(EPoverview.workData.cellNames);
            [tempData]=ep_reorderData(EPoverview.workData,'cells',sortOrder);
            if ~isempty(tempData)
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=tempData;
            end;
        elseif sortSpec < length(specList)-1
            theSpec=find(strcmp(specList{sortSpec},EPoverview.workData.trialSpecNames));
            if isempty(theSpec)
                disp('Sort spec not present.')
            else
                if sortSpec < length(specList)
                    sortArray=cellfun(@num2str,EPoverview.workData.trialSpecs(:,theSpec,1),'UniformOutput',false);
                    numArray=cellfun(@str2double,sortArray,'UniformOutput',false);
                    numArray=cell2mat(numArray);
                    if ~any(isnan(numArray))
                        sortArray=numArray;
                    end;
                    [B, sortOrder]=sort(sortArray);
                    
                    [tempData]=ep_reorderData(EPoverview.workData,'cells',sortOrder);
                    if ~isempty(tempData)
                        EPoverview.lastData=EPoverview.workData;
                        EPoverview.workData=tempData;
                    end;
                end;
            end;
        end;
        
        ep_editData('start');
        
    case 'sortTrials'
        
        sortSpec=get(EPoverview.handles.trials.sort,'Value');
        specList=get(EPoverview.handles.trials.sort,'String');
        
        if sortSpec == length(specList)-2
            [B, sortOrder]=sort(EPoverview.workData.cellNames);
            [tempData]=ep_reorderData(EPoverview.workData,'cells',sortOrder);
            if ~isempty(tempData)
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=tempData;
            end;
        elseif sortSpec == length(specList)-1
            [B, sortOrder]=sort(EPoverview.workData.trialNames);
            [tempData]=ep_reorderData(EPoverview.workData,'cells',sortOrder);
            if ~isempty(tempData)
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=tempData;
            end;
        elseif sortSpec < length(specList)-2
            theSpec=find(strcmp(specList{sortSpec},EPoverview.workData.trialSpecNames));
            if isempty(theSpec)
                disp('Sort spec not present.')
            else
                if sortSpec < length(specList)
                    sortArray=cellfun(@num2str,EPoverview.workData.trialSpecs(:,theSpec,1),'UniformOutput',false);
                    numArray=cellfun(@str2double,sortArray,'UniformOutput',false);
                    numArray=cell2mat(numArray);
                    if ~any(isnan(numArray))
                        sortArray=numArray;
                    end;
                    [B, sortOrder]=sort(sortArray);
                    [tempData]=ep_reorderData(EPoverview.workData,'cells',sortOrder);
                    if ~isempty(tempData)
                        EPoverview.lastData=EPoverview.workData;
                        EPoverview.workData=tempData;
                    end;
                end;
            end;
        end;
        
        ep_editData('start');
        
    case 'changeStartSamp'
        
        newStartSamp=str2double(get(EPoverview.handles.samples.startSamp,'string'));
        if isnan(newStartSamp)
            set(EPoverview.handles.samples.startSamp,'string',sprintf('%d', 1))
            msg{1}='The new start sample must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartSamp < 1
            set(EPoverview.handles.samples.startSamp,'string',sprintf('%d', 1))
            msg{1}='The new start sample cannot be earlier than the existing start sample time.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartSamp == 1
            msg{1}='The start sample was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartSamp > length(EPoverview.workData.timeNames)
            msg{1}='The new start sample cannot be larger than the number of total samples.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        trimSamples=newStartSamp-1;
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[trimSamples+1:length(EPoverview.workData.timeNames)],[],[],[],[]});
        
        ep_editData('start');
        
    case 'changeEndSamp'
        
        newEndSamp=str2double(get(EPoverview.handles.samples.endSamp,'string'));
        if isnan(newEndSamp)
            set(EPoverview.handles.samples.endSamp,'string',sprintf('%d', (length(EPoverview.workData.timeNames))))
            msg{1}='The new end sample must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndSamp > length(EPoverview.workData.timeNames)
            set(EPoverview.handles.samples.endSamp,'string',sprintf('%d', (length(EPoverview.workData.timeNames))))
            msg{1}='The new end sample cannot be later than the existing end sample time.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndSamp == length(EPoverview.workData.timeNames)
            msg{1}='The end sample was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndSamp < 1
            msg{1}='The new end sample cannot be smaller than one.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        

        trimSamples=length(EPoverview.workData.timeNames)-newEndSamp;
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[1:length(EPoverview.workData.timeNames)-trimSamples],[],[],[],[]});

        ep_editData('start');

    case 'changeStartMs'
        
        newStartMs=str2double(get(EPoverview.handles.samples.startMs,'string'));
        if isnan(newStartMs)
            set(EPoverview.handles.samples.startMs,'string',sprintf('%d', 1))
            msg{1}='The new start ms must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartMs < EPoverview.workData.timeNames(1)
            set(EPoverview.handles.samples.startMs,'string',sprintf('%d', 1))
            msg{1}='The new start ms cannot be earlier than the existing start ms.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartMs == EPoverview.workData.timeNames(1)
            msg{1}='The start ms was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartMs > EPoverview.workData.timeNames(end)
            msg{1}='The new start ms cannot be larger than the end ms.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        trimSamples=length(find(newStartMs > EPoverview.workData.timeNames));
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[trimSamples+1:length(EPoverview.workData.timeNames)],[],[],[],[]});
        
        ep_editData('start');
        
    case 'changeEndMs'
        
        newEndMs=str2double(get(EPoverview.handles.samples.endMs,'string'));
        if isnan(newEndMs)
            set(EPoverview.handles.samples.endMs,'string',sprintf('%d', EPoverview.workData.timeNames(end)+(1000/EPoverview.workData.Fs)));
            msg{1}='The new end ms must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndMs > EPoverview.workData.timeNames(end)+(1000/EPoverview.workData.Fs)
            set(EPoverview.handles.samples.endMs,'string',sprintf('%d', EPoverview.workData.timeNames(end)+(1000/EPoverview.workData.Fs)));
            msg{1}='The new end ms cannot be later than the existing end ms.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndMs == EPoverview.workData.timeNames(end)
            msg{1}='The end ms was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndMs < EPoverview.workData.timeNames(1)
            set(EPoverview.handles.samples.endMs,'string',sprintf('%d', EPoverview.workData.timeNames(end)+(1000/EPoverview.workData.Fs)));
            msg{1}='The new end ms cannot be smaller than the start ms.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        

        trimSamples=length(find(newEndMs < (EPoverview.workData.timeNames+(1000/EPoverview.workData.Fs))));
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[1:length(EPoverview.workData.timeNames)-trimSamples],[],[],[],[]});

        ep_editData('start');

    case 'halfRate'
        
        if mod(length(EPoverview.workData.timeNames),2) ~= 0
            msg{1}='The total number of samples must be an even number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        EPoverview.workData.baseline=EPoverview.workData.baseline/2;
        EPoverview.workData.Fs=EPoverview.workData.Fs/2;
        EPoverview.workData.recTime=EPoverview.workData.recTime/2;
        
        EPoverview.workData.timeNames=EPoverview.workData.timeNames(1:2:end);
        
        if ~isempty(EPoverview.workData.facVecT)
            newFacVecT=zeros(numPoints/2,numFacs);
            for i=1:length(EPoverview.workData.timeNames)
                newFacVecT(i,:)=(EPoverview.workData.facVecT((i-1)*2+1,:)+EPoverview.workData.facVecT(i*2,:))./2;
            end;
            EPoverview.workData.facVecT=newFacVecT;
        else
            newData=zeros(numChans,numPoints/2,numWaves,numSubs);
            for i=1:length(EPoverview.workData.timeNames)
                newData(:,i,:,:)=(EPoverview.workData.data(:,(i-1)*2+1,:,:,:,:,:)+EPoverview.workData.data(:,i*2,:,:,:,:,:))./2;
            end;
            EPoverview.workData.data=newData;
        end;

        if ~isempty(EPoverview.workData.noise)
            newData=zeros(numChans,numPoints/2,numWaves,numSubs);
            for i=1:length(EPoverview.workData.timeNames)
                newData(:,i,:,:)=(EPoverview.workData.noise(:,(i-1)*2+1,:,:,:)+EPoverview.workData.noise(:,i*2,:,:,:))./2;
            end;
            
            EPoverview.workData.noise=newData;
        end;
        
        if ~isempty(EPoverview.workData.std)
            newData=zeros(numChans,numPoints/2,numWaves,numSubs);
            for i=1:length(EPoverview.workData.timeNames)
                newData(:,i,:,:)=(EPoverview.workData.std(:,(i-1)*2+1,:,:,:,:)+EPoverview.workData.std(:,i*2,:,:,:,:))./2;
            end;
            
            EPoverview.workData.std=newData;
        end;
        
        if ~isempty(EPoverview.workData.stdCM)
            newData=zeros(numChans,numPoints/2,numWaves,numSubs);
            for i=1:length(EPoverview.workData.timeNames)
                newData(:,i,:,:)=(EPoverview.workData.stdCM(:,(i-1)*2+1,:,:,:,:)+EPoverview.workData.stdCM(:,i*2,:,:,:,:))./2;
            end;
            
            EPoverview.workData.stdCM=newData;
        end;
        
        if ~isempty(EPoverview.workData.facData)
            newData=zeros(numChans,numPoints/2,numWaves,numSubs);
            for i=1:length(EPoverview.workData.timeNames)
                newData(:,i,:,:)=(EPoverview.workData.facData(:,(i-1)*2+1,:,:,:,:,:)+EPoverview.workData.facData(:,i*2,:,:,:,:,:))./2;
            end;
            
            EPoverview.workData.facData=newData;
        end;
        
        for subject = 1:numSubs
            for wave = 1:numWaves
                for event = 1:length(EPoverview.workData.events{subject,wave})
                    EPoverview.workData.events{subject,wave}(event).sample=(EPoverview.workData.events{subject,wave}(event).sample)/2;
                end;
            end;
        end;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'twiceRate'
        %interpolate between timepoints to double the sampling rate.  The final interpolated time point will be at the very end of the epoch 
        %and will therefore only have a real timepoint in front of it.
        
        EPoverview.lastData=EPoverview.workData;
        
        EPoverview.workData.baseline=EPoverview.workData.baseline*2;
        EPoverview.workData.Fs=EPoverview.workData.Fs*2;
        EPoverview.workData.recTime=EPoverview.workData.recTime*2;
        
        sampleSize=median(diff(EPoverview.workData.timeNames))/2;
        
        newtimeNames=[EPoverview.workData.timeNames(1):sampleSize:EPoverview.workData.timeNames(end)+sampleSize]';
        
        if ~isempty(EPoverview.workData.facVecT)
            EPoverview.workData.facVecT=interp1(EPoverview.workData.timeNames,EPoverview.workData.facVecT,newtimeNames,'linear','extrap');
        else
            newData=zeros(size(EPoverview.workData.data));
            newData(1,2*numPoints,1,1,1,1,1)=0;
            for iCell=1:size(EPoverview.workData.data,3)
                for iSub=1:size(EPoverview.workData.data,4)
                    for iFac=1:size(EPoverview.workData.data,5)
                        for iFreq=1:size(EPoverview.workData.data,6)
                            for iRel=1:size(EPoverview.workData.data,7)
                                newData(:,:,iCell,iSub,iFac,iFreq,iRel)=interp1(EPoverview.workData.timeNames,squeeze(EPoverview.workData.data(:,:,iCell,iSub,iFac,iFreq,iRel))',newtimeNames,'linear','extrap')';
                            end;
                        end;
                    end;
                end;
            end;
            EPoverview.workData.data=newData;
        end;
        
        if ~isempty(EPoverview.workData.noise)
            newData=zeros(size(EPoverview.workData.noise));
            newData(1,2*numPoints,1,1,1,1,1)=0;
            for iCell=1:size(EPoverview.workData.data,3)
                for iSub=1:size(EPoverview.workData.data,4)
                    for iFac=1:size(EPoverview.workData.data,5)
                        for iFreq=1:size(EPoverview.workData.data,6)
                            for iRel=1:size(EPoverview.workData.data,7)
                                newData(:,:,iCell,iSub,iFac,iFreq,iRel)=interp1(EPoverview.workData.timeNames,squeeze(EPoverview.workData.noise(:,:,iCell,iSub,iFac,iFreq,iRel))',newtimeNames,'linear','extrap')';
                            end;
                        end;
                    end;
                end;
            end;
            EPoverview.workData.noise=newData;
        end;
        
        if ~isempty(EPoverview.workData.std)
            newData=zeros(size(EPoverview.workData.std));
            newData(1,2*numPoints,1,1,1,1,1)=0;
            for iCell=1:size(EPoverview.workData.data,3)
                for iSub=1:size(EPoverview.workData.data,4)
                    for iFac=1:size(EPoverview.workData.data,5)
                        for iFreq=1:size(EPoverview.workData.data,6)
                            for iRel=1:size(EPoverview.workData.data,7)
                                newData(:,:,iCell,iSub,iFac,iFreq,iRel)=interp1(EPoverview.workData.timeNames,squeeze(EPoverview.workData.std(:,:,iCell,iSub,iFac,iFreq,iRel))',newtimeNames,'linear','extrap')';
                            end;
                        end;
                    end;
                end;
            end;
            EPoverview.workData.std=newData;
        end;
        
        if ~isempty(EPoverview.workData.stdCM)
            newData=zeros(size(EPoverview.workData.stdCM));
            newData(1,2*numPoints,1,1,1,1,1)=0;
            for iCell=1:size(EPoverview.workData.data,3)
                for iSub=1:size(EPoverview.workData.data,4)
                    for iFac=1:size(EPoverview.workData.data,5)
                        for iFreq=1:size(EPoverview.workData.data,6)
                            for iRel=1:size(EPoverview.workData.data,7)
                                newData(:,:,iCell,iSub,iFac,iFreq,iRel)=interp1(EPoverview.workData.timeNames,squeeze(EPoverview.workData.stdCM(:,:,iCell,iSub,iFac,iFreq,iRel))',newtimeNames,'linear','extrap')';
                            end;
                        end;
                    end;
                end;
            end;
            EPoverview.workData.stdCM=newData;
        end;
        
        if ~isempty(EPoverview.workData.facData)
            newData=zeros(size(EPoverview.workData.facData));
            newData(1,2*numPoints,1,1,1,1,1)=0;
            for iCell=1:size(EPoverview.workData.data,3)
                for iSub=1:size(EPoverview.workData.data,4)
                    for iFac=1:size(EPoverview.workData.data,5)
                        for iFreq=1:size(EPoverview.workData.data,6)
                            for iRel=1:size(EPoverview.workData.data,7)
                                newData(:,:,iCell,iSub,iFac,iFreq,iRel)=interp1(EPoverview.workData.timeNames,squeeze(EPoverview.workData.facData(:,:,iCell,iSub,iFac,iFreq,iRel))',newtimeNames,'linear','extrap')';
                            end;
                        end;
                    end;
                end;
            end;
            EPoverview.workData.facData=newData;
        end;
        
        for subject = 1:numSubs
            for wave = 1:numWaves
                for event = 1:length(EPoverview.workData.events{subject,wave})
                    EPoverview.workData.events{subject,wave}(event).sample=(EPoverview.workData.events{subject,wave}(event).sample)*2-1;
                end;
            end;
        end;
        
        EPoverview.workData.timeNames=newtimeNames;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');    
    
    case 'changeStartBin'
        
        newStartBin=str2double(get(EPoverview.handles.freqs.startBin,'string'));
        if isnan(newStartBin)
            set(EPoverview.handles.freqs.startBin,'string',sprintf('%d', 1))
            msg{1}='The new start frequency bin must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartBin < 1
            set(EPoverview.handles.freqs.startBin,'string',sprintf('%d', 1))
            msg{1}='The new start frequency bin cannot be earlier than the existing start frequency bin.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartBin == 1
            msg{1}='The start frequency bin was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartBin > length(EPoverview.workData.freqNames)
            msg{1}='The new start frequency bin cannot be larger than the largest frequency bin.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        trimFreqs=newStartBin-1;

        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[],[],[],[],[1+trimFreqs:length(EPoverview.workData.freqNames)]});
        
        ep_editData('start');
        
    case 'changeEndBin'
        
        newEndBin=str2double(get(EPoverview.handles.freqs.endBin,'string'));
        if isnan(newEndBin)
            set(EPoverview.handles.freqs.endBin,'string',sprintf('%d', (length(EPoverview.workData.freqNames))))
            msg{1}='The new end frequency bin must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndBin > length(EPoverview.workData.freqNames)
            set(EPoverview.handles.freqs.endBin,'string',sprintf('%d', (length(EPoverview.workData.freqNames))))
            msg{1}='The new end frequency bin cannot be later than the existing end frequency bin.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndBin == length(EPoverview.workData.freqNames)
            msg{1}='The end frequency bin was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndBin < 1
            msg{1}='The new end frequency bin cannot be smaller than one.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        trimFreqs=length(EPoverview.workData.freqNames)-newEndBin;
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[],[],[],[],[1:length(EPoverview.workData.freqNames)-trimFreqs]});
        
        ep_editData('start');
                
    case 'changeStartFreq'
        
        newStartFreq=str2double(get(EPoverview.handles.freqs.startFreq,'string'));
        if isnan(newStartFreq)
            set(EPoverview.handles.freqs.startFreq,'string',sprintf('%d', EPoverview.workData.freqNames(1)))
            msg{1}='The new start frequency must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartFreq < EPoverview.workData.freqNames(1)
            set(EPoverview.handles.freqs.startFreq,'string',sprintf('%d', EPoverview.workData.freqNames(1)))
            msg{1}='The new start frequency cannot be earlier than the existing start frequency.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartFreq == EPoverview.workData.freqNames(1)
            msg{1}='The start frequency bin was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newStartFreq > EPoverview.workData.freqNames(end)
            msg{1}='The new start frequency cannot be larger than the largest frequency.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        trimFreqs=length(find(newStartFreq > EPoverview.workData.freqNames));

        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[],[],[],[],[1+trimFreqs:length(EPoverview.workData.freqNames)]});
        
        ep_editData('start');
        
    case 'changeEndFreq'
        
        newEndFreq=str2double(get(EPoverview.handles.freqs.endFreq,'string'));
        if isnan(newEndFreq)
            set(EPoverview.handles.freqs.endFreq,'string',sprintf('%d', (EPoverview.workData.freqNames(end))))
            msg{1}='The new end frequency must be a number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndFreq > EPoverview.workData.freqNames(end)
            set(EPoverview.handles.freqs.endFreq,'string',sprintf('%d', EPoverview.workData.freqNames(end)))
            msg{1}='The new end frequency cannot be later than the existing end frequency.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndFreq == EPoverview.workData.freqNames(end)
            msg{1}='The end frequency was not changed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if newEndFreq < EPoverview.workData.freqNames(1)
            msg{1}='The new end frequency cannot be smaller than the first one.';
            set(EPoverview.handles.freqs.endFreq,'string',sprintf('%d', EPoverview.workData.freqNames(end)))
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        trimFreqs=length(find(newEndFreq < EPoverview.workData.freqNames));
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[],[],[],[],[1:length(EPoverview.workData.freqNames)-trimFreqs]});
        
        ep_editData('start');
                
    case 'halfBinning'
        
        if mod(length(EPoverview.workData.freqNames),2) ~= 0
            msg{1}='The total number of frequency bins must be an even number.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        EPoverview.workData.freqNames=EPoverview.workData.freqNames(1:2:end);
        
        if ~isempty(EPoverview.workData.facVecF)
            newFacVecF=zeros(numFreqs/2,numFacs);
            for i=1:length(EPoverview.workData.freqNames)
                newFacVecF(i,:)=(EPoverview.workData.facVecF((i-1)*2+1,:)+EPoverview.workData.facVecF(i*2,:))./2;
            end;
        else
            newData=zeros(numChans,numPoints,numWaves,numSubs,numFacs,numFreqs/2,numRels);
            for i=1:length(EPoverview.workData.freqNames)
                %while the data are not spectral density internally, the data amplitude needs to be adjusted 
                %for what it would have been if the Hz bin size was originally larger, so summing the merged bins
                %rather than taking the mean of the merged bins.
                newData(:,:,:,:,:,i,:)=(EPoverview.workData.data(:,:,:,:,:,(i-1)*2+1,:)+EPoverview.workData.data(:,:,:,:,:,i*2,:));
            end;
        end;
        
        EPoverview.workData.data=newData;
        
        if ~isempty(EPoverview.workData.facData)
            newData=zeros(numChans,numPoints,numWaves,numSubs,numFacs,numFreqs/2,numRels);
            for i=1:length(EPoverview.workData.freqNames)
                newData(:,:,:,:,:,i,:)=(EPoverview.workData.facData(:,:,:,:,:,(i-1)*2+1,:)+EPoverview.workData.facData(:,:,:,:,:,i*2,:))./2;
            end;
            
            EPoverview.workData.facData=newData;
        end;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'Binning'
        
        EPoverview.lastData=EPoverview.workData;
        
        newBinning=str2num(get(EPoverview.handles.freqs.hBinning,'String'));
        oldBinning=EPoverview.workData.freqNames(2)-EPoverview.workData.freqNames(1);
        EPoverview.workData.freqNames=(EPoverview.workData.freqNames/oldBinning)*newBinning;
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'editEvents'
        
        EPoverview.lastData=EPoverview.workData;
        
        tableData=get(EPoverview.handles.events.hTable,'Data');
        tableNames=get(EPoverview.handles.events.hTable,'ColumnName');
        
        eventCounter=0;
        for subject=1:numSubs
            for wave=1:numWaves
                for event=1:length(EPoverview.workData.events{subject,wave})
                    eventCounter=eventCounter+1;
                    EPoverview.workData.events{subject,wave}(event).type=tableData{eventCounter,6};
                    EPoverview.workData.events{subject,wave}(event).sample=tableData{eventCounter,7};
                    EPoverview.workData.events{subject,wave}(event).value=tableData{eventCounter,8};
                    EPoverview.workData.events{subject,wave}(event).duration=tableData{eventCounter,10};
                    
                    for iKey=1:length(EPoverview.workData.events{subject,wave}(event).keys)
                        theKeyName=EPoverview.workData.events{subject,wave}(event).keys(iKey).code;
                        EPoverview.workData.events{subject,wave}(event).keys(iKey).data=tableData{eventCounter,find(strcmp(theKeyName,tableNames))};
                    end;
                end;
            end;
        end;
        
        EPoverview.events.select=cell2mat(tableData(:,1));
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
    case 'deleteEvents'
        if sum(EPoverview.events.select) ==0
            msg{1}='No events specified for deletion.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        keepEvents=not(EPoverview.events.select);
        
        eventCounter=0;
        for subject=1:numSubs
            for wave=1:numWaves
                eventList=[];
                for event=1:length(EPoverview.workData.events{subject,wave})
                    eventCounter=eventCounter+1;
                    if keepEvents(eventCounter)
                        eventList=[eventList event];
                    end;
                end;
                EPoverview.workData.events{subject,wave}=EPoverview.workData.events{subject,wave}(eventList);
            end;
        end;
        
        EPoverview.events.select=repmat(false,numEvents,1);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'clearEventWeights'
        
        EPoverview.events.weights=repmat(0,numEvents,1);
        EPoverview.events.select=repmat(false,numEvents,1);
        
        ep_editData('start');
        
    case 'allEvents'
        
        EPoverview.events.weights=repmat(1,numEvents,1);
        EPoverview.events.select=repmat(true,numEvents,1);
        
        ep_editData('start');
        
    case 'exportEvents'
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Events',[EPoverview.workData.dataName EPmain.preferences.general.specSuffix]);
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        ep_writeEventText(EPoverview.workData, [PathName filesep FileName]);
        
        ep_editData('start');
        
    case 'importEvents'
        
        [FileName,PathName,FilterIndex] = uigetfile('*.txt','Import Events',[EPoverview.workData.dataName EPmain.preferences.general.specSuffix]);
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        EPoverview.workData = ep_readEventText(EPoverview.workData, [PathName filesep FileName]);
        
        ep_editData('start');

    case 'editFacs'
        
        EPoverview.lastData=EPoverview.workData;
        
        tableData=get(EPoverview.handles.factors.hTable,'Data');
        
        if any(any(cellfun(@isempty,tableData(:,4))))
            disp('Factor names cannot be empty.');
            ep_editData('start');
            return
        end;

        %factors reordered?
        if any(cell2mat(tableData(:,1))-[1:size(tableData,1)]')
            [tableData,sortOrder] = sortrows(tableData,1);
            [tempData]=ep_reorderData(EPoverview.workData,'factors',sortOrder);
            if ~isempty(tempData)
                EPoverview.lastData=EPoverview.workData;
                EPoverview.workData=tempData;
            end;
            ep_editData('start');
        else %table data edited but not reordered
            EPoverview.factors.select=cell2mat(tableData(:,2));
            EPoverview.factors.weights=cell2mat(tableData(:,3));
            EPoverview.workData.facNames =tableData(:,4);
            EPoverview.workData.facTypes =tableData(:,5);
            set(EPoverview.handles.factors.hWeightSum,'String',sprintf('%f',sum(EPoverview.factors.weights)));
            [err]=ep_checkEPfile(EPoverview.workData);
        end;
        
    case 'deleteFacs'
        if sum(EPoverview.factors.select) ==0
            msg{1}='No factors specified for deletion.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        if sum(EPoverview.factors.select) ==numFacs
            msg{1}='There would be nothing left.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        keepFactors=not(EPoverview.factors.select);
        
        singleFacs=find(ismember(EPoverview.workData.facTypes,'SGL')) ;
        singleFacsNum=length(singleFacs);
        keepSingleFacs=zeros(singleFacsNum,1);
        for i=1:singleFacsNum
            if keepFactors(singleFacs(i))
                keepSingleFacs(i)=1;
            end;
        end;
        
        keepSingleFacs=intersect(singleFacs,find(keepFactors));
        
        [EPoverview.workData]=ep_selectData(EPoverview.workData,{[],[],[],[],keepSingleFacs,[]});
        
        EPoverview.factors.select=repmat(false,length(EPoverview.workData.facNames),1);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'addFacs'
        
        tableData=get(EPoverview.handles.factors.hTable,'Data');
        weights=cell2mat(tableData(:,3));
        
        if ~any(weights)
            msg{1}='No weights specified for forming new factor average.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(find(weights)) ==1
            msg{1}='Need to specify at least two factors to combine.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        facList=find(weights);
        facWeights=weights(facList);
        
        %make sure new channel doesn't duplicate existing channel names
        comb=1;
        while any(strcmp(['combined' num2str(comb)], EPoverview.workData.facNames))
            comb=comb+1;
        end;
 
        tempData=ep_combineData(EPoverview.workData,'factors',facList,facWeights,['combined' num2str(comb)]);
        if ~isempty(tempData)
            EPoverview.lastData=EPoverview.workData;
            EPoverview.workData=tempData;
        end;
        
        EPoverview.factors.weights=repmat(0,length(EPoverview.workData.facNames),1);
        EPoverview.factors.select=repmat(false,length(EPoverview.workData.facNames),1);
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'clearFacWeights'
        
        EPoverview.factors.weights=repmat(0,numFacs,1);
        EPoverview.factors.select=repmat(false,numFacs,1);
        
        ep_editData('start');
        
    case 'allFacs'
        
        EPoverview.factors.weights=repmat(1,numFacs,1);
        EPoverview.factors.select=repmat(true,numFacs,1);
        
        ep_editData('start');
        
    case 'appendFacs'
        
        disp('Making assumption that file format is same as that of the original file.');
        importFormat=EPoverview.workData.fileFormat;
        fileType=EPoverview.workData.dataType;
        
        textPrefs.firstRow=EPmain.preferences.general.firstRow;
        textPrefs.lastRow=EPmain.preferences.general.lastRow;
        textPrefs.firstCol=EPmain.preferences.general.firstCol;
        textPrefs.lastCol=EPmain.preferences.general.lastCol;
        textPrefs.orientation=EPmain.preferences.general.orientation;
        
        inArg{1}='format';
        inArg{2}=importFormat;
        inArg{3}='type';
        inArg{4}=fileType;
        inArg{5}='elecPrefs';
        inArg{6}=EPmain.preferences.general.rotateHead;
        inArg{7}='textPrefs';
        inArg{8}=textPrefs;
        inArg{9}='screenSize';
        inArg{10}=EPmain.scrsz;
        inArg{11}='FontSize';
        inArg{12}=EPmain.fontsize;
        
        [newData]=ep_readData(inArg);
        
        if isempty(newData)
            return
        end;
        
        tableData=get(EPoverview.handles.factors.hTable,'Data');
        
        if ~strcmp(fileType,newData.dataType)
            msg{1}='The file types do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
         
        if length(EPoverview.workData.cellNames) ~= length(newData.cellNames)
            msg{1}='The numbers of cells do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(EPoverview.workData.chanNames) ~= length(newData.chanNames)
            msg{1}='The numbers of channels do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(EPoverview.workData.timeNames) ~= length(newData.timeNames)
            msg{1}='The numbers of time points do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if length(EPoverview.workData.subames) ~= length(newData.subNames)
            msg{1}='The numbers of subjects do not match.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if size(EPoverview.workData.facVecS,2) ~= size(newData.facVecS,2)
            msg{1}='The two datasets are not equally factor compressed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if size(EPoverview.workData.facVecT,2) ~= size(newData.facVecT,2)
            msg{1}='The two datasets are not equally factor compressed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        EPoverview.lastData=EPoverview.workData;
        
        if ~isempty(EPoverview.workData.facVecS)
            EPoverview.workData.facVecS(:,end+1:length(newData.facNames))=newData.facVecS;
        end;
        
        if ~isempty(EPoverview.workData.facVecT)
            EPoverview.workData.facVecT(:,end+1:length(newData.facNames))=newData.facVecT;
        end;
        
        if ~isempty(EPoverview.workData.facVar)
            EPoverview.workData.facVar(:,end+1:length(newData.facNames))=newData.facVar;
        end;
        
        if ~isempty(EPoverview.workData.facVarQ)
            EPoverview.workData.facVarQ(:,end+1:length(newData.facNames))=newData.facVarQ;
        end;
        
        EPoverview.workData.data(:,:,:,:,end+1:end+size(newData.data,5),:,:)=newData.data;
        EPoverview.workData.facData(:,:,:,:,end+1:end+size(newData.facData,5),:,:)=newData.facData;
        EPoverview.workData.facNames(end+1:end+length(newData.facNames),1)=newData.facNames;
        EPoverview.workData.facTypes(end+1:end+length(newData.facNames),1)=newData.facTypes;
        
        EPoverview.factors.select=[EPoverview.factors.select repmat(false,length(newData.facNames),1)];
        EPoverview.factors.weights=[EPoverview.factors.weights repmat(0,length(newData.facNames),1)];
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
        
    case 'exportFacs'
        
        tableData=get(EPoverview.handles.factors.hTable,'Data');
        tableData=tableData(:,4:end);
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Table');
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        [pathstr, name, ext] = fileparts(FileName);
        if ~strcmp(ext,'.txt')
            FileName=[FileName '.txt'];
        end;
        
        fid=fopen([PathName FileName],'w');
        
        for i=1:size(tableData,1)
            for j=1:size(tableData,2)
                theCell=tableData{i,j};
                if isnumeric(theCell)
                    fprintf(fid,'%d', theCell);
                else
                    fprintf(fid,'%s', theCell);
                end;
                if j < size(tableData,2)
                    fprintf(fid,'\t');
                end;
            end;
            if i < size(tableData,1)
                fprintf(fid,'\r');
            end;
        end;
        
        ep_editData('start');
        
    case 'exportQC'
        
        tableData=get(EPoverview.handles.QC.hTable,'Data');
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Table');
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        [pathstr, name, ext] = fileparts(FileName);
        if ~strcmp(ext,'.txt')
            FileName=[FileName '.txt'];
        end;
        
        fid=fopen([PathName FileName],'w');
        
        for i=1:size(tableData,1)
            for j=1:size(tableData,2)
                theCell=tableData{i,j};
                if isnumeric(theCell)
                    fprintf(fid,'%d', theCell);
                else
                    fprintf(fid,'%s', theCell);
                end;
                if j < size(tableData,2)
                    fprintf(fid,'\t');
                end;
            end;
            if i < size(tableData,1)
                fprintf(fid,'\r');
            end;
        end;
        
        ep_editData('start');
        
    case 'exportPCA'
        
        tableData=get(EPoverview.handles.PCA.hTable,'Data');
        
        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Export Table');
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        [pathstr, name, ext] = fileparts(FileName);
        if ~strcmp(ext,'.txt')
            FileName=[FileName '.txt'];
        end;
        
        fid=fopen([PathName FileName],'w');
        
        for i=1:size(tableData,1)
            for j=1:size(tableData,2)
                theCell=tableData{i,j};
                if isnumeric(theCell)
                    fprintf(fid,'%d', theCell);
                else
                    fprintf(fid,'%s', theCell);
                end;
                if j < size(tableData,2)
                    fprintf(fid,'\t');
                end;
            end;
            if i < size(tableData,1)
                fprintf(fid,'\r');
            end;
        end;
        
        ep_editData('start');
        
    case 'importTrials'
        
        [FileName,PathName,FilterIndex] = uigetfile({'*.txt';'*.epoc'},'Cell Names');
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                return
            end;
        end;
        
        fid=fopen([PathName FileName],'r');
        tempVar=fgetl(fid);
        delim='\t';
        if ~isempty(strfind(tempVar,',')) && isempty(strfind(tempVar,'\t'))
            delim=','; %if there are commas and no tabs, assume it is a comma-delimited file.
        elseif ~isempty(strfind(tempVar,' ')) && isempty(strfind(tempVar,'\t'))
            delim=' '; %if there are spaces and no tabs, assume it is a space-delimited file.
        end;
        
        numcols=length(regexp(tempVar,delim))+1; %determine number of columns based on tab markers
        if regexp(tempVar,[delim '$'])
            numcols=numcols-1; %if there is an extra tab at the end of the line, drop it.
        end;
        
        frewind(fid);
        theData=textscan(fid, [repmat('%s',1,numcols)],'Delimiter',delim);
        fclose(fid);
        if size(theData{1},1) ~= numWaves
            msg{1}='Number of names in text file does not match number of trials in dataset.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        EPoverview.workData.cellNames=theData{1};
        
        %rename the combined trial names
        uniqueCells=unique(EPoverview.workData.cellNames,'first');
        numCells=length(uniqueCells);
        trialCount=ones(numCells,1);
        for wave=1:length(EPoverview.workData.cellNames)
            theCell=strcmp(EPoverview.workData.cellNames(wave),uniqueCells);
            EPoverview.workData.trialNames(wave)=trialCount(theCell);
            trialCount(theCell)=trialCount(theCell)+1;
        end;
        %group the trials by cell
        [uniqueCellNames m cellIDs]=unique(EPoverview.workData.cellNames); %get number of trials in each cell
        cellNums=hist(cellIDs,unique(cellIDs));
        trialCounter=zeros(numCells,1);
        cellNames2=cell(size(EPoverview.workData.cellNames));
        cellTypes2=cell(size(EPoverview.workData.cellTypes));
        trialNames2=zeros(size(EPoverview.workData.trialNames));
        trialSpecs2=cell(size(EPoverview.workData.trialSpecs));
        events2=cell(size(EPoverview.workData.events));
        avgNum2=zeros(size(EPoverview.workData.avgNum));
        covNum2=zeros(size(EPoverview.workData.covNum));
        subNum2=zeros(size(EPoverview.workData.subNum));
        data2=zeros(size(EPoverview.workData.data));
        if ~isempty(EPoverview.workData.facData)
            facData2=zeros(size(EPoverview.workData.facData));
        end;
        if ~isempty(EPoverview.workData.noise)
            noise2=zeros(size(EPoverview.workData.noise));
        end;
        if ~isempty(EPoverview.workData.std)
            std2=zeros(size(EPoverview.workData.std));
        end;
        if ~isempty(EPoverview.workData.stdCM)
            stdCM2=zeros(size(EPoverview.workData.stdCM));
        end;
        
        blinkTrial2=zeros(size(EPoverview.workData.analysis.blinkTrial));
        saccadeTrial2=zeros(size(EPoverview.workData.analysis.saccadeTrial));
        saccadeOnset2=zeros(size(EPoverview.workData.analysis.saccadeOnset));
        moveTrial2=zeros(size(EPoverview.workData.analysis.moveTrial));
        badTrials2=zeros(size(EPoverview.workData.analysis.badTrials));
        badChans2=zeros(size(EPoverview.workData.analysis.badChans));
        
        for theTrial = 1:length(cellIDs)
            newCell=cellIDs(theTrial);
            trialCounter(newCell)=trialCounter(newCell)+1;
            data2(:,:,sum(cellNums(1:newCell-1))+trialCounter(newCell),:,:,:)=EPoverview.workData.data(:,:,theTrial,:,:,:,:);
            if ~isempty(EPoverview.workData.facData)
                facData2(:,:,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.facData(:,:,theTrial,:,:,:,:);
            end;
            if ~isempty(EPoverview.workData.noise)
                noise2(:,:,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.noise(:,:,theTrial,:,:);
            end;
            if ~isempty(EPoverview.workData.std)
                std2(:,:,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.std(:,:,theTrial,:,:,:);
            end;
            if ~isempty(EPoverview.workData.stdCM)
                stdCM2(:,:,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.stdCM(:,:,theTrial,:,:,:);
            end;
            
            cellNames2(sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.cellNames(theTrial);
            cellTypes2(sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.cellTypes(theTrial);
            trialNames2(sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.trialNames(theTrial);
            trialSpecs2(sum(cellNums(1:newCell-1))+trialCounter(newCell),:)=EPoverview.workData.trialSpecs(theTrial,:);
            events2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.events(1,theTrial);
            avgNum2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.avgNum(1,theTrial);
            covNum2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.covNum(1,theTrial);
            subNum2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))=EPoverview.workData.subNum(1,theTrial);
            blinkTrial2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))= EPoverview.workData.analysis.blinkTrial(1,theTrial);
            saccadeTrial2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))= EPoverview.workData.analysis.saccadeTrial(1,theTrial);
            saccadeOnset2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))= EPoverview.workData.analysis.saccadeOnset(1,theTrial);
            moveTrial2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))= EPoverview.workData.analysis.moveTrial(1,theTrial);
            badTrials2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell))= EPoverview.workData.analysis.badTrials(1,theTrial);
            badChans2(1,sum(cellNums(1:newCell-1))+trialCounter(newCell),:)= EPoverview.workData.analysis.badChans(1,theTrial,:);
        end;
        
        EPoverview.workData.cellNames=cellNames2;
        EPoverview.workData.cellTypes=cellTypes2;
        EPoverview.workData.trialNames=trialNames2;
        EPoverview.workData.trialSpecs=trialSpecs2;
        EPoverview.workData.events=events2;
        EPoverview.workData.avgNum=avgNum2;
        EPoverview.workData.covNum=covNum2;
        EPoverview.workData.subNum=subNum2;
        EPoverview.workData.data=data2;
        if ~isempty(EPoverview.workData.facData)
            EPoverview.workData.data=facData2;
        end;
        if ~isempty(EPoverview.workData.noise)
            EPoverview.workData.noise=noise2;
        end;
        if ~isempty(EPoverview.workData.std)
            EPoverview.workData.std=std2;
        end;
        if ~isempty(EPoverview.workData.stdCM)
            EPoverview.workData.std=stdCM2;
        end;
        
        EPoverview.workData.analysis.blinkTrial=blinkTrial2;
        EPoverview.workData.analysis.saccadeTrial=saccadeTrial2;
        EPoverview.workData.analysis.saccadeOnset=saccadeOnset2;
        EPoverview.workData.analysis.moveTrial=moveTrial2;
        EPoverview.workData.analysis.badTrials=badTrials2;
        EPoverview.workData.analysis.badChans=badChans2;
        
        EPoverview.cells.select=repmat(false,numCells,1);
        EPoverview.cells.weights=repmat(0,numCells,1);
        
        ep_editData('start');
    case 'done'
        try
            EPver=ver('EP_Toolkit');
        catch
            EPver='unavailable'; %workaround for bug in earlier version of Matlab
        end;
        EPoverview.workData.EPver=EPver;
        EPoverview.workData.ver=ver;
        EPoverview.workData.date=date;
        EPoverview.workData.history(end+1)={'ep_editData'};
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        EPdata=EPoverview.workData;
        
        if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
            msg{1}='The work directory cannot be found.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        set(EPoverview.handles.hDone,'enable','off');
        drawnow;
        
        if ~strcmp(EPoverview.workData.dataName,EPdataset.dataset(EPoverview.dataset).dataName)
            button = questdlg('Did you want this to be a new separate dataset, leaving the original unchanged?','Duplicate dataset?','New','Rename','New');
            if strcmp(button,'New')
                sameName=1;
                suffix=0;
                nameSuffix=EPoverview.workData.dataName;
                while sameName
                    sameName=0;
                    for i=1:length(EPdataset.dataset)
                        if strcmp(EPdataset.dataset(i).dataName,nameSuffix)
                            sameName=1;
                        end;
                    end;
                    if sameName
                        suffix=suffix+1;
                        nameSuffix=[EPoverview.workData.dataName '-' num2str(suffix)];
                    end;
                end;
                EPoverview.workData.dataName=nameSuffix;
                EPdataset=ep_saveEPdataset(EPdataset,EPdata,length(EPdataset.dataset)+1,'no');
            else
                delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(EPoverview.dataset).dataName '.mat']);
                EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],EPoverview.dataset));
                EPdataset=ep_saveEPdataset(EPdataset,EPdata,EPoverview.dataset,'no');
            end;
        else
            delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(EPoverview.dataset).dataName '.mat']);
            EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],EPoverview.dataset));
            EPdataset=ep_saveEPdataset(EPdataset,EPdata,EPoverview.dataset,'no');
        end;
        
        if ~isempty(EPoverview.workData.facNames)
            disp('Note: edits made to a PCA file will not be reflected in the follow-up PCA in a two-step PCA procedure as the PCA information is kept separate from the edited information.');
        end;
        
        close('EditData');
        EPoverview=[];
        EPoverview='done';
        ep('start');
        
    case 'undo'
        
        if isempty(EPoverview.lastData)
            msg{1}='Nothing to undo.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        tempVar=EPoverview.workData;
        EPoverview.workData=EPoverview.lastData;
        EPoverview.lastData=tempVar;
        
        numPoints=length(EPoverview.workData.timeNames);
        numChans=length(EPoverview.workData.chanNames);
        numCells=length(unique(EPoverview.workData.cellNames));
        numWaves=length(EPoverview.workData.cellNames);
        numSubs=length(EPoverview.workData.subNames);
        numSpecs=length(EPoverview.workData.subjectSpecNames);
        numImplicit=length(EPoverview.workData.implicit);
        numFacs=length(EPoverview.workData.facNames);
        
        numEvents=0;
        for subject=1:numSubs
            for wave=1:numWaves
                numEvents=numEvents+length(EPoverview.workData.events{subject,wave});
            end;
        end;
        
        EPoverview.subjects.select=repmat(false,numSubs,1);
        EPoverview.cells.select=repmat(false,numCells,1);
        EPoverview.channels.select=repmat(false,numChans+numImplicit,1);
        EPoverview.subjects.weights=repmat(0,numSubs,1);
        EPoverview.cells.weights=repmat(0,numCells,1);
        EPoverview.channels.weights=repmat(0,numChans+numImplicit,1);
        EPoverview.events.select=repmat(false,numEvents,1);
        EPoverview.events.weights=repmat(0,numEvents,1);
        EPoverview.factors.select=repmat(false,numFacs,1);
        EPoverview.factors.weights=repmat(0,numFacs,1);
        
        ep_editData('start');
        
    case 'cancel'
        close('EditData');
        EPoverview=[];
        
    case 'new'
        
        try
            EPver=ver('EP_Toolkit');
        catch
            EPver='unavailable'; %workaround for bug in earlier version of Matlab
        end;
        EPoverview.workData.EPver=EPver;
        EPoverview.workData.ver=ver;
        EPoverview.workData.date=date;
        EPoverview.workData.history(end+1)={'ep_editData'};
        
        [err]=ep_checkEPfile(EPoverview.workData);
        
        EPdata=EPoverview.workData;
        
        if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
            msg{1}='The work directory cannot be found.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        set(EPoverview.handles.hNew,'enable','off');
        drawnow;
        
        sameName=1;
        suffix=0;
        nameSuffix=EPoverview.workData.dataName;
        while sameName
            sameName=0;
            for i=1:length(EPdataset.dataset)
                if strcmp(EPdataset.dataset(i).dataName,nameSuffix)
                    sameName=1;
                end;
            end;
            if sameName
                suffix=suffix+1;
                nameSuffix=[EPoverview.workData.dataName '-' num2str(suffix)];
            end;
        end;
        EPoverview.workData.dataName=nameSuffix;

        EPdataset=ep_saveEPdataset(EPdataset,EPdata,length(EPdataset.dataset)+1,'no');
         
        if ~isempty(EPoverview.workData.facNames)
            disp('Note: edits made to a PCA file will not be reflected in the follow-up PCA in a two-step PCA procedure as the PCA information is kept separate from the edited information.');
        end;
        
        close('EditData');
        EPoverview=[];
        ep('start');
        
    otherwise
        disp(['Don''t recognize the task ' varargin{1}]);
        
end

if ~strcmp(varargin{1},{'cancel','done'})
    refresh;
    drawnow;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editCells(src,eventdata)
%change the contents of the Cells table

global EPoverview

EPoverview.lastData=EPoverview.workData;

tableData=get(EPoverview.handles.cells.hTable,'Data');

theRow=eventdata.Indices(1);
theCol=eventdata.Indices(2);

numCells=length(unique(EPoverview.workData.cellNames));

switch theCol
    case 1 %reorder cells
        
        [tableData2,index] = sortrows(tableData,1);
        
        counter=0;
        sortOrder=zeros(length(EPoverview.workData.cellNames),1);
        for i=1:numCells %going through new order of reordered cells
            trialIndex=find(strcmp(tableData{index(i),4},EPoverview.workData.cellNames));
            numTrials=length(trialIndex); %number of trials in the cell
            sortOrder(counter+1:counter+numTrials)=trialIndex;
            counter=counter+numTrials;
        end;
        
        [tempData]=ep_reorderData(EPoverview.workData,'cells',sortOrder);
        if ~isempty(tempData)
            EPoverview.lastData=EPoverview.workData;
            EPoverview.workData=tempData;
        end;
        
        ep_editData('start');
        
    case 2 %select checked
        EPoverview.cells.select=cell2mat(tableData(:,2));
    case 3 %weight set
        EPoverview.cells.weights=cell2mat(tableData(:,3));
        ep_editData('start');
    case 4 %cell name changed
        newCellName=eventdata.NewData;
        oldCellName=eventdata.PreviousData;
        if isempty(newCellName)
            disp('Cell names cannot be empty.');
            ep_editData('start');
            return
        end;
        if ((length(unique(EPoverview.workData.cellNames)) ~= length(EPoverview.workData.cellNames))) && strcmp(EPoverview.workData.dataType,'average')
            disp('Cell names cannot be duplicated except in single trial data.');
            ep_editData('start');
            return
        end;
        if strcmp(EPoverview.workData.dataType,'single_trial')
            cellIndex=find(strcmp(oldCellName,EPoverview.workData.cellNames));
            EPoverview.workData.trialNames(cellIndex)=EPoverview.workData.trialNames(cellIndex)+length(find(strcmp(newCellName,EPoverview.workData.cellNames)));
            [EPoverview.workData.cellNames{cellIndex}]=deal(newCellName);
        else
            EPoverview.workData.cellNames{theRow}=newCellName;
        end;
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
    case 5 %cell type changed
        newCellType=eventdata.NewData;
        oldCellType=eventdata.PreviousData;
        if isempty(newCellType)
            disp('Cell types cannot be empty.');
            ep_editData('start');
            return
        end;
        EPoverview.workData.cellTypes{theRow}=newCellType;
        [err]=ep_checkEPfile(EPoverview.workData);
        
        ep_editData('start');
end;



