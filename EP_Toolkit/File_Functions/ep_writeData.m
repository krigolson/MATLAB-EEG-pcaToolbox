function [err]=ep_writeData(EPdata,outFileName,eventSuffix,subjectSpecSuffix,format,adds,prefs);
%  [err]=ep_writeData(EPdata,outFileName,eventSuffix,format,adds,prefs);
%       Writes EP data to output files of various available formats.  The calling routine is responsible for checking
%       for file overwriting.
%
%Inputs:
%  EPdata         : Structured array with the data and accompanying information.  See readData.
%  outFileName  : Name of the file to be imported.  No suffix needed as it will be added.
%                 Will assume active directory.
%  eventSuffix  : The suffix added to event text files (e.g., "_evt.txt").
%  subjectSpecSuffix : The suffix added to subject spec text files (e.g., "_sub.txt").
%
% Optional---------------------
%  format       : File format of the output file.  Can be 'egi_egis', 'egi_egia', 'egi_sbin', 'text', or 'ep_mat'.
%  adds         : Cell string array of types of information to be written out:
%                   SGLchan (single channels, including EEG, MGM, MGA, MGP, ANS, and REF)
%                   REGchan (regional channels)
%                   SGLcell (single cells)
%                   CMBcell (combined cells)
%                   RAW (single trial data)
%                   AVG (averaged data)
%                   GAV (grand average data)
%                   SGLfac (single factors)
%                   CMBfac (combined factors)
%  prefs        : Cell array of preferences
%                   'montage off' means do not add montage information to EGIS format files.

%
%Outputs:
%  err          : Error code.
%
%  Outputs the EP data as a data file.  Cells are sorted alphabetically in order of cell names.  For single_trial
%  type the trials are grouped together by cell.  For average type the subjects are grouped together by cell.
%  Regional channel averages stripped out for EGI formats.
%  Also, in EGI formats, subject adds are saved in a separate file to avoid confusion.
%
%History:
%  by Joseph Dien (3/10/09)
%  jdien07@mac.com
%
%
% bugfix 4/27/09 JD
% Fix for file spec fields being all the same for a trial in EGIS session files.
%
% bugfix 9/14/09 JD
% Crash when saving factor files in simple binary format.
% Out of memory errors when saving factor files in simple binary format.
%
% modified & bugfix 7/22/09 JD
% Eliminated the latency field from the event structure.  Added support for 5D data field.
% Added option to control types of combined data output.  Added support for EP files
% Will only assume active directory is the location for the output file if no path is
% provided in the output file name.  Adds subject IDs to EGIS average output if not already available.
% Handles files with adds.  Consolidated code for factor files.
% EGIS subject specs converted to strings.
%
% bugfix 9/17/09 JD
% Not saving output file when format is text and data type is factors.
%
% bugfix 10/29/09 JD
% Crash when saving simple binary format files and the data file has information in the subject spec fields.
%
% bugfix 11/4/09 JD
% Crash when saving simple binary format files and the data file was not originally an EGIS format file.
% Crash when saving single trial simple binary format file and there are more than one events per trial.
% Single trial simple binary files being saved as ".bin.bin"
% Thanks to Grega Repov.
%
% bugfix 1/31/10 JD
% Fixed not saving event data to single_trial simple binary files and putting events at the very first sample for all simple binary files.
%
% bugfix 2/8/10 JD
% Stopped adding "_ep" suffix in this function so that it can instead be done at the "ep" level where error-checking
% happens.
%
% bugfix & modified 2/28/10 JD
% When writing EGIS files, cell names and experiment names are terminated so that are not padded out with spaces when
% read by some programs.
% Added option to turn off adding montage information to EGIS file format files due to incompatibility issues with some
% versions of NetStation.
% Fixed crash when saving EGIS file format file with a "sex" subject spec field.
%
% modified 3/27/10 JD
% Added ability to save data using .set file format.
%
% bugfix & modified 5/25/10 JD
% When saving EEGlab file formats, cells are now saved as separate .set files.
%
% modified 6/16/10 JD
% Changed the suffix of EP files to ".ept".
%
% modified 8/25/10 JD
% Now supports writing out unsegmented simple binary files.
%
%  modified 10/12/10 JD
%  For continuous files, data now divided into one second epochs and can be artifact rejected in an epochwise fashion
%  in same fashion as segmented data.  Bad trials deleted during file saving to formats other than EP.
%
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
%
% bugfix 2/15/10 JD
% Fixed crash when saving continuous simple binary or EP file with events and with no bad trials.
%
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% bugfix 5/30/12 JD
% Fixed crash when saving continuous simple binary files with bad segments
% and changed so that bad data segments are marked with 1000 microvolts
% rather than just deleted.
%
% modified 2/3/13 JD
% Added writing ERPlab format files.
%
% bugfix 4/5/13 JD
% Fixed ERPlab fields ERP.ntrials.rejected and ERP.ntrials.invalid to be vectors of zeros rather than empty set.
% Due Joseph Orr.
%
% modified 10/21/13 JD
% Added support for writing nTrials fields from ERPlab files.
%
% bugfix 10/24/13 JD
% Fixed warning message when saving file in simple binary format and there are events with fractional durations.
%
% modified 10/28/13 JD
% Dropped setting bad one second epochs of continuous data to 1000mv when marked bad.
%
% bugfix 10/29/13 JD
% Fixed EEGlab .set files not including study name in file name.
% Fixed not checking to see if file name already exists when saving an EEGLab .study file.
%
% modified 3/18/14 JD
% Improved saving data in text format so will check for overwriting and for FFT data will save as freq by chan.
% Added support for saving data in Neuromag FIFF file format.
%
% modified 3/24/14 JD
% Added support for bad channels for Neuromag FIFF file format.
% Added cov field.
%
% modified 4/14/14 JD
% Use covNum field when saving FIFF files rather than avgNum.
%
% modified 4/24/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure to text output.
%
% bugfix 7/24/14 JD
% Fixed crash when writing out simple binary files when there are events with empty fields.
%
% bugfix 8/27/14 JD
% Fixed crash for files where event values are all numbers rather than
% strings.
%
% modified 9/17/14 JD
% Adding event output for files that don't support events.
%
% bugfix 12/22/14 JD
% Fixed crash when saving factor files with more than one factor in simple
% binary file format.
% Fixed not handling event times correctly when saving files in simple
% binary form and event time samples are decimals (warning was being
% displayed by Matlab).
%
% modified 10/26/15 JD
% EGIS cells no longer sorted alphabetically when saved.
%
% bugfix 3/4/16 JD
% Fixed error when saving text files.
% Fixed error when saving separate event file for text files.
%
% modified 11/5/16 JD
% Added support for writing out subject spec text files.
%
% bugfix 12/1/17 JD
% Fixed dropping last part of file name after a period when writing out non-EGIS files.
%
% bugfix 4/25/18 JD
% If trying to save a file in EP format and preferences are set to v6 or v7 and it is over 2GB, then instead of not saving it will now save in v7.3 format.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[err]=ep_checkEPfile(EPdata);
if err
    return;
end;

if ~exist('eventSuffix','var')
    eventSuffix = '_evt.txt';
end;

if ~exist('subjectSpecSuffix','var')
    eventSuffix = '_sub.txt';
end;

if ~exist('format','var')
    format = EPdata.fileFormat;
end;

if ~exist('adds','var')
    adds = {'SGLchan','REGchan','SGLcell','CMBcell','RAW','AVG','GAV','SGLfac','CMBfac'};
end;

if ~exist('prefs','var')
    prefs={};
end;

[EPdata]=ep_stripAdds(EPdata,adds);
if isempty(EPdata.data)
    err=1;
    msg{1}='No data were left after stripping out unwanted data types.';
    [msg]=ep_errorMsg(msg);
    return
end;

type=EPdata.dataType;

[pathstr, fileName, ext]=fileparts(outFileName);
if isempty(pathstr)
    pathstr=pwd;
end;
outFileName=[pathstr filesep fileName ext];

%drop bad epochs if continuous data
% if strcmp(type,'continuous') && ~strcmp(format,'ep_mat')
%     if strcmp(EPdata.dataType,'continuous')
%         numEpochs=floor(size(EPdata.data,2)/EPdata.Fs); %excess time points are tacked onto final epoch
%         if numEpochs == 0
%             numEpochs =1;
%         end;
%     end;
%     counter=0;
%     for epoch=1:numEpochs
%         if EPdata.analysis.badTrials(epoch-counter)
%             if epoch == numEpochs %excess time points are tacked onto final epoch
%                 epochPoints=(epoch-1-counter)*EPdata.Fs+1:(epoch-1-counter)*EPdata.Fs-size(EPdata.data,2);
%             else
%                 epochPoints=(epoch-1-counter)*EPdata.Fs+1:(epoch-counter)*EPdata.Fs;
%             end;
%
%             EPdata.data(:,epochPoints,:,:,:)=1000; %mark bad data with huge microvolts
% %             if ~isempty(EPdata.noise)
% %                 EPdata.noise(:,epochPoints,:,:,:)=[];
% %             end;
% %             if ~isempty(EPdata.std)
% %                 EPdata.std(:,epochPoints,:,:,:)=[];
% %             end;
% %             EPdata.timeNames(epochPoints)=[];
%
% %             droppedEvents=[];
% %             for theEvent=1:length(EPdata.events{1})
% %                 startSample=EPdata.events{1}(theEvent).sample;
% %                 endSample=EPdata.events{1}(theEvent).sample+EPdata.events{1}(theEvent).duration-1;
% %                 if epoch == numEpochs %excess time points are tacked onto final epoch
% %                     eventPoints=(epoch-1)*EPdata.Fs+1:(epoch-1)*EPdata.Fs-size(EPdata.data,2);
% %                 else
% %                     eventPoints=(epoch-1)*EPdata.Fs+1:(epoch)*EPdata.Fs;
% %                 end;
% %
% %                 if (startSample >= min(eventPoints) && (startSample <= max(eventPoints)))
% %                     if endSample <= max(eventPoints)
% %                         droppedEvents(end+1)=theEvent;
% %                     else
% %                         EPdata.events{1}(theEvent).sample=min(epochPoints);
% %                         EPdata.events{1}(theEvent).duration=EPdata.events{1}(theEvent).duration-(max(epochPoints)-startSample+1);
% %                     end;
% %                 elseif (endSample >= min(eventPoints)) && (endSample <= max(eventPoints))
% %                     EPdata.events{1}(theEvent).duration=EPdata.events{1}(theEvent).duration-(endSample-min(epochPoints)+1);
% %                 end;
% %             end;
% %             EPdata.events{1}(droppedEvents)=[];
% %             EPdata.analysis.badChans(:,epoch-counter,:)=[];
% %             EPdata.analysis.moveTrial(:,epoch-counter)=[];
% %             EPdata.analysis.blinkTrial(:,epoch-counter)=[];
% %             EPdata.analysis.saccadeTrial(:,epoch-counter)=[];
% %             EPdata.analysis.saccadeOnset(:,epoch-counter)=[];
% %             EPdata.analysis.badTrials(:,epoch-counter)=[];
%              counter=counter+1;
%         end;
%     end;
%     if counter
%         disp('Note: marking bad one-second segments in data with 1000 µv time points.');
%     end;
% end;

switch format
    case 'ep_mat'
        
        sameName=1;
        theNumber=0;
        [pathstr, fileName, ext] = fileparts(outFileName);
        fileNameStem=fileName;
        while sameName
            sameName=0;
            if exist([pathstr filesep fileName '.ept'],'file')
                sameName=1;
            end;
            if sameName
                theNumber=theNumber+1;
                fileName=[fileNameStem '-' num2str(theNumber)];
            end;
        end;
        
        warning('');
        warning('off','MATLAB:save:sizeTooBigForMATFile')
        save('-mat', [pathstr filesep fileName '.ept'], 'EPdata');
        warning('on')
        if ~isempty(lastwarn)
            [msg,msgID] = lastwarn;
            if strcmp(msgID,'MATLAB:save:sizeTooBigForMATFile')
                disp('The file was too big for v6 or v7 so will save in v7.3 format.')
                save('-mat', [pathstr filesep fileName '.ept'], 'EPdata','-v7.3');
            end;
        end;
        
    case 'egi_egis'
        if ~strcmp(type,'single_trial')
            err=1;
            msg{1}='File format is inconsistent with the file type.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if ~isempty(EPdata.freqNames)
            err=1;
            msg{1}='EGIS format cannot represent frequency data.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        [pathstr, fileName, ext] = fileparts(outFileName);
        if ~isempty(find(ismember(EPdata.subTypes,{'AVG','GAV'}), 1)) && ismember('RAW',EPdata.subTypes)
            adds2=setdiff(adds,'RAW');
            [err]=ep_writeData(EPdata,[pathstr filesep fileName '_AVG' ext],eventSuffix,subjectSpecSuffix,'egi_egia',adds2); %save average data as a separate EGIS AVE file
            adds=setdiff(adds,{'AVG','GAV'});
            [EPdata]=ep_stripAdds(EPdata,adds); %remove the average data which has now been saved.
        end;
        
        if isempty(EPdata.data) %return if the single-trial data had been deleted.
            return
        end;
        
        numChans= length(EPdata.chanNames);
        numPoints = length(EPdata.timeNames);
        numSubs = length(EPdata.subNames);
        numFacs= length(EPdata.facNames);
        numTrials = length(EPdata.trialNames);
        
        ses_hdr_offsets_v;
        
        SampleRate = EPdata.Fs;
        
        expname=EPdata.ename;
        
        if size(expname,1) > size(expname,2)
            expname=expname'; %make sure the experiment name is horizontal
        end;
        
        [cellNames m cellIDs]=unique(EPdata.cellNames,'stable'); %get number of trials in each cell
        numCells=length(cellNames);
        cellNums=hist(cellIDs,unique(cellIDs));
        
        fhdr=[];
        chdr=[];
        fhdr(NChan) = numChans;
        fhdr(NCells) = numCells;
        fhdr(HdrVer) = 3;
        fcom=[];
        ftext=[];
        
        for theSubjectSpec = 1:length(EPdata.subjectSpecNames)
            switch EPdata.subjectSpecNames{theSubjectSpec}
                case 'RunDateMo'
                    fhdr(RunDateMo)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'RunDateDay';
                    fhdr(RunDateDay)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'RunDateYr';
                    fhdr(RunDateYr)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'RunTimeHr';
                    fhdr(RunTimeHr)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'RunTimeMin';
                    fhdr(RunTimeMin)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'RunTimeSec';
                    fhdr(RunTimeSec)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'SubjID';
                    fhdr(SubjID)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'Handed';
                    fhdr(Handed)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'Sex';
                    fhdr(Sex)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'Age';
                    fhdr(Age)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'ExperID';
                    fhdr(ExperID)=str2num(EPdata.subjectSpecs{1,theSubjectSpec});
                case 'Comment';
                    fcom=EPdata.subjectSpecs{1,theSubjectSpec};
                case 'Text';
                    ftext=EPdata.subjectSpecs{1,theSubjectSpec};
                otherwise
            end;
        end;
        
        ename = [expname char(0) blanks(80)];
        ename = ename(1:80);
        fhdr(LComment)=length(fcom);
        fhdr(LText)=length(ftext);
        
        numSubSpecs = size(EPdata.trialSpecs,2);
        
        %Construct cell headers
        chdr = zeros(numCells,(6+numSubSpecs*max(cellNums)));
        for i = 1:numCells
            chdr(i,CellID) = 90 + i;
            chdr(i,NObs) = cellNums(i);
            chdr(i,NPoints) = numPoints;
            chdr(i,SampRate) = SampleRate;
            chdr(i,LSpec) = numSubSpecs*2;
        end;
        
        L = ((chdr(:,LSpec).*chdr(:,NObs))'+90);
        fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))) = L';
        nlcellhdr = fhdr(NCells)*2;
        lhdr = 130+nlcellhdr+ fhdr(NChan) * 4 + sum(L)+fhdr(LComment)+fhdr(LText);
        fhdr(LPad) = (512-rem(lhdr,512));
        lhdr = lhdr+fhdr(LPad);
        fhdr(LHeader) = lhdr;
        fhdr(LData) = sum((chdr(:,NPoints).*chdr(:,NObs))*fhdr(NChan))*2;
        
        if lhdr > 32767
            err=1;
            msg{1}='Header size is too large.  You need to shrink the file or switch to a different file format.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        for cell=1:fhdr(NCells)
            theName=[cellNames{cell}(1,:) char(0) blanks(80)];
            cnames(cell,:)=theName(1:80);
        end;
        czeros=zeros(numChans,1);
        cgains=ones(numChans,1)*2000;
        fhdr(BytOrd) = 16909060; %always write out EGIS files big-endian
        
        if ~isempty(EPdata.facNames) %if factors in file
            numFiles=numFacs;
        else
            numFiles=1;
        end;
        
        for outFile=1:numFiles
            
            if numFiles > 1
                thisOutFileName=[pathstr filesep fileName 'fac' sprintf('%03d',outFile) '.egis'];
                
            else
                thisOutFileName=[pathstr filesep fileName '.egis'];
            end;
            
            if exist(thisOutFileName,'file')
                err=1;
                msg{1}=[thisOutFileName ' already exists!'];
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            outFID = fopen(thisOutFileName, 'w','ieee-be'); %NetStation assumes EGIS files are always big-endian
            if (outFID == -1)
                error('Error creating output file!');
            end;
            
            theData=zeros(numChans,numPoints,numTrials);
            trialCounter=zeros(numCells);
            for theTrial = 1:length(cellIDs)
                newCell=cellIDs(theTrial);
                trialCounter(newCell)=trialCounter(newCell)+1;
                chdr(newCell,6+(trialCounter(newCell)-1)*numSubSpecs:5+trialCounter(newCell)*numSubSpecs)=[EPdata.trialSpecs{theTrial,:}];
                theData(:,:,sum(cellNums(1:newCell-1))+trialCounter(newCell))=ep_expandFacs(EPdata,[],[],theTrial,[],outFile,[]); %rearrange trials so that they are grouped by cell
            end;
            
            theData=reshape(theData,numChans,[]);
            
            %Write out output data
            [status]=wt_ses_hdr_v2(outFID,fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext);	%Write header
            
            %NetStation uses 5bins/microvolt scaling for EGIS session files.
            eval('fwrite(outFID, int16(round(theData*5)), ''int16'',''ieee-be'');');
            
            ST = fclose(outFID);
            if (ST == -1)
                err=1;
                msg{1}='Error closing output file!';
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            try
                if isunix && ismac
                    [err out]=system(['osascript ''' which('changeType.scpt') ''' ''' thisOutFileName ''' EGIS NETs']);
                    if err == 1
                        disp('Unable to set file type and creator.  If you wish to open this file with Netstation, you will need to set this file information manually.  See tutorial.');
                    end;
                    
                    if ~strcmp('no montage',prefs)
                        if isempty(EPdata.montage)
                            EPdata.montage=ep_askForMontage;
                        end;
                        [err out]=ep_addSLAY(thisOutFileName,EPdata.montage);
                        if err == 1
                            disp('Unable to insert montage information into the output EGIS file because:');
                            disp(out);
                        end;
                    end;
                end
            catch
            end
        end
        
        ep_writeEventText(EPdata, [pathstr filesep fileName eventSuffix]);
        ep_writeSubjectText(EPdata, [pathstr filesep fileName subjectSpecSuffix])
        
    case 'egi_egia'
        
        if strcmp(type,'single_trial')
            err=1;
            msg{1}='File format is inconsistent with the file type.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if ~isempty(EPdata.freqNames)
            err=1;
            msg{1}='EGIS format cannot represent frequency data.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        [pathstr, fileName, ext] = fileparts(outFileName);
        if ismember('GAV',EPdata.subTypes) && ismember('AVG',EPdata.subTypes)
            adds2=setdiff(adds,'AVG');
            [err]=ep_writeData(EPdata,[pathstr filesep fileName '_GAV' ext],eventSuffix,subjectSpecSuffix,format,adds2); %save grand average data as a separate EGIS AVE file
            adds=setdiff(adds,'GAV');
            [EPdata]=ep_stripAdds(EPdata,adds); %remove the grand average data which has now been saved.
        end;
        
        if isempty(EPdata.data) %return if no more data to write, as in grand average data but no average data.
            return
        end;
        
        numChans= length(EPdata.chanNames);
        numPoints = length(EPdata.timeNames);
        numSubs = length(EPdata.subNames);
        numFacs= length(EPdata.facNames);
        numTrials = length(EPdata.trialNames);
        
        if ~isempty(EPdata.facNames) %if factors in file
            numFiles=numSubs; %if just GAV data then numSubs will just equal one
        else
            numFiles=1;
        end;
        
        for outFile=1:numFiles %if a factor file then spit out one per subject
            
            if numFiles > 1
                thisOutFileName=[pathstr filesep fileName 'sub' sprintf('%03d',outFile) '.egis'];
            else
                thisOutFileName=[pathstr filesep fileName '.egis'];
            end;
            
            if exist(thisOutFileName,'file')
                err=1;
                msg{1}=[thisOutFileName ' already exists!'];
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            outFID = fopen(thisOutFileName, 'w','ieee-be'); %NetStation assumes EGIS files are always big-endian
            if (outFID == -1)
                err=1;
                msg{1}='Error creating output file!';
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            ave_hdr_offsets_v;
            
            SampleRate = EPdata.Fs;
            
            expname=EPdata.ename;
            
            if size(expname,1) > size(expname,2)
                expname=expname'; %make sure the experiment name is horizontal
            end;
            
            [cellNames m cellIDs]=unique(EPdata.cellNames,'stable'); %get number of trials in each cell
            numCells=length(cellNames);
            cellNums=hist(cellIDs,unique(cellIDs));
            
            fhdr=[];
            chdr=[];
            fhdr(NChan) = numChans;
            fhdr(NCells) = numCells;
            fhdr(HdrVer) = -1;
            fhdr(BaseDur) = EPdata.baseline*(1000/SampleRate);
            fhdr(ScaleBins) = 500;
            fhdr(ScaleCal) = 50;
            fhdr(LastDone) = numSubs;
            
            fcom=[];
            ftext=[];
            for theSubjectSpec = 1:length(EPdata.subjectSpecNames)
                switch EPdata.subjectSpecNames{theSubjectSpec}
                    case 'Comment';
                        fcom=EPdata.subjectSpecs{1,theSubjectSpec};
                    case 'Text';
                        ftext=EPdata.subjectSpecs{1,theSubjectSpec};
                    otherwise
                end;
            end;
            
            fhdr(LComment)=length(fcom);
            fhdr(LText)=length(ftext);
            ename = [expname char(0) blanks(80)];
            ename = ename(1:80);
            
            if isempty(EPdata.facNames)
                numOutSubs=numSubs;
            else
                numOutSubs=numFacs; %if spitting out subject files, the "subjects" in each such file are actually factors.
            end;
            
            %Construct cell headers
            chdr = zeros(numCells,(6+14*numOutSubs));
            for i = 1:numCells
                chdr(i,CellID) = 90 + i;
                chdr(i,NObs) = numOutSubs;
                chdr(i,NPoints) = numPoints;
                chdr(i,SampRate) = SampleRate;
                chdr(i,LSpec) = 14*2;
            end;
            
            %set up the data
            theData=zeros(numChans,numPoints*numCells*numOutSubs);
            for theCell = 1:numCells
                newCell=cellIDs(theCell);
                
                for theOutSub = 1:numOutSubs
                    if isempty(EPdata.facNames) %not factor data
                        theData(:,(newCell-1)*numOutSubs*numPoints+(theOutSub-1)*numPoints+1:(newCell-1)*numOutSubs*numPoints+theOutSub*numPoints)=EPdata.data(:,:,theCell,theOutSub,1); %rearrange data into 2D matrix
                        chdr(newCell,(theOutSub-1)*14+NAvg)=EPdata.avgNum(theOutSub,newCell);
                        for theSubjectSpec = 1:length(EPdata.subjectSpecNames)
                            if ~isempty(EPdata.subjectSpecs{theOutSub,theSubjectSpec})
                                switch EPdata.subjectSpecNames{theSubjectSpec}
                                    case 'RunDateMo'
                                        chdr(newCell,(theOutSub-1)*14+RunDateMo_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'RunDateDay';
                                        chdr(newCell,(theOutSub-1)*14+RunDateDay_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'RunDateYr';
                                        chdr(newCell,(theOutSub-1)*14+RunDateYr_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'RunTimeHr';
                                        chdr(newCell,(theOutSub-1)*14+RunTimeHr_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'RunTimeMin';
                                        chdr(newCell,(theOutSub-1)*14+RunTimeMin_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'RunTimeSec';
                                        chdr(newCell,(theOutSub-1)*14+RunTimeSec_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'SubjID';
                                        chdr(newCell,(theOutSub-1)*14+SubjID_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'Handed';
                                        chdr(newCell,(theOutSub-1)*14+Handed_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'Sex';
                                        chdr(newCell,(theOutSub-1)*14+Sex_Ave)=find(strcmp(EPdata.subjectSpecs{theOutSub,theSubjectSpec}',cellstr(unique(char(EPdata.subjectSpecs{:,theSubjectSpec})))));
                                    case 'Age';
                                        chdr(newCell,(theOutSub-1)*14+Age_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    case 'ExperID';
                                        chdr(newCell,(theOutSub-1)*14+ExperID_Ave)=str2num(EPdata.subjectSpecs{theOutSub,theSubjectSpec});
                                    otherwise
                                end;
                            end;
                        end;
                    else %factor data
                        theData(:,(newCell-1)*numOutSubs*numPoints+(theOutSub-1)*numPoints+1:(newCell-1)*numOutSubs*numPoints+theOutSub*numPoints)=ep_expandFacs(EPdata,[],[],theCell,outFile,theOutSub,[]); %rearrange data into 2D matrix
                    end;
                end;
                if sum(chdr(newCell,SubjID_Ave:14:end)) == 0
                    chdr(newCell,SubjID_Ave:14:end)=1:numOutSubs; %if no subject IDs were available, put some in.
                end;
            end;
            
            L = ((chdr(:,LSpec).*chdr(:,NObs))'+90);
            fhdr(LCellHdr:(LCellHdr+(fhdr(NCells)-1))) = L';
            nlcellhdr = fhdr(NCells)*2;
            lhdr = 130+nlcellhdr+ fhdr(NChan) * 4 + sum(L)+fhdr(LComment)+fhdr(LText);
            fhdr(LPad) = (512-rem(lhdr,512));
            lhdr = lhdr+fhdr(LPad);
            fhdr(LHeader) = lhdr;
            fhdr(LData) = sum((chdr(:,NPoints).*chdr(:,NObs))*fhdr(NChan))*2;
            
            if lhdr > 32767
                err=1;
                msg{1}='Header size is too large.  You need to shrink the file or switch to a different file format.';
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            for cell=1:fhdr(NCells)
                theName=[cellNames{cell}(1,:) char(0) blanks(80)];
                cnames(cell,:)=theName(1:80);
            end;
            fhdr(BytOrd) = 16909060; %always write out EGIS files big-endian
            
            %Write out output data
            [status]=wt_PCAave_hdr_v(outFID,fhdr,chdr,ename,cnames,fcom,ftext);	%Write header
            
            %NetStation uses 500bins/microvolt scaling for EGIS average files.
            eval(['fwrite(outFID, int16(round(theData*500)), ''int16'',''ieee-be'');']);
            
            ST = fclose(outFID);
            if (ST == -1)
                err=1;
                msg{1}='Error closing output file!';
                [msg]=ep_errorMsg(msg);
                return
            end;
            
            try
                if isunix && ismac
                    [err out]=system(['osascript ''' which('changeType.scpt') ''' ''' thisOutFileName ''' EGIA NETs']);
                    if err == 1
                        disp('Unable to set file type and creator.  If you wish to open this file with Netstation, you will need to set this file information manually.  See tutorial.');
                    end;
                    
                    if ~strcmp('no montage',prefs)
                        if isempty(EPdata.montage)
                            EPdata.montage=ep_askForMontage;
                        end;
                        [err out]=ep_addSLAY(thisOutFileName,EPdata.montage);
                        if err == 1
                            disp('Unable to insert montage information into the output EGIS file because:');
                            disp(out);
                        end;
                    end;
                end
            catch
            end
        end;
        
        ep_writeEventText(EPdata, [pathstr filesep fileName eventSuffix]);
        ep_writeSubjectText(EPdata, [pathstr filesep fileName subjectSpecSuffix])
        
    case 'egi_sbin'
        
        if ~isempty(EPdata.freqNames)
            err=1;
            msg{1}='Simple binary format cannot represent frequency data.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        switch type
            case {'single_trial','continuous'}
                
                if ~isempty(find(ismember(EPdata.subTypes,{'AVG','GAV'}), 1)) && ismember('RAW',EPdata.subTypes)
                    adds2=setdiff(adds,'RAW');
                    [err]=ep_writeData(EPdata,[outFileName '_AVG'],eventSuffix,subjectSpecSuffix,format,adds2); %save average data as a separate file
                    adds=setdiff(adds,{'AVG','GAV'});
                    [EPdata]=ep_stripAdds(EPdata,adds); %remove the average data which has now been saved.
                end;
                
                if isempty(EPdata.data) %return if the single-trial data had been deleted.
                    return
                end;
                
                numChans= length(EPdata.chanNames);
                numPoints = length(EPdata.timeNames);
                numSubs = length(EPdata.subNames);
                numFacs= length(EPdata.facNames);
                
                
                SampleRate = EPdata.Fs;
                
                [cellNames m cellIDs]=unique(EPdata.cellNames,'stable'); %get number of trials in each cell
                numCells=length(cellNames);
                cellNums=hist(cellIDs,unique(cellIDs));
                
                if strcmp(type,'single_trial')
                    CateNames=[];
                    for cell = 1:numCells
                        CateNames{cell}=cellNames{cell};
                        CatLengths(cell)=length(CateNames{cell});
                    end;
                    header_array(1) = 5; %version is segmented single-precision
                    numWaves = length(EPdata.trialNames);
                else %Continuous file
                    CateNames=[];
                    CatLengths=[];
                    header_array(1) = 4; %version is unsegmented single-precision
                    numWaves = 1;
                end;
                
                header_array(9) = SampleRate;
                header_array(10) = numChans;
                header_array(11) = 1; %Gain
                header_array(12) = 0; %bits
                header_array(13) = 0; %range
                header_array(14) = numCells; %number of categories
                header_array(15) = numWaves;
                header_array(16) = numPoints; % making assumption that number of samples is same for all cells
                
                for theSubjectSpec = 1:length(EPdata.subjectSpecNames)
                    switch EPdata.subjectSpecNames{theSubjectSpec}
                        case 'RunDateYr'
                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec};
                            if ~isempty(theSpec)
                                header_array(2)=str2num(theSpec);
                            end;
                        case 'RunDateMo';
                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec};
                            if ~isempty(theSpec)
                                header_array(3)=str2num(theSpec);
                            end;
                        case 'RunDateDay';
                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec};
                            if ~isempty(theSpec)
                                header_array(4)=str2num(theSpec);
                            end;
                        case 'RunTimeHr';
                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec};
                            if ~isempty(theSpec)
                                header_array(5)=str2num(theSpec);
                            end;
                        case 'RunTimeMin';
                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec};
                            if ~isempty(theSpec)
                                header_array(6)=str2num(theSpec);
                            end;
                        case 'RunTimeSec';
                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec};
                            if ~isempty(theSpec)
                                header_array(7)=str2num(theSpec);
                            end;
                        case 'RunTimeMsec';
                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec};
                            if ~isempty(theSpec)
                                header_array(8)=str2num(theSpec);
                            end;
                        otherwise
                    end;
                end;
                
                areEvents=0;
                for theSubject=1:size(EPdata.events,1)
                    for theCell=1:size(EPdata.events,2)
                        for theEvent=1:length(EPdata.events{theSubject,theCell})
                            if ~isempty(EPdata.events{theSubject,theCell}(theEvent))
                                areEvents=1;
                            end;
                        end;
                    end;
                end;
                
                if areEvents %if there are events...
                    eventCounter=1;
                    for theSubject=1:size(EPdata.events,1)
                        for theCell=1:size(EPdata.events,2)
                            for theEvent=1:length(EPdata.events{theSubject,theCell})
                                eventValues{eventCounter}=EPdata.events{theSubject,theCell}(theEvent).value;
                                eventContents{eventCounter}=EPdata.events{theSubject,theCell}(theEvent);
                                eventContents{eventCounter}.subject=theSubject;
                                eventContents{eventCounter}.cell=theCell;
                                eventCounter=eventCounter+1;
                            end;
                        end;
                    end;
                    
                    uniqueEvents=unique(cellfun(@num2str,eventValues(find(~cellfun(@isempty,eventValues))),'UniformOutput',false'));
                    for theCode=1:length(uniqueEvents)
                        temp=[uniqueEvents{theCode} blanks(4)];
                        EventCodes(theCode,1:4)=temp(1:4); %events have four-letter labels
                    end;
                    eventData=zeros(length(uniqueEvents),numWaves*numPoints);
                    for event=1:length(eventValues)
                        duration=eventContents{event}.duration;
                        epochSample=ceil(mod(eventContents{event}.sample,numPoints));
                        if epochSample==0
                            epochSample=1;
                        end;
                        if isempty(duration)
                            duration=0;
                        end;
                        theCell=eventContents{event}.cell;
                        theSubject=eventContents{event}.subject;
                        theWave=(theSubject-1)*numCells+theCell;
                        startSample=(theWave-1)*numPoints;
                        duration=ceil(eventContents{event}.duration);
                        if duration ==0
                            duration =1;
                        end;
                        eventData(find(strcmp(uniqueEvents,eventContents{event}.value)),startSample+epochSample:startSample+epochSample+duration-1)=1;
                    end;
                else
                    EventCodes=[];
                    eventData=[];
                    eventContents=[];
                end;
                
                header_array(17) = size(EventCodes,1); %Number of events per segment.
                
                segHdr=zeros(numWaves,2);
                
                if strcmp(type,'single_trial')
                    for theTrial = 1:length(cellIDs) %data are already grouped by cell by readData
                        segHdr(theTrial,1)=find(strcmp(CateNames,EPdata.cellNames{theTrial})); %Figure out which cell each segment belongs to.
                        segHdr(theTrial,2)=(theTrial-1)*numPoints*(1000/SampleRate);
                    end;
                end;
                
                if ~isempty(EPdata.facNames) %if factors in file
                    numFiles=numFacs;
                else
                    numFiles=1;
                end;
                
                [pathstr, fileName, ext]=fileparts(outFileName);
                fileNameStem=fileName;
                for outFile=1:numFiles
                    
                    if numFiles > 1
                        fileNameStemFac=[fileNameStem 'fac' sprintf('%03d',outFile)];
                    else
                        fileNameStemFac=fileNameStem;
                    end;
                    
                    theFileName=fileNameStemFac;
                    
                    sameName=1;
                    theNumber=0;
                    while sameName
                        sameName=0;
                        if exist([pathstr filesep theFileName '.sbin'],'file')
                            sameName=1;
                        end;
                        if sameName
                            theNumber=theNumber+1;
                            theFileName=[fileNameStemFac '-' num2str(theNumber)];
                        end;
                    end;
                    
                    theData=squeeze(ep_expandFacs(EPdata,[],[],[],[],outFile,[],[])); %drop the subject 4th dimension (only one subject for session files).
                    
                    theData=reshape(theData,numChans,[]);
                    
                    [SUCCESS] = Write_eGLY([pathstr filesep theFileName '.sbin'], header_array, CateNames, CatLengths, EventCodes, segHdr, eventData, theData);
                    
                    if (SUCCESS == 0)
                        err=1;
                        msg{1}='Error closing output file!';
                        [msg]=ep_errorMsg(msg);
                        return
                    end;
                    
                    try
                        if isunix && ismac
                            [err out]=system(['osascript ''' which('changeType.scpt') ''' ''' [thisOutFileName '.sbin'] ''' eGLY NETs']);
                            if err == 1
                                disp('Unable to set file type and creator.  If you wish to open this file with Netstation, you will need to set this file information manually.  See tutorial.');
                            end;
                        end
                    catch
                    end
                end;
                
            case 'average'
                
                if ismember('GAV',EPdata.subTypes) && ismember('AVG',EPdata.subTypes)
                    adds2=setdiff(adds,'AVG');
                    [err]=ep_writeData(EPdata,[outFileName '_GAV'],eventSuffix,subjectSpecSuffix,format,adds2); %save grand average data as a separate file
                    adds=setdiff(adds,'GAV');
                    [EPdata]=ep_stripAdds(EPdata,adds); %remove the grand average data which has now been saved.
                end;
                
                if isempty(EPdata.data) %return if no more data to write, as in grand average data but no average data.
                    return
                end;
                
                numChans= length(EPdata.chanNames);
                numPoints = length(EPdata.timeNames);
                numSubs = length(EPdata.subNames);
                numFacs= length(EPdata.facNames);
                numTrials = length(EPdata.trialNames);
                
                if ~isempty(EPdata.facNames) %if factors in file
                    numFiles=numSubs; %if just GAV data then numSubs will just equal one
                else
                    numFiles=1;
                end;
                
                [pathstr, fileName, ext]=fileparts(outFileName);
                fileNameStem=fileName;
                for outFile=1:numFiles %if a factor file then spit out one file (with full set of factors) per subject
                    
                    if numFiles > 1
                        fileNameSubStem=[fileNameStem 'sub' sprintf('%03d',outFile)];
                    else
                        fileNameSubStem=fileNameStem;
                    end;
                    
                    while sameName
                        sameName=0;
                        if exist([pathstr filesep theFileName '.sbin'],'file')
                            sameName=1;
                        end;
                        if sameName
                            theNumber=theNumber+1;
                            theFileName=[fileNameSubStem '-' num2str(theNumber)];
                        end;
                    end;
                    
                    if ~isempty(EPdata.facNames)
                        numWaves = length(EPdata.cellNames)*length(EPdata.facNames);
                    else
                        numWaves = length(EPdata.cellNames)*length(EPdata.subNames);
                    end;
                    SampleRate = EPdata.Fs;
                    
                    %[A index]=sort(eventValues);
                    %[b cellNums n]=unique(A); %get number of trials in each cell
                    %cellNums=[cellNums(1) diff(cellNums)];
                    
                    [cellNames m cellIDs]=unique(EPdata.cellNames,'stable'); %get number of trials in each cell
                    numCells=length(cellNames);
                    cellNums=hist(cellIDs,unique(cellIDs));
                    
                    CateNames=[];
                    for cell = 1:numCells
                        CateNames{cell}=cellNames{cell};
                        CatLengths(cell)=length(CateNames{cell});
                    end;
                    
                    header_array(1) = 5; %version is segmented single-precision
                    header_array(9) = SampleRate;
                    header_array(10) = numChans;
                    header_array(11) = 1; %Gain
                    header_array(12) = 0; %bits
                    header_array(13) = 0; %range
                    header_array(14) = numCells; %number of categories
                    header_array(15) = numWaves;
                    header_array(16) = numPoints; % making assumption that number of samples is same for all cells
                    
                    areEvents=0;
                    for theSubject=1:size(EPdata.events,1)
                        for theCell=1:size(EPdata.events,2)
                            for theEvent=1:length(EPdata.events{theSubject,theCell})
                                if ~isempty(EPdata.events{theSubject,theCell}(theEvent))
                                    areEvents=1;
                                end;
                            end;
                        end;
                    end;
                    
                    if areEvents %if there are events...
                        eventCounter=1;
                        for theSubject=1:size(EPdata.events,1)
                            for theCell=1:size(EPdata.events,2)
                                for theEvent=1:length(EPdata.events{theSubject,theCell})
                                    eventValues{eventCounter}=EPdata.events{theSubject,theCell}(theEvent).value;
                                    eventContents{eventCounter}=EPdata.events{theSubject,theCell}(theEvent);
                                    eventContents{eventCounter}.subject=theSubject;
                                    eventContents{eventCounter}.cell=theCell;
                                    eventCounter=eventCounter+1;
                                end;
                            end;
                        end;
                        
                        uniqueEvents=unique(cellfun(@num2str,eventValues(find(~cellfun(@isempty,eventValues))),'UniformOutput',false'));
                        for theCode=1:length(uniqueEvents)
                            temp=[uniqueEvents{theCode} blanks(4)];
                            EventCodes(theCode,1:4)=temp(1:4); %events have four-letter labels
                        end;
                        eventData=zeros(length(uniqueEvents),numWaves*numPoints);
                        for event=1:length(eventValues)
                            duration=eventContents{event}.duration;
                            epochSample=mod(eventContents{event}.sample,numPoints);
                            if isempty(duration)
                                duration=0;
                            end;
                            theCell=eventContents{event}.cell;
                            theSubject=eventContents{event}.subject;
                            theWave=(theSubject-1)*numCells+theCell;
                            startSample=(theWave-1)*numPoints;
                            duration=eventContents{event}.duration;
                            if duration ==0
                                duration =1;
                            end;
                            eventData(find(strcmp(uniqueEvents,eventContents{event}.value)),startSample+epochSample:startSample+epochSample+duration-1)=1;
                        end;
                    else
                        EventCodes=[];
                        eventData=[];
                        eventContents=[];
                    end;
                    
                    header_array(17) = size(EventCodes,1); %Number of events per segment.
                    
                    if isempty(EPdata.facNames)
                        numOutSubs=numSubs;
                    else
                        numOutSubs=numFacs; %if spitting out subject files, the "subjects" in each such file are actually factors.
                    end;
                    
                    segHdr=zeros(numWaves,2);
                    
                    %set up the data
                    theData=zeros(numChans,numPoints*numCells*numOutSubs); %grouped by cells (all subjects for first cell, then all for second cell, etc.)
                    for theCell = 1:numCells
                        newCell=cellIDs(theCell);
                        
                        for theOutSub = 1:numOutSubs
                            if isempty(EPdata.facNames) %not factor data
                                theData(:,(newCell-1)*numOutSubs*numPoints+(theOutSub-1)*numPoints+1:(newCell-1)*numOutSubs*numPoints+theOutSub*numPoints)=EPdata.data(:,:,theCell,theOutSub,1); %rearrange data into 2D matrix
                                for theSubjectSpec = 1:length(EPdata.subjectSpecNames)
                                    switch EPdata.subjectSpecNames{theSubjectSpec}
                                        case 'RunDateYr'
                                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec}; %simple binary format can accommodate only one subject's specs per file so just use the first one.
                                            if ~isempty(theSpec)
                                                header_array(2)=str2num(theSpec);
                                            end;
                                        case 'RunDateMo';
                                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec}; %simple binary format can accommodate only one subject's specs per file so just use the first one.
                                            if ~isempty(theSpec)
                                                header_array(3)=str2num(theSpec);
                                            end;
                                        case 'RunDateDay';
                                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec}; %simple binary format can accommodate only one subject's specs per file so just use the first one.
                                            if ~isempty(theSpec)
                                                header_array(4)=str2num(theSpec);
                                            end;
                                        case 'RunTimeHr';
                                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec}; %simple binary format can accommodate only one subject's specs per file so just use the first one.
                                            if ~isempty(theSpec)
                                                header_array(5)=str2num(theSpec);
                                            end;
                                        case 'RunTimeMin';
                                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec}; %simple binary format can accommodate only one subject's specs per file so just use the first one.
                                            if ~isempty(theSpec)
                                                header_array(6)=str2num(theSpec);
                                            end;
                                        case 'RunTimeSec';
                                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec}; %simple binary format can accommodate only one subject's specs per file so just use the first one.
                                            if ~isempty(theSpec)
                                                header_array(7)=str2num(theSpec);
                                            end;
                                        case 'RunTimeMsec';
                                            theSpec=EPdata.subjectSpecs{1,theSubjectSpec}; %simple binary format can accommodate only one subject's specs per file so just use the first one.
                                            if ~isempty(theSpec)
                                                header_array(8)=str2num(theSpec);
                                            end;
                                        otherwise
                                    end;
                                end;
                            else %factor data
                                theData(:,(newCell-1)*numOutSubs*numPoints+(theOutSub-1)*numPoints+1:(newCell-1)*numOutSubs*numPoints+theOutSub*numPoints)=ep_expandFacs(EPdata,[],[],theCell,outFile,theOutSub,[]); %rearrange data into 2D matrix
                            end;
                        end;
                    end;
                    
                    for theWave=1:numWaves
                        segHdr(theWave,1)=floor((theWave-1)/numOutSubs)+1;    %cell
                        segHdr(theWave,2)=(theWave-1)*numPoints*(1000/SampleRate);    %time stamp
                    end;
                    
                    [SUCCESS] = Write_eGLY([pathstr filesep theFileName '.sbin'], header_array, CateNames, CatLengths, EventCodes, segHdr, eventData, theData);
                    
                    if (SUCCESS == 0)
                        err=1;
                        msg{1}='Error closing output file!';
                        [msg]=ep_errorMsg(msg);
                        return
                    end;
                    
                    try
                        if isunix && ismac
                            [err out]=system(['osascript ''' which('changeType.scpt') ''' ''' thisOutFileName ''' eGLY NETs']);
                            if err == 1
                                disp('Unable to set file type and creator.  If you wish to open this file with Netstation, you will need to set this file information manually.  See tutorial.');
                            end;
                        end
                    catch
                    end
                end;
                
            otherwise
                err=1;
                msg{1}='Not a supported data type for simple binary files.';
                [msg]=ep_errorMsg(msg);
                return
        end;
        
        ep_writeEventText(EPdata, [outFileName eventSuffix]);
        ep_writeSubjectText(EPdata, [outFileName subjectSpecSuffix])
        
    case 'text'
        numChans= length(EPdata.chanNames);
        numPoints = length(EPdata.timeNames);
        numSubs = length(EPdata.subNames);
        numFacs= length(EPdata.facNames);
        numCells = length(EPdata.cellNames);
        numFreqs = length(EPdata.freqNames);
        
        maxCellName=0;
        for i=1:numCells
            maxCellName=max(maxCellName,length(EPdata.cellNames{i}));
        end;
        
        if length(unique(EPdata.subNames)) == length(EPdata.subNames)
            theSubNames=EPdata.subNames;
        else
            for i=1:numSubs
                theSubNames{i}=sprintf(['S%0' num2str(floor(log10(numSubs))+1) 'd'],i);
            end;
        end;
        
        maxSubName=0;
        for i=1:numSubs
            maxSubName=max(maxSubName,length(theSubNames{i}));
        end;
        
        if ~isempty(EPdata.facNames)
            maxFacName=0;
            for i=1:numFacs
                maxFacName=max(maxFacName,length(EPdata.facNames{i}));
            end;
        end;
        
        if ~isempty(EPdata.freqNames)
            maxFreqName=0;
            for i=1:numFreqs
                maxFreqName=max(maxFreqName,length(num2str(floor(10*EPdata.freqNames(i))))); %number of integer digits plus first decimal
            end;
            maxFreqName=maxFreqName+1; %plus the period
        end;
        
        
        if ~isempty(EPdata.trialNames)
            maxTrialDigits=ceil(log10(max([EPdata.trialNames])));
        end;
        
        if ~isempty(EPdata.freqNames) && isempty(EPdata.timeNames) & isempty(EPdata.relNames) %if FFT data, output as freq x chan rather than time x chan
            numFreqs=1;
            isFFT=1;
        else
            isFFT=0;
        end;
        
        if ~isreal(EPdata.data)
            isComplex=1;
        else
            isComplex=0;
        end;
        
        for theCell=1:numCells
            for theSub=1:numSubs
                for theFactor=1:max(1,numFacs)
                    for theFreq=1:max(1,numFreqs)
                        if isFFT
                            theData=squeeze(ep_expandFacs(EPdata,[],[],theCell,theSub,theFactor,[]))'; %if FFT data, output as freq x chan rather than time x chan
                        else
                            theData=squeeze(ep_expandFacs(EPdata,[],[],theCell,theSub,theFactor,theFreq))';
                        end;
                        
                        if isComplex
                            theDataImag=imag(theData);
                            theData=real(theData);
                        end;
                        
                        facName=[];
                        if ~isempty(EPdata.facNames)
                            facName = ['_' repmat('_',1,maxFacName-length(EPdata.facNames{theFactor})) EPdata.facNames{theFactor}];
                        end;
                        freqName=[];
                        if ~isempty(EPdata.freqNames) && ~isFFT
                            freqName = ['_' sprintf(['%0' num2str(maxFreqName) '.1f'],EPdata.freqNames(theFreq))];
                        end;
                        subName = ['_' repmat('_',1,maxSubName-length(theSubNames{theSub})) theSubNames{theSub}];
                        if ~isempty(EPdata.trialNames)
                            cellName = ['_' repmat('_',1,maxCellName-length(EPdata.cellNames{theCell})) EPdata.cellNames{theCell} '_' sprintf(['%0' num2str(maxTrialDigits) 'd'],EPdata.trialNames(theCell))];
                        else
                            cellName = ['_' repmat('_',1,maxCellName-length(EPdata.cellNames{theCell})) EPdata.cellNames{theCell}];
                        end;
                        
                        if isComplex
                            theName=[outFileName subName cellName facName freqName '_Real.txt'];
                        else
                            theName=[outFileName subName cellName facName freqName '.txt'];
                        end;
                        
                        sameName=1;
                        theNumber=0;
                        [pathstr, name, ext] = fileparts(theName);
                        fileName=name;
                        while sameName
                            sameName=0;
                            if exist([pathstr filesep fileName ext],'file')
                                sameName=1;
                            end;
                            if sameName
                                theNumber=theNumber+1;
                                fileName=[name '-' num2str(theNumber)];
                            end;
                        end;
                        if isComplex
                            dlmwrite([pathstr filesep fileName ext],theData, sprintf('\t'));
                            dlmwrite([pathstr filesep strrep(fileName,'_Real','_Imag') ext],theDataImag, sprintf('\t'));
                        else
                            dlmwrite([pathstr filesep fileName ext],theData, sprintf('\t'));
                        end;
                    end;
                end;
            end;
        end;
        
        ep_writeEventText(EPdata, [pathstr filesep fileName eventSuffix]);
        ep_writeSubjectText(EPdata, [outFileName subjectSpecSuffix])
        
    case 'eeglab_set'
        if ~isempty(EPdata.facNames)
            err=1;
            msg{1}='Error: saving factor files in .set format not currently supported.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if ~isempty(EPdata.freqNames)
            err=1;
            msg{1}='.set format cannot represent frequency data.  EEGlab philosophy is to store the original data and then to generate transforms from it as needed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        [ALLEEG]=ep_ep2alleeg(EPdata);
        sameName=1;
        theNumber=0;
        [pathstr, name, ext] = fileparts(outFileName);
        while sameName
            sameName=0;
            if exist([pathstr filesep fileName '.study'],'file')
                sameName=1;
            end;
            if sameName
                theNumber=theNumber+1;
                fileName=[name '-' num2str(theNumber)];
            end;
        end;
        for dataset=1:length(ALLEEG)
            ALLEEG(dataset).filename=[fileName '-' ALLEEG(dataset).subject '-' ALLEEG(dataset).condition '.set'];
            ALLEEG(dataset).filepath=pathstr;
            pop_saveset(ALLEEG(dataset),'filename',ALLEEG(dataset).filename,'filepath',ALLEEG(dataset).filepath);
        end;
        
        STUDY=[];
        STUDY.name=EPdata.dataName;
        STUDY.task='';
        STUDY.notes='';
        STUDY.filename=[fileName '.study.'];
        STUDY.cluster=[];
        STUDY.history=[];
        for dataset=1:length(ALLEEG)
            STUDY.datasetinfo(dataset).filepath=ALLEEG(dataset).filepath;
            STUDY.datasetinfo(dataset).filename=ALLEEG(dataset).filename;
            STUDY.datasetinfo(dataset).subject=ALLEEG(dataset).subject;
            STUDY.datasetinfo(dataset).session=[];
            STUDY.datasetinfo(dataset).condition=ALLEEG(dataset).condition;
            STUDY.datasetinfo(dataset).group='';
            STUDY.datasetinfo(dataset).index=dataset;
            STUDY.datasetinfo(dataset).comps=[];
        end;
        for chan=1:length(ALLEEG(1).chanlocs)
            STUDY.changrp(chan).name=ALLEEG(1).chanlocs(chan).labels;
            STUDY.changrp(chan).channels{1}=ALLEEG(1).chanlocs(chan).labels;
            for i=1:length(ALLEEG)
                STUDY.changrp(chan).setinds{i,1}=i;
                STUDY.changrp(chan).allinds{i,1}=1;
            end;
            STUDY.changrp(chan).centroid=[];
        end;
        STUDY.filepath=pathstr;
        theSubjects={ALLEEG.subject};
        if ~isempty(theSubjects)
            STUDY.subject=unique(theSubjects);
        else
            STUDY.subject=[];
        end;
        theGroups={ALLEEG.group};
        if ~isempty(theGroups)
            STUDY.group=unique(theGroups);
        else
            STUDY.group=[];
        end;
        theSessions=[ALLEEG.session];
        if ~isempty(theSessions)
            STUDY.session=unique(theSessions);
        else
            STUDY.session=[];
        end;
        theConditions={ALLEEG.condition};
        if ~isempty(theConditions)
            STUDY.condition=unique(theConditions);
        else
            STUDY.condition=[];
        end;
        STUDY.setind=[1:length(ALLEEG)]';
        STUDY.etc=[];
        STUDY.preclust=[];
        STUDY.changrpstatus='all channels present in all datasets';
        STUDY.saved='no';
        
        save('-mat', [pathstr filesep fileName '.study'], 'STUDY');
        
    case 'eeglab_erp'
        if ~strcmp(type,'average')
            err=1;
            msg{1}=['ERPlab files are only for subject average data.'];
            [msg]=ep_errorMsg(msg);
            return
        end;
        if ~isempty(EPdata.facNames)
            err=1;
            msg{1}='Error: saving factor files in .erp format not currently supported.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        if ~isempty(EPdata.freqNames)
            err=1;
            msg{1}='.erp format cannot represent frequency data.  EEGlab philosophy is to store the original data and then to generate transforms from it as needed.';
            [msg]=ep_errorMsg(msg);
            return
        end;
        
        for subject=1:length(EPdata.subNames)
            [EEGALL]=ep_ep2alleeg(ep_selectData(EPdata,{[],[],[],subject,[],[]}));
            ERP=EEGALL(1);
            ERP.bindata=zeros(size(ERP.data,1),size(ERP.data,2),length(EEGALL));
            for theCell=1:length(EEGALL);
                ERP.bindata(:,:,theCell)=EEGALL(theCell).data;
                ERP.bindescr{theCell}=EEGALL(theCell).condition;
            end;
            ERP.nbin=length(EEGALL);
            ERP.ntrials.accepted=squeeze(EPdata.avgNum(subject,:));
            ERP.ntrials.rejected=squeeze(EPdata.analysis.badTrials(subject,:));
            ERP.ntrials.invalid=zeros(size(ERP.ntrials.accepted));
            ERP.erpname=EPdata.subNames{subject};
            ERP.data=[];
            ERP.condition=[];
            ERP.EVENTLIST = [];
            ERP.binerror = [];
            ERP.splinefile = '';
            ERP.ntrials.arflags = zeros(ERP.nbin,8);
            ERP.isfilt=0;
            ERP.ref=[];
            ERP.workfiles=[];
            ERP.nchan=length(EPdata.chanNames);
            ERP.version='4.0.2.1';
            ERP.pexcluded=round(1000*(sum(ERP.ntrials.rejected)/(sum(ERP.ntrials.accepted)+sum(ERP.ntrials.rejected))))/10;
            
            if ~isempty(which('checkERP'))
                checking = checkERP(ERP);
                if ~checking
                    err=1;
                    msg{1}=['Problem with ERPlab file.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end;
            end;
            sameName=1;
            theNumber=0;
            [pathstr, name, ext] = fileparts(outFileName);
            fileName=[name '-' EPdata.subNames{subject}];
            while sameName
                sameName=0;
                if exist([pathstr filesep fileName '.erp'],'file')
                    sameName=1;
                end;
                if sameName
                    theNumber=theNumber+1;
                    fileName=[name '-' EPdata.subNames{subject} '-' num2str(theNumber)];
                end;
            end;
            save('-mat', [pathstr filesep fileName '.erp'], 'ERP');
        end;
        
    case 'neuromag_fif'
        
        %For now, data will be saved as if it is raw data.  Writing evoked data will require additional code writing.
        
        numChans= length(EPdata.chanNames);
        numPoints = length(EPdata.timeNames);
        numSubs = length(EPdata.subNames);
        numFacs= length(EPdata.facNames);
        numCells = length(EPdata.cellNames);
        numFreqs = length(EPdata.freqNames);
        
        if any(strcmp('REGchan', adds)) && any(strcmp('REG',EPdata.chanTypes))
            disp('FIFF files are not really set up for regional channels so not including them.');
            adds=setdiff(adds,'REGchan');
            [EPdata]=ep_stripAdds(EPdata,adds);
        end;
        
        maxCellName=0;
        for i=1:numCells
            maxCellName=max(maxCellName,length(EPdata.cellNames{i}));
        end;
        
        if length(unique(EPdata.subNames)) == length(EPdata.subNames)
            theSubNames=EPdata.subNames;
        else
            for i=1:numSubs
                theSubNames{i}=sprintf(['S%0' num2str(floor(log10(numSubs))+1) 'd'],i);
            end;
        end;
        
        maxSubName=0;
        for i=1:numSubs
            maxSubName=max(maxSubName,length(theSubNames{i}));
        end;
        
        if ~isempty(EPdata.facNames)
            maxFacName=0;
            for i=1:numFacs
                maxFacName=max(maxFacName,length(EPdata.facNames{i}));
            end;
        end;
        
        if ~isempty(EPdata.freqNames)
            maxFreqName=0;
            for i=1:numFreqs
                maxFreqName=max(maxFreqName,length(num2str(floor(10*EPdata.freqNames(i)))))+1; %number of integer digits plus period plus first decimal
            end;
        end;
        
        if ~isempty(EPdata.trialNames)
            maxTrialDigits=ceil(log10(max([EPdata.trialNames])));
        end;
        
        if ~isempty(EPdata.freqNames) && isempty(EPdata.timeNames) %if FFT data, output as freq x chan rather than time x chan
            numFreqs=1;
            isFFT=1;
        else
            isFFT=0;
        end;
        
        [pathstr, fileName, ext]=fileparts(outFileName);
        for theSub=1:numSubs
            for theFactor=1:max(1,numFacs)
                for theFreq=1:max(1,numFreqs)
                    %FIFF has a special data type for FFT data but I don't know enough yet about it to implement it.
                    if isFFT
                        theData=squeeze(ep_expandFacs(EPdata,[],[],[],theSub,theFactor,[])); %if FFT data, output as freq x chan rather than time x chan
                    else
                        theData=squeeze(ep_expandFacs(EPdata,[],[],[],theSub,theFactor,theFreq));
                    end;
                    facName=[];
                    if ~isempty(EPdata.facNames)
                        facName = ['_' repmat('_',1,maxFacName-length(EPdata.facNames{theFactor})) EPdata.facNames{theFactor}];
                    end;
                    freqName=[];
                    if ~isempty(EPdata.freqNames) && ~isFFT
                        freqName = ['_' sprintf(['%0' num2str(maxFreqName) '.1f'],EPdata.freqNames(theFreq))];
                    end;
                    if numSubs > 1
                        subName = ['_' repmat('_',1,maxSubName-length(theSubNames{theSub})) theSubNames{theSub}];
                    else
                        subName='';
                    end;
                    theName=[pathstr filesep fileName subName facName freqName '.fif'];
                    
                    sameName=1;
                    theNumber=0;
                    [pathstr, name, ext] = fileparts(theName);
                    fileName=name;
                    while sameName
                        sameName=0;
                        if exist([pathstr filesep fileName ext],'file')
                            sameName=1;
                        end;
                        if sameName
                            theNumber=theNumber+1;
                            fileName=[name '-' num2str(theNumber)];
                        end;
                    end;
                    
                    %fieldtrip2fiff isn't meeting my needs.  In particular, FieldTrip doesn't seem to be set up to
                    %handle fiducial points that are not data channels and it's critical that they be included in the
                    %file for the MNE software.  Nonetheless, fiffdata is set up to echo the output from FieldTrip in
                    %case it might be helpful later on.
                    
                    fiffdata.file_id.version = NaN;
                    fiffdata.file_id.machid  = [NaN;NaN];
                    fiffdata.file_id.secs    = NaN;
                    fiffdata.file_id.usecs   = NaN;
                    
                    fiffdata.meas_id.version = NaN;
                    fiffdata.meas_id.machid  = [NaN;NaN];
                    fiffdata.meas_id.secs    = NaN;
                    fiffdata.meas_id.usecs   = NaN;
                    
                    fiffdata.meas_date       = [NaN;NaN];
                    fiffdata.nchan    = numChans;
                    fiffdata.sfreq = EPdata.Fs;
                    fiffdata.highpass = NaN;
                    fiffdata.lowpass  = NaN;
                    
                    fiffdata.chs=[];
                    fiffdata.ch_names = EPdata.chanNames;
                    
                    % according to FieldTrip, no strictly necessary, but the inverse functions in MNE works better if
                    % this matrix is present
                    fiffdata.dev_head_t.from = 1;
                    fiffdata.dev_head_t.to = 4;
                    fiffdata.dev_head_t.trans = eye(4);
                    
                    fiffdata.ctf_head_t = [];
                    fiffdata.dev_ctf_t = [];
                    
                    fiffdata.dig=[];
                    
                    fiffdata.bads={};
                    for iChan=1:length(EPdata.chanNames)
                        if strcmp(EPdata.dataType,'average')
                            if any(isnan(EPdata.analysis.badChans(theSub,:,iChan)))
                                fiffdata.bads{1,end+1}=EPdata.chanNames{iChan};
                            end;
                        else
                            if any(EPdata.analysis.badChans(theSub,:,iChan) < 0)
                                fiffdata.bads{1,end+1}=EPdata.chanNames{iChan};
                            end;
                        end;
                    end;
                    
                    fiffdata.projs = struct('kind', {}, 'active', {}, 'desc', {}, 'data', {});
                    fiffdata.comps = struct('ctfkind', {}, 'kind', {}, 'save_calibrated', {}, ...
                        'rowcals', {}, 'colcals', {}, 'data', {});
                    
                    fiffdata.acq_pars = []; % needed by raw
                    fiffdata.acq_stim = []; % needed by raw
                    
                    fiffdata.evoked=[];
                    fiffdata.info=[];
                    
                    fiffdata.vartriallength = 0;
                    fiffdata.isaverage = 0;
                    fiffdata.isepoched = 0;
                    fiffdata.iscontinuous = 0;
                    
                    identCounter=zeros(2,1);
                    for iChan=1:length(EPdata.chanNames)
                        fiffdata.chs(iChan).scanno=iChan;
                        if any(strcmp(EPdata.chanTypes{iChan},{'EEG','REG'}))
                            fiffdata.chs(iChan).kind=2;
                            fiffdata.dig(end+1).kind=3;
                            identCounter(1)=identCounter(1)+1;
                            if ~isempty(EPdata.eloc)
                                fiffdata.chs(iChan).logno=identCounter(1);
                                fiffdata.dig(end).ident=identCounter(1);
                            end;
                        elseif any(strcmp(EPdata.chanTypes{iChan},{'MGM','MGA','MGP'}))
                            fiffdata.chs(iChan).kind=1;
                            fiffdata.dig(end+1).kind=4;
                            identCounter(2)=identCounter(2)+1;
                            if ~isempty(EPdata.eloc)
                                fiffdata.chs(iChan).logno=identCounter(2);
                                fiffdata.dig(end).ident=identCounter(2);
                            end;
                        end
                        fiffdata.chs(iChan).range=1;
                        fiffdata.chs(iChan).cal=1;
                        fiffdata.chs(iChan).coil_type=1;
                        
                        theChan=find(strcmp(EPdata.chanNames{iChan},{EPdata.eloc(:).labels}));
                        if ~isempty(theChan) && any(strcmp(EPdata.chanTypes{iChan},{'EEG','REG'})) && ~isempty(EPdata.eloc)
                            fiffdata.chs(iChan).loc=[EPdata.eloc(theChan).Y; EPdata.eloc(theChan).X; EPdata.eloc(theChan).Z; 0; 0; 0; 0; 1; 0; 0; 0; 1]; %FIFF is Y+ is nose whereas EEGLab (and hence EP) seems to be X+ is the nose
                            fiffdata.chs(iChan).coil_trans=[];
                            fiffdata.chs(iChan).eeg_loc=[EPdata.eloc(theChan).Y 0; EPdata.eloc(theChan).X 0; EPdata.eloc(theChan).Z 0]; %FIFF is Y+ is nose whereas EEGLab (and hence EP) seems to be X+ is the nose
                            
                            fiffdata.dig(end).r(1)=EPdata.eloc(theChan).Y/100;
                            fiffdata.dig(end).r(2)=EPdata.eloc(theChan).X/100;
                            fiffdata.dig(end).r(3)=EPdata.eloc(theChan).Z/100;
                            fiffdata.dig(end).coord_frame=4;
                        else
                            fiffdata.chs(iChan).loc=[0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 1];
                            fiffdata.chs(iChan).coil_trans=[];
                            fiffdata.chs(iChan).eeg_loc=[0 0; 0 0; 0 0];
                            
                        end
                        fiffdata.chs(iChan).coord_frame=4;
                        fiffdata.chs(iChan).unit=107; %volt
                        fiffdata.chs(iChan).unit_mul=-6; %the volts are micro: 10^(-6)
                        fiffdata.chs(iChan).ch_name=EPdata.chanNames{iChan};
                    end
                    if ~isempty(EPdata.eloc)
                        if length(EPdata.implicit)==3 %add three fiducials
                            fiffdata.dig(end+1).kind=1;
                            fiffdata.dig(end).ident=1;
                            fiffdata.dig(end).r(1)=EPdata.implicit(1).Y/100;
                            fiffdata.dig(end).r(2)=EPdata.implicit(1).X/100;
                            fiffdata.dig(end).r(3)=EPdata.implicit(1).Z/100;
                            fiffdata.dig(end).coord_frame=4;
                            fiffdata.dig(end+1).kind=1;
                            fiffdata.dig(end).ident=1;
                            fiffdata.dig(end).r(1)=EPdata.implicit(2).Y/100;
                            fiffdata.dig(end).r(2)=EPdata.implicit(2).X/100;
                            fiffdata.dig(end).r(3)=EPdata.implicit(2).Z/100;
                            fiffdata.dig(end).coord_frame=4;
                            fiffdata.dig(end+1).kind=1;
                            fiffdata.dig(end).ident=1;
                            fiffdata.dig(end).r(1)=EPdata.implicit(3).Y/100;
                            fiffdata.dig(end).r(2)=EPdata.implicit(3).X/100;
                            fiffdata.dig(end).r(3)=EPdata.implicit(3).Z/100;
                            fiffdata.dig(end).coord_frame=4;
                        else
                            disp('There needed to be three implicit sites to include the fiducials.');
                        end;
                    end;
                    
                    fiffdata.info.file_id=fiffdata.file_id;
                    fiffdata.info.meas_id=fiffdata.meas_id;
                    fiffdata.info.meas_date=fiffdata.meas_date;
                    fiffdata.info.nchan=fiffdata.nchan;
                    fiffdata.info.sfreq=fiffdata.sfreq;
                    fiffdata.info.highpass=fiffdata.highpass;
                    fiffdata.info.lowpass=fiffdata.lowpass;
                    fiffdata.info.chs=fiffdata.chs;
                    fiffdata.info.ch_names=fiffdata.ch_names;
                    fiffdata.info.dev_head_t=fiffdata.dev_head_t;
                    fiffdata.info.ctf_head_t=fiffdata.ctf_head_t;
                    fiffdata.info.dev_ctf_t=fiffdata.dev_ctf_t;
                    fiffdata.info.dig=fiffdata.dig;
                    fiffdata.info.bads=fiffdata.bads;
                    fiffdata.info.projs=fiffdata.projs;
                    fiffdata.info.comps=fiffdata.comps;
                    fiffdata.info.acq_pars=fiffdata.acq_pars;
                    fiffdata.info.acq_stim=fiffdata.acq_stim;
                    fiffdata.info.filename=[fileName '.fif'];
                    
                    if any(strcmp(EPdata.dataType,{'average','single_trial'}))
                        fiffdata.isepoched = 1;
                        if strcmp(EPdata.dataType,'average')
                            theKind=100;
                            fiffdata.isaverage = 1;
                        else
                            theKind=102;
                        end;
                        
                        for iCell=1:length(EPdata.cellNames)
                            fiffdata.evoked.aspect_kind=theKind;
                            fiffdata.evoked.is_smsh=0;
                            fiffdata.evoked.nave=EPdata.covNum(theSub,iCell); %Use estimated effective sample size rather than true sample size, for purposes of estimating noise levels.
                            fiffdata.evoked.first=-EPdata.baseline;
                            fiffdata.evoked.last=numPoints-EPdata.baseline-1;
                            fiffdata.evoked.comment=EPdata.cellNames{iCell};
                            fiffdata.evoked.times=EPdata.timeNames;
                            fiffdata.evoked.epochs=squeeze(ep_expandFacs(EPdata,[],[],iCell,theSub,theFactor,theFreq));
                        end;
                        
                        
                        fiff_write_evoked([pathstr filesep fileName ext], fiffdata);
                        
                        if ~isempty(EPdata.cov)
                            disp('Adding cov information to a separate fiff file ending in -cov.fif.');
                            theCov.kind=1;
                            theCov.diag=0;
                            theCov.dim=size(EPdata.cov.covMatrix,2);
                            theCov.names=EPdata.chanNames';
                            theCov.data=squeeze(EPdata.cov.covMatrix(theSub,:,:));
                            theCov.projs=struct([]);
                            theCov.bads=EPdata.chanNames(find(isnan(EPdata.cov.covMatrix(1,:,1))))';
                            theCov.nfree=EPdata.cov.Nq;
                            theCov.eig=[];
                            theCov.eigvec=[];
                            mne_write_cov_file([pathstr filesep fileName '-cov' ext],theCov);
                        end;
                    else
                        
                        %todo - events
                        rmfield(fiffdata,'evoked'); %not for continuous files
                        
                        [outfid, cals] = fiff_start_writing_raw([pathstr filesep fileName ext], fiffdata);
                        fiff_write_raw_buffer(outfid, squeeze(EPdata.data), cals);
                        fiff_finish_writing_raw(outfid);
                    end;
                    
                end;
                
            end;
        end;
        
        ep_writeEventText(EPdata, [pathstr filesep fileName eventSuffix]);
        ep_writeSubjectText(EPdata, [pathstr filesep fileName subjectSpecSuffix])
        
    otherwise
        err=1;
        msg{1}=[format ' is not a supported file format for file writing.'];
        [msg]=ep_errorMsg(msg);
        return
end;
