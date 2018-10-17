function abortEP=ep_initEP;

% ep_initEP - abortEP=ep_initEP -
%          Set up EP Toolkit and check for installed software packages.
%
%Outputs:
%  abortEP       : Equals 1 if need to abort EP, else 0.
%
%History
%  by Joseph Dien (1/31/09)
%  jdien07@mac.com
%
%  bugfix 8/3/10 JD
%  Fixed crash when EP Toolkit folder name has been changed.
%
%  modified 9/14/13 JD
%  Initializes mff at the outset so that global variables bug doesn't crash the toolkit down the line.
%
%  bugfix 3/23/14 JD
%  Adds MNE toolbox to the path to support reading and writing FIFF files.
%
%  modified 8/12/14 JD
%  Adds new EP Toolkit directories automatically.
%
%  modified 10/8/14 JD
%  Checks for function name conflicts and aborts EP if any detected.
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

global EP_COPYRIGHT_NOTE EP_VER

abortEP=0;

try
    
    EPpath=which('ep.m');
    [pathstr, name, ext]=fileparts(EPpath);
    [a b]=strtok(pathstr,filesep);
    b=b(2:end);
    EPver=ver(b);
    
    if isempty(EPver.Version)
        disp(' ');
        disp('**************************************************************');
        disp('It looks as though you may have two versions of the EP Toolkit installed.  You need to delete one of them and then restart.');
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
        
    end;
    
    MATLABver=ver('MATLAB');
    [a b]=strtok(MATLABver.Version,'.');
    b=b(2:end);
    
    if (str2num(a) < 9)
        disp(' ');
        disp('**************************************************************');
        disp('Warning: The EP Toolkit is only supported for versions Matlab 9.00 (2016a) onwards.');
        disp('You are welcome to try it out but please do not ask me for help with any problems you may encounter.');
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
    end;
catch
    EPver.Version='';
    MATLABver='unavailable'; %workaround for bug in earlier version of Matlab
end;

if ~exist('EP_VER','var')
    EP_VER='';
end;

if isempty(EP_COPYRIGHT_NOTE) || ~strcmp(EPver.Version,EP_VER)
    disp(' ');
    disp(['ERP PCA Toolbox ' EPver.Version ' Copyright (C) 1999-2018  Joseph Dien']);
    disp('This program comes with ABSOLUTELY NO WARRANTY.');
    disp('This is free software, and you are welcome to redistribute it');
    disp('under certain conditions; see http://www.gnu.org/licenses/ for details.');
    disp(' ');
    disp('Checking installation locations of FieldTrip and EEGlab:');
    eval('which ft_read_header')
    eval('which runica')
    
    EP_COPYRIGHT_NOTE='Y';
    EP_VER=EPver.Version; %ensure that the copyright notice will be printed again if the version is updated during a session.
end;

if ~exist('ft_read_header','file') && exist('read_header','file')
    disp(' ');
    disp('**************************************************************');
    disp('Warning: Fieldtrip needs to be updated a current version.');
    disp('**************************************************************');
    disp(' ');
    warndlg('See command window for warning message.','!! Warning !!')
end;

if ~exist('ft_read_header','file') || ~exist('ft_read_data','file') || ~exist('ft_read_event','file')
    disp(' ');
    disp('**************************************************************');
    disp(['Warning: Fieldtrip has either not been installed or put in the path.  This needs to be done before running EP Toolkit or many things will not work.']);
    disp('Remember to use ft_defaults to add the necessary subdirectories too, not just the main directory.');
    disp(['See tutorial file for details.']);
    disp('**************************************************************');
    disp(' ');
    warndlg('See command window for warning message.','!! Warning !!')
    abortEP=1;  
else
    fieldTripLoc=which('ft_read_header');
    if ~isempty(findstr(fieldTripLoc,['external' filesep 'fieldtrip'])) || ~isempty(findstr(fieldTripLoc,['external' filesep 'fieldtrip-partial'])) || ~isempty(findstr(fieldTripLoc,['external' filesep 'fileio']))
        disp(' ');
        disp('**************************************************************');
        disp('Warning: You appear to be using the version of Fieldtrip that comes with EEGlab and/or in SPM (in its external directory).');
        disp('If you rely on it alone, it will tend to be an obsolete version, with unpatched bugs causing problems.');
        disp('If you have both a current version of FieldTrip installed and the older version that comes with EEGlab and/or SPM, then they tend to conflict and cause problems.');
        disp('It is recommended that you delete EEGlab and SPM''s copies or otherwise disable them when using the EP Toolkit and instead rely on just having the most recent version of FieldTrip installed separately.');
        disp(['See tutorial file for details.']);
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
    end;
    
    p=path;
    if ~isempty(findstr(p,['compat' filesep 'matlablt2010b']))
        disp(' ');
        disp('**************************************************************');
        disp(['Warning: It appears that the subdirectories of fieldtrip were not added to the path using ft_defaults.']);
        disp('This means that there are subdirectories are on the path that should not be, which can destabilize Matlab.');
        disp(['See tutorial file for details.']);
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
    end;

    fieldTripLoc=which('ft_freqanalysis');
    if ~isempty(findstr(fieldTripLoc,['external' filesep 'fieldtrip']))
        disp(' ');
        disp('**************************************************************');
        disp('Warning: You appear to be using the version of ft_freqanalysis that comes with SPM (in its external directory):');
        disp(which('ft_freqanalysis'));
        disp('It is recommended that you delete it or at least remove it from Matlab''s path list.');
        disp(['See tutorial file for details.']);
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
    end;
    theLoc=findstr(which('ft_read_header'),'fieldtrip');
    if str2num(fieldTripLoc(theLoc+10:theLoc+17)) < 20150924
        disp(' ');
        disp('**************************************************************');
        disp('Warning: You appear to be using a version of FieldTrip that is too old to work properly with the EP Toolkit.');
        disp(['It dates back to: ' datestr([str2num(fieldTripLoc(theLoc+10:theLoc+13)) str2num(fieldTripLoc(theLoc+14:theLoc+15)) str2num(fieldTripLoc(theLoc+16:theLoc+17)) 0 0 0])])
        disp('This may cause problems importing files and performing dipole fitting.');
        disp('It is recommended that you delete it and download the current version of FieldTrip.');
        disp('If you just downloaded it, try looking again.  If the file names were ordered by name rather than date and you took the top one, then that was the oldest version.')
        disp(['See tutorial file for details.']);
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
    end;
    cellfunLoc=which('cellfun');
    if ~isempty(findstr(cellfunLoc,['compat' filesep 'R']))
        disp(' ');
        disp('**************************************************************');
        disp('Warning: There is a copy of cellfun.m in the Fieldtrip folder that needs to be deleted immediately or it will destabilize Matlab.');
        disp('It is located in the comp directory, in folders named R13 and/or R14.');
        disp('It is recommended that you delete it immediately.');
        disp(['See tutorial file for details.']);
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
    end;
end;

if ~exist('traditionaldipfit','file')
    disp(' ');
    disp('**************************************************************');
    disp('Warning: The dipfit plugin appears not to be installed in your copy of EEGlab or EEGlab itself has not been installed).');
    disp('This plug-in is needed for dipole source analyses and any attempt to perform them without this plugin will result in a crash.');
    disp(['See tutorial file for details.']);
    disp('**************************************************************');
    disp(' ');
    warndlg('See command window for warning message.','!! Warning !!')
end;

if exist('icadefs','file')
    icadefs;
    if exist(ICABINARY) == 2
        [pathstr, name, ext]=fileparts(ICABINARY);
        if ismac && ~(strcmp(name,'ica_osx') || strcmp(name,'ica_osxG5'))
            disp(['Warning: You have a binary ICA program installed.  The file is called: ' name]);
            disp('It is not one of the two that is appropriate for OS X.');
            disp('See tutorial for instructions on how to deinstall it.');
            warndlg('See command window for warning message.','!! Warning !!')
        end;
        try
            if isunix && ~ismac && ~(strcmp(name,'ica_linux') || strcmp(name,'ica_bsd'))
                disp(['Warning: You have a binary ICA program installed.  The file is called: ' name]);
                disp('It is not one of the two that is appropriate for UNIX.');
                disp('See tutorial for instructions on how to deinstall it.');
                warndlg('See command window for warning message.','!! Warning !!')
            end;
        catch
        end
        if ispc
            disp(['Warning: You have a binary ICA program installed.  The file is called: ' name]);
            disp('There is not currently one that is appropriate for Windows.');
            warndlg('See command window for warning message.','!! Warning !!')
        end;
        if ~isempty(findstr(which('binica'),['external' filesep 'eeglab']))
            disp(' ');
            disp('**************************************************************');
            disp('Warning: You appear to be using the version of binica that comes with Fieldtrip (in its external directory).');
            disp('It is recommended that you delete it.');
            disp(['See tutorial file for details.']);
            disp('**************************************************************');
            disp(' ');
            warndlg('See command window for warning message.','!! Warning !!')
        end;
    else
        if ~isempty(findstr(which('runica'),['external' filesep 'eeglab']))
            disp(' ');
            disp('**************************************************************');
            disp('Warning: You appear to be using the version of runica that comes with Fieldtrip or SPM (in their external directories):');
            disp(which('runica'));
            disp('It is recommended that you delete it.');
            disp(['See tutorial file for details.']);
            disp('**************************************************************');
            disp(' ');
            warndlg('See command window for warning message.','!! Warning !!')
        end;
    end
    p=path;
    if ~isempty(findstr(p,['functions' filesep 'octavefunc']))
        disp(' ');
        disp('**************************************************************');
        disp(['Warning: It appears that the subdirectories of EEGlab were not added to the path by just starting eeglab.']);
        disp(['Instead it was apparently added by using the "add with subfolders" button in the Set Paths matlab window.']);
        disp('This means that there are subdirectories are on the path that should not be, which can destabilize Matlab.');
        disp(['See tutorial file for details.']);
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
    end;
    
    
else
    disp(' ');
    disp('**************************************************************');
    disp(['Warning: ICA is not installed.  You will need to install it and add it to the path,']);
    disp(['otherwise you will not be able to use the Infomax ICA rotation.  See tutorial.']);
    disp('Remember to add its subdirectories too, not just the main directory.');
    disp('**************************************************************');
    disp(' ');
    warndlg('See command window for warning message.','!! Warning !!')
end;

if exist('biosig/NaN','dir')
    disp(' ');
    disp('**************************************************************');
    disp(['Warning: The Biosig directory "NaN" is in the path and may cause problems.']);
    disp(['It is recommended that it be removed from the path.']);
    disp('**************************************************************');
    disp(' ');
    warndlg('See command window for warning message.','!! Warning !!')
end;

if exist('biosig/maybe-missing','dir')
    disp(' ');
    disp('**************************************************************');
    disp(['Warning: The Biosig directory "maybe-missing" is in the path and may cause problems.']);
    disp(['It is recommended that it be removed from the path.']);
    disp('**************************************************************');
    disp(' ');
    warndlg('See command window for warning message.','!! Warning !!')
end;

scrsz = get(0,'ScreenSize');
if max(scrsz) <2
    msg{1}='Error: Matlab is currently unable to determine the size of the monitor.  Please restart Matlab.';
    [msg]=ep_errorMsg(msg);
    return
end;

if scrsz(4) < 500
    disp(' ');
    disp('**************************************************************');
    disp(['The current vertical screen height of ' num2str(scrsz(4)) ' pixels is too low.  It needs to be at least 500.']);
    disp('**************************************************************');
    disp(' ');
    warndlg('See command window for warning message.','!! Warning !!')
end;

disp('Checking to see if there is a newer version of the EP Toolkit available.')
[sourceForge, status]= urlread('http://sourceforge.net/projects/erppcatoolkit/files/erppcatoolkit/');
if status
    fileLine='<tr title="EP_Toolkit ';
    whichLines=findstr(sourceForge,fileLine);
    EPToolkitFiles=[];
    for i=1:length(whichLines)
        a=str2num(sourceForge(whichLines(i)+length(fileLine):whichLines(i)+length(fileLine)+2));
        if ~isempty(a)
            EPToolkitFiles(i)=a/100;
        end;
    end;
    if max(EPToolkitFiles) > str2num(EP_VER)
        disp(' ');
        disp('**************************************************************');
        disp(['More recent version (' num2str(max(EPToolkitFiles)) ') of EP Toolkit now available for download.']);
        disp('**************************************************************');
        disp(' ');
        warndlg(['More recent version (' num2str(max(EPToolkitFiles)) ') of EP Toolkit now available for download.'],'Alert')
    end;
else
    disp('Could not reach internet to check for updates to EP Toolkit.');
end;
disp('Done checking.')

lastwarn('');
if exist('fieldTripLoc','var')
    %set up mff at the outset so the mff bug doesn't crash the toolkit down the line.
    addpath([fileparts(fieldTripLoc) filesep 'external' filesep 'egi_mff']);
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for Matlab bug resulting in global variables being cleared
    globalTemp=cell(0);
    globalList=whos('global');
    varList=whos;
    for i=1:length(globalList)
        eval(['global ' globalList(i).name ';']);
        eval(['globalTemp{end+1}=' globalList(i).name ';']);
    end;
    %%%%%%%%%%%%%%%%%%%%%%
    
    evalc('mff_setup;');
    
    %%%%%%%%%%%%%%%%%%%%%%
    %workaround for Matlab bug resulting in global variables being cleared
    varNames={varList.name};
    for i=1:length(globalList)
        eval(['global ' globalList(i).name ';']);
        eval([globalList(i).name '=globalTemp{i};']);
        if ~any(strcmp(globalList(i).name,varNames)) %was global variable originally out of scope?
            eval(['clear ' globalList(i).name ';']); %clears link to global variable without affecting it
        end;
    end;
    clear globalTemp globalList varNames varList;
    %%%%%%%%%%%%%%%%%%%%%%
    
    %add the mne toolkit to the path to support reading and writing FIFF files.
    addpath([fileparts(fieldTripLoc) filesep 'external' filesep 'mne']);
end;

%update path with EP Toolkit directories in case a new one has been added
p=path;
EPtoolkit=strfind(p,'EP_Toolkit:');
if ~isempty(EPtoolkit)
    p=p(1:EPtoolkit+length('EP_Toolkit:')-2);
    colonPos=strfind(p,':');
    if ~isempty(colonPos)
        p=p(colonPos(end)+1:end);
    end;
    p2=genpath(p);
    addpath(p2);
end;

[msgstr, msgid] = lastwarn;
if strcmp(msgid,'MATLAB:dispatcher:nameConflict')
    disp('***ALERT**** - You have a function on your path that will override one of Matlab''s own functions, causing unpredictable errors.')
    disp('This can happen, for example, if you add FieldTrip to your path by adding everything rather than via ft_defaults.');
    disp('See the tutorial for full instructions on installing fieldtrip.');
    disp('To avoid problems, EP Toolkit will not run until this problem is solved.');
    abortEP=1;   
end;

