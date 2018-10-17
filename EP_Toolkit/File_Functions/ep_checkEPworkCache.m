function EPdataset=ep_checkEPworkCache(EPdataset);
%  EPdataset=ep_checkEPworkCache(EPdataset);
%       Checks integrity of the EPwork cache and regenerates it from scratch if out-of-date or corrupted.
%
%Inputs:
%  EPdataset      : Structured array with the list of files in the work directory
%     .EPwork     : The path of the work directory.
%     .dataName   : The list of dataset names.
%
%Outputs:
%  EPdataset      : Structured array with the list of files in the work directory

%History:
%  by Joseph Dien (11/1/09)
%  jdien07@mac.com
%
%  modified 1/20/10 JD
%  Added implicit, baseline, and Fs to cache contents.
%
%  modified 5/12/10 JD
%  Added ced to cache contents.
%
%  modified 6/15/10 JD
%  Added saved flag to cache contents.
%
%  modified 1/24/12 JD
%  Added frequencies to cache contents.
%
%  modified 1/11/13 JD
%  Added power to cache contents.
%
%  modified 10/13/13 JD
%  Added trial specs and events to cache contents.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 4/23/14 JD
% Added pca and subjectSpecNames and timeNames fields to the cache to speed up the windowing function.
%
% modified 5/25/14 JD
% Added facVecT and FacVecS and facVecF to the cache to support the sampleTest function.
%
% bugfix 9/9/14 JD
% Fixed crash when dataset in working set had older version of event structure.
%
% modified 9/3/15 JD
% Added recTime to the cache to support the segment function.
%
% modified 1/23/16 JD
% Added support for complex spectral data.
%
% modified 1/25/17 JD
% Eliminated support for .power field.
% Eliminated EPver field.
%
% modified 6/13/17 JD
% Added .timeUnits field.
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

cacheBad=0;

if isempty(EPdataset.dataset)
    EPdataset.dataset='';
end;
if ~isempty(EPdataset.dataset)
%     if ~isfield(EPdataset,'EPver')
%         cacheBad=1;
%     end;
    if ~isfield(EPdataset.dataset,'dataName')
        cacheBad=2;
    end;
    if ~isfield(EPdataset.dataset,'chanNames')
        cacheBad=3;
    end;
    if ~isfield(EPdataset.dataset,'timeNames')
        cacheBad=4;
    end;
    if ~isfield(EPdataset.dataset,'cellNames')
        cacheBad=5;
    end;
    if ~isfield(EPdataset.dataset,'trialNames')
        cacheBad=6;
    end;
    if ~isfield(EPdataset.dataset,'subNames')
        cacheBad=7;
    end;
    if ~isfield(EPdataset.dataset,'facNames')
        cacheBad=8;
    end;
    if ~isfield(EPdataset.dataset,'freqNames')
        cacheBad=9;
    end;
    if ~isfield(EPdataset.dataset,'relNames')
        cacheBad=10;
    end;
    if ~isfield(EPdataset.dataset,'eloc')
        cacheBad=11;
    end;
    if ~isfield(EPdataset.dataset,'dataType')
        cacheBad=12;
    end;
    if ~isfield(EPdataset.dataset,'chanTypes')
        cacheBad=13;
    end;
    if ~isfield(EPdataset.dataset,'subTypes')
        cacheBad=14;
    end;
    if ~isfield(EPdataset.dataset,'cellTypes')
        cacheBad=15;
    end;
    if ~isfield(EPdataset.dataset,'facTypes')
        cacheBad=16;
    end;
    if ~isfield(EPdataset.dataset,'plotMVmin')
        cacheBad=17;
    end;
    if ~isfield(EPdataset.dataset,'plotMVmax')
        cacheBad=18;
    end;
    if ~isfield(EPdataset.dataset,'implicit')
        cacheBad=19;
    end;
    if ~isfield(EPdataset.dataset,'baseline')
        cacheBad=20;
    end;
    if ~isfield(EPdataset.dataset,'Fs')
        cacheBad=1;
    end;
    if ~isfield(EPdataset.dataset,'ced')
        cacheBad=21;
    end;
    if ~isfield(EPdataset.dataset,'saved')
        cacheBad=22;
    end;
    if ~isfield(EPdataset.dataset,'trialSpecNames')
        cacheBad=23;
    end;
    if ~isfield(EPdataset.dataset,'trialSpecs')
        cacheBad=24;
    end;
    if ~isfield(EPdataset.dataset,'pca')
        cacheBad=25;
    end;
    if ~isfield(EPdataset.dataset,'subjectSpecNames')
        cacheBad=26;
    end;
    if ~isfield(EPdataset.dataset,'timeNames')
        cacheBad=27;
    end;
    if ~isfield(EPdataset.dataset,'facVecT')
        cacheBad=28;
    end;
    if ~isfield(EPdataset.dataset,'facVecS')
        cacheBad=29;
    end;
    if ~isfield(EPdataset.dataset,'facVecF')
        cacheBad=30;
    end;
    if ~isfield(EPdataset.dataset,'recTime')
        cacheBad=31;
    end;
    if ~isfield(EPdataset.dataset,'plotMVminAbs')
        cacheBad=32;
    end;
    if ~isfield(EPdataset.dataset,'plotMVmaxAbs')
        cacheBad=33;
    end;
    if ~isfield(EPdataset.dataset,'timeUnits')
        cacheBad=34;
    end;
end;

if cacheBad
    disp('EPwork cache out-of-date or corrupted.  Regenerating EPwork cache.');
    
    newCache=[];
    
    for i=1:length(EPdataset.dataset)
        oldEPdata=ep_loadEPdataset(EPdataset,i);
        newDataset=ep_addToEPworkCache(oldEPdata);
        newCache=[newCache newDataset];
    end;
    EPdataset.dataset=newCache;
%     try
%         EPver=ver('EP_Toolkit');
%     catch
%         EPver='unavailable'; %workaround for bug in earlier version of Matlab
%     end;
%     EPdataset.EPver=EPver;
    
    eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
    
    disp('Done.')
end;