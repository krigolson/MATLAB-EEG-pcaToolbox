function newDataset=ep_addToEPworkCache(EPdata);
%  newDataset=addToEPworkCache(EPdata);
%       Generates the EPwork cache entry for the data file.
%
%Inputs:
%  EPdata      : Structured array with the data and accompanying information.  See readData.
%
%Outputs:
%  newDataset      : Structured array with the list of files in the work directory
%     .EPwork     : The path of the work directory.
%     .dataName   : The list of dataset names.
%

%History:
%  by Joseph Dien (11/1/09)
%  jdien07@mac.com
%
%  modified 1/24/10 JD
%  Added implicit, baseline, and Fs to cache contents.
%  For continuous files, calculates plotMVmin and max for each one second epoch.
%
%  modified 5/12/10 JD
%  Added ced to cache contents.
%
%  modified 6/15/10 JD
%  Added saved flag to cache contents.
%
%  bufix 7/24/10 JD
%  Fixed min and max of factors sometimes not being calculated correctly.
%
%  modified 4/25/11 JD
%  Added transforms to cache contents.
%
%  modified 8/22/11 JD
%  Eliminated transforms from cache contents.
%
%  modified 2/7/12 JD
%  Added frequencies to cache contents.
%
%  modified 1/11/13 JD
%  Added power to cache contents.
% 
% modified 10/8/13 JD
% Restricted min and max voltage to EEG channels.
%
%  modified 10/13/13 JD
%  Added trial specs and events to cache contents.
%
%  bufix 12/23/13 JD
%  Fixed crash when regenerating cache and there is a spatial PCA in the
%  active set.
%
% bufix 3/12/14 JD
% Handles decimal sampling rates gracefully.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 4/23/14 JD
% Added pca and subjectSpecNames and timeNames fields to the cache to speed up the windowing function.
%
% bufix 4/28/14 JD
% Fixed not computing frequency min and max correctly for PCA datasets.
%
% modified 4/30/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure, including complex numbers.
%
% modified 5/25/14 JD
% Added facVecT and FacVecS and facVecF to the cache to support the sampleTest function.
%
% modified 9/3/15 JD
% Added recTime to the cache to support the segment function.
%
% bugfix & modified 1/23/16 JD
% Fixed no support for frequency PCA of continuous data.
% Added support for complex spectral data.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
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

SGLfacs=size(EPdata.data,5);
if ~isempty(EPdata.facData)
    CMBfacs=size(EPdata.facData,5);
end;

EEGchans=find(strcmp('EEG',EPdata.chanTypes));

if strcmp(EPdata.dataType,'continuous')
    numEpochs=ceil(length(EPdata.timeNames)/ceil(EPdata.Fs));
    newDataset.plotMVmin=zeros(numEpochs,1,max(length(EPdata.facNames),1));
    newDataset.plotMVmax=zeros(numEpochs,1,max(length(EPdata.facNames),1));
    for theEpoch=1:numEpochs
        for theFactor=1:SGLfacs
            timeMin=1;
            timeMax=1;
            timeAbsMin=1;
            timeAbsMax=1;
            firstSample=ceil(1+(theEpoch-1)*EPdata.Fs);
            lastSample=ceil(min(theEpoch*EPdata.Fs,length(EPdata.timeNames)));
            if ~isempty(EPdata.facVecT)
                timeMin=min(EPdata.facVecT(:,theFactor));
                timeMax=max(EPdata.facVecT(:,theFactor));
                if ~isempty(EPdata.freqNames)
                    timeAbsMin=min(abs(EPdata.facVecT(:,theFactor)));
                    timeAbsMax=max(abs(EPdata.facVecT(:,theFactor)));
                end;
                firstSample=1;
                lastSample=1;
            end;
            spaceMin=1;
            spaceMax=1;
            spaceAbsMin=1;
            spaceAbsMax=1;
            theChans=EEGchans;
            if ~isempty(EPdata.facVecS)
                spaceMin=min(EPdata.facVecS(EEGchans,theFactor));
                spaceMax=max(EPdata.facVecS(EEGchans,theFactor));
                if ~isempty(EPdata.freqNames)
                    spaceAbsMin=min(abs(EPdata.facVecS(EEGchans,theFactor)));
                    spaceAbsMax=max(abs(EPdata.facVecS(EEGchans,theFactor)));
                end;
                theChans=1;
            end;
            freqMin=1;
            freqMax=1;
            freqAbsMin=1;
            freqAbsMax=1;
            if ~isempty(EPdata.facVecF)
                freqMin=min(EPdata.facVecF(:,theFactor));
                freqMax=max(EPdata.facVecF(:,theFactor));
                if ~isempty(EPdata.freqNames)
                    freqAbsMin=min(abs(EPdata.facVecF(:,theFactor)));
                    freqAbsMax=max(abs(EPdata.facVecF(:,theFactor)));
                end;
            end;
            theRels=1;
            if ~isempty(EPdata.relNames)
                theRels=theChans;
            end;
            if isreal(EPdata.data)
                minFac=min(min(min(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels))));
                maxFac=max(max(max(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels))));
                if ~isempty(EPdata.freqNames)
                    minAbsFac=min(min(min(abs(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels)))));
                    maxAbsFac=max(max(max(abs(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels)))));
                end;
            else
                minFac=min([min(min(min(real(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels))))) min(min(min(imag(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels)))))]);
                maxFac=max([max(max(max(real(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels))))) max(max(max(imag(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels)))))]);
                if ~isempty(EPdata.freqNames)
                    minAbsFac=min(min(min(abs(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels)))));
                    maxAbsFac=max(max(max(abs(EPdata.data(theChans,firstSample:lastSample,:,:,theFactor,:,theRels)))));
                end;
            end;
            if ~isempty(EPdata.freqNames)
                absMinMax=kron(kron(kron([timeAbsMin timeAbsMax],[spaceAbsMin spaceAbsMax]),[freqAbsMin freqAbsMax]),[minAbsFac maxAbsFac]);
            else
                absMinMax=0;
            end;
            minMax=kron(kron(kron([timeMin timeMax],[spaceMin spaceMax]),[freqMin freqMax]),[minFac maxFac]);
            newDataset.plotMVmin(theEpoch,1,theFactor)=min(minMax);
            newDataset.plotMVmax(theEpoch,1,theFactor)=max(minMax);
            newDataset.plotMVminAbs(theEpoch,1,theFactor)=min(absMinMax);
            newDataset.plotMVmaxAbs(theEpoch,1,theFactor)=max(absMinMax);
        end;
        if ~isempty(EPdata.facData)
            for theFactor=1:CMBfacs
                newDataset.plotMVmin(theEpoch,1,theFactor+SGLfacs)=min(min(min(EPdata.facData(EEGchans,firstSample:lastSample,:,:,theFactor,:))));
                newDataset.plotMVmax(theEpoch,1,theFactor+SGLfacs)=max(max(max(EPdata.facData(EEGchans,firstSample:lastSample,:,:,theFactor,:))));
                if ~isempty(EPdata.freqNames)
                    newDataset.plotMVminAbs(theEpoch,1,theFactor+SGLfacs)=min(min(min(abs(EPdata.facData(EEGchans,firstSample:lastSample,:,:,theFactor,:)))));
                    newDataset.plotMVmaxAbs(theEpoch,1,theFactor+SGLfacs)=max(max(max(abs(EPdata.facData(EEGchans,firstSample:lastSample,:,:,theFactor,:)))));
                else
                    newDataset.plotMVminAbs(theEpoch,1,theFactor+SGLfacs)=0;
                    newDataset.plotMVmaxAbs(theEpoch,1,theFactor+SGLfacs)=0;
                end;
            end;
        end;
    end;
else
    newDataset.plotMVmin=zeros(length(EPdata.cellNames),length(EPdata.subNames),max(length(EPdata.facNames),1));
    newDataset.plotMVmax=zeros(length(EPdata.cellNames),length(EPdata.subNames),max(length(EPdata.facNames),1));
    for theCell=1:length(EPdata.cellNames)
        for theSub=1:length(EPdata.subNames)
            for theFactor=1:SGLfacs
                timeMin=1;
                timeMax=1;
                timeAbsMin=1;
                timeAbsMax=1;
                if ~isempty(EPdata.facVecT)
                    timeMin=min(EPdata.facVecT(:,theFactor));
                    timeMax=max(EPdata.facVecT(:,theFactor));
                    if ~isempty(EPdata.freqNames)
                        timeAbsMin=min(abs(EPdata.facVecT(:,theFactor)));
                        timeAbsMax=max(abs(EPdata.facVecT(:,theFactor)));
                    end;
                end;
                spaceMin=1;
                spaceMax=1;
                spaceAbsMin=1;
                spaceAbsMax=1;
                theChans=EEGchans;
                if ~isempty(EPdata.facVecS)
                    spaceMin=min(EPdata.facVecS(EEGchans,theFactor));
                    spaceMax=max(EPdata.facVecS(EEGchans,theFactor));
                    if ~isempty(EPdata.freqNames)
                        spaceAbsMin=min(abs(EPdata.facVecS(EEGchans,theFactor)));
                        spaceAbsMax=max(abs(EPdata.facVecS(EEGchans,theFactor)));
                    end;
                    theChans=1;
                end;
                freqMin=1;
                freqMax=1;
                freqAbsMin=1;
                freqAbsMax=1;
                if ~isempty(EPdata.facVecF)
                    freqMin=min(EPdata.facVecF(:,theFactor));
                    freqMax=max(EPdata.facVecF(:,theFactor));
                    if ~isempty(EPdata.freqNames)
                        freqAbsMin=min(abs(EPdata.facVecF(:,theFactor)));
                        freqAbsMax=max(abs(EPdata.facVecF(:,theFactor)));
                    end;
                end;
                theRels=1;
                if ~isempty(EPdata.relNames)
                    theRels=EEGchans;
                end;
                if isreal(EPdata.data)
                    minFac=min(reshape(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels),1,[])');
                    maxFac=max(reshape(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels),1,[])');
                    if ~isempty(EPdata.freqNames)
                        minAbsFac=min(abs(reshape(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels),1,[])'));
                        maxAbsFac=max(abs(reshape(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels),1,[])'));
                    end;
                else
                    minFac=min([min(reshape(real(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels)),1,[])') min(reshape(imag(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels)),1,[])')]);
                    maxFac=max([max(reshape(real(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels)),1,[])') max(reshape(imag(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels)),1,[])')]);
                    if ~isempty(EPdata.freqNames)
                        minAbsFac=min(abs(reshape(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels),1,[])'));
                        maxAbsFac=max(abs(reshape(EPdata.data(theChans,:,theCell,theSub,theFactor,:,theRels),1,[])'));
                    end;
                end;
                if ~isempty(EPdata.freqNames)
                    absMinMax=kron(kron(kron([timeAbsMin timeAbsMax],[spaceAbsMin spaceAbsMax]),[freqAbsMin freqAbsMax]),[minAbsFac maxAbsFac]);
                else
                    absMinMax=0;
                end;
                minMax=kron(kron(kron([timeMin timeMax],[spaceMin spaceMax]),[freqMin freqMax]),[minFac maxFac]);
                newDataset.plotMVmin(theCell,theSub,theFactor)=min(minMax);
                newDataset.plotMVmax(theCell,theSub,theFactor)=max(minMax);
                newDataset.plotMVminAbs(theCell,theSub,theFactor)=min(absMinMax);
                newDataset.plotMVmaxAbs(theCell,theSub,theFactor)=max(absMinMax);
            end;
            if ~isempty(EPdata.facData)
                for theFactor=1:CMBfacs
                    newDataset.plotMVmin(theCell,theSub,theFactor+SGLfacs)=min(min(min(EPdata.facData(EEGchans,:,theCell,theSub,theFactor,:,theRels))));
                    newDataset.plotMVmax(theCell,theSub,theFactor+SGLfacs)=max(max(max(EPdata.facData(EEGchans,:,theCell,theSub,theFactor,:,theRels))));
                    if ~isempty(EPdata.freqNames)
                        newDataset.plotMVminAbs(theCell,theSub,theFactor+SGLfacs)=min(min(min(abs(EPdata.facData(EEGchans,:,theCell,theSub,theFactor,:,theRels)))));
                        newDataset.plotMVmaxAbs(theCell,theSub,theFactor+SGLfacs)=max(max(max(abs(EPdata.facData(EEGchans,:,theCell,theSub,theFactor,:,theRels)))));
                    else
                        newDataset.plotMVminAbs(theCell,theSub,theFactor+SGLfacs)=0;
                        newDataset.plotMVmaxAbs(theCell,theSub,theFactor+SGLfacs)=0;
                    end;
                end;
            end;
        end;
    end;
end;

newDataset.dataName=EPdata.dataName;
newDataset.chanNames=EPdata.chanNames;
newDataset.timeNames=EPdata.timeNames;
newDataset.cellNames=EPdata.cellNames;
newDataset.trialNames=EPdata.trialNames;
newDataset.subNames=EPdata.subNames;
newDataset.facNames=EPdata.facNames;
newDataset.freqNames=EPdata.freqNames;
newDataset.relNames=EPdata.relNames;
newDataset.eloc=EPdata.eloc;
newDataset.dataType=EPdata.dataType;
newDataset.chanTypes=EPdata.chanTypes;
newDataset.subTypes=EPdata.subTypes;
newDataset.cellTypes=EPdata.cellTypes;
newDataset.facTypes=EPdata.facTypes;
newDataset.implicit=EPdata.implicit;
newDataset.baseline=EPdata.baseline;
newDataset.Fs=EPdata.Fs;
newDataset.ced=EPdata.ced;
newDataset.saved='yes';
newDataset.trialSpecNames=EPdata.trialSpecNames;
newDataset.trialSpecs=EPdata.trialSpecs;
newDataset.events=EPdata.events;
newDataset.pca=EPdata.pca;
newDataset.subjectSpecNames=EPdata.subjectSpecNames;
newDataset.timeNames=EPdata.timeNames;
newDataset.facVecT=EPdata.facVecT;
newDataset.facVecS=EPdata.facVecS;
newDataset.facVecF=EPdata.facVecF;
newDataset.recTime=EPdata.recTime;
newDataset.timeUnits=EPdata.timeUnits;

