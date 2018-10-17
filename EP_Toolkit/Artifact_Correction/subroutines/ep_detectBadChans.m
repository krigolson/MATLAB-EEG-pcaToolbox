function [badChans shortChans outputLog]=ep_detectBadChans(EPdata, badDataCriteria, theSubject);
%  [badChans shortChans outputLog]=ep_detectBadChans(EPdata, badDataCriteria, theSubject);
%       Detects bad channels by identifying ones that cannot be readily
%       predicted by the neighboring channels and by detecting flat data channels.
%       Also notes shorted channels.
%
% Cohen, J., & Cohen, P. (1983). Applied multiple regression/correlation analysis for the behavioral sciences. Hillsdale, NJ: Lawrence Erlbaum Associates. 
%
%Inputs:
%  EPdata         : Structured array with the data and accompanying information.  See readData.
%   badDataCriteria:  Criteria for detecting bad data.
%       .neighbors: number of electrodes considered to be neighbors
%       .badchan:   minimum predictability (multiple R) from neighbors to not be considered globally bad
%       .badtrials: percentage of good trials chan is bad to declare a channel globally bad
%       .saturation: followed by range of acceptable data values.  Time points with a channel outside this range will be excluded.
%   theSubject: which subject of the file is being processed.
%
%Outputs:
%  badChans : List of bad channels.
%  shortChans: List of shorted channels.
%  outputLog: output messages from bad channel detection process

%History:
%  by Joseph Dien (2/8/09)
%  jdien07@mac.com
%
% modified 3/14/09 JD
% Changed to use EP format data to provide more flexibility with I/O functions.
%
% modified 5/15/09 JD
% Treats flat channels as bad data and excludes from bad channel detection routine.
% Dropped eloc as separate input parameter (now part of data).
% Flat channel not bad if it is the reference channel.
%
% modified 9/4/09 JD
% Added support for multiple refChans to deal with mean mastoid data where the presence of the two reference channels (correlated -1)
% was causing ICA problems.
%
% modified 10/28/09 JD
% Added detection of channels perfectly correlated with a reference channel and which were therefore flat prior to rereferencing.
%
% modified 11/12/09 JD
% Correlated channels can be only nearly perfect (e.g., .9999) and still trigger bad channel code, to account for rounding errors etc.
%
% bugfix 11/20/09 JD
% Replaced "union" commands with "unique" commands because certain situations caused the "union" command to crash in
% Matlab 2007.
%
% modified & bugfix 12/3/09 JD
% Detects non-reference channels that are perfectly correlated and identifies them as bad channels as they must be
% shorted together.  Fixed bug where test of correlation with reference only detecting +1 correlation, not -1 correlation.
% Fixed bug where if there is an explicit reference channel and it is flat, then all reference channels marked bad and
% real bad channels are no longer marked bad.
% Don't apply correlated neighbors test to the reference channels as distant reference channels will always be labeled
% bad.
%
% modified & bugfix 2/24/10 JD
% Now works on average files.
% Fixed bug where neighboring channels for determing whether a channel is not correlating with its neighbors sometimes not chosen
% correctly, which could lead to too many channels being dubbed globally bad.
% No longer treating shorted channels as being bad (too conservative).  Instead just displaying a warning message.
% If there are two reference channels (as in mean mastoids), then no longer require that they have a -1 correlation as one may just be bad.
% If there are two reference channels (as in mean mastoids), then they are still marked as bad channels if they are
% flat.
% Added log output.
% When there are shorted channels, prints out the channel pairs.
%
% bugfix 6/7/11 JD
% Fixed crash when bad channels detection in artifact correction routine generated error message.
% 
% modified 1/25/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
% 
% modified 9/22/13 JD
% Restricted bad channel detection to EEG channels.
% 
% modified 5/23/16 JD
% Detects channels that went flat partway through a session (more than 10%) and labels them as global bad channels.
%
% modified 4/26/17 JD
% Excludes time points beyond a certain range from global bad channel detection, blink, and saccade routines.
% Divides the data into sections according to the badtrials parameter and checks the badchan multiple R predictability threshold for each section to detect channels that went bad partway through the session.
%
% bugfix 9/30/17 JD
% Median corrects data prior to saturation check to ensure channels with merely high offsets are not treated as bad data.
% Median based on the section rather than the entire dataset.
%
% bugfix 11/26/17 JD
% Fixed calculation of section median bad channel.
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

badChans{1}=-1;
outputLog={};

if nargin < 3
    theSubject =1;
end;

if strcmp(EPdata.dataType,'factors')
    msg='This function does not support factor files.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end;

elecDistances=ep_closestChans(EPdata.eloc);

numSubs=length(EPdata.subNames);
shortChans=[];
numChans=length(EPdata.chanNames);
if length(EPdata.eloc) ~= numChans
    msg='Error: The number of channels in the electrode coordinates file is different from that of the data.';
    disp(msg);
    outputLog{end+1}=msg;    
    badChans=-1;
    return;
end;

EEGchans=find(strcmp('EEG',EPdata.chanTypes));
testData=reshape(EPdata.data(EEGchans,:,:,theSubject),length(EEGchans),[]);
goodPoints = find((max(testData-repmat(median(testData')',1,size(testData,2))) < badDataCriteria.saturation(2)) & (min(testData-repmat(median(testData')',1,size(testData,2))) > badDataCriteria.saturation(1)));
if isempty(goodPoints) || (length(goodPoints) < (size(testData,2)/10))
    goodPoints=[1:size(testData,2)];
    msg='Too many points exceeded saturation check for bad channel detection so dropping this criterion.';
    disp(msg);
    outputLog{end+1}=msg;    
end;
testData=testData(:,goodPoints);
allData=reshape(EPdata.data(EEGchans,:,:,theSubject),length(EEGchans),[]);
normData=pinv(diag(std(testData,0,2)))*testData;
normData=normData-diag(mean(normData,2))*ones(size(normData));
badChans=find(~std(normData'))'; %flat channels are bad.
goodChans=find(std(normData'))';
elecDistances(badChans,:)=inf;
elecDistances(:,badChans)=inf;
origRefChan=EPdata.reference.original;
currRefChan=EPdata.reference.current;

%global bad channel is flat for more than 10% of the segments
if strcmp(EPdata.dataType,'continuous')
    numTrials=floor(size(EPdata.data,3)/EPdata.Fs); %ignore excess
    flatTrialsChans=zeros(length(goodChans),numTrials);
    for iTrial = 1:numTrials
        flatTrialsChans(:,iTrial)=var(squeeze(EPdata.data(goodChans,(iTrial-1)*EPdata.Fs+1:iTrial*EPdata.Fs,1,theSubject))');
    end;
else
    numTrials=size(EPdata.data,3);
    flatTrialsChans=zeros(length(goodChans),numTrials);
    for iTrial = 1:numTrials
        flatTrialsChans(:,iTrial)=var(squeeze(EPdata.data(goodChans,:,iTrial,theSubject))');
    end;
end;
badChans=[badChans; goodChans(find(sum(flatTrialsChans'==0)>(numTrials/10)))];
goodChans=setdiff(goodChans,badChans);

if isempty(goodChans)
    msg='Error: No good EEG channels left.';
    disp(msg);
    outputLog{end+1}=msg; 
    badChans=-1;
    return;
end;

corrs=corrcoef([testData']);
corrSigns=sign(corrs);
corrs=corrSigns.*ceil(abs(corrs)*1000)/1000;

%look for bad channels that are perfectly correlated with each other and must therefore be shorted together.

goodNonRefChans=setdiff(goodChans,currRefChan);

shortPairs=[];
for chan1= 1:length(goodNonRefChans)
    theChan1= goodNonRefChans (chan1);
    for chan2=chan1+1:length(goodNonRefChans)
        theChan2= goodNonRefChans (chan2);
        if abs(corrs(theChan1, theChan2))== 1
            shortChans=unique([shortChans theChan1, theChan2]);
            shortPairs=[shortPairs EPdata.chanNames{theChan1} '-' EPdata.chanNames{theChan2} '; '];
        end;
    end;
end;
if ~isempty(shortPairs)
    msg=['Warning: shorted channels: ' shortPairs];
    disp(msg);
    outputLog{end+1}=msg;
end;

%look for bad channels that are not correlated with other channels
numRegressors=min(badDataCriteria.neighbors,length(goodChans)-1);
neighbors=zeros(length(EEGchans),numRegressors);
for iChan=1:length(goodChans)
    chan=goodChans(iChan);
    [E IX]=sort(elecDistances(chan,goodChans));
    neighbors(chan,:)=goodChans(IX(2:numRegressors+1));
end;

secLength=floor((badDataCriteria.badtrials/100)*length(goodPoints));
secAllLength=floor((badDataCriteria.badtrials/100)*size(EPdata.data,2));
if secLength >= EPdata.Fs %if section is at least a second long
    numSecs=floor(100/badDataCriteria.badtrials);
else
    numSecs=1;
    secLength=length(goodPoints);
    secAllLength=size(EPdata.data,2);
end;

chanPred=zeros(length(EEGchans),numSecs);
badPred=zeros(length(EEGchans),1);
for iSec=1:numSecs
    secPoints=((iSec-1)*secLength)+1:iSec*secLength;
    secAllPoints=((iSec-1)*secAllLength)+1:iSec*secAllLength;
    for iChan=1:length(goodChans)
        theChan=goodChans(iChan);
        Y=normData(theChan,secPoints)';
        X=normData(neighbors(theChan,:),secPoints)';
        R=corrcoef([Y X]);
        R=R(2:end,1);
        B=[ones(secLength,1) X]\Y;
        chanPred(theChan,iSec)=sqrt(sum(B(2:end).*R)); %Multiple R Cohen & Cohen (1983), p. 86
        if chanPred(theChan,iSec) < badDataCriteria.badchan
            badPred(theChan)=1;
        end;
        if (length(find(abs(allData(theChan,secAllPoints)-median(allData(theChan,secAllPoints)))>500))/secLength) > (badDataCriteria.badtrials/100)
            badPred(theChan)=1;
        end;
    end;
end;

nonRefGoodChans=setdiff(goodChans,currRefChan);
badChans=[badChans; nonRefGoodChans(find(badPred(nonRefGoodChans)))]; %don't check ref chans for being bad by local channel predictability

if length(currRefChan) == 1
    if std(normData(currRefChan,:)')==0
        badChans=setdiff(badChans,currRefChan); %flat channel is not bad if it is the reference channel
    end;
end;

goodChans=setdiff(goodChans,badChans);

if isempty(goodChans)
    msg='Error: No good EEG channels left.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end;

%look for bad channels that are perfectly correlated with the reference channel(s) and were therefore flat prior to
%rereferencing.

if length(origRefChan) > 2
    msg='There are more than two original reference channels indicated.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
elseif length(origRefChan) == 2
    refCorrs=abs(corrs(origRefChan(1),:))';
    badChans=unique([badChans; setdiff(find(refCorrs == 1),origRefChan)]);
    badChans=intersect(badChans,EEGchans);
elseif length(origRefChan) == 1
    refCorrs=abs(corrs(origRefChan(1),:))';
    badChans=unique([badChans; setdiff(find(refCorrs == 1),origRefChan)]);
    badChans=intersect(badChans,EEGchans);
end;

if ~strcmp(EPdata.reference.type,'AVG')
    %If there are two current reference channels and they are not perfectly inversely correlated, there is a problem
    if length(origRefChan) == 2
        if (corrs(currRefChan(1),currRefChan(2)) ~= -1) && isempty(intersect(currRefChan,badChans))
            msg=['The two current reference channels should have a perfect inverse correlation and do not (' num2str(corrs(currRefChan(1),currRefChan(2))) ') so something is wrong.'];
            disp(msg);
            outputLog{end+1}=msg;
            badChans=-1;
            return;
        end;
    end;
end;

goodChans=setdiff(goodChans,badChans);

if isempty(goodChans)
    msg='Error: No good channels left.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end;





