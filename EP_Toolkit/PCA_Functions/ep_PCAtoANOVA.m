function ep_PCAtoANOVA(FactorResults, changrps, inputcells, startwin, endwin, measure, factor, fname, stand);

%  ep_PCAtoANOVA - ep_PCAtoANOVA(FactorResults, changrps, inputcells, startwin, endwin, measure, factor, fname, stand) - generates summary measures for ANOVAs from PCA results
%
%  Outputs text file suitable for ANOVA from PCA results.  Each line is a subject.  As for columns, cells will vary slowest while chan groups will
%  vary fastest.  Chan groups will be generated for spatial as well as temporal PCAs even though chan effects will not occur for spatial PCAs.
%  Results will be in microvolts (or msec) for a specific channel (or mean of group of channels) for a specific window to make results exactly comparable
%  to conventional analyses.
%
%Inputs:
%  As fields in the FactorResults structure
%  .PCAmode	 : Format for output data ('temp': time points as columns; 'spat': channels as columns).
%  .FacScr 	 : is the standardized factor score matrix (observations, factors)
%  .FacPat	 : Factor pattern matrix - produces standardized variables from scores, scaled by communality  (variables, factors)
%  .varSD     : Standard deviations of the variables.
%  .numchan	 : Number of channels.
%  .timePoints: Number of timepoints.
%  .numSubs   : Number of subjects.
%  .Fs          : The sampling frequency in Hz.
%  .baseline    : The number of samples in the baseline.
%
%  changrps : Channel going into each output channel grouping (channel group, channels).  Should be padded with zeros to maximum number of input chans.
%				Example:
%					changrps = [6 8 0 0;...
%					1 2 3 4];
% inputcells : Cells included in analysis.  Should be specified in order desired for output.
% startwin   : Start of time window in samples (1-based).
% endwin	 : End of time window in samples.
% measure	 : Type of measure for windows: mean, minpeak, maxpeak, minlatency, maxlatency, mincentroid, maxcentroid.
% factor	 : Factor to be output
% fname		 : name of output file
% stand      : Whether to scale the data or leave standardized  ('st'=standardized, 'sc'=scaled).  If standardized, different factors can be compared.  If not, it's in microvolts.
%
%  Citation for centroid measure of latency:
%  Dien, J., Spencer, K. M., & Donchin, E. (2004). Parsing the "Late Positive Complex": Mental chronometry and the ERP
%  components that inhabit the neighborhood of the P300. Psychophysiology, 41(5), 665-678.

%History:
%  by Joseph Dien (4/4/01)
%  jdien07@mac.com
%
%  modified 10/4/01 JD
%  Added ability to output standardized factor scores so that different factors can be compared within a single ANOVA.
%
%  modified 11/19/02 JD
%  Added error message for when Endwindow is larger than timepoints.  Latency output in milliseconds rather than samples.
%
%  modified 11/24/02 JD
%  Added centroid option to types of measures.
%
%  bugfix 5/29/05 JD
%  Fixed channel groups feature for spatial PCAs
%
%  modified 11/9/06 JD
%  added error check for stand variable
%
%  modified (2/6/08) JD
%  Copyright notice appears only once per session.  Input factor results are now packaged in a
%  structured variable.  Accommodates two-step PCA results by itself.
%
%  modified (11/7/08) JD
%  Fixed incorrect initialization that was slowing things down dramatically.
%
% modified (2/25/09) JD
% Gets baseline and sampling rate from structured variable.  Assumes bins equals 1.
%
% modified (1/29/14) JD
% Updated field names.
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

if nargin < 14
  stand = 'SC';
end

if ~any(strcmp(stand,{'ST','SC'}))
    error('stand parameter must equal either SC for scaled or ST for standardized');
end;

if ~isfield(FactorResults,'varSD')
    error('No varSD field.');
end
varSD=FactorResults.varSD;
if ~isfield(FactorResults,'numchan')
    error('No numchan field.');
end
numchans=FactorResults.numchan;
if ~isfield(FactorResults,'timepoints')
    error('No timepoints field.');
end
timepoints=FactorResults.timepoints;
if ~isfield(FactorResults,'numSubs')
    error('No numSubs field.');
end
numSubs=FactorResults.numSubs;
if ~isfield(FactorResults,'numCells')
    error('No numCells field.');
end
numCells=FactorResults.numCells;

if ~isfield(FactorResults,'Fs')
    error('No Fs field.');
end
SampleRate=FactorResults.Fs;
dig=1000/SampleRate;

if ~isfield(FactorResults,'baseline')
    error('No baseline field.');
end
baseline=FactorResults.baseline*(1000/SampleRate);

if endwin < startwin
	error('Endwindow sample is smaller than startwindow sample.');
end;

if endwin > timepoints
	error('Endwindow sample is larger than number of time points.');
end;

if isfield(FactorResults,'PCAmode2')  %if this is the result of a two-step PCA
    tempVarSD = ones(1,size(FactorResults.varSDST,2));	%since varSD is not the same for all the factors, set this to one and then multiply the pattern matrix by the appropriate numbers.

    numfacs2 = FactorResults.numFacs2;	%number of factors retained in second PCA step
    numfacs1 = FactorResults.numFacs;	%number of factors retained in first PCA step
    tempFacPatST=zeros(size(FactorResults.FacPatST,1), size(FactorResults.FacPatST,2));
    for i = 1:numfacs1
        for k = 1:numfacs2
            tempFacPatST(:,(i-1)*numfacs2+k)=diag(FactorResults.varSDST(i,:))*(FactorResults.FacPatST(:,(i-1)*numfacs2+k));
        end;
    end;
    PCAmode=FactorResults.PCAmode2;
    if strcmp(PCAmode,'temp')
        tempVar = numchans;
    elseif strcmp(PCAmode,'spat')
        tempVar = timepoints;
    else
        error('PCAmode in factor results structure must be set to either temp or spat');
    end;
    tempFacScrST=zeros(size(FactorResults.FacScrST,1)*tempVar, size(FactorResults.FacScrST,2));
%     for i = 1:numfacs1
%         for j = 1:numfacs2
%             for k = 1:(numSubs*numCells)
%                 tempFacScrST(1+tempVar*(k-1):tempVar*k,(i-1)*numfacs2+j)=diag(FactorResults.varSD)*FactorResults.FacPat(:,i)*FactorResults.FacScrST(k,(i-1)*numfacs2+j);
%             end;
%         end;
%     end;
    
    for i = 1:numfacs1
        for j = 1:numfacs2
            tempFacScrST(:,(i-1)*numfacs2+j)=reshape((diag(FactorResults.varSD)*FactorResults.FacPat(:,i))*FactorResults.FacScrST(:,(i-1)*numfacs2+j)',[],1);
        end;
    end;    
    
    PCAmode=FactorResults.PCAmode2;
    FacScr=tempFacScrST;
    FacPat=tempFacPatST;
    varSD=tempVarSD;
elseif isfield(FactorResults,'PCAmode')
    PCAmode=FactorResults.PCAmode;
    if ~isfield(FactorResults,'FacScr')
        error('No FacScr field.');
    end
    FacScr=FactorResults.FacScr;
    if ~isfield(FactorResults,'FacPat')
        error('No FacPat field.');
    end
    FacPat=FactorResults.FacPat;
else
    error('No PCAmode field.');
end;

if stand == 'SC'
	FacPat=diag(varSD)*FacPat;	%prepare for conversion
elseif stand == 'ST'
	FacPat=ones(size(FacPat,1),size(FacPat,2)); %set to ones so that the scores are output still standardized
else
	error('stand parameter must equal either SC for scaled or ST for standardized');
end


outputdata = zeros(numSubs,(size(changrps,1)*length(inputcells)));

if PCAmode == 'temp'
	for subject = 1:numSubs
		for cell = 1:length(inputcells)
		    for group = 1:size(changrps,1)
				data = 0;
				count = 0;
				for withingrp = 1:size(changrps,2)
					if changrps(group, withingrp) ~= 0
						data = data + FacScr(((inputcells(cell)-1)*numSubs*numchans)+((subject-1)*numchans)+changrps(group, withingrp),factor);
						count = count +1;
					end;
				end;
				data = data /count;
				outputdata(subject,(cell-1)*size(changrps,1)+group) = data;
			end;
		end;
	end;

	switch measure
		case 'mean', outputdata = mean(FacPat(startwin:endwin,factor)).*outputdata;
		case 'minpeak', outputdata = min(FacPat(startwin:endwin,factor)).*outputdata;
		case 'maxpeak', outputdata = max(FacPat(startwin:endwin,factor)).*outputdata;
		case 'minlatency', [dummy, latency] = min(FacPat(startwin:endwin,factor)); outputdata = ones(numSubs,(size(changrps,1)*length(inputcells))).*(latency * dig) + baseline;
		case 'maxlatency', [dummy, latency] = max(FacPat(startwin:endwin,factor)); outputdata = ones(numSubs,(size(changrps,1)*length(inputcells))).*(latency * dig) + baseline;
        case 'mincentroid',
            epoch=FacPat(startwin:endwin,factor);
            maxvalue=max(epoch);
            outputdata = sum([startwin:endwin]'.*(epoch-maxvalue)/sum(epoch-maxvalue)).*outputdata*dig + baseline;
        case 'maxcentroid',
            epoch=FacPat(startwin:endwin,factor);
            minvalue=min(epoch);
            outputdata = sum([startwin:endwin]'.*(epoch-minvalue)/sum(epoch-minvalue)).*outputdata*dig + baseline;
		otherwise, error('Measure should be one of the following options: mean, minpeak, maxpeak, minlatency, maxlatency, mincentroid, maxcentroid');
	end;
        
elseif PCAmode == 'spat'
	for subject = 1:numSubs
		for cell = 1:length(inputcells)
		    for group = 1:size(changrps,1)
				data = 0;
				count = 0;
				for withingrp = 1:size(changrps,2)
					if changrps(group, withingrp) ~= 0
						data = data + FacPat(changrps(group, withingrp), factor);
						count = count +1;
					end;
				end;
				data = data /count;
				startscore = ((inputcells(cell)-1)*numSubs*timepoints)+((subject-1)*timepoints)+startwin;
				endscore = ((inputcells(cell)-1)*numSubs*timepoints)+((subject-1)*timepoints)+endwin;
				switch measure
					case 'mean', outputdata(subject,(cell-1)*size(changrps,1)+group) = mean(FacScr(startscore:endscore,factor))*data;
					case 'minpeak', outputdata(subject,(cell-1)*size(changrps,1)+group) = min(FacScr(startscore:endscore,factor))*data;
					case 'maxpeak', outputdata(subject,(cell-1)*size(changrps,1)+group) = max(FacScr(startscore:endscore,factor))*data;
					case 'minlatency', [dummy, latency] = min(FacScr(startscore:endscore,factor)); outputdata(subject,(cell-1)*size(changrps,1)+group) = (latency * dig) + baseline;
					case 'maxlatency', [dummy, latency] = max(FacScr(startscore:endscore,factor)); outputdata(subject,(cell-1)*size(changrps,1)+group) = (latency * dig) + baseline;
                    case 'mincentroid',
                        epoch=FacScr(startscore:endscore,factor);
                        maxvalue=max(epoch);
                        outputdata(subject,(cell-1)*size(changrps,1)+group) = sum([startwin:endwin]'.*(epoch-maxvalue)/sum(epoch-maxvalue))*dig + baseline;
                    case 'maxcentroid',
                        epoch=FacScr(startscore:endscore,factor);
                        minvalue=min(epoch);
                        outputdata(subject,(cell-1)*size(changrps,1)+group) = sum([startwin:endwin]'.*(epoch-minvalue)/sum(epoch-minvalue))*dig + baseline;
                    otherwise, error('Measure should be one of the following options: mean, minpeak, maxpeak, minlatency, maxlatency, mincentroid, maxcentroid');
                end;
			end;
		end;
	end;
else
	error('PCAmode should be either ''temp'' for temporal or ''spat'' for spatial');
end;

if exist(fname,'file')
    msg{1}=[fname ' already exists!'];
    [msg]=ep_errorMsg(msg);
    return
end;

outFID=fopen(fname,'w');
if (outFID == -1)
    msg{1}='Error creating output file!';
    [msg]=ep_errorMsg(msg);
    return
end;

%file header
fprintf(outFID,'%s\r',fname);
fprintf(outFID,'%d-%d',num2str(startwin),num2str(endwin));
fprintf(outFID,'\r');


fprintf(outFID,'%s\r',measure);


fprintf(outFID,'\r');
fprintf(outFID,'\r');

for theCell = 1:length(inputcells)
    for theArea = 1:size(changrps,1)
        fprintf(outFID,'%s\t',num2str(inputcells(theCell)));
    end
end;

fprintf(outFID,'\r');

for theCell = 1:length(inputcells)
    for theArea = 1:size(changrps,1)
        fprintf(outFID,'%s\t',num2str(theArea));
    end
end;

fprintf(outFID,'\r');

%file data
for theSub=1:numSubs
    for theCell = 1:size(outputdata,2)
        fprintf(outFID,'%f\t',outputdata(theSub,theCell));
    end;
    fprintf(outFID,[num2str(theSub) '\r']);
end;

fclose(outFID);