function [expandedData]=ep_expandFacs(EPdata,chans,points,cells,subs,facs,freqs,rels);
%  [expandedData]=ep_expandFacs(EPdata,chans,points,cells,subs,facs,freqs);
%       Uncompresses factor data by multiplying loadings (in facVecS and facVecT and facVecF) by the factor scores.
%       When there are rels, only provides those corresponding to the channels specified.
%
%Inputs:
%  EPdata         : Structured array with the data and accompanying information in EP file format.  See readData.
%  chans          : Channels to expand.
%  points         : Time points to expand.
%  cells          : Cells to expand.
%  subs           : Subjects to expand.
%  facs           : Factors to expand.
%  freqs          : Freqs to expand
%  rels           : Rels to expand
%
%Outputs:
%  expandedData        : Uncompressed voltage data, as specified by the inputs

%History:
%  by Joseph Dien (7/17/09)
%  jdien07@mac.com

% bugfix 2/3/11 JD
% Fixed only final factor being expanded for spatial PCAs.
%
% modified 2/1/12 JD
% Added support for 6th dimension of frequency
%
% bugfix 5/30/13 JD
% Fixed channel of non-spatial PCAs being set to 1 (when windowing for the most part) rather than the actual channel.
%
% modified 4/24/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% bugfix 7/27/14 JD
% Fixed time point of spatial PCAs and frequency of temporal and spatial PCAs of FFT/TFT data indexing from first sample
% rather than from first specified sample when windowing.

%
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

if nargin < 2
    chans=[1:length(EPdata.chanNames)];
end;

if nargin < 3
    points=[1:length(EPdata.timeNames)];
end;

if nargin < 4
    cells=[1:length(EPdata.cellNames)];
end;

if nargin < 5
    subs=[1:length(EPdata.subNames)];
end;

if nargin < 6
    facs=[1:length(EPdata.facNames)];
end;

if nargin < 7
    facs=[1:length(EPdata.freqNames)];
end;

if nargin < 8
    rels=[1:length(EPdata.relNames)];
end;

if isempty(chans)
    chans=[1:length(EPdata.chanNames)];
end;

if isempty(points)
    if isempty(EPdata.timeNames)
        points=1;
    else
        points=[1:length(EPdata.timeNames)];
    end;
end;

if isempty(cells)
    cells=[1:length(EPdata.cellNames)];
end;

if isempty(subs)
    subs=[1:length(EPdata.subNames)];
end;

if isempty(facs)
    facs=[1:length(EPdata.facNames)];
end;

if isempty(freqs)
    if isempty(EPdata.freqNames)
        freqs=1;
    else
        freqs=[1:length(EPdata.freqNames)];
    end;
end;

if isempty(facs)
    facs=1;
end;

if isempty(rels)
    rels=[1:length(EPdata.relNames)];
    if isempty(rels)
        rels=1;
    end;
end;

SGLfacs=facs(find(facs <= size(EPdata.data,5)));
if ~isempty(EPdata.facData)
    CMBfacs=facs(find(facs > size(EPdata.data,5)))-size(EPdata.data,5);
end;

expandedData=zeros(length(chans),length(points),length(cells),length(subs),length(facs),length(freqs),length(rels));

if isempty(EPdata.facVecS) && isempty(EPdata.facVecT) && isempty(EPdata.facVecF)
    expandedData=EPdata.data(chans,points,cells,subs,SGLfacs,freqs,rels);
else
    for i=1:length(chans)
        theChan=chans(i);
        theComChan=theChan;
        for k=1:length(points)
            thePoint=points(k);
            theComPoint=thePoint;
            for n=1:length(freqs)
                theFreq=freqs(n);
                theComFreq=theFreq;
                for m=1:length(SGLfacs)
                    theFac=SGLfacs(m);
                    theNum=1;
                    if ~isempty(EPdata.facVecS)
                        theNum=theNum*EPdata.facVecS(theChan,theFac);
                        theComChan=1;
                    end;
                    if ~isempty(EPdata.facVecT)
                        theNum=theNum*EPdata.facVecT(thePoint,theFac);
                        theComPoint=1;
                    end;
                    if ~isempty(EPdata.facVecF)
                        theNum=theNum*EPdata.facVecF(theFreq,theFac);
                        theComFreq=1;
                    end;
                    expandedData(i,k,:,:,m,n,:)=EPdata.data(theComChan,theComPoint,cells,subs,theFac,theComFreq,rels)*theNum;
                end;
            end;
        end;
    end;
end;

if ~isempty(EPdata.facData)
    expandedData(:,:,:,:,length(SGLfacs)+1:length(SGLfacs)+length(CMBfacs),:,:)=EPdata.facData(chans,points,cells,subs,CMBfacs,freqs,rels);
end;