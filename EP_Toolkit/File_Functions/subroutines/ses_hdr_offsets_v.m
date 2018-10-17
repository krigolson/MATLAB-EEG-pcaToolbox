%ses_hdr_offsets.m
%
%Declares a set of global constants for reading the fields of the
%EGIS session header vectors returned by function rd_ses_hdr_vecs.
%
%Note: Because of name conflicts this file should not be used with
%the non-vectorized script rd_ses_hdr.
%
%See also:rd_ses_hdr_vecs, wt_ses_hdr_vecs.

%Modification History:
%	Original Version by Brian Rakitin, 3/2/95.
%	06/05/95 PJ Equated NObs with NTrials

%Session File Description from "EGIS Hdr Conv 1.02"
BytOrd=1;
HdrVer=2;
LHeader=3;
LData=4;

%Session File Session Information
RunDate=5;
RunDateMo=5;
RunDateDay=6;
RunDateYr=7;
RunTime=8;
RunTimeHr=8;
RunTimeMin=9;
RunTimeSec=10;
SubjID=11;
Handed=12;
Sex=13;
Age=14;
ExperID=15;
EdVer=16;
CalFlag=17;

%Session File Data Description
NCells=18;
NChan=19;
LComment=20;
LText=21;
LPad=22;
BrdGain=23;
LCellHdr=24;

%Cell Header Offsets
CellID=1;
NTrials=2;
NObs=2;
NPoints=3;
NSamps=3;
NSamp=3;
SampRate=4;
LSpec=5;
TSpecs=6;
