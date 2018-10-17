function ep_contrastANOVA
% ep_contrastANOVA -
% Allows contrasts to be specified for the robust ANOVA module.

%History
%  by Joseph Dien (7/23/09)
%  jdien07@mac.com
%
%  bugfix 8/27/09 JD
%  Location of windows appearing partly off screen on some systems.  Fixed.
%
%  bugfix 1/20/11 JD
%  Fixed levels of ANOVA factors being listed in the wrong order when there are more than one ANOVA factor.
%  Fixed crash when no between factor specified.
%  Fixed error message incorrectly rejecting between or within contrast specified as being "1", which should mean
%  no contrast of that type.
%
%  bugfix 3/24/11 JD
%  Fixed crash when performing contrasts with dataset having no between group factors.
%
%  bugfix 11/1/13 JD
%  Fixes font sizes on Windows.
%
%  bugfix 1/12/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% modified 12/11/15 JD
% Checks if p-value variability exceeds two standard deviations over the threshold.
% Includes warning summary of %age of singular and nearly singular matrices in the output.
%
% bugfix 6/12/17 JD
% Fixed crash when using contrast function.
% Only presents the between contrast controls if there are more than one between levels.
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

global EPmain EPContrast

scrsz = EPmain.scrsz;
windowHeight=scrsz(4);

contrastFigure=findobj('Name', 'Contrast Window');

if ~isempty(contrastFigure)
    close(contrastFigure)
end;

EPContrast.handles.figure=figure('Name', 'Contrast Window', 'NumberTitle', 'off', 'Position',[201 1 700 windowHeight]);
colormap jet;

%initial error checking of settings
for i=1:6
    if xor(isempty(EPmain.anova.data.betweenName{i}),(EPmain.anova.data.between(i)>length(EPmain.anova.data.columnNames)))
        msg{1}='All between-group independent variables need to have a three letter label specified.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if ~isempty(EPmain.anova.data.factor{i}) && isempty(EPmain.anova.data.levels{i})
        msg{1}='There is a within-group factor with a name but no levels.';
        [msg]=ep_errorMsg(msg);
        return
    end
end;

set(EPmain.handles.anova.contrast,'enable','off');
set(EPmain.handles.anova.load,'enable','off');
set(EPmain.handles.anova.view,'enable','off');
set(EPmain.handles.anova.run,'enable','off');
set(EPmain.handles.anova.done,'enable','off');
drawnow

EPContrast.handles.run = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-130 60 30], 'Callback', @ANOVAcontrast);

EPContrast.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-160 60 30], 'Callback', ['close(''Contrast Window'');', 'ep(''start'');']);

%within-group levels table
totalFactorLevels=EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1;

factorLevels=cell(0);
factorNames=cell(0);
for theFactor=1:length(EPmain.anova.data.factor)
    if ~isempty(EPmain.anova.data.factor{theFactor}) && ~isempty(EPmain.anova.data.levels{theFactor})
        factorLevels{end+1}=EPmain.anova.data.levels{theFactor};
        factorNames{end+1}=EPmain.anova.data.factor{theFactor};
    end;
end;

repFactor=1;
repCount=1;
lvlCount=1;
numFactors=length(factorLevels);
for theFactor=numFactors:-1:1
    for theRow=1:totalFactorLevels
        tableDataW{theRow,theFactor+1}=factorLevels{theFactor}(lvlCount);
        repCount=repCount+1;
        if repCount > repFactor
            repCount=1;
            lvlCount=lvlCount+1;
            if lvlCount > length(factorLevels{theFactor})
                lvlCount=1;
            end;
        end;
    end;
    repFactor=repFactor*length(factorLevels{theFactor});
end;

for i=1:totalFactorLevels
    tableDataW{i,1}=EPmain.anova.data.columnNames{EPmain.anova.data.leftColumn+i-1};
    tableDataW{i,2+numFactors}=0;
    tableDataW{i,3+numFactors}=0;
    tableDataW{i,4+numFactors}=0;
    tableDataW{i,5+numFactors}=0;
    tableDataW{i,6+numFactors}=0;
end;

tableNamesW{1}='Cells';
for i=1:numFactors
    tableNamesW{i+1}=factorNames{i};
end;

tableNamesW{2+numFactors}='Con1';
tableNamesW{3+numFactors}='Con2';
tableNamesW{4+numFactors}='Con3';
tableNamesW{5+numFactors}='Con4';
tableNamesW{6+numFactors}='Con5';

columnEditableW =  [repmat(false,1,1+numFactors) repmat(true,1,5)];

columnWidthW{1}=135;
for i=2:length(tableNamesW)
    columnWidthW{i}=50;
end;

try
    EPContrast.handles.withinTable = uitable('Data',tableDataW,'ColumnName',tableNamesW,'FontSize',EPmain.fontsize,...
        'ColumnWidth',columnWidthW,...
        'ColumnEditable', columnEditableW, 'Position',[100 windowHeight-400 569 300]);
catch
    uicontrol('Style','text','HorizontalAlignment','left','String', 'This function does not work with this version of Matlab.','FontSize',EPmain.fontsize,...
        'Position',[100 windowHeight-400 569 300]);
end;

%between-group levels table

betweenFactors=[];
totBetweenLevels=1;
betweenFactorLevels=cell(0);
for theFactor=1:6
    if ~isempty(EPmain.anova.data.betweenName{theFactor})
        theLevels={EPmain.anova.data.betweenLvl{:,EPmain.anova.data.between(theFactor)-EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1}};
        betweenFactors(end+1)=EPmain.anova.data.between(theFactor);
        betweenFactorLevels{end+1}=unique(theLevels);
        totBetweenLevels=totBetweenLevels*length(unique(theLevels));
    end;
end;

if totBetweenLevels>1
    repFactor=1;
    repCount=1;
    lvlCount=1;
    numFactors=length(betweenFactors);
    for theFactor=1:length(betweenFactorLevels)
        for theRow=1:totBetweenLevels
            tableDataB{theRow,theFactor}=betweenFactorLevels{theFactor}{lvlCount};
            repCount=repCount+1;
            if repCount > repFactor
                repCount=1;
                lvlCount=lvlCount+1;
                if lvlCount > length(betweenFactorLevels{theFactor})
                    lvlCount=1;
                end;
            end;
        end;
        repFactor=repFactor*length(betweenFactorLevels{theFactor});
    end;
    
    for i=1:numFactors
        tableNamesB{i}=EPmain.anova.data.betweenName{i};
    end;
    
    for i=1:totBetweenLevels
        tableDataB{i,1+numFactors}=0;
        tableDataB{i,2+numFactors}=0;
        tableDataB{i,3+numFactors}=0;
        tableDataB{i,4+numFactors}=0;
        tableDataB{i,5+numFactors}=0;
    end;
    
    tableNamesB{1+numFactors}='Con1';
    tableNamesB{2+numFactors}='Con2';
    tableNamesB{3+numFactors}='Con3';
    tableNamesB{4+numFactors}='Con4';
    tableNamesB{5+numFactors}='Con5';
    
    columnEditableB =  [repmat(false,1,numFactors) repmat(true,1,5)];
    
    for i=1:length(tableNamesB)
        columnWidthB{i}=50;
    end;
    
    try
        EPContrast.handles.betweenTable = uitable('Data',tableDataB,'ColumnName',tableNamesB,'FontSize',EPmain.fontsize,...
            'ColumnWidth',columnWidthB,...
            'ColumnEditable', columnEditableB, 'Position',[100 windowHeight-720 569 300]);
    catch
        uicontrol('Style','text','HorizontalAlignment','left','String', 'This function does not work with this version of Matlab.','FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-720 569 300]);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ANOVAcontrast(src,eventdata) %run the ANOVA contrasts
global EPmain EPContrast

set(EPContrast.handles.run,'enable','off');
set(EPContrast.handles.done,'enable','off');

WithinTableData=get(EPContrast.handles.withinTable,'Data');
BetweenTableData=get(EPContrast.handles.betweenTable,'Data');

numContrasts=0;
for i=1:5
    withinCon=[WithinTableData{:,size(WithinTableData,2)-5+i}];
    betweenCon=[BetweenTableData{:,size(BetweenTableData,2)-5+i}];
    if any(withinCon) || any(betweenCon)
        numContrasts=numContrasts+1;
        if any(withinCon)
            U{i}=withinCon';
            if (sum(withinCon) ~= 0) && (length(find(withinCon ==1)) ~= length(withinCon))
                msg{1}='Error: within contrast does not sum to zero.';
                [msg]=ep_errorMsg(msg);
                set(EPContrast.handles.run,'enable','on');
                set(EPContrast.handles.done,'enable','on');
                return
            end;
        else
            U{i}=repmat(1,size(WithinTableData,1),1);
        end;
        if any(betweenCon)
            C{i}=betweenCon;
            if (sum(betweenCon) ~= 0) && (length(find(betweenCon ==1)) ~= length(betweenCon))
                msg{1}='Error: between contrast does not sum to zero.';
                [msg]=ep_errorMsg(msg);
                set(EPContrast.handles.run,'enable','on');
                set(EPContrast.handles.done,'enable','on');
                return
            end;
        else
            C{i}=repmat(1,1,size(BetweenTableData,1));
        end;
    end;
end;

[outFileName, pathname] = uiputfile('*.html','ANOVA output:');
if outFileName == 0
    return %user hit cancel on file requestor
end;
outFileName=[pathname outFileName];
if exist(outFileName,'file')
    delete(outFileName); %user must have clicked "yes" to whether to replace existing file
end;
[pathstr, fileName, ext] = fileparts(outFileName);
if ~strcmp(ext,'.html')
    ext='.html';
end;

outfid=fopen([pathstr filesep fileName ext],'w');

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
    set(EPContrast.handles.run,'enable','on');
    set(EPContrast.handles.done,'enable','on');
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

numComps=EPmain.anova.numComps;
if numComps==0
    numComps=1;
end;
alpha.uncorrected=.05; %threshold for declaring statistical significance.
alpha.corrected=alpha.uncorrected/numComps; %with Bonferroni correction

fprintf(outfid,'<font face="Courier">');

for theFile=1:length(ANOVAfiles)
    [pathstr, fileName, ext] = fileparts(ANOVAfiles{theFile});
    infid=fopen([pathname fileName ext],'r');
    if infid == -1
        msg{1}=['There was an error trying to open ' fileName '.'];
        [msg]=ep_errorMsg(msg);
        set(EPContrast.handles.run,'enable','on');
        set(EPContrast.handles.done,'enable','on');
        return
    end
    
    fprintf(outfid,'%s </BR>',fileName);
    
    %load in the ANOVA data
    ANOVA = textscan(infid, '%s','delimiter','\n');
    ANOVAdata.name=deblank(ANOVA{1}{1});
    ANOVAdata.window=deblank(ANOVA{1}{2});
    ANOVAdata.measure=deblank(ANOVA{1}{3});
    ANOVAdata.changrp=deblank(ANOVA{1}{4});
    ANOVAdata.factor=deblank(ANOVA{1}{5});
    remain = ANOVA{1}{6};
    numCols=length(find(strfind(remain,sprintf('\t'))))+1;
    if strcmp(remain(end),sprintf('\t'))
        numCols=numCols-1;
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
            [theNum, remain] = strtok(remain,sprintf('\t'));
            ANOVAdata.data(i,k)=str2num(theNum);
        end
        for k = 1:numSpecs
            [ANOVAdata.betweenLvl{i,k}, remain] = strtok(remain,sprintf('\t'));
        end
    end;
    
    %organize the between groups
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
                set(EPContrast.handles.run,'enable','on');
                set(EPContrast.handles.done,'enable','on');
                return
            end;
            
            if length(unique({ANOVAdata.betweenLvl{:,EPmain.anova.data.between(theFactor)-numDataCols}})) ==1
                msg{1}=['The between group variable ' EPmain.anova.data.betweenName{theFactor} ' in ' fileName ' has only a single level.'];
                [msg]=ep_errorMsg(msg);
                set(EPContrast.handles.run,'enable','on');
                set(EPContrast.handles.done,'enable','on');
                return
            end;
            
            betweenFactors(end+1)=EPmain.anova.data.between(theFactor);
            totBetweenLevels=totBetweenLevels*length(unique(theLevels));
        end;
    end;
    
    if isempty(betweenFactors)
        subjects=size(ANOVAdata.data,1);
        for theContrast=1:numContrasts
            C{theContrast}=1;
        end;
    else
        [sortedBetween,index] = sortrows(ANOVAdata.betweenLvl(:,betweenFactors-numDataCols));
        ANOVAdata.data(index,:)=ANOVAdata.data;
        
        betweenCombos=cellstr(cell2mat(sortedBetween));
        if length(unique(betweenCombos)) ~= totBetweenLevels
            msg{1}='The between group factors are incompletely crossed.  Nested designs are not currently supported.';
            [msg]=ep_errorMsg(msg);
            set(EPContrast.handles.run,'enable','on');
            set(EPContrast.handles.done,'enable','on');
            return
        end;
        
        uniqueCombos=unique(betweenCombos);
        subjects=zeros(totBetweenLevels,1);
        for i=1:totBetweenLevels
            subjects(i)=length(find(strcmp(uniqueCombos{i},betweenCombos)));
        end;
        if any(subjects < 2)
            msg{1}='Each between group cell needs at least two subjects to be able to estimate error variance.';
            [msg]=ep_errorMsg(msg);
            set(EPContrast.handles.run,'enable','on');
            set(EPContrast.handles.done,'enable','on');
            return
        end;
    end;
    
    Y=ANOVAdata.data;
    for theContrast=1:numContrasts
        fprintf(outfid,'Contrast: %d </BR>',theContrast);
        fprintf(outfid,'Within: %s </BR>',num2str(U{theContrast}'));
        fprintf(outfid,'Between: %s </BR>',num2str(C{theContrast}));
        resultsTot=cell(EPmain.preferences.anova.reps,1);
        pList=zeros(EPmain.preferences.anova.reps,1);
        for iRep=1:EPmain.preferences.anova.reps
            [MUHAT, SIGMA, RESULTS]=ep_WJGLMml(Y, subjects, C{theContrast}, U{theContrast}, 1, EPmain.preferences.anova.trimming, 1, EPmain.preferences.anova.bootstrap, EPmain.preferences.anova.seed*iRep, EPmain.preferences.anova.missing, 0, alpha, 1, 0, 0);
            resultsTot{iRep}=RESULTS;
            pList(iRep)=RESULTS(4);
        end;
        [A B]=sort(pList);
        RESULTS=resultsTot{find(B==ceil(length(pList)/2))};
        if isempty(RESULTS)
            return
        end;
        alpha.pvar=std(pList)*2;
        ep_printADFResults(RESULTS, alpha, 0, outfid, EPmain.preferences.anova.bootstrap);
    end;
end;

set(EPContrast.handles.run,'enable','on');
set(EPContrast.handles.done,'enable','on');


