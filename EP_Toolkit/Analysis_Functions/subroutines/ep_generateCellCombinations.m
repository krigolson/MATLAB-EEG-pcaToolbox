function [cellNames, cellCombinations] = ep_generateCellCombinations(factorNames, levelNames);
% ep_generateCellCombinations - [cellNames, cellCombinations] = ep_generateCellCombinations(factorNames, levelNames);
%Inputs
%  factorNames: The names of the within factors.  Each should be three letters long.
%               The order should correspond to the structure of the data with the first varying the slowest.
%               Set to empty set if no within factors.
%  levelNames : The names of the levels within each within factor.  One row for each factor.  Each level name should be a single letter.
%               The order should correspond to the order of the levels in the data.
%
%Outputs
%  cellNames    : cell array with the names of the experimental cells
%  cellCombinations : cell array with the combinations of experimental cells
%
%  Takes the analysis structure of the experiment and generates the combinations of cells that would be of interest
%  for examining the nature of main effects and interactions.  For within-factors only.
%  In order of factors, the largest number is the one that changes most quickly as one goes across the columns

%History
%  by Joseph Dien (2/24/09)
%  jdien07@mac.com
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

cellNames=[];
cellCombinations=[];

if ~iscell(levelNames)
    msg{1}='Error: The levelNames array needs to be a cell array (i.e., delimited by curly brackets).';
    [msg]=ep_errorMsg(msg);
    return
end

facNum = size(factorNames,1);

%how many levels in each factor?
facSize = zeros(facNum,1);
for i = 1:facNum
    facSize(i)=length(levelNames{i});
end;

%which columns are associated with each factor?
totalLevels=sum(facSize);
levelColumns={};
totalColumns=prod(facSize);

levelColumnsCounter=1;
for theFactor = 1:facNum
    for theLevel = 1:facSize(theFactor);
        clumpSize=totalColumns/prod(facSize(1:theFactor)); %factor levels repeat in clumps of this size
        clumpNum=(totalColumns/facSize(theFactor))/clumpSize; %number of clumps
        columnIndex=repmat([zeros(1,(theLevel-1)*clumpSize) ones(1,clumpSize) zeros(1,(facSize(theFactor)-theLevel)*clumpSize)],1,clumpNum);
        levelColumns{theFactor,theLevel}=find(columnIndex);
    end;
end;

factorCombinations=cell(facNum,1);
for theFactorLvl=1:facNum %level of interactions, starting with main effects
    factorCombinations{theFactorLvl}=nchoosek(1:facNum,theFactorLvl); %all the combinations of factors
    for theFactorComb=1:size(factorCombinations{theFactorLvl},1) %each combination of factors
        theComb=factorCombinations{theFactorLvl}(theFactorComb,:); %the combination of factors being worked on
        totalCombLevels=prod(facSize(theComb));
        combNames=cell(totalCombLevels,1);
        combCells=cell(totalCombLevels,1);
        for whichFactor = 1:length(theComb) %which of the combination is being added to the list
            inClumpCounter=1;
            lvlCounter=1;
            clumpSize=totalCombLevels/prod(facSize(theComb(1:whichFactor))); %level name repeats in clumps of this size
            for combLevel=1:totalCombLevels
                %names of the combinations of cells
                theName=levelNames{theComb(whichFactor)}(lvlCounter);
                combNames{combLevel}=[combNames{combLevel} theName];
                %columns going into each combination
                if isempty(combCells{combLevel,1})
                    combCells{combLevel,1}=levelColumns{theComb(whichFactor),lvlCounter};
                else
                    combCells{combLevel,1}= intersect(levelColumns{theComb(whichFactor),lvlCounter}, combCells{combLevel,1});
                end; 
                inClumpCounter=inClumpCounter+1;
                if inClumpCounter > clumpSize
                    inClumpCounter=1;
                    lvlCounter=lvlCounter+1;
                    if lvlCounter > facSize(theComb(whichFactor))
                        lvlCounter=1;
                    end;
                end;
            end;
        end;
        cellNames=[cellNames; combNames];
        cellCombinations=[cellCombinations; combCells];
    end;
end;












