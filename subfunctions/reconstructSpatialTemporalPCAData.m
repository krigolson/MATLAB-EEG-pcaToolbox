function [factorResults] = reconstructSpatialTemporalPCAData(data,factorResults)

    % takes the output of a temporal PCA after a spatial PCA and reconstructs the Factor Scores
    % into scores.

    numberOfConditions = size(data,2);
    numberOfSubjects = size(data,3);

    tempData = factorResults.FacScr;
    tempData = tempData';
    scoreData = [];
    
    counter = 1;

    for subjectCounter = 1:numberOfSubjects 
        for conditionCounter = 1:numberOfConditions
            scoreData(:,conditionCounter,subjectCounter) = tempData(:,counter);
            counter = counter + 1;
        end
    end
    
    factorResults.STPCAScores = scoreData;

end