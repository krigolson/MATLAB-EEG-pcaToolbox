function [factorResults] = reconstructTemporalPCAData(data,factorResults)

    % takes the output of a temporal PCA and reconstructs the Factor Scores
    % into temporal component data (time = components).

    numberOfConditions = size(data,3);
    numberOfSubjects = size(data,4);
    numberOfChannels = size(data,1);

    tempData = factorResults.FacScr;
    virtualData = [];
    
    x1 = 1;
    x2 = numberOfChannels;

    for subjectCounter = 1:numberOfSubjects
        for conditionCounter = 1:numberOfConditions   
            virtualData(:,:,conditionCounter,subjectCounter) = tempData(x1:x2,:);
            x1 = x2 + 1;
            x2 = x2 + numberOfChannels;
        end
    end
    
    factorResults.virtualData = virtualData;

end