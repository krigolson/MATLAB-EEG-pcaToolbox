function [factorResults] = reconstructSpatialPCAData(data,factorResults)

    % takes the output of a spatial PCA and reconstructs the Factor Scores
    % into Virtual ERP data (channels = components).

    numberOfConditions = size(data,3);
    numberOfSubjects = size(data,4);
    numberOfTimepoints = size(data,2);

    tempData = factorResults.FacScr;
    tempData = tempData';
    virtualData = [];
    
    x1 = 1;
    x2 = numberOfTimepoints;

    for subjectCounter = 1:numberOfSubjects
        for conditionCounter = 1:numberOfConditions   
            virtualData(:,:,conditionCounter,subjectCounter) = tempData(:,x1:x2);
            x1 = x2 + 1;
            x2 = x2 + numberOfTimepoints;
        end
    end
    
    factorResults.virtualERPs = virtualData;

end