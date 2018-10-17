function [spatialData] = prepareSpatialData(data)

    % takes a channels x time x conditions x participants data set and makes it
    % channels x (time x conditions x participants) for spatial PCA

    numberOfConditions = size(data,3);
    numberOfSubjects = size(data,4);
    
    spatialData = [];

    % put the data into a format suitable for spatial PCA
    for subjectCounter = 1:numberOfSubjects
        for conditionCounter = 1:numberOfConditions
            spatialData = [spatialData squeeze(data(:,:,conditionCounter,subjectCounter))];
        end
    end
    
    spatialData = spatialData';

end