function [spatialtemporalData] = prepareSpatialTemporalData(data)

    % takes a time x conditions x participants data set and makes it
    % time x (conditions x participants) for spatial PCA after temporal PCA

    numberOfConditions = size(data,2);
    numberOfSubjects = size(data,3);
    
    spatialtemporalData = [];

    % put the data into a format suitable for spatial PCA
    for subjectCounter = 1:numberOfSubjects
        for conditionCounter = 1:numberOfConditions
            spatialtemporalData = [spatialtemporalData squeeze(data(:,conditionCounter,subjectCounter))];
        end
    end
    
    spatialtemporalData = spatialtemporalData';

end