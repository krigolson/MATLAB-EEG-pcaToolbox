function [temporalspatialData] = prepareTemporalSpatialData(data)

    % takes a channels x time x conditions x participants data set and makes it
    % channels x (time x conditions x participants) for spatial PCA

    numberOfConditions = size(data,2);
    numberOfSubjects = size(data,3);
    
    temporalspatialData = [];

    % put the data into a format suitable for spatial PCA
    for subjectCounter = 1:numberOfSubjects
        for conditionCounter = 1:numberOfConditions
            temporalspatialData = [temporalspatialData squeeze(data(:,conditionCounter,subjectCounter))];
        end
    end
    
    temporalspatialData = temporalspatialData';

end