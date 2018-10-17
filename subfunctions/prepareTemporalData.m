function [temporalData] = prepareTemporalData(data)

    % takes a channels x time x conditions x participants data set and makes it
    % time x (channels x conditions x participants) for temporal PCA

    numberOfConditions = size(data,3);
    numberOfSubjects = size(data,4);
    
    temporalData = [];

    % put the data into a format suitable for spatial PCA
    for subjectCounter = 1:numberOfSubjects
        for conditionCounter = 1:numberOfConditions
            temporalData = [temporalData; squeeze(data(:,:,conditionCounter,subjectCounter))];
        end
    end

end