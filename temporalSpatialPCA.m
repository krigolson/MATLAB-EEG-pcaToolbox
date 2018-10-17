function [PCAResults TSPCAResults] = temporalSpatialPCA(data,channelLocations,timeVector,trotation,srotation,minimumVariance)

    % original PCA code from Joe Dien, shell code developed by O. Krigolson
    % supported rotations are 'IMAX' for ICA, 'PMAX' for Promax, and 'VMAX' for varimax rotations
    % assumed data format is a 4D matrix, channels x time x conditions x
    % participants
    % a time vector for plotting is also needed and channel locations
    % finally, you must specify the minimum amount of variance to show

    [temporalPCAData] = prepareTemporalData(data);

    [PCAResults] = ep_doPCA('asis',trotation,'3','COV',10,temporalPCAData,'N','N');

    [PCAResults] = reconstructTemporalPCAData(data,PCAResults);

    plotTemporalLoadings(PCAResults,timeVector,minimumVariance);
    
    selectedTemporalFactor = input('Which temporal factor do you want to run a spatial PCA on ? : ');
    
    spatialData = squeeze(PCAResults.virtualData(:,selectedTemporalFactor,:,:));
    
    [temporalspatialData] = prepareTemporalSpatialData(spatialData);
    
    [TSPCAResults] = ep_doPCA('asis',srotation,'3','COV',10,temporalspatialData,'N','N');
    
    [TSPCAResults] = reconstructTemporalSpatialPCAData(spatialData,TSPCAResults);
    
    plotSpatialLoadings(TSPCAResults,channelLocations,minimumVariance);
    
    selectedSpatialFactor = input('Which spatial factor do you want scores for ? : ');
    
    scoreData = squeeze(TSPCAResults.TSPCAScores(selectedSpatialFactor,:,:));
    
    scoreData = scoreData';
    
    TSPCAResults.scoreData = scoreData;
    
    for counter = 1:size(scoreData,2)
        CIs = makeCIs(scoreData(:,counter));
        yErr(counter) = CIs(4);
        y(counter) = CIs(1);
    end
    
    figure;
    barwitherr(yErr,y);
    title(['TSPCA Scores for the ' num2str(selectedTemporalFactor) ' temporal factor and the ' num2str(selectedSpatialFactor) ' spatial factor']);

end