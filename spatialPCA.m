function [PCAResults] = spatialPCA(data,channelLocations,rotation,minimumVariance)

    % original PCA code from Joe Dien, shell code developed by O. Krigolson
    % supported rotations are 'IMAX' for ICA, 'PMAX' for Promax, and 'VMAX' for varimax rotations
    % assumed data format is a 4D matrix, channels x time x conditions x
    % participants
    % a channel location file is also needed
    % finally, you must specify the minimum amount of variance to show

    [spatialPCAData] = prepareSpatialData(data);

    [PCAResults] = ep_doPCA('asis',rotation,'3','COV',10,spatialPCAData,'N','N');

    [PCAResults] = reconstructSpatialPCAData(data,PCAResults);

    plotSpatialLoadings(PCAResults,channelLocations,minimumVariance);

end