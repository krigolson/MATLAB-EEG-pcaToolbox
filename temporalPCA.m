function [PCAResults] = temporalPCA(data,timeVector,rotation,minimumVariance);

    % original PCA code from Joe Dien, shell code developed by O. Krigolson
    % supported rotations are 'IMAX' for ICA, 'PMAX' for Promax, and 'VMAX' for varimax rotations
    % assumed data format is a 4D matrix, channels x time x conditions x
    % participants
    % a time vector for plotting is also needed
    % finally, you must specify the minimum amount of variance to show

    [temporalPCAData] = prepareTemporalData(data);

    [PCAResults] = ep_doPCA('asis',rotation,'3','COV',10,temporalPCAData,'N','N');

    [PCAResults] = reconstructTemporalPCAData(data,PCAResults);

    plotTemporalLoadings(PCAResults,timeVector,minimumVariance);

end