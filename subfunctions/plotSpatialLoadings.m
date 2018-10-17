function plotSpatialLoadings(PCAResults,chanlocs,percentToPlot)

    % plots spatial factor loadings from a spatial PCA on a topography

    spatialLoadings = PCAResults.FacPat;
    spatialVariance = PCAResults.facVar*100;

    for toposToPlot = 1:size(spatialLoadings,2)

        if spatialVariance(toposToPlot) < percentToPlot
            toposToPlot = toposToPlot - 1;
            break
        end

    end

    figure;
    plotPosition = 1;
    
    if toposToPlot < 5
        prows = 2;
        pcolumns = 2;
    else
        prows = 2;
        pcolumns = 4;
    end

    while 1

        subplot(prows,pcolumns,plotPosition);
        topoplot(spatialLoadings(:,plotPosition),chanlocs,'shrink','on','plotrad',0.6);

        spatialVar = num2str(spatialVariance(plotPosition));
        currentFactor = num2str(plotPosition);
        factorText = strcat('Factor: ', currentFactor); 
        spatialVar = ['Variance Accounted For: ' spatialVar ' %'];
        title({factorText;spatialVar});

        plotPosition = plotPosition + 1;
        if plotPosition > toposToPlot
            break
        end

    end
    
end