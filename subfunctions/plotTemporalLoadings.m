function plotTemporalLoadings(PCAResults,timeVector,percentToPlot);

    temporalLoadings = PCAResults.FacPat;
    temporalVariance = PCAResults.facVar*100;

    for componentsToPlot = 1:size(temporalLoadings,2)

        if temporalVariance(componentsToPlot) < percentToPlot
            componentsToPlot = componentsToPlot - 1;
            break
        end

    end

    figure;

    subplot(3,3,1);
    for counter = 1:componentsToPlot
        plot(timeVector,temporalLoadings(:,counter));
        hold on;
    end
    title('All Temporal Factors');
    hold off;
    
    for counter = 1:componentsToPlot
        
        subplot(3,3,counter+1);
        plot(timeVector,temporalLoadings(:,counter));
        temporalVar = num2str(temporalVariance(counter));
        currentFactor = num2str(counter);
        factorText = strcat('Factor: ', currentFactor); 
        spatialVar = ['Variance Accounted For: ' temporalVar ' %'];
        title({factorText;temporalVar});
        
    end

end