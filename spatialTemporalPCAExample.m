clear all;
close all;
clc;

% demo to show how to run a spatial-temporal PCA. Loads a provided sample
% data set of 250 people who did a standard visual oddball task. A quick
% ERP analysis is provided for comparison with the PCA. Requires EEGLAB as
% well as my Github files. The data is organized as channels (30) x time x
% conditions (2) x participants (250). Use IMAX rotation for spatial PCA
% and PMAX for temporal PCA

% note, you can accomplish the full STPCA with one call but I have broken
% it out here so you can see the steps
% [PCAResults TSPCAResults] = spatialTemporalPCA(data,chanlocs,timeVector,'PMAX','IMAX',5)

% select the factors for the spatial and temporal PCAs. Note, you need to
% run the analysis to actually know which one to look at so use 1 and 1 the
% first time and then rerun as needed
% select the first spatial factor
selectedSpatialFactor = 1;
% select the first temporal factor
selectedTemporalFactor = 1;

% load the data
load('sampleData');
% load a channel locations file
load('locs');

% simple ERP analysis with plot
theChannel = 12; % select Pz for classic P300B
dwData = squeeze(data(:,:,1,:) - data(:,:,2,:));
meanData = mean(data,4);
meanDWData = squeeze(mean(dwData,3));
timeVector = -200:2:800;
subplot(4,3,1);
plot(timeVector,meanData(theChannel,:,1));
hold on;
plot(timeVector,meanData(theChannel,:,2));
title('Conditional Waveforms');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
subplot(4,3,2);
plot(timeVector,meanDWData(theChannel,:));
title('Difference Waveforms');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
subplot(4,3,3);
topoData = squeeze(meanDWData(:,290));
doPlot2DTopo(topoData,chanlocs);
title('Peak Topography');
xlabel('Time (ms)');
ylabel('Voltage (uV)');

% prepare data for spatial PCA
[spatialPCAData] = prepareSpatialData(data);

% run a spatial PCA first, use IMAX rotation as per Joe Dien's
% recommendations
PCAResults = [];
[PCAResults] = ep_doPCA('asis','IMAX','3','COV',10,spatialPCAData,'K','N');

% reconstruct the output to get virtual ERPs
[PCAResults] = reconstructSpatialPCAData(data,PCAResults);

% plot the first three spatial factors
for plotCounter = 4:6
    subplot(4,3,plotCounter);
    topoplot(PCAResults.FacPat(:,(plotCounter-3)),chanlocs,'shrink','on','plotrad',0.6);
    spatialVar = num2str(PCAResults.facVar(plotCounter-3)*100);
    currentFactor = num2str(plotCounter-3);
    factorText = strcat('Factor: ', num2str(plotCounter-3)); 
    spatialVar = ['Variance Accounted For: ' spatialVar ' %'];
    title({factorText;spatialVar});
    set(gcf,'color','w');
    axis([-0.6 0.6 -0.6 0.6]);
end

% arrange the virtual erps for temporal PCA 
temporalData = squeeze(PCAResults.virtualERPs(selectedSpatialFactor,:,:,:));
[spatialTemporalData] = prepareSpatialTemporalData(temporalData);

% to the PCA on the virtual erps, use a VMAX rotation as per Joe Dien's
% recommendations
[STPCAResults] = ep_doPCA('asis','PMAX','3','COV',10,spatialTemporalData,'K','N');

% recontrust the data post PCA
[STPCAResults] = reconstructSpatialTemporalPCAData(temporalData,STPCAResults);

% plot the first three temporal factors
for plotCounter = 4:6
    subplot(4,3,plotCounter+3);
    plot(timeVector,STPCAResults.FacPat(:,(plotCounter-3)));
    temporalVar = num2str(STPCAResults.facVar(plotCounter-3)*100);
    currentFactor = num2str(plotCounter-3);
    factorText = strcat('Factor: ', num2str(plotCounter-3)); 
    temporalVar = ['Variance Accounted For: ' temporalVar ' %'];
    xlabel({factorText;temporalVar});
    set(gcf,'color','w');
end

% extract the virtual ERPs and plot them
virtualERPs = PCAResults.virtualERPs;
meanVirtualERPs = squeeze(mean(virtualERPs,4));
subplot(4,3,10);
plot(timeVector,meanVirtualERPs(1,:,1));
hold on;
plot(timeVector,meanVirtualERPs(1,:,2));
xlabel('Virtual ERPs');

% extract the STPCA scores
scoreData = squeeze(STPCAResults.STPCAScores(selectedTemporalFactor,:,:));
scoreData = scoreData';
STPCAResults.scoreData = scoreData;

subplot(4,3,11);
bar(mean(scoreData,1));
xlabel(['Mean STPCA Scores for the ' num2str(selectedSpatialFactor) ' spatial factor and the ' num2str(selectedTemporalFactor) ' temporal factor']);