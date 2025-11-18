% sandbox_for_holdout_modulation

%% Startup
clear; clc;
startup_SPERRFY;  % path setting

%% Load data
main_import_processed_data;  % import data (.mat)
main_data_preparation;

%% Parameter setting
main_parameter_setting;
% small number of null connectome test
numNullTest = 10;

%% Mouse connectome data analysis
rng(1);
mouseConnectomeDataAnalysis = ConnectomeAnalysisRunner(connectomeMatrix, geneExprLevels,crossInfo, ...
    pairDataGenerationOptions,reconstructionParameters,holdoutParameters_main,"FactoryDescription","Main Data Analysis");
mouseConnectomeDataAnalysis = mouseConnectomeDataAnalysis.runAnalysis(ReconstructionStoreTag=true);

%% PI Similarity
holdoutPISimilarityList = mouseConnectomeDataAnalysis.OverallModelAndResults.calculateHoldoutPISimilarityList();


%%
dataShow = holdoutPISimilarityList;
dataName = "PI SImilarity Score";
yLimVec = [0,1];
%%
dataShow = mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutSummary.TestAllCorrelations;
dataName = "Correlation Coefficients";
yLimVec = [min(dataShow,[],'all')-0.01,1];

%%
dimShow = 50;
figure;
mu = mean(dataShow,1);
sigma = std(dataShow,0,1);
x = 1:size(dataShow,2);
errorbar(x,mu,sigma,'o-','LineWidth', 1.5, 'CapSize',6);
xlabel("Correlation Componet (rank)")
ylabel(dataName)
xlim([0,dimShow+1])
ylim(yLimVec)
title(strcat(dataName," (Holdout)"))

%%
dimShow = 50;
figure;
boxplot(dataShow);
xlabel("Correlation Componet (rank)")
ylabel(dataName)
xlim([0,dimShow+1])
ylim(yLimVec)
title(strcat(dataName," (Holdout)"))

%%
dimShow = 50;
figure;
sz = 6;
c = colorList13(9,:);
x = ones(height(dataShow),1) * (1:dimShow);
h = swarmchart(x,dataShow(:,1:dimShow),sz,c,'filled');
xlabel("Correlation Componet (rank)")
ylabel(dataName)
xlim([0,dimShow+1])
ylim(yLimVec)
title(strcat(dataName," (Holdout)"))

%%
hold on
x = 1:dimShow;
y = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.CCAResultsData.r(1,:);
scatter(x,y,sz*2,'red','filled')


%% Reconstruction Matrix
holdoutPIDistanceList = mouseConnectomeDataAnalysis.OverallModelAndResults.calculateHoldoutPIDistanceList();

%%
figure;
t = tiledlayout(1,3);
nexttile
PIDiffData = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.PIDistanceMatrix;
imagesc(abs(PIDiffData - mean(holdoutPIDistanceList,3)));
title("|Full Data - Holdout Mean|")
axis square
colorbar
nexttile
imagesc(mean(holdoutPIDistanceList,3));
title("Holdout Mean")
axis square
colorbar
nexttile
imagesc(std(holdoutPIDistanceList,0,3));
title("Holdout Std")
axis square
colorbar

%% Make ROC of Holdout Test
f = figure;
%c = [0 0.4470 0.7410 0.05];
c = [0 0 0 0.05];
for i = 1:holdoutParameters_main.NumSplits
    h = gca();
    holdoutReconstResults = mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutUnits(i).TestMetrics.TestPseudeReconstructionResults;
    holdoutReconstResults.plotROC("ParentAxes",h,"TitleOff",true,"LineWidth",0.1,"PlotColor",c);
    hold on
end
[meanFPR,meanTPR] = mouseConnectomeDataAnalysis.OverallModelAndResults.getHoldoutROCmean("TestOrTrain","Test");
meanAUC = mean(mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC,2);
plot(meanFPR,meanTPR,"Color",[0 0.4470 0.7410],"LineStyle","-","LineWidth",2)
titleStr = strcat("ROC plot of holdout reconstruction test (mean AUC = ",num2str(round(meanAUC,3))," )");
title(gca,titleStr)
hold on
distanceReconstResults.plotROC("ParentAxes",h);



%%