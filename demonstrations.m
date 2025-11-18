%% demonstration.m
%% Startup
clear; clc;
startup_SPERRFY;  % path setting

%% Data Preparation
main_import_processed_data;  % import processed data (.mat)
main_data_preparation; % make basic objects for analysis

%%
% Figure 3a
f3a = FigMaker.imageConnectionMatrix_MRColorLabels(connectomeMatrix,crossInfo);

%% Parameter Setting
main_parameter_setting;

%% Mouse connectome data analysis (Figure 3~5)
rng(1);
% --perform CCA and reconstruction--
mouseConnectomeDataAnalysis = ConnectomeAnalysisRunner(connectomeMatrix, geneExprLevels,crossInfo, ...
    pairDataGenerationOptions,reconstructionParameters,holdoutParameters_main,"FactoryDescription","Main Data Analysis");
mouseConnectomeDataAnalysis = mouseConnectomeDataAnalysis.runAnalysis(ReconstructionStoreTag=true);

%%
% Fig. 3b
f3b = FigMaker.plotAllCorrCoeff_FullAndHoldoutTest(mouseConnectomeDataAnalysis);
% Fig. 3c
f3c = FigMaker.scatterCorrelations(mouseConnectomeDataAnalysis);
% Fig. 3d
f3d = FigMaker.imageWiringPIs(mouseConnectomeDataAnalysis,brainSpace,brainInfo);

%% Perform gene analysis
relatedGeneAnalysisResults = mouseConnectomeDataAnalysis.performRelatedGeneAnalysis();
%% 
% Figure 4a
f4a = FigMaker.histogramsGeneDistributionSimilarity(relatedGeneAnalysisResults);
% Figure 4b
f4b = FigMaker.image3DSimilarGeneDistributions(relatedGeneAnalysisResults,brainSpace,brainInfo);

%% Reconstruction by distance matrix
normalizedDistanceMatrix = crossInfo.DistanceMatrix / max(crossInfo.DistanceMatrix,[],'all');
distanceReconstResults = ReconstructionResults(reconstructionParameters);
distanceReconstResults = distanceReconstResults.compute(mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData, ...
    connectomeMatrix.Matrix,crossInfo.DomainDefineMatrix,"ArbitralPIDistanceMatrix",normalizedDistanceMatrix,"StoreMatrixTag",true);
%% 
% figure 5a
f5a = FigMaker.imageWiringPIDiffForReconstruction(mouseConnectomeDataAnalysis); % need colorbar: 
% figure 5b
f5b = FigMaker.plotROCHoldoutAndDistReconst(mouseConnectomeDataAnalysis,distanceReconstResults);
% figure 5c
f5c = FigMaker.swarmHoldoutAUC(mouseConnectomeDataAnalysis,distanceReconstResults);
% figure 5d
f5d = FigMaker.scatterPIDiffvsRegionDistance(mouseConnectomeDataAnalysis);

%% Perform null connectome analysis (Figure 6)
% --make null connectome model--
nullGenerator_global = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"RandomGeneration",0, ...
    "Description","RandomGeneration, Golobally");
nullGenerator_local = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"RandomGeneration",1, ...
    "Description","RandomGeneration, Locally");

% --make sample matrix--
rng(2)
randomConnectomeSample_global = nullGenerator_global.generate;
randomConnectomeSample_local = nullGenerator_local.generate;

% !!!For Demonstration: Small Number of Test Times!!!
numRandomConnectome = 25; % number of random connectome data per each null model
                  % original analysis: numNullTest = 1000
% --run test--
nullBatchRunner_global = NullModelAnalysisBatchRunner(nullGenerator_global,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_global = nullBatchRunner_global.run(mouseConnectomeDataAnalysis);
nullBatchRunner_local = NullModelAnalysisBatchRunner(nullGenerator_local,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_local = nullBatchRunner_local.run(mouseConnectomeDataAnalysis);


%%
% figure 6a
f6a(1) = FigMaker.imageConnectionMatrix_MRColorLabels(ConnectionMatrix(randomConnectomeSample_global),crossInfo);
f6a(2) = FigMaker.imageConnectionMatrix_MRColorLabels(ConnectionMatrix(randomConnectomeSample_local),crossInfo);
% figure 6b
f6b = FigMaker.plotRandomTestCorrCoef(mouseConnectomeDataAnalysis,nullBatchRunner_global,nullBatchRunner_local);
% figure 6c
f6c = FigMaker.swarmRandomTestCorrCoef(mouseConnectomeDataAnalysis,nullBatchRunner_global,nullBatchRunner_local,"SwarmSize",12);
% figure 6d
f6d = FigMaker.swarmRandomTestAUC(mouseConnectomeDataAnalysis,nullBatchRunner_global,nullBatchRunner_local,"SwarmSize",12);
% figure 6e
f6e = FigMaker.swarmRandomTestSimilarity(nullBatchRunner_global,nullBatchRunner_local,"SwarmSize",12);





