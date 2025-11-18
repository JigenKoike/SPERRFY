% main_additional_analysis.m


%% PI Similarity of Holdout Test
holdoutPISimilarityList = mouseConnectomeDataAnalysis.OverallModelAndResults.calculateHoldoutPISimilarityList();

%% Calculate Cosine Similarity between Wiring PIs
wiringPIdata = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
PIMatSource = wiringPIdata.WiringPISource(:,1:dimPI);
PIMatTarget = wiringPIdata.WiringPITarget(:,1:dimPI);

PISim_ss = corr(PIMatSource,PIMatSource);
PISim_tt = corr(PIMatTarget,PIMatTarget);
PISim_st = corr(PIMatSource,PIMatTarget);


%% Calculate Spatial Autocorrelation
regionDistanceMedian = SpatialAutocorrelationCalculator.gethFromMedian(brainRegionDistanceMatrix);
hList_exponential = SpatialAutocorrelationCalculator.makehList_exponentialScale(regionDistanceMedian,2,-7,1);
hList_quantile = SpatialAutocorrelationCalculator.makehList_quantile(brainRegionDistanceMatrix,0.01);

% PI
hList = hList_exponential;
spatialAutocorrelationMatrix_PISource = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    wiringPIdata.WiringPISource(:,:),brainRegionDistanceMatrix,hList);
spatialAutocorrelationMatrix_PITarget = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    wiringPIdata.WiringPITarget(:,:),brainRegionDistanceMatrix,hList);

% genes
spatialAutocorrelationMatrix_geneExpression = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    geneExprLevels.ExpressionMatrix,brainRegionDistanceMatrix,hList);

% PC
spatialAutocorrelationMatrix_PC = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    geneExprLevels.PCAMatrix,brainRegionDistanceMatrix,hList);


%% Reconstruction by distance matrix
normalizedDistanceMatrix = crossInfo.DistanceMatrix / max(crossInfo.DistanceMatrix,[],'all');
distanceReconstResults = ReconstructionResults(reconstructionParameters);
distanceReconstResults = distanceReconstResults.compute(mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData, ...
    connectomeMatrix.Matrix,crossInfo.DomainDefineMatrix,"ArbitralPIDistanceMatrix",normalizedDistanceMatrix,"StoreMatrixTag",true);


%% Make Reciprocal & Unidirectional Index
reciprocalFeatureMasks = ConnectomeFeatureAnalyzer.makeReciprocalUnidirectionalMasks(connectomeMatrix,crossInfo);


%% Distance Preserved with 50um or 200 um Bin Width
rng(4)
nullBatchRunner_3_50um = NullModelAnalysisBatchRunner(nullGenerator_3,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_3_50um = nullBatchRunner_3_50um.run(mouseConnectomeDataAnalysis,"DistanceBinWidth",50);
nullBatchRunner_4_50um = NullModelAnalysisBatchRunner(nullGenerator_4,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_4_50um = nullBatchRunner_4_50um.run(mouseConnectomeDataAnalysis,"DistanceBinWidth",50);

nullBatchRunner_3_200um = NullModelAnalysisBatchRunner(nullGenerator_3,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_3_200um = nullBatchRunner_3_200um.run(mouseConnectomeDataAnalysis,"DistanceBinWidth",200);
nullBatchRunner_4_200um = NullModelAnalysisBatchRunner(nullGenerator_4,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_4_200um = nullBatchRunner_4_200um.run(mouseConnectomeDataAnalysis,"DistanceBinWidth",200);


%% Gene Surrogate Test
main_gene_surrogate
