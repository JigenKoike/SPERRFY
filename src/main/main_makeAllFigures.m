% main_makeAllFigures.m
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
colorList = orderedcolors("gem");
%% Directory
cd(strcat(projectRoot,"/results/figures"));
dt = datetime('now','Format','yyyyMMdd_HHmmss');
dtstr = string(dt);
mkdir(dtstr)
cd(dtstr)

%% Make Object
FM = FigMaker();
%% Figure 3 (Connectome Data & Wiring PI Results)
% Fig. 3a
f3a = FM.imageConnectionMatrix_MRColorLabels(connectomeMatrix,crossInfo);
% Fig. 3b
f3b = FM.plotAllCorrCoeff_FullAndHoldoutTest(mouseConnectomeDataAnalysis);
% Fig. 3c
f3c = FM.scatterCorrelations(mouseConnectomeDataAnalysis);
% Fig. 3d
f3d = FM.imageWiringPIs(mouseConnectomeDataAnalysis,brainSpace,brainInfo);
%%
% display similarity value
PISimBetweenSourceTarget = diag(PISim_st);
PISimBetweenSourceTarget = string(compose('%.3f', PISimBetweenSourceTarget)); 
%disp("Wiring PI Similarity between Source and Target")
%disp(PISimBetweenSourceTarget(1:5).')

%%
exportgraphics(f3a,'Figure3A.tiff');
exportgraphics(f3b,'Figure3B.pdf');
exportgraphics(f3c,'Figure3C.pdf');
exportgraphics(f3d,'Figure3D.pdf');

%% Figure S1 (Major Region)
% Fig. S1a
fs1a = FM.image3DMajorRegions(brainSpace,brainInfo);
% Fig. S1b
fs1b = FM.imageConnectomeProperties(connectomeMatrix,crossInfo);

%%
exportgraphics(fs1a,'FigureS1A.pdf');
exportgraphics(fs1b,'FigureS1B.pdf');

%% Figure S2 (Gene Expression Data)
% Fig. S2a
fs2a = FM.imageGeneExpressionMatrix(geneExprLevels,brainInfo);
% Fig. S2b
fs2b = FM.imagePCMatrix(geneExprLevels,brainInfo);
% Fig. S2c
fs2c = FM.plotPCACumContribution(geneExprLevels);

%%
exportgraphics(fs2a,'FigureS2A.tiff');
exportgraphics(fs2b,'FigureS2B.tiff');
exportgraphics(fs2c,'FigureS2C.pdf');


%% Figure S3 (Procedure of Holdout Test)
% data preparations
holdoutExample = mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutUnits(1);
crossInfoExample_train = CrossRegionInformation(brainInfo,brainInfo);
crossInfoExample_train.setDomainDefineMatrix(holdoutExample.SplitInfo.TrainMask);
crossInfoExample_test = CrossRegionInformation(brainInfo,brainInfo);
crossInfoExample_test.setDomainDefineMatrix(holdoutExample.SplitInfo.TestMask);
% Fig. S3
fs3(1) = FM.imageConnectionMatrix_MRColorLabels(connectomeMatrix,crossInfoExample_train,"DataType","DomainDefine");
fs3(2) = FM.imageConnectionMatrix_MRColorLabels(connectomeMatrix,crossInfoExample_test,"DataType","DomainDefine");
fs3(3) = FM.imageConnectionMatrix_MRColorLabels(connectomeMatrix,crossInfoExample_train);
fs3(4) = FM.imageConnectionMatrix_MRColorLabels(connectomeMatrix,crossInfoExample_test);

%%
exportgraphics(fs3(1),"FigureS3_1.tiff");
exportgraphics(fs3(2),"FiuregS3_2.tiff");
exportgraphics(fs3(3),"FigureS3_3.tiff");
exportgraphics(fs3(4),"FigureS3_4.tiff");

%% Figure S4 (Source-Target Correlation and Spatial Autocorrelation of Wiring PI)
% Fig. S4A
fs4a = FM.plotWiringPISpatialAutocorrelation(spatialAutocorrelationMatrix_PISource(:,1),spatialAutocorrelationMatrix_PITarget(:,1));
% Fig. S4B
fs4b = FM.heatmapWiringPIGradientCorrelation(PISim_ss,PISim_tt,PISim_st);

%%
exportgraphics(fs4a,'FigureS4A.pdf');
exportgraphics(fs4b,'FigureS4B.tiff');


%% Figure S5 (Wiring PI Values per Major Regions)
% Fig. S5A
fs5a = FM.swarmPIValuesByMR(mouseConnectomeDataAnalysis);
% Fig. S5B
fs5b = FM.scatterCorrelations_SpecificMRPairColor(mouseConnectomeDataAnalysis,pairDataGenerationOptions);

%%
exportgraphics(fs5a,'FigureS5A.pdf');
exportgraphics(fs5b,'FigureS5B.pdf');



%% Figure 4
% figure 4a
f4a = FM.histogramsGeneDistributionSimilarity(relatedGeneAnalysisResults);
% figure 4b
f4b = FM.image3DSimilarGeneDistributions(relatedGeneAnalysisResults,brainSpace,brainInfo);
%%
exportgraphics(f4a,'Figure4A.pdf');
exportgraphics(f4b,'Figure4B.pdf');


%% Figure S7 (Gene Cosine Similarity vs Spatial Autocorrelation)
% (Figure S6 is produced by g:Profiler)
% Figure S7
fs7 = FM.scatterGenePISimSpatialAutocorrelation(relatedGeneAnalysisResults,spatialAutocorrelationMatrix_geneExpression(:,1));
%%
exportgraphics(fs7,'FigureS7.pdf');

%% Figure 5
% figure 5a
f5a = FM.imageWiringPIDiffForReconstruction(mouseConnectomeDataAnalysis); % need colorbar: 
% figure 5b
f5b = FM.plotROCHoldoutAndDistReconst(mouseConnectomeDataAnalysis,distanceReconstResults);
% figure 5c
f5c = FM.swarmHoldoutAUC(mouseConnectomeDataAnalysis,distanceReconstResults);
% figure 5d
f5d = FM.scatterPIDiffvsRegionDistance(mouseConnectomeDataAnalysis);

%%
exportgraphics(f5a,'Figure5A.tiff');
exportgraphics(f5b,'Figure5B.pdf');
exportgraphics(f5c,'Figure5C.pdf');
exportgraphics(f5d,'Figure5D.pdf');

%% Figure S8
% Fig. S8a
thresholdList = [0.14,0.18,0.22];
fs8a = FM.imageReconstructedMatrixSamples(mouseConnectomeDataAnalysis,"ThreshList",thresholdList);
% Fig. S8b
fs8b = FM.plotReconstructionROC(mouseConnectomeDataAnalysis,"ThreshList",thresholdList);

%%
for i = 1:length(thresholdList)
    fName = strcat("FigureS8A_",num2str(i),"_",compose('%.2f',thresholdList(i)));
    exportgraphics(fs8a(i),strcat(fName,".tiff"));
end
exportgraphics(fs8b,'FigureS8B.pdf');

%% Figure 6
% figure 6a
f6a(1) = FM.imageConnectionMatrix_MRColorLabels(ConnectionMatrix(randomConnectomeSample_1),crossInfo);
f6a(2) = FM.imageConnectionMatrix_MRColorLabels(ConnectionMatrix(randomConnectomeSample_2),crossInfo);
%%
% figure 6b
f6b = FM.plotRandomTestCorrCoef(mouseConnectomeDataAnalysis,nullBatchRunner_1,nullBatchRunner_2);
% figure 6c
f6c = FM.swarmRandomTestCorrCoef(mouseConnectomeDataAnalysis,nullBatchRunner_1,nullBatchRunner_2);
% figure 6d
f6d = FM.swarmRandomTestAUC(mouseConnectomeDataAnalysis,nullBatchRunner_1,nullBatchRunner_2);
% figure 6e
f6e = FM.swarmRandomTestSimilarity(nullBatchRunner_1,nullBatchRunner_2);
%%
exportgraphics(f6a(1),"Figure6A_1.tiff",'BackgroundColor','none');
exportgraphics(f6a(2),"Figure6A_2.tiff",'BackgroundColor','none');
exportgraphics(f6b,'Figure6B.pdf');
exportgraphics(f6c,'Figure6C.pdf');
exportgraphics(f6d,'Figure6D.pdf');
exportgraphics(f6e,'Figure6E.pdf');

%% Figure S9
% Figure S9A
fs9a(1) = FM.swarmNullModelsCorrCoeff_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_1,nullBatchRunner_3,nullBatchRunner_5);
fs9a(2) = FM.swarmNullModelsAUC_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_1,nullBatchRunner_3,nullBatchRunner_5);
fs9a(3) = FM.swarmNullModelsPISim_3nulls(...
    nullBatchRunner_1,nullBatchRunner_3,nullBatchRunner_5);

% Figure S9B
fs9b(1) = FM.swarmNullModelsCorrCoeff_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_2,nullBatchRunner_4,nullBatchRunner_6);
fs9b(2) = FM.swarmNullModelsAUC_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_2,nullBatchRunner_4,nullBatchRunner_6);
fs9b(3) = FM.swarmNullModelsPISim_3nulls(...
    nullBatchRunner_2,nullBatchRunner_4,nullBatchRunner_6);

%%
exportgraphics(fs9a(1),'FigureS9A_1.pdf');
exportgraphics(fs9a(2),'FigureS9A_2.pdf');
exportgraphics(fs9a(3),'FigureS9A_3.pdf');
exportgraphics(fs9b(1),'FigureS9B_1.pdf');
exportgraphics(fs9b(2),'FigureS9B_2.pdf');
exportgraphics(fs9b(3),'FigureS9B_3.pdf');



%% Figure S10
fs10(1) = FM.swarmNullModelsCorrCoeff_1nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_GeneSurrogate);
fs10(2) = FM.swarmNullModelsAUC_1nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_GeneSurrogate);
fs10(3) = FM.swarmNullModelsPISim_1nulls(...
    nullBatchRunner_GeneSurrogate);

%%
exportgraphics(fs10(1),'FigureS10_1.pdf');
exportgraphics(fs10(2),'FigureS10_2.pdf');
exportgraphics(fs10(3),'FigureS10_3.pdf');


%% Figure S11
fs11(1) = FM.plotCorrVariousPC(variousPCTestResults);
fs11(2) = FM.barAUCVariousPC(variousPCTestResults);
%%
exportgraphics(fs11(1),'FigureS11_1.pdf');
exportgraphics(fs11(2),'FigureS11_2.pdf');

%% Figure S12
fs12a = FM.histogramRegionPairDistance(crossInfo,connectomeMatrix);
fs12b = FM.histogramRegionPairDistance_perMRPair(crossInfo,connectomeMatrix);

%%
exportgraphics(fs12a,'FigureS12A.pdf');
exportgraphics(fs12b,'FigureS12B.pdf');

exportgraphics(fs12a,'FigureS12A.tiff');
exportgraphics(fs12b,'FigureS12B.tiff');

%% Figure S13
% Fig S13A
fs13a(1) = FM.swarmNullModelsCorrCoeff_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_3_200um,nullBatchRunner_3,nullBatchRunner_3_50um);
fs13a(2) = FM.swarmNullModelsAUC_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_3_200um,nullBatchRunner_3,nullBatchRunner_3_50um);
fs13a(3) = FM.swarmNullModelsPISim_3nulls(...
    nullBatchRunner_3_200um,nullBatchRunner_3,nullBatchRunner_3_50um);

% Fig S13B
fs13b(1) = FM.swarmNullModelsCorrCoeff_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_4_200um,nullBatchRunner_4,nullBatchRunner_4_50um);
fs13b(2) = FM.swarmNullModelsAUC_3nulls(mouseConnectomeDataAnalysis,...
    nullBatchRunner_4_200um,nullBatchRunner_4,nullBatchRunner_4_50um);
fs13b(3) = FM.swarmNullModelsPISim_3nulls(...
    nullBatchRunner_4_200um,nullBatchRunner_4,nullBatchRunner_4_50um);

%%
exportgraphics(fs13a(1),'FigureS13A_1.pdf');
exportgraphics(fs13a(2),'FigureS13A_2.pdf');
exportgraphics(fs13a(3),'FigureS13A_3.pdf');
exportgraphics(fs13b(1),'FigureS13B_1.pdf');
exportgraphics(fs13b(2),'FigureS13B_2.pdf');
exportgraphics(fs13b(3),'FigureS13B_3.pdf');

%% Legends
f_lg_MR = FM.makeLegendMRs(brainInfo);
f_lg_con = FM.makeLegendConnectome();
f_lg_rand = FM.makeLegendRandomTest();

%%
colorListGem = orderedcolors("gem");
connNullNames = ["Original data","Random generation","Distance preserving","Network preserving"];
f_lg_connNull = FM.makeLegend_any(connNullNames,colorListGem(1:4,:),"MarkerSize",8);

nullDistNames = ["Original data","200μm bins","100μm bins","50μm bins"];
f_lg_distBin = FM.makeLegend_any(nullDistNames,colorListGem(1:4,:),"MarkerSize",8);

surrNullNames = ["Original data","Spatially randomized gene expression"];
f_lg_surrNull = FM.makeLegend_any(surrNullNames,colorListGem(1:2,:),"FDepth",50,"FWidth",250,"MarkerSize",8);

%%
colorListReef = [orderedcolors("reef")];
colorListReef(6,:) = [0.5,0.5,0.5];
majorMRPairNames = string({
    "Isocortex → Thalamus"; ...
    "Isocortex → Midbrain"; ...
    "Midbrain → Thalamus"; ...
    "Thalamus → Isocortex"; ...
    "Hypothalamus → Thalamus"; ...
    "Others"});
f_lg_majMRPair = FM.makeLegend_any(majorMRPairNames,colorListReef,"FDepth",150,"MarkerSize",8);


%%
exportgraphics(f_lg_MR,'legend_MRs.pdf');
exportgraphics(f_lg_con,'legend_connection.pdf');
exportgraphics(f_lg_rand,'legend_randomTest.pdf');

%%
exportgraphics(f_lg_connNull,'legend_nullModels.pdf');
exportgraphics(f_lg_distBin,'legend_distanceBins.pdf');
exportgraphics(f_lg_surrNull,'legend_surrogates.pdf');
%%
exportgraphics(f_lg_majMRPair,'legend_majorMRPairs.pdf');

%% Close Figures
close all
%cd(projectRoot)


%
%% Tables
[sourceSimGeneTableList,targetSimGeneTableList] =relatedGeneAnalysisResults.makeTopGeneInformationTables_withPQValues(...
    geneCosSimPvalueMat_s,geneCosSimPvalueMat_t,geneCosSimFDRQvalueMat_s,geneCosSimFDRQvalueMat_t);
% make xlsx files
%%
KTopList = 30;
for d = 1:dimPI
    simGeneTableSource = sourceSimGeneTableList{d,1};
    simGeneTableSource = simGeneTableSource(1:KTopList,:);
    simGeneTableTarget = targetSimGeneTableList{d,1};
    simGeneTableTarget = simGeneTableTarget(1:KTopList,:);
    fileNameSource = strcat("similarGeneListSource_",num2str(d),".xlsx");
    fileNameTarget = strcat("similarGeneListTarget_",num2str(d),".xlsx");
    writetable(simGeneTableSource, fileNameSource);
    writetable(simGeneTableTarget, fileNameTarget);
end
% memo: Table 1 is the concatenation of these tables

%}
%% Table S1　& S2 and prepare gene list for enrichment analysis

geneAcronymsAll = string(table2array(geneInformationTable(:,'acronym')));
geneAcronymsSmallp = string([]);

pthresh = 0.05;
ctr = 0;
for d = 1:dimPI
    simGeneTableSource = sourceSimGeneTableList{d,1};
    plist = table2array(simGeneTableSource(:,"P-value"));
    extractedTable = simGeneTableSource(plist < pthresh,:);
    nExtractedGenes = height(extractedTable);
    geneAcronymsSmallp(ctr+1:ctr+nExtractedGenes,1) = string(table2array(extractedTable(:,"Gene acronym")));
    ctr = ctr + nExtractedGenes;
    fileNameSource = strcat("similarGeneListSource_smallPvalue_",num2str(d),".xlsx");
    %writetable(extractedTable, fileNameSource);
end

for d = 1:dimPI
    simGeneTableTarget = targetSimGeneTableList{d,1};
    plist = table2array(simGeneTableTarget(:,"P-value"));
    extractedTable = simGeneTableTarget(plist < pthresh,:);
    nExtractedGenes = height(extractedTable);
    geneAcronymsSmallp(ctr+1:ctr+nExtractedGenes,1) = string(table2array(extractedTable(:,"Gene acronym")));
    ctr = ctr + nExtractedGenes;
    fileNameTarget = strcat("similarGeneListTarget_smallPvalue_",num2str(d),".xlsx");
    %writetable(extractedTable, fileNameTarget);
end

geneAcronymsSmallp = unique(geneAcronymsSmallp);

%% Make files for GO enrichment analysis

writematrix(geneAcronymsSmallp, "allSmallPvalueGenes.xlsx");
writematrix(geneAcronymsAll, "allGeneAcronyms.xlsx");