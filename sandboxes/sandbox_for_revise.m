%% sandbox for revise

%% data preparation
%% Startup
clear; clc;
startup_SPERRFY;  % path setting

%% Load data
main_import_processed_data;  % import data (.mat)
main_data_preparation;
main_parameter_setting;


%% Run main data analysis
rng(1)
mouseConnectomeDataAnalysis = ConnectomeAnalysisRunner(connectomeMatrix, geneExprLevels,crossInfo, ...
    pairDataGenerationOptions,reconstructionParameters,holdoutParameters_main,"FactoryDescription","Main Data Analysis");
mouseConnectomeDataAnalysis = mouseConnectomeDataAnalysis.runAnalysis(ReconstructionStoreTag=true);

%% Perform gene analysis
KTopGene = 30;
relatedGeneAnalysisResults = mouseConnectomeDataAnalysis.performRelatedGeneAnalysis(KTopGene);


%% Make Reciprocal & Unidirectional Index
reciprocalFeatureMasks = ConnectomeFeatureAnalyzer.makeReciprocalUnidirectionalMasks(connectomeMatrix,crossInfo);



%% Reconstruction by distance matrix
normalizedDistanceMatrix = crossInfo.DistanceMatrix / max(crossInfo.DistanceMatrix,[],'all');
distanceExponentialReconstResults = ReconstructionResults(reconstructionParameters);
distanceExponentialReconstResults = distanceExponentialReconstResults.compute(mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData, ...
    connectomeMatrix.Matrix,crossInfo.DomainDefineMatrix,"ArbitralPIDistanceMatrix",normalizedDistanceMatrix,"StoreMatrixTag",true);

%%
figure;
hROC = distanceExponentialReconstResults.plotROC();
figure;
sampleMat = cell2mat(distanceExponentialReconstResults.PredictedMatrices(1,20));
imagesc(sampleMat);


%% See Similarities between PIs

wiringPIdata = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
PIMatSource = wiringPIdata.WiringPISource(:,1:dimPI);
PIMatTarget = wiringPIdata.WiringPITarget(:,1:dimPI);

PISim_ss = corr(PIMatSource,PIMatSource);
PISim_tt = corr(PIMatTarget,PIMatTarget);
PISim_st = corr(PIMatSource,PIMatTarget);

%%
figure;
t = tiledlayout(1,3);
title(t,"PI Similarity")
nexttile
heatmap(PISim_ss)
colormap(redblue_cp)
clim([-1,1])
title("Source - Source")
nexttile
heatmap(PISim_tt)
colormap(redblue_cp)
clim([-1,1])
title("Target - Target")
nexttile
heatmap(PISim_st)
colormap(redblue_cp)
clim([-1,1])
title("Source - Target")

%% Calculate Spatial Autocorrelation
regionDistanceMedian = SpatialAutocorrelationCalculator.gethFromMedian(brainRegionDistanceMatrix);
hList_exponential = SpatialAutocorrelationCalculator.makehList_exponentialScale(regionDistanceMedian,2,-7,1);
hList_quantile = SpatialAutocorrelationCalculator.makehList_quantile(brainRegionDistanceMatrix,0.01);

%%
figure;histogram(brainRegionDistanceMatrix(brainRegionDistanceMatrix>0))
hold on
scatter(hList_exponential,zeros(size(hList_exponential)),'.','r');
hold on
scatter(hList_quantile,zeros(size(hList_quantile)),'.','b');

%% 
% PI
hList = hList_exponential;
spatialAutocorrelationMatrix_PISource = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    wiringPIdata.WiringPISource(:,:),brainRegionDistanceMatrix,hList);
spatialAutocorrelationMatrix_PITarget = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    wiringPIdata.WiringPITarget(:,:),brainRegionDistanceMatrix,hList);

% genes
spatialAutocorrelationMatrix_geneExpression = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    geneExprLevels.ExpressionMatrix,brainRegionDistanceMatrix,hList);

%% PC
spatialAutocorrelationMatrix_PC = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    geneExprLevels.PCAMatrix,brainRegionDistanceMatrix,hList);

%% 　See Distribution of Spatial Autocorrelation
% show only lowest h
% plot about Wiring PI
figure;
plot([1:50],spatialAutocorrelationMatrix_PISource(:,1),'r','LineStyle','-','Marker','.','MarkerSize',12);
hold on
plot([1:50],spatialAutocorrelationMatrix_PITarget(:,1),'b','LineStyle','-','Marker','.','MarkerSize',12);
ylim([0,1])
legend(["Source","Target"])
title("Spatial Autocorrelation of Wiring PI")

%%
figure;
X = [spatialAutocorrelationMatrix_PISource(:,1) spatialAutocorrelationMatrix_PITarget(:,1)];
bar(X)

%%
% plot of PC
figure;
%plot([1:height(spatialAutocorrelationMatrix_PC)],spatialAutocorrelationMatrix_PC(:,1),'-','Marker','.');
plot([1:50],spatialAutocorrelationMatrix_PC(1:50,1),'k','LineStyle','-','Marker','.');
ylim([0,1])

%% plot PI & PC

figure;
plot([1:50],spatialAutocorrelationMatrix_PISource(:,1),'r','LineStyle','-','Marker','.');
hold on
plot([1:50],spatialAutocorrelationMatrix_PITarget(:,1),'b','LineStyle','-','Marker','.');
hold on
plot([1:50],spatialAutocorrelationMatrix_PC(1:50,1),'k','LineStyle','-','Marker','.');
ylim([0,1])
legend(["Source PI","Target PI","PC"])
title("Spatial Autocorrelation of Wiring PI & PC")

%% histogram of gene
figure;
histogram(spatialAutocorrelationMatrix_geneExpression(:,1),50);
xlim([0,1])
title("Spatial Autocorrelation of Gene Expression")

%% PI Similarity vs Autocorrelation: about Genes
figure;
t = tiledlayout(2,dimPI);
title(t,"PI Similarity Abs (y axis) vs Spatial Autocorrelation (x axis)")
x = spatialAutocorrelationMatrix_geneExpression(:,1);
for d = 1:dimPI
    nexttile
    y = abs(relatedGeneAnalysisResults.CorrCoefMatrixSource(:,d));
    scatter(x,y,'.');
    r = corr(x,y);
    hold on
    spAuCor_PI = spatialAutocorrelationMatrix_PISource(d,1);
    plot([spAuCor_PI,spAuCor_PI],[0,1]);
    titleStr = strcat("{PI_s}^{",num2str(d),"} : r = ",num2str(round(r,3)));
    title(titleStr);
    xlim([0,1])
    ylim([0,1])
end
for d = 1:dimPI
    nexttile
    y = abs(relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,d));
    scatter(x,y,'.');
    r = corr(x,y);
    hold on
    spAuCor_PI = spatialAutocorrelationMatrix_PITarget(d,1);
    plot([spAuCor_PI,spAuCor_PI],[0,1]);
    titleStr = strcat("{PI_t}^{",num2str(d),"} : r = ",num2str(round(r,3)));
    title(titleStr);
    xlim([0,1])
    ylim([0,1])
end

%% Autocorrelation diff
figure;
t = tiledlayout(2,dimPI);
title(t,"PI Similarity Abs (y axis) vs Spatial Autocorrelation (x axis)")
x = spatialAutocorrelationMatrix_geneExpression(:,1);
for d = 1:dimPI
    nexttile
    x = abs(spatialAutocorrelationMatrix_geneExpression(:,1) - spatialAutocorrelationMatrix_PISource(d,1));
    y = abs(relatedGeneAnalysisResults.CorrCoefMatrixSource(:,d));
    scatter(x,y,'.');
    r = corr(x,y);
    titleStr = strcat("{PI_s}^{",num2str(d),"} : r = ",num2str(round(r,3)));
    title(titleStr);
    xlim([0,1])
    ylim([0,1])
end
for d = 1:dimPI
    nexttile
    x = abs(spatialAutocorrelationMatrix_geneExpression(:,1) - spatialAutocorrelationMatrix_PITarget(d,1));
    y = abs(relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,d));
    scatter(x,y,'.');
    r = corr(x,y);
    titleStr = strcat("{PI_t}^{",num2str(d),"} : r = ",num2str(round(r,3)));
    title(titleStr);
    xlim([0,1])
    ylim([0,1])
end

%% Reconstruction by Distance Matrix (isotoropic, exponential)
distanceMatrix = brainRegionDistanceMatrix;
LList = length(hList);
distAUCList = zeros(LList,1);
for l = 1:LList
    exponentialWeightMatrix = SpatialAutocorrelationCalculator.distance2exponentialWeight( ...
        distanceMatrix,hList(l),0);
    arbitralMatrix = ones(size(exponentialWeightMatrix)) - exponentialWeightMatrix;
    distanceExponentialReconstResults = ReconstructionResults(reconstructionParameters);
    distanceExponentialReconstResults = distanceExponentialReconstResults.compute(mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData, ...
        connectomeMatrix.Matrix,crossInfo.DomainDefineMatrix,"ArbitralPIDistanceMatrix",arbitralMatrix,"StoreMatrixTag",true);
    distAUCList(l) = distanceExponentialReconstResults.AUC_ROC;
end
%
figure;
plot(hList,distAUCList,"LineStyle",'-',"Marker",".");
set(gca, 'XScale', 'log')
box off
ylim([0.5,1])

%% Histogram of connection distance
distanceList_connections = distanceMatrix(connectomeMatrix.Matrix .* crossInfo.DomainDefineMatrix == 1);
distanceList_regionPairs = distanceMatrix(crossInfo.DomainDefineMatrix == 1);
binEdges = 0:100:13000;
figure;
histogram(distanceList_regionPairs,'BinEdges',binEdges);
hold on
histogram(distanceList_connections,'BinEdges',binEdges);
legend(["All region pairs","Connections"])

%% Per MR Pair
MRNameList = unique(brainInfo.getMajorRegions,'stable');

MRPairMask = crossInfo.makeSpecificMRPairMask(MRNameList(1),MRNameList(7));
distanceList_connections = distanceMatrix(connectomeMatrix.Matrix .* MRPairMask == 1);
distanceList_regionPairs = distanceMatrix(MRPairMask == 1);
figure;
histogram(distanceList_regionPairs,'BinEdges',binEdges);
hold on
histogram(distanceList_connections,'BinEdges',binEdges);
legend(["All region pairs","Connections"])
titleStr = strcat(MRNameList(1)," -> ",MRNameList(7));
title(titleStr)
%%
figure;
t = tiledlayout(13,13);
ctr = 0;
for sMRInd = 1:13
    sMR = MRNameList(sMRInd);
    for tMRInd = 1:13
        tMR = MRNameList(tMRInd);
        ctr = ctr + 1;
        MRPairMask = crossInfo.makeSpecificMRPairMask(sMR,tMR);
        if sMR == tMR

        else
            nexttile(ctr)
            distanceList_connections = distanceMatrix(connectomeMatrix.Matrix .* MRPairMask == 1);
            distanceList_regionPairs = distanceMatrix(MRPairMask == 1);
            histogram(distanceList_regionPairs,'BinEdges',binEdges, ...
                'EdgeColor', 'none', ...              
                'FaceAlpha', 0.5);              
            hold on
            if ~isempty(distanceList_connections)
                histogram(distanceList_connections,'BinEdges',binEdges, ...
                'EdgeColor', 'none', ...              
                'FaceAlpha', 1);
                titleStr = strcat("#Connections: ",num2str(sum(connectomeMatrix.Matrix .* MRPairMask == 1,"all")));
                title(titleStr)
            else
                title("No Connections")
            end
            %legend(["All region pairs","Connections"])
            ylim([0,100])
            box off
            %titleStr = strcat(sMR," -> ",tMR);
                      
            if tMRInd == 1
                ylabel(sMR)
            end
            if sMRInd == 13
                xlabel(tMR)
            end
        end
    end
end

%% Connection Distance vs PI Diff
PIDiffMat = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.PIDistanceMatrix;

sz = 3;
figure;
t = tiledlayout(1,3);
x_all = distanceMatrix(crossInfo.DomainDefineMatrix == 1);
y_all = PIDiffMat(crossInfo.DomainDefineMatrix == 1);
r_all = corr(x_all,y_all);
x_conn = distanceMatrix(connectomeMatrix.Matrix .* crossInfo.DomainDefineMatrix == 1);
y_conn = PIDiffMat(connectomeMatrix.Matrix .* crossInfo.DomainDefineMatrix == 1);
r_conn = corr(x_conn,y_conn);
x_noConn = distanceMatrix((~connectomeMatrix.Matrix) .* crossInfo.DomainDefineMatrix == 1);
y_noConn = PIDiffMat((~connectomeMatrix.Matrix) .* crossInfo.DomainDefineMatrix == 1);
r_noConn = corr(x_noConn,y_noConn);
nexttile
scatter(x_noConn,y_noConn,sz,'k','filled','MarkerFaceAlpha',1);
hold on
scatter(x_conn,y_conn,sz,'red','filled');
xlabel("Region Distance")
ylabel("PI Distance")
titleStr = strcat("All Region Pair: r = ",num2str(round(r_all,4)));
title(titleStr)
xlim([0,12500])
ylim([0,1])
nexttile
scatter(x_conn,y_conn,sz,'red','filled');
xlabel("Region Distance")
ylabel("PI Distance")
titleStr = strcat("Conection: r = ",num2str(round(r_conn,4)));
title(titleStr)
xlim([0,12500])
ylim([0,1])
nexttile
scatter(x_noConn,y_noConn,sz,'k','filled');
xlabel("Region Distance")
titleStr = strcat("No Conection: r = ",num2str(round(r_noConn,4)));
title(titleStr)
xlim([0,12500])
ylim([0,1])


%% See Reciprocal & Unidirectional

x_rec = distanceMatrix(reciprocalFeatureMasks.reciprocal == 1);
y_rec = PIDiffMat(reciprocalFeatureMasks.reciprocal == 1);
x_uni = distanceMatrix(reciprocalFeatureMasks.unidirectional == 1);
y_uni = PIDiffMat(reciprocalFeatureMasks.unidirectional == 1);
x_cou = distanceMatrix(reciprocalFeatureMasks.unidirectionalCounter == 1);
y_cou = PIDiffMat(reciprocalFeatureMasks.unidirectionalCounter == 1);
x_oth = distanceMatrix(reciprocalFeatureMasks.others == 1);
y_oth = PIDiffMat(reciprocalFeatureMasks.others == 1);

figure;
sz = 3;
h1 = scatter(x_oth,y_oth,sz,'k','filled','MarkerFaceAlpha',1);
hold on
h2 = scatter(x_cou,y_cou,sz,'green','filled');
hold on
h3 = scatter(x_uni,y_uni,sz,'blue','filled','MarkerFaceAlpha',1);
hold on
h4 = scatter(x_rec,y_rec,sz,'red','filled');
xlabel("Region Distance")
ylabel("PI Distance")
titleStr = strcat("All Region Pair: r = ",num2str(round(r_all,4)));
title(titleStr)
legend([h4 h3 h2 h1],["Reciprocal","Unidirectional","Uniderectional Counter","Others"]);
xlim([0,12500])
ylim([0,1])

%% Histogram
figure;
binEdges = 0:100:13000;
alpha = 0.3;
h1 = histogram(x_oth,'BinEdges',binEdges,'Normalization','probability',"FaceColor","black","FaceAlpha",alpha);
hold on
h2 = histogram(x_cou,'BinEdges',binEdges,'Normalization','probability',"FaceColor","green","FaceAlpha", alpha);
hold on
h3 = histogram(x_uni,'BinEdges',binEdges,'Normalization','probability',"FaceColor","blue","FaceAlpha",alpha);
hold on
h4 = histogram(x_rec,'BinEdges',binEdges,'Normalization','probability',"FaceColor","red","FaceAlpha",alpha);

titleStr = strcat("Histogram of Region Distance for Each Region Pair");
title(titleStr)
legend([h4 h3 h2 h1],["Reciprocal","Unidirectional","Uniderectional Counter","Others"]);
xlabel("Region Distance")
ylabel("Proportion")

%%
figure;
binEdges = 0:0.01:1;
h1 = histogram(y_oth,'BinEdges',binEdges,'Normalization','probability',"FaceColor","black","FaceAlpha",alpha);
hold on
h2 = histogram(y_cou,'BinEdges',binEdges,'Normalization','probability',"FaceColor","green","FaceAlpha",alpha);
hold on
h3 = histogram(y_uni,'BinEdges',binEdges,'Normalization','probability',"FaceColor","blue","FaceAlpha",alpha);
hold on
h4 = histogram(y_rec,'BinEdges',binEdges,'Normalization','probability',"FaceColor","red","FaceAlpha",alpha);
titleStr = strcat("Histogram of Wiring PI Distance for Each Region Pair");
title(titleStr)
legend([h4 h3 h2 h1],["Reciprocal","Unidirectional","Uniderectional Counter","Others"]);
xlabel("PI Difference (Normalized)")
ylabel("Proportion")

%% per PI Component
figure;
tiledlayout(2,5);

binEdges = 0:0.01:1;
alpha = 0.3;
for d = 1:10
    nexttile
    diffMat_PIComponent = abs(wiringPIdata.WiringPISource(:,d) - wiringPIdata.WiringPITarget(:,d)');
    diffMat_PIComponent = diffMat_PIComponent / max(diffMat_PIComponent,[],'all');
    diff_rec = diffMat_PIComponent(reciprocalFeatureMasks.reciprocal == 1);
    diff_uni = diffMat_PIComponent(reciprocalFeatureMasks.unidirectional == 1);
    diff_cou = diffMat_PIComponent(reciprocalFeatureMasks.unidirectionalCounter == 1);
    diff_oth = diffMat_PIComponent(reciprocalFeatureMasks.others == 1);
    h1 = histogram(diff_oth,'BinEdges',binEdges,'Normalization','probability',"FaceColor","black","FaceAlpha",alpha);
    hold on
    h2 = histogram(diff_cou,'BinEdges',binEdges,'Normalization','probability',"FaceColor","green","FaceAlpha",alpha);
    hold on
    h3 = histogram(diff_uni,'BinEdges',binEdges,'Normalization','probability',"FaceColor","blue","FaceAlpha",alpha);
    hold on
    h4 = histogram(diff_rec,'BinEdges',binEdges,'Normalization','probability',"FaceColor","red","FaceAlpha",alpha);
    titleStr = strcat("PI^{(",num2str(d),")}");
    title(titleStr)
    legend([h4 h3 h2 h1],["Reciprocal","Unidirectional","Uniderectional Counter","Others"]);
    xlabel("PI Difference (Normalized)")
    ylabel("Proportion")
end


%% Gene List of High Autocorrelation
[~,topKGeneIndexLists_autocorr] = maxk(spatialAutocorrelationMatrix_geneExpression(:,1),KTopGene,1,"ComparisonMethod","abs");
geneAcronyms_all = geneInformationTable(:,'acronym');
geneNames_all = geneInformationTable(:,'name');
topAutocorrGenes_acronyms = geneAcronyms_all(topKGeneIndexLists_autocorr,:);
topAutocorrGenes_names = geneNames_all(topKGeneIndexLists_autocorr,:);

%% Make Gene Name List (All)
geneAcronyms_all = geneInformationTable(:,'acronym');

%writetable(geneAcronyms_all,'GeneList_Acronym.csv','WriteVariableNames',false);

%% check gene distribution by GO
% run "sandbox_gProfiler.m before runnning this part

goColumnID = 2; % Axon Guidance, etc
goIndexList = L_any(:,goColumnID);
figure;
t = tiledlayout(2,dimPI);
title(t,Legend_sel.TermName(goColumnID))
x_go = spatialAutocorrelationMatrix_geneExpression(goIndexList,1);
x_except = spatialAutocorrelationMatrix_geneExpression(~goIndexList,1);

for d = 1:dimPI
    nexttile
    y_go = abs(relatedGeneAnalysisResults.CorrCoefMatrixSource(goIndexList,d));
    y_except = abs(relatedGeneAnalysisResults.CorrCoefMatrixSource(~goIndexList,d));
    scatter(x_except,y_except,'.','k');
    hold on
    scatter(x_go,y_go,'.','r');
    r_go = corr(x_go,y_go);
    r_except = corr(x_except,y_except);
    titleStr = strcat("{PI_s}^{",num2str(d),"} : r = ",num2str(round(r_go,3)));
    title(titleStr);
    xlim([0,1])
    ylim([0,1])
end

for d = 1:dimPI
    nexttile
    y_go = abs(relatedGeneAnalysisResults.CorrCoefMatrixTarget(goIndexList,d));
    y_except = abs(relatedGeneAnalysisResults.CorrCoefMatrixTarget(~goIndexList,d));
    scatter(x_except,y_except,'.','k');
    hold on
    scatter(x_go,y_go,'.','r');
    r_go = corr(x_go,y_go);
    r_except = corr(x_except,y_except);
    titleStr = strcat("{PI_s}^{",num2str(d),"} : r = ",num2str(round(r_go,3)));
    title(titleStr);
    xlim([0,1])
    ylim([0,1])
end

%% Make GO Index List for Each PI
topKGeneGOLogicalMat_source = zeros(KTopGene,dimPI);
topKGeneGOLogicalMat_target = zeros(KTopGene,dimPI);

for d = 1:dimPI
    topKGeneGOLogicalMat_source(:,d) = goIndexList(relatedGeneAnalysisResults.TopKGeneIndexListsSource(:,d));
    topKGeneGOLogicalMat_target(:,d) = goIndexList(relatedGeneAnalysisResults.TopKGeneIndexListsTarget(:,d));
end


figure;
t = tiledlayout(1,2);
title(t,Legend_sel.TermName(goColumnID));
nexttile
imagesc(topKGeneGOLogicalMat_source);
title("Source")
yticks(1:30)
nexttile
imagesc(topKGeneGOLogicalMat_target);
title("target")
yticks(1:30)

%% Make MR Index of Connection Pair

focusingMRPairList = string({
    "Isocortex","Thalamus"; ...
    "Isocortex","Midbrain"; ...
    "Midbrain","Thalamus"; ...
    "Thalamus","Isocortex"; ...
    "Hypothalamus","Thalamus"});


speceficMRPairIndexList = zeros(height(mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.CCAResultsData.U),1);

for k = 1:height(focusingMRPairList)
    [spcfcMRPairMask,~,~] = mouseConnectomeDataAnalysis.Factory.getSpecificMRPairMask( ...
        focusingMRPairList(k,1),focusingMRPairList(k,2),pairDataGenerationOptions);
    speceficMRPairIndexList = speceficMRPairIndexList + spcfcMRPairMask * k;
end

%%

colorListTmp = [0.2,0.2,0.2;orderedcolors("reef")];
legendStr = string({
    "Others"; ...
    "Isocortex -> Thalamus"; ...
    "Isocortex -> Midbrain"; ...
    "Midbrain -> Thalamus"; ...
    "Thalamus -> Isocortex"; ...
    "Hypothalamus -> Thalamus"});

%%
sz = 12;
U = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.CCAResultsData.U;
V = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.CCAResultsData.V;
figure;
t = tiledlayout(1,5);
for d = 1:5
    nexttile
    for k = 0:5
        u = U(speceficMRPairIndexList==k,d);
        v = V(speceficMRPairIndexList==k,d);
        scatter(u,v,sz,colorListTmp(k+1,:),'filled');
        hold on
    end
    box off
    axis square
end
legend(legendStr)



%% Gene Expression Shuffling per Major Region

% genes
aaa_geneSA_data = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    geneExprLevels.ExpressionMatrix,brainRegionDistanceMatrix,hList(1));

aaa_shuffledGeneMatrix = geneExprLevels.makeRegionShuffledGeneExprLevels_withinMR(brainInfo);

aaa_geneSA_shuffled = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    aaa_shuffledGeneMatrix,brainRegionDistanceMatrix,hList(1));

%%
rng(5)
NShuffle = 10;
aaa_geneSA_shuffled_List = zeros(height(aaa_geneSA_data),NShuffle);
for n = 1:NShuffle
    aaa_shuffledGeneMatrix = geneExprLevels.makeRegionShuffledGeneExprLevels_withinMR(brainInfo);
    aaa_geneSA_shuffled_List(:,n) = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
        aaa_shuffledGeneMatrix,brainRegionDistanceMatrix,hList(1));
end

%%
figure;
binEdges = 0:0.05:1;
histogram(aaa_geneSA_data,binEdges);
hold on
histogram(aaa_geneSA_shuffled,binEdges)
legend(["data","shuffled"])
title("spatial autocorrelation of gene expression")

%%
mean(aaa_geneSA_shuffled,"all")
mean(aaa_geneSA_data)


%% See Sample of BrainSmash Results

load(strcat(projectRoot,"/sandboxes/brainsmashData/surrogates.mat"));
%%
aaa_gene1 = geneExpressionMatrix(:,1);
aaa_gene1_surrogatedMatrix = surrogates.';

%%
aaa_gene1_SA = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    aaa_gene1,brainRegionDistanceMatrix,hList(1));

%%
aaa_gene1_surrogated_SA = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    aaa_gene1_surrogatedMatrix,brainRegionDistanceMatrix,hList(1));

%%
figure;
x = 1;
h1 = scatter(x,aaa_gene1_surrogated_SA,20,'k','filled');
hold on
h2 = scatter(x,aaa_gene1_SA,36,'r','filled');
ylim([0,1])
ylabel("Moran's I")
legend([h2,h1],["gene1","surrogates"])

%% See Full Data of Surrogates
%load(strcat(projectRoot,"/sandboxes/brainsmashData/surrogate_genes_org.mat"));
load(strcat(projectRoot,"/sandboxes/brainsmashData/surrogate_genes.mat"));

%%
figure;
for d = 1:10
aaa_gene_d = geneExpressionMatrix(:,d);
aaa_gene_d_surrogatedMatrix = permute(squeeze(surrogate_genes(:,:,d)),[2,1]);
%
aaa_gene_d_SA = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    aaa_gene_d,brainRegionDistanceMatrix,hList(1));

aaa_gene_d_surrogated_SA = SpatialAutocorrelationCalculator.calculateMoransIMatrix( ...
    aaa_gene_d_surrogatedMatrix,brainRegionDistanceMatrix,hList(1));

%

x = d*ones(size(aaa_gene_d_surrogated_SA));
h1 = swarmchart(x,aaa_gene_d_surrogated_SA,12,'k','filled');
hold on
h2 = scatter(d,aaa_gene_d_SA,36,'r','filled');
ylim([0,1])
end
%%
ylabel("Moran's I")
title("Samples of Spatial Autocorrelation of Surrogates (surrogate\_genes)")
xlabel("Gene Index")
%%
xlim([0,11])
legend([h2(1),h1(1)],["gene_data","surrogates"])

%%
aaa_surrogateGenes = permute(surrogate_genes,[2,3,1]);
aaa_MoransIMat_sample = SpatialAutocorrelationCalculator.calculateMoransIMatrix_forSurrogates(...
    aaa_surrogateGenes(:,1:10,1:50),brainRegionDistanceMatrix,hList(1));

%%