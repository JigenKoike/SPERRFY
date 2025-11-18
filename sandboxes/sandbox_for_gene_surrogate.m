% sandbox_gene_surrogate.m

load(strcat(projectRoot,"/data/processed/surrogate_genes_org.mat"));

surrogateGeneMatrices = permute(surrogate_genes,[2,3,1]);
clear surrogate_genes

%% Calculate Moran's I

%{
surrogateGeneMoransIMatrix = SpatialAutocorrelationCalculator.calculateMoransIMatrix_forSurrogates(...
    surrogateGeneMatrices(:,:,:),brainRegionDistanceMatrix,hList(1),"PrintProgress",true,"CounterInterval",10);

%%
save("surrogateGeneMoransIMatrix_org.mat","surrogateGeneMoransIMatrix");
%}
load(strcat(projectRoot,"/data/processed/surrogateGeneMoransIMatrix_org.mat"))





%% Random Test with Surrogate Genes
rng(6)
nullGenerator_dummy = nullGenerator_1;
nullBatchRunner_GeneSurrogate = NullModelAnalysisBatchRunner(nullGenerator_dummy,numRandomConnectome,holdoutParameters_nulls);
[nullBatchRunner_GeneSurrogate,surrogateWPIMatrices_s,surrogateWPIMatrices_t] = nullBatchRunner_GeneSurrogate.run_geneSurrogateNullModel(surrogateGeneMatrices(:,:,:),mouseConnectomeDataAnalysis);

%% Calculate Cosine Similarity between Surrogate Genes and Wiring PI from Surrogate
[~,dimGene,NSurrogate] = size(surrogateGeneMatrices);
surrogateGenePICosSimMatrices_s = zeros([NSurrogate,dimGene,dimPI]);
surrogateGenePICosSimMatrices_t = zeros([NSurrogate,dimGene,dimPI]);

for nsr = 1:NSurrogate
    surroGeneMat = squeeze(surrogateGeneMatrices(:,:,nsr));
    surroWPIMat_s = squeeze(surrogateWPIMatrices_s(:,:,nsr));
    surroWPIMat_t = squeeze(surrogateWPIMatrices_t(:,:,nsr));
    for dg = 1:dimGene
        surroGeneVec = surroGeneMat(:,dg);
        for dpi = 1:dimPI
            surroWPIVec_s = surroWPIMat_s(:,dpi);
            surroWPIVec_t = surroWPIMat_t(:,dpi);
            cosSim_s = corr(surroGeneVec,surroWPIVec_s);
            cosSim_t = corr(surroGeneVec,surroWPIVec_t);
            surrogateGenePICosSimMatrices_s(nsr,dg,dpi) = cosSim_s;
            surrogateGenePICosSimMatrices_t(nsr,dg,dpi) = cosSim_t;
        end
    end
end

%
% extract Top Similarity Value
surrogateTopSimilarites_s = squeeze(max(abs(surrogateGenePICosSimMatrices_s),[],2));
surrogateTopSimilarites_t = squeeze(max(abs(surrogateGenePICosSimMatrices_t),[],2));


%% Calculate P-value for each gene
geneCosSimPvalueMat_s = zeros(dimGene,dimPI);
geneCosSimPvalueMat_t = zeros(dimGene,dimPI);
for dg = 1:dimGene
    for dpi = 1:dimPI
        % source
        empiricalDist = squeeze(surrogateGenePICosSimMatrices_s(:,dg,dpi));
        cosSim_data = relatedGeneAnalysisResults.CorrCoefMatrixSource(dg,dpi);
        pval = empiricalPvalue_abs(empiricalDist,cosSim_data);
        geneCosSimPvalueMat_s(dg,dpi) = pval;
        % target
        empiricalDist = squeeze(surrogateGenePICosSimMatrices_t(:,dg,dpi));
        cosSim_data = relatedGeneAnalysisResults.CorrCoefMatrixTarget(dg,dpi);
        pval = empiricalPvalue_abs(empiricalDist,cosSim_data);
        geneCosSimPvalueMat_t(dg,dpi) = pval;
    end
end

%% FDR Correction by Benjamin-Hochberg
geneCosSimFDRQvalueMat_s = zeros(dimGene,dimPI);
geneCosSimFDRQvalueMat_t = zeros(dimGene,dimPI);
for dpi = 1:dimPI
    % source
    plist = geneCosSimPvalueMat_s(:,dpi);
    qvalList = fdr_bh(plist);
    geneCosSimFDRQvalueMat_s(:,dpi) = qvalList;
    % target
    plist = geneCosSimPvalueMat_t(:,dpi);
    qvalList = fdr_bh(plist);
    geneCosSimFDRQvalueMat_t(:,dpi) = qvalList;
end




%%
% see distribution of 1 gene
dg = 1;
dpi = 1;
figure;
binEdge = 0:0.025:1;
empiricalDist = squeeze(surrogateGenePICosSimMatrices_s(:,dg,dpi));
cosSim_data = relatedGeneAnalysisResults.CorrCoefMatrixSource(dg,dpi);
histogram(abs(empiricalDist),binEdge);
hold on 
scatter(abs(cosSim_data),0);
pval = empiricalPvalue_abs(empiricalDist,cosSim_data);
title(strcat("p-value = ",num2str(round(pval,3))))






%%

%{
%% Calculate Similarity Score of Surrogates Compared to Wiring PI Gradient
wiringPIPairs_data = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
reshapedMat = reshape(surrogateGeneMatrices,size(surrogateGeneMatrices,1),size(surrogateGeneMatrices,2)*size(surrogateGeneMatrices,3));
allSurrogatedGeneExprLevels = GeneExpressionLevels(reshapedMat,0,"dummy","dummy");
%
[surrogateCosSim_WPIs,surrogateCosSim_WPIt] = ConnectomeAnalysisRunner.calculateWiringPIGeneCorrelations(wiringPIPairs_data,allSurrogatedGeneExprLevels);

%%
save("surrogateCosSim_WPIs_org.mat","surrogateCosSim_WPIs")
save("surrogateCosSim_WPIt_org.mat","surrogateCosSim_WPIt")
%}
%load(strcat(projectRoot,"/data/prpcessed/surrogateCosSim_WPIt_org.mat"))
%load(strcat(projectRoot,"/data/prpcessed/surrogateCosSim_WPIt_org.mat"))
%%




%% See Results (Tentative)
%% See Distribution of Moran's I

figure;
binEdge = 0:0.025:1;
h1 = histogram(surrogateGeneMoransIMatrix(:),binEdge,'Normalization','probability');
hold on
h2 = histogram(spatialAutocorrelationMatrix_geneExpression(:,1),binEdge,'Normalization','probability');
legend([h2,h1],["data","surrogate"])
xlabel("Moran's I")
xlim([0,1])
ylabel("probability")
title("Moran's I of surrogate genes")

%
mean(spatialAutocorrelationMatrix_geneExpression(:,1))
mean(surrogateGeneMoransIMatrix(:))


%% Null Models
%% Compare　Correlation 

DimShow = 5;
r_data = mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutSummary.TestAllCorrelations(:,1:DimShow);
r_null1 = nullBatchRunner_1.HoldoutTestCorrMeans(:,1:DimShow);
r_null2 = nullBatchRunner_2.HoldoutTestCorrMeans(:,1:DimShow);
r_null3 = nullBatchRunner_3.HoldoutTestCorrMeans(:,1:DimShow);
r_null4 = nullBatchRunner_4.HoldoutTestCorrMeans(:,1:DimShow);
r_null5 = nullBatchRunner_5.HoldoutTestCorrMeans(:,1:DimShow);
r_null6 = nullBatchRunner_6.HoldoutTestCorrMeans(:,1:DimShow);
r_null7 = nullBatchRunner_GeneSurrogate.HoldoutTestCorrMeans(:,1:DimShow);

intervalWidth = 6;
x_data = ones(height(r_data),1) * (1:intervalWidth:1+intervalWidth*(DimShow-1));
x_sur = ones(height(r_null7),1) * (2:intervalWidth:2+intervalWidth*(DimShow-1));
x_gen = ones(height(r_null1),1) * (3:intervalWidth:3+intervalWidth*(DimShow-1));
x_dis = ones(height(r_null1),1) * (4:intervalWidth:4+intervalWidth*(DimShow-1));
x_net = ones(height(r_null1),1) * (5:intervalWidth:5+intervalWidth*(DimShow-1));

legendsStr = ["Original Connectome","Gene Surrogate","Random Generation","Distance Preserved","Network Preserved"]; 
xticksVec = [intervalWidth/2:intervalWidth:intervalWidth/2+intervalWidth*(DimShow-1)];
xticklabelsVec = [1:5];
xlimVec = [0,intervalWidth*DimShow];
ylimVec = [-0.1,1.001];


sz = 3;
szData = 6;
colorListTmp = colorList13([1,9,3,7,4],:);

%%
figure;
h1 = swarmchart(x_data,r_data,szData,colorListTmp(1,:),'filled');
hold on
h2 = swarmchart(x_sur,r_null7,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_gen,r_null2,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_dis,r_null4,sz,colorListTmp(4,:),'filled');
hold on
h5 = swarmchart(x_net,r_null6,sz,colorListTmp(5,:),'filled');

legend([h1(1),h2(1),h3(1),h4(1),h5(1)],legendsStr);
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
xticks(xticksVec)
xticklabels(xticklabelsVec)
ylim(ylimVec)
ylabel("Correlation Coefficient")
title("Randomized Model (Locally)")

%% AUC
%% AUC
AUC_data = mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC.';
AUC_null1 = nullBatchRunner_1.HoldoutTestAUCMeans;
AUC_null2 = nullBatchRunner_2.HoldoutTestAUCMeans;
AUC_null3 = nullBatchRunner_3.HoldoutTestAUCMeans;
AUC_null4 = nullBatchRunner_4.HoldoutTestAUCMeans;
AUC_null5 = nullBatchRunner_5.HoldoutTestAUCMeans;
AUC_null6 = nullBatchRunner_6.HoldoutTestAUCMeans;
AUC_null7 = nullBatchRunner_GeneSurrogate.HoldoutTestAUCMeans;

%%
x_data = ones(height(AUC_data),1) * 1;
x_sur = ones(height(AUC_null7),1) * 2;
x_gen = ones(height(AUC_null1),1) * 3;
x_dis = ones(height(AUC_null1),1) * 4;
x_net = ones(height(AUC_null1),1) * 5;

xlimVec = [0,6];
ylimVec = [0.4,1];

figure;
h1 = swarmchart(x_data,AUC_data,szData,colorListTmp(1,:),'filled');
hold on
h2 = swarmchart(x_sur,AUC_null7,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_gen,AUC_null2,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_dis,AUC_null4,sz,colorListTmp(4,:),'filled');
hold on
h5 = swarmchart(x_net,AUC_null6,sz,colorListTmp(5,:),'filled');
legend([h1,h2,h3,h4,h5],legendsStr);
xlim(xlimVec)
xlabel("Model")
xticks([])
ylim(ylimVec)
yticks([0.5:0.1:1])
ylabel("AUC")
title("AUC of Randomized Model")

%% PI Similarity

DimShow = 5;

PISim_null1 = nullBatchRunner_1.WirigPISimilarities(:,1:DimShow);
PISim_null2 = nullBatchRunner_2.WirigPISimilarities(:,1:DimShow);
PISim_null3 = nullBatchRunner_3.WirigPISimilarities(:,1:DimShow);
PISim_null4 = nullBatchRunner_4.WirigPISimilarities(:,1:DimShow);
PISim_null5 = nullBatchRunner_5.WirigPISimilarities(:,1:DimShow);
PISim_null6 = nullBatchRunner_6.WirigPISimilarities(:,1:DimShow);
PISim_null7 = nullBatchRunner_GeneSurrogate.WirigPISimilarities(:,1:DimShow);


%%
intervalWidth = 5;
x_sur = ones(height(PISim_null7),1) * (1:intervalWidth:1+intervalWidth*(DimShow-1));
x_gen = ones(height(PISim_null1),1) * (2:intervalWidth:2+intervalWidth*(DimShow-1));
x_dis = ones(height(PISim_null1),1) * (3:intervalWidth:3+intervalWidth*(DimShow-1));
x_net = ones(height(PISim_null1),1) * (4:intervalWidth:4+intervalWidth*(DimShow-1));
xticklabelsVec = [1:5];
xticksVec = [intervalWidth/2:intervalWidth:2.5+intervalWidth*(DimShow-1)];
xlimVec = [0,intervalWidth*DimShow];
ylimVec = [0,1];

figure;
h2 = swarmchart(x_sur,PISim_null7,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_gen,PISim_null2,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_dis,PISim_null4,sz,colorListTmp(4,:),'filled');
hold on
h5 = swarmchart(x_net,PISim_null6,sz,colorListTmp(5,:),'filled');
legend([h2(1),h3(1),h4(1),h5(1)],legendsStr(2:5));
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
xticks(xticksVec)
xticklabels(xticklabelsVec)
ylim(ylimVec)
ylabel("PI Similarity Score")
title("Randomized Model (Locally)")


%% See Distribution of Cosine Similarity Compared to Wiring PI Gradients
figure;
t = tiledlayout(2,5);
binEdge = -1:0.025:1;
for d = 1:5
nexttile
cosSimDist = squeeze(surrogateGenePICosSimMatrices_s(:,:,d));
h1 = histogram(cosSimDist(:),binEdge,'Normalization','probability');
hold on
h2 = histogram(relatedGeneAnalysisResults.CorrCoefMatrixSource(:,d),binEdge,'Normalization','probability');
title(strcat("PI_s^{(",num2str(d),")}"));
xlim([-1,1])
end
legend([h2,h1],["Data","Surrogate Test"])
for d = 1:5
nexttile
cosSimDist = squeeze(surrogateGenePICosSimMatrices_t(:,:,d));
histogram(cosSimDist(:),binEdge,'Normalization','probability');
hold on
histogram(relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,d),binEdge,'Normalization','probability');
title(strcat("PI_t^{(",num2str(d),")}"));
xlim([-1,1])
end
xlabel(t,"Cosine Similarity")
title(t,"Distribution of Cosine Similarity (Surrogate & Data)")

%%
figure;
t = tiledlayout(2,5);
binEdge = 0:0.025:1;
for d = 1:5
nexttile
cosSimDist = squeeze(surrogateGenePICosSimMatrices_s(:,:,d));
h1 = histogram(abs(cosSimDist(:)),binEdge,'Normalization','probability');
hold on
h2 = histogram(abs(relatedGeneAnalysisResults.CorrCoefMatrixSource(:,d)),binEdge,'Normalization','probability');
title(strcat("PI_s^{(",num2str(d),")}"));
xlim([0,1])
end
legend([h2,h1],["Data","Surrogate Test"])
for d = 1:5
nexttile
cosSimDist = squeeze(surrogateGenePICosSimMatrices_t(:,:,d));
histogram(abs(cosSimDist(:)),binEdge,'Normalization','probability');
hold on
histogram(abs(relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,d)),binEdge,'Normalization','probability');
title(strcat("PI_t^{(",num2str(d),")}"));
xlim([0,1])
end
xlabel(t,"Cosine Similarity (Abs)")
title(t,"Distribution of Cosine Similarity (Abs) (Surrogate & Data)")


%% Distribution of Top Similarity Score 
KTop = 10;
figure;
t = tiledlayout(2,5);
binEdge = 0:0.025:1;
for d = 1:5
nexttile
h1 = histogram(surrogateTopSimilarites_s(:,d),binEdge,'Normalization','probability');
hold on
ktopSimData = maxk(relatedGeneAnalysisResults.CorrCoefMatrixSource(:,d),KTop,"ComparisonMethod","abs");
h2 = scatter(abs(ktopSimData),zeros(size(ktopSimData)),'filled');
title(strcat("PI_s^{(",num2str(d),")}"));
xlim([0,1])
end
legend([h2,h1],["Top 10 Genes of Data","Surrogate Test"])
for d = 1:5
nexttile
histogram(surrogateTopSimilarites_t(:,d),binEdge,'Normalization','probability');
hold on
ktopSimData = maxk(relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,d),KTop,"ComparisonMethod","abs");
h2 = scatter(abs(ktopSimData),zeros(size(ktopSimData)),'filled');
title(strcat("PI_t^{(",num2str(d),")}"));
xlim([0,1])
end
xlabel(t,"Cosine Similarity (Abs)")
title(t,"Distribution of Top Similarity Score (Surrogate)")

%% plot cosSim vs pVal

figure;
t = tiledlayout(2,dimPI);
for dpi = 1:dimPI
    nexttile
    cosSimList = relatedGeneAnalysisResults.CorrCoefMatrixSource(:,dpi);
    pValList = geneCosSimPvalueMat_s(:,dpi);
    scatter(abs(cosSimList),pValList,'.');
    r = corr(pValList,abs(cosSimList));
    title(strcat("PI_s^{(",num2str(dpi),")}"))
    ylabel("p-value")
    xlabel("Cosine Similarity (abs)")
    set(gca, 'YScale', 'log') 
    xlim([0,1])
    ylim([10^(-4),1])
end
for dpi = 1:dimPI
    nexttile
    cosSimList = relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,dpi);
    pValList = geneCosSimPvalueMat_t(:,dpi);
    scatter(abs(cosSimList),pValList,'.');
    r = corr(pValList,abs(cosSimList));
    title(strcat("PI_t^{(",num2str(dpi),")}"))
    ylabel("p-value")
    xlabel("Cosine Similarity (abs)")
    set(gca, 'YScale', 'log') 
    xlim([0,1])
    ylim([10^(-4),1])
end


%% plot cosSim vs q-value

figure;
t = tiledlayout(2,dimPI);
qthresh = 0.05;
for dpi = 1:dimPI
    nexttile
    cosSimList = relatedGeneAnalysisResults.CorrCoefMatrixSource(:,dpi);
    pValList = geneCosSimFDRQvalueMat_s(:,dpi);
    scatter(abs(cosSimList),pValList,'.');
    
    title(strcat("PI_s^{(",num2str(dpi),")}"))
    ylabel("FDR q-value")
    xlabel("Cosine Similarity (abs)")
    set(gca, 'YScale', 'log') 
    hold on
    plot([0,1],[qthresh,qthresh])
    xlim([0,1])
    %ylim([0,1])
    ylim([10^(-2),1])
end
for dpi = 1:dimPI
    nexttile
    cosSimList = relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,dpi);
    pValList = geneCosSimFDRQvalueMat_t(:,dpi);
    scatter(abs(cosSimList),pValList,'.');
    
    title(strcat("PI_t^{(",num2str(dpi),")}"))
    ylabel("FDR q-value")
    xlabel("Cosine Similarity (abs)")
    set(gca, 'YScale', 'log') 
    hold on
    plot([0,1],[qthresh,qthresh])
    xlim([0,1])
    %ylim([0,1])
    ylim([10^(-2),1])
end
