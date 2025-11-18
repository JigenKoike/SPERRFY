% main_gene_surrogate.m

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
%load(strcat(projectRoot,"/data/processed/surrogateGeneMoransIMatrix_org.mat"))


%% Calculate Similarity between Surrogate Genes and Original Wiring PI
[~,dimGene,NSurrogate] = size(surrogateGeneMatrices);

surrogateGeneCosSimMatrices_dataPIs = zeros([NSurrogate,dimGene,dimPI]);
surrogateGeneCosSimMatrices_dataPIt = zeros([NSurrogate,dimGene,dimPI]);

wiringPIdata = mouseConnectomeDataAnalysis.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
wPIMat_s =  wiringPIdata.WiringPISource(:,1:dimPI);
wPIMat_t =  wiringPIdata.WiringPITarget(:,1:dimPI);



for dg = 1:dimGene
    surroGeneMat = squeeze(surrogateGeneMatrices(:,dg,:));
    for nsr = 1:NSurrogate
        surroGeneVec = surroGeneMat(:,nsr);
        for dpi = 1:dimPI
            wPIVec_s = wPIMat_s(:,dpi);
            wPIVec_t = wPIMat_t(:,dpi);
            cosSim_s = corr(surroGeneVec,wPIVec_s);
            cosSim_t = corr(surroGeneVec,wPIVec_t);
            surrogateGeneCosSimMatrices_dataPIs(nsr,dg,dpi) = cosSim_s;
            surrogateGeneCosSimMatrices_dataPIt(nsr,dg,dpi) = cosSim_t;
        end
    end
end



%% Calculate P-value for each gene
geneCosSimPvalueMat_s = zeros(dimGene,dimPI);
geneCosSimPvalueMat_t = zeros(dimGene,dimPI);
for dg = 1:dimGene
    for dpi = 1:dimPI
        % source
        empiricalDist = squeeze(surrogateGeneCosSimMatrices_dataPIs(:,dg,dpi));
        cosSim_data = relatedGeneAnalysisResults.CorrCoefMatrixSource(dg,dpi);
        pval = empiricalPvalue_abs(empiricalDist,cosSim_data);
        geneCosSimPvalueMat_s(dg,dpi) = pval;
        % target
        empiricalDist = squeeze(surrogateGeneCosSimMatrices_dataPIt(:,dg,dpi));
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



%% Random Test with Surrogate Genes
rng(6)
nullGenerator_dummy = nullGenerator_1;
nullBatchRunner_GeneSurrogate = NullModelAnalysisBatchRunner(nullGenerator_dummy,numRandomConnectome,holdoutParameters_nulls);
[nullBatchRunner_GeneSurrogate,surrogateWPIMatrices_s,surrogateWPIMatrices_t] = nullBatchRunner_GeneSurrogate.run_geneSurrogateNullModel(surrogateGeneMatrices(:,:,:),mouseConnectomeDataAnalysis);













%%






