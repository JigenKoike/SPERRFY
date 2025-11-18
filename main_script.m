% main_script.m

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
relatedGeneAnalysisResults = mouseConnectomeDataAnalysis.performRelatedGeneAnalysis();

%% Perform null connectome analysis

nullGenerator_1 = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"RandomGeneration",0, ...
    "Description","RandomGeneration, Golobally");
nullGenerator_2 = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"RandomGeneration",1, ...
    "Description","RandomGeneration, Locally");
nullGenerator_3 = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"DistancePreserve",0, ...
    "Description","RandomGeneration, Golobally");
nullGenerator_4 = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"DistancePreserve",1, ...
    "Description","RandomGeneration, Locally");
nullGenerator_5 = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"NetworkPreserve",0, ...
    "Description","RandomGeneration, Golobally");
nullGenerator_6 = NullConnectomeGeneratorModel(connectomeMatrix,crossInfo,"NetworkPreserve",1, ...
    "Description","RandomGeneration, Locally");
% sample matrix
rng(2)
randomConnectomeSample_1 = nullGenerator_1.generate;
randomConnectomeSample_2 = nullGenerator_2.generate;
randomConnectomeSample_3 = nullGenerator_3.generate;
randomConnectomeSample_4 = nullGenerator_4.generate;
randomConnectomeSample_5 = nullGenerator_5.generate;
randomConnectomeSample_6 = nullGenerator_6.generate;
% run
rng(3)
nullBatchRunner_1 = NullModelAnalysisBatchRunner(nullGenerator_1,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_1 = nullBatchRunner_1.run(mouseConnectomeDataAnalysis);
nullBatchRunner_2 = NullModelAnalysisBatchRunner(nullGenerator_2,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_2 = nullBatchRunner_2.run(mouseConnectomeDataAnalysis);
nullBatchRunner_3 = NullModelAnalysisBatchRunner(nullGenerator_3,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_3 = nullBatchRunner_3.run(mouseConnectomeDataAnalysis);
nullBatchRunner_4 = NullModelAnalysisBatchRunner(nullGenerator_4,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_4 = nullBatchRunner_4.run(mouseConnectomeDataAnalysis);
nullBatchRunner_5 = NullModelAnalysisBatchRunner(nullGenerator_5,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_5 = nullBatchRunner_5.run(mouseConnectomeDataAnalysis);
nullBatchRunner_6 = NullModelAnalysisBatchRunner(nullGenerator_6,numRandomConnectome,holdoutParameters_nulls);
nullBatchRunner_6 = nullBatchRunner_6.run(mouseConnectomeDataAnalysis);

%%

%% Test for Various PC
PCDimList = [10,25,50,75,100];
L_PCList = length(PCDimList);
rngSeed = 4;
rng(rngSeed)
clear variousPCTestResults

for l = 1:L_PCList
    dimPC = PCDimList(l);
    if dimPC ~= 50
        variousPC_pairDataGeneratioOptions = PairDataGenerationOptions(noiseLevel,dataType,dimPC);
        variousPC_AnalysisModel = ConnectomeAnalysisRunner(connectomeMatrix, geneExprLevels,crossInfo, ...
            variousPC_pairDataGeneratioOptions,reconstructionParameters,holdoutParameters_main,"FactoryDescription","Various PC");
        variousPC_AnalysisModel = variousPC_AnalysisModel.runAnalysis;
        variousPCTestResults(l) = variousPC_AnalysisModel.OverallModelAndResults;
    elseif dimPC == 50
        variousPCTestResults(l) = mouseConnectomeDataAnalysis.OverallModelAndResults;
    end
end

%%
main_additional_analysis;

%% Save data
cd(strcat(projectRoot,"/results/calculatedVariables"));
dt = datetime('now','Format','yyyyMMdd_HHmmss');
dtstr = string(dt);
flName = strcat("results_",dtstr,".mat");
save(flName, '-v7.3');
cd(projectRoot);


%% Make Figures

% main_makeAllFigures;