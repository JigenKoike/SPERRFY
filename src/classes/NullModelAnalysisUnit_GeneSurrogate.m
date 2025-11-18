classdef NullModelAnalysisUnit_GeneSurrogate
    % Holds and executes a single null model analysis, leveraging ConnectomeAnalysisRunner

    properties
        AnalysisRunner ConnectomeAnalysisRunner
        WiringPISimilarity
    end

    methods
        function obj = NullModelAnalysisUnit_GeneSurrogate(surrogateGeneMatrix,randHoldoutOptions,origRunner)
            arguments
                surrogateGeneMatrix (:,:) double
                randHoldoutOptions HoldoutParameters
                origRunner ConnectomeAnalysisRunner
            end
            origFactory = origRunner.Factory;
            origGeneExpr = origFactory.GeneExpressionLevels;
            origConnMat = ConnectionMatrix(origFactory.ConnectionMatrix.Matrix);
            origCrossInfo = copyCrossRegionInfo(origFactory.CrossRegionInformation);
            newGeneExprLevels = GeneExpressionLevels(surrogateGeneMatrix, ...
                origGeneExpr.AnnotationIndices,origGeneExpr.GeneAcronyms,origGeneExpr.GeneNames);
            newGeneExprLevels.performPCA;
            newGeneExprLevels.performGeneClustering;
            nullConnRunner = ConnectomeAnalysisRunner(origConnMat,newGeneExprLevels,origCrossInfo, ...
                origRunner.PairGenerationOptions,origRunner.ReconstructionOptions,randHoldoutOptions);
            obj.AnalysisRunner = nullConnRunner;            
        end
        
        %
        function obj = run(obj,options)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
                options.RngSeed = [];
            end
            if ~isempty(options.RngSeed)
                rng(options.RngSeed)
                disp("rng is used for running null model analysis unit")
            end
            % Set up and execute analysis using runner
            obj.AnalysisRunner = obj.AnalysisRunner.runAnalysis("ReconstructionStoreTag",false);
        end

        function obj = calculatePISimilarityScore(obj,refWiringPIPairs)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
                refWiringPIPairs WiringPIPairs
            end
            unitWiringPI = obj.AnalysisRunner.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
            simScoreList = calculatePISimilarities(unitWiringPI,refWiringPIPairs);
            obj.WiringPISimilarity = simScoreList;
        end

        function fullDataCorrs = getFullDataCorrs(obj)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
            end
            fullDataCorrs = obj.AnalysisRunner.OverallModelAndResults.FullAnalysisModel.CCAResultsData.r;
        end

        function fullDataAUCs = getFullDataAUC(obj)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
            end
            fullDataAUCs = obj.AnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.AUC_ROC;
        end

        function holdoutTrainCorrs = getHoldoutTrainMeanCorrs(obj)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
            end
            holdoutTrainCorrs = obj.AnalysisRunner.OverallModelAndResults.HoldoutSummary.TrainMeanCorrelation;
        end

        function holdoutTrainAUCs = getHoldoutTrainMeanAUC(obj)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
            end
            holdoutTrainAUCs = obj.AnalysisRunner.OverallModelAndResults.HoldoutSummary.TrainMeanAUC_ROC;
        end

        function holdoutTestCorrs = getHoldoutTestMeanCorrs(obj)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
            end
            holdoutTestCorrs = obj.AnalysisRunner.OverallModelAndResults.HoldoutSummary.TestMeanCorrelation;
        end

        function holdoutTestAUCs = getHoldoutTestMeanAUC(obj)
            arguments
                obj NullModelAnalysisUnit_GeneSurrogate
            end
            holdoutTestAUCs = obj.AnalysisRunner.OverallModelAndResults.HoldoutSummary.TestMeanAUC_ROC;
        end

        %}
        



    end
end