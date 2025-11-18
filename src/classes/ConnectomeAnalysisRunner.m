classdef ConnectomeAnalysisRunner
    % ConnectomeAnalysisRunner
    % Executes full and holdout analyses for a given connectome and shared data.

    properties
        Factory AnalysisDataFactory
        PairGenerationOptions PairDataGenerationOptions
        ReconstructionOptions ReconstructionParameters
        HoldoutOptions HoldoutParameters
        OverallModelAndResults OverallAnalysisGroup
    end

    methods
        function obj = ConnectomeAnalysisRunner(connectomeMatrix, geneExprLevels,crossRegionInfo, pairDataGenerationOptions,reconstParams,holdoutParams,options)
            arguments
                connectomeMatrix ConnectionMatrix
                geneExprLevels GeneExpressionLevels
                crossRegionInfo CrossRegionInformation
                pairDataGenerationOptions PairDataGenerationOptions
                reconstParams ReconstructionParameters
                holdoutParams HoldoutParameters
                options.FactoryDescription = '';
            end
            obj.Factory = AnalysisDataFactory(connectomeMatrix,geneExprLevels,crossRegionInfo,options.FactoryDescription);
            obj.PairGenerationOptions = pairDataGenerationOptions;
            obj.ReconstructionOptions = reconstParams;  
            obj.HoldoutOptions = holdoutParams;
        end

        function obj = runAnalysis(obj,options)
            arguments
                obj ConnectomeAnalysisRunner
                options.ReconstructionStoreTag = false;
            end
            % Full analysis
            fullPairData = obj.Factory.generateFullPairData(obj.PairGenerationOptions);
            fullModel = PairDataCCAAnalysisModel();
            fullModel = fullModel.runCCAAnalysis(fullPairData,obj.Factory.GeneExpressionLevels,obj.PairGenerationOptions);
            % %%% how manage store matrix tag ??
            fullModel = fullModel.runReconstructionAnalysis(obj.ReconstructionOptions, ...
                obj.Factory.ConnectionMatrix.Matrix,obj.Factory.CrossRegionInformation.DomainDefineMatrix,"StoreMatrixTag",options.ReconstructionStoreTag);

            % Holdout analyses
            if ~isempty(obj.HoldoutOptions)
                holdoutUnitArray = obj.performHoldout();
            else
                holdoutUnitArray = [];
            end

            % Overall result
            overallGroup = OverallAnalysisGroup(fullModel);
            overallGroup = overallGroup.addHoldoutUnitArray(holdoutUnitArray);
            overallGroup = overallGroup.computeHoldoutSummary();

            obj.OverallModelAndResults = overallGroup;
        end

        function holdoutUnitArray = performHoldout(obj)
            numSplits = obj.HoldoutOptions.NumSplits;
            holdoutUnitCells = cell(numSplits,1);
            crossInfo = obj.Factory.CrossRegionInformation;
            trainRatio = obj.HoldoutOptions.TrainRatio;
            geneExprLevels = obj.Factory.GeneExpressionLevels;
            trueConnMat = obj.Factory.ConnectionMatrix.Matrix;
            for i = 1:numSplits
                [trainMask,testMask] = generateRandomHoldoutMasks(crossInfo,trainRatio,[]);
                holdoutSplit = HoldoutSplit(trainMask, testMask, "random split");
                [trainSet, testSet] = obj.Factory.generateHoldoutPairData(holdoutSplit,obj.PairGenerationOptions);
                trainModel = PairDataCCAAnalysisModel();
                trainModel = trainModel.runCCAAnalysis(trainSet,geneExprLevels,obj.PairGenerationOptions);
                trainModel = trainModel.runReconstructionAnalysis(obj.ReconstructionOptions,trueConnMat,trainSet.DomainDefineMask,"StoreMatrixTag",false);    

                holdoutUnit = HoldoutAnalysisUnit(trainModel,holdoutSplit);
                holdoutUnit = holdoutUnit.runTestMetrics(testSet,trueConnMat);
                holdoutUnitCells{i,1} = holdoutUnit;
            end
            holdoutUnitArray = [holdoutUnitCells{:}].';
        end

        function [relatedGeneAnalysisResults] = performRelatedGeneAnalysis(obj,KTopGene)
            arguments(Input)
                obj ConnectomeAnalysisRunner
                KTopGene = [];
            end
            arguments(Output)
                relatedGeneAnalysisResults RelatedGeneAnalysisResults
            end
            if isempty(KTopGene) == 1;
                KTopGene = width(obj.Factory.GeneExpressionLevels.ExpressionMatrix);
            end
            wiringPIPairs = obj.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
            geneExprLevels = obj.Factory.GeneExpressionLevels;
            [corrCoeffMatrixSource,corrCoeffMatrixTarget] = obj.calculateWiringPIGeneCorrelations(wiringPIPairs,geneExprLevels);
            [~,topKGeneIndexListsSource] = maxk(corrCoeffMatrixSource,KTopGene,1,"ComparisonMethod","abs");
            [~,topKGeneIndexListsTarget] = maxk(corrCoeffMatrixTarget,KTopGene,1,"ComparisonMethod","abs");
            relatedGeneAnalysisResults = RelatedGeneAnalysisResults(geneExprLevels,corrCoeffMatrixSource,corrCoeffMatrixTarget, ...
                topKGeneIndexListsSource,topKGeneIndexListsTarget);
        end

    end

    methods(Static)
        function [corrCoeffMatSource,corrCoeffMatTarget] = calculateWiringPIGeneCorrelations(wiringPIPairs,geneExprLevels)
            arguments(Input)
                wiringPIPairs WiringPIPairs
                geneExprLevels GeneExpressionLevels
            end
            sourcePIs = wiringPIPairs.WiringPISource;
            targetPIs = wiringPIPairs.WiringPITarget;
            geneExprMatrix = geneExprLevels.ExpressionMatrix;
            DimPI = width(sourcePIs);
            NGene = width(geneExprMatrix);
            if height(sourcePIs) ~= height(geneExprMatrix)
                error("region size is different between wiring PI ande genes")
            end
            corrCoeffMatSource = zeros(NGene,DimPI);
            corrCoeffMatTarget = zeros(NGene,DimPI);
            for ng = 1:NGene
                geneExprVec = geneExprMatrix(:,ng);
                for dpi = 1:DimPI
                    sPIVec = sourcePIs(:,dpi);
                    tPIVec = targetPIs(:,dpi);
                    corrCoeffMatSource(ng,dpi) = corr(geneExprVec,sPIVec);
                    corrCoeffMatTarget(ng,dpi) = corr(geneExprVec,tPIVec);
                end
            end
        end


    end


end
