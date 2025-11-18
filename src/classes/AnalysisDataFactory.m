classdef AnalysisDataFactory
    properties
        ConnectionMatrix ConnectionMatrix         
        GeneExpressionLevels GeneExpressionLevels   
        CrossRegionInformation CrossRegionInformation  

        Description string = ""   % explanation
    end

    methods
        function obj = AnalysisDataFactory(connMatrix, geneExprLevels, crossRegionInfo, description)
            arguments
                connMatrix ConnectionMatrix
                geneExprLevels GeneExpressionLevels
                crossRegionInfo CrossRegionInformation
                description string = "";
            end

            obj.ConnectionMatrix = connMatrix;
            obj.GeneExpressionLevels = geneExprLevels;
            obj.CrossRegionInformation = crossRegionInfo;
            obj.Description = description;
        end

        function pairDataSet = generatePairData(obj, regionMask, pairDataGenerationOptions,options)
            % make connection index pairs with region mask
            arguments
                obj AnalysisDataFactory
                regionMask % [] or (:,:) logical
                pairDataGenerationOptions PairDataGenerationOptions
                options.Description string = "";
            end
            % 1. get source/target index
            if isempty(regionMask)
                regionMask = obj.CrossRegionInformation.DomainDefineMatrix;
            end
            maskedConnMat = obj.ConnectionMatrix.Matrix .* regionMask;
            [srcIdx, tgtIdx] = find(maskedConnMat);  % [source,target] index for each pair
            % 2. extract gene expression levels
            switch pairDataGenerationOptions.DataType
                case "PC"
                    expMatrix = obj.GeneExpressionLevels.PCAMatrix(:,1:pairDataGenerationOptions.DataDimension);
                case "gene"
                    expMatrix = obj.GeneExpressionLevels.ExpressionMatrix;
                    pairDataGenerationOptions.DataDimension = width(expMatrix);
            end
            srcExpr = expMatrix(srcIdx,:);
            tgtExpr = expMatrix(tgtIdx,:);
            % 3. add noise
            if ~isempty(pairDataGenerationOptions.RngSeed)
                rng(pairDataGenerationOptions.RngSeed);  % Optional random seed
                disp("rng is used for pair data generation")
            end
            noiseLevel = pairDataGenerationOptions.NoiseLevel;
            srcExpr = srcExpr + noiseLevel * randn(size(srcExpr));
            tgtExpr = tgtExpr + noiseLevel * randn(size(tgtExpr));
            % 4. record GenerationInfo
            genInfo = struct( ...
                'FactoryClass', class(obj), ...
                'PairDataGenerationOptions', pairDataGenerationOptions ...
            );
            % 5. make ConnectomeGenePairDataSet object
            pairDataSet = ConnectomeGenePairDataSet( ...
                srcIdx,tgtIdx,srcExpr,tgtExpr,regionMask, ...
                GenerationInfo = genInfo, ...
                Description = options.Description);
        end

        function [fullPairData] = generateFullPairData(obj,pairDataGenerationOptions,options);
            arguments
                obj AnalysisDataFactory
                pairDataGenerationOptions PairDataGenerationOptions
                options.Description = "";
            end
            fullPairData = obj.generatePairData(obj.CrossRegionInformation.DomainDefineMatrix,pairDataGenerationOptions, ...
                "Description",options.Description);
        end

        function [trainSet, testSet] = generateHoldoutPairData(obj, holdoutSplit, pairDataGenerationOptions, options)
            % generate Train/Test pair dataset using HoldourtSplit
            arguments
                obj
                holdoutSplit HoldoutSplit
                pairDataGenerationOptions PairDataGenerationOptions;
                options.Description string = ""
            end

            trainSet = obj.generatePairData(holdoutSplit.TrainMask, pairDataGenerationOptions, ...
                Description = "Train - " + holdoutSplit.Description);

            testSet = obj.generatePairData(holdoutSplit.TestMask, pairDataGenerationOptions, ...
                Description = "Test - " + holdoutSplit.Description);
        end

        function [srcIndList,tgtIndList] = getRegionIndPair(obj, regionMask, pairDataGenerationOptions)
            arguments
                obj AnalysisDataFactory
                regionMask % [] or (:,:) logical
                pairDataGenerationOptions PairDataGenerationOptions
            end
            pairDataSet = obj.generatePairData(regionMask, pairDataGenerationOptions);
            srcIndList = pairDataSet.SourceIDs;
            tgtIndList = pairDataSet.TargetIDs;
        end

        function [speceficMRPairListMask,srcMask,tgtMask] = getSpecificMRPairMask(obj,srcMRName,tgtMRName,pairDataGenerationOptions)
            arguments
                obj AnalysisDataFactory
                srcMRName string
                tgtMRName string
                pairDataGenerationOptions PairDataGenerationOptions
            end
            [srcIndList,tgtIndList] = obj.getRegionIndPair([],pairDataGenerationOptions);
            crossInfo = obj.CrossRegionInformation;
            matchingRegionIndicesSource = find(crossInfo.SourceRegionInfo.RegionTable.major_region == srcMRName);
            matchingRegionIndicesTarget = find(crossInfo.SourceRegionInfo.RegionTable.major_region == tgtMRName);
            srcMask = ismember(srcIndList,matchingRegionIndicesSource);
            tgtMask = ismember(tgtIndList,matchingRegionIndicesTarget);
            speceficMRPairListMask = srcMask .* tgtMask;
        end

    end
end