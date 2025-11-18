classdef CrossRegionInformation < handle
    % CrossRegionInformation
    % Manages information between source and target brain regions, such as distance and domain definition.

    properties
        SourceRegionInfo BrainRegionInformation
        TargetRegionInfo BrainRegionInformation
        DistanceMatrix double = []            % [N_source x N_target] Euclidean distances between regions
        DomainDefineMatrix logical = []       % [N_source x N_target] Mask of meaningful region pairs
    end

    methods
        function obj = CrossRegionInformation(srcInfo, tgtInfo)
            arguments
                srcInfo BrainRegionInformation
                tgtInfo BrainRegionInformation
            end
            obj.SourceRegionInfo = srcInfo;
            obj.TargetRegionInfo = tgtInfo;
        end

        function setDistanceMatrix(obj, D)
            arguments
                obj CrossRegionInformation
                D double
            end
            obj.DistanceMatrix = D;
        end
        
        function setDomainDefineMatrix(obj, M)
            arguments
                obj CrossRegionInformation
                M logical
            end
            obj.DomainDefineMatrix = M;
        end

        function copyObj = copyCrossRegionInfo(obj)
            arguments
                obj CrossRegionInformation
            end
            copyObj = CrossRegionInformation(obj.SourceRegionInfo,obj.TargetRegionInfo);
            copyObj.setDomainDefineMatrix(obj.DomainDefineMatrix);
            copyObj.setDistanceMatrix(obj.DistanceMatrix);
        end

        function M = makeDifferentMajorRegionDomainDefineMatrix(obj)
            % Define domain by matching Major Region between source and target
            srcMajor = obj.SourceRegionInfo.getMajorRegions();
            tgtMajor = obj.TargetRegionInfo.getMajorRegions();
            Nsrc = numel(srcMajor);
            Ntgt = numel(tgtMajor);
            M = false(Nsrc, Ntgt);
            for i = 1:Nsrc
                for j = 1:Ntgt
                    M(i,j) = strcmp(srcMajor{i}, tgtMajor{j});
                end
            end
            M = ~M;
        end

        function M = makeDomainDefineMatrixByDistance(obj, threshold)
            % Define domain by distance threshold (only regions closer than threshold are valid)
            if isempty(obj.DistanceMatrix)
                error("DistanceMatrix is empty")
            end
            M = obj.DistanceMatrix <= threshold;
        end

        function M = getDomainDefineMatrix(obj)
            % Getter for domain mask
            M = obj.DomainDefineMatrix;
        end

        function D = getDistanceMatrix(obj)
            % Getter for distance matrix
            D = obj.DistanceMatrix;
        end

        function mask = getDomainDefineMask(obj)
            domDefMat = obj.DomainDefineMatrix;
            mask = ones(size(domDefMat));
            mask(~domDefMat) = NaN;
        end

        function maskedMatrix = makeMaskedMatrix(obj,matrix)
            if isempty(obj.DomainDefineMatrix) == 0
                domDefask = obj.getDomainDefineMask;
                maskedMatrix = matrix .* domDefask;
            else
                maskedMatrix = matrix;
            end
        end

        function [trainMask, testMask] = generateRandomHoldoutMasks(obj, trainRatio, rngSeed)
            % Generates random train/test masks based on DomainDefineMatrix.
            % Inputs:
            %   - trainRatio: proportion of valid pairs to include in the training set (0~1)
            %   - rngSeed: optional random seed for reproducibility
            % Outputs:
            %   - trainMask: logical matrix (same size as DomainDefineMatrix)
            %   - testMask: logical matrix (same size as DomainDefineMatrix)
            arguments
                obj CrossRegionInformation
                trainRatio (1,1) {mustBeNumeric} 
                rngSeed = [];
            end
            domainMask = obj.DomainDefineMatrix;
            [rowIdx, colIdx] = find(domainMask);    
            numValid = numel(rowIdx);
            numTrain = round(trainRatio * numValid);
            % Optional random seed
            if ~isempty(rngSeed)
                rng(rngSeed);
                disp("rng is used for making holdout mask")
            end
            % Randomly permute indices
            permIdx = randperm(numValid);
            % Split indices
            trainIdx = permIdx(1:numTrain);
            testIdx = permIdx(numTrain+1:end);
            % Initialize masks
            sz = size(domainMask);
            trainMask = false(sz);
            testMask = false(sz);
            % Set true for selected indices
            linIdx = sub2ind(sz, rowIdx, colIdx);
            trainMask(linIdx(trainIdx)) = true;
            testMask(linIdx(testIdx)) = true;
        end


        function [MRPairMask] = makeSpecificMRPairMask(obj, sourceMRName, targetMRName)
            arguments(Input)
                obj CrossRegionInformation
                sourceMRName string
                targetMRName string
            end
            sourceMRLogical = (obj.SourceRegionInfo.RegionTable.major_region == sourceMRName);
            targetMRLogical = (obj.TargetRegionInfo.RegionTable.major_region == targetMRName);
            MRPairMask = sourceMRLogical * targetMRLogical' ;           
        end



    end
end
