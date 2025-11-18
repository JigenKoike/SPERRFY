classdef MajorRegionBlockRandomizer
    % Applies randomization method to each MajorRegion pair block, or to the whole matrix

    properties 
        RandomizationMethod function_handle
        UseBlockwiseLogical (1,1) logical
    end

    methods
        function obj = MajorRegionBlockRandomizer(randomizationMethod, useBlockwise)
            arguments          
                randomizationMethod function_handle
                useBlockwise logical
            end     
            obj.RandomizationMethod = randomizationMethod;
            obj.UseBlockwiseLogical = useBlockwise;
        end

        function [randomizedMatrix,newDomDef] = apply(obj,connMat,crossInfo,options)
            arguments 
                obj MajorRegionBlockRandomizer 
                connMat ConnectionMatrix
                crossInfo CrossRegionInformation
                options.DistanceBinWidth = [];
            end    
            origMat = connMat.Matrix;
            domDefMat = crossInfo.DomainDefineMatrix;
            distMat = crossInfo.DistanceMatrix;
            newMat = zeros(size(origMat));
            newDomDef = zeros(size(domDefMat));
            if isequal(obj.RandomizationMethod,@ConnectomeRandomizationMethods.networkPreserveGeneration)
                if obj.UseBlockwiseLogical
                    MRCategoricalSource = crossInfo.SourceRegionInfo.getMajorRegions;
                    MRCategoricalTarget = crossInfo.TargetRegionInfo.getMajorRegions;
                    [newMat,newDomDef] = ConnectomeRandomizationMethods.networkPreserveGeneration(origMat,domDefMat,[], ...
                        "MRCategoricalSrc",MRCategoricalSource,"MRCategoricalTgt",MRCategoricalTarget,"SameRegionTag",true);
                else
                    [newMat,newDomDef] = ConnectomeRandomizationMethods.networkPreserveGeneration(origMat,domDefMat,[]);
                end
            else
                if obj.UseBlockwiseLogical
                    sourceMRList = unique(crossInfo.SourceRegionInfo.getMajorRegions,"stable");
                    targetMRList = unique(crossInfo.TargetRegionInfo.getMajorRegions,"stable");
                    for i = 1:numel(sourceMRList)
                        srcMajor = sourceMRList(i,1);
                        for j = 1:numel(targetMRList)
                            tgtMajor = targetMRList(j,1);
        
                            srcIdx = find(crossInfo.SourceRegionInfo.getMajorRegions == srcMajor);
                            tgtIdx = find(crossInfo.TargetRegionInfo.getMajorRegions == tgtMajor);
        
                            subMat = origMat(srcIdx, tgtIdx);
                            subDomDef = domDefMat(srcIdx,tgtIdx);
        
                            if ~isempty(distMat)
                                subDist = distMat(srcIdx, tgtIdx);
                            else
                                subDist = [];
                            end
                            if isempty(options.DistanceBinWidth) == 0
                                [randomizedBlock,domDefBlock] = ConnectomeRandomizationMethods.distancePreservedGeneration(subMat, subDomDef, subDist,"BinWidth",options.DistanceBinWidth);
                            else
                                [randomizedBlock,domDefBlock] = obj.RandomizationMethod(subMat, subDomDef, subDist);
                            end
                            newMat(srcIdx, tgtIdx) = randomizedBlock;
                            newDomDef(srcIdx,tgtIdx) = domDefBlock;
                        end
                    end
                else
                    [randomizedAll,domDefAll] = obj.RandomizationMethod(origMat, domDefMat, distMat);
                    newMat = randomizedAll;
                    newDomDef = domDefAll;
                end
            end
            randomizedMatrix = newMat;
        end
    end
end
