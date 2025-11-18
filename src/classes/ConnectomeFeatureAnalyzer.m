classdef ConnectomeFeatureAnalyzer

    methods (Static)


        function [reciprocalFeatureMasks] = makeReciprocalUnidirectionalMasks(connMat,crossInfo)
            arguments(Input)
                connMat ConnectionMatrix
                crossInfo CrossRegionInformation
            end
            arguments(Output)
                reciprocalFeatureMasks
            end
            binariMat = connMat.Matrix .* crossInfo.DomainDefineMatrix;
            symmetrical = binariMat + binariMat';
            asymmetrical = binariMat - binariMat';
            reciprocal = (symmetrical == 2);
            unidirectional = (asymmetrical == 1);
            unidirCounter = (asymmetrical == -1);
            others = logical((symmetrical == 0) .* crossInfo.DomainDefineMatrix);
            reciprocalFeatureMasks.reciprocal = reciprocal;
            reciprocalFeatureMasks.unidirectional = unidirectional;
            reciprocalFeatureMasks.unidirectionalCounter = unidirCounter;
            reciprocalFeatureMasks.others = others;
        end

    end

end
