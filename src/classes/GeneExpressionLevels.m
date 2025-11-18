classdef GeneExpressionLevels < handle
    % GeneExpressionLevels (PCA-specialized)
    % Manages brain region gene expression levels and their PCA representation.

    properties
        ExpressionMatrix double          % [N_regions x N_genes] original expression matrix
        AnnotationIndices double         % [N_regions x 1] annotation indices for regions
        GeneAcronyms string              % [1 x N_genes] gene acronyms
        GeneNames string                 % [1 x N_genes] gene full names

        PCAMatrix double = []            % [N_regions x N_components] PCA-transformed data
        PCAInfo PCAResults = []       % Struct holding PCA details: Loadings, ExplainedVariance, MeanExpression
        ClusteringResults struct = []
    end

    methods
        function obj = GeneExpressionLevels(exprMatrix, annotationIdx, geneAcronyms, geneNames)
            arguments
                exprMatrix double
                annotationIdx double
                geneAcronyms string
                geneNames string
            end

            obj.ExpressionMatrix = exprMatrix;
            obj.AnnotationIndices = annotationIdx;
            obj.GeneAcronyms = geneAcronyms;
            obj.GeneNames = geneNames;
        end

        function performGeneClustering(obj,options)
            arguments
                obj GeneExpressionLevels
                options.pdistDistance = 'correlation';
                options.linkageMethod = 'average';
            end
            geneExpressionMatrix = obj.ExpressionMatrix;
            data_T = geneExpressionMatrix';
            D = pdist(data_T,options.pdistDistance);
            Z = linkage(D,options.linkageMethod);
            leafOrder = optimalleaforder(Z, D);
            sortedGeneExpression = geneExpressionMatrix(:, leafOrder);
            obj.ClusteringResults = struct("SortedGeneExpression",sortedGeneExpression, ...
                "D",D,"Z",Z,"LeafOrder",leafOrder,"Options",options);
        end

        function performPCA(obj)
            arguments
                obj GeneExpressionLevels
            end
            pcaResults = PCAResults(obj.ExpressionMatrix);

            obj.PCAMatrix = pcaResults.score;
            obj.PCAInfo = pcaResults;
        end

        function pcs = getPCs(obj)
            pcs = obj.PCAMatrix;
        end

        function info = getPCAInfo(obj)
            info = obj.PCAInfo;
        end

        function [shuffledExprMatrix] = makeRegionShuffledGeneExprLevels_withinMR(obj,brainInfo)
            arguments
                obj GeneExpressionLevels
                brainInfo BrainRegionInformation
            end
            shuffledRegionIndex = brainInfo.getShuffledRegionIndexWithinMRs();
            shuffledExprMatrix = obj.ExpressionMatrix(shuffledRegionIndex,:);
        end


        

    end
end
