classdef ReconstructionResults
    % ReconstructionResults
    % Computes connection reconstruction accuracy using wiring PI.

    properties
        Parameters ReconstructionParameters
        PIDistanceMatrix        % Matrix of |PIs(xs) - PIt(xt)| differences
        Thresholds            % Threshold values used
        TPR                   % True positive rate per threshold
        FPR                   % False positive rate per threshold
        Precision             % Precision per threshold
        Recall                % Recall per threshold (same as TPR)
        AUC_ROC               % Area under ROC curve
        AUC_PR                % Area under PR curve
        ROCPoints             % [FPR, TPR] pairs
        PRPoints              % [Recall, Precision] pairs
        PredictedMatrices     % Cell array of predicted binary connection matrices per threshold (optional)
    end

    methods
        function obj = ReconstructionResults(reconstructionParameters)
            arguments
                reconstructionParameters ReconstructionParameters
            end
            obj.Parameters = reconstructionParameters;
        end

        function obj = compute(obj,wiringPIPairs,trueConnMat,domainMask,options)
            arguments
                obj ReconstructionResults
                wiringPIPairs WiringPIPairs
                trueConnMat (:,:) logical
                domainMask (:,:) logical
                options.StoreMatrixTag (1,1) logical = 0;
                options.ArbitralPIDistanceMatrix = [];
            end
            dimPI = obj.Parameters.DimPI;
            thresholds = obj.Parameters.getThresholds();
            if ~isempty(options.ArbitralPIDistanceMatrix)
                obj.PIDistanceMatrix = options.ArbitralPIDistanceMatrix;
            else
                sourcePI = wiringPIPairs.WiringPISource;
                targetPI = wiringPIPairs.WiringPITarget;
                obj.PIDistanceMatrix = obj.computePIDistanceMatrix(sourcePI, targetPI, dimPI);
            end
            [TPRs, FPRs, Precisions, Recalls, predMats] = obj.computeROC(obj.PIDistanceMatrix, thresholds, trueConnMat, domainMask);

            obj.TPR = TPRs;
            obj.FPR = FPRs;
            obj.Precision = Precisions;
            obj.Recall = Recalls;
            obj.Thresholds = thresholds;
            % ROC
            [auc_ROC, FPRSorted, TPRSorted] = obj.computeAUC(FPRs,TPRs);
            obj.AUC_ROC = auc_ROC;
            obj.ROCPoints = [FPRSorted, TPRSorted];
            % Precision-Recall
            [auc_PR, RecallSorted, PrecisionSorted] = obj.computeAUC(Recalls,Precisions);
            obj.AUC_PR = auc_PR;
            obj.PRPoints = [RecallSorted, PrecisionSorted];
            % optional: save predicted matrix
            if options.StoreMatrixTag
                obj.PredictedMatrices = predMats;
            else
                obj.PredictedMatrices = {};
            end
        end
        
        function hAx = plotROC(obj, options)
            arguments
                obj ReconstructionResults
                options.ParentAxes = [];
                options.PlotThresholds = [];
                options.PlotColor = [0 0.4470 0.7410];
                options.LineWidth = 1;
                options.LineStyle = "-";
                options.ScatterMarker = ".";
                options.ScatterColor = [0.8500 0.3250 0.0980];
                options.ScatterSize = 24;
                options.TitleOff = false;
                options.ReturnAx = 0;
            end
            % Get axes
            hAx = getOrCreateAxes(options.ParentAxes);
            % plot ROC
            ax = plot(hAx,obj.FPR, obj.TPR, 'Color',options.PlotColor,'LineStyle',options.LineStyle, ...
                'LineWidth',options.LineWidth,'Marker',"none");
            xlabel(hAx, 'FPR');
            ylabel(hAx, 'TPR');
            if options.TitleOff == false
                title(hAx, sprintf('ROC Curve (AUC = %.3f)', obj.AUC_ROC));
            end
            box on
            axis square
            xlim([0,1])
            ylim([0,1])
            xticks([0,0.5,1])
            yticks([0,0.5,1])
            if ~isempty(options.PlotThresholds)
                thresholdIndices = obj.getThresholdIndices(options.PlotThresholds);
                plotFPR = obj.FPR(thresholdIndices);
                plotTPR = obj.TPR(thresholdIndices);
                hold on
                scatter(hAx,plotFPR,plotTPR,options.ScatterSize,options.ScatterColor,'filled');
            end
            if options.ReturnAx == 1
                hAx = ax;
            end
        end

        function plotPR(obj, options)
            arguments
                obj ReconstructionResults
                options.ParentAxes = [];
                options.PlotThresholds = [];
                options.PlotColor = [0 0.4470 0.7410];
                options.LineWidth = 1;
                options.LineStyle = "-";
                options.ScatterMarker = ".";
                options.ScatterColor = [0.8500 0.3250 0.0980];
                options.ScatterSize = 24;
            end
            % Get axes
            hAx = getOrCreateAxes(options.ParentAxes);
            % plot ROC
            plot(hAx,obj.Recall, obj.Precision, 'Color',options.PlotColor,'LineStyle',options.LineStyle, ...
                'LineWidth',options.LineWidth,'Marker',"none");
            xlabel(hAx, 'Recall');
            ylabel(hAx, 'Precision');
            title(hAx, sprintf('ROC Curve (AUC = %.3f)', obj.AUC_ROC));
            box on
            axis square
            xlim([0,1])
            ylim([0,1])
            xticks([0,0.5,1])
            yticks([0,0.5,1])
            if ~isempty(options.PlotThresholds)
                thresholdIndices = obj.getThresholdIndices(options.PlotThresholds);
                plotRecall = obj.Recall(thresholdIndices);
                plotPrecision = obj.Precision(thresholdIndices);
                hold on
                scatter(hAx,plotRecall,plotPrecision,options.ScatterSize,options.ScatterColor,'filled');
            end
        end

        function threshholdIndices = getThresholdIndices(obj,threshholdsList)
            arguments
                obj ReconstructionResults
                threshholdsList 
            end
            nThreshold = numel(threshholdsList);
            threshholdIndices = zeros(size(threshholdsList));
            for n = 1:nThreshold
                id = find(obj.Thresholds == threshholdsList(n));
                if isempty(id)
                    error("the selected threshold was not found")
                end
                threshholdIndices(n) = id;
            end
        end
    end

    methods (Static)
        function diffMatrix = computePIDistanceMatrix(sourcePI, targetPI, D)
            diffMatrix = zeros(height(sourcePI),height(targetPI));
            for d = 1:D
                diffMatrix = diffMatrix + abs(sourcePI(:, d) - targetPI(:, d)');
            end
            diffMatrix = diffMatrix / max(diffMatrix,[],'all');
        end

        function [TPRs, FPRs, Precisions, Recalls, predMats] = computeROC(piDiffMatrix, thresholds, trueConnMatrix, domainMask)
            arguments
                piDiffMatrix (:,:) double
                thresholds (:,1) double
                trueConnMatrix (:,:) logical
                domainMask (:,:) logical
            end
            trueFlat = trueConnMatrix(domainMask);
            numThresholds = numel(thresholds);
            TPRs = zeros(1, numThresholds);
            FPRs = zeros(1, numThresholds);
            Precisions = zeros(1, numThresholds);
            Recalls = zeros(1, numThresholds);
            predMats = cell(1, numThresholds);

            for i = 1:numThresholds
                th = thresholds(i);
                predMatrix = piDiffMatrix <= th;
                predFlat = predMatrix(domainMask);

                TP = sum(predFlat & trueFlat);
                FP = sum(predFlat & ~trueFlat);
                FN = sum(~predFlat & trueFlat);
                TN = sum(~predFlat & ~trueFlat);

                TPRs(i) = TP / (TP + FN);
                FPRs(i) = FP / (FP + TN);
                if TP + FP == 0
                    eps = 0.001;
                else 
                    eps = 0;
                end
                Precisions(i) = TP / (TP + FP + eps);  % eps to avoid division by zero
                Recalls(i) = TPRs(i);
                predMats{i} = predMatrix .* domainMask;
            end
        end

        function [auc, xSorted, ySorted] = computeAUC(x, y)
            [xSorted, idx] = sort(x);
            ySorted = y(idx);
            % Ensure left end (0)
            if xSorted(1) > 0
                xSorted = [0, xSorted];
                ySorted = [ySorted(1), ySorted];
            end
            % Ensure right end (1)
            if xSorted(end) < 1
                xSorted = [xSorted, 1];
                ySorted = [ySorted, ySorted(end)];
            end
            auc = trapz(xSorted, ySorted);
        end

       

    end
end
