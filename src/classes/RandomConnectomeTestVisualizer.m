classdef RandomConnectomeTestVisualizer
    methods (Static)

        function [hAx,hPl] = plotRandomTestCorrelations(nullBatchRunner,options)
            arguments
                nullBatchRunner NullModelAnalysisBatchRunner
                options.ParentAxes = [];
                options.StdErrorBar = "off"; % or "on"
                options.DataType = "Full"; % or "Train" or "Test"
                options.Color = [0 0.4470 0.7410];
                options.MarkerSize = 6;
                options.LineWidth = 1;
            end
            % Get axes
            hAx = getOrCreateAxes(options.ParentAxes);
            % data select
            ylimVec = [-0.1,1.0001];
            switch options.DataType
                case "Full"
                    plotDataList = nullBatchRunner.FullDataCorrs;
                case "Train"
                    plotDataList = nullBatchRunner.HoldoutTrainCorrMeans;
                case "Test"
                    plotDataList = nullBatchRunner.HoldoutTestCorrMeans;
                case "PISimilarity"
                    plotDataList = nullBatchRunner.WirigPISimilarities;
                    ylimVec = [0,1.0001];
            end
            meanList = mean(plotDataList,1);
            stdList = std(plotDataList,1);
            switch options.StdErrorBar
                case "off"
                    hPl = plot(meanList,'.-','MarkerSize',options.MarkerSize,'Color',options.Color,'LineWidth',options.LineWidth);
                case "on"
                    capSz = 1;
                    hPl = errorbar(1:width(stdList),meanList,stdList,'Color',options.Color,'LineWidth',options.LineWidth,'CapSize',capSz);
                    hPl.LineWidth = 1;
                    hold on
                    hPl = plot(meanList,'-','Color',options.Color,'LineWidth',options.LineWidth);
            end
            ylim(ylimVec)
            xlim([0,width(meanList)])
            yticks([0:0.2:1])
            box off
        end

        function [hAx,hPl] = swarmRandomTestCorrelations(nullBatchRunner,DimShow,options)
            arguments
                nullBatchRunner NullModelAnalysisBatchRunner
                DimShow (1,1) {mustBeInteger}
                options.ParentAxes = [];
                options.DataType = "Full"; % or "Train" or "Test"
                options.Color = [0 0.4470 0.7410];
                options.SwarmSize = 1.5;
            end
            % Get axes
            hAx = getOrCreateAxes(options.ParentAxes);
            % data select
            ylimVec = [-0.1,1.0001];
            switch options.DataType
                case "Full"
                    plotDataList = nullBatchRunner.FullDataCorrs(:,1:DimShow);
                case "Train"
                    plotDataList = nullBatchRunner.HoldoutTrainCorrMeans(:,1:DimShow);
                case "Test"
                    plotDataList = nullBatchRunner.HoldoutTestCorrMeans(:,1:DimShow);
                case "PISimilarity"
                    plotDataList = nullBatchRunner.WirigPISimilarities(:,1:DimShow);
                    ylimVec = [0,1.0001];
            end
            x = ones(height(plotDataList),1) * (1:DimShow);
            hPl = swarmchart(x,plotDataList,options.SwarmSize,options.Color,'filled');
            ylim(ylimVec)
            xlim([0.5,DimShow+0.5])
            xticks([1:DimShow])
            yticks([0:0.2:1])
            box off
        end

        function [hAx] = compareCorrPlot(connAnalysisRunner,globalRunner,localRunner,options)
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                globalRunner NullModelAnalysisBatchRunner
                localRunner NullModelAnalysisBatchRunner
                options.DataType = "Full"; % or "Train" or "Test"
                options.StdErrorBar = "on";
                options.MarkerSize = 6;
                options.LineWidth = 1;
                options.ParentAxes = [];
                options.ColorList = [];
            end
            % Get axes
            hAx = getOrCreateAxes(options.ParentAxes);
            if isempty(options.ColorList)
                colorList = orderedcolors("gem");
            else
                colorList = options.ColorList;
            end
            switch options.DataType
                case "Full"
                    plotDataOriginal = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.CCAResultsData.r;
                case "Train"
                    plotDataOriginal = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TrainMeanCorrelation;
                case "Test"
                    plotDataOriginal = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestMeanCorrelation;
            end
            h(1) = plot(plotDataOriginal,'.-','Color',colorList(1,:),"MarkerSize",options.MarkerSize,"LineWidth",options.LineWidth);
            hold on
            [~,h(2)] = RandomConnectomeTestVisualizer.plotRandomTestCorrelations(globalRunner, ...
                "Color",colorList(2,:),"DataType",options.DataType,"StdErrorBar",options.StdErrorBar,"ParentAxes",hAx, ...
                "MarkerSize",options.MarkerSize,"LineWidth",options.LineWidth);
            hold on
            [~,h(3)] = RandomConnectomeTestVisualizer.plotRandomTestCorrelations(localRunner, ...
                "Color",colorList(3,:),"DataType",options.DataType,"StdErrorBar",options.StdErrorBar,"ParentAxes",hAx, ...
                "MarkerSize",options.MarkerSize,"LineWidth",options.LineWidth);
            box off
            ylim([-0.1,1.0001])
            ylabel("Correlation Coefficient")
            xlabel("Correlation Component (Rank)")
            legendLabels = ["Original connectome data",strcat("Globally randomized model"),strcat("Locally randomized model")];
            legend(h,legendLabels);
        end

        function [hAx] = compareCorrSwarm(connAnalysisRunner,globalRunner,localRunner,DimShow,options)
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                globalRunner NullModelAnalysisBatchRunner
                localRunner NullModelAnalysisBatchRunner
                DimShow (:,:) {mustBeInteger}
                options.DataType = "Full"; % or "Train" or "Test"
                options.SwarmSize = 6;
                options.ParentAxes = [];
                options.ColorList = [];
            end
            % Get axes
            hAx = getOrCreateAxes(options.ParentAxes);
            if isempty(options.ColorList)
                colorList = orderedcolors("gem");
            else
                colorList = options.ColorList;
            end
            switch options.DataType
                case "Full"
                    plotDataOriginal = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.CCAResultsData.r(:,1:DimShow);
                case "Train"
                    plotDataOriginal = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TrainMeanCorrelation;
                case "Test"
                    %plotDataOriginal = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestMeanCorrelation;
                    plotDataOriginal = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestAllCorrelations(:,1:DimShow);
            end
            x = ones(height(plotDataOriginal),1) * (1:DimShow);
            swarmchart(x,plotDataOriginal,options.SwarmSize*2,colorList(1,:),'filled');
            hold on
            RandomConnectomeTestVisualizer.swarmRandomTestCorrelations(globalRunner,DimShow, ...
                "Color",colorList(2,:),"DataType",options.DataType,"ParentAxes",hAx, ...
                "SwarmSize",options.SwarmSize);
            hold on
            RandomConnectomeTestVisualizer.swarmRandomTestCorrelations(localRunner,DimShow, ...
                "Color",colorList(3,:),"DataType",options.DataType,"ParentAxes",hAx, ...
                "SwarmSize",options.SwarmSize);
            ylim([-0.15,1.0001])
            ylabel("Correlation coefficient")
            xlabel("Correlation component (Rank)")
            box off
        end

        function [hAx] = compareAUC(connAnalysisRunner,globalRunner,localRunner,options)
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                globalRunner NullModelAnalysisBatchRunner
                localRunner NullModelAnalysisBatchRunner
                options.DataType = "Full"; % or "Train" or "Test"
                options.SwarmSize = 1.5;
                options.PlotSize = 24;
                options.ParentAxes = [];
                options.ColorList = [];
            end
            % Get axes
            hAx = getOrCreateAxes(options.ParentAxes);
            if isempty(options.ColorList)
                colorList = orderedcolors("gem");
            else
                colorList = options.ColorList;
            end
            switch options.DataType
                case "Full"
                    originalAUC = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.AUC_ROC;
                    globalAUC = globalRunner.FullDataAUCs;
                    localAUC = localRunner.FullDataAUCs;
                case "Train"
                    originalAUC = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TrainMeanAUC_ROC;
                    globalAUC = globalRunner.HoldoutTrainAUCMeans;
                    localAUC = localRunner.HoldoutTrainAUCMeans;
                case "Test"
                    %originalAUC = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestMeanAUC_ROC;
                    originalAUC = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC.';
                    originalAUC_mean = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestMeanAUC_ROC;
                    globalAUC = globalRunner.HoldoutTestAUCMeans;
                    localAUC = localRunner.HoldoutTestAUCMeans;
            end
            x = ones(size(originalAUC));
            swarmchart(x,originalAUC,options.SwarmSize*2,colorList(1,:),'filled');
            hold on 
            x = 2*ones(size(globalAUC));
            swarmchart(x,globalAUC,options.SwarmSize,colorList(2,:),'filled');
            hold on
            x = 3*ones(size(localAUC));
            swarmchart(x,localAUC,options.SwarmSize,colorList(3,:),'filled');
            hold on
            %plot([0.5,3.5],[originalAUC_mean,originalAUC_mean],'k')
            xticks(1:3)
            xticklabels(["Original",strcat("Global"),strcat("Local")])
            ylim([0.45,1])
            yticks([0.5:0.1:1])
            %axis square         
        end

        function h = intervalSwarm(swarmData,intervalWidth,mLocation,options)
            arguments
                swarmData (:,:) double % [NData, DimData]
                intervalWidth (1,1) {mustBeInteger}
                mLocation (1,1) {mustBeInteger}
                options.SwarmSize = 1.5;
                options.DimSHow = 5;
                options.Color = [0 0 0];
            end
            DimShow = options.DimSHow;
            swarmData = swarmData(:,1:DimShow);
            xswarm = ones(height(swarmData),1) * (mLocation:intervalWidth:mLocation+intervalWidth*(DimShow-1));
            h = swarmchart(xswarm,swarmData,options.SwarmSize,options.Color,"filled");
        end

        


    end

end