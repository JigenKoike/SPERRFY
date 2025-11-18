classdef PairDataCCAResultsVisualizer
    % PairDataCCAResultsVisualizer
    % Visualizer for CCA analysis results, including correlation plots
    % and 3D visualization of wiring PI components.

    methods (Static)

        function ax = scatterCCACorrelation(ccaResults, componentIdx, options)
            % Scatter plot of source vs. target projections for a given component.
            arguments
                ccaResults CCAResults
                componentIdx (1,1) {mustBeInteger, mustBePositive} % Index of CCA component to plot
                options.MarkerSize = 1;
                options.Color = [0 0.4470 0.7410];
                options.ParentAxes = [];% must be 'matlab.graphics.axis.Axes' or empty
                options.Axis0Center = true;
                options.SignificantDigits = 2;
                options.PIAxisLabel = "on";
                options.FontSize = 12;
            end
            ax = getOrCreateAxes(options.ParentAxes);
            sourcePIValues = ccaResults.U(:, componentIdx);
            targetPIValues = ccaResults.V(:, componentIdx);
            scatter(ax, sourcePIValues,targetPIValues, ...
                options.MarkerSize,options.Color,'.');
            r = ccaResults.r(componentIdx);
            axis equal
            axis square
            if options.Axis0Center == true
                maxAbsValue = max(abs([sourcePIValues;targetPIValues]));
                xlim([-maxAbsValue,maxAbsValue])
                ylim([-maxAbsValue,maxAbsValue])
            end
            fmt = sprintf('%%.%dg',options.SignificantDigits);
            titleStr = strcat("r_{",num2str(componentIdx),"} = ",compose(fmt,r));
            title(titleStr,"FontSize",options.FontSize);
            box('off');
            if options.PIAxisLabel == "on"
                xlabel(strcat("{{PI}_{s}}^{(",num2str(componentIdx),")}"),"FontSize",options.FontSize)
                ylabel(strcat("{{PI}_{t}}^{(",num2str(componentIdx),")}"),"FontSize",options.FontSize)
            end
        end

        function t = tiledScatterCorrelation(ccaResults,maxComponentId,options)
            arguments
                ccaResults CCAResults
                maxComponentId (1,1) {mustBeInteger, mustBePositive} % Index of CCA component to plot
                options.MarkerSize = [];
                options.Color = [0 0.4470 0.7410];
                options.Axis0Center = true;
            end           
            t = tiledlayout(1,maxComponentId);
            for dPI = 1:maxComponentId
                ax = nexttile;
                PairDataCCAResultsVisualizer.scatterCCACorrelation(ccaResults, dPI, ...
                    "Axis0Center",options.Axis0Center,"Color",options.Color,"MarkerSize",options.MarkerSize, ...
                    "ParentAxes",ax);
            end
        end


        function plotWiringPI3D(wiringPIPairs, stTag, componentIdx, brain3D, brainInfo, options)
            % Plot 3D brain map of wiring PI for a specific component.
            arguments
                wiringPIPairs WiringPIPairs
                stTag string % "s" or "t" ; source/target
                componentIdx (1,1) {mustBeInteger, mustBePositive}
                brain3D BrainSpace3D 
                brainInfo BrainRegionInformation
                options.Clim = [];
                options.ViewAzEl (1,2) {mustBeNumeric} = [70,20];
                options.PlotSize (1,1) {mustBeNumeric} = 8;
                options.ColorList = redblue_cp;
                options.CBarPos = "bottom";
                options.ParentAxes = [];
                options.FontSize = 12;
            end
            hAx = getOrCreateAxes(options.ParentAxes);
            % Extract PI values for plotting
            switch stTag
                case "s"
                    piValues = wiringPIPairs.WiringPISource(:,componentIdx);
                case "t"
                    piValues = wiringPIPairs.WiringPITarget(:,componentIdx);
            end
            if isempty(options.Clim)
                options.Clim = [min(piValues),max(piValues)];
            end
            % Use Brain3DStructureVisualizer for plotting
            Brain3DStructureVisualizer.imageBrainRegionalValuesIn3D(piValues, brain3D, brainInfo, ...
                "Clim",options.Clim,"ColorList",options.ColorList,"CBarPos",options.CBarPos, ...
                "ViewAzEl",options.ViewAzEl,"PlotSize",options.PlotSize, ...
                "ParentAxes",hAx);
            titleStr = strcat("{PI_",stTag,"}^{(",num2str(componentIdx),")}");
            title(hAx, titleStr,"FontSize",options.FontSize);
        end

        function t = tiledPlotWiringPIs3D(wiringPIPairs, dimShow, brain3D, brainInfo, options)
            arguments
                wiringPIPairs WiringPIPairs
                dimShow (1,1) {mustBeInteger, mustBePositive}
                brain3D BrainSpace3D 
                brainInfo BrainRegionInformation
                options.Clim = [];
                options.ViewAzEl (1,2) {mustBeNumeric} = [70,20];
                options.PlotSize (1,1) {mustBeNumeric} = 8;
                options.ColorList = redblue_cp;
                options.ParentAxes = [];
            end
            t = tiledlayout(2,dimShow);
            for dpi = 1:dimShow
                ax = nexttile;
                PairDataCCAResultsVisualizer.plotWiringPI3D(wiringPIPairs,"s",dpi,brain3D,brainInfo, ...
                    "Clim",options.Clim,"ViewAzEl",options.ViewAzEl,"PlotSize",options.PlotSize, ...
                    "ColorList",options.ColorList,"CBarPos","none","ParentAxes",ax);
            end
            for dpi = 1:dimShow
                ax = nexttile;
                PairDataCCAResultsVisualizer.plotWiringPI3D(wiringPIPairs,"t",dpi,brain3D,brainInfo, ...
                    "Clim",options.Clim,"ViewAzEl",options.ViewAzEl,"PlotSize",options.PlotSize, ...
                    "ColorList",options.ColorList,"CBarPos","bottom","ParentAxes",ax);             
            end
            
        end




        function [t] = tiledSwarmPIValuesByMR(wiringPIPairs,crossInfo,DimShow,options)
            arguments
                wiringPIPairs WiringPIPairs
                crossInfo CrossRegionInformation
                DimShow (1,1) {mustBeInteger}
                options.ParentAxes = [];
                options.SwarmSize = 16;
                options.SrcLabelOff logical = 1;
            end
            hAx = getOrCreateAxes(options.ParentAxes);
            load("colorList13.mat")
            sz = options.SwarmSize;
            [srcMRNames,~,~] = getMajorRegionInfo(crossInfo.SourceRegionInfo);
            [tgtMRNames,~,~] = getMajorRegionInfo(crossInfo.TargetRegionInfo);
            srcPIValuesPerRegion = wiringPIPairs.WiringPISource;
            tgtPIValuesPerRegion = wiringPIPairs.WiringPITarget;
            majorRegionNameCategoricalListSource = getMajorRegions(crossInfo.SourceRegionInfo);
            majorRegionNameCategoricalListTarget = getMajorRegions(crossInfo.TargetRegionInfo);
            xSwarmSource = categorical(majorRegionNameCategoricalListSource,srcMRNames,'Ordinal',true);
            xSwarmTarget = categorical(majorRegionNameCategoricalListTarget,tgtMRNames,'Ordinal',true);
            cs = makeMRColorListPerRegion(crossInfo.SourceRegionInfo);
            ct = makeMRColorListPerRegion(crossInfo.TargetRegionInfo);
            ylimList = zeros(DimShow,2);
            for ds = 1:DimShow
                srcPI = srcPIValuesPerRegion(:,ds);
                tgtPI = tgtPIValuesPerRegion(:,ds);
                ylimList(ds,1) = min([srcPI;tgtPI]);
                ylimList(ds,2) = max([srcPI;tgtPI]);
            end
            t = tiledlayout(2,DimShow);
            for ds = 1:DimShow
                nexttile
                ySwarmSource = srcPIValuesPerRegion(:,ds);
                swarmchart(xSwarmSource,ySwarmSource,sz,cs,'filled')
                title(strcat("{PI_s}^{(",num2str(ds),")}"))
                ylim(ylimList(ds,:));
                if options.SrcLabelOff == 1;
                    xticklabels([]);
                end
            end
            for ds = 1:DimShow
                nexttile
                ySwarmTarget = tgtPIValuesPerRegion(:,ds);
                swarmchart(xSwarmTarget,ySwarmTarget,sz,ct,'filled')
                title(strcat("{PI_t}^{(",num2str(ds),")}"))
                ylim(ylimList(ds,:));
            end
                
        end

    end
end


