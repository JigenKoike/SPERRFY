classdef FigMaker

    properties
        ColorListGem (:,3) double
        ColorList13 (:,3) double
    end

    methods
        function obj = FigMaker()
            obj.ColorListGem = orderedcolors("gem");
            load("colorList13.mat") 
            obj.ColorList13 = colorList13;
        end
    end


    methods (Static)
        function posiVec = figPosition(fwidth,fdepth)
            xleft = 55;
            ybottom = 5;
            posiVec = [xleft ybottom xleft+fwidth ybottom+fdepth];
        end

        % Figure 3
        function f = imageConnectionMatrix_MRColorLabels(connMat,crossInfo,options)
            % for f3a, fs3(1:4), f6a(1:2)
            arguments(Input)
                connMat ConnectionMatrix
                crossInfo CrossRegionInformation
                options.DataType = "Connections";
                options.FWidth = 390;
                options.FDepth = 440;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            ConnectomeVisualizer.plotConnectomeMatrix_withMRColorLabels(connMat,crossInfo,"DataType",options.DataType);
            set(f,'Color','w')
        end

        function f = plotAllCorrCoeff_FullAndHoldoutTest(connAnalysisRunner,options)
            % f3b
            arguments(Input)
                connAnalysisRunner ConnectomeAnalysisRunner
                options.FWidth = 1000;
                options.FDepth = 440;
                options.DimShow = 50;
                options.SzScatter = 6;
                options.SzFullDataRatio = 3;

            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            colorList = orderedcolors("gem");
            c = colorList(3,:);
            dimShow = options.DimShow;
            sz = options.SzScatter;
            % holdout test
            testCorr = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestAllCorrelations;
            x = ones(height(testCorr),1) * (1:dimShow);
            h1 = swarmchart(x,testCorr(:,1:dimShow),sz,c*1.05,'filled','AlphaData',0.5);
            xlabel("Correlation Componet (rank)")
            ylabel("Correlation Coefficients")
            xlim([0,dimShow+1])
            yLimVec = [min(testCorr,[],'all')-0.005,1];
            ylim(yLimVec)
            %title(strcat("Correlation Coefficient"))
            hold on
            % full data
            x = 1:dimShow;
            y = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.CCAResultsData.r(1,:);
            h2 = plot(x,y,'LineStyle','-','Color',colorList(2,:),'Marker','.','LineWidth',1,'MarkerSize',sz*options.SzFullDataRatio);
            %h2 = scatter(x,y,sz*options.SzFullDataRatio,'red','filled');
            hold on
            y = mean(testCorr,1);
            h3 = plot(x,y,'LineStyle','-','Color',c*0.9,'Marker','.','LineWidth',1,'MarkerSize',sz*options.SzFullDataRatio);
            %h2 = scatter(x,y,sz*options.SzFullDataRatio,'red','filled');
            hold on
            x = [5.5,5.5];
            y = yLimVec;
            h4 = plot(x,y,'Color','k','LineStyle','--','LineWidth',0.5);
            legend([h2,h3],["Full Data","Holdout Test"],'FontSize',12);
            set(f,'Color','w')
        end

        function f = scatterCorrelations(connAnalysisRunner,options)
            % f3c
            arguments(Input)
                connAnalysisRunner ConnectomeAnalysisRunner
                options.DimShow = 5;
                options.FWidthUnit = 200;
                options.FDepth = 200;
                options.Color = [0.05 0.4670 0.7610];%[0 0.4470 0.7410];
            end
            DimShow = options.DimShow;
            mkSz = 1;
            fwidth = options.FWidthUnit*DimShow;
            fdepth = options.FDepth;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            ccaResults = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.CCAResultsData;
            clist = orderedcolors("gem");
            t = PairDataCCAResultsVisualizer.tiledScatterCorrelation(ccaResults,DimShow,"MarkerSize",mkSz,"Color",options.Color);
            set(f,'Color','w')
        end

        function f = imageWiringPIs(connAnalysisRunner,brainSpace,brainInfo,options)
            % f3d
            arguments(Input)
                connAnalysisRunner ConnectomeAnalysisRunner
                brainSpace BrainSpace3D
                brainInfo BrainRegionInformation
                options.DimShow = 5;
                options.FWidthUnit = 250;
                options.FDepthUnit = 225;
            end
            climVec = [-2.5,2.5];
            DimShow = options.DimShow;
            fwidth = options.FWidthUnit*DimShow;
            fdepth = options.FDepthUnit*2;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            wiringPIPairs = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
            PairDataCCAResultsVisualizer.tiledPlotWiringPIs3D(wiringPIPairs,DimShow,brainSpace,brainInfo, ...
                "Clim",climVec);
            set(f,'Color','w')
        end

        % Figure S1
        function f = image3DMajorRegions(brainSpace,brainInfo,options)
            % fs1a
            arguments(Input)
                brainSpace BrainSpace3D
                brainInfo BrainRegionInformation
                options.FWidth = 500;
                options.FDepth = 400;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            Brain3DStructureVisualizer.image3DBrainMajorRegion(brainSpace,brainInfo);
            set(f,'Color','w')
        end

        
        function f = imageConnectomeProperties(connMat,crossInfo,options)
            % fs1b
            arguments(Input)
                connMat ConnectionMatrix
                crossInfo CrossRegionInformation
                options.FWidth = 950;
                options.FDepth = 450;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            t = tiledlayout(1,2);
            ax = nexttile;
            ConnectomeVisualizer.heatmapProjectionFeaturesBetweenMRs(connMat,crossInfo,"Feature","counts","ParentAxes",ax);
            colorbar off
            ax = nexttile;
            ConnectomeVisualizer.heatmapProjectionFeaturesBetweenMRs(connMat,crossInfo,"Feature","density","ParentAxes",ax);
            colorbar
            h.FontName = 'Arial';
            set(f,'Color','w')
        end

        % Figure S2
        function f = imageGeneExpressionMatrix(geneExprLevels,brainInfo,options)
            % fs2a
            arguments(Input)
                geneExprLevels GeneExpressionLevels
                brainInfo BrainRegionInformation
                options.FWidth = 415;
                options.FDepth = 400;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            geneExpressionMatrix = geneExprLevels.ClusteringResults.SortedGeneExpression;
            mainAx = axes('Position',[0.1 0.1 0.8 0.8],'OuterPosition',[0 0 1 1]);
            hIm = imagesc(geneExpressionMatrix);
            mainAx.PlotBoxAspectRatio = [1 1 1];
            % colorbar setting
            cb = colorbar;
            expStd = mean(std(geneExpressionMatrix,1,1));
            cLimVec = [-2*expStd,2*expStd];
            clim(cLimVec);
            cb.Ticks = cb.Limits;
            cb.TickLabels = {'Low','High'};
            % label setting
            [~,MRInitialIndexList,~] = getMajorRegionInfo(brainInfo);
            yticks(MRInitialIndexList - 0.5);
            yticklabels([])
            NGenes = width(geneExpressionMatrix);
            xlabelString = strcat(num2str(NGenes)," genes");
            xlabel(xlabelString)
            xticks([])
            % MR color setting
            gapWidth = 0.01;
            barWidth = 0.01;
            hold on
            leftAx = setMRColorLabel(brainInfo);
            leftAx.Position = [mainAx.Position(1) - gapWidth, ...
                                   mainAx.Position(2), ...
                                   barWidth, ...
                                   mainAx.Position(4)];
            linkaxes([mainAx,leftAx],'y');
            set(f,'Color','w')
        end

        function f = imagePCMatrix(geneExprLevels,brainInfo,DimPC,options)
            % fs2b
            arguments(Input)
                geneExprLevels GeneExpressionLevels
                brainInfo BrainRegionInformation
                DimPC {mustBeInteger} = 50;
                options.FWidth = 415;
                options.FDepth = 400;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            expressionMatrix = geneExprLevels.PCAMatrix;
            expressionMatrix = expressionMatrix(:,1:DimPC);
            mainAx = axes('Position',[0.1 0.1 0.8 0.8],'OuterPosition',[0 0 1 1]);
            im = imagesc(expressionMatrix);
            mainAx.PlotBoxAspectRatio = [1 1 1];
            % colorbar setting
            cb = colorbar;
            expStd = mean(std(expressionMatrix,1,1));
            cLimVec = [-2*expStd,2*expStd];
            clim(cLimVec);
            cb.Ticks = cb.Limits;
            cb.TickLabels = {'Low','High'};
            % label setting
            [~,MRInitialIndexList,~] = getMajorRegionInfo(brainInfo);
            yticks(MRInitialIndexList - 0.5);
            yticklabels([])
            NGenes = width(expressionMatrix);
            xlabelString = strcat("#PC");
            xlabel(xlabelString)
            xticks([])
            % MR color setting
            gapWidth = 0.01;
            barWidth = 0.01;
            hold on
            leftAx = setMRColorLabel(brainInfo);
            leftAx.Position = [mainAx.Position(1) - gapWidth, ...
                                   mainAx.Position(2), ...
                                   barWidth, ...
                                   mainAx.Position(4)];
            linkaxes([mainAx,leftAx],'y');
            set(f,'Color','w')
        end

        function f = plotPCACumContribution(geneExprLevels,options)
            % fs2c
            arguments
                geneExprLevels GeneExpressionLevels
                options.DimPCA = 50;
                options.FWidth = 400;
                options.FDepth = 300;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            geneExprLevels.PCAInfo.plotCumulativeContribution("PlotDim",options.DimPCA);
            set(f,'Color','w')          
        end

        % Figure S5
        function f = swarmPIValuesByMR(connAnalysisRunner,options)
            % fs5a
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                options.DimShow = 5;
                options.FWidthUnit = 300;
                options.FDepthUnit = 250;
            end
            DimShow = options.DimShow;
            fwidth = options.FWidthUnit*DimShow;
            fdepth = options.FDepthUnit*2;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            wiringPIPairs = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.WiringPIPairsData;
            crossInfo = connAnalysisRunner.Factory.CrossRegionInformation;
            PairDataCCAResultsVisualizer.tiledSwarmPIValuesByMR(wiringPIPairs,crossInfo,DimShow);
            set(f,'Color','w')
        end

        function f = scatterCorrelations_SpecificMRPairColor(connAnalysisRunner,pairDataGenerationOptions,options)
            % fs5b
            arguments(Input)
                connAnalysisRunner ConnectomeAnalysisRunner
                pairDataGenerationOptions PairDataGenerationOptions
                options.DimShow = 5;
                options.MkSize = 12;
                options.FWidthUnit = 400;
                options.FDepth = 450;
                options.PIAxisLabel = "on";
                options.IsLegend = 0;
                options.IsAxis0Center = 0;
            end
            DimShow = options.DimShow;
            mkSz = options.MkSize;
            fwidth = options.FWidthUnit*DimShow;
            fdepth = options.FDepth;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            U = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.CCAResultsData.U;
            V = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.CCAResultsData.V;
            % make focusing MR pair list
            focusingMRPairList = string({
                "Isocortex","Thalamus"; ...
                "Isocortex","Midbrain"; ...
                "Midbrain","Thalamus"; ...
                "Thalamus","Isocortex"; ...
                "Hypothalamus","Thalamus"});
            speceficMRPairIndexList = zeros(height(U),1);
            for k = 1:height(focusingMRPairList)
                [spcfcMRPairMask,~,~] = connAnalysisRunner.Factory.getSpecificMRPairMask( ...
                    focusingMRPairList(k,1),focusingMRPairList(k,2),pairDataGenerationOptions);
                speceficMRPairIndexList = speceficMRPairIndexList + spcfcMRPairMask * k;
            end
            % make legend string
            colorListTmp = [0.5,0.5,0.5;orderedcolors("reef")];
            legendStr = string({
                "Isocortex → Thalamus"; ...
                "Isocortex → Midbrain"; ...
                "Midbrain → Thalamus"; ...
                "Thalamus → Isocortex"; ...
                "Hypothalamus → Thalamus"; ...
                "Others"});
            % make plot
            t = tiledlayout(1,5);
            for d = 1:5
                nexttile
                limVec = [min(min(U(:,d)),min(V(:,d))),max(max(U(:,d)),max(V(:,d)))];
                if options.IsAxis0Center == 1
                    limVec = [-max(abs(limVec)),max(abs(limVec))];
                end
                for k = 0:5
                    u = U(speceficMRPairIndexList==k,d);
                    v = V(speceficMRPairIndexList==k,d);
                    if k == 0
                        h(k+1) = scatter(u,v,mkSz*0.8,colorListTmp(k+1,:),'filled',"MarkerFaceAlpha",1);
                    else
                        h(k+1) = scatter(u,v,mkSz,colorListTmp(k+1,:),'filled');
                    end
                    hold on
                end
                if options.PIAxisLabel == "on"
                xlabel(strcat("{{PI}_{s}}^{(",num2str(d),")}"))
                ylabel(strcat("{{PI}_{t}}^{(",num2str(d),")}"))
                end
                box off
                axis square
                xlim(limVec)
                ylim(limVec)
                if d == 1
                    if options.IsLegend == 1;
                        legend([h(2),h(3),h(4),h(5),h(6),h(1)],legendStr,'Location','northwest');
                    end
                end
            end
            set(f,'Color','w')
        end

        % Figure S4
        function f = plotWiringPISpatialAutocorrelation(saVec_PIs,saVec_PIt,options)
            % fs4a
            arguments
                saVec_PIs (:,1) double
                saVec_PIt (:,1) double
                options.FWidth = 1000;
                options.FDepth = 400;
                options.LineWidth = 2;
                options.FontSize = 14;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            dimPI = numel(saVec_PIs);
            clist = orderedcolors("gem");
            plot([1:dimPI],saVec_PIs.','Color',clist(1,:),'LineStyle','-','Marker','.','MarkerSize',16,'LineWidth',options.LineWidth);
            hold on
            plot([1:dimPI],saVec_PIt.','Color',clist(2,:),'LineStyle','-','Marker','.','MarkerSize',16,'LineWidth',options.LineWidth);
            ylim([0,1])
            legend(["Source","Target"],'FontSize',options.FontSize)
            %title("Spatial autocorrelation of wiring PI")
            ylabel("Moran's I",'FontSize',options.FontSize)
            xlabel("Correlation component (Rank)",'FontSize',options.FontSize)
            box off
            set(f,'Color','w')
        end

        function f = heatmapWiringPIGradientCorrelation(PISim_ss,PISim_tt,PISim_st,options)
            % fs4a
            arguments
                PISim_ss (:,:) double
                PISim_tt (:,:) double
                PISim_st (:,:) double
                options.FWidth = 975;
                options.FDepth = 300;
                options.FontSize = 14;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            t = tiledlayout(1,3);
            %title(t,"Correlation between each wiring PI Gradient")
            nexttile
            heatmap(PISim_ss,"CellLabelFormat",'%.2f');
            colormap(redblue_cp)
            clim([-1,1])
            title("Source - Source")
            xlabel("Source")
            ylabel("Source")
            colorbar off
            nexttile
            heatmap(PISim_tt,"CellLabelFormat",'%.2f');
            colormap(redblue_cp)
            clim([-1,1])
            title("Target - Target")
            xlabel("Target")
            ylabel("Target")
            colorbar off
            nexttile
            heatmap(PISim_st,"CellLabelFormat",'%.2f');
            colormap(redblue_cp)
            clim([-1,1])
            title("Source - Target")
            xlabel("Target")
            ylabel("Source")
            set(f,'Color','w')
        end

        % Figure 4
        function f = histogramsGeneDistributionSimilarity(relatedGeneAnalysisResults,options)
            % f4a
            arguments
                relatedGeneAnalysisResults RelatedGeneAnalysisResults
                options.IsAbs = 0;
                options.DimShow = 5;
                options.FWidthUnit = 225;
                options.FDepthUnit = 225;
                options.Color = [0 0.4270 0.6610];
            end
            DimShow = options.DimShow;
            fwidth = options.FWidthUnit*DimShow;
            fdepth = options.FDepthUnit*2;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            t = relatedGeneAnalysisResults.tiledHistogramsGeneCorrCoef(DimShow,"IsAbs",options.IsAbs,"Color",options.Color);
            set(f,'Color','w')
        end

        function f = image3DSimilarGeneDistributions(relatedGeneAnalysisResults,brainSpace,brainInfo,options)
            arguments
                relatedGeneAnalysisResults RelatedGeneAnalysisResults
                brainSpace BrainSpace3D
                brainInfo BrainRegionInformation
                options.DimShow = 5;
                options.FWidthUnit = 250;
                options.FDepthUnit = 225;
            end
            DimShow = options.DimShow;
            fwidth = options.FWidthUnit*DimShow;
            fdepth = options.FDepthUnit*2;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            t = relatedGeneAnalysisResults.tiledImage3DMostCorrelatedGenes(brainSpace,brainInfo,DimShow);
            set(f,'Color','w')
        end

        % Figure S7
        function f = scatterGenePISimSpatialAutocorrelation(relatedGeneAnalysisResults,geneSAVec,options)
            % fs7
            arguments
                relatedGeneAnalysisResults RelatedGeneAnalysisResults
                geneSAVec (:,1) double
                %options.IsAbs = 0;
                options.DimShow = 5;
                options.FWidthUnit = 225;
                options.FDepthUnit = 225;
            end
            DimShow = options.DimShow;
            fwidth = options.FWidthUnit*DimShow;
            fdepth = options.FDepthUnit*2;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));            
            t = tiledlayout(2,DimShow);
            %title(t,"PI similarity vs Spatial autocorrelation")
            x = geneSAVec;
            for d = 1:DimShow
                nexttile
                y = abs(relatedGeneAnalysisResults.CorrCoefMatrixSource(:,d));
                scatter(x,y,'.');
                r = corr(x,y);
                titleStr = strcat("{PI_s}^{(",num2str(d),")} : r = ",num2str(round(r,3)));
                title(titleStr);
                xlim([0,1])
                ylim([0,1])
                axis square
            end
            for d = 1:DimShow
                nexttile
                y = abs(relatedGeneAnalysisResults.CorrCoefMatrixTarget(:,d));
                scatter(x,y,'.');
                r = corr(x,y);
                hold on
                titleStr = strcat("{PI_t}^{(",num2str(d),")} : r = ",num2str(round(r,3)));
                title(titleStr);
                xlim([0,1])
                ylim([0,1])
                axis square
            end
            xlabel(t,"Spatial autocorrelation",'FontName','Arial')
            ylabel(t,"Cosine similarity (Absolute value)",'FontName','Arial') 
            set(f,'Color','w')
        end

        % Figure 5
        function f = imageWiringPIDiffForReconstruction(connAnalysisRunner,options)
            % f5a
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                options.FWidth = 390;
                options.FDepth = 440;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            piDiffMat = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.PIDistanceMatrix;
            connMat = connAnalysisRunner.Factory.ConnectionMatrix;
            crossInfo = connAnalysisRunner.Factory.CrossRegionInformation;
            ConnectomeVisualizer.plotConnectomeMatrix_withMRColorLabels(connMat,crossInfo,"AnyMatrix",piDiffMat);
            set(f,'Color','w')
        end

        function f = plotROCHoldoutAndDistReconst(connAnalysisRunner,distanceReconstResults,options)
            % f5b
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                distanceReconstResults ReconstructionResults
                options.FWidth = 300;
                options.FDepth = 300;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));            
            colorlist = orderedcolors("gem");  
            c_holdout = [colorlist(3,:)*1.05 0.075];
            for i = 1:connAnalysisRunner.HoldoutOptions.NumSplits
                h = gca();
                holdoutReconstResults = connAnalysisRunner.OverallModelAndResults.HoldoutUnits(i).TestMetrics.TestPseudeReconstructionResults;
                h1(i) = holdoutReconstResults.plotROC("ParentAxes",h,"TitleOff",true,"LineWidth",0.1,"PlotColor",c_holdout,"ReturnAx",1);
                hold on
            end
            [meanFPR,meanTPR] = connAnalysisRunner.OverallModelAndResults.getHoldoutROCmean("TestOrTrain","Test");
            meanAUC = mean(connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC,2);
            h2 = plot(meanFPR,meanTPR,"Color",colorlist(3,:)*0.9,"LineStyle","-","LineWidth",2);
            %titleStr = strcat("ROC plot of connection reconstruction");
            %title(gca,titleStr)
            hold on
            c_distance = [0 0 0];
            h3 = distanceReconstResults.plotROC("ParentAxes",h,"TitleOff",true,"LineWidth",1,"PlotColor",c_distance,"ReturnAx",1);            
            h4 = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.plotROC(...
                "ParentAxes",h,"TitleOff",true,"LineWidth",1,"PlotColor",colorlist(2,:),"ReturnAx",1);
            %distAUC = distanceReconstResults.AUC_ROC;
            legendStr2 = strcat("Wiring PI (Holdout)");
            legendStr3 = strcat("Physical region distance");
            legendStr4 = strcat("Wiring PI (Full data)");
            legendStr = [legendStr2,legendStr4,legendStr3];
            legend([h2 h4(1) h3(1)],legendStr,"Location","southeast","FontSize",12);
            set(f,'Color','w')
        end

        function f = swarmHoldoutAUC(connAnalysisRunner,distanceReconstResults,options)
            % f5c
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                distanceReconstResults ReconstructionResults
                options.FWidth = 200;
                options.FDepth = 300;
                options.SzSwarm = 6;
                options.SzDataRation = 6;
            end
            clist = orderedcolors("gem");
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            fullAUC = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.AUC_ROC;
            testAUCs = connAnalysisRunner.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC;
            meanTestAUC = mean(testAUCs);
            distanceAUC = distanceReconstResults.AUC_ROC;
            xswarm = zeros(size(testAUCs));
            h1 = swarmchart(xswarm,testAUCs,options.SzSwarm,clist(3,:)*1.05,"filled");
            hold on
            h2 = scatter(0,meanTestAUC,options.SzSwarm*options.SzDataRation,clist(3,:)*0.9,'filled');
            hold on 
            h3 = scatter(0,fullAUC,options.SzSwarm*options.SzDataRation,clist(2,:),'filled');
            hold on
            h4 = scatter(0,distanceAUC,options.SzSwarm*options.SzDataRation,[0.2,0.2,0.2],'filled');
            xlim([-0.6,0.6])
            ylim([0.45,1])
            xticks([])
            ylabel("AUC")
            legendStr2 = strcat("Wiring PI (Holdout)");
            legendStr3 = strcat("Wiring PI (Full data)");
            legendStr4 = strcat("Physical region distance");
            legendStr = [legendStr2,legendStr3,legendStr4];
            legend([h2 h3 h4],legendStr,"Location","south","FontSize",12)

            set(f,'Color','w')
        end

        function f = scatterPIDiffvsRegionDistance(connAnalysisRunner,options)
            % f5d
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                options.SzScatter = 3;
                options.FWidth = 325;
                options.FDepth = 325;
                options.AlphaBlack = 0.6;
                options.AlphaRed = 0.8;
                options.FontSize = 14;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            sz = options.SzScatter;
            PIDiffMat = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.PIDistanceMatrix;
            crossInfo = connAnalysisRunner.Factory.CrossRegionInformation;
            connectomeMatrix = connAnalysisRunner.Factory.ConnectionMatrix;
            distanceMatrix = crossInfo.DistanceMatrix;
            x_all = distanceMatrix(crossInfo.DomainDefineMatrix == 1);
            y_all = PIDiffMat(crossInfo.DomainDefineMatrix == 1);
            r_all = corr(x_all,y_all);
            x_conn = distanceMatrix(connectomeMatrix.Matrix .* crossInfo.DomainDefineMatrix == 1);
            y_conn = PIDiffMat(connectomeMatrix.Matrix .* crossInfo.DomainDefineMatrix == 1);
            r_conn = corr(x_conn,y_conn);
            x_noConn = distanceMatrix((~connectomeMatrix.Matrix) .* crossInfo.DomainDefineMatrix == 1);
            y_noConn = PIDiffMat((~connectomeMatrix.Matrix) .* crossInfo.DomainDefineMatrix == 1);
            r_noConn = corr(x_noConn,y_noConn);
            nexttile
            h1 = scatter(x_noConn,y_noConn,sz,'k','filled','MarkerFaceAlpha',options.AlphaBlack);
            hold on
            h2 = scatter(x_conn,y_conn,sz,'red','filled','MarkerFaceAlpha',options.AlphaRed);
            legend([h2, h1],["Connected","Unconnected"],"Location","northeast","FontSize",options.FontSize*0.8);
            titleStr = strcat("All region pair: r = ",num2str(round(r_all,2)));
            title(titleStr,"FontSize",options.FontSize)
            xlim([0,12500])
            ylim([0,1])
            xlabel("Region distance (μm)","FontSize",options.FontSize)
            axis square                
            ylabel("Wiring PI difference (d_{PI})","FontSize",options.FontSize)
            set(f,'Color','w')
        end


        function f = scatterPIDiffvsRegionDistance_tiled(connAnalysisRunner,options)
            % f5d_old
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                options.SzScatter = 3;
                options.FWidth = 1000;
                options.FDepth = 325;
                options.AlphaBlack = 0.6;
                options.AlphaRed = 0.8;
                options.FontSize = 14;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            sz = options.SzScatter;
            PIDiffMat = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData.PIDistanceMatrix;
            crossInfo = connAnalysisRunner.Factory.CrossRegionInformation;
            connectomeMatrix = connAnalysisRunner.Factory.ConnectionMatrix;
            distanceMatrix = crossInfo.DistanceMatrix;
            t = tiledlayout(1,3);
            x_all = distanceMatrix(crossInfo.DomainDefineMatrix == 1);
            y_all = PIDiffMat(crossInfo.DomainDefineMatrix == 1);
            r_all = corr(x_all,y_all);
            x_conn = distanceMatrix(connectomeMatrix.Matrix .* crossInfo.DomainDefineMatrix == 1);
            y_conn = PIDiffMat(connectomeMatrix.Matrix .* crossInfo.DomainDefineMatrix == 1);
            r_conn = corr(x_conn,y_conn);
            x_noConn = distanceMatrix((~connectomeMatrix.Matrix) .* crossInfo.DomainDefineMatrix == 1);
            y_noConn = PIDiffMat((~connectomeMatrix.Matrix) .* crossInfo.DomainDefineMatrix == 1);
            r_noConn = corr(x_noConn,y_noConn);
            nexttile
            h1 = scatter(x_noConn,y_noConn,sz,'k','filled','MarkerFaceAlpha',options.AlphaBlack);
            hold on
            h2 = scatter(x_conn,y_conn,sz,'red','filled','MarkerFaceAlpha',options.AlphaRed);
            legend([h2, h1],["Connected","Unconnected"],"Location","northeast","FontSize",options.FontSize*0.75);
            titleStr = strcat("All region pair: r = ",num2str(round(r_all,2)));
            title(titleStr,"FontSize",options.FontSize)
            xlim([0,12500])
            ylim([0,1])
            xlabel("Region distance","FontSize",options.FontSize)
            axis square
            nexttile
            scatter(x_conn,y_conn,sz,'red','filled','MarkerFaceAlpha',options.AlphaBlack);
            titleStr = strcat("Connected: r = ",num2str(round(r_conn,2)));
            title(titleStr,"FontSize",options.FontSize)
            xlim([0,12500])
            ylim([0,1])
            xlabel("Region distance","FontSize",options.FontSize)
            axis square
            nexttile
            scatter(x_noConn,y_noConn,sz,'k','filled','MarkerFaceAlpha',options.AlphaRed);
            titleStr = strcat("Unconnected: r = ",num2str(round(r_noConn,2)));
            title(titleStr,"FontSize",options.FontSize)
            xlim([0,12500])
            ylim([0,1])
            xlabel("Region Distance","FontSize",options.FontSize)
            axis square            
            ylabel(t,"Wiring PI Difference","FontSize",options.FontSize)
            set(f,'Color','w')
        end

        % Figure S8
        function f = imageReconstructedMatrixSamples(connAnalysisRunner,options)
            % f8sa
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                options.ThreshList = [0.14,0.18,0.22];
                options.FWidth = 390;
                options.FDepth = 440;
            end      
            reconstResults = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData;
            threshIndices = reconstResults.getThresholdIndices(options.ThreshList);
            crossInfo = connAnalysisRunner.Factory.CrossRegionInformation;
            for i = 1:numel(options.ThreshList)
                f(i) = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
                predMat = cell2mat(reconstResults.PredictedMatrices(threshIndices(i)));
                reconstMat = ConnectionMatrix(predMat);
                ConnectomeVisualizer.plotConnectomeMatrix_withMRColorLabels(reconstMat,crossInfo,"DataType","Connections");
                set(f,'Color','w')
            end
        end

        function f = plotReconstructionROC(connAnalysisRunner,options)
            % f8sb
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                options.ThreshList = [0.14,0.18,0.22];
                options.FWidth = 300;
                options.FDepth = 300;
            end      
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            reconstResults = connAnalysisRunner.OverallModelAndResults.FullAnalysisModel.ReconstructionResultsData;
            reconstResults.plotROC("PlotThresholds",options.ThreshList);
            set(f,'Color','w')
        end

        % Figure 6
        function f = plotRandomTestCorrCoef(connAnalysisRunner,globalRunner,localRunner,options)
            % f6b
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                globalRunner NullModelAnalysisBatchRunner
                localRunner NullModelAnalysisBatchRunner
                options.DataType = "Test"; % or "Train" or "Test"
                options.StdErrorBar = "off";
                options.MarkerSize = 12;
                options.LineWidth = 1.5;
                options.FWidth = 450;
                options.FDepth = 300;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            RandomConnectomeTestVisualizer.compareCorrPlot(connAnalysisRunner,globalRunner,localRunner,...
                "DataType",options.DataType,"StdErrorBar",options.StdErrorBar,"MarkerSize",options.MarkerSize);
            set(f,'Color','w')
        end

        function f = swarmRandomTestCorrCoef(connAnalysisRunner,globalRunner,localRunner,options)
            % f6d
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                globalRunner NullModelAnalysisBatchRunner
                localRunner NullModelAnalysisBatchRunner
                options.DataType = "Test"; % or "Train" or "Test"
                options.DimShow = 5;
                options.SwarmSize = 1.5;
                options.FWidth = 300;
                options.FDepth = 300;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            RandomConnectomeTestVisualizer.compareCorrSwarm(connAnalysisRunner,globalRunner,localRunner,options.DimShow,...
                "SwarmSize",options.SwarmSize,"DataType","Test");
            box off
            set(f,'Color','w')
        end

        function f = swarmRandomTestAUC(connAnalysisRunner,globalRunner,localRunner,options)
            % f6d
            arguments
                connAnalysisRunner ConnectomeAnalysisRunner
                globalRunner NullModelAnalysisBatchRunner
                localRunner NullModelAnalysisBatchRunner
                options.DataType = "Test"; % or "Train" or "Test"
                options.SwarmSize = 1.5;
                options.FWidth = 150;
                options.FDepth = 300;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            RandomConnectomeTestVisualizer.compareAUC(connAnalysisRunner,globalRunner,localRunner,...
                "SwarmSize",options.SwarmSize,"DataType",options.DataType);
            ylabel("AUC")
            box off
            set(f,'Color','w')
        end

        function f = swarmRandomTestSimilarity(globalRunner,localRunner,options)
            % f6d
            arguments
                globalRunner NullModelAnalysisBatchRunner
                localRunner NullModelAnalysisBatchRunner
                options.SwarmSize = 1.5;
                options.FWidth = 450;
                options.FDepth = 300;
                options.DimShow = 5;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            globalSimilarities = globalRunner.WirigPISimilarities(:,1:options.DimShow);
            localSimilarities = localRunner.WirigPISimilarities(:,1:options.DimShow);
            colorlist = orderedcolors("gem");
            h1 = RandomConnectomeTestVisualizer.intervalSwarm(globalSimilarities,3,1,"Color",colorlist(2,:),"SwarmSize",options.SwarmSize);
            hold on
            h2 = RandomConnectomeTestVisualizer.intervalSwarm(localSimilarities,3,2,"Color",colorlist(3,:),"SwarmSize",options.SwarmSize);
            hold on
            ylabel("Similarity score")
            yticks([0:0.2:1])
            xticks([1.5:3:1.5+3*(options.DimShow-1)])
            xticklabels([1:5])
            xlabel("Correlation components (Rank)")
            box off
            set(f,'Color','w')
        end

        % Figure S9
        function f = swarmNullModelsCorrCoeff_3nulls(origRunner,nullRunner1,nullRunner2,nullRunner3,options)
            % fs9a(1)
            arguments
                origRunner ConnectomeAnalysisRunner
                nullRunner1 NullModelAnalysisBatchRunner
                nullRunner2 NullModelAnalysisBatchRunner
                nullRunner3 NullModelAnalysisBatchRunner
                options.SwarmSize = 1.5;
                options.FWidth = 450;
                options.FDepth = 300;
                options.DimShow = 5;
                options.ColorList = orderedcolors("gem");
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            origCorrCoeff = origRunner.OverallModelAndResults.HoldoutSummary.TestAllCorrelations(:,1:options.DimShow);
            null1CorrCoeff = nullRunner1.HoldoutTestCorrMeans(:,1:options.DimShow);
            null2CorrCoeff = nullRunner2.HoldoutTestCorrMeans(:,1:options.DimShow);
            null3CorrCoeff = nullRunner3.HoldoutTestCorrMeans(:,1:options.DimShow);
            intervalWidth = 5;
            colorlist = options.ColorList;
            h1 = RandomConnectomeTestVisualizer.intervalSwarm(origCorrCoeff,intervalWidth,1,"Color",colorlist(1,:),"SwarmSize",options.SwarmSize);
            hold on
            h2 = RandomConnectomeTestVisualizer.intervalSwarm(null1CorrCoeff,intervalWidth,2,"Color",colorlist(2,:));
            hold on
            h3 = RandomConnectomeTestVisualizer.intervalSwarm(null2CorrCoeff,intervalWidth,3,"Color",colorlist(3,:));
            hold on
            h4 = RandomConnectomeTestVisualizer.intervalSwarm(null3CorrCoeff,intervalWidth,4,"Color",colorlist(4,:));
            hold on
            ylabel("Correlation coefficient")
            yticks([0:0.2:1])
            ylim([-0.1,1])
            tickStart = (intervalWidth)/2;
            xticks([tickStart:intervalWidth:tickStart+intervalWidth*(options.DimShow-1)])
            xticklabels([1:options.DimShow])
            xlabel("Correlation Components (Rank)")
            box off
            set(f,'Color','w')
        end

        function f = swarmNullModelsAUC_3nulls(origRunner,nullRunner1,nullRunner2,nullRunner3,options)
            % fs9a(2)
            arguments
                origRunner ConnectomeAnalysisRunner
                nullRunner1 NullModelAnalysisBatchRunner
                nullRunner2 NullModelAnalysisBatchRunner
                nullRunner3 NullModelAnalysisBatchRunner
                options.SwarmSize = 1.5;
                options.FWidth = 200;
                options.FDepth = 300;
                options.DimShow = 1;
                options.ColorList = orderedcolors("gem");
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            origCorrCoeff = origRunner.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC.';
            null1CorrCoeff = nullRunner1.HoldoutTestAUCMeans;
            null2CorrCoeff = nullRunner2.HoldoutTestAUCMeans;
            null3CorrCoeff = nullRunner3.HoldoutTestAUCMeans;
            intervalWidth = 4;
            colorlist = options.ColorList;
            h1 = RandomConnectomeTestVisualizer.intervalSwarm(origCorrCoeff,intervalWidth,1,"Color",colorlist(1,:),"DimSHow",1);
            hold on
            h2 = RandomConnectomeTestVisualizer.intervalSwarm(null1CorrCoeff,intervalWidth,2,"Color",colorlist(2,:),"DimSHow",1);
            hold on
            h3 = RandomConnectomeTestVisualizer.intervalSwarm(null2CorrCoeff,intervalWidth,3,"Color",colorlist(3,:),"DimSHow",1);
            hold on
            h4 = RandomConnectomeTestVisualizer.intervalSwarm(null3CorrCoeff,intervalWidth,4,"Color",colorlist(4,:),"DimSHow",1);
            hold on
            ylabel("AUC")
            yticks([0.5:0.1:1])
            ylim([0.45,1])
            xticks([])
            xticklabels([])
            xlabel("")
            box off
            set(f,'Color','w')
        end

        function f = swarmNullModelsPISim_3nulls(nullRunner1,nullRunner2,nullRunner3,options)
            % fs9a(3)
            arguments
                nullRunner1 NullModelAnalysisBatchRunner
                nullRunner2 NullModelAnalysisBatchRunner
                nullRunner3 NullModelAnalysisBatchRunner
                options.SwarmSize = 1.5;
                options.FWidth = 400;
                options.FDepth = 300;
                options.DimShow = 5;
                options.ColorList = orderedcolors("gem");
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            null1CorrCoeff = nullRunner1.WirigPISimilarities(:,1:options.DimShow);
            null2CorrCoeff = nullRunner2.WirigPISimilarities(:,1:options.DimShow);
            null3CorrCoeff = nullRunner3.WirigPISimilarities(:,1:options.DimShow);
            intervalWidth = 4;
            colorlist = options.ColorList;
            h2 = RandomConnectomeTestVisualizer.intervalSwarm(null1CorrCoeff,intervalWidth,1,"Color",colorlist(2,:));
            hold on
            h3 = RandomConnectomeTestVisualizer.intervalSwarm(null2CorrCoeff,intervalWidth,2,"Color",colorlist(3,:));
            hold on
            h4 = RandomConnectomeTestVisualizer.intervalSwarm(null3CorrCoeff,intervalWidth,3,"Color",colorlist(4,:));
            hold on
            ylabel("Wiring PI similarity")
            yticks([0:0.2:1])
            ylim([0,1])
            tickStart = (intervalWidth)/2;
            xticks([tickStart:intervalWidth:tickStart+intervalWidth*(options.DimShow-1)])
            xticklabels([1:options.DimShow])
            xlabel("Correlation Components (Rank)")
            box off
            set(f,'Color','w')
        end

        

        % Figure S12
        function f = histogramRegionPairDistance(crossInfo,connMat,options)
            % fs12a
            arguments
                crossInfo CrossRegionInformation
                connMat ConnectionMatrix
                options.binEdges = 0:100:12500;
                options.FWidth = 600;
                options.FDepth = 400;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            distanceMatrix = crossInfo.DistanceMatrix;
            distanceList_connections = distanceMatrix(connMat.Matrix .* crossInfo.DomainDefineMatrix == 1);
            distanceList_regionPairs = distanceMatrix(crossInfo.DomainDefineMatrix == 1);
            binEdges = options.binEdges;
            histogram(distanceList_regionPairs,'BinEdges',binEdges,"LineWidth",0.01);
            hold on
            histogram(distanceList_connections,'BinEdges',binEdges);
            xlabel("Distance (μm)")
            xlim([0,12500])
            legend(["All region pairs","Connections"])
            box off
            set(f,'Color','w')
        end


        function f = histogramRegionPairDistance_perMRPair(crossInfo,connMat,options)
            % fs12b
            arguments
                crossInfo CrossRegionInformation
                connMat ConnectionMatrix
                options.binEdges = 0:100:12500;
                options.FWidth = 2000;
                options.FDepth = 1500;
                options.MRLebel = 0;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            t = tiledlayout(13,13);
            distanceMatrix = crossInfo.DistanceMatrix;
            binEdges = options.binEdges;
            MRNameList = unique(crossInfo.SourceRegionInfo.getMajorRegions,'stable');
            ctr = 0;
            for sMRInd = 1:13
                sMR = MRNameList(sMRInd);
                for tMRInd = 1:13
                    tMR = MRNameList(tMRInd);
                    ctr = ctr + 1;
                    MRPairMask = crossInfo.makeSpecificMRPairMask(sMR,tMR);
                    if sMR == tMR
                        ax = nexttile(ctr);
                        set(ax,'XTick',[])
                        set(ax, 'XColor', [1 1 1])  
                        set(ax,'YTick',[])
                        set(ax, 'YColor', [1 1 1])  
                        %axis off

            
                    else
                        ax = nexttile(ctr);
                        distanceList_connections = distanceMatrix(connMat.Matrix .* MRPairMask == 1);
                        distanceList_regionPairs = distanceMatrix(MRPairMask == 1);
                        histogram(distanceList_regionPairs,'BinEdges',binEdges, ...
                            'EdgeColor', 'none', ...              
                            'FaceAlpha', 0.5);              
                        hold on
                        if ~isempty(distanceList_connections)
                            histogram(distanceList_connections,'BinEdges',binEdges, ...
                            'EdgeColor', 'none', ...              
                            'FaceAlpha', 1);
                            titleStr = strcat("#Connections: ",num2str(sum(connMat.Matrix .* MRPairMask == 1,"all")));
                            title(titleStr)
                        else
                            title("No Connections")
                        end
                        box off
                        %titleStr = strcat(sMR," -> ",tMR);
                        xlim([0,12500])
                                  
                    end
                    if options.MRLebel == 1
                        if tMRInd == 1
                            ylabel(ax,sMR,'Color','k')
                        end
                        if sMRInd == 13
                            xlabel(ax,tMR,'Color','k')
                        end
                    end
                end
            end
            if options.MRLebel == 1
                xlabel(t,"Target")
                ylabel(t,"Source")
            end
            set(f,'Color','w')
        end


        % Figure S13: reuse S9

        % Figure S10
        function f = swarmNullModelsCorrCoeff_1nulls(origRunner,nullRunner1,options)
            % fs10(1)
            arguments
                origRunner ConnectomeAnalysisRunner
                nullRunner1 NullModelAnalysisBatchRunner
                options.SwarmSize = 1.5;
                options.FWidth = 450;
                options.FDepth = 300;
                options.DimShow = 5;
                options.ColorList = orderedcolors("gem");
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            origCorrCoeff = origRunner.OverallModelAndResults.HoldoutSummary.TestAllCorrelations(:,1:options.DimShow);
            null1CorrCoeff = nullRunner1.HoldoutTestCorrMeans(:,1:options.DimShow);
            intervalWidth = 3;
            colorlist = options.ColorList;
            h1 = RandomConnectomeTestVisualizer.intervalSwarm(origCorrCoeff,intervalWidth,1,"Color",colorlist(1,:),"SwarmSize",options.SwarmSize);
            hold on
            h2 = RandomConnectomeTestVisualizer.intervalSwarm(null1CorrCoeff,intervalWidth,2,"Color",colorlist(2,:));
            ylabel("Correlation coefficient")
            yticks([0:0.2:1])
            ylim([-0.1,1])
            tickStart = (intervalWidth)/2;
            xticks([tickStart:intervalWidth:tickStart+intervalWidth*(options.DimShow-1)])
            xticklabels([1:options.DimShow])
            xlabel("Correlation Components (Rank)")
            box off
            set(f,'Color','w')
        end

        function f = swarmNullModelsAUC_1nulls(origRunner,nullRunner1,options)
            % fs10(2)
            arguments
                origRunner ConnectomeAnalysisRunner
                nullRunner1 NullModelAnalysisBatchRunner
                
                options.SwarmSize = 1.5;
                options.FWidth = 200;
                options.FDepth = 300;
                options.DimShow = 1;
                options.ColorList = orderedcolors("gem");
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            origCorrCoeff = origRunner.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC.';
            null1CorrCoeff = nullRunner1.HoldoutTestAUCMeans;
            intervalWidth = 2;
            colorlist = options.ColorList;
            h1 = RandomConnectomeTestVisualizer.intervalSwarm(origCorrCoeff,intervalWidth,1,"Color",colorlist(1,:),"DimSHow",1);
            hold on
            h2 = RandomConnectomeTestVisualizer.intervalSwarm(null1CorrCoeff,intervalWidth,2,"Color",colorlist(2,:),"DimSHow",1);
            ylabel("AUC")
            yticks([0.5:0.1:1])
            ylim([0.45,1])
            xticks([])
            xticklabels([])
            xlabel("")
            box off
            set(f,'Color','w')
        end

        function f = swarmNullModelsPISim_1nulls(nullRunner1,options)
            % fs10(3)
            arguments
                nullRunner1 NullModelAnalysisBatchRunner
                options.SwarmSize = 1.5;
                options.FWidth = 400;
                options.FDepth = 300;
                options.DimShow = 5;
                options.ColorList = orderedcolors("gem");
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            null1CorrCoeff = nullRunner1.WirigPISimilarities(:,1:options.DimShow);
            intervalWidth = 1;
            colorlist = options.ColorList;
            h2 = RandomConnectomeTestVisualizer.intervalSwarm(null1CorrCoeff,intervalWidth,1,"Color",colorlist(2,:));
            hold on
            ylabel("Wiring PI similarity")
            yticks([0:0.2:1])
            ylim([0,1])
            xticks(1:5)
            xticklabels([1:options.DimShow])
            xlabel("Correlation Components (Rank)")
            box off
            set(f,'Color','w')
        end

        % Figure S11
        function f = plotCorrVariousPC(variousPCTestResults,options)
            % fs11(1)
            arguments
                variousPCTestResults (1,:) OverallAnalysisGroup
                options.LineWidth = 2;
                options.MarkerSize = 12;
                options.PlotDim = 10;
                options.FWidth = 800;
                options.FDepth = 300;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            t = tiledlayout(1,3);
            lineWidth = options.LineWidth;
            mkSz = options.MarkerSize;
            pltDim = options.PlotDim;
            nexttile
            for l = 1:length(variousPCTestResults)
                l_model = variousPCTestResults(l);
                rList  = l_model.FullAnalysisModel.CCAResultsData.r(1,1:pltDim);
                h = plot(rList,'.-','LineWidth',lineWidth,'MarkerSize',mkSz);
                ylim([-0.1,1.0001])
                xlim([1,10])
                hold on
            end
            axis square
            box off
            xlabel("Correlaton Component (Rank)")
            ylabel("Correlation coefficient")
            nexttile
            for l = 1:length(variousPCTestResults)
                l_model = variousPCTestResults(l);
                rList  = l_model.HoldoutSummary.TrainMeanCorrelation(1,1:pltDim);
                h = plot(rList,'.-','LineWidth',lineWidth,'MarkerSize',mkSz);
                ylim([-0.1,1.0001])
                xlim([1,10])
                hold on
            end
            axis square
            box off
            xlabel("Correlaton Component (Rank)")
            ylabel("Correlation coefficient")
            nexttile
            for l = 1:length(variousPCTestResults)
                l_model = variousPCTestResults(l);
                rList  = l_model.HoldoutSummary.TestMeanCorrelation(1,1:pltDim);
                rList = rList(1:10);
                h = plot(rList,'.-','LineWidth',lineWidth,'MarkerSize',mkSz);
                ylim([-0.1,1.0001])
                xlim([1,10])
                hold on
            end
            axis square
            box off
            xlabel("Correlaton Component (Rank)")
            ylabel("Correlation coefficient")
            % legend
            legendlabel = ["10","25","50","75","100"];
            legend(legendlabel)
            set(f,'Color','w')
        end

        function f = barAUCVariousPC(variousPCTestResults,options)
            % fs11(2)
            arguments
                variousPCTestResults (1,:) OverallAnalysisGroup
                options.FWidth = 800;
                options.FDepth = 300;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            t = tiledlayout(1,3);
            colorList = orderedcolors("gem");
            NVariousPC = length(variousPCTestResults);
            AUCListWhole = zeros(NVariousPC,1);
            AUCListTrain = zeros(NVariousPC,1);
            AUCListTest = zeros(NVariousPC,1);
            for l = 1:length(variousPCTestResults)
                AUCListWhole(l,1) = variousPCTestResults(l).FullAnalysisModel.ReconstructionResultsData.AUC_ROC;
                AUCListTrain(l,:) = variousPCTestResults(l).HoldoutSummary.TrainMeanAUC_ROC;
                AUCListTest(l,:) = variousPCTestResults(l).HoldoutSummary.TestMeanAUC_ROC;
            end
            nexttile
            b = bar(AUCListWhole);
            b.FaceColor = 'flat';
            b.CData = colorList(1:NVariousPC,:);
            axis square
            box off
            ylim([0.5 1])
            xticklabels([10 25 50 75 100])
            xlabel("Number of PC")
            ylabel("AUC")
            nexttile
            b = bar(mean(AUCListTrain,2));
            b.FaceColor = 'flat';
            b.CData = colorList(1:NVariousPC,:);
            ylim([0.5 1])
            xticklabels([10 25 50 75 100])
            axis square
            box off
            xlabel("Number of PC")
            ylabel("AUC")
            nexttile
            b = bar(mean(AUCListTest,2));
            b.FaceColor = 'flat';
            b.CData = colorList(1:NVariousPC,:);
            ylim([0.5 1])
            xticklabels([10 25 50 75 100])
            axis square
            box off
            xlabel("Number of PC")
            ylabel("AUC")
            set(f,'Color','w')
        end


        % Legends only
        function f = makeLegendMRs(brainInfo)
            arguments
                brainInfo BrainRegionInformation
            end
            MRNameList = brainInfo.getMajorRegionInfo;
            load('colorList13.mat')
            fwidth = 200;
            fdepth = 250;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            makeLegendsByDummyPlot(MRNameList,colorList13,'Marker','square');
            set(f,'Color','w')
        end

        function f = makeLegendConnectome()
            labelNameList = ["Connected","Not connected","Excluded from analysis"];
            colorList = zeros(3,3);
            cmap = parula;
            colorList(1,:) = cmap(256,:);
            colorList(2,:) = cmap(1,:);
            colorList(3,:) = [0 0 0];
            fwidth = 200;
            fdepth = 100;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            makeLegendsByDummyPlot(labelNameList,colorList,'Marker','square');
            set(f,'Color','w')
        end

        function f = makeLegendRandomTest()
            labelNameList = ["Original connectome data","Golobally randomized model","Locally randomized model"];
            colorList = orderedcolors("gem");
            colorList = colorList(1:3,:);
            fwidth = 200;
            fdepth = 100;
            f = figure('Position',FigMaker.figPosition(fwidth,fdepth));
            makeLegendsByDummyPlot(labelNameList,colorList,'Marker','o',"MarkerSize",8);
            set(f,'Color','w')
        end

        function f = makeLegend_any(labelNameList,colorList,options)
            arguments
                labelNameList
                colorList (:,3) double
                options.FWidth = 200;
                options.FDepth = 100;
                options.MarkerSize = 12;
            end
            f = figure('Position',FigMaker.figPosition(options.FWidth,options.FDepth));
            makeLegendsByDummyPlot(labelNameList,colorList,'Marker','o',"MarkerSize",options.MarkerSize);
            set(f,'Color','w')
        end



    end
end
                    