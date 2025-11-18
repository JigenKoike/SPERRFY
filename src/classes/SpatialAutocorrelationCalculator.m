classdef SpatialAutocorrelationCalculator

    methods (Static)

        function moransIMatrix = calculateMoransIMatrix(valueMatrix,distanceMatrix,hList)
            arguments(Input)
                valueMatrix (:,:) double % [NRegion, DimComponents]
                distanceMatrix (:,:) double % [NRegion, NRegion]
                hList (:,1) double % [LhList,1]
            end
            arguments(Output)
                moransIMatrix (:,:) double % [DimComponents, LhList]
            end
            [~,DimComponents] = size(valueMatrix);
            LhList = length(hList);
            moransIMatrix = zeros(DimComponents,LhList);
            for d = 1:DimComponents
                valueVec = valueMatrix(:,d);
                moransIList = SpatialAutocorrelationCalculator.makeMoransIList(valueVec,distanceMatrix,hList);
                moransIMatrix(d,:) = moransIList';
            end
        end

        function moransIMatrix = calculateMoransIMatrix_forSurrogates(surrogatedValueMatrices,distanceMatrix,h_single,options)
            arguments(Input)
                surrogatedValueMatrices  (:,:,:) double % [NRegion, DimComponents, NSample]
                distanceMatrix (:,:) double % [NRegion, NRegion]
                h_single (1,1) double % 
                options.PrintProgress (1,1) logical = false;
                options.CounterInterval = 100;
            end
            arguments(Output)
                moransIMatrix (:,:) double % [DimComponents,NSample]
            end
            [~,DimComponents,NSample] = size(surrogatedValueMatrices);
            moransIMatrix = zeros(DimComponents,NSample);
            for d = 1:DimComponents
                if options.PrintProgress == true
                    if mod(d,options.CounterInterval) == 1;
                        disp(strcat("d = ",num2str(d)));
                    end
                end
                for ns = 1:NSample
                    valueVec = surrogatedValueMatrices(:,d,ns);
                    moransI = SpatialAutocorrelationCalculator.makeMoransIList(valueVec,distanceMatrix,h_single);
                    moransIMatrix(d,ns) = moransI;
                end
            end
        end


        function moransIList = makeMoransIList(valueVec,distanceMatrix,hList)
            moransIList = zeros(size(hList));          
            for nList = 1:length(hList)
                h = hList(nList);
                W = SpatialAutocorrelationCalculator.distance2weight_exponential(distanceMatrix,h);
                moransI = SpatialAutocorrelationCalculator.calculateIfromW(valueVec,W);
                if isnumeric(moransI) == 1
                    moransIList(nList) = SpatialAutocorrelationCalculator.calculateIfromW(valueVec,W);
                else
                    moransIList(nList) = NaN;
                    dispStr = strcat("#List = ",num2str(nList),", h = ",num2str(round(h,4))," : unstable computaion");
                    disp(dispStr);
                end 
            end
        end

        function moransI = calculateIfromW(x,W)
            if isnumeric(W) == 1
                N = width(W);
                x_mean = mean(x);
                z = x - x_mean;
                W_sum = sum(W,'all');
                numerator = z' * W * z;
                denominator = z' * z;
                moransI = (N/W_sum) * (numerator/denominator);
            else
                moransI = "unstable computation";
            end
        end
            

        function weightMat = distance2weight_exponential(D,h)
            epsilon = 10^(-6);
            N = width(D);
            f_exp = @(x) exp(-x/h);
            weightMat_beforeNormalize = arrayfun(f_exp,D);
            weightMat_beforeNormalize = weightMat_beforeNormalize .* (ones(N)-eye(N));
            weightMatSum = sum(weightMat_beforeNormalize);
            if min(weightMatSum) > epsilon
                weightMat = weightMat_beforeNormalize ./ repmat(weightMatSum,size(weightMatSum.'));
            else
                weightMat = "unstable computation";
            end
        end

        function distMedian = gethFromMedian(D)
            % std,mean,median,q25,q75
            DList = D(D>0);
            distMedian = median(DList);
        end

        function hList = makehList_exponentialScale(h,base,exponent_min,exponent_max)
            NList = (exponent_max - exponent_min) + 1;
            hList = zeros(NList,1);
            for n = 1:NList
                hList(n) = h * base^(n + 1 + exponent_min);
            end
        end
            
        function hList = makehList_quantile(D,quantileWidth)
            DList = D(D>0);
            quantileList = quantileWidth:quantileWidth:1;
            NList = length(quantileList);
            hList = zeros(NList,1);
            for n = 1:NList
                p = quantileList(n);
                hList(n) = quantile(DList,p);
            end
                
        end

        function exponentialWeightMat = distance2exponentialWeight(D,h,normalizeByRowTag)
            epsilon = 10^(-6);
            N = width(D);
            f_exp = @(x) exp(-x/h);
            weightMat_beforeNormalize = arrayfun(f_exp,D);
            weightMat_beforeNormalize = weightMat_beforeNormalize .* (ones(N)-eye(N));
            weightMatSum = sum(weightMat_beforeNormalize);
            if min(weightMatSum) > epsilon
                if normalizeByRowTag == 1
                    weightMat = weightMat_beforeNormalize ./ repmat(weightMatSum,size(weightMatSum.'));
                elseif normalizeByRowTag == 0
                    weightMat = weightMat_beforeNormalize;
                end
                weightMatMax = max(weightMat,[],'all');
                exponentialWeightMat = weightMat/ weightMatMax;
            else
                exponentialWeightMat = "unstable computation";
            end
        end


    end

end
