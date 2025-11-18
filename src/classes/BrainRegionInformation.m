classdef BrainRegionInformation < handle
    % BrainRegionInformation
    % Handles brain region metadata and keeps AnnotationIndex for linkage to BrainSpace3D

    properties
        RegionTable table
    end

    methods
        function obj = BrainRegionInformation(tbl)
            arguments
                tbl table
            end

            % Check required fields
            requiredFields = {'atlas_id', 'acronym', 'major_region', 'name'};
            if ~all(ismember(requiredFields, tbl.Properties.VariableNames))
                error('Input table must contain atlas_id, acronym, major_region, and name.');
            end

            % Assign AnnotationIndex if missing
            if ~ismember('AnnotationIndex', tbl.Properties.VariableNames)
                tbl.AnnotationIndex = (1:height(tbl))';
            end

            obj.RegionTable = tbl;
        end

        function ids = getAtlasIDs(obj)
            ids = obj.RegionTable.atlas_id;
        end

        function acronyms = getRegionAcronyms(obj)
            acronyms = obj.RegionTable.acronym;
        end

        function majorRegions = getMajorRegions(obj)
            majorRegions = obj.RegionTable.major_region;
        end

        function names = getRegionNames(obj)
            names = obj.RegionTable.name;
        end

        function annotationIdx = getAnnotationIndices(obj)
            annotationIdx = obj.RegionTable.AnnotationIndex;
        end

        function [MRNames,MRInitialIndexList,MRIndexList] = getMajorRegionInfo(obj)
            [MRNames,MRInitialIndexList,MRIndexList] = unique(obj.RegionTable.major_region,'stable');
        end


        function newObj = subsetByAnnotationIndex(obj, idxList)
            mask = ismember(obj.RegionTable.AnnotationIndex, idxList);
            newTbl = obj.RegionTable(mask, :);
            newObj = BrainRegionInformation(newTbl);
        end

        function newObj = subsetByMajorRegion(obj, majorRegionName)
            mask = strcmp(obj.RegionTable.major_region, majorRegionName);
            newTbl = obj.RegionTable(mask, :);
            newObj = BrainRegionInformation(newTbl);
        end

        function [cListPerRegion] = makeMRColorListPerRegion(obj)
            load("colorList13.mat");
            [~,~,MRIndexList] = getMajorRegionInfo(obj);
            NRegion = height(MRIndexList);
            cListPerRegion = zeros(NRegion,3);
            for nr = 1:NRegion
                MRIndex = MRIndexList(nr);
                cListPerRegion(nr,:) = colorList13(MRIndex,:);
            end
        end


        function [ax] = setMRColorLabel(obj,options)
            arguments
                obj BrainRegionInformation
                %ax 
                options.Vertical = 'on';
            end
            load('colorList13.mat');
            categoryColorList = colorList13;
            [~,~,MRIndexList] = getMajorRegionInfo(obj);
            switch options.Vertical
                case 'off'
                    MRIndexList = MRIndexList.';
            end
            ax = axes;
            imagesc(ax,MRIndexList)
            colormap(ax,categoryColorList)
            axis off
        end

        function [shuffledRegionIndex] = getShuffledRegionIndexWithinMRs(obj)
            arguments
                obj BrainRegionInformation
            end
            [~,~,MRIndexList] = getMajorRegionInfo(obj);
            NRegion = numel(MRIndexList);
            shuffledRegionIndex = zeros(NRegion,1);
            [MRIDs,~,~] = unique(MRIndexList);
            for nmr = 1:length(MRIDs)
                idx = find(MRIndexList == nmr);
                shuffledRegionIndex(idx) = idx(randperm(length(idx)));
            end
        end
    end
end