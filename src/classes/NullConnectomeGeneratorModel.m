classdef NullConnectomeGeneratorModel
    % Generates randomized connectome matrices using a specified method
    properties
        OriginalConnectionMatrix ConnectionMatrix
        CrossRegionInfo CrossRegionInformation
        RandomizationMethod function_handle
        UseBlockwise logical
        Randomizer MajorRegionBlockRandomizer
        Description 
    end

    methods
        function obj = NullConnectomeGeneratorModel(connMatrix, crossInfo, methodName, useBlockwise,options)
            arguments
                connMatrix ConnectionMatrix
                crossInfo CrossRegionInformation
                methodName string
                useBlockwise logical
                options.Description = "";
            end

            obj.OriginalConnectionMatrix = connMatrix;
            obj.CrossRegionInfo = crossInfo;
            obj.RandomizationMethod = ConnectomeRandomizationMethods.randMethodName2handle(methodName);
            obj.UseBlockwise = useBlockwise;
            obj.Randomizer = MajorRegionBlockRandomizer(obj.RandomizationMethod,obj.UseBlockwise);
            obj.Description = options.Description;
        end

        function [randomizedConn,newDomDef] = generate(obj,options)
            arguments
                obj NullConnectomeGeneratorModel
                options.DistanceBinWidth = [];
            end
            [randomizedConn,newDomDef] = obj.Randomizer.apply(obj.OriginalConnectionMatrix,obj.CrossRegionInfo,"DistanceBinWidth",options.DistanceBinWidth);
        end
    end
end
