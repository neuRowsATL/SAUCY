%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BC - April 17, 2018
%%%%% GUI for SAUCY algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SAUCY_GUI
    %SAUCY_GUI GUI to use SAUCY algorithm for spike sorting
    %   SAUCY allows a user to detect and cluster spikes using an
    %   interactive framework.
    
    properties
        saucy_objs = []; % list of SAUCY objects
    end
    
    methods
        function obj = untitled2(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

