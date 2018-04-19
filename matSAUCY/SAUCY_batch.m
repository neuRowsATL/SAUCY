classdef SAUCY_batch < handle
    %meta_SAUCY Manages multiple SAUCY objects for analysis
    %   Large datasets can be managed across multiple data files using
    %   meta_SAUCY
    
    properties
        % mode      What mode to run analysis?
        %   batch       Run across multiple files simultaneously
        %   track       Track analysis results across multiple files
        %   check       Check analysis results
        %   add_latent  Add latent values to analysis results?
        mode = 'batch'
    end
    
    methods
    end
    
end

