classdef SAUCY < handle
    % SAUCY Object class to manage analysis of single file
    %   SAUCY allows analysis for a single file and can be used for
    %   threshold detection, PCA classification, and optimization of false
    %   positives and false negatives.
    
    properties
        experiment_name = '' % base name of experiment
        
        fname_intan = '' % full file name of intan file
        fname_cbin = '' % full file name of cbin file
        fname_neuralnot = '' % full file name of neural_not file
        
        % chan          Channel number(s) for analysis
        %   as INT      Single channel will be used
        %   as TUPLE    Channel subtraction will be used (a b) = a - b
        chan = []
        
        n_clusters = 2 % number of clusters to detect
        
        Fs = 30000 % sampling rate of recording
        
        data = {} % Struct to store data
    end
    
    methods
        function load_data()
           
        end
    end
    
end

