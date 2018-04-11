classdef SAUCY < handle
    % SAUCY Object class to manage analysis of single file
    %   SAUCY allows analysis for a single file and can be used for
    %   threshold detection, PCA classification, and optimization of false
    %   positives and false negatives.
    
    properties
        experiment_name = '' % base name of experiment
        
        data_path = '' % path to data files
        fname_intan = '' % full file name of intan file
        fname_cbin = '' % full file name of cbin file
        fname_neuralnot = '' % full file name of neural_not file
        
        % chan          Channel number(s) for analysis
        %   as INT      Single channel will be used
        %   as TUPLE    Channel subtraction will be used (a b) = a - b
        chan = []
        
        n_clusters = 2 % number of clusters to detect
        
        Fs = 30000 % sampling rate of recording
        F_low = 300 % low cut off
        F_high = 7500 % hi cut off
        
        raw_data = {} % Struct for raw data 
        data = {} % Struct for results
    end
    
    methods
        % Initialize SAUCY object
        function S = SAUCY(experiment_name)
            S.experiment_name = experiment_name;
        end
        
        % Load data into SAUCY object
        function load_data(S, filename)
            if strfind(filename, '.rhd') > 0
               [t_amplifier, t_board_adc, amplifier_data, board_adc_data, frequency_parameters] = read_Intan_RHD2000_nongui_saucy(RHD_name);
               S.raw_data.t_amplifier = t_amplifier;
               S.raw_data.t_board_adc = t_board_adc;
               S.raw_data.amplifier_data = amplifier_data;
               S.raw_data.board_adc_data = board_adc_data;
               
               S.Fs = frequency_parameters.amplifier_sample_rate;
            elseif strfind(filename, '.cbin' > 0
                % use cbin script
            else
                % use default script
            end
        end
    end
    
end

