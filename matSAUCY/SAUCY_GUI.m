%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BC - April 17, 2018
%%%%% GUI for SAUCY algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SAUCY_GUI
    %SAUCY_GUI GUI to use SAUCY algorithm for spike sorting
    %   SAUCY allows a user to detect and cluster spikes using an
    %   interactive framework.
    
    properties
        gui_name % name of gui instance
        
        master_window % handle to master figure window
        
        tab_files % tab to manage data files
        tab_wavArray % tab for showing array data
        tab_spikeDetection % tab for showing spike detection
        tab_spikeClustering % tab for showing spike clustering
        tab_validation % tab to show validation plots
        
        saucy_objs = []; % list of SAUCY objects
    end
    
    methods
        function gui = SAUCY_GUI(gui_namDe)
            gui.gui_name = gui_name;
            fig = uifigure('Visible', 'Off', 'Position', [50, 50, 1500, 900]);
            fig.Name = gui_name;
            set(fig, 'MenuBar', 'none');
            set(fig, 'ToolBar', 'none');
            
            tgroup = uitabgroup('Parent', fig);
            
            gui.tab_files = uitab('Parent', tgroup, 'Title', 'Data Files');
            gui.tab_wavArray = uitab('Parent', tgroup, 'Title', 'Array Data');
            gui.tab_spikeDetection = uitab('Parent', tgroup, 'Title', 'Detection');
            gui.tab_spikeClustering = uitab('Parent', tgroup, 'Title', 'Clustering');
            gui.tab_validation = uitab('Parent', tgroup, 'Title', 'Validation');
            
            gui.master_window = fig;
            movegui(fig, 'center');
            fig.Visible = 'on';
        end
    end
end

