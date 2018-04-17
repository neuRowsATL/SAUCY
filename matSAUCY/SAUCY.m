classdef SAUCY < handle
    % SAUCY Object class to manage analysis of single file
    %   SAUCY allows analysis for a single file and can be used for
    %   threshold detection, PCA classification, and optimization of false
    %   positives and false negatives.
    
    properties
        verbose = 0 % level of verbose output
        warning_str = [];
        
        experiment_name = ''; % base name of experiment
        
        data_path = ''; % path to data files
        fname_mat = ''; % full file name of .mat file
        fname_intan = ''; % full file name of .rhd (Intan) file
        fname_cbin = ''; % full file name of .cbin file
        fname_neuralnot = ''; % full file name of .neuralnot.mat file
        
        col_vec='rgbc';
        col_mat=[1 0 0;0 1 0;0 0 1;1 1 0];
        
        % chan          Channel number(s) for analysis
        %   as INT      Single channel will be used
        %   as TUPLE    Channel subtraction will be used [a b] = a - b
        chan = [];
        
        n_clusters = 2; % number of clusters to detect
        PCA = {};
        kMeans = {};
        
        Fs % sampling rate of recording
        F_low % low cut off
        F_high % hi cut off
        filter_type % type of filter to apply
        
        raw_data = {}; % struct for raw data 
        data = {}; % struct for results
        
        alignment % how are spike times aligned? peak or threshold?
        threshold = []; % Spike detection threshold
        spike_ids = [];
        nsamp_wave = [30 90]; % how many indices to include in spike_mat
        
        % Bookkeeping variables
        TH_original % for plotting purposes only
        inv_sign = 1; % indicates whether waveform is inverted
        inv_waveform = 1; % indicates whether waveform has been inverted or not
        
        spike_mat % matrix of spike waveforms
    end
    
    methods
        % ----- ----- ----- ----- -----
        % Initialize SAUCY object
        % Sets variables:
        %       experiment_name
        %       chan
        function S = SAUCY(experiment_name, chan)
            S.experiment_name = experiment_name;
            S.chan = chan;
            
            if S.verbose > 2
                disp(['Initiated SAUCY for ', experiment_name]);
                if length(chan) == 1
                    disp(['Using unipolar channel: ', num2str(chan)]);
                elseif length(chan) == 2
                    disp(['Using subtracted channels: ', num2str(chan(1)), ' - ' num2str(chan(2))]);
                else
                    error('Invalid channel number or subtraction specification.');
                end
            end
        end
        
        % ----- ----- ----- ----- -----
        % Get 1 second of sample data
        % Returns 1 second of sample data [t, wav]
        function [output] = get_second_samp(S)
            dat = S.data.amplifier_data_filt;
            Fs = S.Fs;
            
            % BC ----->> MIGHT NEED TO RECALIBRATE THIS IF SPIKES ARE SPARSE
            % IN FILE
            
            % Selects 1 second of data from the middle of the file
            % If data is shorter than 1 sec, uses all of data
            mid_dat=round(length(dat)/2)-Fs/2;
            if length(dat)/Fs > 1
                id_for_samp=mid_dat:(mid_dat+Fs-1);    
            else
                id_for_samp=1:length(dat);    
            end

            t_dat = S.raw_data.t_amplifier(id_for_samp);
            wav_dat = S.data.amplifier_data_filt(id_for_samp);
            
            output.ts = t_dat;
            output.wav = wav_dat;
            output.ixs = id_for_samp;
        end
        
        % ----- ----- ----- ----- -----
        % Load data into SAUCY object
        % Sets variables:
        %       raw_data = {t_amplifier, t_board_adc, amplifier_data,
        %           board_adc_data}
        %       data = {filt}
        %       Fs
        %       fname_intan
        function load_data(S, filename)
            if strfind(filename, '.rhd') > 0
                %%%%% ADD CODE TO SAVE .mat FILE FOR QUICK ACCESS
                basename = split(filename, '.');
                S.fname_intan = filename;
                S.fname_mat = strcat(basename{1}, '.mat');
                
                if ~exist(S.fname_mat, 'file')
                    disp(['No .mat file found.', newline]);
                    disp(['Converting .rhd file to .mat file', newline]);
                    read_Intan_nongui(filename);
                end
                
                % Load data from .rhd file
%                 [t_amplifier, t_board_adc, amplifier_data, board_adc_data, frequency_parameters] = ...
%                    read_Intan_RHD2000_nongui_saucy(filename);
                disp(['Loading data...',newline]);
                dat_tmp = load(S.fname_mat, 't_amplifier', 't_board_adc', 'amplifier_data', 'board_adc_data', 'frequency_parameters');
                S.raw_data.t_amplifier = dat_tmp.t_amplifier;
                S.raw_data.t_board_adc = dat_tmp.t_board_adc;
                S.raw_data.amplifier_data = dat_tmp.amplifier_data;
                S.raw_data.board_adc_data = dat_tmp.board_adc_data;
                
                % Set default data for analysis
                S.data.amplifier_data_filt = S.raw_data.amplifier_data(S.chan,:);
               
                S.Fs = dat_tmp.frequency_parameters.amplifier_sample_rate;
               
               
                if S.verbose > 1
                    disp(['Successfully loaded Intan file: ', filename]);
                end
            elseif strfind(filename, '.cbin') > 0
                % Load .cbin file
                
                if S.verbose > 1
                    disp(['Not prepared to load a cbin file! (', filename, ')']);
                end
            elseif strfind(filename, '.mat') > 0
                % Load data from .mat file
                
                if S.verbose > 1
                    disp(['Not equipped to open a mat file! (', filename, ')']);
                end
            else
                error(["I don't know what to do with your file type... (", filename, ")"]);
            end
        end
        
        % ----- ----- ----- ----- -----
        % Filter data if specified
        % Sets variables:
        %       F_low
        %       F_high
        %       filter_type
        function filter_data(S, f_specs, filter_type)
            if nargin == 1
                disp('!! USING DEFAULT BANDPASS LIMITS [350 7000] !!');
                disp('!! USING DEFAULT FILTER TYPE: hanningfir !!');
                S.F_low = 350;
                S.F_high = 7000;
                S.filter_type = 'hanningfir';
            elseif nargin == 2
                if length(f_specs) == 2
                    S.F_low = f_specs(1);
                    S.F_high = f_specs(2);
                end
                disp('!! USING DEFAULT FILTER TYPE: hanningfir !!');
                S.filter_type = 'hanningfir';
            elseif nargin == 3
                if length(f_specs) == 2
                    S.F_low = f_specs(1);
                    S.F_high = f_specs(2);
                end
                S.filter_type = filter_type;
            else
                error('Unable to figure out filter parameters.');    
            end
            
            disp(['Bandpass low cut off: ',num2str(S.F_low),' Hz']);
            disp(['Bandpass high cut off: ',num2str(S.F_high),' Hz']);
            disp(['Using sampling frequency: ',num2str(S.Fs),' s/s']);
            disp(['Using filter type: ', S.filter_type]);
            
            if isfield(S.raw_data, 'amplifier_data')
                S.data.amplifier_data_filt = ...
                    bandpass_filtfilt(S.raw_data.amplifier_data(S.chan, :), ...
                    S.Fs, S.F_low, S.F_high, S.filter_type);
            else
                disp('ERROR: Check that data has been loaded into SAUCY.');
            end
            
            disp(['Filtering complete...',newline]);
        end
        
        % ----- ----- ----- ----- -----
        % Set threshold for spike detection using SAUCY interface
        % Sets variables:  
        %       threshold, TH_original
        %       inv_sign
        %       inv_waveform
        function set_user_threshold(S, data_source, std_th)
            % data_source       filt        Uses filtered data
            %                   raw         Uses raw data
            
            % Set default S.D. for threshold
            if nargin == 1
                data_source = '';
                n_SD_for_TH = 2;
            elseif nargin == 2
                n_SD_for_TH = 2;
            elseif nargin == 3
                n_SD_for_TH = std_th;
            end
            
            if strcmp(data_source, 'filt') & ~strcmp(data_source, 'raw') & ~strcmp(data_source, '')
                if isfield(S.data, 'amplifier_data_filt')
                    disp('Using filtered data to find threshold.');
                    dat = S.data.amplifier_data_filt;
                else
                    disp('ERROR: Check that filter was applied to data.');
                end
            elseif strcmp(data_source, 'raw') | strcmp(data_source, '')
                if isfield(S.raw_data, 'amplifier_data')
                    disp('Using raw data to find threshold.');
                    dat = S.raw_data.amplifier_data(S.chan,:);
                else
                    disp('ERROR: Check that data was loaded into SAUCY.');
                end
            else
                disp('ERROR: Invalid data source format.');
            end
            
            % Plot waveform sample for user selection of threshold
            samp_dat = S.get_second_samp();
            ts = samp_dat.ts;
            wav_samp = samp_dat.wav;
            
            figure();
            plot(ts, wav_samp, 'k');
            set(gca,'xlim',[min(ts) max(ts)]);

            % Prompt user to select threshold
            t_str{1}='Left-click twice between noise and large spike peaks';
            t_str{2}=['or hit return to use mean - ' num2str(n_SD_for_TH) ' S.D.'];
            t_str{3}=['or right-click once to select a threshold by hand ' num2str(n_SD_for_TH) ' S.D.'];
            title(t_str,'fontweight','bold','fontsize',12)
            [tmp, TH_1, BUTTON] = ginput(1);
            
            if isempty(TH_1)
                TH = TH_1;
            elseif BUTTON==1    % Left-click
                % Plot threshold from first click
                xl_tmp = get(gca,'xlim');
                hold on;
                plot(xl_tmp,TH_1*[1 1], 'r:');
                
                [tmp, TH_2] = ginput(1);
                TH_lower=min([TH_1 TH_2]);
                TH_upper=max([TH_1 TH_2]);
                
                % pull waveforms between user-selected thresholds
                dat_clipped_tmp=dat(find(dat > TH_lower & dat < TH_upper));
                mean_above_clipped=mean(dat(find(dat > TH_upper)));
                mean_below_clipped=mean(dat(find(dat < TH_lower)));
                
                % if upward spikes
                if mean_above_clipped-mean(dat_clipped_tmp) > mean(dat_clipped_tmp)-mean_below_clipped
                    TH_sign=1;
                    TH_sign_str='+';
                else   % if downward spikes
                    TH_sign=-1;
                    TH_sign_str='-';
                end
                
                % Typical method: mean(wav) +/- 2*std(wav)
                TH_oldmethod = mean(dat) + TH_sign*n_SD_for_TH*std(dat);
                
                % New method: mean(btwn thresh) +/- 2*std(btwn thresh)
                % BC ----->> CONSIDER USING ADJUSTED STD based on Rossant 2016?
                TH = mean(dat_clipped_tmp)+TH_sign*n_SD_for_TH*std(dat_clipped_tmp);
                
                disp(['Adjusting threshold from MEAN +/- 2*SD = ' num2str(TH_oldmethod) ' to ' num2str(TH)])
                disp(['Computed TH = ' num2str(TH) ' (MEAN ' TH_sign_str ' ' num2str(n_SD_for_TH) '*SD)'])
                TH_set_by_user = TH;
            elseif BUTTON==3    % Right-click
                TH = TH_1;
                disp(['Hand-set TH = ' num2str(TH)])
                TH_set_by_user = TH;
            end

            if isempty(TH)
                TH = mean(dat)-n_SD_for_TH*std(dat);
                TH_set_by_user = TH;
                disp(['Computed TH = ' num2str(TH) ' (MEAN - ' num2str(n_SD_for_TH) '*SD)'])
            end
            
            S.TH_original = TH; % holds track of original threshold value
            S.threshold = TH;
            
            % BC ----->> WHY INVERT AND THEN UN-INVERT?!?
            % problems here
            if  TH < 0
                S.data.amplifier_data_filt = -S.data.amplifier_data_filt;
                S.inv_sign = -1;
                S.threshold = -TH;
                S.inv_waveform = 1;
                
                disp('Inverted waveform...')

                % ----->> POSSIBLE BUG NEEDS TO BE FIXED
                %     if batchmode_condition % this doesn't seem to solve problem consistently - screws up with alternate reloading
                %     vert_spike_lims=fliplr(-vert_spike_lims); disp('test - inverting spike lims too')
                %     end
            else
                S.inv_sign = 1;
                S.inv_waveform = 0;
            end
            
            close(gcf);
            
            figure();
            plot(S.raw_data.t_amplifier, dat);
            hold on;
            plot(xlim, [TH TH], 'r-');
            title('Press [ENTER] to continue...');
            
            x = input('Press [ENTER] to continue...', 's');
            close(gcf);

            % ----->> POSSIBLE BUG NEEDS TO BE FIXED
            % dat and TH are consistent when reloaded...  but looks like vertical lims
            % are bad
        end
        
        % ----- ----- ----- ----- -----
        % Extract spike waveforms
        % Sets variables:
        %       alignment
        function [spike_mat id_spikes] = set_spike_mat(S, alignment, spike_ids, vert_spike_lims)
            % Identify spike peaks (copied from neural_and_song.m)
            %   alignment           Specify how to align spikes
            %       peaks           Align spikes by peak
            %       threshold       Align spikes by threshold crossing
            %   spike_ids           IDs of spikes
            %   vert_spike_lims     Carried from original SAUCY
            
            % If no alignment mode is specified, ALIGN SPIKES TO PEAKS
            if nargin == 1
                S.alignment = 'peak';
            end

            dat = S.data.amplifier_data_filt;
            TH = S.threshold;
            Fs = S.Fs;
            nsamp_wave = S.nsamp_wave;
            
            invert_sign = S.inv_sign;
            
            if S.verbose > 1
                show_timing = true;
            else
                show_timing = false;
            end

            if nargin == 2 | (nargin == 4 & length(spike_ids) == 0)
                S.alignment = alignment;
                % Identify peaks of spikes
                % -----
                id = find(dat>TH);                % id of samples that are above threshold
                sdx = diff(sign(diff(dat)));    % this will be -2 one point before upward-pointing peak, +2 one point before downward-pointing peak, zero else

                id2 = find(sdx<0)+1;            % id2 is ids of upward-pointing peaks
                first_der = diff(dat);

                not_constant = find(first_der);   % first deriv not equal to zero - use this to prevent multiple spikes being assigned to truncated spikes

                % so to be identified as a spike,  sample must be:
                % (1) above threshold
                % (2) an upward-pointing peak and
                % (3) changing in value
                id_peaks = intersect(id,id2);     
                id_peaks = intersect(id_peaks,not_constant);

                if strcmp(S.alignment, 'threshold')
                    % aligns at threshold crossing, not peak
                    id_crossings = find(diff(dat>TH)==1)+1;           % id of samples that first to cross threshold

                    % # of peaks will be slightly higher than # crossings.  Below is a vector
                    % of id_peaks that fall immediately after each threshold crossing.  This is
                    % <<===== FOR PLOTTING ONLY =====>>
                    for x=1:length(id_crossings)
                        id_peaks_after_crossings(x,1)=min(id_peaks(find(id_peaks>id_crossings(x))));
                    end

                    id_peaks = id_crossings;
                end
            elseif nargin == 3
                id_peaks = spike_ids;
            end
                  
            % =====>> This section taken from generate_spike_mat.m
            % first, disinclude peaks too close to start or end of file
            id_peaks = id_peaks(find(id_peaks>nsamp_wave(1) & id_peaks<(length(dat)-nsamp_wave(2))));

            % id_before_peaks is a vector of the indexes of datapoints before the
            % peaks
            id_before_peaks = zeros(nsamp_wave(1),length(id_peaks));  % do memory allocation all at once
            id_before_peaks(1,:) = id_peaks-1;
            for x=2:nsamp_wave(1)
                id_before_peaks(x,:) = id_peaks-x;
            end
            
            % here, we orgainze these into a matrix:
            % each row is a sample number before peak
            % each column in for one entry in id_peaks
            id_before_peaks_mat=reshape(id_before_peaks,nsamp_wave(1),length(id_peaks));

            % id_before_peaks is a vector of the indexes of all datapoints within
            % nsampe_wave of the peaks
            id_whole_waveform=zeros(sum(nsamp_wave)+1,length(id_peaks));% do memory allocation all at once - speeds things up
            id_whole_waveform(1,:)=id_peaks-nsamp_wave(1);
            for x=[-1*nsamp_wave(1)+1:nsamp_wave(2)]
                id_whole_waveform(x+nsamp_wave(1)+1,:)=id_peaks+x;
            end

            % at least one point before peak must be below threshold - this finds
            % ids of the waveforms that contain no violations
            id_before_peaks_mat_AND_below_TH=find(sum(dat(id_before_peaks_mat)*invert_sign<TH));
            if vert_spike_lims(1)~=0
                id_no_vert_limit_violation=find(~sum(dat(id_whole_waveform)*invert_sign<vert_spike_lims(1)) & ~sum(dat(id_whole_waveform)*invert_sign>vert_spike_lims(2)));
            else % (if vertical limits not set, include all)
                [tmp,c]=size(id_before_peaks_mat);
                id_no_vert_limit_violation=1:c;
            end

            % ids of waveforms that pass both tests
            ids_all_crit=intersect(id_peaks(id_before_peaks_mat_AND_below_TH),id_peaks(id_no_vert_limit_violation));
            
            % assemble waveforms into a matrix
            id_waveforms_included=zeros(length(ids_all_crit),sum(nsamp_wave)+1);% do memory allocation all at once - speeds things up
            id_waveforms_included(:,1)=dat(ids_all_crit-nsamp_wave(1));
            for x=[-1*nsamp_wave(1)+1:nsamp_wave(2)]
                id_waveforms_included(:,x+nsamp_wave(1)+1)=dat(ids_all_crit+x);
            end

            % each row is a sample number before peak
            % each column in for one entry in id_peaks
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BC Updated: 2/21/2018
            %spike_mat=reshape(id_waveforms_included,length(id_waveforms_included),sum(nsamp_wave)+1);
            spike_mat = reshape(id_waveforms_included, size(id_waveforms_included,1), sum(nsamp_wave)+1);
            id_spikes = ids_all_crit';
        end
        
        % ===== ===== ===== ===== =====
        % Do clustering using kmeans
        function do_clustering(S, n_clusters, spike_ids)
            %  Un-inverting waveform
            if S.inv_sign == -1
                dat = -S.data.amplifier_data_filt;
            else
                dat = S.data.amplifier_data_filt;
            end

            TH = S.threshold;
            Fs = S.Fs;
            invert_sign = S.inv_sign;
            
            if nargin == 1
                n_clusters = S.n_clusters;
                spike_ids = S.spike_ids;
            elseif nargin == 2
                spike_ids = S.spike_ids;
            elseif nargin == 3
                disp(['Setting user-defined spike times',newline]);
            else
                error("Too many arguments in do_clustering()");
            end
                
                
            % Extract spike waveforms and establish vertical and horizontal
            % limits
            % -----
            col_vec = S.col_vec;
            col_mat = S.col_mat;

            if S.verbose > 1
                show_timing = true;
            else
                show_timing = false;
            end
            
            finished = false;
            
            vert_spike_lims_approved = false;
            horiz_spike_lims_approved = false;
            
            redo_peaks = true;
            replot = true;
            first_plotting = true;
            
            spike_ids = [];
            
            % [# samples before peak        # after peak] to subject to PCA (these are
            % defaults, user will be able to adjust these below)
            nsamp_wave = S.nsamp_wave;
            vert_spike_lims=[0 0]; % use vertical limits to cut off mvmt artifact (default is not to use these)
            id_peaks_save = [];        
            
            while ~finished
                if redo_peaks
                    % =====>> WHAT DOES THIS LINE DO?
                    % id_peaks = id_peaks(find(id_peaks>(nsamp_wave(1)+1)  & id_peaks<(length(dat)-nsamp_wave(2))));
                    
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end

                    if ~exist('spike_mat')  % a new trick - only run generate_spike_mat once
                        % vert_spike_lims are [0 0] when this is called
                        % When re-called in batchmode, vert limits are as they
                        % appear on disply - i think the negative of how they are
                        % actually used with the inverted waveform
                        disp([newline, 'Hack to preserve sign of vert_spike_lims when running in batchmode', newline])
                        [spike_mat id_peaks_save] = S.set_spike_mat(S.alignment, spike_ids, vert_spike_lims);
                    end
                    
                    save spike_mat_tmp spike_mat
                    
                    mean_wave = mean(spike_mat);
                    
                    if show_timing
                        disp(['Elapsed time generating spike_mat= ' num2str(toc)]);
                    end

                    % PRINCOMP(spike_mat) performs PCA on the N-by-P matrix of spikes  and returns the principal component coefficients
                    % Each row of spike_mat is a waveform   COEFF is a P-by-P matrix, each column containing coefficients
                    % for one principal component.  The columns are in order of decreasing   component variance.
                    %
                    % SCORE is the representation of X in the principal component space
                    % Rows of SCORE correspond to observations, columns to components
                    
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end
                    
                    % BC: 4/12/2018
                    % Modified: UPDATED FUNCTION CALL
                    % [COEFF, SCORE,LATENT] = princomp(spike_mat);
                    [COEFF, SCORE, LATENT] = pca(spike_mat);
                    S.PCA.coeff = COEFF;
                    S.PCA.score = SCORE;
                    S.PCA.latent = LATENT;
                    
                    if show_timing
                        disp(['Elapsed time running PCA= ' num2str(toc)]);
                    end

                    % %run this on the file used to set prefs to show that it works
                    % [n,p] = size(spike_mat);
                    % %subtract mean from each column of spike_mat
                    % spike_mat_x0 = spike_mat- repmat(mean(spike_mat,1),n,1);tmp_SCORE=spike_mat_x0*Data.COEFF_matrix;
                    % tmp_SCORE(1:5,1:5);SCORE(1:5,1:5)
            
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end
                    
                    % cluster points in component space, using only the first two
                    % components
                    
                    [IDX,C] = kmeans(SCORE(:,1:2), n_clusters);
                    S.kMeans.idx = IDX;
                    S.kMeans.C = C;
                    
                    if show_timing
                        disp(['Elapsed time grouping points= ' num2str(toc)]);
                    end

                    for x=1:n_clusters
                        wave_tmp=mean_wave+C(x,:)*COEFF(:,1:2)';
                        max_abs_wave(x)=max(abs(wave_tmp));
                    end
                    
                    % determine relative spike sizes to differentiate noise from spike
                    % waves
                    [tmp,wave_id_by_size]=sort(max_abs_wave);

                    %%%%%%%%% reorder clusters by spike size
                    % so noise_spikes are in cluster 1, and spikes are in clusters 2-n
                    IDX_old=IDX;
                    C_old=C;
                    IDX_ORIGINAL=IDX_old;
                    
                    clear IDX C
                    
                    for x=1:n_clusters
                        IDX(find(IDX_old==wave_id_by_size(x)))=x;
                        C(x,:)=C_old(wave_id_by_size(x),:);
                    end
                    
                    for x=1:n_clusters
                        % CHECK if this is right - component{x} will drift as mean centers
                        % do
                        % component{x} is the waveform reconstructed using only the
                        % mean of the first two principle components (cluster centers).
                        component{x} = mean_wave+C(x,:)*COEFF(:,1:2)';
                    end

                    %%%%%%%%%%%%%%%%%%%
                    for x=1:n_clusters
                        wave_tmp=mean_wave+C(x,:)*COEFF(:,1:2)';
                        max_abs_wave(x)=max(abs(wave_tmp));
                    end

                    redo_peaks = false;
                end
                % END OF "redo_peaks" CONDITIONTAL

                if replot
                    cla;
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end
                    
                    % Plot waveforms
                    n_waves_to_show=300;
                    
                    clear spike_mat_tmp
                    
                    % If vertical and horizontal limits are set, create
                    % final plot to show waveforms.
                    %
                    % If both limits are not set, plot sample of waveforms
                    % and prompt user to set vert/horiz limits.
                    % -----
                    if vert_spike_lims_approved & horiz_spike_lims_approved
                        close(gcf);
                        figure();
                        subplot(2,3,1);
                        
                        for z=1:n_clusters
                            id_tmp=find(IDX==z);
                            %                 id_tmp=id_tmp(1:min([n_waves_to_show length(id_tmp)]));
                            if n_waves_to_show >= length(id_tmp)
                                id_tmp=id_tmp(1: length(id_tmp));
                            else
                                [tmp,id_tmp_2]=sort(rand(1,length(id_tmp)));
                                id_tmp=id_tmp(id_tmp_2(1: length(id_tmp)));
                            end
                            id_tmp=id_tmp(1:min([n_waves_to_show length(id_tmp)]));

                            spike_mat_tmp{z}=spike_mat(id_tmp,:);
                            save_example_ids{z}=id_tmp;
                        end
                        
                        smt_color=col_mat(1:n_clusters,:)/2;
                        t_str{2}=['# waveforms - ' num2str(length(spike_mat)) ', only plotting ' num2str(n_waves_to_show) ' of each cluster'];
                    else
                        n_waves_to_show=1000;
                        smt_color=col_mat(1:n_clusters,:);
                        spike_mat_tmp{1}=spike_mat;
                        
                        % For very large files, plotting all spike waveforms will slow
                        % down MATLAB considerably.  This will select 1000 random
                        % waveforms and plot them.  Then, as the vertical and
                        % horizontal limits are restricted, it will continue plotting
                        % only waveforms that remain from that original 1000.
                        if length(spike_mat)>1000
                            if ~exist('id_peaks_save_subset')
                                a=1:length(spike_mat);b=randn(size(a));
                                [tmp,id]=sort(b);
                                id_waveform_subset=a(id(1:1000));
                                id_peaks_save_subset=id_peaks_save(id_waveform_subset);
                                spike_mat_tmp{1}=spike_mat(id_waveform_subset,:);
                            else
                                [tmp,id_waveform_subset_remaining,tmp]=intersect(id_peaks_save,id_peaks_save_subset);
                                spike_mat_tmp{1}=spike_mat(id_waveform_subset_remaining,:);
                            end
                        end
                        smt_color=[0 0 0];

                        t_str{2}=['# waveforms - ' num2str(length(spike_mat)) ' only showing ' num2str(length(spike_mat_tmp{1})) ', chosen randomly'];
                    end
                    
                    
                    xdat=[-nsamp_wave(1):nsamp_wave(2)]/(Fs./1000);% time in msec
                    
                    hold on;
                    for x=1:length(spike_mat_tmp)
                        plot(xdat,spike_mat_tmp{x}','color',smt_color(x,:));    % PLOTS SPIKE WAVEFORMS
                    end
                    
                    xlabel('Time (msec)');
                    if length(S.fname_intan) > 0
                        t_str{1}=RemoveUnderScore(S.fname_intan);
                    elseif length(S.fname_cbin) > 0
                        t_str{1}=RemoveUnderScore(S.fname_cbin);
                    elseif length(S.fname_neuralnot) > 0
                        t_str{1}=RemoveUnderScore(S.fname_neuralnot);
                    else
                        t_str{1}=RemoveUnderScore(S.experiment_name);
                    end
                    title(t_str)
                    
                    if first_plotting
                        orig_xlims=[min(xdat) max(xdat)];
                        first_plotting = false;
                    end
                    
                    set(gca,'xlim',orig_xlims,'xtick',[-1:.25:2]);xl=get(gca,'xlim');
                    yl=get(gca,'ylim');
                    
                    for x=1:2
                        plot(xl,vert_spike_lims(x)*[1 1],'k:');
                    end
                    
                    plot(xl,TH*invert_sign*[1 1],'k--','linew',2)

                    for x=1:n_clusters
                        plot(xdat,component{x},col_vec(x),'linew',2)
                    end
                    
                    replot = false;
                    
                    if show_timing
                        disp(['Elapsed time plotting waveforms= ' num2str(toc)]);
                    end
                    
                end
                % END OF "replot" CONDITIONAL

                % GET SPIKE VERT/HORIZ LIMITS
                if ~vert_spike_lims_approved | ~horiz_spike_lims_approved
                    if ~vert_spike_lims_approved
                        if vert_spike_lims(1)==0;
                            yl=get(gca,'ylim');
                            set(gca,'ylim',1.5*yl);
                        else
                            set(gca,'ylim',1.1*vert_spike_lims);
                        end
                        if 'y'==input('Add a vertical window for mvmt artifacts?  (y for yes, return for no)','s')
                            title('Click above and below waveforms to set window','fontweight','bold','fontsize',16)
                            [tmp,vert_spike_lims(1)]=ginput(1);plot(xl,vert_spike_lims(1)*[1 1],'k:','linew',2);
                            [tmp,vert_spike_lims(2)]=ginput(1);plot(xl,vert_spike_lims(2)*[1 1],'k:','linew',2);
                            vert_spike_lims=sort(vert_spike_lims);
                            cla;
                            redo_peaks = true;
                            replot = true;
                        else
                            vert_spike_lims_approved = true;
                        end

                        if redo_peaks
                            % The below is more efficient than re-generating spike_mat and
                            % id_peaks_save
                            %
                            %
                            %   id of rows of spike_mat that DONT contain a violation of the spike limits
                            %            max(M,[],2) is max across rows matrix M
                            id=find(min(spike_mat,[],2)>vert_spike_lims(1) & max(spike_mat,[],2)<vert_spike_lims(2));
                            spike_mat=spike_mat(id,:);
                            id_peaks_save=id_peaks_save(id);
                        end

                    else
                        if ~horiz_spike_lims_approved
                            yl_tmp=get(gca,'ylim');%plot(-.25*[1 1],yl_tmp,'r:');plot(1*[1 1],yl_tmp,'r:')
                            tmp=input('Set horizontal spike limits by hand?  ("y" for yes, return for no, "d" for default of  [-0.25 +1.00] msec)','s');
                            nsamp_wave_old=nsamp_wave;
                            if strcmp('y',tmp)
                                title('Click to left and right of waveforms to set limits','fontweight','bold','fontsize',16)
                                [horiz_spike_lims(1),tmp]=ginput(1);plot(horiz_spike_lims(1)*[1 1],yl_tmp,'k:','linew',2);
                                [horiz_spike_lims(2),tmp]=ginput(1);plot(horiz_spike_lims(2)*[1 1],yl_tmp,'k:','linew',2);
                                nsamp_wave=[floor(min(horiz_spike_lims*Fs./1000))*-1 ceil(max(horiz_spike_lims*Fs./1000))];
                                redo_peaks = true;
                                replot = true;
                            elseif strcmp('d',tmp)
                                horiz_spike_lims=[-0.25 1];
                                nsamp_wave=[floor(min(horiz_spike_lims*Fs./1000))*-1 ceil(max(horiz_spike_lims*Fs./1000))];
                                redo_peaks = true;
                            else
                                horiz_spike_lims_approved = true;
                            end
                        end
                        replot = true;
                        if redo_peaks
                            disp([newline, 'redo peaks', newline]);
                            % The below is more efficient than re-generating spike_mat and
                            % id_peaks_save
                            %
                            %
                            %   id of rows of spike_mat that DONT contain a violation of the spike limits
                            %            max(M,[],2) is max across rows matrix M
                            clip_start=nsamp_wave_old(1)-nsamp_wave(1);
                            clip_end=nsamp_wave_old(2)-nsamp_wave(2);
                            [r_sp,c_sp]=size(spike_mat);
                            
                            % clip out columns (time points) according to nsamp_wave
                            spike_mat_clipped=spike_mat(:,1+clip_start:c_sp-clip_end);

                            % check to make sure there is a point below threshold
                            % before peak.  this mimics the line
                            %                sum(potential_spike(1:nsamp_wave(1))*invert_sign<TH)
                            % in generate_spike_mat below
                            spike_mat_pre_peak=spike_mat_clipped(:,1:nsamp_wave(1));
                            
                            % TH is positive, and so are peaks
                            spike_mat_pre_peak=spike_mat_pre_peak*invert_sign;
                            
                            %            min(M,[],2) is min across rows matrix M
                            % id is rows of spike_mat_pre_peak that contain at least
                            % one value less than TH
                            id=find(min(spike_mat_pre_peak,[],2)<TH);
                            spike_mat=spike_mat_clipped(id,:);
                            id_peaks_save=id_peaks_save(id);

                        end
                    end                    
                else
                    finished = true;
                    
                    S.spike_mat = spike_mat;
                    S.spike_ids = id_peaks_save;
                    
                    set(gca,'xlim',[min(xdat) max(xdat)]);
                end
                disp([newline, 'End While Loop', newline]);
            end
        end
        
        % ===== ===== ===== ===== =====
        % Optimize clusters
        function optimize_clusters(S, top_right_plot_code, threshold_strategy)
            % top_right_plot_code=1  PCA based on amplitude alone
            % top_right_plot_code=2 determine optimal threshold
            if nargin == 1
                top_right_plot_code = 2;
                threshold_strategy = 1;
            elseif nargin == 2
                threshold_strategy = 1;
            end
            
            % CASES FOR top_right_plot_code = 2
            %    threshold_strategy = 1
            %         Finds largest (furthest from zero) theshold such that the false
            %         positive rate is less than 1%.  That is, the threshold is set
            %         such that if a spike exceeds the threshold, there is a 99% chance
            %         that it came from the cluster with the largets mean amplitude.
            %         As a consequence, if there is a lot of amplitude overlap between
            %         clusters, a large number of spikes will be excluded.  However,
            %         you can have high confidence (99%) that the included spikes were
            %         a part of the biggest-spike cluster
            %
            %        Then, finds mean - 2.5 S.D. of amplitudes for largest cluster.
            %         TH_empirical will the the LARGER (furthest from zero) of
            %         these two measures.
            %
            %     threshold_strategy = 2
            %          Finds the threshold such that the false positive and false
            %          negative rates are equal (or as close to equal as possible).
            %          That is, the threshold is set such that the rate at which spikes
            %          from non-biggest cluster are beyond threshold (fales
            %          positives)
            %          is equal to the rate at which spikes from the biggest-spike
            %          cluster are less than the threshold (false negatives).  Compared
            %          the the other strategy, this will miss fewer spikes, but will
            %          include a much greater number of noise waveforms.
            %
            %  Note that for very large units, these two strategies will produce nearly
            %  identical answers.
            
            if S.verbose > 1
                show_timing = true;
            else
                show_timing = false;
            end
            
            n_clusters = S.n_clusters;
            id_peaks_save = S.spike_ids';
            IDX = S.kMeans.idx';
            C = S.kMeans.C;
            
            dat = S.data.amplifier_data_filt';
            Fs = S.Fs;
            
            TH_original = S.TH_original;
            
            SCORE = S.PCA.score;
            spike_mat = S.spike_mat;
            invert_sign = S.inv_sign;
            nsamp_wave = S.nsamp_wave;
            
            col_vec = S.col_vec;
            col_mat = S.col_mat;
            
            subplot(2,9,4:7);
            hold on;

            if show_timing
                tic;
                disp('Starting timer...');
            end
            
            samp_dat = S.get_second_samp();
            id_for_samp = samp_dat.ixs;
            t_lims=id_for_samp([1 end])/Fs;
            
            for x=1:n_clusters
                id_tmp=find(IDX==x);
                spike_ids_by_pca{x}=id_peaks_save(id_tmp);
                spike_times{x}=spike_ids_by_pca{x}/Fs;              % all spike times - BASED ON PEAKS, NOT CROSSINGS
                all_mags{x}=dat(spike_ids_by_pca{x});
                mean_peaks_later(x)=mean(dat(spike_ids_by_pca{x}));
                
                spike_ids_tmp=spike_times{x}(find(spike_times{x}>t_lims(1) & spike_times{x}<t_lims(2)));
                plot(spike_ids_tmp-t_lims(1)+1/Fs,dat(round(spike_ids_tmp*Fs)),'o','color',col_vec(x),'markerfacecolor',col_vec(x));
            end
            
            set(gca,'xlim',[0 length(id_for_samp)/Fs])

            if show_timing
                disp(['Elapsed time plotting dots on example waveform= ' num2str(toc)]);
            end
            
            plot([1:length(id_for_samp)]/Fs,dat(id_for_samp),'k')
            %    plot([1:Fs]/Fs,dat(id_for_samp),'k')
            set(gca,'ylim',[min(min(spike_mat)) max(max(spike_mat))])
            yl=get(gca,'ylim');
                
            % PCA BASED ON AMPLITUDE ALONE
            if top_right_plot_code==1
                %if plot_flag;subplot(2,10,19:20);cla;hold on;end
                
                subplot(4,9,8:9);
                hold on;
                
                % invert spike mat if downward spike - now all are upward spikes
                spike_mat_tmp=spike_mat*invert_sign;
                
                %    spike_mat_tmp=spike_mat_for_plotting_and_amp*invert_sign;
                % column of amplitudes of minima
                peak_vec=spike_mat_tmp(:,nsamp_wave(1)+1);
                for x=1:length(peak_vec)
                    first_der=diff(spike_mat_tmp(x,:));
                    id=min(intersect(find(first_der>=0),nsamp_wave(1)+1:sum(nsamp_wave)+1));
                    if isempty(id)
                        id=sum(nsamp_wave)+1;
                        disp('Ending spike at end of waveform')
                    end
                    trough_vec(x)=spike_mat_tmp(x,id);
                end
                
                spike_mag_vec=peak_vec-trough_vec';
                
                % demo to show this works
                %     figure(100);clf;hold on;plot(spike_mat_tmp(x,:),'r');id,[peak_vec(x)  trough_vec(x)],return
                [IDX_mags,C_mags] = kmeans(spike_mag_vec, n_clusters);
                [tmp,spike_id_by_size]=sort(C_mags);

                %%%%%%%%% reorder clusters by spike size
                % so noise_spikes are in cluster 1, and spikes are in clusters 2-n
                IDX_mags_old=IDX_mags;
                C_mags_old=C_mags;
                clear IDX_mags C_mags
                for x=1:n_clusters
                    IDX_mags(find(IDX_mags_old==spike_id_by_size(x)))=x;
                    C_mags(x,:)=C_mags_old(spike_id_by_size(x),:);
                end
                for x=1:n_clusters
                    id_tmp=find(IDX_mags==x);
                    n_hist_bins=round(30/n_clusters);
                    [height,ctrs]=hist(spike_mag_vec(id_tmp),n_hist_bins);
                    
                    a=bar(ctrs,height);
                    set(a,'facecolor',col_vec(x))

                    mean_vec(x)=mean(spike_mag_vec(id_tmp));
                    std_vec(x)=std(spike_mag_vec(id_tmp));
                end

                subplot(4,9,17:18);
                hold on;
                
                for x=1:n_clusters
                    n_gaussfit_bins=round(std_vec(x))*5;
                    g=gaussian_lf(std_vec(x),n_gaussfit_bins);
                    plot([1:length(g)]-length(g)/2+mean_vec(x),g,col_vec(x),'linew',2)
                end
                n_points_1D=1000;
                sim_cluster_id_1D=[];
                out_all_1D=[];
                for x=1:n_clusters
                    % compute and plot covariance matrices for each computed cluster
                    id_tmp=find(IDX_mags==x);

                    % generate n_points of data generated from 1D gaussian with
                    % the computed variance and mean
                    normnoise=randn(n_points_1D,1)*std_vec(x)+mean_vec(x);

                    % record cluster number from which each point was generated
                    sim_cluster_id_1D=[sim_cluster_id_1D; x*ones(n_points_1D,1)];

                    % combine synthetic clusters to be sorted according to KMEANS
                    out_all_1D=[ out_all_1D; normnoise];
                end
                [IDX_mags_simulated,C_mags_simulated] = kmeans(out_all_1D,n_clusters);
                [tmp,spike_id_by_size]=sort(C_mags_simulated);
                %%%%%%%%% reorder clusters by spike size
                % so noise_spikes are in cluster 1, and spikes are in clusters 2-n
                IDX_mags_old=IDX_mags_simulated;C_mags_old=C_mags_simulated;
                clear IDX_mags C_mags
                for x=1:n_clusters
                    IDX_mags_simulated(find(IDX_mags_old==spike_id_by_size(x)))=x;
                    C_mags_simulated(x,:)=C_mags_old(spike_id_by_size(x),:);
                end

                disagreements_1D=find(IDX_mags_simulated ~= sim_cluster_id_1D);

                title(['Total: ' num2str(length(disagreements_1D)) ' errors out of ' num2str(length(IDX_mags_simulated))])
                xl=get(gca,'xlim');yl=get(gca,'ylim');y_range=diff(yl);
                new_yl=[yl(1) yl(2)+1.15*y_range];
                set(gca,'ylim',new_yl,'xlim',xl);
                text(xl(1)+.2*diff(xl),new_yl(2)-.2*y_range,['Of ' num2str(n_points_1D*n_clusters) ' simulated points: '],'color','k','fontweight','bold')

                error_str=[];
                for x=1:n_clusters
                    IDX_tmp=IDX_mags_simulated;
                    sim_cluster_id_tmp=sim_cluster_id_1D;
                    IDX_tmp(find(IDX_tmp~=x))=9999;
                    sim_cluster_id_tmp(find(sim_cluster_id_tmp~=x))=9999;
                    n_disagreements_by_cluster(x)=length(find(IDX_tmp ~= sim_cluster_id_tmp));
                    pct_error(x)=n_disagreements_by_cluster(x)/(n_points_1D*n_clusters);
                    pct_error_str{x}=num2str(round(10000*pct_error(x))/100);
                    %    disp([num2str(n_disagreements_by_cluster(x)) ' errors (' pct_error_str{x} '%) between cluster ' num2str(x) ' (color = ' col_vec(x) ') and all others' ])
                    
                    text(xl(1)+.2*diff(xl),new_yl(2)-.2*y_range*(x+1),[num2str(n_disagreements_by_cluster(x)) ' errors (' pct_error_str{x} '%) between cluster ' num2str(x) ' & others' ],'color',col_vec(x),'fontweight','bold')
                    
                    error_str=[error_str pct_error_str{x} '%          '];
                end

            % FIND THE BEST THRESHOLD
            elseif top_right_plot_code==2

                subplot(2,9,8:9);
                cla;
                hold on;

                if show_timing
                    tic;
                    disp('Starting timer...');
                end
                
                concat_mags=[];
                concat_N=[];
                concat_spiketimes=[];
                
                % plot histogram of amplitudes and save magnitudes
                for x=1:n_clusters
                    [N_occurences{x},bin_cov_centerss{x}]=hist(all_mags{x},20);
                    
                    plot(N_occurences{x},bin_cov_centerss{x},[col_vec(x) '-o'],'linew',2,'markersize',4);
                    
                    % concat_mags=[concat_mags all_mags{x}];
                    concat_mags=[concat_mags all_mags{x}' ];
                    concat_N=[concat_N N_occurences{x} ];
                end
                % all spiketimes, in order of time
                [concat_spiketimes,tmp]=sort(cell2mat(spike_times));
                % all magnitudes, in order of time
                mags_resorted=concat_mags(tmp);

                % this will put one potential threshold between each spike magnitude.
                th_sweep_vec=[unique(floor(unique(concat_mags))) ceil(max(concat_mags))];

                if TH_original>0
                    th_sweep_vec=fliplr(th_sweep_vec);
                    disp('Flipping th_sweep_vec');
                end

                % combing spike amplitudes of all clusters except for the biggest
                all_noise_mags=[];
                for x=1:n_clusters-1
                    all_noise_mags=[all_noise_mags; all_mags{x}];
                end
                if max(abs(all_noise_mags))>max(abs(all_mags{end}))
                    S.warning_str{end+1}='Biggest noise spike is bigger than biggest signal spike';
                    disp(' ');
                    disp(S.warning_str{end});
                    disp(' ');
                end

                % sweep through the potential threshold values and record number of
                % each type of error
                for x=1:length(th_sweep_vec)
                    n_spikes_as_noise_vec(x)=length(find(abs(all_mags{end})<abs(th_sweep_vec(x))));%last entry is biggest spike
                    n_spikes_as_spikes_vec(x)=length(find(abs(all_mags{end})>abs(th_sweep_vec(x))));% last entry is biggest spike
                    n_noise_as_spikes_vec(x)=length(find(abs(all_noise_mags)>abs(th_sweep_vec(x))));   %noise
                    total_n_beyond_threshold(x)=length(find(abs(concat_mags)>abs(th_sweep_vec(x))));
                end

                if show_timing
                    disp(['Elapsed time running amp for-loop= ' num2str(toc)]);
                    tic;
                    disp('Starting timer...');
                end

                if threshold_strategy==1
                    warning off % the next line generates a 'divide by zero' warning
                    pct_of_positives_that_are_true_positives=n_spikes_as_spikes_vec./total_n_beyond_threshold;
                    warning on
                    id_tmp=max(find(pct_of_positives_that_are_true_positives>.95));
                    if isempty(id_tmp)
                        [tmp,id_tmp]=max(pct_of_positives_that_are_true_positives);
                        S.warning_str{end+1}=['pct_of_positives_that_are_true_positives never gets to criterion - using peak = ' num2str(tmp)];
                        disp(S.warning_str{end})
                    end
                    pct_of_positives_that_are_true_positives_at_thresh=pct_of_positives_that_are_true_positives(id_tmp);
                    pct_of_spikes_that_are_above_thresh=n_spikes_as_spikes_vec(id_tmp)/length(all_mags{end});
                    if pct_of_spikes_that_are_above_thresh<.1
                        disp('Redoing...')
                        id_tmp=min(find(n_spikes_as_spikes_vec/length(all_mags{end})>.25));
                        pct_of_positives_that_are_true_positives_at_thresh=pct_of_positives_that_are_true_positives(id_tmp);
                        pct_of_spikes_that_are_above_thresh=n_spikes_as_spikes_vec(id_tmp)/length(all_mags{end});
                        S.warning_str{end+1}='Less than 10% of main-cluster spikes beyond threshold, resetting thresh so that 25% are within';
                    end
                    TH_empirical_by_type_1_error=th_sweep_vec(id_tmp);
                    % % figure showing how this works
                    %         ww=gcf;figure(101);clf;subplot(1,2,1);hold on;plot(th_sweep_vec,total_n_beyond_threshold,'b');
                    %         plot(th_sweep_vec,n_spikes_as_spikes_vec,'g');xlabel('blue is total N beyone thresh, green is N spikes as spikes')
                    %         subplot(1,2,2);hold on;plot(th_sweep_vec,pct_of_positives_that_are_true_positives)
                    %         yl=get(gca,'ylim');plot(TH_empirical*[1 1],yl,'k--');figure(ww)

                    % NEW - if no overlap between spike mags in biggest-mag cluster and
                    % others, set empirical threshold at min of biggest-mag cluster
                    if min(abs(all_mags{end}))>max(abs(all_mags{end-1}))
                        id_tmp=find(    abs(all_mags{end})  ==    min(abs(all_mags{end})));
                        TH_empirical_by_STD=all_mags{end}(id_tmp);
                        disp('NOTE:  No overlap between biggest and next-biggest cluster - setting SD-based threshold to min of big cluster')
                    else
                        TH_empirical_by_STD=mean(all_mags{end})-2.5*sign(mean(all_mags{end}))*std(all_mags{end});
                    end

                    if abs(TH_empirical_by_STD)>abs(TH_empirical_by_type_1_error)
                        TH_empirical=TH_empirical_by_STD;
                        disp('Using mean - 2.5 SD of amplitudes as threshold')
                    else
                        TH_empirical=TH_empirical_by_type_1_error;disp('Using type 1 error criterion as threshold')

                          %  TH_empirical=-5000;disp('HACK SETTING TH_empirical to -5000')
                    end
                elseif threshold_strategy==2
                    n_in_spike_cluster=length(all_mags{end});% last entry is biggest spike
                    n_in_noise_cluster=length(all_noise_mags);

                    pct_1_vec=n_noise_as_spikes_vec./(n_spikes_as_spikes_vec+n_noise_as_spikes_vec);
                    pct_2_vec=n_spikes_as_noise_vec./n_in_spike_cluster;

                    abs_diff_vec=abs(pct_1_vec-pct_2_vec);

                    if length(abs_diff_vec)==1;whos;end

                    TH_empirical_id=find(abs_diff_vec==min(abs_diff_vec));
                    TH_empirical=th_sweep_vec(TH_empirical_id);
                    if length(TH_empirical)>1
                        for x=1:length(TH_empirical)
                            whos TH_empirical th_sweep_vec
                            TH_id(x)=find(TH_empirical(x)==th_sweep_vec);
                        end
                        if prod(diff(TH_id))==1
                            TH_empirical=mean(TH_empirical);
                            TH_empirical_id=round(mean(TH_empirical_id));
                        else
                            S.warning_str{end+1}='Two discontinuous possible locations for optimal threshold';
                            disp(S.warning_str{end});
                            TH_empirical=0;
                        end
                    end
                    n_noise_as_spikes=n_noise_as_spikes_vec(TH_empirical_id);
                    n_spikes_as_noise=n_spikes_as_noise_vec(TH_empirical_id);
                    n_spikes_as_spikes=n_spikes_as_spikes_vec(TH_empirical_id);
                end

                title('Spike amplitude and recommended thresh');
                ylabel('Max amp of spike');
                xl=[min(concat_N) max(concat_N)];
                set(gca,'xlim',xl);
                yl=get(gca,'ylim');
                plot(xl,TH_empirical*[1 1],'m-o','linew',3,'markerfacecolor','m')
                if threshold_strategy==1
                    pct_txt=num2str(round(pct_of_positives_that_are_true_positives_at_thresh*1000)/10);
                    pct_txt_2=num2str(round(pct_of_spikes_that_are_above_thresh*1000)/10);
                    xlab_txt{1}='N';
                    xlab_txt{2}=['              Cyan line - Mean - 2.5 S.D. error'];
                    xlab_txt{3}=[ '          Black line - Type I error criterion'];
                    xlab_txt{4}=[pct_txt '% of thresholded spikes are from primary cluster'];
                    xlab_txt{5}=[pct_txt_2 '% of spikes from primary cluster are past thresh'];
                    xlabel(xlab_txt);

                    plot(xl,TH_empirical_by_STD*[1 1],'c','linew',1)
                    plot(xl,TH_empirical_by_type_1_error*[1 1],'k','linew',1)


                elseif threshold_strategy==2
                    xlabel('N');
                    pct_error_noise_as_spikes=round(10000*n_noise_as_spikes/(n_spikes_as_spikes+n_noise_as_spikes))/100;
                    pct_error_spikes_as_noise=round(10000*n_spikes_as_noise/n_in_spike_cluster)/100;
                    text(xl(1)+.1*diff(xl),yl(2)-.1*diff(yl),[num2str(pct_error_noise_as_spikes) '% thresholded spikes actually from other cluster(s)' ])
                    text(xl(1)+.1*diff(xl),yl(2)-.2*diff(yl),[num2str(pct_error_spikes_as_noise) '% from biggest spike cluster excluded by threshold' ])
                end
                subplot(2,9,4:7);xl=get(gca,'xlim');
                plot(xl,TH_empirical*[1 1],'m--','linew',2)

                ids_all_past_TH_empirical=find(abs(mags_resorted)>abs(TH_empirical));
                spiketimes_from_recommended_TH=concat_spiketimes(ids_all_past_TH_empirical);

                if show_timing
                    disp(['Elapsed time computing threshold = ' num2str(toc)]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            %if plot_flag;subplot(2,2,3);hold on;end
            subplot(2,10,11:14);
            cla;
            hold on;

            if show_timing;
                tic;
                disp('Starting timer...');
            end

            only_plot_subset=0;
            %pct_to_plot_PCA=.2;  % percent of each cluster to plot
            n_to_plot_PCA=1000;  % number from each cluster to plot

            for x=1:n_clusters
                id_tmp=find(IDX==x);
                n_points_in_cluster(x)=length(id_tmp);
                cov_centers{x}=mean(SCORE(id_tmp,1:2));
                cov_matrices{x}=cov(SCORE(id_tmp,1),SCORE(id_tmp,2));
                
                ellipsoid=plot_cov(cov_centers{x},cov_matrices{x},col_vec(x),.0455);
                if only_plot_subset
                    %             disp(['ONLY PLOTTING ' num2str(pct_to_plot_PCA*100) ' % OF PCA POINTS IN LOWER LEFT'])
                    %             id_tmp_part=id_tmp(1:round(length(id_tmp)*pct_to_plot_PCA));
                    tstr{1}=['ONLY PLOTTING ' num2str(n_to_plot_PCA) ' OF PCA POINTS FROM EACH CLUSTER'];disp([tstr{1} ' IN LOWER LEFT'])
                    id_tmp_part=save_example_ids{x};
                    plot(SCORE(id_tmp_part,1),SCORE(id_tmp_part,2),[col_vec(x) 'o' ],'markersize',2,'markerfacecolor',col_vec(x))
                else
                    plot(SCORE(id_tmp,1),SCORE(id_tmp,2),[col_vec(x) 'o' ],'markersize',2,'markerfacecolor',col_vec(x))
                end
                plot(C(x,1),C(x,2),[col_vec(x) 'x'] ,'linew',2)

                if 0% threshold_strategy==1% circles the points beyond thresh
                    % these are the ids of spikes from rec thresh
                    [tmp,tmp,id_tmp2]=intersect(Data.spiketimes_from_recommended_TH,spike_times{x});
                    plot(SCORE(id_tmp(id_tmp2),1),SCORE(id_tmp(id_tmp2),2),['mo' ])
                end

            end
            tstr{2}='Ellipses are 95.4% conf intervals (2 S.D.)';
            
            title(tstr);
            xl=get(gca,'xlim');
            yl=get(gca,'xlim');
            set(gca,'xlim',xl,'ylim',yl);
            
            if show_timing,disp(['Elapsed time plotting PCA clusters= ' num2str(toc)]);end

            %if plot_flag;subplot(2,2,4);hold on;end
            subplot(2,10,15:18);
            cla;
            hold on;
            
            if show_timing;
                tic;
                disp('Starting timer...');
            end

            out_all=[];
            sim_cluster_id=[];
            % will simulate a total of 10000 points, with number of points in each
            % cluster proportional to the size of each cluster
            n_points_vec=round(n_points_in_cluster/sum(n_points_in_cluster)*10000);

            for x=1:n_clusters
                % compute and plot covariance matrices for each computed cluster
                id_tmp=find(IDX==x);

                % generate n_points of data generated from 2D gaussian with
                % the computed covariance and center
                n_points=n_points_vec(x);
                normnoise=randn(n_points,2);

                out_newcluster{x}=normnoise*sqrtm(cov_matrices{x});
                out_newcluster{x}(:,1)=out_newcluster{x}(:,1)+cov_centers{x}(1);
                out_newcluster{x}(:,2)=out_newcluster{x}(:,2)+cov_centers{x}(2);

                % plot points generated with these clusters
                
                plot_cov(cov_centers{x},cov_matrices{x},col_vec(x),.0455);
                plot(out_newcluster{x}(:,1),out_newcluster{x}(:,2),[col_vec(x)  'o'],'markersize',2,'markerfacecolor',col_vec(x))


                % record cluster number from which each point was generated
                sim_cluster_id=[sim_cluster_id; x*ones(n_points,1)];

                % combine synthetic clusters to be sorted according to KMEANS
                out_all=[ out_all; out_newcluster{x}];
            end

            L_mat=length(out_all);
            n_blocks=ceil(L_mat/1000);  % do this 1000 points at a time
            dist_mat_by_data_centers=[];
            % divide this up into three blocks to avoid running out of memory
            for x=1:n_blocks
                ind=1+(x-1)*1000:min([x*1000 L_mat]);
                MAT_FOR_DISTANCES=[C;out_all(ind,:)];
                % a matrix of distances between all these points
                % NOTE:  function pdist isnt really the right one, since it calculated
                % distance between ALL points, not just the C matrix.  this is sort of
                % wasteful in terms of memory
                DIST_MAT_tmp=squareform(pdist(MAT_FOR_DISTANCES));
                % clip off distance columns corresponding to  C matrix
                DIST_MAT_tmp=DIST_MAT_tmp(:,n_clusters+1:end);
                % Make an n_clusters x n_simulated_points matrix of distances
                for j=1:n_clusters
                    dist_mat_by_data_centers(j,ind)=DIST_MAT_tmp(j,:);
                end
            end
            clear DIST_MAT_tmp MAT_FOR_DISTANCES

            % find minima of these (IDX_simulated is id of the row of
            % dist_mat_by_data_centers that has minimum distance
            [tmp,IDX_simulated]=min(dist_mat_by_data_centers);
            IDX_simulated=IDX_simulated';

            % cases where the classification of the simulated data (IDX_simulated)
            % doesnt match the classification that generated it (sim_cluster_id)
            disagreements=find(IDX_simulated ~= sim_cluster_id);
            
            for x=1:n_clusters
                id_tmp=intersect(find(IDX_simulated==x),disagreements);
                plot(out_all(id_tmp,1),out_all(id_tmp,2),[col_vec(x) 'o' ])
            end
            set(gca,'xlim',xl,'ylim',yl);
            y_range=diff(yl);
            text(xl(1)+.2*diff(xl),yl(2),['Of ' num2str(sum(n_points_vec)) ' simulated points: '],'color','k','fontweight','bold')


            error_str=[];
            for x=1:n_clusters

                classed_as_x_AND_from_x=length(find(IDX_simulated==x & sim_cluster_id==x));
                NOT_classed_as_x_AND_from_x=length(find(IDX_simulated~=x & sim_cluster_id==x));
                classed_as_x_AND_from_NOT_x=length(find(IDX_simulated==x & sim_cluster_id~=x));
                classed_as_x=length(find(IDX_simulated==x));
                from_x=length(find(sim_cluster_id==x));
                NOT_from_x=length(find(sim_cluster_id~=x));

                FALSE_POSITIVE_n(x)=classed_as_x_AND_from_NOT_x;
                FALSE_NEGATIVE_n(x)=NOT_classed_as_x_AND_from_x;
                TOTAL_ERROR_n(x)=NOT_classed_as_x_AND_from_x+classed_as_x_AND_from_NOT_x;

                FALSE_POSITIVE_RATE(x)=classed_as_x_AND_from_NOT_x/NOT_from_x;
                FALSE_NEGATIVE_RATE(x)=NOT_classed_as_x_AND_from_x/from_x;
                pct_error(x)=(NOT_classed_as_x_AND_from_x+classed_as_x_AND_from_NOT_x)/(from_x + NOT_from_x);

                pct_error_str{x}=num2str(round(10000*pct_error(x))/100);
                
                text(xl(1)+.2*diff(xl),.9*yl(2)-.05*y_range*(x-1),[num2str(TOTAL_ERROR_n(x)) ' errors (' pct_error_str{x} '%) between cluster ' num2str(x) ' & others' ],'color',col_vec(x),'fontweight','bold')
                
                error_str=[error_str pct_error_str{x} '%          '];
            end
            disp([error_str '     '  num2str(TH_empirical) '    ' S.experiment_name])
            if show_timing,disp(['Elapsed time simulating distributions= ' num2str(toc)]);end

            subplot(2,10,19:20);
            cla;
            hold on;
            
            for x=1:n_clusters+1
                id_tmp=find(IDX==x);
                %     spike_ids_by_pca{x}=id_peaks_save(id_tmp);
                %     spike_times{x}=spike_ids_by_pca{x}/32000;
                ISI_vec=0:.0005:.1;ID_1ms=find(ISI_vec<=.001);
                col_tmp='rgbmky';
                max_y=0;
                if x<=n_clusters
                    y=hist(diff(spike_times{x}),ISI_vec);
                else
                    y=hist(diff(spiketimes_from_recommended_TH),ISI_vec);
                end
                pct_Leq_1ms(x)=sum(y(ID_1ms))/sum(y);
                
                if x<=n_clusters
                    plot(ISI_vec(1:end-1),y(1:end-1)/sum(y),[col_tmp(x)],'linew',2)
                else
                    plot(ISI_vec(1:end-1),y(1:end-1)/sum(y),[col_tmp(x-1)],'linew',1)
                end
                yl(x,:)=get(gca,'ylim');
                xl=[0 .1];
                set(gca,'xlim',xl,'ylim',[0 max(yl(:,2))])
                plot([.001 .001],[0 1],'k--')

            end
            
            yl=min(yl,[],1);
            for x=1:n_clusters+1
                tmp=num2str(round(10000*pct_Leq_1ms(x))/100);
                if x<=n_clusters
                    text(xl(1)+.1*diff(xl),yl(2)-.075*(x)*diff(yl),['% ISIs <= 1 ms - ' num2str(tmp) '%'],'fontweight','bold','color',col_tmp(x))
                else
                    text(xl(1)+.1*diff(xl),yl(2)-.075*(x)*diff(yl),['% ISIs <= 1 ms - ' num2str(tmp) '% from thresh'],'color',col_tmp(x-1))
                end
            end
        end
        
        function do_saucy(S)
            % ----->> NEED TO SET UP DEFAULT ARGUMENTS FOR FUNCTIONS
            filter_data(S);
            set_user_threshold(S);
            do_clustering(S);
            optimize_clusters(S);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NEED TO REASSIGN OUTPUT DATA STRUCTURE OR WRITE FUNCTIONS THAT ALLOW
%%%%% ACCESS TO THESE VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% clear spike_mat_tmp
% Data.filename=RHD_name;
% Data.n_units=n_clusters-1;  % number of non-noise waveforms
% Data.vert_spike_lims=vert_spike_lims;
% Data.nsamp_wave=nsamp_wave;
% Data.run_bandpass=run_bandpass;
% Data.bandpass_freqs=bandpass_freqs;
% Data.TH=TH_set_by_user;
% % for x=1:n_clusters-1
% %     Data.spiketimes{x}=spike_ids_by_pca{x+1}/Fs;  % 2 is spike
% % end
% for x=1:n_clusters
%     Data.spiketimes{x}=spike_ids_by_pca{x}/Fs;  % 2 is spike
% end
% Data.recommended_TH=TH_empirical;
% Data.spiketimes_from_recommended_TH=spiketimes_from_recommended_TH;
% waveform_temp=[];
% for z=1:n_clusters
%     id_tmp=find(IDX==z);
%     spike_mat_tmp=spike_mat(id_tmp,:);
%     Data.mean_waveform{z}=mean(spike_mat_tmp);
%     Data.std_waveform{z}=std(spike_mat_tmp);
%     waveform_temp=[waveform_temp; spike_mat_tmp];
% end
% Data.waveforms=waveform_temp;
% %Data.components=component;
% Data.components=SCORE;
% Data.COEFF_matrix=COEFF;
% Data.cov_matrices=cov_matrices;
% Data.cov_centers=cov_centers;
% Data.pct_error=pct_error;
% Data.total_n_samples=length(dat);
% if exist('warning_str')
%     Data.warning_str=warning_str;
% else
%     Data.warning_str=[];
% end
% Data.Fs=Fs;
% Data.inverted_waveform=inverted_waveform;
% Data.pct_ISIs_leq_1ms=pct_Leq_1ms;
% if batchmode_condition==0
%     Data.LATENT = LATENT ;
% end
% if length(ch_num)==1
%     eval(sprintf('save %s.neuralnot_CH%s.mat Data',RHD_name,num2str(ch_num)))
% elseif length(ch_num)==2
%     eval(sprintf('save %s.neuralnot_CH%s_minus_CH%s.mat Data',RHD_name,num2str(ch_num(1)),num2str(ch_num(2))))
% end