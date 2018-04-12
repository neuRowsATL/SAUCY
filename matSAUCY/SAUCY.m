classdef SAUCY < handle
    % SAUCY Object class to manage analysis of single file
    %   SAUCY allows analysis for a single file and can be used for
    %   threshold detection, PCA classification, and optimization of false
    %   positives and false negatives.
    
    properties
        verbose = 0 % level of verbose output
        
        experiment_name = '' % base name of experiment
        
        data_path = '' % path to data files
        fname_intan = '' % full file name of intan file
        fname_cbin = '' % full file name of cbin file
        fname_neuralnot = '' % full file name of neural_not file
        
        % chan          Channel number(s) for analysis
        %   as INT      Single channel will be used
        %   as TUPLE    Channel subtraction will be used [a b] = a - b
        chan = []
        
        n_clusters = 2 % number of clusters to detect
        
        Fs = 30000 % sampling rate of recording
        F_low = 350 % low cut off
        F_high = 7000 % hi cut off
        filter_type = 'hanningfir' % type of filter to apply
        
        raw_data = {} % struct for raw data 
        data = {} % struct for results
        
        inv_sign % indicates whether waveform is inverted
        inv_waveform % indicates whether waveform has been inverted or not
        
        threshold = [] % Spike detection threshold
        use_crossing % indicates whether to align to peaks or threshold
        nsamp_wave = [30 90]
        
        spike_mat % matrix of spike waveforms
        id_peaks_save % id of peaks... ?
    end
    
    methods
        % ----- ----- ----- ----- -----
        % Initialize SAUCY object
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
        % Load data into SAUCY object
        function load_data(S, filename)
            if strfind(filename, '.rhd') > 0
                % Load data from .rhd file
                [t_amplifier, t_board_adc, amplifier_data, board_adc_data, frequency_parameters] = ...
                   read_Intan_RHD2000_nongui_saucy(filename);
                S.raw_data.t_amplifier = t_amplifier;
                S.raw_data.t_board_adc = t_board_adc;
                S.raw_data.amplifier_data = amplifier_data;
                S.raw_data.board_adc_data = board_adc_data;
                
                % Set default data for analysis
                S.data.amplifier_data_filt = S.raw_data.amplifier_data;
               
                S.Fs = frequency_parameters.amplifier_sample_rate;
               
                S.fname_intan = filename;
               
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
        function filter_data(S, f_specs)
            if nargin == 1
                disp('!! USING DEFAULT BANDPASS LIMITS [350 7000] !!');
            elseif nargin == 2
                if length(f_specs) == 2
                    S.F_low = f_specs(1);
                    S.F_high = f_specs(2);
                end
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
        % Set threshold for spike detection
        function set_threshold(S, data_source, std_th)
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
                    dat = S.raw_data.amplifier_data;
                else
                    disp('ERROR: Check that data was loaded into SAUCY.');
                end
            else
                disp('ERROR: Invalid data source format.');
            end
                
            Fs = S.Fs;
            
            % ----->> MIGHT NEED TO RECALIBRATE THIS IF SPIKES ARE SPARSE
            % IN FILE
            
            % Selects 1 second of data from the middle of the file
            % If data is shorter than 1 sec, uses all of data
            mid_dat=round(length(dat)/2)-Fs/2;
            if length(dat)/Fs > 1
                id_for_samp=mid_dat:(mid_dat+Fs-1);    
            else
                id_for_samp=1:length(dat);    
            end

            % Plot waveform sample for user selection of threshold
            figure();
            plot([1:length(id_for_samp)]/Fs, dat(id_for_samp), 'k');
            set(gca,'xlim',[0 length(id_for_samp)/Fs]);

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
                % ----->> CONSIDER USING ADJUSTED STD based on Rossant 2016?
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
            
            TH_original = TH;
            S.threshold = TH;
            
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

            % ----->> POSSIBLE BUG NEEDS TO BE FIXED
            % dat and TH are consistent when reloaded...  but looks like vertical lims
            % are bad
        end
        
        % ----- ----- ----- ----- -----
        % Extract spike waveforms aligned by peak or threshold crossing
        function set_spike_mat(S, alignment)
            % Identify spike peaks (copied from neural_and_song.m)
            % alignment             Specify how to align spikes
            %       peaks           Align spikes by peak
            %       threshold       Align spikes by threshold crossing
            
            % If no alignment mode is specified, ALIGN SPIKES TO PEAKS
            if nargin == 1
                S.use_crossing = false;
            elseif nargin == 2
                if strcmp(alignment, 'peaks') % use peaks to align spikes
                    S.use_crossing = false;
                elseif strcmp(alignment, 'threshold') % use threshold crossing to align spikes
                    S.use_crossing = true;
                end
            end
            
            use_crossing = S.use_crossing;
            
            if S.verbose > 1
                show_timing = true;
            else
                show_timing = false;
            end

            % Identify peaks of spikes
            % -----
            id = find(dat > TH); % id of samples that are above threshold
            sdx = diff(sign(diff(dat))); % this will be -2 one point before upward-pointing peak, +2 one point before downward-pointing peak, zero else
            id2 = find(sdx < 0) + 1; % id2 is ids of UPWARD-POINTING PEAKS
            
            % first deriv not equal to zero - use this to prevent multiple spikes being assigned to truncated spikes
            % so to be identified as a spike,  sample must be (1) above threshold (2) an upward-pointing peak and
            % (3)  changing in value
            first_der=diff(dat);
            not_constant=find(first_der);
            id_peaks=intersect(id,id2);     
            id_peaks=intersect(id_peaks,not_constant);
            
            % Identify threshold crossings
            % -----
            id_crossings=find(diff(dat>TH)==1)+1;           % id of samples that first to cross threshold

            % # of peaks will be slightly higher than # crossings.  Below is a vector
            % of id_peaks that fall immediately after each threshold crossing.  This is
            % for plotting only
            if use_crossing
                for x=1:length(id_crossings)
                    id_peaks_after_crossings(x,1) = min(id_peaks(find(id_peaks>id_crossings(x))));
                end
            end

            %  Un-inverting waveform
            if S.inv_sign == -1
                dat = -S.data.amplifier_data_filt;
            else
                dat = S.data.amplifier_data_filt;
            end

            TH = S.threshold;
            invert_sign = S.inv_sign;
            n_clusters = S.n_clusters;
            
            % Extract spike waveforms and establish vertical and horizontal
            % limits
            % -----
            col_vec='rgbc';
            col_mat=[1 0 0;0 1 0;0 0 1;1 1 0];

            finished = false;
            
            vert_spike_lims_approved = false;
            horiz_spike_lims_approved = false;
            
            redo_peaks = true;
            replot = true;
            first_plotting = true;
            
            % [# samples before peak        # after peak] to subject to PCA (these are
            % defaults, user will be able to adjust these below)
            nsamp_wave = S.nsamp_wave;
            vert_spike_lims=[0 0]; % use vertical limits to cut off mvmt artifact (default is not to use these)
                    
            while ~finished
                if redo_peaks
                    id_peaks = id_peaks(find(id_peaks>(nsamp_wave(1)+1)  & id_peaks<(length(dat)-nsamp_wave(2))));
                    
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end

                    if ~exist('spike_mat')  % a new trick - only run generate_spike_mat once
                        if ~use_crossing
                            % vert_spike_lims are [0 0] when this is called
                            % When re-called in batchmode, vert limits are as they
                            % appear on disply - i think the negative of how they are
                            % actually used with the inverted waveform
                            disp('Hack to preserve sign of vert_spike_lims when running in batchmode')
                            [spike_mat,id_peaks_save]=generate_spike_mat(dat,TH,sort(invert_sign*vert_spike_lims),id_peaks,nsamp_wave,invert_sign);
                        else
                            [spike_mat,id_peaks_save]=generate_spike_mat(dat,TH,vert_spike_lims,id_crossings,nsamp_wave,invert_sign);
                        end
                    end
                    
                    save spike_mat_tmp spike_mat
                    
                    mean_wave = mean(spike_mat);
                    
                    if show_timing
                        disp(['Elapsed time generating spike_mat= ' num2str(toc)]);
                    end

                    %     PRINCOMP(spike_mat) performs PCA on the N-by-P matrix of spikes  and returns the principal component coefficients
                    %     Each row of spike_mat is a waveform   COEFF is a P-by-P matrix, each column containing coefficients
                    %     for one principal component.  The columns are in order of decreasing   component variance.
                    %
                    % SCORE is the representation of X in the principal component space
                    %  Rows of SCORE correspond to observations, columns to components
                    
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end
                    
                    % BC: 4/12/2018
                    % Modified: UPDATED FUNCTION CALL
                    % [COEFF, SCORE,LATENT] = princomp(spike_mat);
                    [COEFF, SCORE,LATENT] = pca(spike_mat);
                    
                    
                    if show_timing
                        disp(['Elapsed time running PCA= ' num2str(toc)]);
                    end

                    %             % run this on the file used to set prefs to show that it works
                    %                         [n,p] = size(spike_mat);
                    %                         % subtract mean from each column of spike_mat
                    %                         spike_mat_x0 = spike_mat- repmat(mean(spike_mat,1),n,1);tmp_SCORE=spike_mat_x0*Data.COEFF_matrix;
                    %                         tmp_SCORE(1:5,1:5);SCORE(1:5,1:5)
            
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end
                    
                    % cluster points in component space, using only the first two
                    % components
                    
                    [IDX,C] = kmeans(SCORE(:,1:2), n_clusters);
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
                        component{x}=mean_wave+C(x,:)*COEFF(:,1:2)';
                    end

                    %%%%%%%%%%%%%%%%%%%
                    for x=1:n_clusters
                        wave_tmp=mean_wave+C(x,:)*COEFF(:,1:2)';
                        max_abs_wave(x)=max(abs(wave_tmp));
                    end

                    redo_peaks = false;
                end

                if replot
                    cla;
                    if show_timing
                        tic;
                        disp('Starting timer...');
                    end
                    
                    % Plot waveforms
                    n_waves_to_show=300;
                    
                    clear spike_mat_tmp
                    
                    if vert_spike_lims_approved & horiz_spike_lims_approved
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
                    t_str{1}=RemoveUnderScore(RHD_name);
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
                    set(gca,'xlim',[min(xdat) max(xdat)]);
                end
            end
        end
    end
    
end

