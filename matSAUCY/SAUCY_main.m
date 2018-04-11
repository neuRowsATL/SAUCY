function SAUCY(RHD_name,ch_num,subtract,RHD_already_loaded,n_clusters)

if nargin==2,n_clusters=2;
    RHD_already_loaded=0;
    analyze(RHD_name,ch_num,subtract,RHD_already_loaded,n_clusters,[],0);
elseif nargin==3
    n_clusters=2;
    analyze(RHD_name,ch_num,subtract,RHD_already_loaded,n_clusters,[],0);
elseif strcmp(n_clusters,'batchmode') | strcmp(n_clusters,'trackmode') | strcmp(n_clusters,'checkmode') | strcmp(n_clusters,'addlatentmode') % | strcmp(n_clusters,'searchmode')
    if length(ch_num)==1
        file_for_parameters=[RHD_name '.neuralnot_CH'  num2str(ch_num) '.mat'];
    elseif length(ch_num)==2 && subtract==1
        file_for_parameters=[RHD_name '.neuralnot_CH'  num2str(ch_num(1)) '_minus_CH' num2str(ch_num(2)) '.mat'];        
    end
    disp(['Loading parameters from file ' file_for_parameters '...'])
    load(file_for_parameters);
    
    if ~strcmp(n_clusters,'trackmode') & ~strcmp(n_clusters,'checkmode') & ~strcmp(n_clusters,'addlatentmode')% addlatentmode is NEW
        if isunix;
            !rm batchfile
            !ls -1 *.rhd > batchfile
        else
            if ~isempty(strfind(RHD_name,'.rhd'))
                !dir /B **.rhd > batchfile
            else
                !dir /B **.mat > batchfile
            end
        end
    elseif strcmp(n_clusters,'checkmode')
        if length(ch_num)==1
        	load([RHD_name '.neuralnot_CH' num2str(ch_num) '.mat']);
        elseif length(ch_num)==2 && subtract==1
            load([RHD_name '.neuralnot_CH' num2str(ch_num(1)) '_minus_CH' num2str(ch_num(2)) '.mat']);  
        end
        
        % 2/20/08 - last input was 3 - prob changed when setting up addlatentmode. put back to 2
        analyze(RHD_name,ch_num,subtract,RHD_already_loaded,Data.n_units+1,Data,2);
        return
    else strcmp(n_clusters,'trackmode') |  strcmp(n_clusters,'addlatentmode') % addlatentmode is NEW
        if length(ch_num)==1
            if isunix
                eval(sprintf('!ls  *neuralnot_CH%s* > batchfile',num2str(ch_num)))
            else
                eval(sprintf('!dir /B *neuralnot_CH%s* > batchfile',num2str(ch_num)))
            end
        elseif length(ch_num)==2
            if isunix
                eval(sprintf('!ls  *neuralnot_CH%s_minus_CH%s* > batchfile',num2str(ch_num(1)),num2str(ch_num(2))))
            else
                eval(sprintf('!dir /B *neuralnot_CH%s_minus_CH%s* > batchfile',num2str(ch_num(1)),num2str(ch_num(2))))
            end
        end
    end
    fid=fopen('batchfile','r');ct=1;
    while 1
        RHD_name=fgetl(fid);
        try
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BC Updated 2/21/2018
            %if isempty(strfind(RHD_name,'neuralnot')) && isempty(strfind(RHD_name,'wav.not.mat')) && isempty(strfind(RHD_name,'coherency')) && isempty(strfind(RHD_name,'mat.not.mat')) && isempty(strfind(RHD_name,'calibration')) && isempty(strfind(RHD_name,'FILTERED')) && isempty(strfind(RHD_name,'spikoclust')) && isempty(strfind(RHD_name,'spike_mat_tmp'))
            if isempty(strfind(RHD_name,'wav.not.mat')) && isempty(strfind(RHD_name,'coherency')) && isempty(strfind(RHD_name,'mat.not.mat')) && isempty(strfind(RHD_name,'calibration')) && isempty(strfind(RHD_name,'FILTERED')) && isempty(strfind(RHD_name,'spikoclust')) && isempty(strfind(RHD_name,'spike_mat_tmp'))

                if (~ischar(RHD_name));fclose(fid);Data=tmp;break;end
                if strcmp(n_clusters,'batchmode') | strcmp(n_clusters,'addlatentmode')
                    %        if strcmp(n_clusters,'batchmode')
                    if ct==1
                        tmp_str=[];
                        for x=1:Data.n_units+1
                            tmp_str=[tmp_str 'Cluster ' num2str(x) '    '];
                        end
                        tmp_str=[tmp_str 'TH recommended   '];
                        tmp_str=[tmp_str 'Filename'];disp(tmp_str)
                    end
                    if strcmp(n_clusters,'batchmode')
                        tmp(ct)=analyze(RHD_name,ch_num,subtract,RHD_already_loaded,Data.n_units+1,Data,1,file_for_parameters);
                    elseif strcmp(n_clusters,'addlatentmode')
                        %"addlatentmode" needs to be passed both the
                        disp(['adding latent values for ' RHD_name])
                        load(RHD_name) % this is actually the name of the current neuralnot file - to pass right structure "Data"
                        RHD_name_shortened=RHD_name(1:end-18); % acutal RHD name
                        tmp(ct)=analyze(RHD_name_shortened,ch_num,subtract,RHD_already_loaded,Data.n_units+1,Data,3);
                        %                tmp(ct)=analyze(RHD_name,ch_num,Data.n_units+1,Data,3);
                    end
                elseif strcmp(n_clusters,'trackmode')
                    display(['Analyzing file ' RHD_name]);load(RHD_name);tmp(ct)=Data;
                end
                ct=ct+1;
            end
        catch
            ct = ct + 1;
        end
    end
elseif isstr(n_clusters)
    error('bad string input for third input')
else
    Data=analyze(RHD_name,ch_num,subtract,RHD_already_loaded,n_clusters,[],0);
end

if  strcmp(n_clusters,'trackmode')% strcmp(n_clusters,'batchmode')
    show_batch_clusterdrift(Data);
end