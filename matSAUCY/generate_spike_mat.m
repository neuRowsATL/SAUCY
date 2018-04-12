function [spike_mat,id_peaks_save]=generate_spike_mat(dat,TH,vert_spike_lims,id_peaks,nsamp_wave,invert_sign)

    % first, disinclude peaks too close to start or end of file
    id_peaks=id_peaks(find(id_peaks>nsamp_wave(1) & id_peaks<(length(dat)-nsamp_wave(2))));

    % id_before_peaks is a vector of the indexes of datapoints before the
    % peaks
    id_before_peaks=zeros(nsamp_wave(1),length(id_peaks));  % do memory allocation all at once
    id_before_peaks(1,:)=id_peaks-1;
    for x=2:nsamp_wave(1)
        id_before_peaks(x,:)=id_peaks-x;
    end
    % here, we orgainze these into a matrix:
    % each row is a sample number before peak
    % each column in for one entry in id_peaks
    id_before_peaks_mat=reshape(id_before_peaks,nsamp_wave(1),length(id_peaks));

    % id_before_peaks is a vector of the indexes of all datapoints withing
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
    spike_mat=reshape(id_waveforms_included, size(id_waveforms_included,1), sum(nsamp_wave)+1);

    id_peaks_save=ids_all_crit';
end
