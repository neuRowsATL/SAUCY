function unipolar_analysis_temp(fname,display_filter,plotflag)
%displays unipolar data in one of five ways. If display_filter is a
%number, then it compares every unipolar channel to that channel number. If
%display_filter is the string 'unipolar', then all 16 channels are
%displayed without any subtraction. If display_filter='mean_subtract', then the mean of all unipolar signals is subtracted from each signal
%If display_filter=='bipolar', the 16
%subtractions that normally occur in the bipolar amplifer are shown (they
%can be changed below if desired). If display_filter=='laplace', the
%channels are subtracted using a 2D laplace filter. If
%display_filter=='d_laplace', diagonals are included in the laplace filter.
if nargin==2
    plotflag=1;
end

n_chan=16;
k=factor(n_chan);

if length(k)==1
    r=k;
    c=1;
elseif length(k)==2
    r=k(2);
    c=k(1);
else
    r=k(ceil(end/2):end);
    c=k(1:floor(end/2));
end
r=4;
c=4;
% time_ids=1:400000;
if ~isempty(strfind(fname,'FILTERED'))
    already_filtered=1;
    load(fname)
    if ~exist('t_amplifier','var')
        t_amplifier=t_board_adc;
    end
    FS=round(1/(t_amplifier(2)-t_amplifier(1)));
else
    already_filtered=0;
    id=strfind(fname,'.mat');
    filtname=[fname(1:id-1) '_WITH_FILTERED_300_7500.mat'];
    if ~exist(filtname,'file')
        load(fname)
        if ~exist('t_amplifier','var')
            t_amplifier=t_board_adc;
        end
        FS=round(1/(t_amplifier(2)-t_amplifier(1)));
        n_sec=60;
        sample_id=1:FS*n_sec;
        for x=1:n_chan
%             filt_data_matrix(x,:)=bandpass_filtfilt(amplifier_data(x,:),FS,300,7500);disp(['filtering channel '  num2str(x)])
            filt_data_matrix(x,:)=bandpass_filtfilt(amplifier_data(x,sample_id),FS,300,7500);disp(['filtering channel '  num2str(x)])
        end
        save(filtname,'filt_data_matrix','t_amplifier','sample_id','t_board_adc')
    else
        load(filtname)
    end
    FS=round(1/(t_amplifier(2)-t_amplifier(1)));
end

t_amplifier=t_amplifier(sample_id);
t_board_adc=t_board_adc(sample_id);

t_amplifier=t_amplifier-t_amplifier(1);
t_board_adc=t_board_adc-t_board_adc(1);

figure
if isnumeric(display_filter) % plot channel specified by "channel_to_compare" against the other 15
    channel_to_compare=display_filter;
    for x=1:n_chan
        ax(x)=subplot(r,c,x);
        if x==channel_to_compare
            plot(t_amplifier,filt_data_matrix(channel_to_compare,:),'r');
            title('Unipolar, others subtracted from this one')
        else
            plot(t_amplifier,filt_data_matrix(channel_to_compare,:)-filt_data_matrix(x,:),'k')
            title(['Ch ' num2str(channel_to_compare) ' - Ch ' num2str(x)])
        end
    end
elseif strcmp(display_filter,'unipolar') % show all unipolar recordings
    for x=1:n_chan
        ax(x)=subplot(4,4,x);
        plot(t_amplifier,filt_data_matrix(x,:),'r');title(['Unipolar channel ' num2str(x)])
    end
elseif strcmp(display_filter,'mean_subtract') % show all unipolar recordings
    for x=1:n_chan
        ax(x)=subplot(r,c,x);
        mean_data=mean(filt_data_matrix);
        plot(t_amplifier,filt_data_matrix(x,:)-mean_data,'r');title(['Mean Subtracted Unipolar channel ' num2str(x)])
    end
elseif strcmp(display_filter,'bipolar')
    % look at a bunch of paired channels
%    pairs_to_use=[1 2;3 4;5 6;7 8; 9 10; 11 12;13 14;15 16];
    pairs_to_use=[1 2;3 4;5 6;7 8;9 10;11 12;13 14;15 16;1 13;5 9;3 10;7 15;4 6;11 16;2 14;8 12];
%     pairs_to_use=[1 3;5 4;13 9;10 7;15 6;11 16;8 12;14 2];
%     pairs_to_use=[1 2;16 3;7 4;13 10;8 5;14 11;9 6;15 12;4 2;6 5;7 1;9 8;10 3;12 11;13 16;15 14];
%     pairs_to_use=[2 3;16 1;10 4;13 7;11 5;14 8;12 6;15 9;5 2;6 4;8 1;9 7;11 3;12 10;14 16;15 13];
%     pairs_to_use=[16 15;1 14;10 13;4 7;9 12;3 6;8 11;2 5;13 15;11 12;10 16;8 9;7 14;5 6;4 1;2 3];
%     pairs_to_use=[15 12;1 2;5 6;13 7;11 5;14 8;12 6;15 9;5 2;6 4;8 1;9 7;11 3;12 10;14 16;15 13];
%     pairs_to_use=[9 7;3 4;4 2;5 7;7 6;6 8;9 10;10 11;11 12;13 15;15 16;16 14; 16 10; 1 5; 1 14; 10 7];
%     pairs_to_use=[1 5; 5 9;9 13;3 7;7 10; 10 15;4 6; 6 11;11 16;2 8;8 12;12 14;1 7;4 8;9 15; 11 14];
%     pairs_to_use=[1 5; 1 9;1 13;5 9;5 13; 9 13;2 6; 2 10;2 14;6 10;6 14;10 14;1 4;5 8;9 12; 13 16];
    for x=1:length(pairs_to_use)
        ax(x)=subplot(4,4,x);
            plot(t_amplifier,filt_data_matrix(pairs_to_use(x,1),:)-filt_data_matrix(pairs_to_use(x,2),:),'k')
            title(['Ch ' num2str(pairs_to_use(x,1)) ' - Ch ' num2str(pairs_to_use(x,2))])
    end
elseif strcmp(display_filter,'laplace')
        laplace_data=zeros(4,length(t_amplifier));
        laplace_data(1,:)=filt_data_matrix(2,:)+filt_data_matrix(5,:)+filt_data_matrix(7,:)+filt_data_matrix(10,:)-4*filt_data_matrix(6,:);
        laplace_data(2,:)=filt_data_matrix(3,:)+filt_data_matrix(6,:)+filt_data_matrix(8,:)+filt_data_matrix(11,:)-4*filt_data_matrix(7,:);
        laplace_data(3,:)=filt_data_matrix(6,:)+filt_data_matrix(9,:)+filt_data_matrix(11,:)+filt_data_matrix(14,:)-4*filt_data_matrix(10,:);
        laplace_data(4,:)=filt_data_matrix(7,:)+filt_data_matrix(10,:)+filt_data_matrix(12,:)+filt_data_matrix(15,:)-4*filt_data_matrix(11,:);
        if already_filtered
            save(fname,'laplace_data','-append')
        else
            save(filtname,'laplace_data','-append')
        end
    if plotflag
        ax(1)=subplot(2,2,1);
        plot(t_amplifier,laplace_data(1,:),'k')
        title(['Ch6 2D Laplace filter'])
        ax(2)=subplot(2,2,2);
        plot(t_amplifier,laplace_data(2,:),'k')
        title(['Ch7 2D Laplace filter'])
        ax(3)=subplot(2,2,3);
        plot(t_amplifier,laplace_data(3,:),'k')
        title(['Ch10 2D Laplace filter'])
        ax(4)=subplot(2,2,4);
        plot(t_amplifier,laplace_data(4,:),'k')
        title(['Ch11 2D Laplace filter'])
    end
elseif strcmp(display_filter,'d_laplace')
    laplace_data=zeros(4,length(t_amplifier));
    laplace_data(1,:)=0.5*filt_data_matrix(1,:)+0.5*filt_data_matrix(3,:)+0.5*filt_data_matrix(9,:)+0.5*filt_data_matrix(11,:)+filt_data_matrix(2,:)+filt_data_matrix(5,:)+filt_data_matrix(7,:)+filt_data_matrix(10,:)-6*filt_data_matrix(6,:);
    laplace_data(2,:)=0.5*filt_data_matrix(2,:)+0.5*filt_data_matrix(4,:)+0.5*filt_data_matrix(10,:)+0.5*filt_data_matrix(12,:)+filt_data_matrix(3,:)+filt_data_matrix(6,:)+filt_data_matrix(8,:)+filt_data_matrix(11,:)-6*filt_data_matrix(7,:);
    laplace_data(3,:)=0.5*filt_data_matrix(5,:)+0.5*filt_data_matrix(7,:)+0.5*filt_data_matrix(13,:)+0.5*filt_data_matrix(15,:)+filt_data_matrix(6,:)+filt_data_matrix(9,:)+filt_data_matrix(11,:)+filt_data_matrix(14,:)-6*filt_data_matrix(10,:);
    laplace_data(4,:)=0.5*filt_data_matrix(6,:)+0.5*filt_data_matrix(8,:)+0.5*filt_data_matrix(14,:)+0.5*filt_data_matrix(16,:)+filt_data_matrix(7,:)+filt_data_matrix(10,:)+filt_data_matrix(12,:)+filt_data_matrix(15,:)-6*filt_data_matrix(11,:);
    if already_filtered
        save(fname,'laplace_data','-append')
    else
        save(filtname,'laplace_data','-append')
    end
    if plotflag
        ax(1)=subplot(2,2,1);
        plot(t_amplifier,laplace_data(1,:),'k')
        title(['Ch6 2D Diagonal Laplace filter'])
        ax(2)=subplot(2,2,2);
        plot(t_amplifier,laplace_data(2,:),'k')
        title(['Ch7 2D  Diagonal Laplace filter'])
        ax(3)=subplot(2,2,3);
        plot(t_amplifier,laplace_data(3,:),'k')
        title(['Ch10 2D Diagonal Laplace filter'])
        ax(4)=subplot(2,2,4);
        plot(t_amplifier,laplace_data(4,:),'k')
        title(['Ch11 2D Diagonal Laplace filter'])
    end
end
if plotflag
    linkaxes(ax,'x')
end

%
% % subtract mean of other channels(not great)
% for x=1:16
% ax(x)=subplot(4,4,x)
% id_not_x=find([1:16]~=x);
% plot(t_amplifier,filt_data_matrix(1,:)-mean(filt_data_matrix(x,id_not_x)),'k')
% end
% linkaxes(ax,'x')
