function  show_batch_cluster_drift(Data);
subplot(2,2,1);hold on;
col_vec='rgbc';
col_mat=[1 0 0;0 1 0;0 0 1;1 1 0];
N_files=length(Data);
const_vect=.5:.5/(N_files-1):1;
n_clusters=length(Data(1).cov_matrices);
%
% put this in later - a list (batchfile?) of filenames that passed or are
% borderline
%
% passed_filenames_arr=[];
% passed_filenames_pct=[];
% boarderline_filenames_arr=[];
% boarderline_filenames_pct=[];
% failed_filenames_arr=[];
% failed_filenames_pct=[];
for y=1:length(Data(1).cov_matrices)
    for x=1:length(Data)
        cov_cent_matrix_arr{y}(x,:)=Data(x).cov_centers{y};
        ellipse_handle=plotcov_SAUCY(Data(x).cov_centers{y},Data(x).cov_matrices{y},col_mat(y,:)*const_vect(x),.0455);
        % if error is less than 1%;plot as filled circle
        if Data(x).pct_error(y)<.01
            plot(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),'ko','markerfacecolor',col_mat(y,:)*const_vect(x))
        elseif .02>Data(x).pct_error(y)>.01 % if marginal - plot gray
            plot(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),'ko','markerfacecolor',[.7 .7 .7])
        end
        ydat=get(ellipse_handle,'YData');xdat=get(ellipse_handle,'XData');
        id=find(ydat==min(ydat));
        if x==1
            text(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),'  (First file)','color',col_mat(y,:)*const_vect(x),'fontsize',7)
        elseif x==length(Data)
            text(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),'  (Last file)','color',col_mat(y,:)*const_vect(x),'fontsize',7)
        end
        plot(cov_cent_matrix_arr{y}(:,1),cov_cent_matrix_arr{y}(:,2),'color',col_mat(y,:))
        subplot(2,n_clusters,n_clusters+y);hold on;
        if Data(x).pct_error(y)<.01
            plot(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),'ko','markerfacecolor',col_mat(y,:)*const_vect(x))
        elseif .02>Data(x).pct_error(y)>.01 % if marginal - plot gray
            plot(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),'ko','markerfacecolor',[.7 .7 .7])
        else  % if not sig - plot white
            plot(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),'o','color',col_mat(y,:)*const_vect(x))
        end
        text(Data(x).cov_centers{y}(1),Data(x).cov_centers{y}(2),RemoveUnderScore(Data(x).filename),'color',col_mat(y,:)*const_vect(x),'fontsize',7,'horizontalalignment','center','verticalalignment','bottom')
        plot(cov_cent_matrix_arr{y}(:,1),cov_cent_matrix_arr{y}(:,2),'color',col_mat(y,:))
        subplot(2,2,1);hold on;
    end
end
title('Center filled in ellipse color - less than 1% error      Filled in gray - 1-2% error','fontweight','bold','fontsize',12)

subplot(2,2,2);hold on;
for x=1:length(Data)
    plot(Data(x).pct_error(end)*100,Data(x).pct_ISIs_leq_1ms(end-1)*100,'o')
end
for x=1:length(Data)
    plot(Data(x).pct_error(end)*100,Data(x).pct_ISIs_leq_1ms(end)*100,'ro')
end
xlabel('% overlap of largest cluster')
ylabel('% ISI violations of largest cluster')
t_str{1}='Blue - spikes from cluster';
t_str{2}='Red - spikes from suggested threshold';
title(t_str)



%
% put this in later - a list (batchfile?) of filenames that passed or are
% borderline
%
passed_filenames_arr=[];
passed_filenames_pct=[];
borderline_filenames_arr=[];
borderline_filenames_pct=[];
failed_filenames_arr=[];
failed_filenames_pct=[];
for x=1:length(Data)
    cov_cent_matrix_arr{1}(x,:)=Data(x).cov_centers{end};
    % if error is less than 1%;plot as filled circle
    if Data(x).pct_error(end)<.01
        passed_filenames_arr{end+1}=Data(x).filename;
        passed_filenames_pct(end+1)=Data(x).pct_error(end);
    elseif .02>Data(x).pct_error(y)>.01 % if marginal - plot gray
        borderline_filenames_arr{end+1}=Data(x).filename;
        borderline_filenames_pct(end+1)=Data(x).pct_error(end);
    else
        failed_filenames_arr{end+1}=Data(x).filename;
        failed_filenames_pct(end+1)=Data(x).pct_error(end);
    end
end
whos passed_filenames_pct passed_filenames_arr Data
clear id
[tmp,id{1}]=sort(passed_filenames_pct);
passed_filenames_pct=passed_filenames_pct(id{1});
[tmp,id{2}]=sort(borderline_filenames_pct);
borderline_filenames_pct=borderline_filenames_pct(id{2});
[tmp,id{3}]=sort(failed_filenames_pct);
failed_filenames_pct=failed_filenames_pct(id{3});
for x=1:length(passed_filenames_arr)
    passed_filenames_arr_new{x}=passed_filenames_arr{id{1}(x)};
end
for x=1:length(borderline_filenames_arr)
    borderline_filenames_arr_new{x}=borderline_filenames_arr{id{2}(x)};
end
for x=1:length(failed_filenames_arr)
    failed_filenames_arr_new{x}=failed_filenames_arr{id{3}(x)};
end

disp('Passing files:')
for x=1:length(passed_filenames_arr)
    disp([passed_filenames_arr_new{x} '     ' num2str(passed_filenames_pct(x))])
end
disp(' ');disp('Borderline files:')
for x=1:length(borderline_filenames_arr)
    disp([borderline_filenames_arr_new{x} '     ' num2str(borderline_filenames_pct(x))])
end
disp(' ');disp('Failed files:')
for x=1:length(failed_filenames_arr)
    disp([failed_filenames_arr_new{x} '     ' num2str(failed_filenames_pct(x))])
end