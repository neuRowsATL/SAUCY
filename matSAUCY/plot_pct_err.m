function [pct_errors] = plot_pct_err( sample_neuralnotmat )
%goes through neuralnotmats matching input and plots pct error in order

id=strfind(sample_neuralnotmat,'.neuralnot');
neuralnotmat_ending=sample_neuralnotmat(id:end);
pct_errors=[];
ISI_viols=[];
ISI_TH_viols=[];
num_spikes=[];
    
!dir /B *neuralnot* > batchfile
    fid=fopen('batchfile','r');
    while 1
        fn=fgetl(fid);
        if ~ischar(fn);break;end
        if ~isempty(strfind(fn,neuralnotmat_ending))
            disp(fn);
            load(fn);
            pct_errors=[pct_errors Data.pct_error(end)];
            ISI_viols=[ISI_viols Data.pct_ISIs_leq_1ms(end-1)];
            ISI_TH_viols=[ISI_TH_viols Data.pct_ISIs_leq_1ms(end)];
            num_spikes=[num_spikes length(Data.spiketimes_from_recommended_TH)];
        end
    end
    
plot(1:length(pct_errors),pct_errors*100,'bo')
xlabel('File number')
ylabel('Percent Error (%)')
% mean_error=mean(pct_errors*100);
figure
plot(1:length(ISI_viols),ISI_viols*100,'bo',1:length(ISI_TH_viols),ISI_TH_viols*100,'ro')
xlabel('File number')
ylabel('ISI 1ms violations (%)')
legend('PCA-based','TH-based')

figure
plot(1:length(num_spikes),num_spikes,'bo')
xlabel('File number')
ylabel('num_spikes')

            


end

