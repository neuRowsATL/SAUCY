%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BC - April 17, 2018
%%%%% This is a test script that uses data recorded from an experiment that
%%%%% was named "bl21lb21_171218" with the data channel of interest #8. The
%%%%% raw data file recorded using Intan was saved to
%%%%% "bl21lb21_171218_134533.rhd" and is available in the working path of
%%%%% MATLAB.
%%%%%
%%%%% Because the data file is larger than github's allowable filesize,
%%%%% please contact the author of this script at bpchung[at]emory.edu.
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Initiate SAUCY object for experiment: bl21lb21_171218
mySaucy = SAUCY('bl21lb21_171218');
mySaucy.chan = 8;
load_data(mySaucy, 'bl21lb21_171218_134533.rhd');
filter_data(mySaucy);

%%%%% Run full SAUCY default analysis
%%%%% set_threshold(S);
%%%%% do_clustering(S);
%%%%% optimize_clusters(S);
% do_saucy(mySaucy);

%%%%% Load spike times and then run SAUCY clustering
load('bl21lb21_171218_CH8_spikes-th.mat');
mySaucy.chan = 8;
mySaucy.threshold = 80;
do_clustering(mySaucy, 2, round(ts*30000));
optimize_clusters(mySaucy);