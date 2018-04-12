% Initiate SAUCY object for experiment: bl21lb21_171218
mySaucy = SAUCY('bl21lb21_171218', 8);

load_data(mySaucy, 'bl21lb21_171218_134533.rhd');
filter_data(mySaucy);
set_threshold(mySaucy, 'filt');

set_spike_mat(mySaucy, 'threshold');
optimize_clusters(mySaucy);