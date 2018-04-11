%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------
% S.A.U.C.Y. - Sam's Ad-hoc Unit ClassYfier for use with INTAN recorded
% data
% ---------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PROCEEDURE:
%   (1).  First, set parameters for PCA by running this program on a single
%   file - e.g. SAUCY_beta('name.cbin',ch_num,n_clusters)
%
%   (2).  Run PCA-based analysis on all *bin files in the current
%   directory:  e.g. SAUCY_beta('name.cbin',ch_num,'batchmode')
%
%   (3).  Check how well algoritm did across all files:
%   SAUCY_beta('name.cbin',ch_num,'trackmode')
%   or for a single file:
%   SAUCY_beta('name.cbin',ch_num,'checkmode')
%
%   (4)  If syllables are labelled (in .not.mat files), look at PSTHs:
%       neural_by_syl_seq_PCA('[sequence of labels]',ch_num)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax:     SAUCY_beta(INPUT_1,INPUT_2,INPUT_3,INPUT_4)
%
% INPUT 1: filename: either an .rhd file or a .neuralnot_chX.mat file
%
% INPUT 2:  ch_num - this is channel number from the beginning, starting with zero.  so if 3 channels were recorded,
% this input must be the number 3 rather than "obs0r". This can also use
% the difference between two channels. If so, ch_num is a vector with
% length 2, where the second channel is subtracted from the first.
% 
% INPUT 3: If the RHD file has previously been converted to a matlab file, set to 1, otherwise 0. 
%
%           INPUT 4 CAN BE EITHER NUMERIC (for running this on a single
%           file) OR A TEXT STRING (for running in batchmode)
%
%   FOR A SINGLE FILE:
% INPUT 4:  n_clusters  - number of clusters to group in principal
% components space.  Since one cluster will be the 'noise' waveforms,
% n_clusters=(#_presumptive_spikes + 1). Running the program this way saves a file called
% RHD_name.neuralnot_CHX.mat, where X is the channel number specified by
% INPUT 2.  See below for a list of what is saved in this file.  If only
% two inputs are give, n_clusters will default to 2 (i.e. one spike cluster
% and one noise cluster).
%
%   IN BATCHMODE:
% INPUT 4:  one of three different text strings specifiying a batchmode
%
%       OPTION 1:   'batchmode' : will run through all *bin files in the
%       current directory and run the program, using the settings saved in
%       RHD_name.neuralnot_CHX.mat.  (So you have to have already run thism
%       program on a single file before running it in batchmode).
%
%       OPTION 2:   'checkmode': allows you to visualize the program output
%       for a single file (named in INPUT 1) that ALREADY HAS a .neuralnot_CHX.mat file.
%
%       OPTION 3:  'trackmode' - plots the contours of clusters across all
%       .neuralnot_CHX.mat files in the current directory
%
%       OPTION 4:  'addlatentmode' - uses saved settings to add the field LATENT_RECONSTRUCTED
%        (see below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The file RHD_name.neuralnot_CHX.mat contains a structure array called
%   'Data' with the following fields:
%
% Data.filename
%        Name of the analyzed file
% Data.n_units
%       Number of non-noise clusters.  This is equal to (#_clusters - 1)
% Data.vert_spike_lims
%       Vertical limit for spike waveforms.  Any waveforms that exceed
%       these will be disincluded.  Used to screen out mvmt artifacts,
% Data.nsamp_wave
%       Sets the horizontal limits on spike waveforms.  A 1x2 vector,
%       [#_samples_before_peak #_samples_after_peak]
% Data.run_bandpass
%       Whether or not to run bandpass filter on raw neural trace.  Used to
%       screen out noise and movement artifact.
% Data.bandpass_freqs
%       If Data.run_bandpass==1, this is a 1x2 vector specifying
%       [low_limit high_limit] for the bandpass
% Data.TH
%       Threshold for waveform selection.  Waveforms are selected based on
%       whether they cross this waveform within the horizontal window
%       specified by DATA.nsame_wave and do not exceed the vertical limits
%       set by Data.vert_spike_lims
% Data.spiketimes
%       A 1xN array of spiketimes (in seconds), where N is the total number of clusters.
%       These are ordered so that Data.spiketimes{1} is the cluster with
%       the smallest mean amplitude (i.e. the noise cluster), and
%       Data.spiketimes{end} is the cluster with the largets spikes.
% Data.recommended_TH
%       Threshold for spiketimes recommended based on an analysis of the
%       amplitude distribution of the spikes in each cluster.  See
%       threshold_strategy in code below for how recommended_TH is
%       computed.
% Data.spiketimes_from_recommended_TH
%       Spiketimes resulting from above.  Note that these are spiketimes
%       for the largets unit only.
% Data.mean_waveform
%       A 1xN array of mean waveforms in each cluster, where N is the
%       total number of clusters.
% Data.std_waveform
%       A 1xN array of standard deviations of the waveforms in each
%       cluster, where N is the total number of clusters.
% Data.COEFF_matrix
%       The matrix of coefficients determined by PCA (implemented with
%       PRINCOMP.m).  Each columns is a basis vector for the matrix of
%       spike waveforms in the new coordinate system, and columns are
%       arranged in order of decreasing component variance.
% Data.SCORE
%       The projection of the data in principal component space.
%       SCORE=spike_mat_x0*COEFF, where spike_mat_x0 is the mean-subtracted
%       matrix of spike waveforms.
% Data.cov_matrices
%       These are the 2x2 covariance matrices of the N clusters.  These are
%       the cov matrices of the first two components, NOT the cov mats of
%       the original data.
% Data.cov_centers
%       Centers of the N clusters.
% Data.pct_error
%       Percent overlap of each cluster with the other clusters.
% Data.total_n_samples
%       Total number of samples in the file.
% Data.Fs
%       Sampling frequence (Hz) of the file
% Data.warning_str
%       Text strings containing messages about anything funny that happened
%       when analyzing the file.
% Data.inverted_waveform
%       True (=1) for downward for downward-pointing spikes.  Used to tell
%       program to flip waveform at certain points.
% Data.pct_ISIs_leq_1ms
%       Percentage of interspike intervals (ISIs) that are less than or
%       equal to 1 ms.  This is an 1x(N+1) vector.  The first N entries
%       correspond to the spiketimes from the N clusters.  The N+1st entry
%       is the ISI violation rate for Data.spiketimes_from_recommended_TH
%    Data.LATENT, Data.LATENT_RECONSTRUCTED
%         These are both vectors of how much variance is accounted for by each principal component.
%         Data.LATENT is the output from PCA run when not in ANY kind of batchmode, Data.LATENT_RECONSTRUCTED
%         is the same thing reconstructed (using code copied from princomp.m) in addlatentmode.  These yeild the same
%         values
