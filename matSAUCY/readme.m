%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------
% S.A.U.C.Y. - Sam's Ad-hoc Unit ClassYfier for use with INTAN recorded
% data
% Last modified: April 17, 2018
% Last modified by: Bryce Chung <bpchung@emory.edu>
% ---------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PROCEDURE:
%   ** Procedure is also outlined in test_SAUCY.m **
%
%   SAUCY()
%   (1)  Create SAUCY object for a single file:
%        mySaucy = SAUCY('experiment_name', ch_num);
%
%   (2)  Load data into program and filter raw data:
%        load_data(mySaucy, 'data_file_name.rhd');
%
%   (3)  Run PCA-based analysis on single file to find parameters for analysis:
%        do_saucy(mySaucy);
%
%   SAUCY_batch()
%   (1)  Run PCA-based analysis on all *bin files in the current
%   directory:  e.g. SAUCY_beta('name.cbin',ch_num,'batchmode')
%
%   (2)  Check how well algoritm did across all files:
%   SAUCY_beta('name.cbin',ch_num,'trackmode')
%   or for a single file:
%   SAUCY_beta('name.cbin',ch_num,'checkmode')
%
%   (3)  If syllables are labelled (in .not.mat files), look at PSTHs:
%       neural_by_syl_seq_PCA('[sequence of labels]',ch_num)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
