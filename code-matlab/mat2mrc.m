warning('off');
addpath(genpath('../../ASPIRE'))

vol = load("/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM volumes/results_L6_pinv/volume_curr10.mat");
vol = vol.volume_curr;
WriteMRC(vol, 1, "/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM volumes/volume_curr10_pinv.mrc");
