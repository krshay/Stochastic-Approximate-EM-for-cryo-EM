warning('off');
addpath(genpath('../../ASPIRE'))

vol = load("/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM/code-matlab/volume_curr28.mat");
vol = vol.volume_curr;
WriteMRC(vol, 1, "/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM volumes/volume_curr28_pinv.mrc");
