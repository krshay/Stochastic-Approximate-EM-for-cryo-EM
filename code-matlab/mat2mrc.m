warning('off');
addpath(genpath('../../ASPIRE'))

vol = load("/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM/code-matlab-server/volume_curr1.mat");
vol = vol.volume_curr;
WriteMRC(vol, 1, "/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM volumes/volume_curr1_pinv.mrc");
