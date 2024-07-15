warning('off');
addpath(genpath('../../ASPIRE'))

vol1 = load("/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM volumes/results_L6_explicit_method/volume_curr_82.mat");
vol1 = vol1.volume_curr;
vol2 = load("/Users/shaykreymer/My Drive/PhD/Code/Stochastic-Approximate-EM-for-cryo-EM volumes/results_L6_pinv/volume_curr81.mat");
vol2 = vol2.volume_curr;
% vol1 =w vol1 / norm(vol1, "fro");
% vol2 = vol2 / norm(vol2, "fro");

% Align rotation
[err, min_error_rotation, vol2_rotated] = calc_volume_relative_error_rotations(vol1, vol2, 50);
[resA, fighandle] = plotFSC(vol1, vol2_rotated, 0.5, 40/11);
