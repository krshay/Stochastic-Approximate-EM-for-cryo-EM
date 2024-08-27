%% Script to generate mrc files from mat files.

warning('off');
addpath(genpath('../../ASPIRE'))

vol = load("../data-mat/volume_curr3.mat");
vol = vol.volume_curr;
WriteMRC(vol, 1, "../data-mrc/volume_curr3.mrc");
