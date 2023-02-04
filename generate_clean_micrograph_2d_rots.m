function [I, locations] ...
    = generate_clean_micrograph_2d_rots(volume, W, N, m, seed)
%     Form an N*N matrix containing target images at random locations and
% rotations.
%     Args:
%         volume: the L * L * L volume.
%         W: separation between images; L for arbitrary spacing
% distribution, 2*L-1 for the well-separated case.
%         N: the width and height of the required micrograph.
%         m: wanted number of images to place.
%         K: number of rotations.
%         seed: random seed (default is 1).
%
%     Returns:
%         I: N*N matrix containing target images at random
% locations and rotations.
%         locations: list of 2-D locations of the target images in I.
if nargin <= 4
    seed = 1;
end
rng(seed);
L = size(volume, 1);
m = round(m);
mask = zeros(N + W, N + W);
% The locations table records the chosen signal locations.
locations = {};
% This counter records how many signals we successfully placed.
placed = 0;
max_trials = 5*m;
I = zeros(N, N);
rotations = rand_rots(m);
for ii=1:max_trials
    % Pick a candidate location for the upper-left corner of the signal
    candidate = randi(N - W, 2, 1);
    % Check if there is enough room, taking the separation rule into
    % account. That is, a square of size WxW with upper-left corner
    % specified by the candidate must be entirely free.
    if ~any(mask(candidate(1):candidate(1)+2*W-1-1, ...
            candidate(2):candidate(2)+2*W-1-1), 'all')
        % Record the successful candidate
        locations{placed + 1} = candidate;
        % Mark the area as reserved
        mask(candidate(1):candidate(1)+W-1, candidate(2):candidate(2)+W-1) = 1;
        placed = placed + 1;
        projection = vol_project(volume, rotations( :, :, randi(placed)));
        I(candidate(1):candidate(1)+L-1, ...
            candidate(2):candidate(2)+L-1) = projection;
        % Stop if we placed sufficiently many signals successfully.
        if placed >= m
            break
        end
    end
end