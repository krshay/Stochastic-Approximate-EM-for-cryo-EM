function [Y, placed, projection] = generate_clean_micrograph_2d_rots_new(volume, W, N, m, seed)

% Inputs:
%   X: signal of size LxL
%   parameters : the parameters define the problem
%
% Outputs:
%   Y: micrograph of size NxN (clean)
%   placed: actual number of occurrences of the signal X in Y.

% Since we attempt placement at random, we need the end-result to be
% quite sparse to have decent chances of finding m spots.

%if N^2 < 10*m*W^2
%    warning('BigMRA:gendata2D', ...
%        'It may be difficult to get this many repetitions...');
%end
rng(seed);


projection = 0;
rangeW = 0:(W-1);

% The mask has the same size as the micrograph. Each pixel in the mask
% corresponds to a possible location for the upper-left corner of the
% signal, to be placed in the micrograph. A value of zero means the
% spot is free, while a value of 1 means the spot is forbidden (either
% because a signal will already be there, or because it is in the area
% of separation of another signal.)
mask = zeros(N, N);
% The locations table recors the chosen signal locations.
locations = zeros(m, 2);
% This counter records how many signals we successfully placed.
placed = 0;
% Since placement is random, there is a chance that we will get a lot
% of failed candidates. To make up for it, we allow for more than m
% trials. But we still put a hard cap on it to avoid an infinite loop.
%max_trials = 2*m;
max_trials = 4*m;

rotations = rand_rots(max_trials);
K = size(rotations, 3);

for counter = 1 : max_trials
    
    % Pick a candidate location for the upper-left corner of the signal
    candidate = randi(N-W+1, 1, 2);
    
    % Check if there is enough room, taking the separation rule into
    % account. That is, a square of size WxW with upper-left corner
    % specified by the candidate must be entirely free.
    if all(mask(candidate(1)+rangeW, candidate(2)+rangeW) == 0)
        
        % Record the successful candidate
        placed = placed + 1;
        locations(placed, :) = candidate;
        
        % Mark the area as reserved
        mask(candidate(1)+rangeW, candidate(2)+rangeW) = 1;
        
        % Stop if we placed sufficiently many signals successfully.
        if placed >= m
            break;
        end
        
    end
    
end

% Now that we have a list of locations, actually go ahead and build the
% micrograph.
L = size(volume, 1);
assert(size(volume, 2) == L, 'X must be square.');
rangeL = 0 : (L-1);
Y = zeros(N);
%    Y = zeros(N, N,'gpuArray');
projections = vol_project(volume, rotations( :, :, randi(K, [1 placed])));
for k = 1 : placed
%              projection = vol_project(volume, rotations( :, :, randi(K)));
%         I(candidate(1):candidate(1)+L-1, ...
%             candidate(2):candidate(2)+L-1) = projection;
    Y(locations(k, 1) + rangeL, locations(k, 2) + rangeL) = projections( :, :, k);
    
end

end