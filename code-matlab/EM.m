function [x_curr, rho_curr] = EM(patches_all, x_init, rho_init, L, ...
    rotations, sigma2, jball, Psilms, ell_max, s_lens, shifts, ...
    NUM_ITERS, ATA, ATI)
% An approximate expectation-maximization (EM) algorithm for
% 3D reconstruction of a volume directly from cryo-EM micrographs.
%
% Inputs:
%
% patches_all: a 3D array of 2D patches
% x_init: an array of 3D coefficients of the initial volume estimate
% rho_init: a 2D array of the initial estimates of the probabilities of
% shifts in the patches
% L: the size of the volume along each dimension
% rotations: a 3D array of rotation matrices
% sigma2: the variance of the Gaussian noise
% jball: a vector of spherical Bessel function indices
% Psilms: a cell array of spherical harmonic functions
% ell_max: the maximum spherical harmonic degree used in the reconstruction
% s_lens: a vector of the sizes of s indices
% shifts: a list of the possible 2D shifts
% NUM_ITERS: number of iterations to perform
% ATA: a precomputation of the multiplication between A^T and A
% A_I: a precomputation of the multipication between A^T and I
%
% Outputs:
%
% x_curr: an array of 3D coefficients of the final volume estimate
% rho_curr: a 2D array of the final estimates of the probabilities
% of shifts in the patches

% Get the number of patches
NUMBER_OF_PATCHES = size(patches_all, 1);

% Calculate the number of patches to use in each iteration
NUMBER_OF_PATCHES_PER_ITERATION = NUMBER_OF_PATCHES;
% floor(STOCHASTIC_FACTOR * ...
%     NUMBER_OF_PATCHES);

% Get the number of rotations
K = size(rotations, 3);

% Set the current x and rho to the initial values
x_curr = x_init;
rho_curr = rho_init;

% Initialize the count to 0, and set the number of iterations to 1000
count = 0;

% Get the total number of shifts
NUM_SHIFTS = size(shifts, 1);

% Start an infinite loop
while count <= NUM_ITERS
    % Increment the count
    count = count + 1;

    % Get the current volume by computing the inverse Fourier transform of
    % the expansion of x with respect to Psilms
    volume_curr = real(icfftn(expand_vol_psilms(x_curr, floor(L / 2), ...
        jball, Psilms, L)));
    norm(volume_curr, "fro")
    save("volume_curr" + num2str(count-1) + ".mat", "volume_curr")

    % Select a random subset of patches to use in the current iteration
    % SAMPLES = randsample(NUMBER_OF_PATCHES, ...
    % NUMBER_OF_PATCHES_PER_ITERATION);
    % patches_all_sampled = patches_all(SAMPLES, :, :);
    patches_all_sampled = patches_all;
    Nd = NUMBER_OF_PATCHES_PER_ITERATION;
    % patches_reshaped = permute(reshape(patches_all_sampled, Nd, L^2), [2, 3, 4, 5, 1]);

    %% E-step
    % TODO: Parallelize.
    tic
    pI_curr = calc_probabilities(L, patches_all_sampled, volume_curr, ...
        rotations, sigma2);
    toc
    likelihood_func_l_rot = pI_curr .* permute(rho_curr, [3, 1, 2]);
    pl_rot_curr = likelihood_func_l_rot ./ sum(likelihood_func_l_rot, ...
        [1, 2, 3]);
    pl_rot_curr_reshaped = reshape(pl_rot_curr, 1, 1, K, (2 * L) ^ 2, Nd);

    %% M-step
    %% Update rho
    rho = squeeze(sum(pl_rot_curr, [1, 4]));
    rho = rho / sum(rho, "all");
    upsilon = rho(L + 1, L + 1); % Density of empty patches.
    % We assume all other shifts are equally probable.
    rho_updated = zeros(2 * L, 2 * L);
    for l=1:NUM_SHIFTS
        if (shifts(l, 1) == L) && (shifts(l, 2) == L)
            rho_updated(shifts(l, 1) + 1, shifts(l, 2) + 1) = upsilon;
        else
            rho_updated(shifts(l, 1) + 1, ...
                shifts(l, 2) + 1) = (1 - upsilon) / (NUM_SHIFTS ...
                - 1 - 2 * (2 * L - 1));
        end
    end

    %% Update x
    tic
    A_k = sum(pl_rot_curr_reshaped .* ATA, [3, 4, 5]);
    y_k = sum(pl_rot_curr_reshaped .* ATI, [3, 4, 5]);
    x_updated = A_k \ y_k;
    toc

    % %% Calculate log-likelihood
    % Ax = pagemtimes(PR_matrix_all_rotations_all_shifts, vec_cell(x_curr));
    % norm_diff = sum(abs(patches_reshaped - Ax) .^ 2, [1, 2]);
    % Q = (1 / Nd) * sum(pl_rot_curr_reshaped .* norm_diff, "all");
    % display(Q)

    %% Update for next iteration
    x_curr = cell_vec(x_updated, x_curr, ell_max, s_lens);
    rho_curr = rho_updated;

end
end
