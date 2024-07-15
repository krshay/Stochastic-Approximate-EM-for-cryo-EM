function x_curr = EM_stochastic_parallel_ls(patches_all, ...
    x_init, rho_init, L, rotations, sigma2, jball, Psilms, ell_max, ...
    s_lens, n_list, lms_list, shifts, psi_lsNn, STOCHASTIC_FACTOR, ...
    gs, NUM_ITERS, ls_A_dagger_patch)
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
% n_list: a list of the n indices
% lms_list: a list of the spherical harmonic degree and order
% Ls: a list of the possible 2D shifts
% psi_lsNn: a cell array of the basis functions
% STOCHASTIC_FACTOR: the stochastic factor S
% gs: the precomputed g function
% NUM_ITERS: number of iterations to perform
%
% Outputs:
%
% x_curr: an array of 3D coefficients of the final volume estimate
% rho_curr: a 2D array of the final estimates of the probabilities
% of shifts in the patches

% Get the number of patches
NUMBER_OF_PATCHES = size(patches_all, 1);

% Calculate the number of patches to use in each iteration
NUMBER_OF_PATCHES_PER_ITERATION = floor(STOCHASTIC_FACTOR * ...
    NUMBER_OF_PATCHES);

% Get the length of the initial x
length_x = length(vec_cell(x_init));

% Get the number of rotations
K = size(rotations, 3);

% Set the current x and rho to the initial values
x_curr = x_init;

% Initialize the count to 0, and set the number of iterations to 1000
count = 0;

% Get the total number of shifts
num_shifts = size(shifts, 1);

% Start an infinite loop
while count <= NUM_ITERS
    % Increment the count
    count = count + 1;

    % Get the current volume by computing the inverse Fourier transform of
    % the expansion of x with respect to Psilms
    volume_curr = real(icfftn(expand_vol_psilms(x_curr, floor(L / 2), ...
        jball, Psilms, L)));

    % Select a random subset of patches to use in the current iteration
    SAMPLES = randsample(NUMBER_OF_PATCHES, ...
        NUMBER_OF_PATCHES_PER_ITERATION);
    patches_all_sampled = patches_all(SAMPLES, :, :);
    Nd = NUMBER_OF_PATCHES_PER_ITERATION;

    %% Calculate weights
    % TODO: Parallelize.
    tic
    pI_curr = calc_patch_likelihood_function(L, patches_all_sampled, ...
        volume_curr, sigma2, rotations);
    toc

    likelihood_func_l_rot = pI_curr .* permute(rho_init, [3, 1, 2]);
    pl_rot_curr = likelihood_func_l_rot ./ sum(likelihood_func_l_rot, ...
        [1, 2, 3]);
    clear likelihood_func_l_rot;
    pl_rot_curr_reshaped = reshape(pl_rot_curr, 1, 1, K, (2 * L) ^ 2, Nd);

    %% Update x
    x_updated = (1 / Nd) * sum(pl_rot_curr_reshaped ...
        .* ls_A_dagger_patch, [3, 4, 5]);

    x_curr = cell_vec(x_updated, x_curr, ell_max, s_lens);
end
end
