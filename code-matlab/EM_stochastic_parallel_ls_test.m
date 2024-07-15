function x_curr = EM_stochastic_parallel_ls_test(patches_all, ...
    x_init, rho_init, L, rotations, sigma2, jball, Psilms, ell_max, ...
    s_lens, n_list, lms_list, Ls, psi_lsNn, STOCHASTIC_FACTOR, ...
    gs, NUM_ITERS, PR_matrix_all_rotations_all_shifts, pl_rot_curr_reshaped)
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
num_shifts = size(Ls, 1);

% Start an infinite loop
while count <= NUM_ITERS
    % Increment the count
    count = count + 1;

    % Get the current volume by computing the inverse Fourier transform of
    % the expansion of x with respect to Psilms
    volume_curr = real(icfftn(expand_vol_psilms(x_curr, floor(L / 2), ...
        jball, Psilms, L)));

    save("volume_curr" + num2str(count-1) + ".mat", "volume_curr")

    % Select a random subset of patches to use in the current iteration
    SAMPLES = randsample(NUMBER_OF_PATCHES, ...
        NUMBER_OF_PATCHES_PER_ITERATION);
    patches_all_sampled = patches_all(SAMPLES, :, :);
    patches_all_sampled = patches_all;
    Nd = NUMBER_OF_PATCHES_PER_ITERATION;

    % %%
    % % TODO: Parallelize.
    % tic
    % Mask = zeros(2*L, 2*L);
    % Mask(1:L, 1:L) = 1;
    % Mask_k = fft2(Mask);
    %
    % patches_norm_squared = (squeeze(pagenorm(permute(patches_all_sampled, [2, 3, 1])))).^2;
    % Z_patches_k = fft2(permute(patches_all_sampled, [2, 3, 1]), 2*L, 2*L);
    %
    % Pfrots = vol_project(volume_curr, rotations);
    % Z_Pfrots_k = fft2(Pfrots, 2*L, 2*L);
    % Z_Pfrots2_k = fft2(Pfrots.^2, 2*L, 2*L);
    % C = real(ifft2(conj(Z_Pfrots_k) .* permute(Z_patches_k, [1, 2, 4, 3])));
    % patches_estimates_norm_squared = real(ifft2(Mask_k .* conj(Z_Pfrots2_k)));
    %
    % norm_squared_value = permute(patches_norm_squared, [2, 3, 4, 1]) - 2 * C + patches_estimates_norm_squared;
    % norm_squared_value_normalized = norm_squared_value - min(norm_squared_value, [], [1, 2, 3]);
    %
    % pI_curr = exp(-norm_squared_value_normalized / (2 * sigma2));
    % pI_curr = pI_curr ./ sum(pI_curr, [1, 2, 3]);
    % pI_curr = permute(pI_curr, [3, 1, 2, 4]);
    % toc
    %
    % likelihood_func_l_rot = pI_curr .* permute(rho_init, [3, 1, 2]);
    % pl_rot_curr = likelihood_func_l_rot ./ sum(likelihood_func_l_rot, ...
    %     [1, 2, 3]);
    % clear likelihood_func_l_rot;
    % pl_rot_curr_reshaped = reshape(pl_rot_curr, 1, 1, K, (2 * L) ^ 2, Nd);
    %% Update x
    patches_reshaped = permute(reshape(patches_all_sampled, Nd, L^2), [2, 3, 4, 5, 1]);
    % [ls_A_dagger_patch, rCond] = pagemldivide(PR_matrix_all_rotations_all_shifts, patches_reshaped);
    % %
    % x_updated1 = (1 / Nd) * sum(pl_rot_curr_reshaped .* ls_A_dagger_patch, [3, 4, 5]);
    % PINV_PR_matrix_all_rotations_all_shifts = permute(zeros(size(PR_matrix_all_rotations_all_shifts)), [2, 1, 3, 4]);
    Gamma = 1;
    % for s=1:size(Ls, 1)
    %     for k=1:K
    %         % [~, S, ~] = svd(PR_matrix_all_rotations_all_shifts( :, :, k, s));
    %         % lambda_max = max(S(:));
    %         A = PR_matrix_all_rotations_all_shifts( :, :, k, s);
    %         PINV_PR_matrix_all_rotations_all_shifts( :, :, k, s) = (A' * A + (Gamma ^ 2) * eye(size(PR_matrix_all_rotations_all_shifts, 2))) \ A';
    %         % PINV_PR_matrix_all_rotations_all_shifts( :, :, k, s) = pinv(PR_matrix_all_rotations_all_shifts( :, :, k, s), lambda_max / 100);
    %     end
    % end
    PINV_PR_matrix_all_rotations_all_shifts = pagemldivide(pagemtimes(pagectranspose(PR_matrix_all_rotations_all_shifts), PR_matrix_all_rotations_all_shifts) + (Gamma ^ 2) * eye(size(PR_matrix_all_rotations_all_shifts, 2)), pagectranspose(PR_matrix_all_rotations_all_shifts));
    ls_A_pinv_patch = pagemtimes(PINV_PR_matrix_all_rotations_all_shifts, patches_reshaped);
        x_updated = (1 / Nd) * sum(pl_rot_curr_reshaped .* ls_A_pinv_patch, [3, 4, 5]);
    x_curr = cell_vec(x_updated, x_curr, ell_max, s_lens);

    Ax = pagemtimes(PR_matrix_all_rotations_all_shifts, vec_cell(x_curr));
    norm_diff = sum(abs(patches_reshaped - Ax) .^ 2, [1, 2]);
    Q = sum(pl_rot_curr_reshaped .* norm_diff, "all");
    Q


    % Ax = pagemtimes(PR_matrix_all_rotations_all_shifts, vec_cell(x_curr));
    % norm_diff = sum(abs(patches_reshaped - Ax) .^ 2, [1, 2]);
    % Q = sum(pl_rot_curr_reshaped .* norm_diff, "all");
    % Q
end
end
