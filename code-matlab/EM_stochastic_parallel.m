function [x_curr, rho_curr] = EM_stochastic_parallel(patches_all, ...
    x_init, rho_init, L, rotations, sigma2, jball, Psilms, ell_max, ...
    s_lens, n_list, lms_list, Ls, psi_lsNn, STOCHASTIC_FACTOR, ...
    gs, NUM_ITERS)
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
rho_curr = rho_init;

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

    % Select a random subset of patches to use in the current iteration
    SAMPLES = randsample(NUMBER_OF_PATCHES, ...
        NUMBER_OF_PATCHES_PER_ITERATION);
    patches_all_sampled = patches_all(SAMPLES, :, :);
    Nd = NUMBER_OF_PATCHES_PER_ITERATION;

    % % Preprocess the I * psi values for the sampled patches
    % Ipsis_cell_sampled = preprocess_Ipsi(L, Ls, psi_lsNn, ell_max, ...
    %     s_lens, n_list, patches_all_sampled);
    %%
    tic
    % Calculate the probability function
    Ss = cell(num_shifts, 1);
    for l=1:num_shifts
        S = zeros(K, Nd);
        for k=1:K
            S(k, :) = pI_l_rot_x(patches_all_sampled, Ls(l, :), ...
                rotations( :, :, k), volume_curr);
        end
        Ss{l} = S;
    end

    S = zeros(K, 2*L, 2*L, Nd);
    for l=1:num_shifts
        S(:, Ls(l, 1) + 1, Ls(l, 2) + 1, :) = Ss{l};
    end
    % clear Ss

    % Reduce the size of S to only include the non-zero elements, and
    % normalize the values in S to be between 0 and 1
    S_reduced = zeros(K, num_shifts, Nd);
    for ll=1:num_shifts
        S_reduced(:, ll, :) = squeeze(S( :, Ls(ll, 1) + 1, ...
            Ls(ll, 2) + 1, :));
    end
    S_reduced(S_reduced == 0) = NaN;
    S_normalized = S - permute(min(S_reduced, [], [1, 2]), [1, 2, 4, 3]);
    % clear S
    % clear S_reduced

    S_normalized = max(S_normalized, 0);
    pI_curr_old = exp(-S_normalized / (2 * sigma2));
    % clear S_normalized
    pI_curr_old = permute(pI_curr_old, [4, 2, 3, 1]);
    pI_curr_old = pI_curr_old ./ sum(pI_curr_old, [2, 3, 4]);
    pI_curr_old = permute(pI_curr_old, [4, 2, 3, 1]);
    % for p=1:Nd
    %     sum_p = 0;
    %     for ll=1:num_shifts
    %         sum_p_ll = sum(pI_curr_old(:, Ls(ll, 1) + 1, Ls(ll, 2) + 1, p));
    %         sum_p = sum_p + sum_p_ll;
    %
    %     end
    %     pI_curr_old(:, :, :, p) = pI_curr_old(:, :, :, p) / sum_p;
    % end
    toc

    %%
    % tic
    % norm_value = zeros(NUMBER_OF_PATCHES_PER_ITERATION, num_shifts, K);
    % 
    % Mask = zeros(2*L, 2*L);
    % Mask(1:L, 1:L) = 1;
    % Mask_k = fft2(Mask);
    % Mask_k_k = fft2(Mask_k);
    % Z_patches_k = fft2(permute(patches_all_sampled, [2, 3, 1]), 2*L, 2*L);
    % 
    % precompute_MZ_patches = Mask_k_k .* ...
    %     fft2(Z_patches_k);
    % [u, v] = meshgrid(0:2*L-1, 0:2*L-1);
    % Ls_u = Ls(:, 1);
    % Ls_v = Ls(:, 2);
    % Pfrots = vol_project(volume_curr, rotations);
    % Z_Pfrots_k = fft2(Pfrots, 2*L, 2*L);
    % for l=1:num_shifts
    %     coeffs = exp(-1j * 2 * pi ...
    %         * (u * Ls_u(l)...
    %         + v * Ls_v(l)) / (2*L));
    %     norm_value( :, l, :) = squeeze(pagenorm( ...
    %         precompute_MZ_patches - Mask_k_k .* fft2(permute(coeffs .* Z_Pfrots_k, [1, 2, 4, 3]) ...
    %         ), "fro"));
    % end
    % norm_squared_value_old = norm_value.^2 / ((2 * L)^2 * num_shifts^3);
    % 
    % norm_squared_value_old_normalized = norm_squared_value_old - min(norm_squared_value_old, [], [2, 3]);
    % % pI_curr = exp(-norm_squared_value_normalized / (2 * sigma2));
    % %
    % % pI_curr = pI_curr ./ sum(pI_curr, [2, 3]);
    % % pI_curr = permute(reshape(pI_curr, [NUMBER_OF_PATCHES_PER_ITERATION, 2 * L, 2 * L, K]), [4, 3, 2, 1]);
    % toc

    tic

    Mask = zeros(2*L, 2*L);
    Mask(1:L, 1:L) = 1;
    Mask_k = fft2(Mask);

    patches_norm_squared = (squeeze(pagenorm(permute(patches_all_sampled, [2, 3, 1])))).^2;
    Z_patches_k = fft2(permute(patches_all_sampled, [2, 3, 1]), 2*L, 2*L);

    Pfrots = vol_project(volume_curr, rotations);
    Z_Pfrots_k = fft2(Pfrots, 2*L, 2*L);
    Z_Pfrots2_k = fft2(Pfrots.^2, 2*L, 2*L);
    C = real(ifft2(conj(Z_Pfrots_k) .* permute(Z_patches_k, [1, 2, 4, 3])));
    patches_estimates_norm_squared = real(ifft2(Mask_k .* conj(Z_Pfrots2_k)));
    % tmps1 = zeros(22^2, 1);
    % tmps2 = zeros(22^2, NUMBER_OF_PATCHES_PER_ITERATION);
    % for l=1:22^2
    % [~, tmp] = pI_l_rot_x(patches_all_sampled, Ls(l, :), ...
    %             rotations( :, :, 1), volume_curr);
    % tmps1(l) = norm(squeeze(tmp), "fro")^2;
    % tmps2(l, :) = sum(tmp .* patches_all_sampled, [2, 3]);
    % end
    % tmps1 = reshape(tmps1, 22, 22);
    % est_1 = patches_estimates_norm_squared(:, :, 1);
    % tmps2 = reshape(tmps2, 22, 22, 289);
    % C1 = squeeze(C( :, :, 1, :));

    norm_squared_value = permute(patches_norm_squared, [2, 3, 4, 1]) - 2 * C + patches_estimates_norm_squared;
    norm_squared_value_normalized = norm_squared_value - min(norm_squared_value, [], [1, 2, 3]);
    pI_curr = exp(-norm_squared_value_normalized / (2 * sigma2));

    pI_curr = pI_curr ./ sum(pI_curr, [1, 2, 3]);
    pI_curr = permute(pI_curr, [3, 1, 2, 4]);
    toc

    likelihood_func_l_rot = pI_curr .* permute(rho_curr, [3, 1, 2]);
    pl_rot_curr = likelihood_func_l_rot ./ sum(likelihood_func_l_rot, ...
        [1, 2, 3]);
    clear likelihood_func_l_rot;
    %%
    % Update rho
    rho_updated = update_rho(pl_rot_curr, Nd);

    pl_rot_curr_cell = cell(num_shifts, 1);

    for l=1:num_shifts
        pl_rot_curr_cell{l} = squeeze(pl_rot_curr(:, Ls(l, 1) + 1, ...
            Ls(l, 2) + 1, :));
    end

    % Update x
    x_updated = update_x(Ls, rotations, ell_max, ...
        s_lens, length_x, lms_list, gs, Ipsis_cell_sampled, ...
        pl_rot_curr_cell);

    x_curr = cell_vec(x_updated, x_curr, ell_max, s_lens);
    rho_curr = rho_updated;
end
end