function [x_curr, rho_curr] = EM_3d_parallel_shifts_sampling(patches, x_init, rho_init, ...
    L, rotations, sigma2, jball, Psilms, ell_max, s_lens, n_list, lms_list, vol_true, Ls, Ls_sampled, psipsi_cell, Ipsis_cell)

length_x = length(vec_cell(x_init));
Nd = size(patches, 1);
K = size(rotations, 3);

x_curr = x_init;
rho_curr = rho_init;

log_likelihood_prev = -inf;
count = 0;

NUM_ITERS = 100;

errs = zeros(NUM_ITERS, 1);
log_likelihoods = zeros(NUM_ITERS, 1);
cond_nums = zeros(NUM_ITERS, 1);

sampling_factor = size(Ls_sampled, 1);

while true
    count = count + 1;
    volume_curr = real(icfftn(expand_vol_psilms(x_curr, floor(L / 2), ...
        jball, Psilms, L)));
    save("volume_curr_"+num2str(count)+".mat", "volume_curr");
    [err, ~, volume_curr_rotated] = ...
        calc_volume_relative_error_rotations(vol_true, volume_curr, 20);
    errs(count) = err;
    Ss = cell(sampling_factor, 1);
    parfor l=1:sampling_factor
        S = zeros(K, Nd);
        for k=1:K
            S(k, :) = pI_l_rot_x(patches, Ls_sampled(l, :), ...
                rotations( :, :, k), volume_curr);
        end
        Ss{l} = S;
    end
    S = zeros(K, 2*L, 2*L, Nd);
    for l=1:sampling_factor
        S(:, Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1, :) = Ss{l};
    end
    S_reduced = zeros(K, sampling_factor, Nd);
    for ll=1:sampling_factor
        S_reduced(:, ll, :) = squeeze(S( :, Ls_sampled(ll, 1) + 1, Ls_sampled(ll, 2) + 1, :));
    end
    S_normalized = S - permute(min(S_reduced, [], [1, 2]), [1, 2, 4, 3]);
    S_normalized = max(S_normalized, 0);
    pI_curr = exp(-S_normalized / (2 * sigma2));
    for p=1:Nd
        sum_p = 0;
        for ll=1:sampling_factor
            sum_p = sum_p + sum(pI_curr(:, Ls_sampled(ll, 1) + 1, Ls_sampled(ll, 2) + 1, p));
        end
        pI_curr(:, :, :, p) = pI_curr(:, :, :, p) / sum_p;
    end
    likelihood_func_l_rot = pI_curr .* permute(rho_curr, [3, 1, 2]);
    pl_rot_curr = likelihood_func_l_rot ./ sum(likelihood_func_l_rot, ...
        [1, 2, 3]);
    pl_rot_curr_cell_1 = cell(sampling_factor, 1);
    pl_rot_curr_cell_2 = cell(sampling_factor, 1);
    pl_rot_curr_cell_3 = cell(sampling_factor, 1);

    for l=1:sampling_factor
        pl_rot_curr_cell_1{l} = squeeze(pl_rot_curr(1:50, Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1, :));
        pl_rot_curr_cell_2{l} = squeeze(pl_rot_curr(51:101, Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1, :));
        pl_rot_curr_cell_3{l} = squeeze(pl_rot_curr(102:152, Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1, :));


    end


    pI_curr_for_likelihood = exp(-S / (2*sigma2));
    log_likelihood = sum(squeeze(log(sum(pI_curr_for_likelihood .* ...
        permute(rho_curr, [3, 1, 2]), [1, 2, 3]))));

    log_likelihood
    log_likelihoods(count) = log_likelihood;

    rho_updated = rho_step(pl_rot_curr, Nd);


    tic;
    [x_updated, A, y] = x_step_parallel_all(Ls_sampled, rotations, ell_max, ...
        s_lens, n_list, length_x, lms_list, psipsi_cell, Ipsis_cell, pl_rot_curr_cell_1, pl_rot_curr_cell_2, pl_rot_curr_cell_3);
    toc
    cond_num = cond(A);
    cond_nums(count) = cond_num;

    x_curr = cell_vec(x_updated, x_curr, ell_max, s_lens);
    rho_curr = rho_updated;

    if count > NUM_ITERS
        break;
    end

    if log_likelihood < log_likelihood_prev
        'ERROR'
    end

    log_likelihood_prev = log_likelihood;
end
end