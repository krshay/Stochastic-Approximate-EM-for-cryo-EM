function [x_curr, rho_curr] = EM_3d(patches, x_init, rho_init, ...
    L, rotations, sigma2, jball, Psilms, ell_max, psi_lsNn, s_lens, n_list, psipsi, lms_list, vol_true)

length_x = length(vec_cell(x_init));
Nd = size(patches, 1);
K = size(rotations, 3);
Ls = calc_shifts(L);
x_curr = x_init;
rho_curr = rho_init;

log_likelihood_prev = -inf;
count = 0;

NUM_ITERS = 50;

errs = zeros(NUM_ITERS, 1);
log_likelihoods = zeros(NUM_ITERS, 1);
cond_nums = zeros(NUM_ITERS, 1);


S = zeros(K, 2*L, 2*L, Nd);
while true
    count = count + 1;
    volume_curr = real(icfftn(expand_vol_psilms(x_curr, floor(L / 2), ...
        jball, Psilms, L)));
    err = calc_volume_relative_error_rotations(vol_true, volume_curr, 20);
    errs(count) = err;
    tic
    for k=1:K
        for l=1:size(Ls, 1)
            S(k, Ls(l, 1) + 1, Ls(l, 2) + 1, :) = pI_l_rot_x(patches, ...
                Ls(l, :), rotations( :, :, k), volume_curr);
        end
    end
    toc
    S_normalized = S - min(S, [], [1, 2, 3]);
    pI_curr = exp(-S_normalized / (2 * sigma2));
    pI_curr = pI_curr ./ sum(pI_curr, [1, 2, 3]);
    likelihood_func_l_rot  = pI_curr .* permute(rho_curr, [3, 1, 2]);
    pl_rot_curr = likelihood_func_l_rot ./ sum(likelihood_func_l_rot, ...
        [1, 2, 3]);

    pI_curr_for_likelihood = exp(-S / (2*sigma2));
    log_likelihood = sum(squeeze(log(sum(pI_curr_for_likelihood .* ...
        permute(rho_curr, [3, 1, 2]), [1, 2, 3]))));
    log_likelihoods(count) = log_likelihood;

    rho_updated = rho_step(pl_rot_curr, Nd);
%     tic;
%     [x_check, A_check, y_check] = x_step(pl_rot_curr, patches, Ls, rotations, L, ell_max, ...
%         psi_lsNn, s_lens, n_list, length_x);
%     toc;

    tic;
    [x_updated, A, y] = x_step_accelerated(pl_rot_curr, patches, Ls, rotations, L, ell_max, ...
        psi_lsNn, s_lens, n_list, length_x, psipsi, lms_list);
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