function [x_curr, rho_curr] = EM_3d_parallel_nonempty_reduced_shifts(patches_nonempty, relevant_patches, x_init, rho_init, ...
    L, rotations, sigma2, jball, Psilms, ell_max, s_lens, n_list, lms_list, ...
    vol_true, Ls, Ls_sampled, psipsi_cell, number_in_series, Ipsis_cell_all, meaningful_shifts, nonempty_patches)


length_x = length(vec_cell(x_init));
K = size(rotations, 3);

x_curr = x_init;
rho_curr = rho_init;

log_likelihood_prev = -inf;
count = 0;

NUM_ITERS = 100;

log_likelihoods = zeros(NUM_ITERS, 1);
cond_nums = zeros(NUM_ITERS, 1);

sampling_factor = size(Ls_sampled, 1);

while true
    rho_temp = zeros(size(rho_init));
    count = count + 1;
    volume_curr = real(icfftn(expand_vol_psilms(x_curr, floor(L / 2), ...
        jball, Psilms, L)));
    save("volume_curr_"+num2str(count)+".mat", "volume_curr");

    pl_rot_curr_cell_total = cell(number_in_series, 1);

    for nn=1:number_in_series
        Nd = size(patches_nonempty{nn}, 1);

        Ss = cell(sampling_factor, 1);
        parfor l=1:sampling_factor-1
            S = zeros(K, Nd);
            for k=1:K
                S(k, meaningful_shifts{nn, l}) = pI_l_rot_x(relevant_patches{nn, l}, Ls_sampled(l, :), ...
                    rotations( :, :, k), volume_curr);
            end
            Ss{l} = S;
        end
        l = sampling_factor;
        S = zeros(K, Nd);
        for k=1:K
            S(k, :) = pI_l_rot_x(patches_nonempty{nn}, Ls_sampled(l, :), ...
                rotations( :, :, k), volume_curr);
        end
        Ss{l} = S;
        S = zeros(K, 2*L, 2*L, Nd);
        for l=1:sampling_factor
            S(:, Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1, :) = Ss{l};
        end
        clear Ss
        S_reduced = zeros(K, sampling_factor, Nd);
        for ll=1:sampling_factor
            S_reduced(:, ll, :) = squeeze(S( :, Ls_sampled(ll, 1) + 1, Ls_sampled(ll, 2) + 1, :));
        end
        S_reduced(S_reduced == 0) = NaN;
        S_normalized = S - permute(min(S_reduced, [], [1, 2]), [1, 2, 4, 3]);
        clear S_reduced
        clear S
        S_normalized = max(S_normalized, 0);
        pI_curr = exp(-S_normalized / (2 * sigma2));
        clear S_normalized
        for p=1:Nd
            sum_p = 0;
            for ll=1:sampling_factor
                sum_p_ll = sum(pI_curr(:, Ls_sampled(ll, 1) + 1, Ls_sampled(ll, 2) + 1, p));
                if sum_p_ll ~= K
                    sum_p = sum_p + sum_p_ll;
                else
                    pI_curr(:, Ls_sampled(ll, 1) + 1, Ls_sampled(ll, 2) + 1, p) = 0;
                end
            end
            pI_curr(:, :, :, p) = pI_curr(:, :, :, p) / sum_p;
        end
        likelihood_func_l_rot = pI_curr .* permute(rho_curr, [3, 1, 2]);
        pl_rot_curr = likelihood_func_l_rot ./ sum(likelihood_func_l_rot, ...
            [1, 2, 3]);

        rho_temp = rho_temp + rho_step(pl_rot_curr, Nd)/number_in_series;

        pl_rot_curr_cell = cell(sampling_factor, 1);

        for l=1:sampling_factor-1
            pl_rot_curr_cell{l} = squeeze(pl_rot_curr(:, Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1, meaningful_shifts{nn, l}));
        end
        l = sampling_factor;
        pl_rot_curr_cell{l} = squeeze(pl_rot_curr(:, Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1, :));

        pl_rot_curr_cell_total{nn} = pl_rot_curr_cell;
    end

    rho_updated = rho_temp;

    tic;
    [x_updated, A, y] = x_step_parallel_nonempty_reduced_shifts(Ls_sampled, rotations, ell_max, ...
        s_lens, n_list, length_x, lms_list, psipsi_cell, number_in_series, Ipsis_cell_all, pl_rot_curr_cell_total);
    toc
    cond_num = cond(A);
    cond_nums(count) = cond_num;

    x_curr = cell_vec(x_updated, x_curr, ell_max, s_lens);
    rho_curr = rho_updated;

    if count > NUM_ITERS
        break;
    end

end
end