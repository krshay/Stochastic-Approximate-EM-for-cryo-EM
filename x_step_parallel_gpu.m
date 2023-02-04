function [x_updated, A, y] = x_step_parallel_gpu(pI_rot_curr, patches, Ls, rotations, L, ...
    ell_max, psi_lsNn, s_lens, n_list, length_x, psipsi, lms_list)

K = size(rotations, 3);


As = cell(size(Ls, 1), 1);
ys = cell(size(Ls, 1), 1);
patches = gpuArray(patches);
psipsi = gpuArray(psipsi);
parfor l=1:size(Ls, 1)
    A = zeros(length_x, length_x);
    y = zeros(length_x, 1);
    for k=1:K
        [yaw, pitch, roll] = dcm2angle(rotations( :, :, k), 'ZYZ');
        roll = roll + pi/2;
        omega = [yaw, pitch, roll];
        pI_l_k = pI_rot_curr(k, Ls(l, 1) + 1, Ls(l, 2) + 1, :);
        sum_pI_l_k = sum(pI_l_k, 'all');
        for ell=0:ell_max
            D = wignerd(ell, omega);
            for m=0:ell
                for s=1:s_lens(ell + 1)
                    lms = find(sum(lms_list == [ell, m, s], 2) == 3);
                    for N=-ell:ell
                        for n=1:n_list(abs(N) + 1)
                            y(lms) = y(lms) ...
                                + sum( ...
                                squeeze(pI_l_k) ...
                                .* sum( ...
                                patches * D(-N + ell + 1, -m + ell + 1) ...
                                .* permute(CTZ(psi_lsNn{ell + 1} ...
                                {N + ell + 1}( :, :, n, s), ...
                                Ls(l, :), L), [3, 1, 2]), ...
                                [2, 3]));
                            for ell_tilde=0:ell
                                D_tilde = wignerd(ell_tilde, omega);
                                for m_tilde=-ell_tilde:ell_tilde
                                    for s_tilde=1:s_lens(ell_tilde + 1)
                                        lms_tilde = find(sum(lms_list == [ell_tilde, m_tilde, s_tilde], 2) == 3);
                                        for N_tilde=-ell_tilde:ell_tilde
                                            for n_tilde=1:n_list(abs(N_tilde) + 1)
                                                A(lms, lms_tilde) = A(lms, lms_tilde) + sum_pI_l_k ...
                                                    * D(-N + ell + 1, -m + ell + 1) * D_tilde(-N_tilde + ell_tilde + 1, -m_tilde + ell_tilde + 1) * psipsi(l, ell+1, N+ell+1, ell_tilde+1, N_tilde+ell_tilde+1, n, s, n_tilde, s_tilde);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        ys{l} = y;
        As{l} = A;
    end
end
A = zeros(length_x, length_x);
y = zeros(length_x, 1);
for l=1:size(Ls, 1)
    A = A + As{l};
    y = y + ys{l};
end
for ell=0:ell_max
    for m=-ell:0
        for s=1:s_lens(ell + 1)
            y(sum(lms_list == [ell, m, s], 2) == 3) = (-1)^(ell - m) * conj(y(lms_list(:, 1) == ell & lms_list(:, 2) == -m & lms_list(:, 3) == s));
            for ell_tilde=0:ell
                for m_tilde=-ell_tilde:ell_tilde
                    for s_tilde=1:s_lens(ell_tilde + 1)
                        A(sum(lms_list == [ell, m, s], 2) == 3, sum(lms_list == [ell_tilde, m_tilde, s_tilde], 2) == 3) = (-1)^(ell - m) * (-1)^(ell_tilde - m_tilde) * conj(A(lms_list(:, 1) == ell & lms_list(:, 2) == -m & lms_list(:, 3) == s, lms_list(:, 1) == ell_tilde & lms_list(:, 2) == -m_tilde & lms_list(:, 3) == s_tilde));
                    end
                end
            end
        end
    end
end
A = tril(A) + (tril(A, -1)).';
x_updated = A \ y;
end