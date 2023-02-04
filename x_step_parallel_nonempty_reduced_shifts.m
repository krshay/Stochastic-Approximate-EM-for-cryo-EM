function [x_updated, A, y] = x_step_parallel_nonempty_reduced_shifts(Ls, rotations, ...
    ell_max, s_lens, n_list, length_x, lms_list, psipsi_cell, number_in_series, ...
    Ipsis_cell_all, pl_rot_curr_total)

K = size(rotations, 3);

As = cell(size(Ls, 1), 1);
ys = cell(size(Ls, 1), 1);


[yaw, pitch, roll] = dcm2angle(rotations, 'ZYZ');
roll = roll + pi/2;
omegas = [yaw, pitch, roll];

for nn=1:number_in_series
    pl_rot_curr = pl_rot_curr_total{nn};
    Ipsis_cell = Ipsis_cell_all{nn};
    parfor l=1:size(Ls, 1)
        psipsi_l = psipsi_cell{l};
        Ipsi_l = Ipsis_cell{l};
        A = zeros(length_x, length_x);
        y = zeros(length_x, 1);

        pI_l = pl_rot_curr{l};
        sum_pI_l = sum(pI_l, 2);
        for ell=0:ell_max
            D = calc_wignerd(ell, omegas);
            for m=0:ell
                for s=1:s_lens(ell + 1)
                    lms = find(sum(lms_list == [ell, m, s], 2) == 3);
                    y(lms) = y(lms) + sum(squeeze(sum(squeeze(Ipsi_l(:, ell + 1, 1:2 * ell + 1, s)) .* permute(flip(D(:, -m + ell + 1, :), 1), [2, 1, 3]), 2)) .* pI_l.', "all");
                    for N=-ell:ell
                        for n=1:n_list(abs(N) + 1)
                            for ell_tilde=0:ell
                                D_tilde = calc_wignerd(ell_tilde, omegas);
                                for m_tilde=-ell_tilde:ell_tilde
                                    for s_tilde=1:s_lens(ell_tilde + 1)
                                        lms_tilde = find(sum(lms_list == [ell_tilde, m_tilde, s_tilde], 2) == 3);
                                        for N_tilde=-ell_tilde:ell_tilde
                                            for n_tilde=1:n_list(abs(N_tilde) + 1)
                                                A(lms, lms_tilde) = A(lms, lms_tilde) + (squeeze(D(-N + ell + 1, -m + ell + 1, :)) .* squeeze(D_tilde(-N_tilde + ell_tilde + 1, -m_tilde + ell_tilde + 1, :))).' * sum_pI_l * psipsi_l(ell+1, N+ell+1, ell_tilde+1, N_tilde+ell_tilde+1, n, s, n_tilde, s_tilde);
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
        if nn > 1
            ys{l} = ys{l} + y;
            As{l} = As{l} + A;
        else
            ys{l} = y;
            As{l} = A;
        end
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