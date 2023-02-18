function [x_updated, A, y] = update_x(Ls, rotations, ...
    ell_max, s_lens, length_x, lms_list, gs, ...
    Ipsis_cell, pl_rot_curr)

As = cell(size(Ls, 1), 1);
ys = cell(size(Ls, 1), 1);


[yaw, pitch, roll] = dcm2angle(rotations, 'ZYZ');
roll = roll + pi/2;
omegas = [yaw, pitch, roll];

num_shifts = size(Ls, 1);


parfor l=1:num_shifts
    g = gs{l};
    Ipsi_l = Ipsis_cell{l};
    pI_l = pl_rot_curr{l};

    pI_l_sum = sum(pI_l, 2);

    A = squeeze(sum(g .* pI_l_sum, 1));
    y = zeros(length_x, 1);

    for lms=1:size(lms_list, 1)
        ell = lms_list(lms, 1);
        m = lms_list(lms, 2);
        s = lms_list(lms, 3);
        if m >= 0
            D = calc_wignerd(ell, omegas);
            y(lms) = sum(squeeze(sum(squeeze(Ipsi_l(:, ell + 1, 1:2 * ell + 1, s)) .* permute(flip(D(:, -m + ell + 1, :), 1), [2, 1, 3]), 2)) .* pI_l.', "all");
        end
    end

    ys{l} = y;
    As{l} = A;
end

A = zeros(length_x, length_x);
y = zeros(length_x, 1);
for l=1:num_shifts
    A = A + As{l};
    y = y + ys{l};
end
for ell=0:ell_max
    for m=-ell:-1
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