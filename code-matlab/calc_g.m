function gs = calc_g(Ls, L, rotations, psi_lsNn, lms_list, n_list)
% This function calculates the g term of Eq. (35).
K = size(rotations, 3);
[yaw, pitch, roll] = dcm2angle(rotations, 'ZYZ');
roll = roll + pi/2;
omegas = [yaw, pitch, roll];

num_shifts = size(Ls, 1);

gs = cell(num_shifts, 1);
parfor l=1:num_shifts
    l
    g = zeros(K, size(lms_list, 1), size(lms_list, 1));
    for lms=1:size(lms_list, 1)
        ell = lms_list(lms, 1);

        m = lms_list(lms, 2);
        if m >= 0
            s = lms_list(lms, 3);

            D = calc_wignerd(ell, omegas);

            for lms_tilde = 1:size(lms_list, 1)
                ell_tilde = lms_list(lms_tilde, 1);
                if ell_tilde <= ell
                    m_tilde = lms_list(lms_tilde, 2);
                    s_tilde = lms_list(lms_tilde, 3);
                    D_tilde = calc_wignerd(ell_tilde, omegas);
                    for N=-ell:ell
                        for n=1:n_list(abs(N) + 1)
                            for N_tilde=-ell_tilde:ell_tilde
                                for n_tilde=1:n_list(abs(N_tilde) + 1)
                                    g( :, lms, lms_tilde)  = g( :, lms, lms_tilde) + (squeeze(D(-N + ell + 1, -m + ell + 1, :)) .* squeeze(D_tilde(-N_tilde + ell_tilde + 1, -m_tilde + ell_tilde + 1, :))) * ...
                                        sum(CTZ(psi_lsNn{ell + 1}{N + ell + 1}( :, :, n, s), Ls(l, :), L) ...
                                        .* CTZ(psi_lsNn{ell_tilde + 1}{N_tilde + ell_tilde + 1}( :, :, n_tilde,s_tilde), Ls(l, :), L), "all");
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    gs{l} = g;
end
end
