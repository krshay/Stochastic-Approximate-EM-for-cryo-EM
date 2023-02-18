function [x_updated, A, y] = x_step(pI_rot_curr, patches, Ls, rotations, L, ...
    ell_max, psi_lsNn, s_lens, n_list, length_x)

K = size(rotations, 3);

y = zeros(length_x, 1);
A = zeros(length_x, length_x);
for k=1:K
    [yaw, pitch, roll] = dcm2angle(rotations( :, :, k), 'ZYZ');
    roll = roll + pi/2;
    omega = [yaw, pitch, roll];
    for l=1:size(Ls, 1)
        lms = 0;
        for ell=0:ell_max
            D = wignerd(ell, omega);
            for m=-ell:ell
                for s=1:s_lens(ell + 1)
                    lms = lms + 1;
                    for N=-ell:ell
                        for n=1:n_list(abs(N) + 1)
                            y(lms) = y(lms) ...
                                + sum( ...
                                squeeze(pI_rot_curr(k, Ls(l, 1) + 1, ...
                                Ls(l, 2) + 1, :)) ...
                                .* sum( ...
                                patches * D(-N + ell + 1, -m + ell + 1) ...
                                .* permute(CTZ(psi_lsNn{ell + 1} ...
                                {N + ell + 1}( :, :, n, s), ...
                                Ls(l, :), L), [3, 1, 2]), ...
                                [2, 3]));
                            lms_tilde = 0;
                            for ell_tilde=0:ell_max
                                D_tilde = wignerd(ell_tilde, omega);
                                for m_tilde=-ell_tilde:ell_tilde
                                    for s_tilde=1:s_lens(ell_tilde + 1)
                                        lms_tilde = lms_tilde + 1;
                                        for N_tilde=-ell_tilde:ell_tilde
                                            for n_tilde=1:n_list(abs(N_tilde) + 1)
                                                A(lms, lms_tilde) = A(lms, lms_tilde) + sum(pI_rot_curr(k, Ls(l, 1) + 1, Ls(l, 2) + 1, :), 'all') ...
                                                    * sum(D(-N + ell + 1, -m + ell + 1) * CTZ(psi_lsNn{ell + 1}{N + ell + 1}( :, :, n, s), Ls(l, :), L) ...
                                                    * D_tilde(-N_tilde + ell_tilde + 1, -m_tilde + ell_tilde + 1) .* CTZ(psi_lsNn{ell_tilde + 1}{N_tilde + ell_tilde + 1}( :, :, n_tilde, s_tilde), Ls(l, :), L), 'all');
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
    end
end
x_updated = A \ y;
end