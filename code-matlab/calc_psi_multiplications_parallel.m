function psipsi = calc_psi_multiplications_parallel(psi_lsNn, L, ell_max, n_list, s_lens, Ls)

psipsis = cell(size(Ls, 1), 1);
parfor l=1:size(Ls, 1)
    psipsi = zeros(ell_max + 1, 2*ell_max+1, ell_max + 1, 2*ell_max+1, ...
        max(n_list), max(s_lens), max(n_list), max(s_lens));
    for ell=0:ell_max
        for s=1:s_lens(ell + 1)
            for N=-ell:ell
                for n=1:n_list(abs(N) + 1)
                    for ell_tilde=0:ell_max
                        for s_tilde=1:s_lens(ell_tilde + 1)
                            for N_tilde=-ell_tilde:ell_tilde
                                for n_tilde=1:n_list(abs(N_tilde) + 1)
                                    psipsi(ell+1, N+ell+1, ell_tilde+1, N_tilde+ell_tilde+1, n, s, n_tilde, s_tilde) ...
                                        = sum(CTZ(psi_lsNn{ell + 1}{N + ell + 1}( :, :, n, s), Ls(l, :), L) ...
                                        .* CTZ(psi_lsNn{ell_tilde + 1}{N_tilde + ell_tilde + 1}( :, :, n_tilde,s_tilde), Ls(l, :), L), "all");
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    psipsis{l} = psipsi;
end
psipsi = zeros(size(Ls, 1), ell_max + 1, 2*ell_max+1, ell_max + 1, 2*ell_max+1, ...
    max(n_list), max(s_lens), max(n_list), max(s_lens));
for l=1:size(Ls, 1)
    psipsi(l, :, :, :, :, :, :, :, :) = psipsis{l};
end
end