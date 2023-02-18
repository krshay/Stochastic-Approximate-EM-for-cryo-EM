function psi_lNs = calc_psi_lNs(ell_max, psi_Nn, betas, L, s_lens, n_list)
psi_lNs = cell(ell_max+1,1);
for ell=0:ell_max
    psi_lNs{ell+1} = cell(2*ell+1, 1);
    for N = -ell:ell
        psi_lNs{ell+1}{N+ell+1} = psi_Nn{N+ell_max+1} .* ...
            permute(betas{ell+1}{abs(N)+1}(1:s_lens(ell + 1), :), [3, 2, 1]);
        if N < 0
            psi_lNs{ell+1}{N+ell+1} = (-1)^(N) * psi_lNs{ell+1}{N+ell+1};
        end
        psi_lNs{ell+1}{N+ell+1} = reshape(psi_lNs{ell+1}{N+ell+1}, ...
            L, L, n_list(abs(N) + 1), s_lens(ell + 1));
        psi_lNs{ell+1}{N+ell+1} = icfft2(psi_lNs{ell+1}{N+ell+1});
    end
end
end