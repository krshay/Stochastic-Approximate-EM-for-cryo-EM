function [b_Nn, b_minusNn] = compute_bNn(ell_max, n_list, x_lms, betas, s_lens, omega)
b_Nn = zeros(ell_max+1, max(n_list));
b_minusNn = zeros(ell_max + 1, max(n_list));
for N=0:ell_max
    for n=0:n_list(N+1)-1
        for ell=N:ell_max
            D_ell = wignerd(ell, omega);
            b_Nn(N+1, n+1) = b_Nn(N+1, n+1) + ...
                (x_lms{ell+1}( :, :) * D_ell(N + ell + 1, :).').' ...
                * betas{ell+1}{N + 1}(1:s_lens(ell+1), n+1);
            b_minusNn(N+1, n+1) = b_minusNn(N+1, n+1) + ...
                (x_lms{ell+1}( :, :) * D_ell(-N + ell + 1, :).').' ...
                * (-1)^N * conj(betas{ell+1}{N + 1}(1:s_lens(ell+1), n+1));
        end
    end
end
end
