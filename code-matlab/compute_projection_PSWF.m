function proj = compute_projection_PSWF(x_lms, ell_max, s_lens, psi_lsNn, ...
    omega)
proj = 0;
for ell=0:ell_max
    D = wignerd(ell, omega);
    for m=-ell:ell
            for N=-ell:ell
                proj = proj + sum(permute(x_lms{ell + 1}(1:s_lens(ell+1), m + ell + 1), [3, 4, 2, 1]) ...
                * D(-N + ell + 1, -m + ell + 1) ...
                .* psi_lsNn{ell + 1}{N + ell + 1}, [3, 4]);
            end
    end
end
proj = real(proj);
end