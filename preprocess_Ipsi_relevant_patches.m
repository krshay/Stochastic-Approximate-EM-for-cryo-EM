function Ipsis = preprocess_Ipsi_relevant_patches(L, Ls, psi_lsNn, ell_max, ...
    s_lens, n_list, relevant_patches, number_of_micrograph, patches_nonempty)

Ipsis = cell(size(Ls, 1), 1);
parfor l=1:size(Ls, 1) - 1
    Ipsi = zeros(size(relevant_patches{number_of_micrograph, l}, 1), ell_max+1, 2*ell_max + 1, max(s_lens));
    for ell=0:ell_max
        for s=1:s_lens(ell + 1)
            for N=-ell:ell
                Ipsi( :, ell + 1, N + ell + 1, s) ...
                    = sum(relevant_patches{number_of_micrograph, l} .* permute(CTZ(psi_lsNn{ell + 1} ...
                    {N + ell + 1}( :, :, 1:n_list(abs(N) + 1), s), ...
                    Ls(l, :), L), [4, 1, 2, 3]), ...
                    [2, 3, 4]);
            end
        end
    end
    Ipsis{l} = Ipsi;
end
l = size(Ls, 1);
    Ipsi = zeros(size(patches_nonempty{number_of_micrograph}, 1), ell_max+1, 2*ell_max + 1, max(s_lens));
    for ell=0:ell_max
        for s=1:s_lens(ell + 1)
            for N=-ell:ell
                Ipsi( :, ell + 1, N + ell + 1, s) ...
                    = sum(patches_nonempty{number_of_micrograph} .* permute(CTZ(psi_lsNn{ell + 1} ...
                    {N + ell + 1}( :, :, 1:n_list(abs(N) + 1), s), ...
                    Ls(l, :), L), [4, 1, 2, 3]), ...
                    [2, 3, 4]);
            end
        end
    end
    Ipsis{l} = Ipsi;
end