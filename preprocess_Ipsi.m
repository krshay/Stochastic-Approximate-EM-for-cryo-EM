function Ipsis = preprocess_Ipsi(L, Ls, psi_lsNn, ell_max, s_lens, n_list, patches)

Ipsis = cell(size(Ls, 1), 1);
parfor l=1:size(Ls, 1)
    Ipsi = zeros(size(patches, 1), ell_max+1, 2*ell_max + 1, max(s_lens));
    for ell=0:ell_max
        for s=1:s_lens(ell + 1)
            for N=-ell:ell
                Ipsi( :, ell + 1, N + ell + 1, s) ...
                    = sum(patches .* permute(CTZ(psi_lsNn{ell + 1} ...
                    {N + ell + 1}( :, :, 1:n_list(abs(N) + 1), s), ...
                    Ls(l, :), L), [4, 1, 2, 3]), ...
                    [2, 3, 4]);
            end
        end
    end
    Ipsis{l} = Ipsi;
end
% Ipsi = zeros(size(Ls, 1), size(patches, 1), ell_max+1, 2*ell_max + 1, max(s_lens));
% for l=1:size(Ls, 1)
%     Ipsi_l = Ipsis{l};
%     save("/data/shaykreymer/vonneuman/Ipsi_"+num2str(l)+".mat", "Ipsi_l")
% end
end