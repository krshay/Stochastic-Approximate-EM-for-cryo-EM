function x = fill_negative_ms(x_positive_ms, ell_max, s_lens)
x = cell(ell_max + 1, 1);
idx_positive_ms = 0;
for ell=0:ell_max
    x{ell + 1} = zeros(s_lens(ell + 1), 2 * ell + 1);
    for m=0:ell
        for s=1:s_lens(ell + 1)
            idx_positive_ms = idx_positive_ms + 1;
            x{ell + 1}(s, m + ell + 1) = x_positive_ms(idx_positive_ms);
            x{ell + 1}(s, -m + ell + 1) = (-1)^(ell + m) * conj(x_positive_ms(idx_positive_ms));
        end
    end
end


end
