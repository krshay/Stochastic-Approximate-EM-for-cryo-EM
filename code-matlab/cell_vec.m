function x_cell = cell_vec(x, x_cell, ell_max, s_lens)
lms = 0;
for ell=0:ell_max
    for m=-ell:ell
        for s=1:s_lens(ell + 1)
            lms = lms + 1;
            x_cell{ell + 1}(s, m + ell + 1) = x(lms);
        end
    end
end
end