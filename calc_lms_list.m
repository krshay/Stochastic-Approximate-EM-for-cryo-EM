function lms_list = calc_lms_list(ell_max, s_lens)
lms_list = [];
lms = 0;
for ell=0:ell_max
    for m=-ell:ell
        for s=1:s_lens(ell + 1)
            lms = lms + 1;
            lms_list = [lms_list; ell, m, s];
        end
    end
end
end