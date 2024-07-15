function [x_lms, vol_rec] = expand_vol_spherical_basis(vol, N, ell_max, ...
    L, Psilms, jball)

x_lms = evaluate_xlms(vol, ell_max, Psilms, jball);
vol_rec = expand_vol_psilms(x_lms, N, jball, Psilms, L);
end
