function D = calc_wignerd(ell, omegas)
K = size(omegas, 1);
D = zeros(2 * ell + 1, 2 * ell + 1, K);
for k=1:K
    D( :, :, k) = wignerd(ell, omegas(k, :));
end
end