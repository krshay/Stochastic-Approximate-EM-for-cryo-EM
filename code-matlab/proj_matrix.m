function A = proj_matrix(L, psi_lsNn, lms_list, n_list, rotations)
NUM_ROTATIONS = size(rotations, 3);
A = zeros(L^2, size(lms_list, 1), NUM_ROTATIONS);

[yaw, pitch, roll] = dcm2angle(rotations, 'ZYZ');
roll = roll + pi/2;
omegas = [yaw, pitch, roll];

for k=1:size(lms_list, 1)
    A_k = zeros(L, L, NUM_ROTATIONS);
    ell = lms_list(k, 1);
    D = calc_wignerd(ell, omegas);

    m = lms_list(k, 2);
    s = lms_list(k, 3);
    for N=-ell:ell
        D_N_m_ell = D(-N + ell + 1, -m + ell + 1, :);
        tmp = sum(psi_lsNn{ell + 1}{N + ell + 1}( :, :, 1:n_list(abs(N) + 1), s), 3);
        A_k = A_k + D_N_m_ell .* tmp;
    end
    A( :, k, :) = reshape(A_k, [L^2, NUM_ROTATIONS]);
end

end