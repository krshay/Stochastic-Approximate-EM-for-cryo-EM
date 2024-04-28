function A = proj_matrix(L, psi_lsNn, lms_list, n_list, rotations)
NUM_ROTATIONS = size(rotations, 3);
A = zeros(L, L, size(lms_list, 1), NUM_ROTATIONS);

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
    % A_k_zero_pad = zeros(2 * L, 2 * L, size(A_k, 3));
    % A_k_zero_pad(1:L, 1:L, :) = A_k;
    % A_k_zero_pad_shifted = circshift(circshift(A_k_zero_pad, 7, 1), 14, 2);
    % A_k_zero_pad_shifted_cropped = A_k_zero_pad_shifted(1:L, 1:L, :);
    A( :, :, k, :) = A_k;
end

end