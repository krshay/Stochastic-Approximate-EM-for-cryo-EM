function A = shifted_proj_matrix(proj_matrix_all_rotations, L, shifts)

A = zeros(L, L, size(proj_matrix_all_rotations, 3), ...
    size(proj_matrix_all_rotations, 4), size(shifts, 1));
A_zero_pad = zeros(2 * L, 2 * L, size(proj_matrix_all_rotations, 3), ...
    size(proj_matrix_all_rotations, 4));

A_zero_pad(1:L, 1:L, :, :) = proj_matrix_all_rotations;

for s=1:size(shifts, 1)
    A_zero_pad_shifted = circshift(circshift(A_zero_pad, ...
        shifts(s, 1), 1), shifts(s, 2), 2);
    A_zero_pad_shifted_cropped = A_zero_pad_shifted(1:L, 1:L, :, :);
    A( :, :, :, :, s) = A_zero_pad_shifted_cropped;
end

end