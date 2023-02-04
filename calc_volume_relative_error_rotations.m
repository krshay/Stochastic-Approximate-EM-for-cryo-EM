function [err, min_error_rotation, vol_est_rotated] ...
    = calc_volume_relative_error_rotations(vol_true, vol_est, K)

L = size(vol_true, 1);
rot_matrices = genRotationsGrid(K);
number_of_rotation_matrices = size(rot_matrices, 3);

vols_est_rotated = zeros(L, L, L, number_of_rotation_matrices);
for i=1:number_of_rotation_matrices
    vols_est_rotated( :, :, :, i) = fastrotate3d(vol_est, ...
        rot_matrices( :, :, i));
end
relative_errors = squeeze(sum((vol_true - vols_est_rotated).^2, [1, 2, 3])) ...
    / sum(vol_true.^2, 'all');
[err, min_rot_index] = min(relative_errors);
min_error_rotation = rot_matrices( :, :, min_rot_index);
vol_est_rotated = vols_est_rotated( :, :, :, min_rot_index);
end