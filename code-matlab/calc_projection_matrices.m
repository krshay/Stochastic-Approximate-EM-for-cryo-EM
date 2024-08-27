function PR_matrix_all_rotations_all_shifts = calc_projection_matrices(L, psi_lsNn, lms_list, n_list, rotations, shifts)

NUM_LMS = size(lms_list, 1);
NUM_ROTATIONS = size(rotations, 3);
NUM_SHIFTS = size(shifts, 1);

PR_matrix_all_rotations = proj_matrix(L, psi_lsNn, lms_list, n_list, rotations);

PR_matrix_all_rotations_all_shifts = shifted_proj_matrix(PR_matrix_all_rotations, L, shifts);

PR_matrix_all_rotations_all_shifts = reshape(PR_matrix_all_rotations_all_shifts, L^2, NUM_LMS, NUM_ROTATIONS, NUM_SHIFTS);

end