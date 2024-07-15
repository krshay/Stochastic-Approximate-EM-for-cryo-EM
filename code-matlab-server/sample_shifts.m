function sampled_shifts = sample_shifts(shifts, L, sampling_factor)
shifts(shifts( :, 1) == L, :) = [];
shifts(shifts( :, 2) == L, :) = [];
shifts_taken = randperm(size(shifts, 1), sampling_factor); 
sampled_shifts = shifts(shifts_taken, :);

shifts_1d = 0:2:2*L-1;
shifts = cartprod(shifts_1d, shifts_1d);
shifts(shifts( :, 1) == L, :) = [];
shifts(shifts( :, 2) == L, :) = [];
sampled_shifts = shifts;
sampled_shifts(length(sampled_shifts) + 1, :) = [L, L];
end