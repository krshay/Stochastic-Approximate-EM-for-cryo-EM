function cov_matrix = calc_cov_matrix(noise_stack, L)
cov_matrix = zeros(L^2, L^2);
noise_variables = reshape(noise_stack, L^2, size(noise_stack, 3));
for i=1:L^2
    for j=1:L^2
        cov_matrix(i, j) = mean((noise_variables(i, :) - mean(noise_variables(i, :))) .* (noise_variables(j, :) - mean(noise_variables(j, :))));
    end
end
end
