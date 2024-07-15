function pI_curr = calc_patch_likelihood_function(L, patches, volume_curr, sigma2, rotations)
    Mask = zeros(2*L, 2*L);
    Mask(1:L, 1:L) = 1;
    Mask_k = fft2(Mask);

    patches_norm_squared = (squeeze(pagenorm(permute(patches, [2, 3, 1])))).^2;
    Z_patches_k = fft2(permute(patches, [2, 3, 1]), 2*L, 2*L);

    Pfrots = vol_project(volume_curr, rotations);
    Z_Pfrots_k = fft2(Pfrots, 2*L, 2*L);
    Z_Pfrots2_k = fft2(Pfrots.^2, 2*L, 2*L);
    C = real(ifft2(conj(Z_Pfrots_k) .* permute(Z_patches_k, [1, 2, 4, 3])));
    patches_estimates_norm_squared = real(ifft2(Mask_k .* conj(Z_Pfrots2_k)));

    norm_squared_value = permute(patches_norm_squared, [2, 3, 4, 1]) - 2 * C + patches_estimates_norm_squared;
    norm_squared_value_normalized = norm_squared_value - min(norm_squared_value, [], [1, 2, 3]);
    
    pI_curr = exp(-norm_squared_value_normalized / (2 * sigma2));
    pI_curr = pI_curr ./ sum(pI_curr, [1, 2, 3]);
    pI_curr = permute(pI_curr, [3, 1, 2, 4]);
end
