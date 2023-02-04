function betas = sph_Bessel_to_2D_PSWF_factors(ell_max, n_list, s_max, L)
% Compute factors proportional to the integral of 2D PSWFs and spherical
% Bessel functions, used to evaluate the bispectrum.
% 
% Inputs: 
%   * maxL: cutoff for spherical harmonics expansion
%   * n_list: list of number of radial frequencies for each angular
%   frequency of 2D PSWFs
%   * maxS: maximum radial frequency for spherical Bessel functions
%   * L: length of volume or projection
% 
% Outputs:
%   * betas: a cell array indexed (order of spherical harmonics) x (angular
%   frequency for 2D PSWFs) containing the factors.
% 
% Eitan Levin, June 2018

c = pi*L;
[k, w] = lgwt(200*n_list(1), 0, 1);
w = w.*k;

Y_l = YN2YL(getSH(ell_max, [0, pi/2], 'complex'));
j_l = generate_spherical_bessel_basis(ell_max, s_max*ones(ell_max+1,1), ...
    1/2, k/2);

betas = cell(ell_max+1,1);
for ell = 0:ell_max
    betas{ell+1} = cell(ell+1, 1);
    for N=0:ell % only compute for positive N
        [R_Nn, alpha_Nn_2D] = PSWF_radial_2D(abs(N), n_list(abs(N)+1) ...
            - 1, c, k); % generate 2D radial prolates
        R_Nn = bsxfun(@times, R_Nn, 4./alpha_Nn_2D(:).^2.');
        j_l_curr = sqrt(2*pi)*Y_l{ell+1}(ell+1+N)*j_l{ell+1};
        
        betas{ell+1}{N+1} = j_l_curr.'*diag(w)*R_Nn; % s x n
    end
end
end
