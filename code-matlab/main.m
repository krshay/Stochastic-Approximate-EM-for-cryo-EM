close all;
clear variables;

warning('off');
addpath(genpath('./ASPIRE'))
addpath('./easyspin-5.2.33/easyspin')
addpath('./SphericalHarmonics')

rng(10);

% Length of each dimension of the volume.
L = 11;

%% Volume generation
vol = load("emd_2660_cropped.mat");
vol = vol.emd_2660_cropped;

LL = size(vol, 1);

vol = LL * vol / sqrt(sum(vol.^2, 'all'));

[X, Y, Z] = meshgrid(-LL/2:LL/2-1, -LL/2:LL/2-1, -LL/2:LL/2-1);
tmp = sqrt((X.^2 + Y.^2 + Z.^2)) <= LL/2;
vol_cropped = zeros(size(vol));
for i=1:LL
    for j=1:LL
        for k=1:LL
            if tmp(i, j, k)
                vol_cropped(i, j, k) = vol(i, j, k);
            end
        end
    end
end
vol = vol_cropped;


% %% Volume expansion in 3-D Fourier-Bessel basis
% r_cut = 1 / 2;
% rad_size = floor(L / 2);
% % Generate 3-D Fourier-Bessel basis
% [Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
%     r_cut, ell_max, L);
% % Expand volume in 3-D Fourier-Bessel basis
% [~, vol_true_trunc] = expand_vol_spherical_basis(vol, rad_size, ell_max, ...
%     L, Psilms, jball);
% vol_true_trunc = real(icfftn(vol_true_trunc));
% % WriteMRC(vol_true_trunc, 1, 'volume.mrc');
% x_lms = expand_vol_spherical_basis(vol_true_trunc, rad_size, ell_max, L, ...
%     Psilms, jball);


%% Micrograph generation
% Density of projection images in the measurement
gamma = 0.4;
% Size of the micrograph.
N = 6000 * LL / L;
N = N - mod(N, LL);

W = 2 * LL - 1;

rotations = genRotationsGrid(15);
average_norm_squared = calc_average_norm_squared_projection(vol, ...
    rotations);

SNR = 1;

I_clean = generate_clean_micrograph_2d_rots_new(vol, ...
    W, N, round(gamma * (N ^ 2 / LL ^ 2)), 1);
sigma2 = average_norm_squared / (LL^2 * SNR);
sigma = sqrt(sigma2);
I = I_clean + normrnd(0, sigma, size(I_clean));
N_downsampled = round(N * L / LL);
I_downsampled = cryo_downsample(I, N_downsampled);
N_downsampled_correct_factorial = N_downsampled - mod(N_downsampled, L);
I_downsampled = I_downsampled(1:N_downsampled_correct_factorial, ...
    1:N_downsampled_correct_factorial);
patches = micrograph2patches(I_downsampled, L);
Nd = size(patches, 1);

%% Initializations
% Volume expansion in 3-D Fourier-Bessel basis
ell_max = 10;
r_cut = 1 / 2;
rad_size = floor(L / 2);
% Generate 3-D Fourier-Bessel basis
[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis( ...
    rad_size, r_cut, ell_max, L);

vol_init = normrnd(0, 1, [L, L, L]);
vol_init = 50 * vol_init / sqrt(sum(vol_init.^2, 'all'));
[~, vol_init_trunc] = expand_vol_spherical_basis(vol_init, rad_size, ...
    ell_max, L, Psilms, jball);
vol_init_trunc = real(icfftn(vol_init_trunc));

x_init = expand_vol_spherical_basis(vol_init_trunc, rad_size, ell_max, ...
    L, Psilms, jball);

Ls = calc_shifts(L);
NUMBER_OF_SHIFTS = size(Ls, 1);

beta0 = 0.55;
rho_init = zeros(2*L, 2*L);
for l=1:NUMBER_OF_SHIFTS
    if (Ls(l, 1) == L) && (Ls(l, 2) == L)
        rho_init(Ls(l, 1) + 1, Ls(l, 2) + 1) = beta0;
    else
        rho_init(Ls(l, 1) + 1, ...
            Ls(l, 2) + 1) = (1 - beta0) / (NUMBER_OF_SHIFTS - 1);
    end
end

% Projection expansion.
s_lens = gen_s_list(ell_max, r_cut, rad_size);

% PSWF parameters
beta = 1;       % Bandlimit ratio (between 0 and 1) -
% smaller values stand for greater oversampling
T = 1e-1;       % Truncation parameter
realFlag = 1;   % Flag indicating whether the data is assumed to
% be real-valued

% Obtain coefficients while evaluating PSWFs.
projection = vol_project(vol_true_trunc, eye(3));
[PSWF_coeffs, PSWF_Nn_p] = pswf_t_f(projection, rad_size, beta, T, ...
    realFlag, []);

num_coeffs = find(PSWF_Nn_p.ang_freq <= ell_max, 1, 'last');
PSWF_Nn_p.n_list = PSWF_Nn_p.n_list(1:ell_max+1);
PSWF_Nn_p.alpha_Nn = PSWF_Nn_p.alpha_Nn(1:num_coeffs);
PSWF_Nn_p.ang_freq = PSWF_Nn_p.ang_freq(1:num_coeffs);
PSWF_Nn_p.rad_freq = PSWF_Nn_p.rad_freq(1:num_coeffs);
PSWF_Nn_p.samples = PSWF_Nn_p.samples( :, 1:num_coeffs);
n_list = PSWF_Nn_p.n_list;
betas = sph_Bessel_to_2D_PSWF_factors(ell_max, n_list, ...
    s_lens(1), rad_size);

psi_Nn = PSWF_2D_full_cart(PSWF_Nn_p.n_list, L, beta);

psi_lsNn = calc_psi_lNs(ell_max, psi_Nn, betas, L, s_lens, n_list);

lms_list = calc_lms_list(ell_max, s_lens);

gs = calc_g(Ls, L, rotations, psi_lsNn, lms_list, n_list);

STOCHASTIC_FACTOR = 0.05;
rotations = genRotationsGrid(25);
NUM_ITERS = 10;

%% Stochastic approximate expectation-maximization
[x_est, rho_est] = EM_stochastic_parallel(patches, x_init, rho_init, ...
    L, rotations, sigma2, jball, Psilms, ell_max, s_lens, n_list, ...
    lms_list, Ls, psi_lsNn, STOCHASTIC_FACTOR, gs, NUM_ITERS);
