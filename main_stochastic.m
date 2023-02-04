close all;
clear variables;

warning('off');
addpath(genpath('/data/shaykreymer/PhD/ASPIRE'))
addpath('/data/shaykreymer/PhD/easyspin-5.2.33/easyspin')
addpath('/data/shaykreymer/PhD/SphericalHarmonics')

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
gamma = 0.4;
N = 143 * LL / L;
N = N - mod(N, LL);

Nd = (N / LL)^2;

W = 2 * LL - 1;

rotations = genRotationsGrid(15);
average_norm_squared = calc_average_norm_squared_projection(vol, rotations);

SNR = 1;

I_clean = generate_clean_micrograph_2d_rots_new(vol, ...
    W, N, round(gamma * (N ^ 2 / LL ^ 2)), 1);
sigma2 =  average_norm_squared / (LL^2 * SNR);
sigma = sqrt(sigma2);
I = I_clean + normrnd(0, sigma, size(I_clean));
N_downsampled = round(N * L / LL);
I_downsampled = cryo_downsample(I, N_downsampled);
N_downsampled_correct_factorial = N_downsampled - mod(N_downsampled, L);
I_downsampled = I_downsampled(1:N_downsampled_correct_factorial, 1:N_downsampled_correct_factorial);
patches = micrograph2patches(I_downsampled, L);
% patches_clean = micrograph2patches(I_clean, L);
Nd = size(patches, 1);

clear I;
clear I_clean;

'Finished micrograph generation.'

%% Initializations


%% ell_max = 4
% Volume expansion in 3-D Fourier-Bessel basis
ell_max = 10;
r_cut = 1 / 2;
rad_size = floor(L / 2);
% Generate 3-D Fourier-Bessel basis
[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
    r_cut, ell_max, L);

vol_true_trunc = cryo_downsample(vol, 11);

vol_init = normrnd(0, 1, [L, L, L]); vol_init = 50 * vol_init / sqrt(sum(vol_init.^2, 'all'));
vol_init = load("./L8_24122022/volume_curr_16.mat");
vol_init = vol_init.volume_curr;
[~, vol_init_trunc] = expand_vol_spherical_basis(vol_init, rad_size, ...
    ell_max, L, Psilms, jball);
vol_init_trunc = real(icfftn(vol_init_trunc));

err = calc_volume_relative_error_rotations(vol_true_trunc, vol_init_trunc, 20);
display(["Current estimation_error = ", err])

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
        rho_init(Ls(l, 1) + 1, Ls(l, 2) + 1) = (1 - beta0) / (NUMBER_OF_SHIFTS - 1);
    end
end

rho_init = load("./L8_24122022/rho_curr_16.mat");
rho_init = rho_init.rho_curr;

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


STOCHASTIC_FACTOR = 0.05;

rotations = genRotationsGrid(25);
% tic;
% tmp = load("gs_K552_ellmax6_13122022.mat");
% 
% gs = tmp.gs;
% toc;
tic;
gs = calc_g(Ls, L, rotations, psi_lsNn, lms_list, n_list);
toc;


[x_est, rho_est] = EM_3d_parallel_stochastic_noprecomputation_Ipsis(patches, x_init, ...
    rho_init, L, rotations, sigma2, jball, Psilms, ell_max, ...
    s_lens, n_list, lms_list, vol_true_trunc, Ls, ...
    Ls, psi_lsNn, STOCHASTIC_FACTOR, gs);












































% close all;
% clear variables;
% 
% warning('off');
% addpath(genpath('/data/shaykreymer/PhD/ASPIRE'))
% addpath('/data/shaykreymer/PhD/easyspin-5.2.33/easyspin')
% addpath('/data/shaykreymer/PhD/SphericalHarmonics')
% 
% rng(10);
% 
% % Length of each dimension of the volume.
% L = 11;
% 
% LL = 192;
% % Truncation parameter
% ell_max = 6;
% 
% %% Volume generation
% vol = load("TRPV1.mat");
% vol = vol.TRPV1_vol;
% 
% vol = LL * vol / sqrt(sum(vol.^2, 'all'));
% 
% [X, Y, Z] = meshgrid(-LL/2:LL/2-1, -LL/2:LL/2-1, -LL/2:LL/2-1);
% tmp = sqrt((X.^2 + Y.^2 + Z.^2)) <= LL/2;
% vol_cropped = zeros(size(vol));
% for i=1:LL
%     for j=1:LL
%         for k=1:LL
%             if tmp(i, j, k)
%                 vol_cropped(i, j, k) = vol(i, j, k);
%             end
%         end
%     end
% end
% vol = vol_cropped;
% 
% 
% %% Volume expansion in 3-D Fourier-Bessel basis
% r_cut = 1 / 2;
% rad_size = floor(L / 2);
% % Generate 3-D Fourier-Bessel basis
% [Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
%     r_cut, ell_max, L);
% % Expand volume in 3-D Fourier-Bessel basis
% % [~, vol_true_trunc] = expand_vol_spherical_basis(vol, rad_size, ell_max, ...
% %     L, Psilms, jball);
% % vol_true_trunc = vol;real(icfftn(vol_true_trunc));
% % WriteMRC(vol_true_trunc, 1, 'volume.mrc');
% % x_lms = expand_vol_spherical_basis(vol_true_trunc, rad_size, ell_max, L, ...
% %     Psilms, jball);
% 
% 
% %% Micrograph generation
% gamma = 0.4;
% N = 1000 * LL / L;
% N = N - mod(N, LL);
% 
% SNR = 1;
% 
% Nd = (N / LL)^2;
% 
% rotations = genRotationsGrid(15);
% average_norm_squared = calc_average_norm_squared_projection(vol, rotations);
% 
% W = 2 * LL - 1;
% 
% sigma2 = average_norm_squared / (SNR * LL^2);
% sigma = sqrt(sigma2);
% 
% number_of_micrographs = 1;
% 
% patches = cell(number_of_micrographs, 1);
% Is = cell(number_of_micrographs, 1);
% 
% for mic=1:number_of_micrographs
%     [I_clean_mic, ~, proj1] = generate_clean_micrograph_2d_rots_new(vol, ...
%         W, N, round(gamma * (N ^ 2 / LL ^ 2)), mic);
%     I_mic_original = I_clean_mic + normrnd(0, sigma, size(I_clean_mic));
%     I_mic = cryo_downsample(I_mic_original, round(N * L / LL));
%     patches_mic = micrograph2patches(I_mic, L);
%     clear I_mic_original;
%     clear I_clean_mic;
%     patches{mic} = patches_mic;
%     Is{mic} = I_mic;
%     Nd_eachpatch = size(patches_mic, 1);
% end
% 
% 'Finished micrograph generation.'
% 
% %% Initializations
% 
% ell_max = 6;
% 
% patches_all = zeros(Nd_eachpatch*number_of_micrographs, L, L);
% for nn=1:number_of_micrographs
%     patches_all((nn - 1) * Nd_eachpatch + 1: nn*Nd_eachpatch, :, :) = patches{nn};
% end
% 
% %% Volume expansion in 3-D Fourier-Bessel basis
% r_cut = 1 / 2;
% rad_size = floor(L / 2);
% % Generate 3-D Fourier-Bessel basis
% [Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
%     r_cut, ell_max, L);
% 
% vol_init = normrnd(0, 1, size(cryo_downsample(vol, L))); vol_init = vol_init / sqrt(sum(vol_init.^2, 'all'));
% 
% 
% % vol_init = load('/data/shaykreymer/wigner/MTD-3D-EM_L11_BPTI_v21112022_no_expansion_downsampled/Results_L8/volume_curr_7.mat');
% % vol_init = vol_init.volume_curr;
% [~, vol_init_trunc] = expand_vol_spherical_basis(vol_init, rad_size, ...
%     ell_max, L, Psilms, jball);
% vol_init_trunc = real(icfftn(vol_init_trunc)); 
% vol_init_trunc = 50 * vol_init_trunc / norm(vol_init_trunc, "fro");
% % err    = calc_volume_relative_error_rotations(vol_true_trunc, vol_init_trunc, 20)
% 
% % vol_init_trunc = load("/data/shaykreymer/PhD/MTD-3D-EM_V20072022_L11_TRPV1/L4_results_26072022/volume_curr_12.mat");
% % vol_init_trunc = vol_init_trunc.volume_curr;
% 
% % clos e
% 
% x_init = expand_vol_spherical_basis(vol_init_trunc, rad_size, ell_max, ...
%     L, Psilms, jball);
% 
% 
% SAMPLING_FACTOR = 150;
% Ls = calc_shifts(L);
% Ls_sampled = Ls; % sample_shifts(Ls, L, SAMPLING_FACTOR);
% SAMPLING_FACTOR = size(Ls_sampled, 1);
% 
% beta0 = 0.55;
% rho_init = zeros(2*L, 2*L);
% for l=1:SAMPLING_FACTOR
%     if (Ls_sampled(l, 1) == L) && (Ls_sampled(l, 2) == L)
%         rho_init(Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1) = beta0;
%     else
%         rho_init(Ls_sampled(l, 1) + 1, Ls_sampled(l, 2) + 1) = (1 - beta0) / (SAMPLING_FACTOR - 1);
%     end
% end
% 
% % rho_init = load('/data/shaykreymer/wigner/MTD-3D-EM_L11_BPTI_v21112022_no_expansion_downsampled/Results_L8/rho_curr_7.mat');
% % rho_init = rho_init.rho_curr;
% 
% % rho_init = load("/data/shaykreymer/PhD/MTD-3D-EM_V20072022_L11_TRPV1/L4_results_26072022/rho_curr_12.mat");
% % rho_init = rho_init.rho_curr;
% 
% % Projection expansion.
% 
% s_lens = gen_s_list(ell_max, r_cut, rad_size);
% 
% % PSWF parameters
% beta = 1;       % Bandlimit ratio (between 0 and 1) -
% % smaller values stand for greater oversampling
% T = 1e-1;       % Truncation parameter
% realFlag = 1;   % Flag indicating whether the data is assumed to
% % be real-valued
% 
% % Obtain coefficients while evaluating PSWFs.
% projection = vol_project(cryo_downsample(vol, L), eye(3));
% [PSWF_coeffs, PSWF_Nn_p] = pswf_t_f(projection, rad_size, beta, T, ...
%     realFlag, []);
% 
% num_coeffs = find(PSWF_Nn_p.ang_freq <= ell_max, 1, 'last');
% PSWF_Nn_p.n_list = PSWF_Nn_p.n_list(1:ell_max+1);
% PSWF_Nn_p.alpha_Nn = PSWF_Nn_p.alpha_Nn(1:num_coeffs);
% PSWF_Nn_p.ang_freq = PSWF_Nn_p.ang_freq(1:num_coeffs);
% PSWF_Nn_p.rad_freq = PSWF_Nn_p.rad_freq(1:num_coeffs);
% PSWF_Nn_p.samples = PSWF_Nn_p.samples( :, 1:num_coeffs);
% n_list = PSWF_Nn_p.n_list;
% betas = sph_Bessel_to_2D_PSWF_factors(ell_max, n_list, ...
%     s_lens(1), rad_size);
% 
% psi_Nn = PSWF_2D_full_cart(PSWF_Nn_p.n_list, L, beta);
% 
% psi_lsNn = calc_psi_lNs(ell_max, psi_Nn, betas, L, s_lens, n_list);
% 
% lms_list = calc_lms_list(ell_max, s_lens);
% 
% % Precomputations.
% % tic
% % Ipsis_cell_all = preprocess_Ipsi(L, Ls_sampled, psi_lsNn, ell_max, s_lens, n_list, patches_all);
% % toc
% 
% 
% tic
% psipsi = calc_psi_multiplications_parallel(psi_lsNn, L, ell_max , n_list, s_lens, Ls_sampled);
% toc
% psipsi_cell = cell(SAMPLING_FACTOR, 1);
% for l=1:SAMPLING_FACTOR
%     psipsi_cell{l} = squeeze(psipsi(l, :, :, :, :, :, :, :, :));
% end
% whos psipsi
% clear psipsi
% 
% rotations = genRotationsGrid(25);
% 
% STOCHASTIC_FACTOR = 0.05;
% 
% % [x_est, rho_est] = EM_3d_parallel_stochastic(patches_all, x_init, ...
% %     rho_init, L, rotations, sigma2, jball, Psilms, ell_max, ...
% %     s_lens, n_list, lms_list, vol_true_trunc, Ls, ...
% %     Ls_sampled, psipsi_cell, Ipsis_cell_all, STOCHASTIC_FACTOR);
% 
% [x_est, rho_est] = EM_3d_parallel_stochastic_noprecomputation_Ipsis(patches_all, x_init, ...
%     rho_init, L, rotations, sigma2, jball, Psilms, ell_max, ...
%     s_lens, n_list, lms_list, cryo_downsample(vol, L), Ls, ...
%     Ls_sampled, psipsi_cell, psi_lsNn, STOCHASTIC_FACTOR);
% 
% %% Analysis
% vol_est = real(icfftn(expand_vol_psilms(x_est, floor(L / 2), jb1all, ...
%     Psilms, L)));
% [err, min_error_rotation, vol_est_rotated] ...
%     = calc_volume_relative_error_rotations(vol_true_trunc, vol_est, 40);
% display(['Estimation error = ', num2str(err)])
% [resA, fighandle] = plotFSC(vol, vol_est_rotated, 0.5, 40/13);
