close all;
clear variables;

warning('off');
addpath(genpath('../../ASPIRE'))
addpath('../../easyspin-5.2.33/easyspin')
addpath('../../SphericalHarmonics')

rng(10);

% Length of each dimension of the volume.
LL = 240;

% Length of each dimension of the downsampled volume.
L = 11;

% directory = '../EMPIAR_10028_data/';
% files = dir(fullfile(directory, '*.mrcs'));
% 
% patches = [];
% for i=[1, 7, 9] %, 11, 12, 14]%, 15, 16, 17, 18, 19, 20]
%     %% Micrographs load
%     micrograph = ReadMRC(fullfile(files(i).folder, files(i).name));
%     micrograph = micrograph - mean(micrograph, "all");
%     micrograph = micrograph / norm(micrograph, "fro");
% 
%     micrograph_downsampled = cryo_downsample(micrograph, ...
%         round(size(micrograph, 1) * L / LL));
% 
%     micrograph_downsampled = micrograph_downsampled / norm(micrograph_downsampled, "fro");
% 
%     N = size(micrograph_downsampled, 1);
%     N = N - mod(N, L);
%     micrograph_downsampled = micrograph_downsampled(1:N, 1:N);
%     if i==1
%         noise_stack = zeros(L, L, 30);
%         noise_stack( :, :, 1) = micrograph_downsampled(31:41, 1:11);
%         noise_stack( :, :, 2) = micrograph_downsampled(42:52, 1:11);
%         noise_stack( :, :, 3) = micrograph_downsampled(53:63, 1:11);
%         noise_stack( :, :, 4) = micrograph_downsampled(111:121, 1:11);
%         noise_stack( :, :, 5) = micrograph_downsampled(122:132, 1:11);
%         noise_stack( :, :, 6) = micrograph_downsampled(133:143, 1:11);
%         noise_stack( :, :, 7) = micrograph_downsampled(159:169, 5:15);
%         noise_stack( :, :, 8) = micrograph_downsampled(170:180, 5:15);
%         noise_stack( :, :, 9) = micrograph_downsampled(8:18, 65:75);
%         noise_stack( :, :, 10) = micrograph_downsampled(31:41, 70:80);
%         noise_stack( :, :, 11) = micrograph_downsampled(104:114, 62:72);
%         noise_stack( :, :, 12) = micrograph_downsampled(115:125, 62:72);
%         noise_stack( :, :, 13) = micrograph_downsampled(169:179, 43:53);
%         noise_stack( :, :, 14) = micrograph_downsampled(36:46, 115:125);
%         noise_stack( :, :, 15) = micrograph_downsampled(47:57, 115:125);
%         noise_stack( :, :, 16) = micrograph_downsampled(58:68, 115:125);
%         noise_stack( :, :, 17) = micrograph_downsampled(69:79, 115:125);
%         noise_stack( :, :, 18) = micrograph_downsampled(36:46, 126:136);
%         noise_stack( :, :, 19) = micrograph_downsampled(47:57, 126:136);
%         noise_stack( :, :, 20) = micrograph_downsampled(58:68, 126:136);
%         noise_stack( :, :, 21) = micrograph_downsampled(69:79, 126:136);
%         noise_stack( :, :, 22) = micrograph_downsampled(104:114, 97:107);
%         noise_stack( :, :, 23) = micrograph_downsampled(104:114, 108:118);
%         noise_stack( :, :, 24) = micrograph_downsampled(164:174, 113:123);
%         noise_stack( :, :, 25) = micrograph_downsampled(164:174, 124:134);
%         noise_stack( :, :, 26) = micrograph_downsampled(39:49, 153:163);
%         noise_stack( :, :, 27) = micrograph_downsampled(50:60, 153:163);
%         noise_stack( :, :, 28) = micrograph_downsampled(61:71, 153:163);
%         noise_stack( :, :, 29) = micrograph_downsampled(61:71, 175:185);
%         noise_stack( :, :, 30) = micrograph_downsampled(72:82, 175:185);
% 
% 
%         cov_matrix = calc_cov_matrix(noise_stack, L);
%     end
%     sigma2 = mean(noise_stack.^2, "all");
%     patches = cat(1, patches, micrograph2patches(micrograph_downsampled, L));
% end
% Nd = size(patches, 1);


%% Volume generation
vol = load("TRPV1.mat");
vol = vol.TRPV1_vol;

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

%% Micrograph generation
gamma = 0.4;
N = 70;
N = N - mod(N, L);
Nd = (N / LL)^2;

W = 2 * LL - 1;

W = 2 * L - 1;

% rotations = genRotationsGrid(15);
% average_norm_squared = calc_average_norm_squared_projection(vol, rotations);

SNR = 1;

N = 11;
%% Initializations


%% ell_max = 6
% Volume expansion in 3-D Fourier-Bessel basis
ell_max = 6;
r_cut = 1 / 2;
rad_size = floor(L / 2);
% Generate 3-D Fourier-Bessel basis
[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
    r_cut, ell_max, L);

vol_true_downsampled = cryo_downsample(vol, 11);
vol_true_downsampled = vol_true_downsampled / norm(vol_true_downsampled, "fro");

vol_init = normrnd(0, 1, [L, L, L]); vol_init = vol_init / norm(vol_init, "fro");
% vol_init = load("./L8_results_24122022/volume_curr_9.mat");
% vol_init = vol_init.volume_curr;
[~, vol_init_trunc] = expand_vol_spherical_basis(vol_init, rad_size, ...
    ell_max, L, Psilms, jball);
vol_init_trunc = real(icfftn(vol_init_trunc));
vol_init_trunc = vol_init_trunc / norm(vol_init_trunc, "fro");
err = calc_volume_relative_error_rotations(vol_true_downsampled, vol_init_trunc, 20);
display(["Current estimation_error = ", err])

x_init = expand_vol_spherical_basis(vol_init_trunc, rad_size, ell_max, ...
    L, Psilms, jball);

shifts = calc_shifts(L);
NUM_SHIFTS = size(shifts, 1);

beta0 = 0.55;
rho_init = zeros(2*L, 2*L);
for l=1:NUM_SHIFTS
    if (shifts(l, 1) == L) && (shifts(l, 2) == L)
        rho_init(shifts(l, 1) + 1, shifts(l, 2) + 1) = beta0;
    else
        rho_init(shifts(l, 1) + 1, shifts(l, 2) + 1) = (1 - beta0) / (NUM_SHIFTS - 1);
    end
end

% rho_init = load("./L8_results_24122022/rho_curr_9.mat");
% rho_init = rho_init.rho_curr;


% Projection expansion.

s_lens = gen_s_list(ell_max, r_cut, rad_size);

% PSWF parameters
beta = 1;       % Bandlimit ratio (between 0 and 1) -
% smaller values stand for greater oversampling
T = 1e-1;       % Truncation parameter
realFlag = 1;   % Flag indicating whether the data is assumed to
% be real-valued

% Obtain coefficients while evaluating PSWFs.
projection = vol_project(vol_init_trunc, eye(3));
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

NUM_LMS = size(lms_list, 1);

STOCHASTIC_FACTOR = 1;

rotations = genRotationsGrid(15);
rotations = eye(3);
NUM_ROTATIONS = size(rotations, 3);

%% Check projection matrix
tic
PR_matrix_all_rotations = proj_matrix(L, psi_lsNn, lms_list, n_list, rotations);
toc
tic
PR_matrix_all_rotations_all_shifts = shifted_proj_matrix(PR_matrix_all_rotations, L, shifts);
toc
PR_matrix_all_rotations_all_shifts = reshape(PR_matrix_all_rotations_all_shifts, L^2, NUM_LMS, NUM_ROTATIONS, NUM_SHIFTS);

vol_true_downsampled = cryo_downsample(vol, 11);

[~, vol_downsampled_trunc] = expand_vol_spherical_basis(vol_true_downsampled, rad_size, ...
    ell_max, L, Psilms, jball);
vol_downsampled_trunc = real(icfftn(vol_downsampled_trunc));
vol_downsampled_trunc = vol_downsampled_trunc / norm(vol_downsampled_trunc, "fro");

x_true = expand_vol_spherical_basis(vol_downsampled_trunc, rad_size, ell_max, ...
    L, Psilms, jball);
x_true_vec = zeros(size(lms_list, 1), 1);
for i=1:length(lms_list)
    ell = lms_list(i, 1);
    m = lms_list(i, 2);
    s = lms_list(i, 3);
    x_true_vec(i) = x_true{ell + 1}(s, m + ell + 1);
end

I_clean = generate_clean_micrograph_2d_discrete_rots_projection_matrix(x_true_vec, ...
    W, N, round(gamma * (N ^ 2 / L ^ 2)), rotations, PR_matrix_all_rotations_all_shifts, L, 1);
% sigma2 =  average_norm_squared / (L^2 * SNR);
I_clean = real(reshape(PR_matrix_all_rotations_all_shifts( :, :, 1, 1) * x_true_vec, L, L));
sigma2 = 0.001;
sigma = sqrt(sigma2);
I = I_clean + normrnd(0, sigma, size(I_clean));
% N_downsampled = round(N * L / LL);
% I_downsampled = cryo_downsample(I, N_downsampled);
% N_downsampled_correct_factorial = N_downsampled - mod(N_downsampled, L);
% I_downsampled = I_downsampled(1:N_downsampled_correct_factorial, 1:N_downsampled_correct_factorial);
patches = micrograph2patches(I, L);
% patches_clean = micrograph2patches(I_clean, L);
Nd = size(patches, 1);

% clear I;
% clear I_clean;

'Finished micrographs loading.'



x_init_vec = zeros(size(lms_list, 1), 1);
for i=1:length(lms_list)
    ell = lms_list(i, 1);
    m = lms_list(i, 2);
    s = lms_list(i, 3);
    x_init_vec(i) = x_init{ell + 1}(s, m + ell + 1);
end

patches_reshaped = permute(reshape(patches, Nd, L^2), [2, 3, 4, 5, 1]);
% ls_A_dagger_patch = pagemldivide(PR_matrix_all_rotations_all_shifts, patches_reshaped);
% tic;
% gs = calc_g(Ls, L, rotations, psi_lsNn, lms_list, n_list);
% toc;
gs = 0;

NUM_ITERS = 1000;

% %% Stochastic approximate expectation-maximization
% [x_est, rho_est] = EM_stochastic_parallel(patches, x_init, rho_init, ...
%     L, rotations, sigma2, jball, Psilms, ell_max, s_lens, n_list, ...
%     lms_list, shifts, psi_lsNn, STOCHASTIC_FACTOR, gs, NUM_ITERS);

%% Stochastic approximate expectation-maximization using least-squares solution
x_est = EM_stochastic_parallel_ls(patches, x_init, rho_init, ...
    L, rotations, sigma2, jball, Psilms, ell_max, s_lens, n_list, ...
    lms_list, shifts, psi_lsNn, STOCHASTIC_FACTOR, gs, NUM_ITERS, PR_matrix_all_rotations_all_shifts);
