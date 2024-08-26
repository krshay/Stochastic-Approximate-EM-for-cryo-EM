close all;
clear variables;

warning('off');
addpath(genpath('../../ASPIRE'))
addpath('../../easyspin-5.2.33/easyspin')
addpath('../../SphericalHarmonics')


rng(1);


% Length of each dimension of the volume.
LL = 240;

% Length of each dimension of the downsampled volume.
L = 7;

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

% vol = load("shepp_logan_11.mat");
% vol = vol.shepp_logan_11;


%% Initializations


%% ell_max = 6
% Volume expansion in 3-D Fourier-Bessel basis
ell_max = 4;
r_cut = 1 / 2;
rad_size = floor(L / 2);
% Generate 3-D Fourier-Bessel basis
[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
    r_cut, ell_max, L);

vol_true_downsampled = cryo_downsample(vol, L);
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

upsilon0 = 0.55; % 0.4971;
rho_init = zeros(2*L, 2*L);
for l=1:NUM_SHIFTS
    if (shifts(l, 1) == L) && (shifts(l, 2) == L)
        rho_init(shifts(l, 1) + 1, shifts(l, 2) + 1) = upsilon0;
    else
        rho_init(shifts(l, 1) + 1, shifts(l, 2) + 1) = (1 - upsilon0) / (NUM_SHIFTS - 1);
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


STOCHASTIC_FACTOR = 1;

rotations = genRotationsGrid(20);


%% Calculate projection matrices
tic
PR_matrix_all_rotations_all_shifts = calc_projection_matrices(L, psi_lsNn, lms_list, n_list, rotations, shifts);
toc
PR_matrix_all_rotations_all_shifts_transpose = pagetranspose(PR_matrix_all_rotations_all_shifts);

tic
ATA = pagemtimes(PR_matrix_all_rotations_all_shifts_transpose, PR_matrix_all_rotations_all_shifts);
toc



vol_true_downsampled = cryo_downsample(vol, L);
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

%% Micrograph generation
% rotations = genRotationsGrid(15);
% average_norm_squared = calc_average_norm_squared_projection(vol, rotations);

SNR = 1;

N = 500;
N = round(N / L) * L;
W = 2 * L - 1;
gamma = 0.4;
I_clean = generate_clean_micrograph_2d_rots_new(vol_downsampled_trunc, ...
    W, N, round(gamma * (N ^ 2 / L ^ 2)), 1);
sigma2 = 0.0001;%0.001;
sigma = sqrt(sigma2);
I = I_clean + normrnd(0, sigma, size(I_clean));

patches = micrograph2patches(I, L);
patches_clean = micrograph2patches(I_clean, L);
patches = patches(1:250, :, :);
Nd = size(patches, 1);
patches_reshaped = permute(reshape(patches, Nd, L^2), [2, 3, 4, 5, 1]);

tic
ATI = pagemtimes(PR_matrix_all_rotations_all_shifts_transpose, ...
    patches_reshaped);
toc

'Finished micrographs generation.'


x_init_vec = zeros(size(lms_list, 1), 1);
for i=1:length(lms_list)
    ell = lms_list(i, 1);
    m = lms_list(i, 2);
    s = lms_list(i, 3);
    x_init_vec(i) = x_init{ell + 1}(s, m + ell + 1);
end

NUM_ITERS = 100;

%% Approximate expectation-maximization using least-squares solution
x_est = EM(patches, x_init, rho_init, ...
    L, rotations, sigma2, jball, Psilms, ell_max, s_lens, shifts, ...
    NUM_ITERS, ATA, ATI);
