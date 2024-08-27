function Psilms = generate_psilms(jball, jjorigin, Y, dems_Y, dels_Y, ...
    r_cut, N, r_select_ratio, ell, ell_max, L)

fname='../data-mat/Besselj_L200_S500.mat';
load(fname, 'B');
 
% Choose ell, s, and roots by Nyquist criterion.
B = B(B(:, 3)< 2 * pi * N * r_cut * r_select_ratio & B(:,1) <= ell_max, :);
ell_list = B(:, 1);
s_list = B(:, 2);
roots_ls = B(:, 3);
 
if (mod(L,2) ~= 0)
    [x, y, z] = meshgrid(-N:N, -N:N, -N:N);
else
    [x, y, z] = meshgrid(-N:N-1, -N:N-1, -N:N-1);
end
 
% Convert to r, theta, phi.
r = sqrt( x.^2 + y.^2 + z.^2 ) / N / 2;
siz_lms = length(find(ell_list == ell))*(2*ell+1);
siz_grid = numel(jball);
Psilms = zeros(siz_grid, siz_lms);

% Spherical bessel function generation.
jl_func = @(l, x) besselj(l+1/2, x).*sqrt(pi./(2*x)); % x > 0

% j_{l, s} computation.
ind_ls = find(ell_list == ell);
siz_ls = numel(ind_ls);
j_ls = zeros(siz_grid, siz_ls); 
for ss=1:siz_ls
    root_ls = roots_ls(ind_ls(ss)); % s-th root of j_l(r)
    normalcont_ls = sqrt(2) / abs(r_cut^(3/2) * jl_func(ell+1, root_ls));
    j_ls(:,ss) = jl_func(ell, root_ls * (r(jball) / r_cut)) * normalcont_ls;
    % be careful with the origin point
    if ell == 0
        j_ls(jjorigin, ss) = normalcont_ls;
    else
        j_ls(jjorigin, ss) = 0;
    end
end
 
% Y_{l, m} computation.
Y_lm = Y(dels_Y == ell, :).';

% Fourier-Bessel basis functions computation.
indcol = 0;
for ii=1:(2 * ell + 1)
    Psilms(:, indcol+1:indcol+siz_ls) = j_ls .* Y_lm(:, ii);
    indcol = indcol + siz_ls;
end
end
