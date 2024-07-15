function image = projection_pswf(x_lms, betas, omega, PSWF_Nn_p, ...
    ell_max, s_lens)
img = 0;
for ell=0:ell_max
    D = wignerd(ell, omega);
    for N=-ell:ell
        for n=1:PSWF_Nn_p.n_list(abs(N)+1)
%             addition = sum(x_lms{ell+1}( :, :) .* betas{ell+1} ...
%                 {abs(N) + 1}(1:s_lens(ell+1), n), 'all') ...
%                 * PSWF_Nn_p.samples( :, logical((PSWF_Nn_p.ang_freq == abs(N)) ...
%                 .* (PSWF_Nn_p.rad_freq == n)));
%             if N == 0
%                 img = img + addition;
%             elseif mod(N, 2) == 0
%                 img = img + 2*real(addition);
%             elseif mod(N, 2) == 1
%                 img = img + 2*imag(addition);
%             end
    Nn_idx = find(logical((PSWF_Nn_p.ang_freq == abs(N)) ...
                .* (PSWF_Nn_p.rad_freq == n)) == 1);
            if N >= 0
                img = img + (x_lms{ell+1}( :, :) * D(N + ell + 1, :).').'...
                    * betas{ell+1}{abs(N) + 1}(1:s_lens(ell+1), n) ...
                * PSWF_Nn_p.alpha_Nn(Nn_idx) * PSWF_Nn_p.samples( :, Nn_idx);
            elseif N < 0
                img = img + (x_lms{ell+1}( :, :) * D(N + ell + 1, :).').'...
                * (-1)^N * conj(betas{ell+1}{abs(N) + 1}(1:s_lens(ell+1), n)) ...
                * conj(PSWF_Nn_p.alpha_Nn(Nn_idx) * PSWF_Nn_p.samples( :, Nn_idx));
            end
        end
    end
end
N = PSWF_Nn_p.L;
x_1d_grid = -N:1:N;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= N);

image = zeros(2*N+1, 2*N+1);
image(points_inside_the_circle) = img;
end