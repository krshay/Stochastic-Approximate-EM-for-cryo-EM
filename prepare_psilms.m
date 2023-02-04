function [jball, jjorigin, Y, dems_Y, dels_Y, xyplane, jball_2D] = ...
    prepare_psilms(N, r_cut, ell_max, L)

if (mod(L,2) ~= 0 )
    [x, y, z] = meshgrid(-N:N, -N:N, -N:N);
else
    [x, y, z] = meshgrid(-N:N-1, -N:N-1, -N:N-1);
end

x = x / N / 2;
y = y / N / 2;
z = z / N / 2;
% Convert to r, theta, phi
r = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z ./ r);
phi = acos(x ./ (r .* sin(theta)));
j1 = find(y < 0);
phi(j1) = 2*pi - phi(j1) ;
MACHINEEPS = 1e-15;
phi(theta < MACHINEEPS | theta > pi - MACHINEEPS) = 0; % z = 1 or -1, x=y=0
phi = real( phi); %enforce real numbers
jorigin = find( r < MACHINEEPS ); % be careful at origin point
theta(jorigin) = 0;
phi(jorigin) = 0;
jball = find(r < r_cut); % compact support on r < r_cut
jjorigin = find(r(jball) < MACHINEEPS );
xyplane = find(z(jball) == 0);
 
Y = getSH(ell_max, [phi(jball), theta(jball)], 'complex')';
dems_Y = zeros((ell_max+1)^2, 1); 
dels_Y = dems_Y;
ind_Y = 1;
for ell=0:ell_max
    for m=-ell:ell
        dems_Y(ind_Y) = m;
        dels_Y(ind_Y) = ell;
        ind_Y = ind_Y + 1;
    end
end
 
% Generate 2D jball:
if (mod(L,2) ~= 0)
    [x, y] = meshgrid(-N:N, -N:N);
else
    [x, y] = meshgrid(-N:N-1, -N:N-1);
end
jball_2D = find(sqrt(x.^2 + y.^2)/N/2 < r_cut);
end
