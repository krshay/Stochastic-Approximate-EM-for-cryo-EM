function recovered_hat_v = expand_vol_psilms(x_lms, N, jball, Psilms, L)

if mod(L,2) ~= 0
    recovered_hat_v = zeros(2*N+1, 2*N+1, 2*N+1);
else
    recovered_hat_v = zeros(2*N, 2*N, 2*N);
end

for ii=1:length(Psilms)
    recovered_hat_v(jball) = recovered_hat_v(jball) + Psilms{ii} * x_lms{ii}(:);
end
end
