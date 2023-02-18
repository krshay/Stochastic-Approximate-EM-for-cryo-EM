function rho_updated = rho_step(pl_rot_curr, Nd)

rho_updated = squeeze(sum(pl_rot_curr, [1, 4])) / Nd;
end