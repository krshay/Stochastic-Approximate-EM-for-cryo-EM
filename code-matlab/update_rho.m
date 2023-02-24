function rho_updated = update_rho(pl_rot_curr, Nd)
% This function updates rho according to Eq. (37).
%
% Inputs:
%
% pl_rot_curr: a 4D array that represents the orobability function
% Nd: the number of patches in the current iteration
%
% Outputs:
%
% rho_updated: the updated rho

rho_updated = squeeze(sum(pl_rot_curr, [1, 4])) / Nd;
end