function S = pI_l_rot_x(patches, l, rot, volume)
L = size(volume, 1);
Pfrot = vol_project(volume, rot);

CTZPfrot = CTZ(Pfrot, l, L);
% tmp = permute(CTZPfrot, [3, 1, 2]);
S = squeeze(sum((patches - permute(CTZPfrot, [3, 1, 2])) .^ 2, [2, 3]));
end