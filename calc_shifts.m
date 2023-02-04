function shifts = calc_shifts(L)
shifts_1d = 0:2*L-1;
shifts = cartprod(shifts_1d, shifts_1d);

end