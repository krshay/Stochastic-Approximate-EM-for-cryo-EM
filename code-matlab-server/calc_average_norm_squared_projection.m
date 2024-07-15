function average_norm_squared = calc_average_norm_squared_projection(volume, rotations)

sum_norm = 0;
for i=1:size(rotations, 3)
    projection = vol_project(volume, rotations( :, :, i));
    sum_norm = sum_norm + norm(projection)^2;
end
average_norm_squared = sum_norm / size(rotations, 3);
end