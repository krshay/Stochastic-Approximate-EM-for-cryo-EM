function patches = micrograph2patches(I, L)
N = size(I, 1);
Nd = (N / L)^2;
patches_cell = mat2cell(I, L*ones(N / L, 1), L*ones(N / L, 1));
patches = zeros(Nd, L, L);
count = 0;
for i=1:sqrt(Nd)
    for j=1:sqrt(Nd)
        count = count + 1;
        patches( count, :, :) = patches_cell{i, j};
    end
end
end
