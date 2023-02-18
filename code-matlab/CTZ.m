function CTZim = CTZ(im, l, L)
Zim = zeros(2*L, 2*L, size(im, 3), size(im, 4));
Zim(1:L, 1:L, :, :) = im;
TZim = circshift(circshift(Zim, l(1), 1), l(2), 2);
CTZim = TZim(1:L, 1:L, :, :);
end