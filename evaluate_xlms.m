function [x_lms, C, Psilms, hatv] = evaluate_xlms(vol, ell_max, Psilms, jball)

hatv = fftshift(fftn(ifftshift(vol)));
data = hatv(jball);
x_lms = cell(ell_max+1, 1);
C = cell(ell_max+1, 1);
for ii=0:ell_max
    if mod(ii, 2)==0
        coeff = ((Psilms{ii+1})'*Psilms{ii+1})\((Psilms{ii+1})' ...
            * real(data));
    else
        coeff = ((Psilms{ii+1})'*Psilms{ii+1})\((Psilms{ii+1})' ...
            * imag(data));
        coeff = 1i * coeff;
    end
    x_lms{ii+1} = reshape(coeff, length(coeff) / (2*ii+1), 2*ii+1);
    C{ii+1} = x_lms{ii+1} * x_lms{ii+1}';
end
end
