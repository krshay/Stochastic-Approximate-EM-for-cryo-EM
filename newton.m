
%--------------------------------------------------------------------------
function x = newton(fun,deriv,nu,x0)
% Newton-Raphson for Bessel function FUN or FUN' with initial guess X0
x = x0;
c = 8;
for t = 1:10
    % Newton-Raphson step
    f = fun(nu,x);
    g = nu*f./x;
    df = fun(nu-1,x) - g;
    if deriv == 0
        h = -f./df;
    else
        ddf = (nu*g - df)./x - f;
        h = -df./ddf;
    end
    x = x + h;
    % Convergence criteria
    if all(abs(h) < c*eps(x))
        break
    end
    % Relax convergence criteria
    if t >= 7
        c = 2*c;
    end
end
% Check convergence
if t == 10
    warning('ZEROBESS:Newton','No convergence for Newton-Raphson.')
end
% fprintf('%2d%8d\n',t,numel(x))