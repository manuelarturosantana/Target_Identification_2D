function vals = erf_hwin(z,A)
% An infinite half window using the error function, which is approximantly zero at A /2,
% and 1 at approximant A. For complex values of z, numerical integration with matlab's
% integral function is used. For this reason the implementation is not robust for |imag(z)| %
% large.
% Inputs:
%   z   : A vector of values to evaluate the function at
%   A   : Window is approximantly zero before A/2, and one after A

    rho   = 5.805; % A/2 approx 1.1x10^-16.
    alpha = 0.5;   % For real z the function will be zero before z < \alpha H

    if isreal(z)
        vals = 0.5 + 0.5 * erf(-rho + 2 * rho * ((z - alpha * A) / ((1 - alpha) * A)));
    else
        vals = [];
        f = @(t) exp(-t.^2);
        for ii = 1:length(z)
            % integral function is replacing the calculation of int_0^z e^(-t^2)
            vals(ii) = 0.5 + 0.5 * (2 / sqrt(pi)) *...
                 integral(f,0,-rho + 2 * rho * ((z(ii) - alpha * A) / ((1 - alpha) * A)),'AbsTol',1e-15,'Reltol',1e-15);
        end
    end

end