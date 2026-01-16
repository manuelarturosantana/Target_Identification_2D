function [sigmas,wlims] = gauss_picker(w_0,dist,tol)
    % Given w_0 calculate value of sigmas so that a gaussian
    % e^(-(w - w_0)^2/sigmas)  is tol at w_0 +- dist
    if (nargin < 3)
        tol = 1e-16;
    end
    sigmas = -dist^2 / log(tol);
    wlims = [w_0 - dist, w_0 + dist];
end