function [sigmas,wlims] = gauss_picker(w_0,dist)
    % Given w_0 calculate value of sigmas so that a gaussian
    % e^(-(w - w_0)^2/sigmas)  is 1e-16 at w_0 +- dist
    sigmas = -dist^2 / log(1e-16);
    wlims = [w_0 - dist, w_0 + dist];
end