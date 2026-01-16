function densities = comp_dens(ps,lp)
% Function which computes the frequency domain densities
% Corresponds to step F3 in the paper.
% Inputs:
%   ps  : The problem stucture
%   lp  : The layer potential object
%
% Output:
%   densities a matrix with densities for each frequency in each column

    densities = [];

    if ps.is_open_curve
        bndpts = [lp.curve.X, lp.curve.Y]';
    else
        bndpts = lp.curve.xs;
    end

    %   - For each omega solve the BIE with the incident plane wave. 
    parfor wind = 1:ps.numw
    % for wind = 1:ps.numw
        rhs = -exp(1i * ps.ws(wind) * (ps.kappa * bndpts)).';
        den = lp.bie_mat(ps.ws(wind)) \ rhs;
        densities = [densities, den];
    end

end