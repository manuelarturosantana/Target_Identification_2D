function densities = comp_dens(ps,lp, L, U)
% Function which computes the frequency domain densities
% Corresponds to step F3 in the paper.
% Inputs:
%   ps  : The problem stucture
%   lp  : The layer potential object
%   L,U :  optional Cell arrays of LU factors for each frequency
%
% Output:
%   densities a matrix with densities for each frequency in each column


    use_LU = nargin == 4;


    if ps.is_open_curve
        bndpts = [lp.curve.X, lp.curve.Y]';
    else
        bndpts = lp.curve.xs;
    end

    %   - For each omega solve the BIE with the incident plane wave. 
    ndof = length(bndpts);
    densities = zeros(ndof, ps.numw);
    
    parfor wind = 1:ps.numw
    % for wind = 1:ps.numw
        rhs = -exp(1i * ps.ws(wind) * (ps.kappa * bndpts)).';
        if use_LU == 4
            den = U{wind} \ (L{wind} \ rhs);
        else
            den = lp.bie_mat(ps.ws(wind)) \ rhs;
        end
        
        densities(:, wind) = den;
    end

end