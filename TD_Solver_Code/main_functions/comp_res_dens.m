function [res_dens] = comp_res_dens(ps, lp, pols, L, U)
    % Function which computes the residues of each pole.
    % Inputs: 
    %     ps  : The problem structure
    %     lp  : The layer potential object.
    %     pols: The computed poles from the AAA method.
    %     L,U :  optional Cell arrays cell arrays of LU factors for each frequency  
    %
    % Outputs: 
    %  res_dens : A N x n_res_w x num_pols tensore containing the densities for the contour 
    %              Integration points around each residue
    %  bk_res   : A n_res_w x num_poles x num_k tensor which contains the non-planewave part
    %             of the density
    % Some notes:
    % Convergence studies suggest 2-6 frequency domain points is sufficient for convergence.
    % In specific we are computing the residues associated to poles of U_k^slow (see equation
    % 22 in the paper), which means we can not worry about the highly oscillatory term e^{iskw}

    res_dens  = zeros(ps.N,ps.n_res_w,length(pols));
    
    use_LU = nargin == 5;

    if ps.is_open_curve
        bndpts = [lp.curve.X, lp.curve.Y]';
    else
        bndpts = lp.curve.xs;
    end

    % For loop indicies must be variables, note in a structure for par for.
    n_res_w = ps.n_res_w;
    parfor pind = 1:length(pols)
    % for pind = 1:length(pols)
        % Contour around the pole
        circ = pols(pind) + ps.cont_rad * exp(1i * 2*pi *(1:ps.n_res_w) / ps.n_res_w);
        for wind = 1:n_res_w
            % Create the RHS plane wave modulated by the input signal
            rhs = -exp(1i * circ(wind) * (ps.kappa * bndpts)).';
            % Now compute the density
            if use_LU
                 U_loc = U{pind}{wind}; L_loc = L{pind}{wind}
                 res_dens(:,wind,pind) = U_loc \ (L_loc \ rhs);
            else
                res_dens(:,wind,pind) = lp.bie_mat(circ(wind)) \ rhs;
            end
        end
    end
end