function r_res = comp_r_res(ps,lp,pols,res_dens)
    % Computes the residue at a given point r = (x,y) of the function evalued from the
    % density.
    % Inputs: 
    %   ps       : The problem structure
    %   lp       : The layer potential object
    %   pols     : The computed pols
    %   res_dens : The densitites precomputed for computing the contour integration densities.
    % Outputs:
    %   r_res : a numx x numy x num pols matrix containing the residue computation
    %               for each x,y,pol

    r_res = zeros(ps.numx, ps.numy, length(pols));
    % Jacobian for contour integration.
    eit = exp(1i * 2*pi *(1:ps.n_res_w) / ps.n_res_w);
    dt  = 2 * pi / ps.n_res_w;

    % For loop indicies must be variables, note in a structure for par for.
    numx = ps.numx; numy = ps.numy;
    parfor pind = 1:length(pols)
    % for pind = 1:length(pols)
        % Create the circlular contour around each pole
        circ = pols(pind) + ps.cont_rad * exp(1i * 2*pi *(1:ps.n_res_w) / ps.n_res_w);
        % Function to cancel the error in the first integral
        f2 = @(z) 1 ./ (z - pols(pind));
        % Compute the value of the density at every point
        for xind = 1:numx
            for yind = 1:numy
                if ps.is_open_curve
                    tt = Test_Distance(lp.curve,ps.xs(xind),ps.xs(yind));
                else
                    tt = lp.curve.test_distance(ps.xs(xind),ps.xs(yind));
                end

                if (tt == 1)
                fvals = [];
                % Evaluate the k independent solution
                for cind = 1:length(circ)
                    if ps.is_open_curve
                        fvals(cind) = Scattered_Field(res_dens(:,cind,pind),[],ps.xs(xind),ps.ys(yind),...
                                lp.curve,circ(cind),'Dirichlet', 'None');
                    else
                        fvals(cind) = lp.sol_rep_mat(circ(cind),[ps.xs(xind);ps.ys(yind)]) * res_dens(:,cind,pind);
                    end
                end
                fvals_k = fvals(:);
                % Contour integration computation of the residue. We use the second kind 
                % formulation only because it seems to help some near machine percision.
                r_res(xind,yind,pind) = 1i *dt * sum(fvals_k .* ps.cont_rad .* eit(:))...
                        / (1i *dt * sum(f2(circ) .* ps.cont_rad .* eit));
                end % if the point is outside the curve.
            end % ys
        end % xs
    end % Poles
end