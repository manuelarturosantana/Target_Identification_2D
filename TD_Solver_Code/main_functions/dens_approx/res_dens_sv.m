function [pols, res_dens] = res_dens_sv(ps, Rs, pol_c)
    % Inputs: 
    %     ps   : The problem structure
    %     Rs   : The rational approximant of the densities cell array from comp_res_sv
    %     pol_c: The cell array of computed poles corresponding to Rs, and comp_res_sv
    %  res_dens : A N x n_res_w x num_pols tensore containing the densities for the contour 
    %              Integration points around each residue
    %     pols: An array containing all the computed poles.
    %

    
    % Grab all computed poles
    pols = [];
    for rind = 1:length(Rs)
        pols = [pols; pol_c{rind}];
    end

    res_dens = zeros(ps.N,ps.n_res_w,length(pols));

    % For loop indicies must be variables, note in a structure for par for.
    n_res_w = ps.n_res_w; 
    % Index immediately before where the current pole aligns
    pstart = 0;
    for rind = 1:length(Rs)
        R = Rs{rind};
        pols_loc = pol_c{rind};
        parfor (pind = 1:length(pols_loc), ps.nworkers)
            % Contour around the pole
            circ = pols_loc(pind) + ps.cont_rad * exp(1i * 2*pi *(1:ps.n_res_w) / ps.n_res_w);
            for wind = 1:n_res_w
                % Now compute the density
                res_dens(:,wind,pstart + pind) = R(circ(wind));
                % Now compute the window dependent bkslow term.
            end
        end
        pstart = pstart + length(pols_loc);
    end
end