function [Ls,Us] = comp_pole_LU(ps,lp, poles)
    % Function which compute the LU of the boundary integral mat for
    % computing the residue of the poles
    % Inputs:
    %   ps: The problem structure
    %   lp: The layer potential object
    %   poles: The computed poles
    % Outputs:
    %   Ls, Us - Cell array of cell arrays containing the output

    Ls = cell(length(poles),1);
    Us = cell(length(poles),1);

    parfor pind = 1:length(poles)
        circ = poles(pind) + ps.cont_rad * exp(1i * 2*pi *(1:ps.n_res_w) / ps.n_res_w);
        [L_loc,U_loc] = comp_bie_LU(circ,lp);
        Ls{pind} = L_loc; Us{pind} = U_loc;
    end
end