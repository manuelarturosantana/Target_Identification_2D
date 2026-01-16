function f_sols = comp_freq_sols(ps,lp, densities, Bkslow, pols, r_res)
    % Make a solution for each x,y and omega value needed.
    % Inputs: 
    %   ps  : The problem structure
    %   lp  : The layer potential object
    %   densities : The densities precomputed for each frequency point
    %   Bkslow    : The slowly varying frequency term
    %   pols      : Optional: The computed poles
    %   r_res     : Optional: The residue at each pole
    % Outputs:
    %   f_sols   : a numx x numy x numw x numk tensor containing the solution at each point.
    %
    % Usage examples
    %  f_sols = comp_freq_sols(ps,lp, densities, Bkslow) when no poles are computed.
    %  f_sols = comp_freq_sols(ps,lp, densities, Bkslow, pols, r_res) when poles and residues are computed;
    %  Any other usage is undefined.

    has_pols_res = (nargin >= 5);

    f_sols = zeros(ps.numx,ps.numy,ps.numw,ps.numk);
    
    % Define r_res and pols so the parfor loop knows what is happening.
    if ~has_pols_res
        r_res = 0; pols = 0;
    end

    % For each point in r
    % For loop indicies must be variables, note in a structure for par for.
     numk = ps.numk; numy = ps.numy; numw = ps.numw;
     % parfor xind = 1:ps.numx
    parfor xind = 1:ps.numx
        xx=ps.xs(xind);
        for yind = 1:numy
            yy=ps.ys(yind);
            if ps.is_open_curve
                tt = Test_Distance(lp.curve,xx,yy);
            else
                tt = lp.curve.test_distance(xx,yy);
            end

            if (tt == 1)
            for wind = 1:numw
                if ps.is_open_curve
                    solw = Scattered_Field(densities(:,wind),[],xx,yy,lp.curve,...
                            ps.ws(wind),'Dirichlet', 'None');
                else
                    solw = lp.sol_rep_mat(ps.ws(wind),[xx;yy]);
                    solw = solw * densities(:,wind);
                end
                
                % Subtract the k independent residue out of the frequency
                % solution before windowing
                if has_pols_res 
                    for pind = 1:length(pols)
                        solw = solw - (r_res(xind,yind,pind) ./ (ps.ws(wind) - pols(pind)));
                    end
                end
                
                for kind = 1:numk
                    f_sols(xind, yind, wind,kind) = Bkslow(kind,wind) * solw;
                end % ks
            end %ws
            end % if tt==1 
        end % ys
    end % xs
end % function