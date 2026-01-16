function full_SEM = sum_SEM(ps,pols,r_res,xinds,yinds)
    % Compute the singularity expansion for all times at the given indicies
    % in x and y
    % Inputs:
    %   ps    : The problem perameters
    %   pols  : The computed poles
    %   r_res : The numx x numy x num_pols array of spatial residues
    %   xind  : A vector of x indicies to compute the expanion at
    %   yind  : A vector y indicies to compute the expanion at
    % Outputs :
    %   fullSEM: A length(xind) x length(yind) x ps.numt matrix containing
    %            all the singularity expansion terms.

    
    full_SEM = zeros(length(xinds),length(yinds),ps.numt);

    % Compute the residue contribution from the incidient field.
    bk_slow_pols = zeros(1,length(pols));
    if strcmp(ps.inc_field,'gaussian')
        for pind = 1:length(pols)
            bk_slow_pols(1,pind) = ps.B_w(pols(pind));
        end
    else
        ps2 = ps;
        ps2.ws = pols; ps2.numw = length(pols);

        bk_slow_pols = comp_bkslow(ps2,true);
        % Multiply by the recentering term to compute the correct residue
        bk_slow_pols = bk_slow_pols .* exp(1i * ps.sk.' * (pols(:).'));
        % We get the residue we need the contribution from B_k as it limits to
        % the pole, and then add all the windows the windows together.
        % Make it a column vector
        bk_slow_pols = sum(bk_slow_pols,1).';
    end
    
    % The r independent complex exponential time term
    % Transposed to have time in the rows, and poles in the columns
    exp_pol = exp(-1i *pols * (ps.ts - ps.t_0)).';

    % parfor loop variables
    numt = ps.numt; numy = length(yinds);
    has_zero_freq = ps.has_zero_freq;
    parfor xx = 1:length(xinds)
        for yy = 1:numy
            % The residue  is a column vector of size num_pols
            res = bk_slow_pols.'.*squeeze(r_res(xinds(xx),yinds(yy),:));
            % The exp_pol is a matrix with time in rows and poles in
            % columns, making this matvec mult a vector of the values for
            % all times
            vals = -1i * exp_pol * res;
            % If there is zero frequency multiply this by 2 to account for the half of the fourier
            % transform
            if has_zero_freq
                vals = vals * 2;
            end
            full_SEM(xx,yy,1:numt) = vals;
        end
    end
end