function bkslow = comp_bkslow(ps, is_only_FC_win)
    % Function to compute the slowly oscillatory part of the transform of the incident signal.
    % Corresponds to Step F2 in the paper
    % Inputs:
    %   ps    : The problem structure
    %   is_only_FC_win : Optional input for the case for problems windowing around zero frequency.
    %           If true, but assume only the FC window is used, which is true for the residue computation,
    %           ,so that the residues are correct. default: false
    %           
    % Outputs:
    %   Bkslow: matrix of size num_windows x num_freq_solutions  

    if nargin == 1
        is_only_FC_win = false;
    end
    bkslow = zeros(ps.numk, ps.numw);

    % Create t_values for Bkslow, and take off the end point due to periodicity.
    % Also M must be even for this to work! (This M just gives us enough t
    % values to do the fourier transform for bkslow and is not connected to t
    % otherwise).
    M  = 500;
    ms = (-M/2:M/2-1);
    dx = (2*ps.H) / M;
    t_vals = ms * dx;

    % For loop indices must be variables, not in a structure
    numk = ps.numk;
    parfor wind = 1:ps.numw
    % for wind = 1:ps.numw
        for kind = 1:numk
            % TODO: Everything is coded up so we can do the windowing strategy,  
            % except in how the extra integral term is computed when adding back the residue.
            % However 
            % for simplicity we will start with just using the exact fourier transform 
            % and therefore only needing one window.
            if strcmp(ps.inc_field,'gaussian')
                bkslow(kind,wind) = ps.B_w(ps.ws(wind));
            else
                % TODO: Code for windowing strategy
                % Following equations 14 and 15 the window function is
                % defined as w_k = w(t - sk) and the change which we then
                % evaluate at w_k(t + sk) = w(t)
                ak = ps.a_t(t_vals + ps.sk(kind)) .* ps.pou(t_vals, ps.H);
                bkslow(kind, wind) = ft(ak,-ps.H,ps.H,ps.ws(wind),false);

                % % TESTING: check the integrals are computed correctly.
                % akt = @(t) ps.a_t(t + ps.sk(kind)) .* ps.pou(t,ps.H);
                % test = integral(@(t) akt(t) .* exp(1i * ps.ws(wind) *t ),-ps.H,ps.H, 'AbsTol',1e-14);
                % abs(test - bkslow(kind,wind))
            end
        end
        % Test case that everything adds up!
        % bkslow2 = bkslow(:,wind).* exp(1i * ps.sk.' * (ps.ws(wind).'));
        % % at = @(t) ps.a_t(t);
        % % test = integral(@(t) at(t) .* exp(1i * ps.ws(wind) *t ),ps.ts(1),ps.ts(end), 'AbsTol',1e-14);
        % abs(ps.B_w(ps.ws(wind)) - sum(bkslow2))
    end

    % Multiply by the window for both integrals.
    if ps.has_zero_freq && ps.do_win_zero
        % Case where we are evaluating the window at the poles, so we only use the FC window
        % 
        if is_only_FC_win
            win_FC = ps.win_FC(ps.ws);
            % Matrix .* row vector multiplies each element of each row by the row vector.
            win_FC = win_FC(:).';
            bkslow = bkslow .* win_FC; 
        else
            win_z  = ps.win_zero(ps.ws(1:ps.gmend));   win_z = win_z(:).';
            win_FC = ps.win_FC(ps.ws(ps.gmend+1:end)); win_FC = win_FC(:).';

            bkslow(:,1:ps.gmend)     = bkslow(:,1:ps.gmend) .* win_z;
            bkslow(:,ps.gmend+1:end) = bkslow(:,ps.gmend+1:end) .* win_FC;
        end
     
    end
end