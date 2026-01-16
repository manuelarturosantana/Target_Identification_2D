function [usol, scatsol, smoothint, zeroint] = comp_time_sol(ps,lp,f_sols,pols, r_res)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% A poor solution for now, but bk_slow for the residues is computed in this %%%%%%%%%
    %%% file in order to then compute the exp_ints.The SEM values of BK slow will be %%%%%
    %%% Computed seperately                                                         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Using the frequency domain solution this function uses a fourier-series exact integration
    % method to compute the time domain solution. This is steps T0-T5 of the algorithm
    % Inputs:
    %   ps     : The problem structure
    %   lp     : The layer potential object
    %   f_sols : The frequency domain solutions from comp_freq_sols
    %   pols   : The poles of the problem if solving a near singular integration problem
    %   r_res  : The computed residues from comp_r_res. Pass in empty if not using.
    %
    % Outputs:
    %   usol     : A num_spat_points x numt tensor containing the total
    %              field at computed times. If ps.is_far_field, then the incident
    %              field is not added in, and scatsol is the same as usol
    %   scatsol  : A num_spat_points x numt tensor containing just the scattered field at computed times.
    %   smoothint: A num_spat_points x numt tensor containing the integral with the residues subtracted out
    %              if residues are used, otherwise returns empty. Returns just the FC window
    %              portion if zero frequency windowing is used.
    %   zeroint : A num_spat_points x numt tensor containing the zero frequency windowed integral
    %              if zero frequency windowing is used. otherwise returns empty.

    % Check to see if the residuals were passed in for integration. In this implementation 
    % we leave the possibility open that the pols will be needed for some other algorithm.
    use_res = nargin == 5;
    smoothint = [];
    zeroint  = [];

    if ps.has_zero_freq

        [fc_coeffs,k,P] = comp_f_coeffs_fc(ps,f_sols);

        % Compute the weights, which in particular depend on the
        % recentering!
        weights = zeros(ps.npatch, ps.numt, ps.chebN + 1,ps.numk);
        for kind = 1:ps.numk
            weights(:,:,:,kind) = comp_fcc_weights(ps.ts - ps.sk(kind),...
                            ps.hs,ps.npatch,ps.chebN,ps.M_asym, true);
        end

        if use_res
            eint = comp_exp_ints(ps,pols, ps.n_eint);
            [usol, scatsol, smoothint, zeroint] = freq_to_time_zero(ps,lp,f_sols,k,P,fc_coeffs,weights,eint,r_res);
        else
            [usol, scatsol, smoothint] = freq_to_time_zero(ps,lp,f_sols,k,P,fc_coeffs,weights);
        end   
    else
        f_coeffs = comp_f_coeffs(ps,f_sols);
        if use_res
            % Compute the extra needed integral
            % TODO: 3000 is a magic number for now. Find a better way to integrate
            eint = comp_exp_ints(ps,pols, ps.n_eint);
            [usol,scatsol, smoothint] = freq_to_time(ps,lp,f_coeffs,eint,r_res);
        else
            [usol,scatsol] = freq_to_time(ps,lp,f_coeffs);
        end
    end
end

function eints = comp_exp_ints(ps,pols,n)
    % Function which computes all the integrals e^{-1iwt}/(w - p) for adding back in the subtracted out term
    % Inputs:
    %   ps   : The problem structure
    %   pols : The pols of the problem
    % 
    % Outputs:
    %   eints : A numpols x numt x numk matrix containing the computed integrals.
    
    [xj,wj] = fejer(n);
    eints = zeros(length(pols),ps.numt,ps.numk);
    % change integration interval from a,b to [-1,1]
    
    % Subtract out the half residue integral we are integrating. Note it 
    % is not singular at 0.
    if ps.has_zero_freq
        a = 0; b = ps.wlims(2);
    else
        a = ps.wlims(1); b = ps.wlims(2); 
    end
    ws = ((b-a) * xj + a + b) / 2;
    ps2 = ps;
    ps2.ws = ws; ps2.numw = n;
    % bk_slow at the integration points
    bk_slow_full = comp_bkslow(ps2,true); 

    % Compute bk_slow at the poles for the second residue subtraction.
    ps2 = ps; ps2.ws = pols; ps2.numw = length(pols);
    bk_slow_pols = comp_bkslow(ps2, true);

    % TODO : The inner most loop may be able to be vectorized.
    % Set these values so we can do a parfor loop
    numt = ps.numt; numk = ps.numk;
    parfor pind = 1:length(pols)
        for tind = 1:numt
            for kind = 1:numk
                % Don't forget to use the recentering of the fourier transform.
                % ps.t_0 is used for the Gaussian, and should be set to 0 otherwise.
                fs = exp(-1i * ws * (ps.ts(tind) - ps.sk(kind) - ps.t_0));
                BK = bk_slow_full(kind,:); BK = BK.';
                % Residue for this pole
                r_t = bk_slow_pols(kind,pind) * exp(-1i * pols(pind) * (ps.ts(tind) - ps.sk(kind) - ps.t_0));
                
                % Smoothed Integral
                I1 = ((b-a) / 2) * sum(((BK .* fs - r_t) ./ (ws - pols(pind))) .* wj);
                % Exact second integral
                I2 = r_t * (log(b - pols(pind)) - log(a - pols(pind)));

                eints(pind,tind,kind) =  I1 + I2;
            end
        end
    end
    
end
