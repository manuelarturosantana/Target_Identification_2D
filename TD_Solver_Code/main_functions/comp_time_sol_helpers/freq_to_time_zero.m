function [usol, scatsol, smoothint, zeroint] = freq_to_time_zero(ps,lp,f_sols,k,P,fc_coeffs,weights,eint,r_res)
    % Performs the integration for each point to go from the frequency to the the time domain
    % when there is zero frequency content
    % Inputs:
    %     ps: The problem structure
    %     lp: The layer potential object
    %     fsols : The frequency domain solutions
    %     k,P: The frequencies and period from fourier continuation.
    %     fc_coeffs : The frequency fourier coefficients of the spacial solution.
    %     weights   : The weights for the Filon-Clenshaw-Curtis rule
    %     eint      : Optional: The computed exponential integrals to subtract out. 
    %     r_res     : Optional: The residuals at each point r.
    % Outputs:
    %   same as comp_time_sol.m
    % 
    % Usage:
    %  [usol, scatsol] = freq_to_time(ps,lp,f_coeffs) No residue subtraction case without zero freqeuncy windowing
    %  [usol, scatsol, smoothint] = freq_to_time(ps,lp,f_coeffs,eint,r_res) Residue subtraction case without zero freqeuncy windowing

    use_res = (nargin >=8);
    usol = zeros(ps.num_spat_pts,ps.numt);
    scatsol = zeros(ps.num_spat_pts,ps.numt);
    if use_res
        nump = size(eint,1);
        % Edge case for if AAA_SV is not resolved enough to find any poles.
        % Then we must set the parameters so that the parfor loop still
        % works.
        if nump == 0
            nump = 1; eint=zeros(nump,ps.numt,ps.numk); r_res=zeros(ps.num_spat_pts,nump);
        end

        smoothint = zeros(ps.num_spat_pts,ps.numt);
        zeroint = zeros(ps.num_spat_pts,ps.numt);
    else % These initializations are to make the parfor parser happy
        nump = 1;
        smoothint = zeros(ps.num_spat_pts,ps.numt);
        zeroint = zeros(ps.num_spat_pts,ps.numt);
        eint=zeros(nump,ps.numt,ps.numk); r_res=zeros(ps.num_spat_pts,nump);
    end

    % Endpoints for inverse fc integration
    fc_a = ps.ws(ps.gmend + 1); fc_b = ps.ws(end);

    % For loop indicies must be variables, not in a structure for par for.
    numt = ps.numt;  numk = ps.numk; 
   
    parfor sind=1:ps.num_spat_pts
            % We test the distance so points inside the curve don't get evaluated.
            if ps.is_far_field
                tt = 1
            elseif ps.is_open_curve
                tt = Test_Distance(lp.curve,ps.xs(sind),ps.ys(sind));
            else
                tt = lp.curve.test_distance(ps.xs(sind),ps.ys(sind));
            end
            if (tt == 1)

            % Compute the fourier continuation integrals, add in the exponential integrals
            for tind=1:numt
                for kind = 1:numk
                    % Compute the integrals for the non-zero frequency
                    % content
                    cms  = fc_coeffs(sind,:,kind);
                    cms = cms(:).';
                    % Here per equation 23 we evaluate at the t talue t - sk
                    ksol_fc = ft_fc(cms,k,P,fc_a,fc_b,ps.ts(tind) - ps.sk(kind));
                    
                    % Compute the zero frequency integrals.
                    zer_freq_intsk = gfcc_int(f_sols(sind,1:ps.gmend,kind),...
                        ps.ts(tind) - ps.sk(kind),ps.ccps,ps.hs, ps.cs,weights(:,tind,:,kind),ps.npatch);
                    ksol = ksol_fc + zer_freq_intsk;
                    
                   
                    if use_res
                        if ps.do_win_zero
                            zeroint(sind,tind)   = zeroint(sind,tind) + zer_freq_intsk;
                            smoothint(sind,tind) = smoothint(sind,tind) + ksol_fc;
                        else
                            smoothint(sind,tind) = smoothint(sind,tind) + ksol;
                        end

                       
                        % Add in the exp integral for each term
                        for pind = 1:nump
                            ksol = ksol + r_res(sind,pind) * eint(pind,tind,kind);
                        end
                    end
                    % Don't forget the fourier transform normalization
                    % term after adding back the residues!
                    % Also multiply by 2 because it is the inverse fourier transform of
                    % a real function.
                    ksol = 1/(pi) * ksol;
                    usol(sind,tind) = usol(sind,tind) + ksol;
                end % ks
                % Saving just the scattered solution.
                scatsol(sind, tind) = usol(sind,tind);
                % Add in the incident wave
                if ~ps.is_far_field
                    usol(sind,tind)   = usol(sind, tind) + ...
                    ps.uinc(ps.xs(sind),ps.ys(sind),ps.ts(tind));
                end
            end %ts
            
            end % If the point is outside the curve.
    end % sind
end % Function.