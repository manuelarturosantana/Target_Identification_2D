function [usol, scatsol, smoothint] = freq_to_time(ps,lp,f_coeffs,eint,r_res)
    % Performs the integration for each point to go from the frequency to the the time domain.
    % Inputs:
    %     ps: The problem structure
    %     lp: The layer potential object
    %     f_coeffs: The frequency fourier coefficients of the spacial solution.
    %     eint    : Optional: The computed exponential integrals to subtract out. 
    %     r_res   : Optional: The residuals at each point r.
    % Outputs:
    %   same as comp_time_sol.m
    % 
    % Usage:
    %  [usol, scatsol] = freq_to_time(ps,lp,f_coeffs) No residue subtraction case
    %  [usol, scatsol, smoothint] = freq_to_time(ps,lp,f_coeffs,eint,r_res) Residue subtraction case

    use_res = (nargin >=4);
    usol = zeros(ps.numx,ps.numy,ps.numt);
    scatsol = zeros(ps.numx,ps.numy,ps.numt);
    if use_res
        nump = size(eint,1);
        % Edge case for if AAA_SV is not resolved enough to find any poles.
        % Then we must set the parameters so that the parfor loop still
        % works.
        if nump == 0
            nump = 1; eint=zeros(nump,ps.numt,ps.numk); r_res=zeros(ps.numx,ps.numy,nump);
        end
        smoothint = zeros(ps.numx,ps.numy,ps.numt);
    else
        % The parfor loop will attempt to classify each variable by looking
        % at endpoints of the loop, so we must allocate something of small
        % enough size.
        nump = 1;
        smoothint = zeros(ps.numx,ps.numy,ps.numt);
        eint=zeros(nump,ps.numt,ps.numk); r_res=zeros(ps.numx,ps.numy,nump);
    end

    % For loop indicies must be variables, note in a structure for par for.
    numt = ps.numt;  numk = ps.numk; numy = ps.numy; 

    parfor xind=1:ps.numx
    % for xind=1:ps.numx
        for yind=1:numy
            % We test the distance so points inside the curve don't get evaluated.
            if ps.is_open_curve
                tt = Test_Distance(lp.curve,ps.xs(xind),ps.xs(yind));
            else
                tt = lp.curve.test_distance(ps.xs(xind),ps.ys(yind));
            end
            if (tt == 1)
            for tind=1:numt
                    for kind = 1:numk
                        cms  = f_coeffs(xind,yind,:,kind);
                        cms = cms(:).';
                        % Here per equation 23 we evaluate at the t talue t - sk
                        % Also we evaulate at the shifted time for the incident field.
                        ksol = ft(cms,ps.wlims(1), ps.wlims(2),...
                                ps.ts(tind) - ps.sk(kind) - ps.t_0, true, true);
                        if use_res
                            smoothint(xind,yind,tind) = smoothint(xind,yind,tind) + ksol;
                            % Then add in the exp integral for each term
                            for pind = 1:nump
                                ksol = ksol + r_res(xind,yind,pind) * eint(pind,tind,kind);
                            end
                        end
                        % Don't forget the fourier transform normalization
                        % term after adding back the residues!
                        ksol = 1/(2*pi) * ksol;
                        usol(xind,yind,tind) = usol(xind,yind,tind) + ksol;
                    end % ks
                    % Saving just the scattered solution.
                    scatsol(xind, yind, tind) = usol(xind,yind,tind);
                    % Add in the incident wave
                    usol(xind, yind, tind)   = usol(xind, yind,tind) + ...
                    ps.uinc(ps.xs(xind),ps.ys(yind),ps.ts(tind));
            end %ts
            end % If the point is outside the curve.
        end %y ind
    end % xind
end % Function.