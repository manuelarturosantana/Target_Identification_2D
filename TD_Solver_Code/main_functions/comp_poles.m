function [pols,errs,nfevals] = comp_poles(ps, lp, plot_poles)
    % Function which computes the poles using the rational approximation method.
    % Inputs:
    %     ps : The problem structure
    %     lp : The layer potential object.
    %     plot_poles: (Optional) if true, plot the poles. Default false
    % Outputs:
    %     pols: The computed poles by the AAA/Secant Algorithm
    %     errs: Vector of errs for the computed poles. Will be [] if secant method is not used.
    %     nfevals: The number of function evaluations used to find the poles.
    %
    % Warning: Currently the second guess for secant method is pol + 1e-5. This may or may not work for your situation.

    if nargin == 2
        plot_poles = false;
    end

    % Initialize scalarized resolvent.
    N = size(lp.bie_mat(1),1);
    u = rand(1,N,'like',1i); v = rand(N,1,'like',1i);
    f = @(w) u * (lp.bie_mat(w) \ v);
    % f = @(w) sr(lp,u,v,w);

    nfevals = nan;
    % If zero is in the limits we only compute the poles along the left
    % half.
    if ps.wlims(1) < 0 && ps.wlims(2) > 0
        rlims = [0,ps.wlims(2)];
    else
        rlims = [ps.wlims(1), ps.wlims(2)];
    end

    rp = rec_params('xlims',rlims,'ylims',[ps.wimag,0],'ntb',ps.p_numy,'nlr',ps.p_numx,"plotls",plot_poles);
    pols = aaa_recursive(f,rp);
    

    if plot_poles
        figure()
        clf
        hold on
        plot(pols,'*')
        hold off
    end
        
    % If computing power allows, refine each pol approximation with a secant method or localized AAA
    if ps.use_secant
        [pols, errs, nfevals_s] = sec_ref(f,pols,ps.sec_its,ps.sec_tol);
        if plot_poles
            hold on
            plot(pols,'^')
            hold off
        end
        nfevals = nfevals + nfevals_s;
    else
        errs = [];
    end

    
end