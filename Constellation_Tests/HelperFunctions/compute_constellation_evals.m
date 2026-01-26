function [poles, lp, mc] = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation)
    % Function which builds a constellation and computes its poles.
    % 
    % Inputs:
    %   curve  : The base curve to build the constellation from
    %   xshift : For each curve beyond the first, shift in x for the
    %            location of the new curve.
    %   yshift : Same as xshift but for y
    %   rp     : Parameters for the aaa recursive calculation
    %   draw_constellation : If true plot the constellation before
    %   computing
    %
    %  Outputs:
    %       poles: The computed poles
    %       lp   : The layer potential object
    %       mc   : The multiple curves object

    [lp, mc] = setup_constellation(curve, xshift, yshift, draw_constellation);


    N_total = (length(xshift) + 1) * curve.N;

    v = randn(N_total,1); u = rand(1,N_total);

    f = @(k) u * (lp.bie_mat(k) \ v);

    poles = aaa_recursive(f,rp);

end