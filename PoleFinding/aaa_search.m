function [pole, num_evals] = aaa_search(bd, f_vals, lp, vec1, vec2, tol, max_ref_steps, ncirc_p, x0)
    % Function to do the search local search using AAA-R or AAA-S. Requires aaa installed 
    % from the chebfun package.
    % Inputs:
    %   bd            : Structure containing the bounds from comp_bounds.m
    %   f_vals/num_xy : Function values from the grid todo the inital AAA search on. Alternatively
    %                   pass in a positive scalar to generate new points in the bounds.
    %   lp            : Layer potential for passing into eval_grid
    %   vec1,vec2     : Random vectors for ri_eval.m
    %   tol           : Stopping tolerance. Evaluated by number of digits for convergence  
    %                   of AAA, or  stopping tolerance on step change of secant method
    %   num_xy        : Number of points in the x and y direction to do grid search for 
    %                   AAA: Try 3 or 4
    %   max_ref_steps : Max number of refinement steps if doing the circle refinement. Pass
    %                   in [] if using secant method Try 3 or 4
    %   ncirc_p       : Number of points to evaluate on the circle for refinement step; Try 3 
    %   x0            : If passed in use the secant method with this as one of the inital guesses.
    %                   usually the coarse grid local minima is used.
    % 
    % Outputs:
    %   pole      : The computed pole
    %   num_evals : The total number of function evaluations done.
    %
    %   Examples: 
    %       pole = aaa_search(bd, f_vals, lp, vec1, vec2, tol, [] 1, 3) use AAA refinement with the grid points.
    %       pole = aaa_search(bd,3, lp, vec1, vec2, tol, [], [], x0) use secant method to converge, with a 3x3 grid of xy points


    % Check if values were passed in for num_xy
    if isscalar(f_vals)
        num_xy = f_vals;
    else
        % There will always be a 3x3 grid around the local minimum on the coarse grid
        num_xy = 3;
    end

    x_c   = linspace(bd.xbds(1), bd.xbds(2), num_xy);
    y_c   = linspace(bd.ybds(1), bd.ybds(2), num_xy); 

    if isscalar(f_vals)
        [~, gridc] = eval_grid(x_c,y_c, lp, vec1, vec2); igridc = gridc.^-1;
    else
        igridc = f_vals.^-1;
    end
    
    xsc   = repmat(x_c,num_xy,1); ysc = repmat(y_c,num_xy,1).';
    [pole, err, poles]  = aaa_pole(igridc, xsc, ysc, bd);

    if abs(err) < tol && isscalar(f_vals)
        num_evals = num_xy^2;
        return
    end

    f = @(x) ri_eval_wrapper(x, lp, vec1, vec2);

    if ~isempty(max_ref_steps)
        delta_x = 2 * pi / ncirc_p;
        t_vals = 0:delta_x:2*pi - delta_x; 
        circ_vals = circ_fun(t_vals);

        ii = 0;
        xsc = xsc(:).';
        ysc = ysc(:).';
        igridc = igridc(:).';
        
        for ii = 1:max_ref_steps
            % Recursive step to see if we can get more digits
            % Add points in a circle around the pole approximant on the order 
            % of accuracy of the pole error.
            
            point = [real(pole); imag(pole)];
            points = point + err * circ_vals;

            xsf = [xsc, real(pole)]; ysf = [ysc, imag(pole)]; fsf = [igridc, f(pole)];
            % TODO: If we are really slow we can parallelize here!
            for jj = 1:ncirc_p
                x = points(1,jj); y = points(2,jj); f_j = f(x + 1i * y);
                xsf = [xsf, x]; ysf = [ysf, y]; fsf = [fsf, f_j];
            end

            % We search for the poles of Ainv, not the zeros of(1 / Ainv)
            [pole2, err] = aaa_pole(fsf.^-1, xsf, ysf, bd);
            if abs(err) < tol
                pole = pole2;
                break
            end
            pole = pole2;
        end
        % The 1 is from evaluating f at the pole
        num_evals = (1 + ncirc_p) * ii;
        if isscalar(f_vals)
            num_evals = num_evals + num_xy^2;
        end
    else
        [pole, ~, its] = secant_method(f, x0, pole, 100, tol, 1);
        if isscalar(f_vals)
            num_evals = its + num_xy^2;
        else
            num_evals = its;
        end
    end
end

function val =  ri_eval_wrapper(x, lp, vec1, vec2)
    % Function to wrape around ri_eval for the zero/pole finding
    [~, val] = ri_eval(real(x), imag(x),lp, vec1, vec2);

end

function [pole, err, poles] = aaa_pole(igrid, xs, ys, bd)
    % Function which computes the single pole of the function given
    % Inputs:
    %   igrid: A grid from eig_grid .^-1 so that it has poles
    %   xs, ys: real and imaginary part of grid data. Size must match igrid.
    %   bd    : The bounds structure for this root
    % Output: 
    %   pole  : The best approximation for the pole using the aaa algorithm.
    %   err   : The difference between the last 
    %   poles : All poles computed from the AAA approximation

    % Here we make a big assumption: If there is only one pole of the actual function the 
    % rational function approximation of degree one has one pole which is an approximation 
    % of the actual pole of the function.
    [~, pole] = aaa(igrid, xs + 1i * ys, 'degree', 1);
    sp = [pole];
    for deg = 1:10
        [~, poles] = aaa(igrid, xs + 1i * ys, 'degree', deg);
        % Make sure we are getting the pole with in the bounds of the coarse grid.
        poles = poles(isinbounds(poles, bd));
        [~,I] = min(abs(sp(end) - poles));
        % Check for convergence, and if it happens break before finishing 
        % recording the new pole so we can get a convergence tolerance.
        if deg > 1 &&  abs(poles(I) - sp(end)) < 1e-12
            sp = [sp; poles(I)];
            break
        end
        sp = [sp; poles(I)];
    end
    
    
    pole = sp(end);
    
    if deg == 2
        err = abs(poles(I) - pole);
    else
        % Because we record that convergence happened, we need to compare
        % the step before convergence.
        err  = abs(sp(end) - sp(end - 2));
    end
end

function bool = isinbounds(x, bd)
    % Given a bound object return true if x in the bounds
    bool1 = real(x) > bd.xbds(1) & real(x) < bd.xbds(2);
    bool2 = imag(x) > bd.ybds(1) & imag(x) < bd.ybds(2);

    bool = bool1 & bool2;
end

function vals = circ_fun(t)
    vals = [cos(t); sin(t)];
end

