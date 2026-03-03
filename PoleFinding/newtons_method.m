function [x_new, f_val, exit_flag,  iters] = newtons_method(x0, func, dfunc, tol, max_iters, xlims)
    % Function to do newtons method for finding the minimum of a 1-d
    % function.
    % Inputs: 
    %   x0   : An array of the starting guesses.
    %   func : The function to find the zero of
    %   dfunc: The derivative of the func
    %   xlims: A vector with 2 elements representing the brackets of the root
    %   tol  : Algorithm terminates if func(x) < tol
    %   max_iters: Algorithm terminates after max_iters iterations
    %  Outputs:
    %   x_new: The best approximation for the root. Will return nan if the initial guess caused newton's method to leave the bounds.
    %   f_val: The function value, will return nan if the initial guess caused newton's method to leave the bounds.
    %   iters: How many iterations the algorithm was run for
    %   exit_flag: 0 If a root was found, 1 if no root was found
   

    % values to store the approximations of the minimums
    minimums = []; f_vals = []; run_iters = [];

    % Loop through the grid of initial guesses
    for ii = 1:length(x0)
    
        iters = 0;
        f_val = inf;
        x_old = x0(ii);
    
        while (iters < max_iters) && (abs(f_val) > tol)
            x_new = x_old - func(x_old) / dfunc(x_old);
    
            % If we have left the bounds, exit.
            if (x_new < xlims(1)) || (x_new > xlims(2))
                x_new = nan;
                f_val = nan;
                break
            end
    
            iters = iters + 1;
            
            x_old = x_new;
            f_val = func(x_new);
        end
    
        % Save the iterations
        minimums  = [minimums, x_new];
        f_vals    = [f_vals, f_val];
        run_iters = [run_iters, iters];

    

    end

    % Choose the best guess for the root Note that nan will not show up in
    % the min unless everything is nan.
    [~,I] = min(f_vals);
    f_val = f_vals(I);
    x_new = minimums(I);
    iters = run_iters(I);

    exit_flag = isnan(f_val);

    if iters == max_iters
        warning("Newtons Method reached max iterations")
    end

end




