function [x,h] = fc_grid(x_a,x_b,n)
    % Function to generate the grid for the fourier continuation
    % Inputs:
    %   x_a, x_b end points of the interval
    %   n: Number of points to use
    % Outputs:
    %   x: The grid points.
    %   h: The step size of the grid.
       h = (x_b - x_a)/(n-1);
       x = linspace(x_a, x_b, n).';
end