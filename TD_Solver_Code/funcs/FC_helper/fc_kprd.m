function [k,prd] = fc_kprd(n,h)
    % Function which computes the wave vectors and the period of the fc series with n points.
    % Written to save these outside of par for loop
    % Inputs:
    %   n   : Length of input series to FC
    %   h   : Distance between points in fc
    % Outputs:
    %   k     : The vector of wave numbers.
    %   per   : The period of the function

    C = 27; % number of continuation points.
    fourPts = n + C; % Number of points in the extended grid
    prd = fourPts*h; % Extended period
    if (mod(fourPts, 2) == 0) % wave vector
        k = transpose([0:fourPts/2, -fourPts/2+1:-1]);
    else
        k = transpose([0:(fourPts-1)/2, -(fourPts-1)/2:-1]);
    end
end