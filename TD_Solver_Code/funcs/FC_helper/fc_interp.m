function [fc_coeffs, k, prd] = fc_interp(fx, h, d,A,Q,Q_tilde)
    % Function which performs the FC interpolation
    % Inputs:
    %    fx, h should be created based off the following:
    %    n = 200; % number of points in the initial domain
    %    x_a = -0.6; x_b = 1; % The beginning and end of the Cartesian grid
    %    h = (x_b - x_a)/(n-1);
    %    x = linspace(x_a, x_b, n).';
    %    fx = f(x)
    %
    %    d: is the number of Gram polynomial interpolation points.
    %       Default is 10. Any number other than 5,10,12 will require precomps to happen. 
    %
    %    A,Q,Q_tilde come from the commented out code. It is faster and more
    %    robust to load them before.
    %
    % Outputs:
    %   coeffs: The coefficients on the FC
    %   k     : The vector of wave numbers.
    %   per   : The period of the function
    %   f_vals: The function values

    
    n = length(fx);

    if nargin == 2
        d =  10; % Number of Gram polynomial interpolation points
    end

    C = 27; % Number of continuation points
    
    Z = 12;
    E = C;
    n_over = 20;
    modes_to_reduce = 4;
    num_digits = 256;
    
    fourPts = n + C; % Number of points in the extended grid
    prd = fourPts*h; % Extended period
    if (mod(fourPts, 2) == 0) % wave vector
        k = transpose([0:fourPts/2, -fourPts/2+1:-1]);
    else
        k = transpose([0:(fourPts-1)/2, -(fourPts-1)/2:-1]);
    end
    % 
    % load(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat']);
    % load(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat']);
    % load(['FC_data/Q_tilde_d',num2str(d),'_C', num2str(C), '.mat']);
    % 
    % % [ArQr, AlQl, ArQ_tilder, AlQ_tildel] = build_cont_mat(A, Q, Q_tilde);
    % A = double(A);
    % Q = double(Q);
    % Q_tilde = double(Q_tilde);
    
    % Derivative spectral coefficients
    der_coeffs0 = 1;
    filter_coeffs = 1;
      
    BC0 = [0, 0; 0, 0];
    
    [~, fc_coeffs] = fc_der(fx, der_coeffs0, filter_coeffs, d, C, A, Q, Q, BC0, h);
end