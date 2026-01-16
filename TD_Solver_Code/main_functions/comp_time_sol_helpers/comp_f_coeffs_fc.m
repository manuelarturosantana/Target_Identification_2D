function [f_coeffs,k,P] = comp_f_coeffs_fc(ps, f_sols)
    % Computes the fc coefficients of the frequency domain solutions.
    % Inputs:
    %   ps  :   The problem structure
    %   f_sols: The values of the f solution
    % Outputs :
    %   f_coeffs: The fourier continuation coefficients of size numx x numy x numw x numk
    %   k     : vector of wave numbers from fc
    %   P     : Period of fc sum
    %
    % WARNING: This implementation is not general! It only works for
    % for the integral from w_c to W in the zero frequency case.

    NCONT_POINTS = 27;
    % The distance in the FC grid.
    h = ps.ws(ps.gmend+2) - ps.ws(ps.gmend+1);

    f_coeffs = zeros(ps.num_spat_res, ps.numwFC + NCONT_POINTS, ps.numk);
    numk = ps.numk;
    
    % Was running into a strange bug with the FC code, so now the matrices
    % are loaded before the loop, this saves time as these matrices are
    % small.
    d = 10; C = 27;
    load(['FC_data/A_d',num2str(d),'_C', num2str(C), '.mat'],'A');
    load(['FC_data/Q_d',num2str(d),'_C', num2str(C), '.mat'],'Q');
    load(['FC_data/Q_tilde_d',num2str(d),'_C', num2str(C), '.mat'],'Q_tilde');

    A = double(A);
    Q = double(Q);
    Q_tilde = double(Q_tilde);
    
    parfor sind=1:ps.num_spat_res
        for kind = 1:numk
                % Only grab the points after the graded mesh frequencies
                omega_vals = f_sols(sind,ps.gmend+1:end,kind);
                % Smoothen by subtracting the residues
                omega_vals = omega_vals(:);
                [fc_coeffs, ~, ~] = fc_interp(omega_vals,h,10,A,Q,Q_tilde);
                f_coeffs(sind,:,kind) = fc_coeffs;
        end
    end
    % Compute the wave vector and the period for the fist FC integral.
    [k, P] = fc_kprd(length(f_sols(1,1,ps.gmend+1:end,1)),h);
end