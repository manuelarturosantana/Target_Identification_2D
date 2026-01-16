function f_coeffs = comp_f_coeffs(ps, f_sols)
    % Computes the fourier coefficients of the frequency domain solution.
    % Inputs:
    %   ps  :   The problem structure
    %   f_sols: The values of the f solution
    %   pols  : The computed poles
    % Outputs :
    %   f_coeffs: The fourier series coefficients
    f_coeffs = zeros(size(f_sols));
    numk = ps.numk;
    numw = ps.numw;
    
    parfor sind=1:ps.num_spat_pts
        for kind = 1:numk
            omega_vals = f_sols(sind,:,kind);
            f_coeffs(sind,:,kind) = (1/numw) * fftshift(fft(ifftshift(omega_vals)));
        end
    end
end