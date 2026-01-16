function f_coeffs = comp_f_coeffs(ps, f_sols)
    % Computes the fourier coefficients of the frequency domain solution.
    % Inputs:
    %   ps  :   The problem structure
    %   f_sols: The values of the f solution
    %   pols  : The computed poles
    % Outputs :
    %   f_coeffs: The fourier series coefficients
    f_coeffs = zeros(size(f_sols));
    numk = ps.numk; numy = ps.numy; 
    
    parfor xind=1:ps.numx
    % for xind=1:ps.numx
        for kind = 1:numk
            for yind=1:numy
                omega_vals = f_sols(xind,yind,:,kind);
                f_coeffs(xind,yind,:,kind) = (1/ps.numw) * fftshift(fft(ifftshift(omega_vals)));
            end
        end
    end
end