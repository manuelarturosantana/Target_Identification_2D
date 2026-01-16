function vals = fc_eval(z,fc_coeffs,k,prd,x_a)
    % Function to evaluate the fourier continuation

    % Force Z to be a row vector.
    z = z(:); z = z.';
    
    vals = real(exp(2*pi*1i*(z - x_a).'*k.'/prd) * ...
            fc_coeffs(:));
end