function int = ft_fc(fc_coeffs,k,P,a,b,t)
    % Function which computes the inverse fourier transform using the sinc method and
    % fourier continuation.
    % Inputs:
    %   fc_coeffs : The fc_coeffs computed between a and b
    %   k         : The vector of wave numbers from fc
    %   P         : The period of the fourier continuation
    %   a,b       : The end points of integration
    %   t         : The time to evaulate the fourier transform at
    A     = (b - a) / 2;
    delta = (b + a) / 2;

    % Ensure inputs didn't get switched around.
    fc_coeffs = fc_coeffs(:); k = k(:);
    
    alpha = P / (2 * pi);
    % This constant is modified so we can call sinc, which includes the pi
    const = (2 * A);
    int = fc_coeffs .* const .* sinc((2 * A / P) * (alpha * t - k));
    % We multiply by the shifting factor since FC starts from 0.
    int = sum(int.* exp(-1i *2*pi* k * (-A) / P));
    int = int * exp(-1i * delta * t);
end