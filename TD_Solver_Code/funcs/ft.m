% TODO: Doesn't work with non-symmetric intervals.
function val = ft(fs,a,b,t, is_inverse, is_coeffs)
    % Performs the inverse fourier transform using the sinc method.
    % Inputs:
    %   fs  : the function values F(omega), evaluated on an equally spaced
    %         omega grid. Warning! Must be an even number of these.
    %   a,b : The limits of the integral
    %   t   : The t value to evaluate the integral at
    %   is_inverse: If true use e^(-iwt) otherwise e^(iwt) Default true
    %   is_coeffs: If true assumes fs are the correct fourier coefficients. Default false
    % Outputs : 
    %   val : The inverse fourier transform at a point

    if nargin < 5
        is_inverse = true;
    end
    if nargin < 6
        is_coeffs = false;
    end

    A     = (b - a) / 2;
    delta = (b + a) / 2;
    M     = length(fs);
    ms    = -M/2:(M/2-1);
    % The fft assumes that the signal is periodic in [0,2*pi]. Thus we take
    % our sample periodic in [a,b], with a negative, and shift it to be
    % periodic in [0, a + b] by calling fftshift. 
    if is_coeffs
        cms = fs;
    else
        cms   = (1/M) * fft(ifftshift(fs));
        % This shifts the coefficients to the form in the paper.
        cms   = fftshift(cms);
    end

    P     = 2 * A; alpha = P / (2 *pi);
    if is_inverse
        % See doc sinc for why this works, even for sinc(0).
        val   = P .* sinc((alpha *t - ms));
        val   = exp(-1i * delta  * t) * sum(cms .* val);
    else
        val   = P .* sinc((-alpha *t - ms));
        val   = exp(1i * delta  * t) * sum(cms .* val);
    end
end



% Sanity Check code if it is needed.
% test_func = @(w) cos(w) .* exp(-(w-1).^2);
% % test_func = @(t) cos(t)
% % ws = linspace(-10,10,100);
% % plot(ws, test_func(ws))
% exact = @(t) exp(-t^2/4) / sqrt(2);
% % a=-pi; b = pi; t = 3;
% a=-4; b = 4; t = 4;
% M = 100; 
% % ms = (0:M-1);
% ms = (-M/2:M/2-1);
% dx = (b - a) / M;
% xn = ms * dx;
% fs = test_func(xn);
% 
% fs = fftshift(fs);
% plot(fs)
% A = (b - a) / 2;
% delta = (b + a) / 2;
% cms = (1/M) * fft(fs);
% cms = fftshift(cms);
% plot(real(cms))
% P = 2 * A;
% is_inverse = false;
% alpha = P / (2 *pi);
% val = P .* sin(((pi * A * 2)/P)* (alpha *t - ms)) ;
% val = val ./ (pi * (alpha * t - ms));
% val = sum(cms .* val);
% % val = sum(cms .* 2*A .* sinc((A/pi) * t - ms))
% val = ft(test_func(xn),a,b,t,is_inverse)
% true_val = integral(@(w) test_func(w) .* exp(1i * w * t),-4,4)
% % exact(t) * sqrt(2 * pi)