% Sanity Check code if it is needed.
test_func = @(w) exp(-(w-10).^2);
ws = linspace(-20,20,100);
figure(1)
plot(ws, test_func(ws))
% test_func = @(w) exp(-(w).^2);
a=4; b = 16; t = -1;
% We need to incorporate delta into this to sample the points at the
% correct spots
delta = (b + a) / 2;
M = 50; 
ms = (-M/2:M/2-1);
dx = (b - a) / M;
xn = ms * dx;
fs = test_func(xn + delta);

fs = fftshift(fs);
% plot(fs)
A = (b - a) / 2;
delta = (b + a) / 2;
cms = (1/M) * fft(fs);
cms = fftshift(cms);


clc
val = ft(cms,a,b,t,true,true)
% val = ft(test_func(xn),a,b,t,true)
true_val = integral(@(w) test_func(w) .* exp(-1i * w * t),a,b)

%%
% test_func = @(w) exp(-(w).^2);
% ws = linspace(-10,10,100);
% figure(1)
% plot(ws, test_func(ws))
% exact = @(t) exp(-t^2/4) / sqrt(2);
% % a=-pi; b = pi; t = 3;
% a=-4; b = 4; t = 0;
% M = 201; 
% ms = (-M/2:M/2-1);
% dx = (b - a) / M;
% xn = ms * dx;
% fs = test_func(xn);
% 
% fs = fftshift(fs);
% % plot(fs)
% A = (b - a) / 2;
% delta = (b + a) / 2;
% cms = (1/M) * fft(fs);
% cms = fftshift(cms);
% % plot(real(cms))
% P = 2 * A;
% is_forward = false;
% alpha = P / (2 *pi);
% val = P .* sin(((pi * A * 2)/P)* (alpha *t - ms)) ;
% val = val ./ (pi * (alpha * t - ms));
% val = sum(cms .* val)
% % val = sum(cms .* 2*A .* sinc((A/pi) * t - ms))
% val = ft(test_func(xn),a,b,t,is_forward)
% true_val = integral(@(w) test_func(w) .* exp(-1i * w * t),a,b,'AbsTol',1e-10)
% % exact(t) * sqrt(2 * pi)


% a_t = @(t) 5 * exp(-(t - 6).^2 / 2);
% a=-10; b = 10;
% M = 41; 
% % ms = (0:M-1);
% ms = (-M/2:M/2-1);
% dx = (b - a) / M;
% xn = ms * dx;
% t_vals = xn;
% 
% 
% ak = a_t(t_vals + s_k(ss)) .* window(t_vals + s_k(ss), H);
% figure(2)
% plot(linspace(-H,H,length(ak)),ak)
% 
% Bkslow(ss, ww) = ft(ak,-H,H,0,false);
% % TESTING: check the integrals are computed correctly.
% % TODO: There is a bug here it seems!
% Bkslow(ss,ww)
% akt = @(t) a_t(t + s_k(ss)) .* window(t + s_k(ss),H);
% integral(@(t) akt(t) .* exp(1i * 0 *t ),-H,H)





















