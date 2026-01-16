function [xk, wf1] =fejer(n)
% Weights of Fejer's first quadature rule
% x_k are the sample points
% wf1 are the weights.
% Taken from Fast Construction ofthe Fejerand Clenshaw-Curtis Quadrature Rules
% Jorg Waldvogel 2006
ns = (1:n) - (1/2);
xk = cos(ns*pi/n);
xk = xk(:);
N= [1:2:n-1]'; l=length(N); m=n-l; K= [0:m-1]';
v0=[2*exp(1i*pi*K/n)./(1-4*K.^2);zeros(l+1,1)];
v1=v0(1:end-1)+conj(v0(end:-1:2));wf1=ifft(v1);
end
