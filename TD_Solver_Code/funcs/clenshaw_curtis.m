function [xk, wcc] = clenshaw_curtis(n)
% clenshaw_curtis weights and nodes in [-1,1] (closed closed)
% Taken from Fast Construction ofthe Fejerand Clenshaw-Curtis Quadrature Rules
% Jorg Waldvogel 2006

ks = 0:n;
xk = cos(ks * pi / n);

N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
v0=[2./N./(N-2); 1/N(end); zeros(m,1)]; v2=-v0(1:end-1)-v0(end:-1:2);
g0=-ones(n,1); g0(1+l)=g0(1+l)+n;g0(1+m)=g0(1+m) +n; 
g=g0/(n^2-1+mod(n,2));wcc=ifft(v2+g);

% wcc_n = w_cc = 0 = (1 / (n^2 - 1 + mod(n,2)) See the paper
wcc = [wcc;wcc(1)];
end