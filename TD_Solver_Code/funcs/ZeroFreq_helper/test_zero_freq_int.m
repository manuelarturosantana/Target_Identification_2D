w_c = 1;          % Integrate between 0 and w_c
ts = 1000;
% ts  = linspace(32,40,1); % T values for ift integration
f = @(z) besselh(0,z);  % Input function        
N   = 8;         % Degree of Tchebyshev interpolant, (so uses N+1 points)
q   = N+1.1;           % Parameter in grated mesh. Must be > N + 1
npatch  = 20;     % Number of patches in grated mesh, \mathcal{M} in the paper
M   = 80;         % Parameter >= N for asymptotic expansion of the higher weights
is_inverse =true; % Reminder to take sign on Fourier transform correctly...

%% Number of patches grows.
% This follows Table 1 in the paper
npatches = [8,16,32,64]; 
g = @(z) f(z) .* exp(-1i * ts * z);
true_sol = integral(g,0,w_c,"AbsTol",1e-13,"RelTol",1e-13);
err = [];
for npind = npatches
    npat  = npind; 
    % q   = N + 1.1; 
    [intws,ccps,cs,hs] = grated_mesh_grid(w_c,npat,N,q);
    weights = comp_fcc_weights(ts,hs,npat,N,M, is_inverse);
    fvals = f(intws);
    sol = gfcc_int(fvals,ts,ccps,hs,cs,weights,npat, is_inverse);
    err = [err,abs(sol - true_sol)];
end
err





%% Number of points per patch grows
N   = 40; 
q   = N + 1.1; 
[intws,ccps,cs,hs] = grated_mesh_grid(w_c,npatch,N,q);
weights = comp_fcc_weights(ts,hs,npatch,N,M, is_inverse);
fvals = f(intws);
conv_sol = gfcc_int(fvals,ts,ccps,hs,cs,weights,npatch, is_inverse);
%% point convergence
ns = 2:1:38;
err = [];
for npind = ns
    N   = npind; 
    % q   = N + 1.1; 
    [intws,ccps,cs,hs] = grated_mesh_grid(w_c,npatch,N,q);
    weights = comp_fcc_weights(ts,hs,npatch,N,M, is_inverse);
    fvals = f(intws);
    sol = gfcc_int(fvals,ts,ccps,hs,cs,weights,npatch, is_inverse);
    err = [err,abs(sol - conv_sol)];
end
err
% To get convergence fix the number of patches, then increase N. Also fix
% q!













