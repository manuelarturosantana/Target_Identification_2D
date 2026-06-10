% Script to compute the RSC of a plane
w_0 = 50.18;
params = build_params( ...
    'k', w_0, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 160, ...
    'n', 20, ...
    'p', 3, ...
    'p_edge', 2);

gparams = build_twinjet_params('preset', 'plane4');
geom = twinjet_plane_geometry(gparams);

C = build_curve(geom, params);
C = refine_curve_wavenumber(C, 2*pi/params.k, 1);

for jj = 1:1
    C = bisect_all_patches(C);
    C = build_global_indexing(C);
    C = get_close_points(C);
end

lp = RP2LP(params, geom);
lp = lp.update_RP_Curve(C);

%%
A = lp.bie_mat(w_0);
[L,U,P] = lu(A);
bndpts = [lp.curve.X, lp.curve.Y]';

n_ff = 500;
[x, y] = ff_points_equally_spaced(0, 2*pi, n_ff);
ff_vals = [];


for ii = 1:length(x)
    ii
    xhat = [x(ii), y(ii)];
    
    rhs = -exp(1i * w_0 * (-xhat * bndpts)).';
    den = U \ (L \ (P*rhs));
    
    ff_vals(ii) = eval_far_field(lp, w_0,den,xhat.');
end

%%

rcs = abs(ff_vals).^2;
angles = linspace(0, 2*pi, n_ff);
figure(1)
clf
polarplot(angles,rcs)

