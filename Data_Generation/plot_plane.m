% paths
addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/chebfun"))

% where to save
save_data_dir = "/home/vhojas/Code/Target_Identification_2D";

% parameters
params = build_params( ...
    'k', 55, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

% build geometry and curve
geom_params = build_quadjet_params();
geom = quadjet_plane_geometry(geom_params);
Curve = build_curve(geom, params);
Curve = refine_curve_wavenumber(Curve, 2*pi/params.k, 1);

fig = figure;
plot_curve(Curve);
saveas(fig, 'quadjet_plot.png');
