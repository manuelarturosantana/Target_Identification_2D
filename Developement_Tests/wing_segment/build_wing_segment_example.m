clear; close all; clc;

% Geometry parameters.
gparams = build_twinjet_params();
gparams.separation_lead = 0;         % leading edge separation
gparams.separation_trail = 0;         % trail edge separation

% Build geometry.
geom = wing_section_geometry(gparams);

% Build curve parameters.
params = struct();
params.n = 16;
params.p_edge = 2;
params.delta = 0.05;

% Build curve.
Curve = build_curve(geom, params);
Curve = bisect_all_patches(Curve);

% Plot options.
opts = struct();
opts.nplot = 300;
opts.show_nodes = true;
opts.show_labels = true;
opts.node_size = 20;
opts.line_width = 1.5;
opts.same_color = false;
opts.plot_open_closed = true;
opts.plot_normals = true;
opts.normal_stride = 2;
opts.normal_scale = 0.15;

% Plot.
plot_curve(Curve, opts);