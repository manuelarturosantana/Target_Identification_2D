params = build_params( ...
    'k', 60, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

names = ["plane1","plane2","plane3","plane4","plane5"];
% names = ["plane1"];
for ii = 1:numel(names)
    gparams = build_twinjet_params('preset', names{ii});
    geom = twinjet_plane_geometry(gparams);
    C0 = build_curve(geom, params);
    C0 = refine_curve_wavenumber(C0, 2*pi/params.k, 1);
    
    figure(ii)
    clf
    tiledlayout(1,2)
    nexttile
    plot_curve(C0)
    title(names{ii} + " ppw = 1")
    
    C1 = bisect_all_patches(C0);
    C1 = build_global_indexing(C1);
    C1 = get_close_points(C1);
    
    nexttile
    plot_curve(C1)
end

