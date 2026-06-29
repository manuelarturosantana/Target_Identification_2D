% test_plane_scatter_self_convergence_field.m
% Could it be that the true poles are lower in the complex plane?
%
clear; clc;
addpath(genpath('/home/msantana/Progamming/rp2d_matlab/src'))
addpath(genpath('/home/msantana/Progamming/Target_Identification_2D'))
addpath(genpath('/home/msantana/Progamming/MATLAB_PACKAGES'))

rng(1);
fprintf("==============================================\n");
fprintf(" Wing Segment Convergence test (poles in user range)\n");
fprintf("==============================================\n\n");

% --------------------------------------------------
% plane parameters
% --------------------------------------------------
params = build_params();
params.n = 20;
params.p_edge = 2;
params.p = 3; 
params.delta = 0.05;
params.k = 60;
params.npolar = 160;
params.solver = 'direct';

% --------------------------------------------------
% User-specified frequency range for pole finding
% --------------------------------------------------
% This is 18/20 for the open circle.
wlims = [48,60];   % <-- edit these bounds
wimag = -0.5;
numw  = 400;

% --------------------------------------------------
% refinement levels
% --------------------------------------------------
NREF = 2;   % comparisons; finest is reference

% --------------------------------------------------
% base geometry
% --------------------------------------------------
% gparams = build_twinjet_params();
% geom = twinjet_plane_geometry(gparams);

gparams = build_twinjet_params();
geom = wing_section_geometry(gparams);

% geom = circularcavity(1, pi/8);
% geom = box_geometry(1, 1);
%%
poles  = {};
all_svd_vals = {};

C0 = build_curve(geom, params);
C0 = refine_curve_wavenumber(C0, 2*pi/params.k, 1);

Curves = cell(NREF + 1, 1);

Curves{1} = C0;

lp = RP2LP(params, geom, 1);

fprintf("Building nested curves + finding poles...\n");
for lev = 1:(NREF+1)

    % build LP object (needed for comp_poles)
    lp = lp.update_RP_Curve(Curves{lev});
    N = lp.RP_Curve.global.N;

    %subdivide frequency range into sub-intervals (same rule as gen script)
    w_len       = wlims(2) - wlims(1);
    num_sub_ints = ceil(w_len / 2.0);
    jump        = w_len / num_sub_ints;

    pols_wp = [];
    for jj = 1:num_sub_ints
        wlims_loc = [wlims(1) + (jj-1)*jump, ...
                     wlims(1) +  jj   *jump];

        psp = problem_data( ...
            'wlims',  wlims_loc, ...
            'w_0',    (wlims_loc(2) + wlims_loc(1)) / 2.0, ...
            'numw',   numw, ...
            'xs',     0, ...
            'ys',     0, ...
            'wimag',  wimag, ...
            'N',      N);
        % pols_wp_loc = 1;
        pols_wp_loc = comp_poles(psp, lp);
        pols_wp_loc = sort(pols_wp_loc, 'comparisonMethod', 'real');
        pols_wp = [pols_wp; pols_wp_loc(:)];
    end
    % u = rand(1,N,'like',1i); v = rand(N,1,'like',1i);
    % f = @(z) u * (lp.bie_mat(z) \ v);
    % 
    % [pols_wp, ~] = aaa_recursive1d(f, wlims(1), wlims(2), 200);

    svd_vals = [];
    for ii = 1:length(pols_wp)
        S = svd(lp.bie_mat(pols_wp(ii)));
        svd_vals(end+1) = S(end);
    end

    poles{end+1} = pols_wp;
    all_svd_vals{end+1} = svd_vals;

    fprintf("  Level %d: patches = %4d, N = %6d, poles found = %d\n", ...
        lev-1, numel(lp.RP_Curve.patches), N, numel(pols_wp));

    % refine for next level (ensure bookkeeping after bisection)
    if lev <= NREF
        Cn = bisect_all_patches(Curves{lev});
        Cn = build_global_indexing(Cn);
        Cn = get_close_points(Cn);
        Curves{lev + 1} = Cn;
    end
end

% --------------------------------------------------
% Report pole convergence across levels
% --------------------------------------------------
save('poles_convergence.mat', 'poles', 'all_svd_vals', 'wlims', 'wimag', 'numw', 'params');
fprintf("\nPoles saved to poles_convergence.mat\n");


%%
% gparams = build_twinjet_params();
% geom = wing_section_geometry(gparams);

C0 = build_curve(geom, params);
C0 = refine_curve_wavenumber(C0, 2*pi/params.k, 1);

figure(1)
clf
plot_curve(C0)

C1 = bisect_all_patches(C1);
C1 = build_global_indexing(C1);
C1 = get_close_points(C1);

figure(2)
clf
plot_curve(C1)

