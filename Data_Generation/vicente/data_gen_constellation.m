% generate_two_cavity_constellation_data.m
clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/chebfun"))

save_data_dir = "/scratch/vhojas/data_two_cavity_constellations_2d";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

w_0 = 10;
freq_gauss_width = 5;
gauss_picker_tol = 1e-10;
numw = 100;
t_end = 400;
t_offset = 30;
num_receivers_per_obs = 100;
wimag = -0.5;

[sigmas, wlims] = gauss_picker( ...
    w_0, freq_gauss_width, gauss_picker_tol);

lambda = (2 * pi) / wlims(2);
tlims = [0, t_end];
numt = ceil((tlims(2) / lambda) * 15);

params = build_params( ...
    'k', w_0 + freq_gauss_width, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

obj_names = [
    "circular"
    "elliptic"
];

obj_tags = [
    "circ_o050"
    "ell_o075"
];

obj_aa = [1.0; 1.6];
obj_bb = [1.0; 0.8];
obj_openings = [0.50; 0.75];

right_opening = 0;
left_opening = pi;
bottom_opening = 3*pi/2;

single_center = [0.0, 0.0];

pair_centers_aligned = [
    -2.25, 0.0
     2.25, 0.0
];

pair_centers_offset = [
    -2.25,  0.5
     2.25, -0.5
];

configs = struct([]);
cfg_id = 0;

for iobj = 1:2
    cfg_id = cfg_id + 1;

    configs(cfg_id).name = "single_" + obj_tags(iobj);
    configs(cfg_id).description = ...
        "Single object at origin, bottom-facing.";
    configs(cfg_id).placement = "single_center";
    configs(cfg_id).object_indices = iobj;
    configs(cfg_id).object_names = obj_names(iobj);

    configs(cfg_id).aa = obj_aa(iobj);
    configs(cfg_id).bb = obj_bb(iobj);
    configs(cfg_id).openings = obj_openings(iobj);
    configs(cfg_id).centers = single_center;
    configs(cfg_id).opening_angles = bottom_opening;
end

case_id = 0;
cases = struct([]);

same_obj_pairs = [
    1, 1
    2, 2
];

diff_obj_pairs = [
    1, 2
    2, 1
];

for ipair = 1:size(same_obj_pairs, 1)
    case_id = case_id + 1;
    ids = same_obj_pairs(ipair, :).';

    cases(case_id).family = "same_dir_same_obj";
    cases(case_id).description = ...
        "Two equal objects, both bottom-facing.";
    cases(case_id).object_indices = ids;
    cases(case_id).opening_angles = [
        bottom_opening
        bottom_opening
    ];
end

for ipair = 1:size(diff_obj_pairs, 1)
    case_id = case_id + 1;
    ids = diff_obj_pairs(ipair, :).';

    cases(case_id).family = "same_dir_diff_obj";
    cases(case_id).description = ...
        "Two different objects, both bottom-facing.";
    cases(case_id).object_indices = ids;
    cases(case_id).opening_angles = [
        bottom_opening
        bottom_opening
    ];
end

for ipair = 1:size(same_obj_pairs, 1)
    case_id = case_id + 1;
    ids = same_obj_pairs(ipair, :).';

    cases(case_id).family = "facing_same_obj";
    cases(case_id).description = ...
        "Two equal objects facing each other.";
    cases(case_id).object_indices = ids;
    cases(case_id).opening_angles = [
        right_opening
        left_opening
    ];
end

for ipair = 1:size(diff_obj_pairs, 1)
    case_id = case_id + 1;
    ids = diff_obj_pairs(ipair, :).';

    cases(case_id).family = "facing_diff_obj";
    cases(case_id).description = ...
        "Two different objects facing each other.";
    cases(case_id).object_indices = ids;
    cases(case_id).opening_angles = [
        right_opening
        left_opening
    ];
end

for icase = 1:numel(cases)
    ids = cases(icase).object_indices;

    for placement_id = 1:2
        if placement_id == 1
            placement = "aligned";
            centers = pair_centers_aligned;
        else
            placement = "offset";
            centers = pair_centers_offset;
        end

        cfg_id = cfg_id + 1;

        configs(cfg_id).name = ...
            cases(icase).family + "_" + placement + "_" ...
            + obj_tags(ids(1)) + "_" + obj_tags(ids(2));

        configs(cfg_id).description = cases(icase).description;
        configs(cfg_id).placement = placement;
        configs(cfg_id).object_indices = ids;
        configs(cfg_id).object_names = obj_names(ids);

        configs(cfg_id).aa = obj_aa(ids);
        configs(cfg_id).bb = obj_bb(ids);
        configs(cfg_id).openings = obj_openings(ids);
        configs(cfg_id).centers = centers;
        configs(cfg_id).opening_angles = cases(icase).opening_angles;
    end
end

fprintf("Number of configurations: %d\n", numel(configs));

% Self-convergence check on one representative configuration.
tol_selfconv = 1e-7;
max_manual_refs = 5;
selfconv_config_id = numel(configs);

[num_manual_refs, selfconv_info] = choose_manual_refinements( ...
    params, configs(selfconv_config_id), tol_selfconv, max_manual_refs);

fprintf("\nUsing %d manual refinements for data generation.\n", ...
    num_manual_refs);

receiver_angles = linspace(0, 2*pi, num_receivers_per_obs + 1);
receiver_angles(end) = [];

xs_base = cos(receiver_angles);
ys_base = sin(receiver_angles);

start = tic;

for ii = 1:numel(configs)
    config = configs(ii);

    fprintf("\n==================================================\n");
    fprintf("Running configuration %d / %d\n", ii, numel(configs));
    fprintf("%s\n", config.name);
    fprintf("==================================================\n");

    geom = cavity_constellation_geometry( ...
        config.aa, ...
        config.bb, ...
        config.openings, ...
        config.centers, ...
        config.opening_angles);

    lp = RP2LP(params, geom);
    lp = refine_lp_manually(lp, num_manual_refs);

    N = lp.RP_Curve.global.N;

    w_len = wlims(2) - wlims(1);
    num_sub_ints = ceil(w_len / 2.0);
    jump = w_len / num_sub_ints;

    pols_wp = [];
    pstart = tic;

    for jj = 1:num_sub_ints
        wlims_loc = [
            wlims(1) + (jj - 1) * jump, ...
            wlims(1) + jj * jump
        ];

        psp = problem_data( ...
            'wlims', wlims_loc, ...
            'w_0', mean(wlims_loc), ...
            'numw', numw, ...
            'xs', 0, ...
            'ys', 0, ...
            'wimag', wimag, ...
            'N', N);

        pols_wp_loc = comp_poles(psp, lp);
        pols_wp_loc = sort(pols_wp_loc, "comparisonMethod", "real");
        pols_wp = [pols_wp; pols_wp_loc(:)];
    end

    pol_comp_time = toc(pstart)

    psp = problem_data( ...
        'wlims', wlims, ...
        'w_0', w_0, ...
        'numw', numw, ...
        'xs', 0, ...
        'ys', 0, ...
        'wimag', wimag, ...
        'N', N);

    tic
    [L_rl, U_rl] = comp_bie_LU(psp.ws, lp);
    LU_time_real_line = toc

    tic
    [L_pole, U_pole] = comp_pole_LU(psp, lp, pols_wp);
    LU_time_pols = toc

    xs = xs_base;
    ys = ys_base;

    uff_all = [];
    r_res_all = [];

    for rind = 1:num_receivers_per_obs
        x = xs(rind);
        y = ys(rind);

        psp = problem_data( ...
            'wlims', wlims, ...
            'N', N, ...
            'Tlims', tlims, ...
            'kappa', [-x, -y], ...
            'numt', numt, ...
            'numw', numw, ...
            'w_0', w_0, ...
            't_0', t_offset, ...
            'sigmas', sigmas, ...
            'wimag', wimag, ...
            'is_far_field', true, ...
            'xs', x, ...
            'ys', y);

        tic
        [uff, ~, r_res] = cts_wpols_LU( ...
            psp, lp, pols_wp, L_rl, U_rl, L_pole, U_pole);
        single_receiver_comp_time = toc

        uff_all = [uff_all; uff(:).'];
        r_res_all = [r_res_all; r_res(:).'];
    end

    ts = psp.ts;
    curveX = lp.curve.X;
    curveY = lp.curve.Y;
    pols = pols_wp;

    config_name = config.name;
    config_info = config;
    solver_params = params;

    fname = sprintf("%03d_%s.mat", ii, config.name);

    save(fullfile(save_data_dir, fname), ...
        "r_res_all", ...
        "uff_all", ...
        "pols", ...
        "ts", ...
        "xs", ...
        "ys", ...
        "receiver_angles", ...
        "curveX", ...
        "curveY", ...
        "config_name", ...
        "config_info", ...
        "solver_params", ...
        "selfconv_info", ...
        "num_manual_refs", ...
        "pol_comp_time", ...
        "LU_time_real_line", ...
        "LU_time_pols");
end

total_time = toc(start)

function [num_manual_refs, info] = choose_manual_refinements( ...
    params, config, tol, max_refs)
% Choose manual refinements using far-field convergence of the solved BIE.

fprintf("\n==================================================\n");
fprintf("Self-convergence calibration\n");
fprintf("Config: %s\n", config.name);
fprintf("==================================================\n");

geom = cavity_constellation_geometry( ...
    config.aa, ...
    config.bb, ...
    config.openings, ...
    config.centers, ...
    config.opening_angles);

xhat_inc = [1/sqrt(2); 1/sqrt(2)];

nobs = 64;
theta = linspace(0, 2*pi, nobs + 1);
theta(end) = [];
xhat_obs = [
    cos(theta)
    sin(theta)
];

Lps = cell(max_refs + 1, 1);
Phis = cell(max_refs + 1, 1);
Ns = zeros(max_refs + 1, 1);

lp = RP2LP(params, geom);

for lev = 1:(max_refs + 1)
    phi = solve_planewave_scattering( ...
        lp.RP_Curve, xhat_inc, params);

    Lps{lev} = lp;
    Phis{lev} = phi;
    Ns(lev) = lp.RP_Curve.global.N;

    fprintf("  Level %d: N = %d\n", lev - 1, Ns(lev));

    if lev <= max_refs
        lp = refine_lp_once(lp);
    end
end

u_ref = Lps{end}.eval_far_field(params.k, Phis{end}, xhat_obs);

errs = zeros(max_refs, 1);

fprintf("\nlev   N        max far-field error\n");
fprintf("-----------------------------------\n");

for lev = 1:max_refs
    u = Lps{lev}.eval_far_field(params.k, Phis{lev}, xhat_obs);
    errs(lev) = max(abs(u(:) - u_ref(:)));

    fprintf("%3d  %6d   %.3e\n", lev - 1, Ns(lev), errs(lev));
end

idx = find(errs <= tol, 1, "first");

if isempty(idx)
    error("Self-convergence failed: no level reached %.3e.", tol);
end

num_manual_refs = idx - 1;

info = struct();
info.config_name = config.name;
info.tol = tol;
info.max_refs = max_refs;
info.num_manual_refs = num_manual_refs;
info.Ns = Ns;
info.errs = errs;

fprintf("-----------------------------------\n");
fprintf("Selected manual refinements: %d\n", num_manual_refs);
fprintf("Selected error: %.3e\n", errs(idx));

end

function lp = refine_lp_manually(lp, nrefs)
% Refine an RP2LP object without changing the RP2LP class.

for j = 1:nrefs
    lp = refine_lp_once(lp);
end

end

function lp = refine_lp_once(lp)
% One manual bisection refinement of lp.RP_Curve and lp.curve.

C = bisect_all_patches(lp.RP_Curve);

lp.RP_Curve = C;

[xb, yb] = get_global_xy(C);
lp.curve.X = xb;
lp.curve.Y = yb;
lp.curve.Flag = "vicente open";

end