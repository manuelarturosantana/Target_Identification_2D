% generate_two_cavity_constellation_data.m
%
% Generate data for one-object and two-object cavity configurations.
%
% Geometry is built only through:
%
%   geom = cavity_constellation_geometry( ...
%       aa, bb, openings, centers, opening_angles);

clear
clc

% Paths.
addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/chebfun"))

% Output directory.
save_data_dir = "/scratch/vhojas/data_two_cavity_constellations_2d";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

% Experiment params.
w_0 = 10;
freq_gauss_width = 5;
gauss_picker_tol = 1e-10;
numw = 100;
t_end = 400;
t_offset = 30;
num_receivers_per_obs = 100;
wimag = -0.5;

% Time-domain and AAA derived params.
[sigmas, wlims] = gauss_picker( ...
    w_0, freq_gauss_width, gauss_picker_tol);

lambda = (2 * pi) / wlims(2);
tlims = [0, t_end];
numt = ceil((tlims(2) / lambda) * 15);

% Solver params.
params = build_params( ...
    'k', w_0 + freq_gauss_width, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

% Object library.
%
% Object 1: circular cavity.
% Object 2: elliptic cavity.
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

% Opening-angle convention for the geometry builder.
right_opening = 0;
left_opening = pi;
bottom_opening = 3*pi/2;

% Single-object center.
single_center = [0.0, 0.0];

% Pair centers.
%
% The aligned case has both objects at the same y-level.
% The offset case shifts the objects by about half an object length in y.
pair_centers_aligned = [
    -2.25, 0.0
     2.25, 0.0
];

pair_centers_offset = [
    -2.25,  0.5
     2.25, -0.5
];

% Build configurations.
configs = struct([]);
cfg_id = 0;

% Single-object configurations.
for iobj = 1:2
    cfg_id = cfg_id + 1;

    configs(cfg_id).name = "single_" + obj_tags(iobj);
    configs(cfg_id).description = ...
        "Single object at the origin, bottom-facing opening.";

    configs(cfg_id).placement = "single_center";
    configs(cfg_id).object_indices = iobj;
    configs(cfg_id).object_names = obj_names(iobj);

    configs(cfg_id).aa = obj_aa(iobj);
    configs(cfg_id).bb = obj_bb(iobj);
    configs(cfg_id).openings = obj_openings(iobj);
    configs(cfg_id).centers = single_center;
    configs(cfg_id).opening_angles = bottom_opening;
end

% Base two-object logical cases.
case_id = 0;
cases = struct([]);

% Same direction, same object.
same_obj_pairs = [
    1, 1
    2, 2
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

% Same direction, different objects.
diff_obj_pairs = [
    1, 2
    2, 1
];

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

% Facing each other, same object.
for ipair = 1:size(same_obj_pairs, 1)
    case_id = case_id + 1;
    ids = same_obj_pairs(ipair, :).';

    cases(case_id).family = "facing_same_obj";
    cases(case_id).description = ...
        "Two equal objects facing each other horizontally.";
    cases(case_id).object_indices = ids;
    cases(case_id).opening_angles = [
        right_opening
        left_opening
    ];
end

% Facing each other, different objects.
for ipair = 1:size(diff_obj_pairs, 1)
    case_id = case_id + 1;
    ids = diff_obj_pairs(ipair, :).';

    cases(case_id).family = "facing_diff_obj";
    cases(case_id).description = ...
        "Two different objects facing each other horizontally.";
    cases(case_id).object_indices = ids;
    cases(case_id).opening_angles = [
        right_opening
        left_opening
    ];
end

% For each two-object logical case, create aligned and offset versions.
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

% Uniform receivers on the unit circle.
receiver_angles = linspace(0, 2*pi, num_receivers_per_obs + 1);
receiver_angles(end) = [];

xs_base = cos(receiver_angles);
ys_base = sin(receiver_angles);

% Run experiment.
start = tic;

for ii = 1:numel(configs)
    config = configs(ii);

    fprintf("\n==================================================\n");
    fprintf("Running configuration %d / %d\n", ii, numel(configs));
    fprintf("%s\n", config.name);
    fprintf("==================================================\n");

    aa = config.aa;
    bb = config.bb;
    openings = config.openings;
    centers = config.centers;
    opening_angles = config.opening_angles;

    geom = cavity_constellation_geometry( ...
        aa, bb, openings, centers, opening_angles);

    lp = RP2LP(params, geom);

    % Total number of points in the curve.
    N = lp.RP_Curve.global.N;

    % Rule of thumb for computing complex resonances.
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

        % Dummy problem data for pole computation.
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

    % Dummy problem for LU decompositions.
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
        "pol_comp_time", ...
        "LU_time_real_line", ...
        "LU_time_pols");
end

total_time = toc(start)
