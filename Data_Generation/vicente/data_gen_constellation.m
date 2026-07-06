% data_gen_constellation.m
clear
clc

% cluster paths
addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

% laptop (windows) paths
% addpath(genpath("C:\Users\vicen\Code\rp2d_matlab\src"));
% addpath(genpath("C:\Users\vicen\Code\Target_Identification_2D"));

save_data_dir = "/scratch/vhojas/data_same_object_constellations_2d_fast";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

% ---- user-defined parameters ----

w_0 = 10;
freq_gauss_width = 5;
gauss_picker_tol = 1e-10;
numw = 100;
t_end = 400;
t_offset = 30;
num_receivers_per_obs = 100;
wimag = -0.5;

max_num_objects = 4;

same_direction_angle = 3*pi/2;
local_ellipse_opening_angle = 3*pi/2;
random_orientation_seed = 1;
center_spacing_factor = 2.6;

% Curve refinement, matching the faster plane script.
points_per_wavelength = 1;
num_extra_bisections = 1;

% ---------------------------------

rng(random_orientation_seed, "twister");

[sigmas, wlims] = gauss_picker( ...
    w_0, freq_gauss_width, gauss_picker_tol);

lambda = (2*pi) / wlims(2);
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

orientation_families = [
    "same_direction"
    "facing_away_origin"
    "random_orientation"
    "facing_origin"
];

configs = struct([]);
cfg_id = 0;

for iobj = 1:numel(obj_names)
    for nobj = 1:max_num_objects

        max_object_radius = max(obj_aa(iobj), obj_bb(iobj));

        centers = make_constellation_centers( ...
            nobj, max_object_radius, center_spacing_factor);

        ids = repmat(iobj, nobj, 1);

        if nobj == 1
            families_here = "single";
        else
            families_here = orientation_families;
        end

        for ifam = 1:numel(families_here)
            family = families_here(ifam);

            opening_angles = make_opening_angles( ...
                centers, family, same_direction_angle);

            if obj_names(iobj) == "elliptic"
                object_rotations = opening_angles_to_ellipse_rotations( ...
                    opening_angles, local_ellipse_opening_angle);
            else
                object_rotations = opening_angles;
            end

            cfg_id = cfg_id + 1;

            configs(cfg_id).name = ...
                sprintf("%s_N%d_%s", family, nobj, obj_tags(iobj));

            configs(cfg_id).description = make_config_description( ...
                obj_names(iobj), nobj, family);

            configs(cfg_id).orientation_family = family;
            configs(cfg_id).num_objects = nobj;
            configs(cfg_id).object_indices = ids;
            configs(cfg_id).object_names = obj_names(ids);

            configs(cfg_id).aa = obj_aa(ids);
            configs(cfg_id).bb = obj_bb(ids);
            configs(cfg_id).openings = obj_openings(ids);
            configs(cfg_id).centers = centers;

            configs(cfg_id).opening_angles = opening_angles;
            configs(cfg_id).object_rotations = object_rotations;

            configs(cfg_id).same_direction_angle = same_direction_angle;
            configs(cfg_id).local_ellipse_opening_angle = local_ellipse_opening_angle;
            configs(cfg_id).random_orientation_seed = random_orientation_seed;
            configs(cfg_id).center_spacing_factor = center_spacing_factor;
        end
    end
end

fprintf("Number of configurations: %d\n", numel(configs));

config_id = str2double(getenv("SLURM_ARRAY_TASK_ID"));

if isnan(config_id)
    config_id = str2double(getenv("CONFIG_ID"));
end

if isnan(config_id)
    error("No config id found. Set SLURM_ARRAY_TASK_ID or CONFIG_ID.");
end

if config_id < 1 || config_id > numel(configs)
    error("Config id %d is outside 1:%d.", config_id, numel(configs));
end

receiver_angles = linspace(0, 2*pi, num_receivers_per_obs + 1);
receiver_angles(end) = [];

xs_base = cos(receiver_angles);
ys_base = sin(receiver_angles);

start = tic;

for ii = config_id
    config = configs(ii);

    fprintf("\n==================================================\n");
    fprintf("Running configuration %d / %d\n", ii, numel(configs));
    fprintf("%s\n", config.name);
    fprintf("==================================================\n");

    geom = build_constellation_geometry_from_config(config);

    % Build a refined curve exactly in the spirit of the faster script:
    % wavelength refinement followed by one uniform bisection.
    C = build_curve(geom, params);
    C = refine_curve_wavenumber(C, 2*pi/params.k, points_per_wavelength);

    for jj = 1:num_extra_bisections
        C = bisect_all_patches(C);
        C = build_global_indexing(C);
        C = get_close_points(C);
    end

    lp = RP2LP(params, geom);

    lp.RP_Curve = C;
    [xb, yb] = get_global_xy(C);
    lp.curve.X = xb;
    lp.curve.Y = yb;
    lp.curve.Flag = "open";

    N = lp.RP_Curve.global.N;
    fprintf("N = %d\n", N);

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
        [uff, ~, r_res, ~, ~, f_sols] = cts_wpols_LU( ...
            psp, lp, pols_wp, L_rl, U_rl, L_pole, U_pole);
        single_receiver_comp_time = toc

        uff_all = [uff_all; uff(:).'];
        r_res_all = [r_res_all; r_res(:).'];

        clear f_sols
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
        "points_per_wavelength", ...
        "num_extra_bisections", ...
        "pol_comp_time", ...
        "LU_time_real_line", ...
        "LU_time_pols");
end

total_time = toc(start)

function centers = make_constellation_centers(nobj, max_object_radius, spacing_factor)

if nobj == 1
    centers = [0.0, 0.0];
    return;
end

min_center_distance = spacing_factor * max_object_radius;
center_radius = min_center_distance / (2 * sin(pi / nobj));

theta = linspace(0, 2*pi, nobj + 1).';
theta(end) = [];

centers = [
    center_radius * cos(theta), ...
    center_radius * sin(theta)
];

end

function opening_angles = make_opening_angles( ...
    centers, family, same_direction_angle)

nobj = size(centers, 1);
opening_angles = zeros(nobj, 1);

switch family
    case {"single", "same_direction"}
        opening_angles(:) = same_direction_angle;

    case "facing_away_origin"
        for j = 1:nobj
            c = centers(j, :);
            if norm(c) == 0
                opening_angles(j) = same_direction_angle;
            else
                opening_angles(j) = atan2(c(2), c(1));
            end
        end

    case "random_orientation"
        opening_angles = 2*pi*rand(nobj, 1);

    case "facing_origin"
        for j = 1:nobj
            c = centers(j, :);
            if norm(c) == 0
                opening_angles(j) = same_direction_angle;
            else
                opening_angles(j) = atan2(-c(2), -c(1));
            end
        end

    otherwise
        error("Unknown orientation family: %s", family);
end

opening_angles = mod(opening_angles, 2*pi);

end

function rotations = opening_angles_to_ellipse_rotations( ...
    opening_angles, local_opening_angle)

rotations = mod(opening_angles - local_opening_angle, 2*pi);

end

function geom = build_constellation_geometry_from_config(config)

if all(config.object_names == "elliptic")
    geom = elliptic_constellation_geometry( ...
        config.aa(1), ...
        config.bb(1), ...
        config.openings(1), ...
        config.centers, ...
        config.object_rotations);

elseif all(config.object_names == "circular")
    geom = cavity_constellation_geometry( ...
        config.aa, ...
        config.bb, ...
        config.openings, ...
        config.centers, ...
        config.opening_angles);

else
    error("Mixed object constellations are not supported by this script.");
end

end

function geom = elliptic_constellation_geometry(a, b, opening, centers, rotations)

geom = struct();
geom.is_closed = false;
geom.pieces = {};

for j = 1:size(centers, 1)
    gj = ellipticcavity(a, b, opening, rotations(j));

    for k = 1:numel(gj.pieces)
        geom.pieces{end + 1, 1} = translate_piece( ...
            gj.pieces{k}, centers(j, 1), centers(j, 2));
    end
end

end

function Q = translate_piece(P, cx, cy)

Q = P;

Q.x = @(t) P.x(t) + cx;
Q.y = @(t) P.y(t) + cy;

if isfield(P, "xj")
    Q.xj = P.xj + cx;
end

if isfield(P, "yj")
    Q.yj = P.yj + cy;
end

end

function description = make_config_description(object_name, nobj, family)

switch family
    case "single"
        family_text = "single canonical bottom-facing opening";

    case "same_direction"
        family_text = "all openings point in the same fixed direction";

    case "facing_away_origin"
        family_text = "all openings point radially away from the origin";

    case "random_orientation"
        family_text = "all openings have independent uniform random orientations";

    case "facing_origin"
        family_text = "all openings point toward the origin";

    otherwise
        family_text = "unknown orientation family";
end

description = sprintf( ...
    "%d same-object %s constellation, %s.", ...
    nobj, object_name, family_text);

end