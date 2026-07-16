% data_gen_deformation_experiment.m
clear
clc

% cluster paths
addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

% laptop (windows) paths
% addpath(genpath("C:\Users\vicen\Code\rp2d_matlab\src"));
% addpath(genpath("C:\Users\vicen\Code\Target_Identification_2D"));

save_data_dir = "/scratch/vhojas/data_deformation_experiment_2d";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

% ---- scattering parameters ----

w_0 = 10;
freq_gauss_width = 5;
gauss_picker_tol = 1e-10;
numw = 100;
t_end = 400;
t_offset = 30;
num_receivers_per_obs = 100;
wimag = -0.5;

points_per_wavelength = 1;
num_extra_bisections = 1;

% ---- deformation experiment ----

object_names = ["circular"; "elliptic"];
object_tags = ["circ_o035"; "ell_o045"];

object_a = [1.0; 1.5];
object_b = [1.0; 0.5];
object_opening = [0.35; 0.45];
object_rotation = [0.0; 0.0];

amplitudes = [2e-2; 1e-2; 5e-3];
nmodes = [2; 4; 6];

num_realizations = 1;
random_seed = 1;

% ---------------------------------

[sigmas, wlims] = gauss_picker( ...
    w_0, freq_gauss_width, gauss_picker_tol);

lambda = (2*pi)/wlims(2);
tlims = [0, t_end];
numt = ceil((tlims(2)/lambda)*15);

params = build_params( ...
    'k', w_0 + freq_gauss_width, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

% Generate one coefficient bank. The same random directions are used for
% both object classes and for every amplitude.
rng(random_seed, "twister");

cos_bank = cell(numel(nmodes), 1);
sin_bank = cell(numel(nmodes), 1);

for imode = 1:numel(nmodes)
    n = nmodes(imode);
    cos_bank{imode} = randn(n, num_realizations);
    sin_bank{imode} = randn(n, num_realizations);
end

configs = build_configs( ...
    object_names, object_tags, object_a, object_b, ...
    object_opening, object_rotation, amplitudes, nmodes, ...
    cos_bank, sin_bank, num_realizations, random_seed);

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
config = configs(config_id);

fprintf("\n==================================================\n");
fprintf("Running configuration %d / %d\n", config_id, numel(configs));
fprintf("%s\n", config.name);
fprintf("==================================================\n");

geom = ellipticcavity_deformed( ...
    config.a, config.b, config.opening, config.rotation, ...
    config.deformation);

% Build and refine the boundary discretization.
C = build_curve(geom, params);
C = refine_curve_wavenumber( ...
    C, 2*pi/params.k, points_per_wavelength);

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

% Compute the poles of the current physical geometry. These are required
% for generating its time-domain response. Classification can later use
% only the two nominal pole sets as its reference dictionary.
w_len = wlims(2) - wlims(1);
num_sub_ints = ceil(w_len/2.0);
jump = w_len/num_sub_ints;

pols_wp = [];
pstart = tic;

for jj = 1:num_sub_ints
    wlims_loc = [ ...
        wlims(1) + (jj - 1)*jump, ...
        wlims(1) + jj*jump ...
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
    pols_wp_loc = sort( ...
        pols_wp_loc, "comparisonMethod", "real");
    pols_wp = [pols_wp; pols_wp_loc(:)];
end

pol_comp_time = toc(pstart);

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
LU_time_real_line = toc;

tic
[L_pole, U_pole] = comp_pole_LU(psp, lp, pols_wp);
LU_time_pols = toc;

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

    [uff, ~, r_res, ~, ~, f_sols] = cts_wpols_LU( ...
        psp, lp, pols_wp, L_rl, U_rl, L_pole, U_pole);

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

fname = sprintf("%03d_%s.mat", config_id, config.name);

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

total_time = toc(start)

function configs = build_configs( ...
    object_names, object_tags, object_a, object_b, ...
    object_opening, object_rotation, amplitudes, nmodes, ...
    cos_bank, sin_bank, num_realizations, random_seed)

config_cells = {};
max_modes = max(nmodes);

for iobject = 1:numel(object_names)
    % One nominal configuration per class.
    deformation = struct();
    deformation.modes = (1:max_modes).';
    deformation.cos_coeffs = zeros(max_modes, 1);
    deformation.sin_coeffs = zeros(max_modes, 1);

    config_cells{end + 1} = make_config( ...
        object_names(iobject), object_tags(iobject), ...
        object_a(iobject), object_b(iobject), ...
        object_opening(iobject), object_rotation(iobject), ...
        0, 0, 0, deformation, random_seed); %#ok<AGROW>

    for imode = 1:numel(nmodes)
        n = nmodes(imode);

        for ireal = 1:num_realizations
            base_cos = cos_bank{imode}(:, ireal);
            base_sin = sin_bank{imode}(:, ireal);

            for iamp = 1:numel(amplitudes)
                amp = amplitudes(iamp);

                deformation = struct();
                deformation.modes = (1:n).';
                deformation.cos_coeffs = amp*base_cos;
                deformation.sin_coeffs = amp*base_sin;

                config_cells{end + 1} = make_config( ...
                    object_names(iobject), object_tags(iobject), ...
                    object_a(iobject), object_b(iobject), ...
                    object_opening(iobject), object_rotation(iobject), ...
                    ireal, n, amp, deformation, ...
                    random_seed); %#ok<AGROW>
            end
        end
    end
end

configs = [config_cells{:}];

end

function config = make_config( ...
    object_name, object_tag, a, b, opening, rotation, ...
    realization_id, num_modes, amplitude, deformation, random_seed)

if amplitude == 0
    config_name = sprintf("%s_nominal", object_tag);
    description = sprintf("Nominal %s cavity.", object_name);
else
    amplitude_tag = strrep(sprintf("%.4f", amplitude), ".", "p");
    config_name = sprintf( ...
        "%s_m%02d_r%03d_e%s", ...
        object_tag, num_modes, realization_id, amplitude_tag);
    description = sprintf( ...
        "%s cavity with %d Fourier modes, realization %d, amplitude %.4g.", ...
        object_name, num_modes, realization_id, amplitude);
end

config = struct();
config.name = config_name;
config.description = description;
config.object_name = object_name;
config.object_tag = object_tag;
config.a = a;
config.b = b;
config.opening = opening;
config.rotation = rotation;
config.realization_id = realization_id;
config.num_modes = num_modes;
config.deformation_amplitude = amplitude;
config.is_nominal = (amplitude == 0);
config.deformation = deformation;
config.random_seed = random_seed;

end
