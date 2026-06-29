% generates data for 5 twinjet and 5 quadjet geometries
% Now we use a refined curve to see if it makes the data better.

% paths
addpath(genpath('/home/msantana/Progamming/rp2d_matlab/src'))
addpath(genpath('/home/msantana/Progamming/Target_Identification_2D'))
addpath(genpath('/home/msantana/Progamming/MATLAB_PACKAGES'))

% ------------------------------------------------------------
% Output directory
% ------------------------------------------------------------
save_data_dir = "./data_gen_planes_check";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

% ------------------------------------------------------------
% Experiment params
% ------------------------------------------------------------
w_0              = 50;            % center frequency of the Gaussian
freq_gauss_width = 5;             % Gaussian frequency width
gauss_picker_tol = 1e-10;         % Gaussian truncation tolerance
numw             = 400;           % number of points for frequency integrals
t_end            = 400;           % final time
t_offset         = 30;            % time offset
num_receivers_per_obs = 100;
are_recivers_random = true;

% receiver angles
angle1 = -90*pi/180;
angle2 =  90*pi/180;

wimag  = -0.3;

% time-domain and AAA derived params
[sigmas, wlims] = gauss_picker(w_0, freq_gauss_width, ...
    gauss_picker_tol);
ws = linspace(wlims(1), wlims(2), 1000);
psp = problem_data( ...
    'wlims', wlims, ...
    'sigmas', sigmas, ...
    'w_0', w_0, ...
    'xs', 0, ...
    'ys', 0, ...
    't_0', t_offset);


% time interval and number of time samples
lambda = (2 * pi) / wlims(2);
tlims = [0, t_end];
numt = ceil((tlims(2) / lambda) * 15);

% ------------------------------------------------------------
% Solver params
% ------------------------------------------------------------
params = build_params( ...
    'k', w_0 + freq_gauss_width, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 160, ...
    'n', 20, ...
    'p', 3, ...
    'p_edge', 2);

% ------------------------------------------------------------
% Build all geometries
% ------------------------------------------------------------
num_presets = 5;

all_geometries = cell(1, 2 * num_presets);
all_geometry_names = cell(1, 2 * num_presets);

idx = 1;

% twinjet presets
for i = 1:num_presets
    preset_name = sprintf("plane%d", i);
    gparams = build_twinjet_params('preset', preset_name);

    all_geometries{idx} = twinjet_plane_geometry(gparams);
    all_geometry_names{idx} = sprintf("twinjet_%s", preset_name);

    idx = idx + 1;
end

% quadjet presets
for i = 1:num_presets
    preset_name = sprintf("plane%d", i);
    gparams = build_quadjet_params('preset', preset_name);

    all_geometries{idx} = quadjet_plane_geometry(gparams);
    all_geometry_names{idx} = sprintf("quadjet_%s", preset_name);

    idx = idx + 1;
end

% ------------------------------------------------------------
% Run experiment
% ------------------------------------------------------------
start = tic;
% for ii = 9:10
for ii = 1:length(all_geometries)
    geom = all_geometries{ii};

    fprintf("\n==================================================\n");
    fprintf("Running object %d / %d: %s\n", ii, ...
        length(all_geometries), all_geometry_names{ii});
    fprintf("==================================================\n");
    
    % Build a refined curve.
    C = build_curve(geom, params);
    C = refine_curve_wavenumber(C, 2*pi/params.k, 1);

    for jj = 1:1
        C = bisect_all_patches(C);
        C = build_global_indexing(C);
        C = get_close_points(C);
    end

    lp = RP2LP(params, geom);
    lp = lp.update_RP_Curve(C);

    % total number of points in the curve
    N = lp.RP_Curve.global.N

    % rule of thumb for computing complex resonances
    w_len = wlims(2) - wlims(1);
    num_sub_ints = ceil(w_len / 2.0);
    jump = w_len / num_sub_ints;

    pols_wp = [];
    pstart = tic;

    for jj = 1:num_sub_ints
        wlims_loc = [wlims(1) + (jj - 1) * jump, ...
                     wlims(1) + jj * jump];

        % dummy problem data for pole computation
        psp = problem_data( ...
            'wlims', wlims_loc, ...
            'w_0', (wlims_loc(2) + wlims_loc(1)) / 2.0, ...
            'numw', numw, ...
            'xs', 0, ...
            'ys', 0, ...
            'wimag', wimag, ...
            'N', N);

        pols_wp_loc = comp_poles(psp, lp);
        pols_wp_loc = sort(pols_wp_loc, ...
            "comparisonMethod", "real");
        pols_wp = [pols_wp; pols_wp_loc(:)];
    end
    pol_comp_time = toc(pstart)
   

   

    % dummy problem for LU decompositions
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

    if are_recivers_random
        [xs, ys] = ff_points_random(angle1, angle2, ...
            num_receivers_per_obs);
    else
        [xs, ys] = ff_points_equally_spaced(angle1, angle2, ...
            num_receivers_per_obs);
    end

    uff_all = [];
    r_res_all = [];
    f_sols_all = [];

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
        [uff, scat_sol, r_res, ~, smoothint, f_sols, zeroint] = cts_wpols_LU(psp, lp, pols_wp, ...
            L_rl, U_rl, L_pole, U_pole);
        single_reciever_comp_time = toc

        uff = uff(:).';
        r_res = r_res(:).';
        f_sols = f_sols(:).';

        uff_all = [uff_all; uff];
        r_res_all = [r_res_all; r_res];
        f_sols_all = [f_sols_all; f_sols];
    end

    ts = psp.ts;
    curveX = lp.curve.X;
    curveY = lp.curve.Y;
    pols = pols_wp;

    save(fullfile(save_data_dir, all_geometry_names{ii} + ".mat"), ...
        "r_res_all", "uff_all", "pols", "ts", "xs", "ys", ...
        "curveX", "curveY","f_sols_all",'psp');
end

total_time = toc(start)
