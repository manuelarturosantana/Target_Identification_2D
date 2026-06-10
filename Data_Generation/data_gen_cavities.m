% paths
addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/chebfun"))

% ------------------------------------------------------------
% Output directory
% ------------------------------------------------------------
save_data_dir = "/scratch/vhojas/data_cavities_2d";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

% ------------------------------------------------------------
% Experiment params
% ------------------------------------------------------------
w_0              = 10;            % center frequency of the Gaussian
freq_gauss_width = 5;             % Gaussian frequency width
gauss_picker_tol = 1e-10;         % Gaussian truncation tolerance
numw             = 100;           % number of points for frequency integrals
t_end            = 400;           % final time
t_offset         = 30;            % time offset
num_receivers_per_obs = 100;
are_recivers_random = true;
angle1 = pi;
angle2 = 2 * pi;
wimag  = -0.5;

% time-domain and AAA derived params
[sigmas, wlims] = gauss_picker(w_0, freq_gauss_width, gauss_picker_tol);
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
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

% ------------------------------------------------------------
% Build all geometries
% ------------------------------------------------------------
openings = 0.125:0.125:2.5;

all_geometries = cell(1, 2 * numel(openings));
all_geometry_names = cell(1, 2 * numel(openings));

idx = 1;

for i = 1:numel(openings)
    opening = openings(i);

    % elliptic cavity
    all_geometries{idx} = ellipticcavity(2, 1, opening);
    all_geometry_names{idx} = sprintf("elliptic_%0.2f", ...
        opening);
    idx = idx + 1;

    % circular cavity
    all_geometries{idx} = circularcavity(1, opening);
    all_geometry_names{idx} = sprintf("circular_%0.2f", ...
        opening);
    idx = idx + 1;
end

% ------------------------------------------------------------
% Run experiment
% ------------------------------------------------------------
start = tic;

for ii = 1:length(all_geometries)
    geom = all_geometries{ii};

    fprintf("\n==================================================\n");
    fprintf("Running object %d / %d: %s\n", ii, length(all_geometries), ...
        all_geometry_names{ii});
    fprintf("==================================================\n");

    lp = RP2LP(params, geom);

    % total number of points in the curve
    N = lp.RP_Curve.global.N;

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
        pols_wp_loc = sort(pols_wp_loc, "comparisonMethod", "real");
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
        [uff, ~, r_res] = cts_wpols_LU(psp, lp, pols_wp, ...
            L_rl, U_rl, L_pole, U_pole);
        single_reciever_comp_time = toc

        uff = uff(:).';
        r_res = r_res(:).';

        uff_all = [uff_all; uff];
        r_res_all = [r_res_all; r_res];
    end

    ts = psp.ts;
    curveX = lp.curve.X;
    curveY = lp.curve.Y;
    pols = pols_wp;

    save(fullfile(save_data_dir, all_geometry_names{ii} + ".mat"), ...
        "r_res_all", "uff_all", "pols", "ts", "xs", "ys", ...
        "curveX", "curveY");
end

total_time = toc(start)