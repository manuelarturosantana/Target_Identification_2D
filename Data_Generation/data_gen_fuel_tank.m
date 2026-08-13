
% generates data for 5 twinjet
% Now we use a refined curve to see if it makes the data better.

% paths
addpath(genpath('/home/msantana/Progamming/rp2d_matlab/src'))
addpath(genpath('/home/msantana/Progamming/Target_Identification_2D'))
addpath(genpath('/home/msantana/Progamming/MATLAB_PACKAGES'))

% ------------------------------------------------------------
% Output directory
% ------------------------------------------------------------
save_data_dir = "./plane_data";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

out_file = fullfile(save_data_dir, "twinjet1_with_fuel_tank.mat");

% ------------------------------------------------------------
% Experiment params
% ------------------------------------------------------------
w_0              = 50;            % center frequency of the Gaussian
freq_gauss_width = 5;             % Gaussian frequency width
gauss_picker_tol = 1e-10;         % Gaussian truncation tolerance
numw             = 400;           % number of points for frequency integrals
t_end            = 400;           % final time
t_offset         = 30;            % time offset
num_receivers_per_obs = 359;
are_recivers_random = false;

% how many receivers to buffer in RAM before flushing to disk
save_every = 20;

% receiver angles
angle1 = 0;
angle2 = 2*pi - 0.0174533; % This ensures that we don't double count the 360 degrees

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
% Build the geometry
% ------------------------------------------------------------
gparams = build_twinjet_fuel_tank_params('preset', 'plane1', ...
    'fuel_tank_loc', 0.45, ...
    'fuel_tank_length', 0.15, ...
    'fuel_tank_width', 0.006, ...
    'fuel_tank_shape1', 0.045, ...
    'fuel_tank_shape2', 0.0, ...
    'fuel_tank_shape3', 0.0, ...
    'fuel_tank_pillness', 0.5);

geom = twinjet_fuel_tank_geometry(gparams);


% ------------------------------------------------------------
% Run experiment
% ------------------------------------------------------------
start = tic;

% Build a refined curve.
C = build_curve(geom, params);
C = refine_curve_wavenumber(C, 2*pi/params.k, 1);

for jj = 1:2
    C = bisect_all_patches(C);
    C = build_global_indexing(C);
    C = get_close_points(C);
end

lp = RP2LP(params, geom);
lp = lp.update_RP_Curve(C);

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

    %dummy problem data for pole computation
    psp = problem_data( ...
        'wlims', wlims_loc, ...
        'w_0', (wlims_loc(2) + wlims_loc(1)) / 2.0, ...
        'numw', numw, ...
        'xs', 0, ...
        'ys', 0, ...
        'wimag', wimag, ...
        'N', N, ...
        'aaa_tol', 1e-8);

    pols_wp_loc = comp_poles(psp, lp);
    pols_wp_loc = sort(pols_wp_loc, ...
        "comparisonMethod", "real");
    pols_wp = [pols_wp; pols_wp_loc(:)];
end
pol_comp_time = toc(pstart)

% ------------------------------------------------------------
% Save poles immediately, before the (expensive/large) LU and
% receiver-loop work below. If anything downstream runs out of
% memory or crashes, the pole computation is not lost.
% ------------------------------------------------------------
pols = pols_wp;
config_name = "twinjet_with_fueltank";
receiver_angles = linspace(angle1, angle2, num_receivers_per_obs);
curveX = lp.curve.X;
curveY = lp.curve.Y;

save(out_file, "pols", "pol_comp_time", "config_name", ...
    "receiver_angles", "curveX", "curveY", "-v7.3");

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

% append receiver coordinates now that they're known
mf = matfile(out_file, "Writable", true);
mf.xs = xs;
mf.ys = ys;

% ------------------------------------------------------------
% Receiver loop with buffered disk flushing every save_every
% receivers, so uff_all/r_res_all/f_sols_all never fully live
% in RAM at once.
% ------------------------------------------------------------
uff_buf    = {};
r_res_buf  = {};
f_sols_buf = {};
rows_written = 0;   % how many receiver rows already flushed to disk
ts_saved = false;

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

    % save ts once (identical across receivers)
    if ~ts_saved
        mf.ts = psp.ts;
        ts_saved = true;
    end

    uff_buf{end+1}    = uff(:).';
    r_res_buf{end+1}  = r_res(:).';
    f_sols_buf{end+1} = f_sols(:).';

    is_last = (rind == num_receivers_per_obs);
    if mod(rind, save_every) == 0 || is_last
        block_uff    = cell2mat(uff_buf(:));
        block_r_res  = cell2mat(r_res_buf(:));
        block_f_sols = cell2mat(f_sols_buf(:));

        row_start = rows_written + 1;
        row_end   = rows_written + size(block_uff, 1);

        mf.uff_all(row_start:row_end, :)    = block_uff;
        mf.r_res_all(row_start:row_end, :)  = block_r_res;
        mf.f_sols_all(row_start:row_end, :) = block_f_sols;

        rows_written = row_end;

        % free the buffers and force MATLAB to reclaim memory
        uff_buf = {}; r_res_buf = {}; f_sols_buf = {};
        clear block_uff block_r_res block_f_sols
    end
end

total_time = toc(start)
