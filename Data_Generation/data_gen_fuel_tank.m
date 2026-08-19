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
% Check for a previous (possibly partial) run to resume from
% ------------------------------------------------------------
resume_run = isfile(out_file);
existing_vars = {};
if resume_run
    mf_probe = matfile(out_file);       % read-only probe
    existing_vars = who(mf_probe);
    fprintf('Found existing output file %s. Attempting to resume.\n', out_file);
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

% ------------------------------------------------------------
% Poles: reuse from disk if this is a resumed run, otherwise compute
% ------------------------------------------------------------
if resume_run && ismember('pols', existing_vars)
    pols_wp = mf_probe.pols;
    pols_wp = pols_wp(:);
    pol_comp_time = NaN;  % not recomputed this run
    fprintf('Loaded %d previously computed poles from disk.\n', numel(pols_wp));
else
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
end

% ------------------------------------------------------------
% Save poles (and other static metadata) immediately, before the
% (expensive/large) LU and receiver-loop work below, but only if
% they aren't already saved from a previous run.
% ------------------------------------------------------------
pols = pols_wp;
config_name = "twinjet_with_fueltank";
receiver_angles = linspace(angle1, angle2, num_receivers_per_obs);
curveX = lp.curve.X;
curveY = lp.curve.Y;

if ~resume_run || ~ismember('pols', existing_vars)
    save(out_file, "pols", "pol_comp_time", "config_name", ...
        "receiver_angles", "curveX", "curveY", "-v7.3");
end

% dummy problem for LU decompositions
psp = problem_data( ...
    'wlims', wlims, ...
    'w_0', w_0, ...
    'numw', numw, ...
    'xs', 0, ...
    'ys', 0, ...
    'wimag', wimag, ...
    'N', N);

% NOTE: L_rl/U_rl/L_pole/U_pole are not persisted to disk (they're
% large and cheap-ish relative to the receiver loop), so these are
% recomputed every run, resumed or not.
tic
[L_rl, U_rl] = comp_bie_LU(psp.ws, lp);
LU_time_real_line = toc

tic
[L_pole, U_pole] = comp_pole_LU(psp, lp, pols_wp);
LU_time_pols = toc

% ------------------------------------------------------------
% Receiver coordinates: reuse from disk if resuming, so the
% receiver set is identical across runs; otherwise generate fresh.
% ------------------------------------------------------------
if resume_run && ismember('xs', existing_vars) && ismember('ys', existing_vars)
    xs = mf_probe.xs;
    ys = mf_probe.ys;
    fprintf('Loaded existing receiver coordinates from disk.\n');
else
    if are_recivers_random
        [xs, ys] = ff_points_random(angle1, angle2, ...
            num_receivers_per_obs);
    else
        [xs, ys] = ff_points_equally_spaced(angle1, angle2, ...
            num_receivers_per_obs);
    end
end

% open (or create) the matfile in writable mode
mf = matfile(out_file, "Writable", true);
mf.xs = xs;
mf.ys = ys;

% ------------------------------------------------------------
% Determine how many receivers are already computed and saved,
% so the loop below can resume from that point.
% ------------------------------------------------------------
rows_written = 0;
ts_saved = false;

if resume_run && ismember('uff_all', existing_vars)
    info = whos(mf_probe, 'uff_all');
    rows_written = info.size(1);
    rows_written = min(rows_written, num_receivers_per_obs);
    fprintf('Resuming: %d of %d receivers already computed.\n', ...
        rows_written, num_receivers_per_obs);
end

if resume_run && ismember('ts', existing_vars)
    ts_saved = true;
end

if rows_written >= num_receivers_per_obs
    fprintf('All %d receivers already computed. Nothing to do.\n', ...
        num_receivers_per_obs);
else

% ------------------------------------------------------------
% Receiver loop with buffered disk flushing every save_every
% receivers, so uff_all/r_res_all/f_sols_all never fully live
% in RAM at once. Starts from rows_written+1 when resuming.
% ------------------------------------------------------------
uff_buf    = {};
r_res_buf  = {};
f_sols_buf = {};

for rind = (rows_written + 1):num_receivers_per_obs
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
    if mod(rind - rows_written, save_every) == 0 || is_last
        block_uff    = cell2mat(uff_buf(:));
        block_r_res  = cell2mat(r_res_buf(:));
        block_f_sols = cell2mat(f_sols_buf(:));

        n_new = size(block_uff, 1);
        row_start = rind - n_new + 1;
        row_end   = rind;

        % Get column dimensions explicitly
        n_cols_uff    = size(block_uff, 2);
        n_cols_r_res  = size(block_r_res, 2);
        n_cols_f_sols = size(block_f_sols, 2);

        % Write using explicit ranges instead of ':'
        mf.uff_all(row_start:row_end, 1:n_cols_uff)       = block_uff;
        mf.r_res_all(row_start:row_end, 1:n_cols_r_res)   = block_r_res;
        mf.f_sols_all(row_start:row_end, 1:n_cols_f_sols) = block_f_sols;

        % Free the buffers and reclaim memory
        uff_buf = {}; r_res_buf = {}; f_sols_buf = {};
        clear block_uff block_r_res block_f_sols
    end
end

end % if rows_written >= num_receivers_per_obs

total_time = toc(start)