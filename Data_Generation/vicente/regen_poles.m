% regen_poles.m
%
% Regenerates a plane's pole file on the CURRENT curve discretization
% (matching plot_pole_eigenfunction.m / the current data_gen_twinjet.m /
% data_gen_quadjet.m exactly: curve_k=55, ppw=1, npolar=160, n=20, p=3,
% p_edge=2, delta=0.05, bisected 2x for twinjet / 1x for quadjet).
%
% Why this exists: the stored pole files under /scratch/vhojas/data_planes_2d
% (e.g. twinjet_plane1.mat, curveX/curveY with only N=552 points) predate
% the current data_gen_twinjet.m (that script's earliest commit postdates
% the .mat file's birth by over a month) and were built on a much coarser
% curve than the N=3520 curve plot_pole_eigenfunction.m rebuilds. A pole of
% a discretized BIE operator is a property of that specific discretization
% -- refining the mesh moves the true singularities -- so none of the
% loaded pole locations are actual singularities of the finer operator.
% Confirmed directly (Data_Generation/vicente/diag_pole.m, job 13112):
% rcond(bie_mat(w)) stays flat (~1e-6/1e-7) across 4 orders of magnitude
% of distance from the claimed pole 1 location instead of collapsing, and
% the raw boundary density doesn't blow up like 1/radius as a real
% residue would. This regenerates poles that are genuine singularities of
% the SAME curve plot_pole_eigenfunction.m/comp_r_res_grid.m use, via the
% same AAA-rational-approximation pole search data_gen_twinjet.m /
% data_gen_quadjet.m use (comp_poles) -- not a local secant refinement of
% a stale guess, which is exactly the kind of search that can converge to
% a spurious root of the random scalarized resolvent instead of a true
% pole when the matrix is only moderately (not exceptionally) ill
% conditioned everywhere nearby, which is what happened here.
%
% Deliberately skips the per-receiver far-field time-domain synthesis
% (359 receivers x cts_wpols_LU) that data_gen_twinjet.m/data_gen_quadjet.m
% also compute -- plot_pole_eigenfunction.m only ever reads "pols" back
% out of these files (S = load(pole_data_file, "pols")), so that synthesis
% is pure unneeded cost here.
%
% Writes to a NEW file (suffix "_regen") rather than overwriting the
% original, so the stale file stays available until the regenerated poles
% are validated against a real eigenfunction plot.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

vhojas_dir   = "/scratch/vhojas/data_planes_2d";
out_dir      = "/scratch/vhojas/data_planes_2d/regen";
if ~exist(out_dir, "dir")
    mkdir(out_dir);
end

configs = struct( ...
    "name", {}, "geometry_type", {}, "preset", {});
for i = 1:5
    configs(i) = struct("name", sprintf("quadjet_plane%d", i), ...
        "geometry_type", "quadjet", "preset", sprintf("plane%d", i));
end
for i = 1:5
    configs(5+i) = struct("name", sprintf("twinjet_plane%d", i), ...
        "geometry_type", "twinjet", "preset", sprintf("plane%d", i));
end

config_id = str2double(getenv("SLURM_ARRAY_TASK_ID"));
if isnan(config_id)
    config_id = str2double(getenv("CONFIG_ID"));
end
if isnan(config_id)
    config_id = 6;   % twinjet_plane1
    fprintf("No SLURM_ARRAY_TASK_ID/CONFIG_ID set; defaulting to config %d (%s).\n", ...
        config_id, configs(config_id).name);
end

config = configs(config_id);
plane_name    = config.name;
geometry_type = config.geometry_type;
plane_preset  = config.preset;

fprintf("Config %d/%d: %s (%s, preset %s)\n", ...
    config_id, numel(configs), plane_name, geometry_type, plane_preset);

%% ------------------------------------------------------------
%  Parallel workers -- same per-run JobStorageLocation pattern as
%  plot_pole_eigenfunction.m / late_time_nacelle_field.m.
% ------------------------------------------------------------

num_workers = str2double(getenv("SLURM_CPUS_PER_TASK"));
if isnan(num_workers) || num_workers < 1
    num_workers = maxNumCompThreads;
end
if isempty(gcp("nocreate"))
    job_tag = getenv("SLURM_JOB_ID");
    if isempty(job_tag)
        job_tag = sprintf("pid%d", feature("getpid"));
    end
    array_tag = getenv("SLURM_ARRAY_TASK_ID");
    if ~isempty(array_tag)
        job_tag = job_tag + "_" + array_tag;
    end
    pool_dir = fullfile(getenv("HOME"), ".matlab", "local_cluster_jobs", "per_run_" + job_tag);
    if ~exist(pool_dir, "dir")
        mkdir(pool_dir);
    end
    cluster = parcluster("local");
    cluster.JobStorageLocation = pool_dir;
    parpool(cluster, num_workers);
end

%% ------------------------------------------------------------
%  Experiment / curve params -- identical to data_gen_twinjet.m /
%  data_gen_quadjet.m's current (committed) values.
% ------------------------------------------------------------

w_0              = 50;
freq_gauss_width = 5;
gauss_picker_tol = 1e-10;
numw             = 400;
wimag            = -0.3;

[sigmas, wlims] = gauss_picker(w_0, freq_gauss_width, gauss_picker_tol);

params = build_params( ...
    'k', w_0 + freq_gauss_width, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 160, ...
    'n', 20, ...
    'p', 3, ...
    'p_edge', 2);

switch geometry_type
    case "twinjet"
        num_bisections = 2;
        gparams = build_twinjet_params('preset', plane_preset);
        geom = twinjet_plane_geometry(gparams);
    case "quadjet"
        num_bisections = 1;
        gparams = build_quadjet_params('preset', plane_preset);
        geom = quadjet_plane_geometry(gparams);
    otherwise
        error("Unknown geometry_type '%s'.", geometry_type);
end

C = build_curve(geom, params);
C = refine_curve_wavenumber(C, 2*pi/params.k, 1);
for jj = 1:num_bisections
    C = bisect_all_patches(C);
    C = build_global_indexing(C);
    C = get_close_points(C);
end

lp = RP2LP(params, geom);
lp = lp.update_RP_Curve(C);

N = lp.RP_Curve.global.N;
fprintf("N = %d boundary points (%d bisections)\n", N, num_bisections);

[xb, yb] = get_global_xy(C);
curveX = xb;
curveY = yb;

%% ------------------------------------------------------------
%  Pole search -- identical sub-interval sweep to data_gen_twinjet.m /
%  data_gen_quadjet.m (comp_poles per sub-interval, AAA-recursive
%  rational approximation + secant refinement -- a global search, not a
%  local refinement of a stale guess, so it can't silently land on a
%  spurious root of the scalarized resolvent the way re-refining an
%  already-wrong starting point can).
% ------------------------------------------------------------

w_len = wlims(2) - wlims(1);
num_sub_ints = ceil(w_len / 2.0);
jump = w_len / num_sub_ints;
fprintf("Frequency range [%.3f, %.3f], %d sub-intervals\n", ...
    wlims(1), wlims(2), num_sub_ints);

pols_wp = [];
pstart = tic;
for jj = 1:num_sub_ints
    wlims_loc = [wlims(1) + (jj - 1) * jump, wlims(1) + jj * jump];

    psp = problem_data( ...
        'wlims', wlims_loc, ...
        'w_0', (wlims_loc(2) + wlims_loc(1)) / 2.0, ...
        'numw', numw, ...
        'xs', 0, ...
        'ys', 0, ...
        'wimag', wimag, ...
        'N', N, ...
        'aaa_tol', 1e-8);

    sub_start = tic;
    pols_wp_loc = comp_poles(psp, lp);
    pols_wp_loc = sort(pols_wp_loc, "comparisonMethod", "real");
    pols_wp = [pols_wp; pols_wp_loc(:)];
    fprintf("  sub-interval %2d/%2d [%.3f, %.3f]: %d poles found (%.1f s elapsed, %.1f s total)\n", ...
        jj, num_sub_ints, wlims_loc(1), wlims_loc(2), numel(pols_wp_loc), ...
        toc(sub_start), toc(pstart));
end
pol_comp_time = toc(pstart);
fprintf("\nTotal pole search time: %.1f s. %d poles found.\n", pol_comp_time, numel(pols_wp));

pols = pols_wp;

%% ------------------------------------------------------------
%  Save -- new file, doesn't touch the original stale one.
% ------------------------------------------------------------

out_path = fullfile(out_dir, plane_name + "_regen.mat");
save(out_path, "pols", "curveX", "curveY", "N", "num_bisections", ...
    "geometry_type", "plane_preset", "plane_name", "wlims", "w_0", ...
    "freq_gauss_width", "wimag", "numw", "pol_comp_time");
fprintf("Saved regenerated poles to: %s\n", out_path);
