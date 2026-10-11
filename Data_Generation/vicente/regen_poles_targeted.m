% regen_poles_targeted.m
%
% regen_poles.m's global sweep over [45,55] (5 windows of width 2,
% default 50x50 AAA sampling grid) found only 2 near-duplicate poles at
% k~45.92-0.283i (Q~81), and ZERO poles anywhere in [47,55] -- the entire
% band the old (stale, N=552-curve) pole file placed all 7 of its poles
% in, including the "dominant" one at k~49.8-0.015i (claimed Q~1663).
%
% Likely cause: the 50x50 grid is fixed resolution regardless of window
% size, so a width-2 window gives ~0.04 real-direction spacing -- coarser
% than the ~0.015 influence radius a Q~1663 pole (Im~-0.015) would have.
% A lower-Q pole (Im~-0.28, wider influence) got caught fine in [45,47];
% a genuinely high-Q one near the real axis can fall between grid
% samples and never trigger aaa_recursive's adaptive refinement.
%
% Fix: search a much NARROWER window directly around where the old
% pole file placed its highest-Q candidates, at the SAME 50x50 grid
% resolution -- this gives much finer effective real/imag spacing for
% free, without touching aaa_tol or p_numx/p_numy, and costs about the
% same as one of the empty width-2 sub-intervals we just measured
% (~1400-1500s), not a full resweep at higher resolution.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

vhojas_dir = "/scratch/vhojas/data_planes_2d";
out_dir    = "/scratch/vhojas/data_planes_2d/regen";
if ~exist(out_dir, "dir")
    mkdir(out_dir);
end

plane_name    = "twinjet_plane1";
geometry_type = "twinjet";
plane_preset  = "plane1";

%% ------------------------------------------------------------
%  Parallel workers -- same profile that ran regen_poles.m cleanly
%  (32 workers / 96G, no OOM).
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
    pool_dir = fullfile(getenv("HOME"), ".matlab", "local_cluster_jobs", "per_run_" + job_tag);
    if ~exist(pool_dir, "dir")
        mkdir(pool_dir);
    end
    cluster = parcluster("local");
    cluster.JobStorageLocation = pool_dir;
    parpool(cluster, num_workers);
end

%% ------------------------------------------------------------
%  Curve -- identical to regen_poles.m / plot_pole_eigenfunction.m.
% ------------------------------------------------------------

params = build_params( ...
    'k', 55, 'formulation', 'sl', 'solver', 'direct', ...
    'delta', 0.05, 'npolar', 160, 'n', 20, 'p', 3, 'p_edge', 2);

num_bisections = 2;
gparams = build_twinjet_params('preset', plane_preset);
geom = twinjet_plane_geometry(gparams);

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
fprintf("N = %d boundary points\n", N);

[xb, yb] = get_global_xy(C);
curveX = xb;
curveY = yb;

%% ------------------------------------------------------------
%  Targeted, narrow-window search around the OLD "dominant" pole's
%  location only (k~49.8-0.015i, the highest-Q/closest-to-real-axis
%  stale guess, and the one the user actually wants an eigenfunction
%  for) -- width 1 in the real direction (vs. regen_poles.m's width 2)
%  and depth restricted to [-0.1, 0] (vs. the full [-0.3, 0] the global
%  sweep used, which is far deeper than a Q~1663 pole needs and buys
%  nothing but coarser effective resolution). Both narrowing changes
%  give a much finer effective grid at the SAME 50x50 resolution,
%  targeted right where this pole's resonance should be, rather than
%  re-searching the other 6 stale poles' locations (mostly Im down to
%  -0.2, well outside this shallow window, and not today's ask).
% ------------------------------------------------------------

S = load(fullfile(vhojas_dir, plane_name + ".mat"), "pols");
[~, dominant_idx] = min(abs(imag(S.pols(:))));
dominant_stale = S.pols(dominant_idx);
fprintf("Dominant stale guess: k = %.6f%+.6fi\n", real(dominant_stale), imag(dominant_stale));

half_width = 0.5;
merged = [real(dominant_stale) - half_width, real(dominant_stale) + half_width];
fprintf("Searching window [%.3f, %.3f]\n", merged(1), merged(2));

wimag_targeted = -0.1;
numw = 400;

pols_wp = [];
pstart = tic;
for jj = 1:size(merged, 1)
    wlims_loc = merged(jj, :);

    psp = problem_data( ...
        'wlims', wlims_loc, ...
        'w_0', mean(wlims_loc), ...
        'numw', numw, ...
        'xs', 0, ...
        'ys', 0, ...
        'wimag', wimag_targeted, ...
        'N', N, ...
        'aaa_tol', 1e-8);

    sub_start = tic;
    pols_wp_loc = comp_poles(psp, lp);
    pols_wp_loc = sort(pols_wp_loc, "comparisonMethod", "real");
    pols_wp = [pols_wp; pols_wp_loc(:)];
    fprintf("  window %d/%d [%.3f, %.3f]: %d poles found (%.1f s elapsed, %.1f s total)\n", ...
        jj, size(merged, 1), wlims_loc(1), wlims_loc(2), numel(pols_wp_loc), ...
        toc(sub_start), toc(pstart));
end
pol_comp_time = toc(pstart);
fprintf("\nTotal targeted search time: %.1f s. %d poles found.\n", pol_comp_time, numel(pols_wp));

pols = pols_wp;

out_path = fullfile(out_dir, plane_name + "_targeted_regen.mat");
save(out_path, "pols", "curveX", "curveY", "N", "num_bisections", ...
    "geometry_type", "plane_preset", "plane_name", "merged", ...
    "wimag_targeted", "numw", "pol_comp_time");
fprintf("Saved to: %s\n", out_path);
