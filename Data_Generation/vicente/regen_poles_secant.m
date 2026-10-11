% regen_poles_secant.m
%
% Cheaper alternative to regen_poles.m's full comp_poles (AAA-recursive)
% global search, which OOM'd on 128 'Processes' workers each holding
% their own N=3520 dense matrix before finishing even the first
% sub-interval (job 13113: 7 oom_kill events at ~6 min).
%
% Instead: secant-refine each of the STALE (N=552-curve) pole guesses
% against the CURRENT N=3520 curve's BIE operator -- same idea
% plot_pole_eigenfunction.m already tried (and which produced the false
% "resolved" result diagnosed in diag_pole.m/job 13112) but with the gap
% that made that attempt unreliable closed:
%
%   1. plot_pole_eigenfunction.m secant-refines from ONE fixed random
%      projection f(w) = u*(bie_mat(w)\v), shared (rng(0)) across all 7
%      poles. If that one random (u,v) pair happens to weakly couple to
%      several poles' true null directions, secant can converge to a
%      spurious stationary point of THAT SPECIFIC scalar function rather
%      than a true pole of the operator -- and can do it for multiple
%      poles at once, which is consistent with all 7 poles failing
%      near-identically. Here each pole is refined from N_SEEDS
%      independent random projections and only accepted if they agree.
%   2. plot_pole_eigenfunction.m trusts a small secant residual
%      (tol=1e-13 on the update step) as proof of convergence to a real
%      pole. But the matrix here is only moderately ill-conditioned
%      everywhere nearby (rcond ~1e-6/1e-7 across 4 orders of magnitude
%      of distance from the claimed pole -- see diag_pole.m), not sharply
%      singular, so a small residual doesn't guarantee a genuine
%      near-singularity. Here every multi-seed-agreeing candidate is
%      additionally validated with the SAME physically-grounded contour-
%      density-blowup ratio check plot_pole_eigenfunction.m already uses
%      for its own sanity check (comparing solve magnitude at cont_rad
%      from the candidate vs. at the candidate's real part) before being
%      accepted, rather than trusting the secant residual alone.
%
% Any pole that fails multi-seed agreement or the ratio check is left
% out and reported, rather than silently kept -- it should fall back to
% a local AAA/rcond scan restricted to a small window (not the full
% regen_poles.m global sweep) in a follow-up, not something this script
% guesses at.

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

plane_name     = "twinjet_plane1";
geometry_type  = "twinjet";
plane_preset   = "plane1";
pole_data_file = fullfile(vhojas_dir, plane_name + ".mat");

%% ------------------------------------------------------------
%  Parallel workers -- modest count, matching the successful diag_pole.sh
%  resource profile (16 workers / 32G handled N=3520 solves fine) rather
%  than the 128-worker request that OOM'd.
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
%  Curve -- identical to plot_pole_eigenfunction.m / regen_poles.m.
% ------------------------------------------------------------

tilt_deg = 5;
kappa = -[cosd(tilt_deg), sind(tilt_deg)];

ppw = 1; npolar = 160; n_quad = 20; p_corner = 3; p_edge = 2; delta_bc = 0.05;
curve_k = 55;
num_bisections = 2;

gparams = build_twinjet_params('preset', plane_preset);
geom = twinjet_plane_geometry(gparams);

params = build_params('k', curve_k, 'formulation', 'sl', 'solver', 'direct', ...
    'delta', delta_bc, 'npolar', npolar, 'n', n_quad, 'p', p_corner, 'p_edge', p_edge);

C = build_curve(geom, params);
C = refine_curve_wavenumber(C, 2*pi/params.k, ppw);
for jj = 1:num_bisections
    C = bisect_all_patches(C);
    C = build_global_indexing(C);
    C = get_close_points(C);
end

lp = RP2LP(params, geom);
lp.RP_Curve = C;
[xb, yb] = get_global_xy(C);
lp.curve.X = xb; lp.curve.Y = yb; lp.curve.Flag = "vicente open";
bndpts = [xb, yb]';
N = lp.RP_Curve.global.N;
fprintf("N = %d boundary points\n", N);

%% ------------------------------------------------------------
%  Load stale (N=552-curve) pole guesses.
% ------------------------------------------------------------

S = load(pole_data_file, "pols");
pols_stale = S.pols(:);
fprintf("Loaded %d stale pole guesses from %s\n", numel(pols_stale), pole_data_file);

%% ------------------------------------------------------------
%  Multi-seed secant refinement + ratio-check validation per pole.
% ------------------------------------------------------------

n_seeds        = 5;
seed_agree_tol = 1e-3;   % max pairwise |location diff| across seeds to accept agreement
num_its        = 10;     % more than plot_pole_eigenfunction.m's default 4, without
                          % running secant_method's full iteration budget needlessly
                          % (each iteration is a full N=3520 dense solve)
sec_tol        = 1e-10;

cont_rad       = 1e-5;
n_res_w        = 10;
ratio_threshold = 1e2;   % same threshold used throughout this codebase

pols_refined = nan(size(pols_stale));
accepted     = false(size(pols_stale));
ratios       = nan(size(pols_stale));
seed_spreads = nan(size(pols_stale));

for j = 1:numel(pols_stale)
    p0 = pols_stale(j);
    fprintf("\nPole %d/%d: stale guess k = %.6f%+.6fi\n", j, numel(pols_stale), real(p0), imag(p0));

    seed_ests = complex(zeros(n_seeds, 1));
    for s = 1:n_seeds
        rng(s);
        u = rand(1, N, 'like', 1i);
        v = rand(N, 1, 'like', 1i);
        f = @(w) u * (lp.bie_mat(w) \ v);
        [x, ~, ~, x_diff] = secant_method(@(z) 1/f(z), p0, p0 + 1e-5, num_its, sec_tol, 1);
        seed_ests(s) = x;
        fprintf("  seed %d: refined to k = %.6f%+.6fi (final step %.3e)\n", ...
            s, real(x), imag(x), abs(x_diff));
    end

    spread = max(abs(seed_ests - mean(seed_ests)));
    seed_spreads(j) = spread;
    if spread > seed_agree_tol
        fprintf("  [REJECTED] seeds disagree (spread = %.3e > %.1e) -- no consistent pole found\n", ...
            spread, seed_agree_tol);
        continue
    end

    candidate = mean(seed_ests);
    if imag(candidate) >= 0
        fprintf("  [REJECTED] candidate k = %.6f%+.6fi is not in the lower half-plane\n", ...
            real(candidate), imag(candidate));
        continue
    end

    % Ratio-check validation -- same physically-grounded test used
    % throughout this codebase (comp_r_res_grid.m's header, plot_pole_
    % eigenfunction.m:266, late_time_nacelle_field.m:566-583): does the
    % boundary density genuinely blow up near this location relative to
    % a generic nearby real-line solve.
    circ = candidate + cont_rad*exp(1i*2*pi*(1:n_res_w)/n_res_w);
    res_dens_j = zeros(N, n_res_w);
    for cind = 1:n_res_w
        rhs = -exp(1i*circ(cind)*(kappa*bndpts)).';
        res_dens_j(:, cind) = lp.bie_mat(circ(cind)) \ rhs;
    end
    rhs_real = -exp(1i*real(candidate)*(kappa*bndpts)).';
    dens_real = lp.bie_mat(real(candidate)) \ rhs_real;
    ratio = max(abs(res_dens_j(:))) / max(abs(dens_real(:)));
    ratios(j) = ratio;

    fprintf("  candidate k = %.6f%+.6fi (seed spread %.3e), contour/real-line ratio = %.3e   %s\n", ...
        real(candidate), imag(candidate), spread, ratio, ...
        ternary(ratio > ratio_threshold, "[ACCEPTED]", "[REJECTED -- not resolved]"));

    if ratio > ratio_threshold
        pols_refined(j) = candidate;
        accepted(j) = true;
    end
end

n_accepted = sum(accepted);
fprintf("\n%d/%d poles refined and validated.\n", n_accepted, numel(pols_stale));
if n_accepted < numel(pols_stale)
    fprintf("Poles NOT resolved this way (need a local AAA/rcond-scan fallback): %s\n", ...
        mat2str(find(~accepted).'));
end

pols = pols_refined(accepted);

out_path = fullfile(out_dir, plane_name + "_secant_regen.mat");
save(out_path, "pols", "pols_stale", "accepted", "ratios", "seed_spreads", ...
    "xb", "yb", "N", "num_bisections", "geometry_type", "plane_preset", "plane_name");
fprintf("Saved to: %s\n", out_path);

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
