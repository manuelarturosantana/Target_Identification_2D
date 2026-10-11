% late_time_nacelle_field.m
%
% Renders a paper-quality snapshot of the *late-time* total field around a
% single plane geometry, illuminated nose-on by a long (narrowband)
% Gaussian pulse tuned to a single trapped nacelle resonance. The image
% is produced by evaluating the pole/residue singularity-expansion
% time-domain synthesis (Anderson-Bruno-Lyon) at exactly one requested
% time -- comp_time_sol.m / freq_to_time.m / sum_SEM.m in TD_Solver_Code
% evaluate each requested time independently (no time-marching), so we
% only ever solve for the single late time we plot, never the times
% before it.
%
% Poles are NOT searched for here -- they are loaded from an existing
% pole-data file already produced by one of the data_gen_* scripts, so
% the (expensive, iterative) AAA pole search is skipped entirely. Only
% the frequency-domain BIE solves and near-field evaluations needed for
% THIS grid of image points and THIS late time are recomputed.
%
% The near-field evaluation itself does not reuse TD_Solver_Code's
% comp_r_res.m/comp_freq_sols.m as-is: those loop over spatial points one
% at a time, which is fine for a handful of far-field receivers but is
% far too slow for an image grid. comp_r_res_grid.m/comp_freq_sols_grid.m
% (same directory) implement identical math but call rp2d_matlab's
% (already point-array-vectorized) evaluate_helmholtz_potential once per
% frequency across the whole grid, exactly as plot_scattering_field.m
% does for a single CW frequency.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

%% ------------------------------------------------------------
%  Resolution -- start low (preview) while iterating, bump up for the
%  final paper-quality render. ppw/num_bisections control the boundary
%  discretization; ng controls the image grid; n/npolar control the
%  per-patch quadrature order.
% ------------------------------------------------------------
% Exactly the curve-discretization parameters
% Data_Generation/data_gen_twinjet.m and data_gen_quadjet.m use
% to find pols_wp in the first place (ppw=1, 1 bisection,
% npolar=160, n=20, p=3) -- NOT a hand-picked "in-between" tier.
% numw is left at preview's value: numw governs the time-domain
% Fourier/Chebyshev reconstruction accuracy, a separate question
% from curve/BIE resolution and not part of this comparison.
ppw            = 1;
num_bisections = 1;
n_quad         = 20;
npolar         = 160;
ng             = 120;
numw           = 40;
p_corner       = 3;
check_convergence = false;

%% ------------------------------------------------------------
%  Target geometry / pole data
%
%  One config per plane in the classification experiment (numbered to
%  match Paper_Figure_Generator's convention: 1-5 quadjet, 6-10
%  twinjet), selected via SLURM_ARRAY_TASK_ID/CONFIG_ID so the same
%  script drives the whole batch as a job array, exactly like
%  data_gen_deformation.m/data_gen_constellation.m do. Falls back to
%  config 6 (twinjet_plane1, the one already validated by hand) when
%  run standalone with neither env var set.
%
%  quadjet_plane5 is missing from vhojas's own
%  /scratch/vhojas/data_planes_2d (only has quadjet 1-4); poles are an
%  intrinsic property of the geometry, not of who generated them, so
%  that one config alone falls back to msantana's complete set.
% ------------------------------------------------------------

vhojas_dir   = "/scratch/vhojas/data_planes_2d";
msantana_dir = "/scratch/msantana/plane_data_360";

configs = struct( ...
    "name", {}, "geometry_type", {}, "preset", {}, "pole_data_file", {});
for i = 1:5
    configs(i) = struct( ...
        "name", sprintf("quadjet_plane%d", i), ...
        "geometry_type", "quadjet", ...
        "preset", sprintf("plane%d", i), ...
        "pole_data_file", fullfile(vhojas_dir, sprintf("quadjet_plane%d.mat", i)));
end
configs(5).pole_data_file = fullfile(msantana_dir, "quadjet_plane5.mat");
for i = 1:5
    configs(5+i) = struct( ...
        "name", sprintf("twinjet_plane%d", i), ...
        "geometry_type", "twinjet", ...
        "preset", sprintf("plane%d", i), ...
        "pole_data_file", fullfile(vhojas_dir, sprintf("twinjet_plane%d.mat", i)));
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
if config_id < 1 || config_id > numel(configs)
    error("Config id %d is outside 1:%d.", config_id, numel(configs));
end

config = configs(config_id);
plane_name      = config.name;            % e.g. "twinjet_plane1", for titles/filenames
geometry_type   = config.geometry_type;    % "twinjet" or "quadjet"
plane_preset    = config.preset;           % e.g. "plane1", for build_*_params
pole_data_file  = config.pole_data_file;
target_pole_idx = [];   % [] -> auto-pick the highest-Q pole in the file

% "pole_match" is only an accurate match to the curve that found pols_wp
% if num_bisections is ALSO matched per geometry type: Data_Generation/
% data_gen_twinjet.m bisects 2x, but data_gen_quadjet.m only bisects 1x
% (confirmed by direct inspection of both scripts) -- everything else
% (ppw, n, npolar, p_corner) is identical between the two. Getting this
% wrong silently breaks pole-residue extraction: comp_r_res_grid.m
% circles the loaded pole location with a cont_rad=1e-5 contour, and a
% pole found on a materially different discretization is generally
% shifted by much more than 1e-5, so the contour integral (which is
% exact for a genuine simple pole) instead integrates a regular,
% non-singular function and returns ~0 -- not an error, just silent
% numerical noise standing in for what should be the dominant resonant
% signal.
if resolution_level == "pole_match"
    switch geometry_type
        case "twinjet"
            num_bisections = 2;
        case "quadjet"
            num_bisections = 1;
    end
end

fprintf("Config %d/%d: %s (%s, preset %s)\n", ...
    config_id, numel(configs), plane_name, geometry_type, plane_preset);
fprintf("Pole data: %s\n", pole_data_file);

tilt_deg = 5;   % incidence tilt off the +x axis, breaks the mirror
                % symmetry that would otherwise leave degenerate
                % nacelle-pair modes unexcited at exact normal incidence

gauss_picker_tol = 1e-10;   % same tolerance used throughout this repo's
                             % gauss_picker calls -- reused below to
                             % define the incident pulse's time-domain
                             % shutoff, not just its frequency window

save_dir = fullfile(fileparts(mfilename('fullpath')), "figures");
if ~exist(save_dir, "dir")
    mkdir(save_dir);
end

%% ------------------------------------------------------------
%  Parallel workers -- the per-frequency/per-pole BIE solves below are
%  each independent and already parfor-ed inside comp_bie_LU.m,
%  comp_pole_LU.m and comp_res_dens.m. Point this at the CPUs SLURM gave
%  us (falls back to all local cores off-cluster).
% ------------------------------------------------------------

num_workers = str2double(getenv("SLURM_CPUS_PER_TASK"));
if isnan(num_workers) || num_workers < 1
    num_workers = maxNumCompThreads;
end
if isempty(gcp("nocreate"))
    % Give this job its own JobStorageLocation: parpool's default lives
    % under $HOME/.matlab/local_cluster_jobs, which is shared (NFS) and
    % identical for every concurrent run of this script -- e.g. two
    % array tasks landing on different nodes at once -- causing
    % intermittent "Failed to open file ... Job*.in.mat" startup
    % failures when they collide on the same job-id file.
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
%  Load already-computed poles for this geometry (skip pole search)
% ------------------------------------------------------------

S = load(pole_data_file, "pols");
pols_wp = S.pols(:);

Q = abs(real(pols_wp))./(2*abs(imag(pols_wp)));
fprintf("Loaded %d poles from %s\n", numel(pols_wp), pole_data_file);
for j = 1:numel(pols_wp)
    fprintf("   pole %2d: k = %8.4f %+8.4fi   Q ~ %10.2f\n", ...
        j, real(pols_wp(j)), imag(pols_wp(j)), Q(j));
end

if isempty(target_pole_idx)
    [~, target_pole_idx] = max(Q);
end
target_pole = pols_wp(target_pole_idx);

fprintf("\nTarget pole %d: k = %.6f %+.6fi  (Q ~ %.1f)\n", ...
    target_pole_idx, real(target_pole), imag(target_pole), ...
    Q(target_pole_idx));

% Excite the target pole selectively: pick the pulse's frequency
% half-width so the envelope has already decayed to gauss_picker_tol at
% the *nearest other pole*, leaving that neighbor essentially unexcited.
other_re = real(pols_wp);
other_re(target_pole_idx) = [];
if isempty(other_re)
    freq_gauss_width = 1.0;
else
    freq_gauss_width = min(abs(other_re-real(target_pole)));
end

% Guard against a pathologically small width if the target pole happens
% to sit right next to a near-degenerate twin (e.g. a symmetric
% nacelle-pair split by a tiny numerical amount) -- an extremely narrow
% width blows up the pulse duration/t_late without meaningfully
% improving mode isolation in that case.
min_freq_gauss_width = 0.05;
if freq_gauss_width < min_freq_gauss_width
    fprintf(['Nearest neighboring pole is only %.4f away; clamping ' ...
        'freq_gauss_width to %.4f.\n'], freq_gauss_width, min_freq_gauss_width);
    freq_gauss_width = min_freq_gauss_width;
end

w_0 = real(target_pole);

fprintf("freq_gauss_width = %.4f (distance to nearest neighboring pole)\n", ...
    freq_gauss_width);

%% ------------------------------------------------------------
%  Incident pulse timing
%
%  Delta below is the pulse's time-domain half-width down to the same
%  gauss_picker_tol used for its frequency window (Gaussians are their
%  own Fourier dual, so this is the exact time-domain mirror of
%  gauss_picker's frequency-domain "dist"). t_offset keeps the same
%  ratio to Delta that the reference plane-generation script
%  (Data_Generation/data_gen_twinjet.m: t_offset=30, freq_gauss_width=5,
%  tol=1e-10 => Delta_ref ~= 9.21) used, so the pulse is negligible at
%  t=0 by the same relative margin regardless of how much longer our
%  narrower-band pulse is. t_late then adds a transit buffer (computed
%  below, once the image grid is known) equal to the time needed for a
%  wavefront to cross the full image domain at unit speed, so the
%  direct/geometric part of the field has fully cleared the plotted
%  region and only the resonant tail remains.
% ------------------------------------------------------------

delta_t = 2*abs(log(gauss_picker_tol))/freq_gauss_width;

t_offset_ref = 30;
freq_width_ref = 5;
delta_t_ref = 2*abs(log(gauss_picker_tol))/freq_width_ref;
startup_ratio = t_offset_ref/delta_t_ref;

t_offset = startup_ratio*delta_t;

fprintf("\nPulse timing: delta_t = %.3f, t_offset = %.3f\n", delta_t, t_offset);

%% ------------------------------------------------------------
%  Geometry + curve
% ------------------------------------------------------------

switch geometry_type
    case "twinjet"
        gparams = build_twinjet_params('preset', plane_preset);
        geom = twinjet_plane_geometry(gparams);
    case "quadjet"
        gparams = build_quadjet_params('preset', plane_preset);
        geom = quadjet_plane_geometry(gparams);
    otherwise
        error("Unknown geometry_type '%s'.", geometry_type);
end

% The curve must be built at exactly the wavenumber
% data_gen_twinjet.m/data_gen_quadjet.m used for refine_curve_wavenumber,
% or the patch layout differs from the one pols_wp was actually found on
% even with ppw/num_bisections/n/npolar/p all matched (confirmed: with
% only those matched, the pole-contour density came out statistically
% indistinguishable from the real-line density -- i.e. still missing the
% true pole). Both data_gen scripts hardcode w_0=50, freq_gauss_width=5
% for this purpose regardless of preset or geometry type, giving k=55 --
% NOT this script's own (much narrower) target-pole-based w_0/
% freq_gauss_width, which exists only to shape the excitation pulse and
% must not leak into the curve construction.
if resolution_level == "pole_match"
    curve_k = 55;
else
    curve_k = w_0+freq_gauss_width;
end

params = build_params( ...
    'k', curve_k, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', npolar, ...
    'n', n_quad, ...
    'p', p_corner, ...
    'p_edge', 2);

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

N = lp.RP_Curve.global.N;
fprintf("\nGeometry: N = %d boundary points (ppw=%d, %d bisections)\n", ...
    N, ppw, num_bisections);

if check_convergence
    check_curve_density_convergence(geom, params, ppw, num_bisections, w_0);
end

%% ------------------------------------------------------------
%  Refine every loaded pole against THIS curve's own BIE operator.
%
%  pols_wp was found on whatever curve data_gen_twinjet.m/
%  data_gen_quadjet.m actually built -- possibly with a different
%  rp2d_matlab code version (data_gen_twinjet.m's addpath points at
%  /home/msantana/Progamming/rp2d_matlab, not this repo's copy) even
%  when every resolution parameter is matched by hand. A direct test
%  (solving lp.bie_mat at the loaded pole vs. cont_rad away vs. a
%  generic frequency, i.e. exactly comp_poles.m's own scalarized
%  resolvent) showed real but insufficient near-pole amplification and
%  a peak not quite where the loaded coordinates say it is -- meaning
%  the loaded pole is only an approximate starting guess for THIS
%  curve's true resonance, not a precise one, and cont_rad=1e-5 residue
%  extraction needs it to be precise.
%
%  Rather than chase bit-identical curve reproduction, snap every
%  loaded pole onto this curve's actual near-singularity: run the same
%  secant refinement comp_poles.m uses (sec_ref/secant_method, 4
%  iterations, tol 1e-13) directly on this lp.bie_mat. This is exact
%  regardless of how the curve was built, so it sidesteps the whole
%  reproducibility question.
fprintf("\nRefining %d loaded poles against this curve's BIE operator...\n", ...
    numel(pols_wp));
rng(0);
u_res = rand(1, N, 'like', 1i);
v_res = rand(N, 1, 'like', 1i);
f_res = @(w) u_res*(lp.bie_mat(w)\v_res);
[pols_wp_refined, sec_errs] = sec_ref(f_res, pols_wp, 4, 1e-13);

% Guard against the secant method jumping off the genuine (decaying,
% Im<0) exterior scattering pole onto a spurious real-axis singularity
% -- pure single-layer ('sl') BIE formulations have extra singularities
% at real wavenumbers matching interior Dirichlet eigenvalues of the
% enclosed region, unrelated to physical resonances, and an unconstrained
% secant step starting this close to the real axis (high-Q poles have
% tiny |Im|) can converge to one of those instead. Observed concretely:
% refining the twinjet_plane1 target pole this way landed on Im~+0e-16,
% Q~4e8 -- nonsense. Reject any refined pole that flips unphysical
% (Im>=0) or jumps further than a small discretization correction should
% ever need, and fall back to the loaded value for that pole instead.
max_pole_shift = 1e-2;
for j = 1:numel(pols_wp)
    shift = abs(pols_wp_refined(j)-pols_wp(j));
    accept = imag(pols_wp_refined(j)) < 0 && shift <= max_pole_shift;
    if accept
        fprintf(['   pole %2d: k = %8.4f %+8.4fi -> %8.4f %+8.4fi   ' ...
            '(shift %.3e, secant resid %.3e) [accepted]\n'], j, ...
            real(pols_wp(j)), imag(pols_wp(j)), ...
            real(pols_wp_refined(j)), imag(pols_wp_refined(j)), ...
            shift, sec_errs(j));
    else
        fprintf(['   pole %2d: k = %8.4f %+8.4fi -> %8.4f %+8.4fi   ' ...
            '(shift %.3e, secant resid %.3e) [REJECTED -- ' ...
            'keeping loaded value]\n'], j, ...
            real(pols_wp(j)), imag(pols_wp(j)), ...
            real(pols_wp_refined(j)), imag(pols_wp_refined(j)), ...
            shift, sec_errs(j));
        pols_wp_refined(j) = pols_wp(j);
    end
end
pols_wp = pols_wp_refined;

% Recompute everything that was derived from the pre-refinement poles:
% Q ranking (target_pole_idx itself is kept fixed -- refinement should
% only shift each pole slightly, not reorder them), the target pole
% value, w_0, freq_gauss_width (distance to nearest *refined* neighbor),
% and the pulse timing derived from it.
Q = abs(real(pols_wp))./(2*abs(imag(pols_wp)));
target_pole = pols_wp(target_pole_idx);
w_0 = real(target_pole);

other_re = real(pols_wp);
other_re(target_pole_idx) = [];
if isempty(other_re)
    freq_gauss_width = 1.0;
else
    freq_gauss_width = min(abs(other_re-real(target_pole)));
end
if freq_gauss_width < min_freq_gauss_width
    fprintf(['Nearest neighboring pole is only %.4f away; clamping ' ...
        'freq_gauss_width to %.4f.\n'], freq_gauss_width, min_freq_gauss_width);
    freq_gauss_width = min_freq_gauss_width;
end

delta_t = 2*abs(log(gauss_picker_tol))/freq_gauss_width;
t_offset = startup_ratio*delta_t;

fprintf(['\nRefined target pole %d: k = %.6f %+.6fi  (Q ~ %.1f)\n' ...
    'Refined freq_gauss_width = %.4f, delta_t = %.3f, t_offset = %.3f\n\n'], ...
    target_pole_idx, real(target_pole), imag(target_pole), Q(target_pole_idx), ...
    freq_gauss_width, delta_t, t_offset);

%% ------------------------------------------------------------
%  Image grid: mask interior (solid) points exactly as
%  plot_scattering_field.m does, so the BIE representation is only
%  evaluated where it's valid (outside the scatterer).
% ------------------------------------------------------------

margin = 0.45;   % trimmed from 0.6 -- less empty padding, more grid
                 % resolution spent near the body for the same ng
padx = margin*max(max(xb)-min(xb), max(yb)-min(yb));
xlim_ = [min(xb)-padx, max(xb)+padx];
ylim_ = [min(yb)-padx, max(yb)+padx];

xv = linspace(xlim_(1), xlim_(2), ng);
yv = linspace(ylim_(1), ylim_(2), ng);
[xx, yy] = meshgrid(xv, yv);
xxv = xx(:); yyv = yy(:);

Solid = enclosed_polyshape(C, [], false);
has_solid = ~isempty(Solid);
if has_solid
    inside = isinterior(Solid, xxv, yyv);
    solid_area = area(Solid);
else
    inside = false(size(xxv));
    solid_area = 0;
end
ext_idx = find(~inside);

% Sanity check: enclosed_polyshape is a planar-graph face tracer and can
% silently mis-trace the main airframe loop into a degenerate,
% near-zero-area face at 3+-way piece junctions (seen concretely on the
% quadjet geometries, where the two-engines-per-side layout doubles the
% number of such junctions vs. twinjet's one-per-side). When that
% happens the whole airframe silhouette is left unmasked, so the BIE
% representation gets evaluated deep inside the solid -- producing large
% non-physical values that dominate the color scale and make the
% genuine (tiny, late-time) exterior field invisible. Compare the traced
% solid area against a crude lower bound (half the fuselage bounding-box
% area) so a repeat of this failure is loud instead of silently wrong.
% Only enforced when the fuselage is supposed to close (closed_front):
% an intentionally open-nose preset (closed_front=false) legitimately
% has little or no main-body solid area (just whatever nacelle loops
% happen to be closed), so the same floor would false-positive there.
fuselage_area_floor = 0.5*(max(xb)-min(xb))*(max(yb)-min(yb));
if gparams.closed_front && solid_area < 0.1*fuselage_area_floor
    warning('late_time_nacelle_field:DegenerateSolidMask', ...
        sprintf([ ...
        'enclosed_polyshape''s traced solid area (%.3g) is suspiciously ' ...
        'small compared to the fuselage bounding-box floor (%.3g), even ' ...
        'though closed_front is true for this preset. This usually ' ...
        'means the main airframe loop was mis-traced (see ' ...
        'enclosed_polyshape.m''s face-tracing) and grid points inside the ' ...
        'solid are NOT being masked out. Do not trust this figure until ' ...
        'this is fixed.'], solid_area, fuselage_area_floor));
end

fprintf("Image grid: %d x %d (%d exterior points of %d total)\n", ...
    ng, ng, numel(ext_idx), numel(xxv));
fprintf("Traced solid area: %.4g (bounding-box floor check: %.4g)\n", ...
    solid_area, fuselage_area_floor);

% t_late must be late enough for BOTH (a) the direct/geometric wavefront
% to have crossed the whole image domain, and (b) the target resonance
% to have visibly rung down -- the residue term's time dependence is
% exp(-|Im(target_pole)|*t) (see comp_exp_ints.m in TD_Solver_Code), so
% for a high-Q pole ((a) alone can leave the field looking like an
% undecayed CW/frequency-domain solution, since |Im(target_pole)| can be
% arbitrarily small. Require the resonant amplitude to have decayed to
% decay_tol relative to its value right after the transit buffer.
transit_buffer = hypot(diff(xlim_), diff(ylim_));
decay_tol = 1e-2;
ring_down_time = -log(decay_tol)/abs(imag(target_pole));
t_late = t_offset+delta_t+max(transit_buffer, ring_down_time);

fprintf("Transit buffer = %.3f, ring-down time (decay to %.0e) = %.3f\n", ...
    transit_buffer, decay_tol, ring_down_time);
fprintf("Resonant amplitude decay factor at t_late: exp(-|Im(pole)|*(t_late-t_offset-delta_t)) = %.3e\n", ...
    exp(-abs(imag(target_pole))*max(transit_buffer, ring_down_time)));
fprintf("t_late = %.3f\n\n", t_late);

%% ------------------------------------------------------------
%  Frequency envelope and problem setup
% ------------------------------------------------------------

[sigmas, wlims] = gauss_picker(w_0, freq_gauss_width, gauss_picker_tol);

tilt = tilt_deg*pi/180;
kappa = -[cos(tilt), sin(tilt)];   % incidence from +x, tilted off-normal

ps = problem_data( ...
    'N', N, ...
    'wlims', wlims, ...
    'numw', numw, ...
    'xs', xxv(ext_idx), ...
    'ys', yyv(ext_idx), ...
    'Tlims', [t_late, t_late], ...
    'numt', 1, ...
    'kappa', kappa, ...
    'inc_field', 'gaussian', ...
    'w_0', w_0, ...
    'sigmas', sigmas, ...
    't_0', t_offset, ...
    'is_open_curve', false, ...
    'is_rp_curve', true, ...
    'is_far_field', false);

%% ------------------------------------------------------------
%  Frequency-domain solves (parfor'd internally) + residues at the
%  already-known poles + vectorized near-field evaluation on the grid
% ------------------------------------------------------------

t_start = tic;

fprintf("Computing real-line LU factorizations (%d frequencies)...\n", ps.numw);
[L_rl, U_rl] = comp_bie_LU(ps.ws, lp);

fprintf("Computing pole-contour LU factorizations (%d poles x %d pts)...\n", ...
    numel(pols_wp), ps.n_res_w);
[L_pole, U_pole] = comp_pole_LU(ps, lp, pols_wp);

fprintf("Computing pole-contour surface densities...\n");
res_dens = comp_res_dens(ps, lp, pols_wp, L_pole, U_pole);

fprintf("Computing real-line surface densities...\n");
densities = comp_dens(ps, lp, L_rl, U_rl);

% Sanity check that the pole-contour solve is actually landing on a
% near-singularity of THIS discretized operator, not silently missing
% it. res_dens is solved exactly cont_rad away from the loaded pole
% location; densities is solved on the real line, which for a
% high-Q target pole already includes points only ~|Im(target_pole)|
% away (the closest a real frequency can get) -- i.e. densities already
% captures whatever near-resonance amplification exists at THIS
% discretization. If the pole-contour density isn't dramatically larger
% than that (orders of magnitude, since cont_rad=1e-5 is ~1e3-1e5x
% closer to the pole than the real line can get), the contour missed
% the true pole of this curve -- e.g. because pols_wp was found on a
% differently-discretized curve -- and comp_r_res_grid.m's residue
% extraction is returning near-zero numerical noise instead of the
% intended resonant signal.
pole_contour_mag = max(abs(res_dens(:, :, target_pole_idx)), [], 'all');
real_line_mag = max(abs(densities), [], 'all');
fprintf(['Pole-contour density max |phi| = %.4e vs real-line max |phi| ' ...
    '= %.4e (ratio %.3e)\n'], pole_contour_mag, real_line_mag, ...
    pole_contour_mag/real_line_mag);
if pole_contour_mag < 1e2*real_line_mag
    warning('late_time_nacelle_field:PoleContourNotSingular', ...
        sprintf([ ...
        'The pole-contour surface density is not dramatically larger ' ...
        'than the real-line density, even though the contour radius ' ...
        '(cont_rad=%.1e) is far closer to the target pole than any ' ...
        'real frequency can get. This suggests the cont_rad contour ' ...
        'is missing the true pole of THIS discretized curve (e.g. ' ...
        'pols_wp was found on a differently-discretized curve) -- the ' ...
        'residue/SEM resonant term is likely numerical noise, not ' ...
        'signal. Do not trust this figure until this is resolved.'], ...
        ps.cont_rad));
end

Bkslow = comp_bkslow(ps);

fprintf("Evaluating near-field residues on the image grid...\n");
r_res = comp_r_res_grid(ps, lp, pols_wp, res_dens);

fprintf("Evaluating near-field frequency solutions on the image grid...\n");
f_sols = comp_freq_sols_grid(ps, lp, densities, Bkslow, pols_wp, r_res);

fprintf("Synthesizing the late-time field...\n");
usol = comp_time_sol(ps, lp, f_sols, pols_wp, r_res);   % total field
                                                          % (incident field
                                                          % is added inside
                                                          % freq_to_time.m
                                                          % since is_far_field
                                                          % is false)

fprintf("Field computation took %.1f s\n\n", toc(t_start));

%% ------------------------------------------------------------
%  Save everything needed to re-render this figure later (different
%  color scheme, axes, line widths, ...) without recomputing anything:
%  the field grid AND the exact Curve/Solid geometry objects used, not
%  just the final image. See plot_late_time_field.m.
% ------------------------------------------------------------

U = nan(ng, ng);
U(ext_idx) = real(usol);
U = reshape(U, ng, ng);

Q_target = Q(target_pole_idx);

fname_base = sprintf("late_time_field_%s_pole%d_%s_ng%d", ...
    plane_name, target_pole_idx, resolution_level, ng);

mat_path = fullfile(save_dir, fname_base+".mat");
save(mat_path, "xv", "yv", "U", "C", "Solid", "t_late", "target_pole", ...
    "target_pole_idx", "Q_target", "pols_wp", "Q", "w_0", ...
    "freq_gauss_width", "sigmas", "tilt", "N", "ppw", "num_bisections", ...
    "ng", "numw", "plane_name", "plane_preset", "geometry_type", "resolution_level");
fprintf("Saved data to:   %s\n", mat_path);

img_path = fullfile(save_dir, fname_base+".png");
plot_late_time_field(mat_path, struct('save_path', img_path));
fprintf("Saved figure to: %s\n", img_path);

%% ------------------------------------------------------------
%  Local helper: curve-resolution convergence check (point 1)
% ------------------------------------------------------------

function check_curve_density_convergence(geom, params, ppw, num_bisections, w_0)
% Solves the surface density at the real frequency w_0 on the requested
% curve and on a curve refined by 2x, and reports the max difference
% (via densities_max_diff.m) so the resolution can be judged converged
% before spending time on the full pole/residue/image pipeline.

paramsC = params;
paramsC.k = w_0;

C1 = build_curve(geom, paramsC);
C1 = refine_curve_wavenumber(C1, 2*pi/paramsC.k, ppw);
for jj = 1:num_bisections
    C1 = bisect_all_patches(C1);
    C1 = build_global_indexing(C1);
    C1 = get_close_points(C1);
end

C2 = build_curve(geom, paramsC);
C2 = refine_curve_wavenumber(C2, 2*pi/paramsC.k, 2*ppw);
for jj = 1:(num_bisections+1)
    C2 = bisect_all_patches(C2);
    C2 = build_global_indexing(C2);
    C2 = get_close_points(C2);
end

phi1 = solve_planewave_scattering(C1, [1;0], paramsC);
phi2 = solve_planewave_scattering(C2, [1;0], paramsC);

dmax = densities_max_diff(C1, phi1, C2, phi2);

fprintf(['Convergence check at k=%.3f: N1=%d, N2=%d, ' ...
    'max density difference = %.3e\n\n'], ...
    w_0, C1.global.N, C2.global.N, dmax);

if dmax > 1e-6
    warning('late_time_nacelle_field:ConvergenceWarning', ...
        sprintf([ ...
        'Surface density changed by %.3e when the curve was refined ' ...
        '2x. Consider increasing ppw/num_bisections.'], dmax));
end

end
