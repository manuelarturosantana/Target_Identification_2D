% plot_pole_eigenfunction.m
%
% Renders the SEM eigenfunction (pole-residue surface density, evaluated
% off the boundary) of a single trapped-mode pole, instead of synthesizing
% a late-time pulse. This sidesteps the whole time-domain pipeline
% (Gaussian bandwidth/duration coupling, t_late timing, etc. -- see
% late_time_nacelle_field.m) and just asks: what does this resonance's
% spatial mode shape actually look like?
%
% Target pole selection: among the loaded poles, walk outward from the
% real axis (smallest |Im(pole)| first, i.e. highest-Q first for poles
% with similar Re) and take the first one whose pole-contour residue is
% actually resolved on THIS curve -- not silent numerical noise standing
% in for a missed singularity. "Resolved" is checked the same way
% late_time_nacelle_field.m already sanity-checks this (comparing the
% cont_rad-contour density against the closest a real frequency can get),
% just applied per-candidate at the boundary-density level (cheap: no
% grid evaluation needed) instead of only for a single pre-chosen pole.
%
% Curve discretization matches the pole-finding scripts exactly
% (Data_Generation/data_gen_twinjet.m / data_gen_quadjet.m: curve_k=55,
% ppw=1, npolar=160, n=20, p=3, p_edge=2, delta=0.05, and bisected 2x for
% twinjet / 1x for quadjet) -- not late_time_nacelle_field.m's
% "resolution_level"-gated approximation of the same thing.
%
% The image grid uses a smaller margin than late_time_nacelle_field.m's
% (points closer to the boundary show more of the mode's near-field
% structure, which is the point here), and the color scale is clipped to
% a high percentile of |field| rather than the hard max: a handful of
% pixels immediately against the curve are dominated by near-singular
% quadrature error (the free-space Green's function blows up as an
% evaluation point approaches the source curve), and on a linear scale
% keyed to the true max those pixels alone wash out everything else.
% Clipping the *color scale* rather than the *grid extent* is deliberate:
% the nacelle ducts here are only ~engine_width=0.01-0.03 wide (see
% build_twinjet_params.m), so a spatial exclusion buffer wide enough to
% hide the boundary-quadrature artifacts would also hide the inside of
% the duct itself -- exactly the region this plot exists to show.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

%% ------------------------------------------------------------
%  Target geometry / pole data (same config table as
%  late_time_nacelle_field.m, same fallback default of config 6 ==
%  twinjet_plane1 -- the plane the noise-floor figure came from).
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

config = configs(config_id);
plane_name     = config.name;
geometry_type  = config.geometry_type;
plane_preset   = config.preset;
pole_data_file = config.pole_data_file;

fprintf("Config %d/%d: %s (%s, preset %s)\n", ...
    config_id, numel(configs), plane_name, geometry_type, plane_preset);
fprintf("Pole data: %s\n", pole_data_file);

tilt_deg = 5;   % same rationale as late_time_nacelle_field.m: an
                % on-axis illumination direction would leave a
                % degenerate symmetric/antisymmetric nacelle-pair mode
                % unexcited by symmetry, which would make its residue
                % vanish for a reason having nothing to do with whether
                % it's "relevant" -- tilting avoids that false negative.
kappa = -[cosd(tilt_deg), sind(tilt_deg)];

ratio_threshold = 1e2;   % same threshold late_time_nacelle_field.m uses
                          % for its own pole-contour-vs-real-line sanity
                          % check (:566-583) -- how much bigger the
                          % cont_rad-contour density must be than the
                          % closest a real frequency can get, to call the
                          % pole "actually resolved" rather than noise.

save_dir = fullfile(fileparts(mfilename('fullpath')), "figures");
if ~exist(save_dir, "dir")
    mkdir(save_dir);
end

%% ------------------------------------------------------------
%  Parallel workers (copied from late_time_nacelle_field.m -- same
%  per-run JobStorageLocation reasoning to avoid SLURM array tasks
%  colliding on the shared default pool directory).
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
%  Geometry + curve, built exactly as the pole-finding scripts build it
%  (Data_Generation/data_gen_twinjet.m:56-64,111-119 /
%  data_gen_quadjet.m, the curve pols_wp was actually found on) -- not
%  late_time_nacelle_field.m's separately-tuned/broken "pole_match" path.
% ------------------------------------------------------------

ppw            = 1;
npolar         = 160;
n_quad         = 20;
p_corner       = 3;
p_edge         = 2;
delta_bc       = 0.05;

pf_w0    = 50;   % pole-finding scripts' own w_0/freq_gauss_width (only
pf_width = 5;    % used here to reproduce curve_k -- see build_params
                 % below; not an incident-pulse parameter in this script
curve_k = pf_w0 + pf_width;

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

params = build_params( ...
    'k', curve_k, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', delta_bc, ...
    'npolar', npolar, ...
    'n', n_quad, ...
    'p', p_corner, ...
    'p_edge', p_edge);

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
fprintf("\nGeometry: N = %d boundary points (ppw=%d, %d bisections, curve_k=%.1f)\n", ...
    N, ppw, num_bisections, curve_k);

%% ------------------------------------------------------------
%  Load poles and refine each against THIS curve's own BIE operator
%  (identical secant-refinement + unphysical-jump guard to
%  late_time_nacelle_field.m:334-394 -- see comments there for why the
%  loaded pole locations are only approximate starting guesses).
% ------------------------------------------------------------

S = load(pole_data_file, "pols");
pols_wp = S.pols(:);
fprintf("Loaded %d poles from %s\n", numel(pols_wp), pole_data_file);

rng(0);
u_res = rand(1, N, 'like', 1i);
v_res = rand(N, 1, 'like', 1i);
f_res = @(w) u_res*(lp.bie_mat(w)\v_res);
[pols_wp_refined, sec_errs] = sec_ref(f_res, pols_wp, 4, 1e-13);

max_pole_shift = 1e-2;
for j = 1:numel(pols_wp)
    shift = abs(pols_wp_refined(j)-pols_wp(j));
    accept = imag(pols_wp_refined(j)) < 0 && shift <= max_pole_shift;
    if accept
        pols_wp(j) = pols_wp_refined(j);
    end
    fprintf("   pole %2d: k = %8.4f %+8.4fi   (secant resid %.3e) [%s]\n", ...
        j, real(pols_wp(j)), imag(pols_wp(j)), sec_errs(j), ...
        ternary(accept, "refined", "kept loaded"));
end

Q = abs(real(pols_wp))./(2*abs(imag(pols_wp)));

%% ------------------------------------------------------------
%  Pick the pole closest to the real axis whose residue is actually
%  resolved on this discretization -- not just the highest raw |Im|
%  rank, since a pole the contour missed would otherwise silently win by
%  "looking" high-Q while its residue is discretization noise.
% ------------------------------------------------------------

cont_rad = 1e-5;   % gen_params.m's own documented default for this
n_res_w  = 10;      % same as gen_params.m's documented default

[~, order] = sort(abs(imag(pols_wp)), 'ascend');

target_idx = [];
target_res_dens = [];
target_ratio = -Inf;
best_idx = [];
best_res_dens = [];
best_ratio = -Inf;
fprintf("\nScreening poles by closeness to the real axis, nearest first:\n");
for oi = 1:numel(order)
    j = order(oi);
    p = pols_wp(j);

    circ = p + cont_rad*exp(1i*2*pi*(1:n_res_w)/n_res_w);
    res_dens_j = zeros(N, n_res_w);
    for cind = 1:n_res_w
        rhs = -exp(1i*circ(cind)*(kappa*bndpts)).';
        res_dens_j(:, cind) = lp.bie_mat(circ(cind)) \ rhs;
    end

    rhs_real = -exp(1i*real(p)*(kappa*bndpts)).';
    dens_real = lp.bie_mat(real(p)) \ rhs_real;

    ratio = max(abs(res_dens_j(:))) / max(abs(dens_real(:)));
    fprintf(['   pole %2d: k = %8.4f %+8.4fi   Q ~ %8.2f   ' ...
        'contour/real-line ratio = %.3e   %s\n'], ...
        j, real(p), imag(p), Q(j), ratio, ...
        ternary(ratio > ratio_threshold, "[RESOLVED -- using this one]", "[not resolved, skipping]"));

    if ratio > best_ratio
        best_idx = j; best_res_dens = res_dens_j; best_ratio = ratio;
    end

    if ratio > ratio_threshold
        target_idx = j;
        target_res_dens = res_dens_j;
        target_ratio = ratio;
        break
    end
end

resolved = ~isempty(target_idx);
if ~resolved
    % None cleared the threshold -- fall back to the best available
    % candidate rather than refusing to plot anything. A ratio well
    % below ratio_threshold means this is NOT a validated eigenfunction
    % (see the RCOND warnings from the secant refinement above: the BIE
    % matrix is landing on near machine-precision-singular at some of
    % these frequencies, which can make the near-pole solve itself
    % numerically meaningless rather than genuinely showing "no
    % resonance here") -- the plot below is labeled accordingly.
    warning('plot_pole_eigenfunction:NoPoleResolved', sprintf([ ...
        'None of the %d loaded poles cleared the ratio_threshold=%.1e ' ...
        'resolved-residue check (best ratio = %.3e, pole %d). Plotting ' ...
        'that best candidate anyway, clearly labeled as unresolved -- ' ...
        'treat this figure as provisional, not a validated eigenfunction.'], ...
        numel(pols_wp), ratio_threshold, best_ratio, best_idx));
    target_idx = best_idx;
    target_res_dens = best_res_dens;
    target_ratio = best_ratio;
end

target_pole = pols_wp(target_idx);
Q_target = Q(target_idx);
fprintf("\nSelected pole %d: k = %.6f %+.6fi  (Q ~ %.1f)  ratio = %.3e  [%s]\n", ...
    target_idx, real(target_pole), imag(target_pole), Q_target, target_ratio, ...
    ternary(resolved, "RESOLVED", "UNRESOLVED -- provisional"));

%% ------------------------------------------------------------
%  Image grid -- smaller margin than late_time_nacelle_field.m (0.15 vs
%  0.45) so the plot is framed tightly around the airframe/nacelles
%  rather than mostly empty far field.
% ------------------------------------------------------------

margin = 0.15;
ng     = 160;
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

fuselage_area_floor = 0.5*(max(xb)-min(xb))*(max(yb)-min(yb));
if gparams.closed_front && solid_area < 0.1*fuselage_area_floor
    warning('plot_pole_eigenfunction:DegenerateSolidMask', ...
        sprintf([ ...
        'enclosed_polyshape''s traced solid area (%.3g) is suspiciously ' ...
        'small vs. the fuselage bounding-box floor (%.3g) even though ' ...
        'closed_front is true -- see late_time_nacelle_field.m''s ' ...
        'identical check for what this usually means. Do not trust this ' ...
        'figure until this is resolved.'], solid_area, fuselage_area_floor));
end

fprintf("\nImage grid: %d x %d (%d exterior points of %d total, margin=%.2f)\n", ...
    ng, ng, numel(ext_idx), numel(xxv), margin);

%% ------------------------------------------------------------
%  Evaluate the eigenfunction (pole-residue field) on the grid.
% ------------------------------------------------------------

ps = problem_data( ...
    'N', N, ...
    'wlims', [curve_k-1, curve_k+1], ...   % dummy -- only gen_params.m's
    'w_0', curve_k, ...                    % gaussian-in-band sanity check
    'xs', xxv(ext_idx), ...                % needs these to be consistent;
    'ys', yyv(ext_idx), ...                % no pulse is actually used here
    'kappa', kappa, ...
    'is_rp_curve', true, ...
    'is_open_curve', false, ...
    'n_res_w', n_res_w, ...
    'cont_rad', cont_rad);

fprintf("\nEvaluating the pole-residue field on the image grid...\n");
r_res = comp_r_res_grid(ps, lp, target_pole, reshape(target_res_dens, N, n_res_w, 1));

Efun = nan(ng, ng);
Efun(ext_idx) = r_res;
Efun = reshape(Efun, ng, ng);

%% ------------------------------------------------------------
%  Plot. Color scale is clipped to a high percentile of |field| rather
%  than the true max -- see the header comment for why a handful of
%  near-boundary quadrature-blowup pixels otherwise dominate a
%  max()-keyed linear scale, and why a spatial exclusion buffer isn't
%  used instead.
% ------------------------------------------------------------

clip_pct = 99;
finite_vals = abs(Efun(~isnan(Efun)));
cmax = prctile(finite_vals, clip_pct);
true_max = max(finite_vals);
fprintf(['\n|field| true max = %.4e, %g-th percentile (color scale cap) ' ...
    '= %.4e (ratio %.2f)\n'], true_max, clip_pct, cmax, true_max/cmax);

fig = figure('Color', 'w');
imagesc(xv, yv, real(Efun));
set(gca, 'YDir', 'normal');
hold on;

if has_solid
    plot(Solid, 'FaceColor', [0.78 0.78 0.78], 'FaceAlpha', 1, ...
        'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 1);
end
for j = 1:numel(C.patches)
    Patch = C.patches{j};
    t = linspace(-1, 1, 100);
    plot(Patch.x(t), Patch.y(t), 'k-', 'LineWidth', 1);
end

axis equal tight;
colormap(diverging_colormap(256));
if cmax > 0
    clim([-cmax, cmax]);
end
cb = colorbar;
cb.Label.String = 'Re(eigenfunction)';
box on;
title(sprintf('%s eigenfunction, pole k=%.4f%+.4fi (Q~%.0f) [%s, ratio=%.2e]', ...
    strrep(plane_name, '_', '\_'), real(target_pole), imag(target_pole), Q_target, ...
    ternary(resolved, "RESOLVED", "UNRESOLVED"), target_ratio));
hold off;

fname_base = sprintf("pole_eigenfunction_%s_pole%d_ng%d_%s", plane_name, target_idx, ng, ...
    ternary(resolved, "resolved", "unresolved"));
img_path = fullfile(save_dir, fname_base+".png");
exportgraphics(fig, img_path, 'Resolution', 300);
fprintf("Saved figure to: %s\n", img_path);

mat_path = fullfile(save_dir, fname_base+".mat");
save(mat_path, "xv", "yv", "Efun", "C", "Solid", "target_pole", ...
    "target_idx", "Q_target", "target_ratio", "resolved", "ratio_threshold", ...
    "pols_wp", "Q", "kappa", "N", "ppw", ...
    "num_bisections", "ng", "margin", "clip_pct", "cmax", "true_max", ...
    "plane_name", "plane_preset", "geometry_type");
fprintf("Saved data to:   %s\n", mat_path);

%% ------------------------------------------------------------
%  Local helpers
%% ------------------------------------------------------------

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
