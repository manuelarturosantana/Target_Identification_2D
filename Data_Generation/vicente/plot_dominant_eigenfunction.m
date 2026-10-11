% plot_dominant_eigenfunction.m
%
% Computes and plots the pole-residue eigenfunction for the dominant
% (closest-to-real-axis) pole of a plane's REGENERATED pole file (see
% regen_poles.m -- the original stale pole files under
% /scratch/vhojas/data_planes_2d were found on a much coarser curve
% and/or an since-changed geometry preset; this only ever reads the
% regenerated ones under .../regen/<plane_name>_regen.mat).
%
% Unlike plot_pole_eigenfunction.m, this doesn't need to walk a list of
% many stale secant-refined candidates looking for one that clears the
% resolved-residue ratio check -- the regenerated poles came from a
% global AAA-rational-approximation search on THIS exact curve (comp_poles,
% same method data_gen_twinjet.m/data_gen_quadjet.m use), so they're
% already genuine singularities of this operator. Still runs the same
% ratio check as a cheap sanity gate before trusting the result, and
% still labels the output RESOLVED/UNRESOLVED accordingly.
%
% ng controls image grid resolution -- pass NG env var (default 40, a
% coarse first-look grid) to sanity check before spending time on a
% finer grid.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

regen_dir = "/scratch/vhojas/data_planes_2d/regen";
save_dir  = fullfile(fileparts(mfilename('fullpath')), "figures");
if ~exist(save_dir, "dir")
    mkdir(save_dir);
end

configs = struct("name", {}, "geometry_type", {}, "preset", {});
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
end

ng = str2double(getenv("NG"));
if isnan(ng)
    ng = 40;   % coarse first-look default
end

config = configs(config_id);
plane_name    = config.name;
geometry_type = config.geometry_type;
plane_preset  = config.preset;
pole_data_file = fullfile(regen_dir, plane_name + "_regen.mat");

fprintf("Config %d/%d: %s (%s, preset %s), ng=%d\n", ...
    config_id, numel(configs), plane_name, geometry_type, plane_preset, ng);
fprintf("Pole data: %s\n", pole_data_file);

%% ------------------------------------------------------------
%  Parallel workers -- same profile that ran the pole regeneration
%  cleanly (no OOM).
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
%  Curve -- identical to regen_poles.m (the curve the pole file's poles
%  are genuine singularities of).
% ------------------------------------------------------------

tilt_deg = 5;
kappa = -[cosd(tilt_deg), sind(tilt_deg)];

ppw = 1; npolar = 160; n_quad = 20; p_corner = 3; p_edge = 2; delta_bc = 0.05;
curve_k = 55;

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
    'k', curve_k, 'formulation', 'sl', 'solver', 'direct', ...
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
%  Dominant pole -- closest to the real axis in the regenerated file.
% ------------------------------------------------------------

S = load(pole_data_file, "pols");
if isempty(S.pols)
    error("No poles in %s -- run regen_poles.m for this config first.", pole_data_file);
end
[~, target_idx] = min(abs(imag(S.pols(:))));
target_pole = S.pols(target_idx);
Q_target = abs(real(target_pole)) / (2*abs(imag(target_pole)));
fprintf("Dominant pole: k = %.6f%+.6fi  (Q ~ %.1f)\n", ...
    real(target_pole), imag(target_pole), Q_target);

%% ------------------------------------------------------------
%  Sanity-check ratio (same physically-grounded test used throughout
%  this codebase) -- and the actual residue density used for the field.
% ------------------------------------------------------------

cont_rad = 1e-5;
n_res_w  = 10;
ratio_threshold = 1e2;

circ = target_pole + cont_rad*exp(1i*2*pi*(1:n_res_w)/n_res_w);
target_res_dens = zeros(N, n_res_w);
for cind = 1:n_res_w
    rhs = -exp(1i*circ(cind)*(kappa*bndpts)).';
    target_res_dens(:, cind) = lp.bie_mat(circ(cind)) \ rhs;
end
rhs_real = -exp(1i*real(target_pole)*(kappa*bndpts)).';
dens_real = lp.bie_mat(real(target_pole)) \ rhs_real;
target_ratio = max(abs(target_res_dens(:))) / max(abs(dens_real(:)));
resolved = target_ratio > ratio_threshold;
fprintf("Contour/real-line ratio = %.3e   [%s]\n", target_ratio, ...
    ternary(resolved, "RESOLVED", "UNRESOLVED -- provisional"));

%% ------------------------------------------------------------
%  Image grid.
% ------------------------------------------------------------

margin = 0.15;
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
else
    inside = false(size(xxv));
end
ext_idx = find(~inside);

fprintf("Image grid: %d x %d (%d exterior points of %d total, margin=%.2f)\n", ...
    ng, ng, numel(ext_idx), numel(xxv), margin);

%% ------------------------------------------------------------
%  Evaluate the eigenfunction on the grid.
% ------------------------------------------------------------

ps = problem_data( ...
    'N', N, ...
    'wlims', [curve_k-1, curve_k+1], ...
    'w_0', curve_k, ...
    'xs', xxv(ext_idx), ...
    'ys', yyv(ext_idx), ...
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
%  Plot.
% ------------------------------------------------------------

clip_pct = 99;
finite_vals = abs(Efun(~isnan(Efun)));
cmax = prctile(finite_vals, clip_pct);
true_max = max(finite_vals);
fprintf("\n|field| true max = %.4e, %g-th percentile (color scale cap) = %.4e (ratio %.2f)\n", ...
    true_max, clip_pct, cmax, true_max/cmax);

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
colormap(parula(256));
if cmax > 0
    clim([-cmax, cmax]);
end
cb = colorbar;
cb.Label.String = 'Re(eigenfunction)';
box on;
title(sprintf('%s eigenfunction, pole k=%.4f%+.4fi (Q~%.0f) [%s, ng=%d]', ...
    strrep(plane_name, '_', '\_'), real(target_pole), imag(target_pole), Q_target, ...
    ternary(resolved, "RESOLVED", "UNRESOLVED"), ng));
hold off;

fname_base = sprintf("dominant_eigenfunction_%s_ng%d_%s", plane_name, ng, ...
    ternary(resolved, "resolved", "unresolved"));
img_path = fullfile(save_dir, fname_base+".png");
exportgraphics(fig, img_path, 'Resolution', 300);
fprintf("Saved figure to: %s\n", img_path);

mat_path = fullfile(save_dir, fname_base+".mat");
save(mat_path, "xv", "yv", "Efun", "C", "Solid", "target_pole", ...
    "target_idx", "Q_target", "target_ratio", "resolved", "ratio_threshold", ...
    "N", "ng", "margin", "clip_pct", "cmax", "true_max", ...
    "plane_name", "plane_preset", "geometry_type");
fprintf("Saved data to:   %s\n", mat_path);

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
