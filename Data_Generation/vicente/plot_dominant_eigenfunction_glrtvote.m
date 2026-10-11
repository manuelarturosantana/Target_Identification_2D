% plot_dominant_eigenfunction_glrtvote.m
%
% Same pipeline as plot_dominant_eigenfunction.m (regenerated poles,
% pole-residue eigenfunction on an image grid), but the target pole is
% chosen differently: instead of the pole closest to the real axis, this
% runs "regular GLRT" (the exact basis/scoring from
% Identification/build_glrt_model.m and classify_glrt_data.m) treating
% each regenerated candidate pole as its own 1-pole model, scores it
% against every one of that plane's 359 real (already-simulated, not
% synthesized here) receiver signals from
% /scratch/msantana/plane_data_360/<plane>.mat, and takes the pole that
% wins (highest score) at the most receivers -- majority vote, not an
% averaged score.
%
% late_time = 40.0 matches run_target_identification_averaged.m.
% svd_tol = 1e-12 matches the build_id_params.m default (unchanged --
% "regular GLRT" means reusing the existing defaults, not new ones).
%
% ng defaults to 128 (exploration resolution, per request) rather than
% plot_dominant_eigenfunction.m's finer default.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

regen_dir   = "/scratch/vhojas/data_planes_2d/regen";
signal_dir  = "/scratch/msantana/plane_data_360";
save_dir    = getenv("SAVE_DIR");
if isempty(save_dir)
    save_dir = fullfile(fileparts(mfilename('fullpath')), "figures");
end
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
    ng = 128;   % exploration resolution
end

config = configs(config_id);
plane_name    = config.name;
geometry_type = config.geometry_type;
plane_preset  = config.preset;
pole_data_file   = fullfile(regen_dir, plane_name + "_regen.mat");
signal_data_file = fullfile(signal_dir, plane_name + ".mat");

fprintf("Config %d/%d: %s (%s, preset %s), ng=%d\n", ...
    config_id, numel(configs), plane_name, geometry_type, plane_preset, ng);
fprintf("Candidate pole data: %s\n", pole_data_file);
fprintf("Receiver signal data: %s\n", signal_data_file);

%% ------------------------------------------------------------
%  Parallel workers.
%
%  Deliberately NOT pre-starting our own 'Processes' pool here (unlike
%  plot_dominant_eigenfunction.m, which this was adapted from) -- this
%  script never calls parfor itself. single_layer_extern.m's own
%  internal parfor (used by evaluate_helmholtz_potential inside
%  comp_r_res_grid) needs a 'Threads'-compatible pool (it broadcasts
%  containers.Map objects) and will start one on demand. Pre-starting a
%  separate 'Processes' pool left both pools alive at once at ng=1024's
%  memory footprint, which crashed with a fatal "malloc(): unaligned
%  tcache chunk detected" (glibc heap corruption) on quadjet_plane2 --
%  did not reproduce at ng=512, so this only bit at the larger grid's
%  memory pressure.
% ------------------------------------------------------------

%% ------------------------------------------------------------
%  Regular GLRT, per candidate pole, majority vote across receivers.
% ------------------------------------------------------------

late_time = 40.0;   % matches run_target_identification_averaged.m
svd_tol   = 1e-12;  % matches Identification/build_id_params.m default

Sp = load(pole_data_file, "pols");
candidate_poles = Sp.pols(:);
npoles = numel(candidate_poles);
if npoles == 0
    error("No candidate poles in %s.", pole_data_file);
end

Ssig = load(signal_data_file, "uff_all", "ts");
ts = Ssig.ts(:);
uff_all = Ssig.uff_all;
nrec = size(uff_all, 1);

late_idx = find(ts >= late_time, 1, "first");
if isempty(late_idx)
    error("late_time = %.1f exceeds the final time sample in %s.", late_time, signal_data_file);
end
tt = ts(late_idx:end);

fprintf("\n%d candidate poles, %d receivers, late_idx=%d (t=%.2f), tt length=%d\n", ...
    npoles, nrec, late_idx, ts(late_idx), numel(tt));

fprintf("Building single-pole GLRT bases...\n");
Ur_all = cell(npoles, 1);
for ip = 1:npoles
    p = candidate_poles(ip);
    decay = imag(p); osc = real(p);
    Zc = exp(tt * decay) .* cos(tt * osc);
    Zs = exp(tt * decay) .* sin(tt * osc);
    Z = [Zc, Zs];
    [U, S, ~] = svd(Z, "econ");
    s = diag(S);
    if isempty(s)
        keep = false(size(s));
    else
        keep = s > svd_tol * max(s);
    end
    Ur_all{ip} = U(:, keep);
end

fprintf("Scoring %d receivers against %d single-pole models...\n", nrec, npoles);
winner = zeros(nrec, 1);
for irec = 1:nrec
    f = real(uff_all(irec, late_idx:end)).';
    scores = zeros(1, npoles);
    for ip = 1:npoles
        Ur = Ur_all{ip};
        scores(ip) = norm(Ur' * f, 2);
    end
    [~, winner(irec)] = max(scores);
end

vote_counts = accumarray(winner, 1, [npoles, 1]);
[max_votes, majority_idx] = max(vote_counts);
target_pole = candidate_poles(majority_idx);
target_idx = majority_idx;

fprintf("\nVote tally (pole: votes/%d receivers):\n", nrec);
[sorted_votes, order] = sort(vote_counts, "descend");
for ii = 1:min(10, npoles)
    jj = order(ii);
    fprintf("  k = %.6f%+.6fi : %d votes (%.1f%%)\n", ...
        real(candidate_poles(jj)), imag(candidate_poles(jj)), ...
        sorted_votes(ii), 100 * sorted_votes(ii) / nrec);
end
fprintf("\nMajority-vote pole: k = %.6f%+.6fi  (%d/%d votes, %.1f%%)\n", ...
    real(target_pole), imag(target_pole), max_votes, nrec, 100 * max_votes / nrec);

%% ------------------------------------------------------------
%  Curve -- identical to regen_poles.m / plot_dominant_eigenfunction.m.
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

Q_target = abs(real(target_pole)) / (2*abs(imag(target_pole)));
fprintf("Target pole: k = %.6f%+.6fi  (Q ~ %.1f)\n", ...
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
title(sprintf('%s eigenfunction (GLRT vote), pole k=%.4f%+.4fi (Q~%.0f) [%s, ng=%d, %d/%d votes]', ...
    strrep(plane_name, '_', '\_'), real(target_pole), imag(target_pole), Q_target, ...
    ternary(resolved, "RESOLVED", "UNRESOLVED"), ng, max_votes, nrec));
hold off;

fname_base = sprintf("dominant_eigenfunction_glrtvote_%s_ng%d_%s", plane_name, ng, ...
    ternary(resolved, "resolved", "unresolved"));
img_path = fullfile(save_dir, fname_base+".png");
exportgraphics(fig, img_path, 'Resolution', 300);
fprintf("Saved figure to: %s\n", img_path);

mat_path = fullfile(save_dir, fname_base+".mat");
save(mat_path, "xv", "yv", "Efun", "C", "Solid", "target_pole", ...
    "target_idx", "Q_target", "target_ratio", "resolved", "ratio_threshold", ...
    "N", "ng", "margin", "clip_pct", "cmax", "true_max", ...
    "plane_name", "plane_preset", "geometry_type", ...
    "candidate_poles", "vote_counts", "winner", "max_votes", "nrec", ...
    "late_time", "svd_tol");
fprintf("Saved data to:   %s\n", mat_path);

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
