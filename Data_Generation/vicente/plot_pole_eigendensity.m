function fig = plot_pole_eigendensity(mat_path, opts)
%PLOT_POLE_EIGENDENSITY  Render |eigenfunction|^2 from a saved pole-eigenfunction .mat.
%
% Loads a .mat file saved by plot_pole_eigenfunction.m (fields xv, yv,
% Efun, C, Solid, target_pole, ...) and plots the eigendensity
% |Efun|.^2 on a sequential colormap, instead of Re(Efun) on the
% diverging one plot_pole_eigenfunction.m itself draws. Sign is
% meaningless for a density, so this never needs the diverging scale.
% Re-renders only -- no re-solve.
%
% opts fields (all optional):
%   colormap    colormap matrix, default parula(256)
%   log_scale   if true, plot log10(density) instead of density itself
%               (default false)
%   log_floor   log10 floor relative to the peak, e.g. -8 clips anything
%               more than 1e8 below it (default -8)
%   clip_pct    percentile of the density used as the color limit cap,
%               matching plot_pole_eigenfunction.m's near-boundary
%               quadrature-blowup rationale (default 99, ignored when
%               log_scale is true)
%   line_width  patch outline / solid edge line width, default 1
%   edge_color  solid outline color, default [0.15 0.15 0.15]
%   face_color  solid fill color, default [0.78 0.78 0.78]
%   title       override the title string, default built from saved metadata
%   save_path   if set, exportgraphics(fig, save_path, 'Resolution', 300)

if nargin < 2
    opts = struct();
end
if ~isfield(opts, 'colormap'), opts.colormap = parula(256); end
if ~isfield(opts, 'log_scale'), opts.log_scale = false; end
if ~isfield(opts, 'log_floor'), opts.log_floor = -8; end
if ~isfield(opts, 'clip_pct'), opts.clip_pct = 99; end
if ~isfield(opts, 'line_width'), opts.line_width = 1; end
if ~isfield(opts, 'edge_color'), opts.edge_color = [0.15, 0.15, 0.15]; end
if ~isfield(opts, 'face_color'), opts.face_color = [0.78, 0.78, 0.78]; end
if ~isfield(opts, 'title'), opts.title = ''; end
if ~isfield(opts, 'save_path'), opts.save_path = ''; end

S = load(mat_path);
density = abs(S.Efun).^2;

if opts.log_scale
    Dplot = log10(density);
else
    Dplot = density;
end

fig = figure('Color', 'w');
imagesc(S.xv, S.yv, Dplot);
set(gca, 'YDir', 'normal');
hold on;

if ~isempty(S.Solid)
    plot(S.Solid, ...
        'FaceColor', opts.face_color, 'FaceAlpha', 1, ...
        'EdgeColor', opts.edge_color, 'LineWidth', opts.line_width);
end

for j = 1:numel(S.C.patches)
    Patch = S.C.patches{j};
    t = linspace(-1, 1, 100);
    plot(Patch.x(t), Patch.y(t), 'k-', 'LineWidth', opts.line_width);
end

axis equal tight;
colormap(opts.colormap);
if opts.log_scale
    pk = max(Dplot(:), [], 'omitnan');
    clim([pk+opts.log_floor, pk]);
    cb = colorbar;
    cb.Label.String = 'log_{10}|eigenfunction|^2';
else
    finite_vals = density(~isnan(density));
    cmax = prctile(finite_vals, opts.clip_pct);
    if cmax > 0
        clim([0, cmax]);
    end
    cb = colorbar;
    cb.Label.String = '|eigenfunction|^2';
end
box on;

if isempty(opts.title)
    title(sprintf('%s eigendensity, pole k=%.4f%+.4fi (Q~%.0f) [%s]', ...
        strrep(S.plane_name, '_', '\_'), real(S.target_pole), imag(S.target_pole), ...
        S.Q_target, ternary(S.resolved, "RESOLVED", "UNRESOLVED")));
else
    title(opts.title);
end
hold off;

if ~isempty(opts.save_path)
    exportgraphics(fig, opts.save_path, 'Resolution', 300);
end

end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
