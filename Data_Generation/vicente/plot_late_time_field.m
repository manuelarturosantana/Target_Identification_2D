function fig = plot_late_time_field(mat_path, opts)
%PLOT_LATE_TIME_FIELD  Render (or re-render) a late-time field snapshot.
%
% Loads a .mat file saved by late_time_nacelle_field.m and draws the
% figure. Everything needed to redraw it -- the field grid, the exact
% Curve/Solid geometry objects used at compute time, and the target-pole
% metadata for the title -- was saved alongside the field itself, so
% this never needs to re-run any solve; it only re-renders. Use this to
% try different color schemes, line widths, or axis choices without
% resubmitting a SLURM job.
%
% opts fields (all optional):
%   colormap     colormap matrix, default parula(256)
%   line_width   patch outline / solid edge line width, default 1
%   edge_color   solid outline color, default [0.15 0.15 0.15]
%   face_color   solid fill color, default [0.78 0.78 0.78]
%   cmax         override the color limit, default max(abs(U))
%   log_scale    if true, plot log10(abs(U)) instead of Re(U) -- a
%                sequential (not diverging) scale, since the sign of
%                Re(U) is meaningless on a log axis. Useful as a first
%                diagnostic: with the near-boundary artifacts still in
%                the picture, this shows whether the field away from
%                them falls off like a sensible decaying field or is
%                just noise (default false)
%   log_floor    log10 floor relative to the max, e.g. -6 clips anything
%                more than 1e6 below the peak (default -8)
%   title        override the title string, default built from saved metadata
%   save_path    if set, exportgraphics(fig, save_path, 'Resolution', 300)

if nargin < 2
    opts = struct();
end
if ~isfield(opts, 'log_scale'), opts.log_scale = false; end
if ~isfield(opts, 'colormap')
    if opts.log_scale
        opts.colormap = parula(256);
    else
        opts.colormap = parula(256);
    end
end
if ~isfield(opts, 'line_width'), opts.line_width = 1; end
if ~isfield(opts, 'edge_color'), opts.edge_color = [0.15, 0.15, 0.15]; end
if ~isfield(opts, 'face_color'), opts.face_color = [0.78, 0.78, 0.78]; end
if ~isfield(opts, 'cmax'), opts.cmax = []; end
if ~isfield(opts, 'log_floor'), opts.log_floor = -8; end
if ~isfield(opts, 'title'), opts.title = ''; end
if ~isfield(opts, 'save_path'), opts.save_path = ''; end

S = load(mat_path);

if opts.log_scale
    Uplot = log10(abs(S.U));
else
    Uplot = S.U;
end

fig = figure('Color', 'w');
imagesc(S.xv, S.yv, Uplot);
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
    pk = max(Uplot(:), [], 'omitnan');
    if isempty(opts.cmax)
        clim([pk+opts.log_floor, pk]);
    else
        clim([opts.cmax(1), pk]);
    end
    cb = colorbar;
    cb.Label.String = 'log_{10}|u_{tot}|';
else
    if isempty(opts.cmax)
        cmax = max(abs(S.U(:)), [], 'omitnan');
    else
        cmax = opts.cmax;
    end
    if cmax > 0
        clim([-cmax, cmax]);
    end
    cb = colorbar;
    cb.Label.String = 'Re(u_{tot})';
end
box on;

if isempty(opts.title)
    if isfield(S, 'plane_name')
        name_str = S.plane_name;   % current field name
    else
        name_str = S.plane_preset; % legacy field name (pre-batch scripts)
    end
    title(sprintf('%s, t = %.1f (pole k=%.3f%+.4fi, Q~%.0f)', ...
        name_str, S.t_late, real(S.target_pole), imag(S.target_pole), ...
        S.Q_target));
else
    title(opts.title);
end
hold off;

if ~isempty(opts.save_path)
    exportgraphics(fig, opts.save_path, 'Resolution', 300);
end

end
