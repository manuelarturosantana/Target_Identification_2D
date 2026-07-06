clear;

clf
tiledlayout(2,5);
fig = figure(1);

for obstacle_idx = 1:10
obstacle_idx

% Plots
curve = build_plane(obstacle_idx);
nexttile

opts = struct(); opts.show_nodes= false; opts.plot_color='black'; opts.line_width=2;
opts.no_gaps = true;
plot_curve2(curve, opts, obstacle_idx);
grid on

end
%%%%%%%%%%%%  Helper Functions %%%%%%%%%%%%%%%%%%

function Curve = build_plane(idx)
    params = build_params( ...
    'k', 10, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 160, ...
    'n', 20, ...
    'p', 3, ...
    'p_edge', 2);

  
    % Look at results.obstacles to see first its the quadjets, then its the
    % twin jets
    if idx <= 5
        preset_name = sprintf("plane%d", idx);
        gparams     = build_quadjet_params('preset', preset_name);
        geometry    = quadjet_plane_geometry(gparams);
    else
        preset_name = sprintf("plane%d", idx - 5);
        gparams     = build_twinjet_params('preset', preset_name);
        geometry    = twinjet_plane_geometry(gparams);
   
    end

    Curve = build_curve(geometry, params);
    Curve = refine_curve_wavenumber(Curve, 2*pi/params.k, 1);

end

function plot_curve2(Curve, opts, idx)
%PLOT_CURVE Plot patches and discretization nodes of a Curve struct.
% TODO: Bad
%
% plot_curve(Curve)
% plot_curve(Curve, opts)
%
% Required fields:
%   Curve.patches{j} with fields x(t), y(t), tj, xj, yj
%
% Optional opts fields:
%   opts.nplot        number of plot points per patch (default 200)
%   opts.show_nodes   true/false (default true)
%   opts.show_labels  true/false (default false)
%   opts.node_size    marker size (default 18)
%   opts.line_width   line width (default 1.25)
%   opts.plot_color   What color to plot the plane (default all patches
%                     different color.)
%   opts.no_gaps      Plot patches from [-1,1] instead of discretization
%                     points. 

if nargin < 2
    opts = struct();
end

if ~isfield(opts, "nplot"), opts.nplot = 200; end
if ~isfield(opts, "show_nodes"), opts.show_nodes = true; end
if ~isfield(opts, "show_labels"), opts.show_labels = false; end
if ~isfield(opts, "node_size"), opts.node_size = 18; end
if ~isfield(opts, "line_width"), opts.line_width = 1.25; end
if ~isfield(opts, "no_gaps"), opts.no_gaps = false; end

patches = Curve.patches;

hold_state = ishold;
hold on;

for j = 1:numel(patches)
    P = patches{j};

    % pick plotting parameter range; assume tj in [-1, 1]
    tmin = min(P.tj);
    tmax = max(P.tj);
    if (opts.no_gaps)
        t = linspace(-1, 1, opts.nplot);
        if (opts.show_nodes)
            warning("plot_curve: Show nodes and no_gaps both set to true in opts. Setting no_gaps=true does not use discretization nodes.")
        end
    else
        t = linspace(tmin, tmax, opts.nplot);
    end
   

    % patch curve
    x = P.x(t);
    y = P.y(t);

    if ~isfield(opts, "plot_color")
        plot(x, y, "-", "LineWidth", opts.line_width);
    else
        plot(x, y, "-", "LineWidth", opts.line_width,'Color',opts.plot_color);
    end

    % discretization nodes
    if opts.show_nodes
        plot(P.xj, P.yj, ".", "MarkerSize", opts.node_size);
    end

    % label near the patch midpoint
    if opts.show_labels
        tm = 0.5 * (tmin + tmax);
        xm = P.x(tm);
        ym = P.y(tm);
        text(xm, ym, sprintf("%d", j), ...
            "VerticalAlignment", "bottom", ...
            "HorizontalAlignment", "center");
    end
end

if (idx <= 5) 
    title("Quadjet " + num2str(idx));
end
   
if (idx > 5)
    title("Twinjet " + num2str(idx - 5));
end



axis equal;
grid on;

if ~hold_state
    hold off;
end

end