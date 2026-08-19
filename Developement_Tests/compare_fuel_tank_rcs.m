% svd_null_density_rcs.m
%
% Compute the SVD nullspace-vector density (right singular vector for the
% smallest singular value of the BIE matrix) for two twinjet geometries
% (with and without fuel tank), each at its own wavenumber k, then plot
% the resulting |far field| vs. observation angle ("RCS-like" curve).

clear;

% paths
addpath(genpath('/home/msantana/Progamming/rp2d_matlab/src'))
addpath(genpath('/home/msantana/Progamming/Target_Identification_2D'))
addpath(genpath('/home/msantana/Progamming/MATLAB_PACKAGES'))

% ------------------------------------------------------------
% Plot mode flag: true -> dB (20*log10|.|), false -> raw |far field|
% ------------------------------------------------------------
use_db = true;

% ------------------------------------------------------------
% User-specified k values (different for each geometry)
% ------------------------------------------------------------
load("/scratch/msantana/plane_data_360/twinjet_plane1.mat");
k_no_fuel_tank = pols; 
load("../Data_Generation/plane_data/twinjet1_with_fuel_tank.mat");
k_fuel_tank    = pols;  

% Curve is refined for the largest k this solver will ever be asked to
% evaluate (per RP2LP's own convention: params.k = "highest frequency").
k_refine = 55;

% ------------------------------------------------------------
% Solver / refinement params (same as file 1 / file 2)
% ------------------------------------------------------------
params = build_params( ...
    'k', k_refine, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 160, ...
    'n', 20, ...
    'p', 3, ...
    'p_edge', 2);

% ------------------------------------------------------------
% Build geometries
% ------------------------------------------------------------

% --- Twinjet 1, no fuel tank (file 1, preset "plane1") ---
gparams_no_ft = build_twinjet_params('preset', 'plane1');
geom_no_ft    = twinjet_plane_geometry(gparams_no_ft);

% --- Twinjet 1, with fuel tank (file 2) ---
gparams_ft = build_twinjet_fuel_tank_params('preset', 'plane1', ...
    'fuel_tank_loc', 0.45, ...
    'fuel_tank_length', 0.15, ...
    'fuel_tank_width', 0.006, ...
    'fuel_tank_shape1', 0.045, ...
    'fuel_tank_shape2', 0.0, ...
    'fuel_tank_shape3', 0.0, ...
    'fuel_tank_pillness', 0.5);
geom_ft = twinjet_fuel_tank_geometry(gparams_ft);

% ------------------------------------------------------------
% Build the refined layer-potential objects (same refine/bisect
% pattern as files 1 and 2)
% ------------------------------------------------------------
disp("Building the Geometries")
lp_no_ft = build_lp(geom_no_ft, params);
lp_ft    = build_lp(geom_ft,    params);

% ------------------------------------------------------------
% SVD nullspace-vector density at each case's own k
% ------------------------------------------------------------
disp("Taking the SVD")
[density_no_ft, sv_no_ft] = svd_null_density(lp_no_ft, k_no_fuel_tank);
[density_ft,    sv_ft]    = svd_null_density(lp_ft,    k_fuel_tank);

fprintf('Smallest singular value (no fuel tank, k=%g): %.4e\n', ...
    k_no_fuel_tank, sv_no_ft);
fprintf('Smallest singular value (fuel tank,    k=%g): %.4e\n', ...
    k_fuel_tank, sv_ft);

% ------------------------------------------------------------
% Far field evaluation: 0 to 360 degrees
% ------------------------------------------------------------
disp("Evaluating the far field")
num_angles = 359;
theta = linspace(0, 2*pi - 2*pi/num_angles, num_angles);
xhat = [cos(theta); sin(theta)];

ff_no_ft = lp_no_ft.eval_far_field(k_no_fuel_tank, density_no_ft, xhat);
ff_ft    = lp_ft.eval_far_field(k_fuel_tank,    density_ft,    xhat);

theta_deg = rad2deg(theta);
%%
% ------------------------------------------------------------
% Convert to plotted quantity based on use_db flag
% ------------------------------------------------------------
if use_db
    val_no_ft = 20 * log10(abs(ff_no_ft));
    val_ft    = 20 * log10(abs(ff_ft));
    val_label = '|Far Field| (dB)';
else
    val_no_ft = abs(ff_no_ft);
    val_ft    = abs(ff_ft);
    val_label = '|Far Field|';
end

% ------------------------------------------------------------
% Plot: side-by-side polar plots
% ------------------------------------------------------------
figure
tiledlayout(1,2);

nexttile
polarplot(theta, val_no_ft, '-k', 'LineWidth', 2);
title(sprintf('No Fuel Tank (k=%g)', k_no_fuel_tank))
ax1 = gca; ax1.FontSize = 15;

nexttile
polarplot(theta, val_ft, '-r', 'LineWidth', 2);
title(sprintf('Fuel Tank (k=%g)', k_fuel_tank))
ax2 = gca; ax2.FontSize = 15;

sgtitle(val_label)


%% ================================================================
% Helper functions
% ================================================================

function lp = build_lp(geom, params)
%BUILD_LP Build a refined curve and wrap it in an RP2LP object.

    C = build_curve(geom, params);
    C = refine_curve_wavenumber(C, 2*pi/params.k, 1);

    for jj = 1:2
        C = bisect_all_patches(C);
        C = build_global_indexing(C);
        C = get_close_points(C);
    end

    lp = RP2LP(params, geom);
    lp = lp.update_RP_Curve(C);
end

function [density, smallest_sv] = svd_null_density(lp, k)
%SVD_NULL_DENSITY Right singular vector for the smallest singular value
%   of the BIE matrix at wavenumber k, used as an (approximate)
%   nullspace/resonance density.

    A = lp.bie_mat(k);
    [~, S, V] = svd(A);
    density = V(:, end);
    smallest_sv = S(end, end);
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

    x = P.x(t);
    y = P.y(t);

    if ~isfield(opts, "plot_color")
        plot(x, y, "-", "LineWidth", opts.line_width);
    else
        plot(x, y, "-", "LineWidth", opts.line_width,'Color',opts.plot_color);
    end

    if opts.show_nodes
        plot(P.xj, P.yj, ".", "MarkerSize", opts.node_size);
    end

    if opts.show_labels
        tm = 0.5 * (tmin + tmax);
        xm = P.x(tm);
        ym = P.y(tm);
        text(xm, ym, sprintf("%d", j), ...
            "VerticalAlignment", "bottom", ...
            "HorizontalAlignment", "center");
    end
end

axis equal;
grid on;

if ~hold_state
    hold off;
end

end