% read_averaged_classification_results.m
%
clear;
results_file_fuel_tank = 'Identification/classification_results/plane_results_averaged_fuel_tank.mat';
results_file_baseline  = 'Identification/classification_results/plane_results_averaged.mat';

% Define angle range half-widths in degrees
angle_bounds = [15, 30, 60, 90];
num_bounds   = length(angle_bounds);

tiledlayout(1,2);

obstacle_idx = 6; % Target obstacle index, used for both results files
show_legend  = true;

results_fuel_tank = load(results_file_fuel_tank);
results_fuel_tank = results_fuel_tank.results;

results_baseline = load(results_file_baseline);
results_baseline = results_baseline.results;

%% Data Extraction
[sigmas_fuel_tank, accuracy_ranges_fuel_tank] = compute_accuracy_ranges( ...
    results_fuel_tank, obstacle_idx, angle_bounds);

[sigmas_baseline, accuracy_ranges_baseline] = compute_accuracy_ranges( ...
    results_baseline, obstacle_idx, angle_bounds);

%% Plot the plane
nexttile

params = build_params('k', 30);
gparams = build_twinjet_fuel_tank_params('preset', 'plane1', ...
    'fuel_tank_loc', 0.45, ...
    'fuel_tank_length', 0.15, ...
    'fuel_tank_width', 0.006, ...
    'fuel_tank_shape1', 0.045, ...
    'fuel_tank_shape2', 0.0, ...
    'fuel_tank_shape3', 0.0, ...
    'fuel_tank_pillness', 0.5);

geom = twinjet_fuel_tank_geometry(gparams);
Curve = build_curve(geom, params);

opts = struct();
opts.show_nodes = false;
opts.plot_color = 'k';
opts.line_width = 1.4;
plot_curve2(Curve, opts);
% title('Twinjet with closed pill-shaped fuel tanks');

%% Plot Results: Accuracy vs Noise Level by Angle Range
nexttile

markers = {'o', 's', '^', 'd'};
colors  = {'b', 'r', 'g', 'm'};
range_labels = { ...
    '[-15^\circ, 15^\circ]', ...
    '[-30^\circ, 30^\circ]', ...
    '[-60^\circ, 60^\circ]', ...
    '[-90^\circ, 90^\circ]'};

% Single semilogx call (hold on/off does not mix well with semilogx):
% baseline (solid) and fuel tank (dotted) curves, paired by color/marker
% per angle bound.
semilogx( ...
    sigmas_baseline,  accuracy_ranges_baseline(:, 1),  ['-' markers{1} colors{1}], ...
    sigmas_fuel_tank, accuracy_ranges_fuel_tank(:, 1), [':' markers{1} colors{1}], ...
    sigmas_baseline,  accuracy_ranges_baseline(:, 2),  ['-' markers{2} colors{2}], ...
    sigmas_fuel_tank, accuracy_ranges_fuel_tank(:, 2), [':' markers{2} colors{2}], ...
    sigmas_baseline,  accuracy_ranges_baseline(:, 3),  ['-' markers{3} colors{3}], ...
    sigmas_fuel_tank, accuracy_ranges_fuel_tank(:, 3), [':' markers{3} colors{3}], ...
    sigmas_baseline,  accuracy_ranges_baseline(:, 4),  ['-' markers{4} colors{4}], ...
    sigmas_fuel_tank, accuracy_ranges_fuel_tank(:, 4), [':' markers{4} colors{4}], ...
    'LineWidth', 2);

legend_entries = strings(1, 2 * num_bounds);
for r = 1:num_bounds
    legend_entries(2*r - 1) = sprintf('%s (No Fuel Tank)', range_labels{r});
    legend_entries(2*r)     = sprintf('%s (Fuel Tank)', range_labels{r});
end

grid on
xlabel("Noise Level")
ylabel("Accuracy")
ax = gca; ax.FontSize = 15;

if show_legend
    legend(legend_entries, 'Location', 'best')
end


function [sigmas, accuracy_ranges] = compute_accuracy_ranges(results, obstacle_idx, angle_bounds)
%COMPUTE_ACCURACY_RANGES Classification accuracy for one obstacle across
%   noise levels (sigma) and receiver-angle bounds.
%
% [sigmas, accuracy_ranges] = compute_accuracy_ranges(results, obstacle_idx, angle_bounds)
%
% accuracy_ranges(k, r) is the accuracy at sigmas(k) restricted to
% receiver angles within +/- angle_bounds(r) degrees of boresight.

trials = results.trials;
sigmas = unique(trials.sigma);
nsigma = numel(sigmas);
num_bounds = numel(angle_bounds);

% Convert receiver angles to degrees and wrap to [-180, 180] range
rec_angles_deg = rad2deg(results.snr.receiver_angle);
angles_wrapped = mod(rec_angles_deg + 180, 360) - 180;

accuracy_ranges = zeros(nsigma, num_bounds);

inds = trials.true_idx == obstacle_idx;
for k = 1:nsigma
    idx = (trials.sigma == sigmas(k)) & inds;
    correct_k = trials.correct(idx);

    for r = 1:num_bounds
        mask = abs(angles_wrapped) <= angle_bounds(r);
        accuracy_ranges(k, r) = mean(correct_k(mask));
    end
end

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

axis equal;
grid on;

if ~hold_state
    hold off;
end

end

return

%% SNR Check for the difference 

clear;
results_file_fuel_tank = 'Identification/classification_results/plane_results_averaged_fuel_tank.mat';
results_file_baseline  = 'Identification/classification_results/plane_results_averaged.mat';
figure(2)
clf
% Define angle range half-widths in degrees
angle_bounds = [15, 30, 60, 90];
num_bounds   = length(angle_bounds);

tiledlayout(2,2);

obstacle_idx = 6; % Same obstacle in both results files (twinjet, scatterer 6)

% Pair each results file with its display label
result_files = {results_file_baseline, results_file_fuel_tank};
result_labels = ["No Fuel Tank", "Fuel Tank"];

loop_cnt = 0; % Track loop iteration for transposing the layout
ax_top = gobjects(1, 2); % Store axes handles for top row to sync limits

for file_idx = 1:2
    loop_cnt = loop_cnt + 1;

    plane_name = result_labels(file_idx);

    snr_mode = 'late'; % 'early' or 'late'
    show_legend = true;
    results = load(result_files{file_idx});
    results = results.results;
    
    if strcmp(snr_mode,'late')
        snr_data = results.snr.snr_late_db;
    else
        snr_data = results.snr.snr_full_db;
    end
    
    trials = results.trials;
    sigmas = unique(trials.sigma);
    nsigma = numel(sigmas);
    
    accuracy     = zeros(nsigma, 1);
    ntrials      = zeros(nsigma, 1);
    nwrong       = zeros(nsigma, 1);
    correct_inds = zeros(nsigma-1, 359);
    
    % Store accuracy across angle ranges:
    % Columns: [-15, 15], [-30, 30], [-60, 60], [-90, 90], All angles
    accuracy_ranges = zeros(nsigma, num_bounds + 1);

    % Convert receiver angles to degrees and wrap to [-180, 180] range
    rec_angles_deg = rad2deg(results.snr.receiver_angle);
    angles_wrapped = mod(rec_angles_deg + 180, 360) - 180;

    % Data Extraction Loop
    inds = trials.true_idx == obstacle_idx;
    for k = 1:nsigma
        idx = (trials.sigma == sigmas(k)) & inds;
        if k ~= 1
            correct_inds(k-1, :) = trials.correct(idx);
        end
        ntrials(k)  = nnz(idx);
        nwrong(k)   = nnz(~trials.correct(idx));
        accuracy(k) = mean(trials.correct(idx));

        correct_k = trials.correct(idx);

        % Calculate accuracy for each angle range (-bound, bound)
        for r = 1:num_bounds
            mask = abs(angles_wrapped) <= angle_bounds(r);
            accuracy_ranges(k, r) = mean(correct_k(mask));
        end
        % Full sweep accuracy
        accuracy_ranges(k, num_bounds + 1) = accuracy(k);
    end

    % Subplot 1 & 2 (Top Row): SNR vs Receiver Angle
    ax_top(loop_cnt) = nexttile(loop_cnt); % Transposed: Tile 1 or Tile 2

    % Predefine explicit colors for consistent matching across plots
    num_snr_curves = length(results.snr.sigma) - 1;
    snr_colors = lines(num_snr_curves);

    legend_data = {};
    hold on
    for ii = 1:num_snr_curves
        % Plot line using explicit color array
        plot(rad2deg(results.snr.receiver_angle), snr_data(ii+1,:,obstacle_idx), ...
             'Color', snr_colors(ii, :), 'LineWidth', 1.5);
        
        legend_data{end+1} = "$\sigma$ = " + num2str(results.snr.sigma(ii+1));
        
        inds_missed = ~correct_inds(ii,:);
        plot(rad2deg(results.snr.receiver_angle(inds_missed)), snr_data(ii+1,inds_missed,obstacle_idx), 'ko');
        
        if sum(~correct_inds(ii,:) ~= 0)
            legend_data{end+1} = "";
        end
    end
    hold off

    if file_idx == 2
        ldg = legend(legend_data, 'Location', 'SouthEast', 'interpreter', 'latex');
        ldg.FontSize = 10;
    end

    xlabel("Scatterer Orientation in Degrees")
    
    % Apply Y-axis label only on Tile 1; remove ticks on Tile 2
    if loop_cnt == 1
        ylabel("Late Time SNR")
    else
        yticklabels([])
    end

    xlim([0, 360])
    title(plane_name)
    ax = gca; ax.FontSize = 15;
    grid on

    % Subplot 3 & 4 (Bottom Row): Accuracy vs Noise Level by Angle Range
    nexttile(loop_cnt + 2) % Transposed: Tile 3 or Tile 4

    plot_styles = {'-o', '-s', '-^', '-d', '--k*'};
    range_labels = { ...
        '[-15^\circ, 15^\circ]', ...
        '[-30^\circ, 30^\circ]', ...
        '[-60^\circ, 60^\circ]', ...
        '[-90^\circ, 90^\circ]', ...
        '[-180^\circ, 180^\circ]'};

    % Single semilogx call to preserve log scaling robustly
    semilogx( ...
        sigmas, accuracy_ranges(:, 1), plot_styles{1}, ...
        sigmas, accuracy_ranges(:, 2), plot_styles{2}, ...
        sigmas, accuracy_ranges(:, 3), plot_styles{3}, ...
        sigmas, accuracy_ranges(:, 4), plot_styles{4}, ...
        sigmas, accuracy_ranges(:, 5), plot_styles{5}, ...
        'LineWidth', 2);

    grid on
    xlabel("Noise Level")
    ylabel("Accuracy")
    ax = gca; ax.FontSize = 15;

    if file_idx == 2
        legend(range_labels, 'Location', 'best')
    end
end

% Link Y-axis limits across top row
linkaxes(ax_top, 'y');