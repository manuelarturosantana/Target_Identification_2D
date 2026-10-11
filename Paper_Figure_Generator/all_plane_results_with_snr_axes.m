% Plot the aggregate plane-classification accuracy with reference SNR axes.
clear;
script_directory = fileparts(mfilename('fullpath'));
project_root = fileparts(script_directory);
results_file = fullfile(project_root, 'Identification', ...
    'classification_results', 'plane_results_averaged.mat');

% These are the four examples selected in
% rp2d_matlab/local/plot_plane_shapes.ipynb (top row, then bottom row).
shape_directory = fullfile(fileparts(project_root), 'rp2d_matlab', ...
    'plane_shapes_mat');
shape_filenames = { ...
    'twinjet_plane2.mat', ...
    'twinjet_plane5.mat', ...
    'quadjet_plane3.mat', ...
    'quadjet_plane5.mat'};

% Define return angle half-widths in degrees to test.
angle_bounds = [15, 30, 60, 90];
num_bounds = length(angle_bounds);

loaded = load(results_file);
results = loaded.results;
trials = results.trials;

sigmas = unique(trials.sigma);
nsigma = numel(sigmas);

all_obstacles = unique(trials.true_idx);
num_obstacles = numel(all_obstacles);

% Convert receiver angles to degrees and wrap to [-180, 180].
rec_angles_deg = rad2deg(results.snr.receiver_angle);
angles_wrapped = mod(rec_angles_deg + 180, 360) - 180;

%% Calculate accuracy across all noise levels, angle ranges, and planes.
overall_accuracy = zeros(nsigma, num_bounds + 1);

for k = 1:nsigma
    for r = 1:num_bounds
        mask_r = abs(angles_wrapped) <= angle_bounds(r);
        pooled_correct_r = [];

        for obs_i = 1:num_obstacles
            obs_idx = all_obstacles(obs_i);
            idx = (trials.sigma == sigmas(k)) & ...
                (trials.true_idx == obs_idx);
            correct_k = trials.correct(idx);
            pooled_correct_r = [pooled_correct_r; correct_k(mask_r)]; %#ok<AGROW>
        end

        overall_accuracy(k, r) = mean(pooled_correct_r);
    end

    pooled_correct_all = [];

    for obs_i = 1:num_obstacles
        obs_idx = all_obstacles(obs_i);
        idx = (trials.sigma == sigmas(k)) & ...
            (trials.true_idx == obs_idx);
        pooled_correct_all = [pooled_correct_all; trials.correct(idx)]; %#ok<AGROW>
    end

    overall_accuracy(k, num_bounds + 1) = mean(pooled_correct_all);
end

%% Find the mean SNR across all planes at the receiver nearest 0 degrees.
receiver_distance = abs(angle(exp(1i * results.snr.receiver_angle(:))));
[~, zero_receiver_idx] = min(receiver_distance);

snr_sigmas = results.snr.sigma(:);
[has_snr, snr_rows] = ismember(sigmas, snr_sigmas);

if ~all(has_snr)
    error('all_plane_results_with_snr_axes:MissingSNR', ...
        'The saved SNR data does not contain every classification sigma.');
end

representative_snr_full = ...
    results.snr.snr_full_mean_db(snr_rows, zero_receiver_idx);
representative_snr_late = ...
    results.snr.snr_late_mean_db(snr_rows, zero_receiver_idx);

% Give sigma = 0 a display-only location on the logarithmic axis.
positive_sigmas = sigmas(sigmas > 0);

if isempty(positive_sigmas)
    zero_plot_sigma = 1e-7;
else
    zero_plot_sigma = min(positive_sigmas) / 10;
end

plot_sigmas = sigmas;
plot_sigmas(sigmas == 0) = zero_plot_sigma;

%% Plot using the formatting from all_plane_results.m.
plot_styles = {'-o', '-s', '-^', '-d', '--k*'};
range_labels = { ...
    '[-15^\circ, 15^\circ]', ...
    '[-30^\circ, 30^\circ]', ...
    '[-60^\circ, 60^\circ]', ...
    '[-90^\circ, 90^\circ]', ...
    '[-180^\circ, 180^\circ]'};

fig1 = figure(101); clf;
fig1.Name = 'Overall Classification Accuracy across All Planes';
ax_sigma = axes(fig1);

semilogx(ax_sigma, ...
    plot_sigmas, overall_accuracy(:, 1), plot_styles{1}, ...
    plot_sigmas, overall_accuracy(:, 2), plot_styles{2}, ...
    plot_sigmas, overall_accuracy(:, 3), plot_styles{3}, ...
    plot_sigmas, overall_accuracy(:, 4), plot_styles{4}, ...
    plot_sigmas, overall_accuracy(:, 5), plot_styles{5}, ...
    'LineWidth', 2.5);

grid(ax_sigma, 'on');
xlabel(ax_sigma, "Noise Level (\sigma)", "FontSize", 14);
ylabel(ax_sigma, "Overall Correct Classification Rate", "FontSize", 14);
legend(ax_sigma, range_labels, 'Location', 'SouthWest', 'FontSize', 12);
ax_sigma.FontSize = 14;
ax_sigma.XTick = plot_sigmas;
ax_sigma.XTickLabel = sigma_tick_labels(sigmas);

% Keep endpoint ticks slightly inside the axes so MATLAB does not clip
% them against the left and right boundaries.
x_padding_factor = 1.15;
ax_sigma.XLim = [min(plot_sigmas) / x_padding_factor, ...
    max(plot_sigmas) * x_padding_factor];

% Shorten only the main axes vertically to make room for the two upper
% reference scales. No plot title is used.
main_position = ax_sigma.Position;
main_position(4) = main_position(4) - 0.22;
ax_sigma.Position = main_position;

ax_full = top_snr_axis(ax_sigma, plot_sigmas, ...
    representative_snr_full, ...
    'Representative full-time SNR at 0^\circ [dB]', 0.03);
ax_late = top_snr_axis(ax_sigma, plot_sigmas, ...
    representative_snr_late, ...
    'Representative late-time SNR at 0^\circ [dB]', 0.16);
linkaxes([ax_sigma, ax_full, ax_late], 'x');

% Keep the three horizontal scales aligned after text layout is resolved.
drawnow;
main_position = ax_sigma.Position;
ax_full.Position([1, 3]) = main_position([1, 3]);
ax_late.Position([1, 3]) = main_position([1, 3]);

%% Add the four paper-example planes as independently movable inset axes.
[plane_shapes, shape_limits] = load_plane_shapes( ...
    shape_directory, shape_filenames);

% Positions are normalized figure coordinates. Each inset is an independent
% axes object, so it can be selected and dragged with MATLAB's plot editor.
inset_width = 0.14;
inset_height = inset_width * fig1.Position(3) / fig1.Position(4);
inset_gap_x = 0.015;
inset_gap_y = 0.015;
inset_left = 0.39;
inset_bottom = 0.22;
plane_axes = gobjects(1, numel(plane_shapes));

for ishape = 1:numel(plane_shapes)
    inset_row = ceil(ishape / 2);
    inset_column = mod(ishape - 1, 2) + 1;
    inset_x = inset_left + (inset_column - 1) * ...
        (inset_width + inset_gap_x);
    inset_y = inset_bottom + (2 - inset_row) * ...
        (inset_height + inset_gap_y);

    plane_axes(ishape) = axes( ...
        'Parent', fig1, ...
        'Units', 'normalized', ...
        'Position', [inset_x, inset_y, inset_width, inset_height], ...
        'Color', 'w', ...
        'Box', 'on', ...
        'XTick', [], ...
        'YTick', [], ...
        'Tag', sprintf('Movable plane inset: %s', shape_filenames{ishape}));
    hold(plane_axes(ishape), 'on');

    for ipatch = 1:numel(plane_shapes{ishape}.X)
        plane_line = plot(plane_axes(ishape), ...
            plane_shapes{ishape}.X{ipatch}, ...
            plane_shapes{ishape}.Y{ipatch}, ...
            'k-', 'LineWidth', 1.25);
        plane_line.HitTest = 'off';
        plane_line.PickableParts = 'none';
    end

    xlim(plane_axes(ishape), shape_limits(1:2));
    ylim(plane_axes(ishape), shape_limits(3:4));
    axis(plane_axes(ishape), 'equal');
    plane_axes(ishape).PlotBoxAspectRatio = [1, 1, 1];
    hold(plane_axes(ishape), 'off');
end

% In the desktop, leave plot-edit mode active so the inset axes can be
% selected by their frames and dragged immediately.
if usejava('desktop')
    plotedit(fig1, 'on');
end


function labels = sigma_tick_labels(sigma_values)
%SIGMA_TICK_LABELS Format fixed logarithmic ticks, including sigma = 0.

labels = strings(size(sigma_values));

for ivalue = 1:numel(sigma_values)
    if sigma_values(ivalue) == 0
        labels(ivalue) = "0";
    else
        labels(ivalue) = sprintf('%.0e', sigma_values(ivalue));
    end
end
end


function labels = snr_tick_labels(snr_values)
%SNR_TICK_LABELS Format SNR values for the two upper axes.

labels = strings(size(snr_values));

for ivalue = 1:numel(snr_values)
    if isinf(snr_values(ivalue))
        labels(ivalue) = "Inf";
    else
        labels(ivalue) = sprintf('%.1f', snr_values(ivalue));
    end
end
end


function ax = top_snr_axis(base_axis, x_ticks, snr_values, label, offset)
%TOP_SNR_AXIS Add a negligible-height x-axis above the main plot.

position = base_axis.Position;
position = [position(1), position(2) + position(4) + offset, ...
    position(3), 1e-12];
ax = axes( ...
    'Parent', base_axis.Parent, ...
    'Position', position, ...
    'Color', 'none', ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'YColor', 'none', ...
    'YTick', [], ...
    'XScale', 'log', ...
    'XLim', base_axis.XLim, ...
    'XTick', x_ticks, ...
    'XTickLabel', snr_tick_labels(snr_values), ...
    'FontSize', base_axis.FontSize, ...
    'Box', 'off');
xlabel(ax, label, 'FontSize', 14);
end


function [plane_shapes, limits] = load_plane_shapes(directory, filenames)
%LOAD_PLANE_SHAPES Load patch coordinates and find one shared square scale.

plane_shapes = cell(size(filenames));
all_x = [];
all_y = [];

for ishape = 1:numel(filenames)
    shape_path = fullfile(directory, filenames{ishape});

    if ~isfile(shape_path)
        error('all_plane_results_with_snr_axes:MissingPlaneShape', ...
            'Plane geometry file not found: %s', shape_path);
    end

    plane_shapes{ishape} = load(shape_path, 'X', 'Y');

    if ~iscell(plane_shapes{ishape}.X) || ...
            ~iscell(plane_shapes{ishape}.Y)
        error('all_plane_results_with_snr_axes:BadPlaneShape', ...
            'Expected X and Y patch cell arrays in %s.', shape_path);
    end

    all_x = [all_x; vertcat(plane_shapes{ishape}.X{:})]; %#ok<AGROW>
    all_y = [all_y; vertcat(plane_shapes{ishape}.Y{:})]; %#ok<AGROW>
end

axis_padding = 0.03;
x_center = 0.5 * (min(all_x) + max(all_x));
y_center = 0.5 * (min(all_y) + max(all_y));
square_span = max(max(all_x) - min(all_x), ...
    max(all_y) - min(all_y)) * (1 + 2 * axis_padding);
half_span = 0.5 * square_span;
limits = [x_center - half_span, x_center + half_span, ...
    y_center - half_span, y_center + half_span];
end
