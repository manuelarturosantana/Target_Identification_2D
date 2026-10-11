% plane_results_with_snr.m
clear;

results_file = ...
    'Identification/classification_results/plane_results_averaged.mat';

% Plot parameters.
obs_indices = [6, 9];
angle_bounds = [15, 30, 60, 90];
opening_bounds = [angle_bounds, 180];
snr_mode = 'late';

sigma_legend_columns = 3; % Set to 0 for one column per sigma
maximum_noise = 1e-1;

outer_x_padding_factor = 1.06;
inner_x_padding_factor = 1.50;

axes_font_size = 16;
label_font_size = 18;
title_font_size = 20;
legend_font_size = 11;
snr_tick_font_size = 13;
snr_label_font_size = 13;
accuracy_marker_size = 10;

loaded = load(results_file);
results = loaded.results;
trials = results.trials;

sigmas = unique(trials.sigma);
nsigma = numel(sigmas);

% Show every positive sigma up to the maximum plotted noise.
noise_ticks = sigmas( ...
    sigmas > 0 & sigmas <= maximum_noise);

noise_tick_labels = compose( ...
    '10^{%d}', round(log10(noise_ticks)));

angles = rad2deg(results.snr.receiver_angle(:));
angles = mod(angles + 180, 360) - 180;
num_receivers = numel(angles);

if strcmp(snr_mode, 'late')
    snr_data = results.snr.snr_late_db;
    snr_label = 'Late-time SNR [dB]';
elseif strcmp(snr_mode, 'full')
    snr_data = results.snr.snr_full_db;
    snr_label = 'Full-time SNR [dB]';
else
    error('snr_mode must be either ''late'' or ''full''.');
end

% Match finite SNR rows to classification noise levels.
snr_sigmas = results.snr.sigma(:);
candidate_rows = find(snr_sigmas > 0);

[matched, classification_rows] = ismember( ...
    snr_sigmas(candidate_rows), sigmas);

snr_rows = candidate_rows(matched);
classification_rows = classification_rows(matched);
num_snr_curves = numel(snr_rows);

if sigma_legend_columns == 0
    num_legend_columns = num_snr_curves;
else
    num_legend_columns = min( ...
        sigma_legend_columns, num_snr_curves);
end

% Find the receiver nearest 0 degrees.
receiver_distance = abs(angle(exp( ...
    1i * results.snr.receiver_angle(:))));

[~, zero_receiver_idx] = min(receiver_distance);

% Extract accuracy and correctness.
num_obstacles = numel(obs_indices);
num_openings = numel(opening_bounds);

accuracy = zeros(nsigma, num_openings, num_obstacles);
correctness = false(nsigma, num_receivers, num_obstacles);

for iobs = 1:num_obstacles
    obstacle_idx = obs_indices(iobs);
    obstacle_mask = trials.true_idx == obstacle_idx;

    for isigma = 1:nsigma
        idx = obstacle_mask & trials.sigma == sigmas(isigma);
        correct = logical(trials.correct(idx));

        correctness(isigma, :, iobs) = correct(:).';

        for iopening = 1:num_openings
            angle_mask = abs(angles) <= opening_bounds(iopening);

            accuracy(isigma, iopening, iobs) = ...
                mean(correct(angle_mask));
        end
    end
end

% Distinct categorical colors for the bottom plots.
aperture_colors = [ ...
    0.835, 0.369, 0.000; ... % orange
    0.000, 0.620, 0.451; ... % bluish green
    0.800, 0.475, 0.655; ... % reddish purple
    0.337, 0.706, 0.914; ... % sky blue
    0.000, 0.000, 0.000];    % full opening

% Compact, approximately square layout.
top_positions = [ ...
    0.090, 0.56, 0.421, 0.37; ...
    0.539, 0.56, 0.421, 0.37];

bottom_positions = [ ...
    0.090, 0.08, 0.421, 0.30; ...
    0.539, 0.08, 0.421, 0.30];

fig = figure(101); clf;
fig.Name = 'Plane Classification Results';
fig.Position(3:4) = [1000, 1000];

ax_top = gobjects(1, num_obstacles);
ax_bottom = gobjects(1, num_obstacles);
ax_full = gobjects(1, num_obstacles);
ax_late = gobjects(1, num_obstacles);

minimum_noise = min(noise_ticks);

for iobs = 1:num_obstacles
    obstacle_idx = obs_indices(iobs);
    plane_name = obstacle_name(obstacle_idx);

    %% Top row: SNR versus orientation.
    ax_top(iobs) = axes( ...
        'Parent', fig, ...
        'Position', top_positions(iobs, :));

    hold(ax_top(iobs), 'on');

    [angles_sorted, sort_idx] = sort(angles);

    noise_colors = lines(num_snr_curves);
    noise_handles = gobjects(num_snr_curves, 1);
    noise_labels = strings(num_snr_curves, 1);

    for icurve = 1:num_snr_curves
        snr_row = snr_rows(icurve);
        sigma_value = snr_sigmas(snr_row);

        curve = reshape( ...
            snr_data(snr_row, sort_idx, obstacle_idx), 1, []);

        missed = ~reshape( ...
            correctness( ...
                classification_rows(icurve), sort_idx, iobs), ...
            1, []);

        noise_handles(icurve) = plot( ...
            ax_top(iobs), ...
            angles_sorted, curve, ...
            'Color', noise_colors(icurve, :), ...
            'LineWidth', 1.7);

        plot(ax_top(iobs), ...
            angles_sorted(missed), ...
            curve(missed), ...
            'ko', ...
            'MarkerSize', 4.5, ...
            'HandleVisibility', 'off');

        noise_labels(icurve) = ...
            sprintf('$\\sigma = %.0e$', sigma_value);
    end

    hold(ax_top(iobs), 'off');
    grid(ax_top(iobs), 'on');
    xlim(ax_top(iobs), [-180, 180]);

    xlabel(ax_top(iobs), ...
        'Scatterer Orientation [degrees]', ...
        'FontSize', label_font_size);

    title(ax_top(iobs), plane_name, ...
        'FontSize', title_font_size);

    ax_top(iobs).FontSize = axes_font_size;

    if iobs == 1
        ylabel(ax_top(iobs), snr_label, ...
            'FontSize', label_font_size);
    else
        ax_top(iobs).YTickLabel = [];
    end

    if iobs == num_obstacles
        noise_legend = legend( ...
            ax_top(iobs), ...
            noise_handles, noise_labels, ...
            'Location', 'southeast', ...
            'Interpreter', 'latex', ...
            'FontSize', legend_font_size, ...
            'NumColumns', num_legend_columns);

        noise_legend.Title.String = 'Noise level';
        noise_legend.Title.FontSize = legend_font_size;
    end

    %% Bottom row: accuracy versus noise level.
    ax_bottom(iobs) = axes( ...
        'Parent', fig, ...
        'Position', bottom_positions(iobs, :));

    ax_bottom(iobs).XScale = 'log';
    hold(ax_bottom(iobs), 'on');

    aperture_handles = gobjects(num_openings, 1);
    aperture_labels = strings(num_openings, 1);

    aperture_markers = {'o', 's', '^', 'd', '*'};
    aperture_line_styles = {'-', '-', '-', '-', '--'};

    for iopening = 1:num_openings
        aperture_handles(iopening) = semilogx( ...
            ax_bottom(iobs), ...
            sigmas, ...
            accuracy(:, iopening, iobs), ...
            'Color', aperture_colors(iopening, :), ...
            'LineStyle', aperture_line_styles{iopening}, ...
            'Marker', aperture_markers{iopening}, ...
            'MarkerSize', accuracy_marker_size, ...
            'LineWidth', 2.2);

        aperture_labels(iopening) = sprintf( ...
            '[-%d^\\circ, %d^\\circ]', ...
            opening_bounds(iopening), ...
            opening_bounds(iopening));
    end

    hold(ax_bottom(iobs), 'off');
    grid(ax_bottom(iobs), 'on');

    % Add more log-scale padding along the inner edge.
    if iobs == 1
        x_limits = [ ...
            minimum_noise / outer_x_padding_factor, ...
            maximum_noise * inner_x_padding_factor];
    else
        x_limits = [ ...
            minimum_noise / inner_x_padding_factor, ...
            maximum_noise * outer_x_padding_factor];
    end

    ax_bottom(iobs).XLim = x_limits;
    ax_bottom(iobs).XTick = noise_ticks;
    ax_bottom(iobs).XTickLabel = noise_tick_labels;
    ax_bottom(iobs).TickLabelInterpreter = 'tex';
    ax_bottom(iobs).XMinorTick = 'on';

    xlabel(ax_bottom(iobs), ...
        'Noise Level (\sigma)', ...
        'FontSize', label_font_size);

    ax_bottom(iobs).FontSize = axes_font_size;

    if iobs == 1
        ylabel(ax_bottom(iobs), 'Accuracy', ...
            'FontSize', label_font_size);
    else
        ax_bottom(iobs).YTickLabel = [];
    end

    if iobs == num_obstacles
        aperture_legend = legend( ...
            ax_bottom(iobs), ...
            aperture_handles, aperture_labels, ...
            'Location', 'southwest', ...
            'FontSize', legend_font_size);

        aperture_legend.Title.String = 'Angular range';
        aperture_legend.Title.FontSize = legend_font_size;
    end

    %% Add compact full-time and late-time SNR axes.
    representative_snr_full = reshape( ...
        results.snr.snr_full_db( ...
            snr_rows, zero_receiver_idx, obstacle_idx), ...
        [], 1);

    representative_snr_late = reshape( ...
        results.snr.snr_late_db( ...
            snr_rows, zero_receiver_idx, obstacle_idx), ...
        [], 1);

    snr_tick_sigmas = snr_sigmas(snr_rows);

    visible_ticks = ...
        snr_tick_sigmas >= minimum_noise & ...
        snr_tick_sigmas <= maximum_noise;

    snr_tick_sigmas = snr_tick_sigmas(visible_ticks);
    representative_snr_full = ...
        representative_snr_full(visible_ticks);
    representative_snr_late = ...
        representative_snr_late(visible_ticks);

    if iobs == 1
        full_label = 'Full-time';
        late_label = 'Late-time';
    else
        full_label = '';
        late_label = '';
    end

    ax_full(iobs) = compact_snr_axis( ...
        ax_bottom(iobs), ...
        snr_tick_sigmas, ...
        representative_snr_full, ...
        full_label, ...
        0.012, ...
        snr_tick_font_size, ...
        snr_label_font_size);

    ax_late(iobs) = compact_snr_axis( ...
        ax_bottom(iobs), ...
        snr_tick_sigmas, ...
        representative_snr_late, ...
        late_label, ...
        0.070, ...
        snr_tick_font_size, ...
        snr_label_font_size);

    linkaxes( ...
        [ax_bottom(iobs), ax_full(iobs), ax_late(iobs)], ...
        'x');
end

linkaxes(ax_top, 'y');
linkaxes(ax_bottom, 'y');


function name = obstacle_name(obstacle_idx)
%OBSTACLE_NAME Return the displayed plane name.

if obstacle_idx < 6
    name = "Quadjet " + num2str(obstacle_idx);
else
    name = "Twinjet " + num2str(obstacle_idx - 5);
end
end


function labels = snr_tick_labels(snr_values)
%SNR_TICK_LABELS Format values for the upper SNR axes.

labels = strings(size(snr_values));

for ivalue = 1:numel(snr_values)
    if isinf(snr_values(ivalue))
        labels(ivalue) = "Inf";
    else
        labels(ivalue) = sprintf('%.1f', snr_values(ivalue));
    end
end
end


function ax = compact_snr_axis( ...
    base_axis, x_ticks, snr_values, label_text, offset, ...
    tick_font_size, label_font_size)
%COMPACT_SNR_AXIS Add a compact SNR scale above an accuracy plot.

position = base_axis.Position;
position = [ ...
    position(1), ...
    position(2) + position(4) + offset, ...
    position(3), ...
    1e-12];

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
    'FontSize', tick_font_size, ...
    'Box', 'off');

% Use a figure annotation so the label survives export.
if ~isempty(label_text)
    label_width = 0.078;
    label_height = 0.028;
    label_gap = 0.006;

    label_position = [ ...
        position(1) - label_width - label_gap, ...
        position(2) - label_height / 2, ...
        label_width, ...
        label_height];

    annotation(base_axis.Parent, ...
        'textbox', label_position, ...
        'String', label_text, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', label_font_size, ...
        'FontWeight', 'normal', ...
        'Margin', 0);
end
end