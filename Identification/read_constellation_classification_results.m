% plot_target_identification_constellation.m
%
% Read and plot four-method GLRT constellation classification results.
%
% Plots:
%   1. Accuracy versus sigma.
%   2. Accuracy versus representative SNR, if compatible SNR data exists.
%   3. Object/constellation shape gallery.

clear
clc

% ---- user-defined parameters ----
results_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'constellation_four_methods_results.mat');

plot_output_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'figures_four_methods');

save_figures = true;

% Shape plotting.
plot_all_shapes = true;

if ~exist(plot_output_folder, 'dir')
    mkdir(plot_output_folder);
end

S = load(results_path, 'results');
results = S.results;

params = results.params;
sigmas = params.noise_sigmas(:);

method_fields = [
    "method_full_naive"
    "method_remove_permutations"
    "method_object_multiset"
    "method_single_poles"
];

method_labels = [
    "1. Full naive"
    "2. Remove permutations"
    "3. Remove geometry"
    "4. Single poles"
];

plot_styles = [
    "o-"
    "s-"
    "^-"
    "d-"
];

% ============================================================
% FIGURE 1: ACCURACY VS SIGMA WITH SNR TOP AXIS
% ============================================================

% Choose which SNR to display on the top x-axis.
% Options:
%   "snr_late_db"
%   "snr_full_db"
snr_field = "snr_late_db";

% Title position control.
% title_pos = [x_center, y_bottom]
%
% Increase title_pos(1) -> move title right.
% Decrease title_pos(1) -> move title left.
% Increase title_pos(2) -> move title up.
% Decrease title_pos(2) -> move title down.
title_pos = [0.50, 0.94];

figure(1)
clf
hold on
box on
grid on

for imethod = 1:numel(method_fields)
    method_result = results.(method_fields(imethod));
    acc = get_accuracy_vector(method_result, sigmas);

    plot_sigma_accuracy( ...
        sigmas, ...
        acc, ...
        plot_styles(imethod), ...
        method_labels(imethod));
end

xlabel('\sigma')
ylabel('Classification accuracy')
legend('Location', 'best')
ylim([0, 1.05])
set(gca, 'XScale', 'log')

ax_bottom = gca;

% Main axes position:
% [left, bottom, width, height]
%
% This controls the plot box, not the title.
% Reduce height if the top SNR axis overlaps the title.
ax_bottom.Position = [0.13, 0.14, 0.78, 0.64];

fix_zero_sigma_axis(ax_bottom, sigmas)

add_snr_top_axis( ...
    ax_bottom, ...
    sigmas, ...
    results.method_full_naive.snr, ...
    snr_field);

% Figure-level title.
% title_pos = [x_center, y_center]
%
% Increase title_pos(1) -> move title right.
% Decrease title_pos(1) -> move title left.
% Increase title_pos(2) -> move title up.
% Decrease title_pos(2) -> move title down.
title_pos = [0.50, 0.90];

title_width = 0.50;
title_height = 0.04;

annotation('textbox', ...
    [title_pos(1) - title_width/2, ...
     title_pos(2) - title_height/2, ...
     title_width, ...
     title_height], ...
    'String', 'GLRT classification accuracy', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold');

if save_figures
    exportgraphics(gcf, fullfile(plot_output_folder, 'accuracy_vs_sigma_snr.png'), ...
        'Resolution', 200);
end

% ============================================================
% FIGURE 3: SHAPE GALLERY
% ============================================================

shape_info = results.shape_info;

if plot_all_shapes
    plot_shape_gallery(shape_info, "All object and constellation shapes");
    fig_name = 'all_shapes.png';
else
    keep = representative_shape_mask(shape_info);
    plot_shape_gallery(shape_info(keep), "Representative object shapes");
    fig_name = 'representative_shapes.png';
end

if save_figures
    exportgraphics(gcf, fullfile(plot_output_folder, fig_name), ...
        'Resolution', 200);
end

fprintf('\nSaved plots to:\n  %s\n', plot_output_folder);

% ============================================================
% LOCAL HELPERS
% ============================================================

function acc = get_accuracy_vector(method_result, sigmas)
    if isfield(method_result, 'accuracy_by_sigma')
        acc = method_result.accuracy_by_sigma(:);
    elseif isfield(method_result, 'trials')
        acc = infer_accuracy_by_sigma(method_result.trials, sigmas);
    else
        error('Method result has neither accuracy_by_sigma nor trials field.');
    end

    if numel(acc) ~= numel(sigmas)
        error('Accuracy vector length does not match number of sigmas.');
    end
end

function acc = infer_accuracy_by_sigma(trials, sigmas)
    acc = zeros(numel(sigmas), 1);

    for isig = 1:numel(sigmas)
        mask = trials.sigma == sigmas(isig);
        acc(isig) = mean(trials.correct(mask));
    end
end

function plot_sigma_accuracy(sigmas, acc, style, label)
%PLOT_SIGMA_ACCURACY
%   Plot accuracy versus sigma, including sigma = 0.
%
%   Since zero cannot be shown on a log scale, sigma = 0 is placed at
%   min_positive_sigma / 10 and relabeled as 0 by fix_zero_sigma_axis.

    sigmas = sigmas(:);
    acc = acc(:);

    pos = sigmas > 0;

    if any(sigmas == 0) && any(pos)
        sigma_min = min(sigmas(pos));
        sigma_plot = sigmas;
        sigma_plot(sigmas == 0) = sigma_min / 10;
    else
        sigma_plot = sigmas;
    end

    plot(sigma_plot, acc, style, ...
        'LineWidth', 1.5, ...
        'MarkerSize', 7, ...
        'DisplayName', label);
end

function fix_zero_sigma_axis(ax, sigmas)
    sigmas = sigmas(:);
    pos = sigmas > 0;

    if ~any(sigmas == 0) || ~any(pos)
        return
    end

    sigma_min = min(sigmas(pos));
    sigma_max = max(sigmas(pos));
    sigma_zero_plot = sigma_min / 10;

    xlim(ax, [sigma_zero_plot / 1.5, sigma_max * 1.5])

    ticks = [sigma_zero_plot; sigmas(pos)];
    tick_labels = strings(size(ticks));

    tick_labels(1) = "0";

    for ii = 2:numel(ticks)
        tick_labels(ii) = sprintf('%.0e', ticks(ii));
    end

    xticks(ax, ticks)
    xticklabels(ax, tick_labels)
end

function snr_vec = get_snr_vector(method_result, sigmas)
%GET_SNR_VECTOR
%   Try to extract a compatible SNR vector from method_result.snr.
%
%   This is intentionally tolerant because compute_representative_snr may
%   return different field names depending on your local version.

    snr_vec = [];

    if ~isfield(method_result, 'snr')
        return
    end

    s = method_result.snr;

    candidate_fields = [
        "late_time_snr"
        "snr_late_time"
        "late_snr"
        "SNR_late_time"
        "snr"
    ];

    if isstruct(s)
        for ii = 1:numel(candidate_fields)
            f = candidate_fields(ii);

            if isfield(s, f)
                candidate = s.(f);
                candidate = candidate(:);

                if numel(candidate) == numel(sigmas)
                    snr_vec = candidate;
                    return
                end
            end
        end
    elseif isnumeric(s)
        candidate = s(:);

        if numel(candidate) == numel(sigmas)
            snr_vec = candidate;
            return
        end
    end
end

function keep = representative_shape_mask(shape_info)
    names = arrayfun(@(s) string(s.name), shape_info);

    keep = false(size(names));

    wanted = [
        "single_circ_o050"
        "single_ell_o075"
        "same_dir_same_obj_aligned_circ_o050_circ_o050"
        "same_dir_same_obj_aligned_ell_o075_ell_o075"
        "same_dir_diff_obj_aligned_circ_o050_ell_o075"
        "facing_diff_obj_aligned_circ_o050_ell_o075"
    ];

    for ii = 1:numel(wanted)
        keep = keep | names == wanted(ii);
    end
end

function plot_shape_gallery(shape_info, fig_title)
    nshape = numel(shape_info);

    if nshape == 0
        error('No shape_info entries to plot.');
    end

    ncols = ceil(sqrt(nshape));
    nrows = ceil(nshape / ncols);

    figure
    clf

    for ii = 1:nshape
        subplot(nrows, ncols, ii)
        hold on
        box on
        axis equal

        x = shape_info(ii).curveX;
        y = shape_info(ii).curveY;

        if isempty(x) || isempty(y)
            title(strrep(string(shape_info(ii).name), '_', '\_'), ...
                'Interpreter', 'tex', ...
                'FontSize', 8);
            text(0.5, 0.5, 'No curve data', ...
                'HorizontalAlignment', 'center');
            axis off
            continue
        end

        plot_curve_data(x, y)

        title(strrep(string(shape_info(ii).name), '_', '\_'), ...
            'Interpreter', 'tex', ...
            'FontSize', 8);

        xlabel('x')
        ylabel('y')
    end

    sgtitle(fig_title)
end

function plot_curve_data(x, y)
%PLOT_CURVE_DATA
%   Plot curve data without fake connections.
%
%   This prevents MATLAB from drawing artificial lines between different
%   objects or across cavity openings.

    if isempty(x) || isempty(y)
        return
    end

    if iscell(x)
        for kk = 1:numel(x)
            plot_curve_vector_segmented(x{kk}, y{kk});
        end
        return
    end

    if isvector(x)
        plot_curve_vector_segmented(x(:), y(:));
        return
    end

    if ismatrix(x)
        for kk = 1:size(x, 2)
            plot_curve_vector_segmented(x(:, kk), y(:, kk));
        end
        return
    end

    error('Unsupported curve data format.');
end

function plot_curve_vector_segmented(x, y)
%PLOT_CURVE_VECTOR_SEGMENTED
%   Plot one curve vector with automatic NaN breaks at large jumps.

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        plot(x, y, '.', 'MarkerSize', 8);
        return
    end

    ds = hypot(diff(x), diff(y));
    positive_ds = ds(ds > 0);

    if isempty(positive_ds)
        plot(x, y, '.', 'MarkerSize', 8);
        return
    end

    h = median(positive_ds);

    jump_factor = 6;
    breaks = ds > jump_factor * h;

    xx = [];
    yy = [];

    start_idx = 1;

    for ii = 1:numel(breaks)
        if breaks(ii)
            xx = [xx; x(start_idx:ii); NaN]; %#ok<AGROW>
            yy = [yy; y(start_idx:ii); NaN]; %#ok<AGROW>
            start_idx = ii + 1;
        end
    end

    xx = [xx; x(start_idx:end)];
    yy = [yy; y(start_idx:end)];

    plot(xx, yy, 'LineWidth', 1.0);
end

function add_snr_top_axis(ax_bottom, sigmas, snr_table, snr_field)
%ADD_SNR_TOP_AXIS  Add representative SNR in dB as top x-axis.
%
%   Bottom axis:
%       sigma
%
%   Top axis:
%       SNR in dB corresponding to the same sigma values.
%
%   The sigma = 0 point is plotted at min_positive_sigma / 10 on the
%   bottom axis and labeled as Inf on the SNR axis.

    if ~istable(snr_table)
        error('Expected method_result.snr to be a table.');
    end

    names = string(snr_table.Properties.VariableNames);

    if ~ismember("sigma", names)
        error('SNR table does not contain field "sigma".');
    end

    if ~ismember(snr_field, names)
        error('SNR table does not contain field "%s".', snr_field);
    end

    snr_sigmas = snr_table.sigma(:);
    snr_db = snr_table.(snr_field)(:);

    if numel(snr_sigmas) ~= numel(sigmas)
        error('SNR sigma vector length does not match params.noise_sigmas.');
    end

    if max(abs(snr_sigmas - sigmas(:))) > 1e-14
        error('SNR sigma grid does not match params.noise_sigmas.');
    end

    sigma_ticks = sigma_plot_positions(sigmas);

    snr_labels = strings(size(snr_db));

    for ii = 1:numel(snr_db)
        if isinf(snr_db(ii))
            snr_labels(ii) = "Inf";
        else
            snr_labels(ii) = sprintf('%.1f', snr_db(ii));
        end
    end

    ax_top = axes( ...
        'Position', ax_bottom.Position, ...
        'XAxisLocation', 'top', ...
        'YAxisLocation', 'right', ...
        'Color', 'none', ...
        'XScale', ax_bottom.XScale, ...
        'XLim', ax_bottom.XLim, ...
        'XTick', sigma_ticks, ...
        'XTickLabel', snr_labels, ...
        'YTick', [], ...
        'Box', 'off');

    xlabel(ax_top, 'Representative late-time SNR (dB)')

    % Keep the bottom axis visually dominant.
    ax_top.YColor = 'none';
end

function sigma_plot = sigma_plot_positions(sigmas)
%SIGMA_PLOT_POSITIONS
%   Return the x locations used for plotting sigma, including sigma = 0.
%
%   Since zero cannot be placed on a log axis, sigma = 0 is plotted at
%   min_positive_sigma / 10.

    sigmas = sigmas(:);
    sigma_plot = sigmas;

    pos = sigmas > 0;

    if any(sigmas == 0) && any(pos)
        sigma_min = min(sigmas(pos));
        sigma_plot(sigmas == 0) = sigma_min / 10;
    end
end
