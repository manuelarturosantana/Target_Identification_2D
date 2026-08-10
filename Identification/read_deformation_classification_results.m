% read_deformation_classification_results_paper.m
%
% Read and summarize deformation-classification results. The classifier uses
% nominal poles as its library and evaluates nominal and deformed cavities.
%
% Figures are divided explicitly into paper and supplementary plots. Paper
% plots are saved as vector PDFs when save_paper_plots is true.

clear
clc
close all

% ---- user-defined parameters ----

results_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'nominal_poles_deformation_results.mat');

data_folder_override = fullfile( ...
    'local', ...
    'data_deformation_experiment_2d');

% Empty means use all sigmas in the detailed supplementary confusion table.
sigma_to_report = [];

% ---- paper plots ----

make_paper_geometry_plot = true;
make_paper_accuracy_plot = true;
make_paper_mismatch_plot = true;
make_paper_pole_plot = true;

save_paper_plots = true;
paper_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'deformation_paper_figures');

% Geometry and pole examples use the largest available Fourier mode count.
% Empty selects the first realization common to the plotted configurations.
paper_example_realization = [];

% Empty selects sigma = 0 when available, otherwise the smallest sigma.
mismatch_sigma = [];

% ---- supplementary plots ----

% These reproduce broad diagnostic summaries and are disabled by default.
make_accuracy_sigma_plot = false;
make_accuracy_snr_plot = false;

% Plot the nominal/deformed boundaries exactly as in the driver.
make_deformation_comparison_plots = false;
save_deformation_comparison_plots = false;
deformation_comparison_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'deformation_geometry_comparisons');

% Plot deformed and matching nominal poles together in the complex plane.
% This is independent of make_config_plots.
make_pole_comparison_plots = false;
save_pole_comparison_plots = false;
pole_comparison_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'deformation_pole_comparisons');

make_config_plots = false;

% This supplementary SNR plot pools realizations and classes. It uses the
% saved representative SNR mapping, not a separate SNR for every trial.
% The sigma = 0 point is omitted because its SNR is infinite.

make_shape_panel = true;
make_pole_panel = true;
make_accuracy_text_panel = true;

% Empty means plot all evaluation configurations.
% These are true_idx values.
configs_to_plot = [];

% Empty means use all sigmas for per-configuration accuracy.
sigmas_for_config_accuracy = [];

dock_figures = true;

% 'separate' gives one docked tab per configuration.
% 'reuse' overwrites the same figure at each configuration.
config_figure_mode = 'separate';

save_config_plots = false;
config_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'deformation_config_diagnostics');

% ---- plotting parameters ----

font_name = 'Helvetica';
tick_font_size = 30;
label_font_size = 40;
title_font_size = 40;
subtitle_font_size = 30;
legend_font_size = 30;
text_font_size = 30;

line_width = 2.5;
axes_line_width = 1.5;
marker_size = 9;
shape_marker_size = 7;
pole_marker_size = 9;
nominal_pole_marker_size = 11;

deformed_shape_color = [0.0000, 0.4470, 0.7410];
nominal_shape_color = [0.15, 0.15, 0.15];
deformed_pole_color = [0.0000, 0.4470, 0.7410];
nominal_pole_color = [0.8500, 0.3250, 0.0980];

sigma_char = char(963);
epsilon_char = char(949);
accuracy_ylabel = 'average classification accuracy';

summary_figure_position = [100, 100, 620, 430];
snr_figure_position = [100, 100, 700, 500];
config_figure_position = [100, 100, 1150, 360];

fig_paper_geometry = 1;
fig_paper_accuracy = 2;
fig_paper_mismatch = 3;
fig_paper_poles = 4;
fig_accuracy_sigma = 20;
fig_accuracy_snr = 21;
fig_config_start = 30;

interpreter = 'none';
grid_state = 'on';

% ---------------------------------

if dock_figures
    set(groot, 'DefaultFigureWindowStyle', 'docked')
else
    set(groot, 'DefaultFigureWindowStyle', 'normal')
end

S = load(results_path, 'results');
results = S.results;
trials = results.trials;

if isfield(results, 'eval_objects')
    eval_objects = results.eval_objects;
else
    eval_objects = results.objects;
end

if isfield(results, 'library_objects')
    library_objects = results.library_objects;
else
    library_objects = results.objects;
end

eval_names = table_strings(eval_objects, {'eval_object_name', 'object_name'});
library_names = table_strings(library_objects, {'object_name'});
true_names = eval_names(trials.true_idx);

num_trials = height(trials);
object_types = strings(num_trials, 1);
config_kinds = strings(num_trials, 1);
num_modes = zeros(num_trials, 1);
realization_ids = zeros(num_trials, 1);
amplitudes = zeros(num_trials, 1);

for ii = 1:num_trials
    metadata = deformation_metadata_from_name(true_names(ii));
    object_types(ii) = metadata.object_type;
    config_kinds(ii) = metadata.config_kind;
    num_modes(ii) = metadata.num_modes;
    realization_ids(ii) = metadata.realization_id;
    amplitudes(ii) = metadata.amplitude;
end

snr_by_sigma = table();
snr_axis_label = '';

if make_accuracy_snr_plot
    [snr_by_sigma, snr_axis_label] = extract_snr_by_sigma(results);
end

fprintf('\nLoaded results from:\n  %s\n', results_path);

if isfield(results, 'method_name')
    fprintf('\nMethod: %s\n', results.method_name);
end

if isfield(results, 'method_description')
    fprintf('%s\n', results.method_description);
end

fprintf('\nOverall accuracy: %.4f\n', mean(trials.correct));

deformed_mask = config_kinds == "deformed";
nominal_mask = config_kinds == "nominal";

if any(deformed_mask)
    fprintf('Deformed-only accuracy: %.4f\n', ...
        mean(trials.correct(deformed_mask)));
end

if any(nominal_mask)
    fprintf('Nominal-only accuracy: %.4f\n', ...
        mean(trials.correct(nominal_mask)));
end

acc_by_sigma = grouped_accuracy( ...
    trials.sigma, trials.correct, 'sigma');


acc_by_config = grouped_accuracy( ...
    trials.true_idx, trials.correct, 'true_idx');
acc_by_config.true_name = eval_names(acc_by_config.true_idx);
acc_by_config = movevars( ...
    acc_by_config, 'true_name', 'After', 'true_idx');

acc_by_object_type = grouped_accuracy( ...
    object_types, trials.correct, 'object_type');

acc_by_config_kind = grouped_accuracy( ...
    config_kinds, trials.correct, 'config_kind');

pred_counts = prediction_counts(trials, library_names);

fprintf('\nAccuracy by sigma:\n')
disp(acc_by_sigma)

fprintf('\nAccuracy by true configuration:\n')
disp(acc_by_config)

fprintf('\nAccuracy by object type:\n')
disp(acc_by_object_type)

fprintf('\nAccuracy by nominal/deformed status:\n')
disp(acc_by_config_kind)

fprintf('\nPrediction counts:\n')
disp(pred_counts)

if isempty(sigma_to_report)
    mask_report = true(num_trials, 1);
else
    mask_report = trials.sigma == sigma_to_report;

    if ~any(mask_report)
        error('Requested sigma_to_report = %g was not found.', ...
            sigma_to_report);
    end
end

fprintf('\nConfusion table')
if ~isempty(sigma_to_report)
    fprintf(' for sigma = %g', sigma_to_report)
end
fprintf(':\n')

disp(confusion_counts( ...
    trials(mask_report, :), eval_names, library_names))

fprintf('\nClass-level accuracy by sigma for deformed targets:\n')
disp(class_accuracy_by_sigma( ...
    trials, object_types, config_kinds, library_names))

% ---- paper plots ----

if save_paper_plots && ~exist(paper_plot_folder, 'dir')
    mkdir(paper_plot_folder)
end

if make_paper_geometry_plot
    plot_paper_geometry_examples( ...
        data_folder_override, ...
        paper_example_realization, ...
        fig_paper_geometry, ...
        dock_figures, ...
        font_name, ...
        tick_font_size, ...
        label_font_size, ...
        title_font_size, ...
        legend_font_size, ...
        axes_line_width, ...
        line_width, ...
        save_paper_plots, ...
        paper_plot_folder)
end

if make_paper_accuracy_plot
    plot_paper_accuracy_by_sigma( ...
        trials, ...
        config_kinds, ...
        num_modes, ...
        amplitudes, ...
        fig_paper_accuracy, ...
        dock_figures, ...
        font_name, ...
        tick_font_size, ...
        label_font_size, ...
        title_font_size, ...
        legend_font_size, ...
        axes_line_width, ...
        line_width, ...
        marker_size, ...
        save_paper_plots, ...
        paper_plot_folder)
end

if make_paper_mismatch_plot
    mismatch_table = build_geometry_mismatch_table( ...
        eval_objects, ...
        eval_names, ...
        library_objects, ...
        library_names, ...
        data_folder_override);

    plot_paper_accuracy_by_mismatch( ...
        trials, ...
        mismatch_table, ...
        object_types, ...
        config_kinds, ...
        mismatch_sigma, ...
        fig_paper_mismatch, ...
        dock_figures, ...
        font_name, ...
        tick_font_size, ...
        label_font_size, ...
        title_font_size, ...
        legend_font_size, ...
        axes_line_width, ...
        marker_size, ...
        save_paper_plots, ...
        paper_plot_folder)
end

if make_paper_pole_plot
    plot_paper_pole_examples( ...
        data_folder_override, ...
        paper_example_realization, ...
        fig_paper_poles, ...
        dock_figures, ...
        font_name, ...
        tick_font_size, ...
        label_font_size, ...
        title_font_size, ...
        legend_font_size, ...
        axes_line_width, ...
        line_width, ...
        pole_marker_size, ...
        nominal_pole_marker_size, ...
        save_paper_plots, ...
        paper_plot_folder)
end

% ---- supplementary plots ----

if make_accuracy_sigma_plot
    plot_accuracy_by_sigma( ...
        acc_by_sigma, ...
        fig_accuracy_sigma, ...
        dock_figures, ...
        summary_figure_position, ...
        font_name, ...
        axes_line_width)
end

if make_accuracy_snr_plot
    plot_accuracy_distribution_by_snr( ...
        trials, ...
        config_kinds, ...
        num_modes, ...
        amplitudes, ...
        realization_ids, ...
        snr_by_sigma, ...
        snr_axis_label, ...
        fig_accuracy_snr, ...
        dock_figures, ...
        snr_figure_position, ...
        font_name, ...
        tick_font_size, ...
        label_font_size, ...
        title_font_size, ...
        axes_line_width, ...
        line_width, ...
        marker_size, ...
        grid_state)
end

if make_deformation_comparison_plots
    plot_deformation_comparisons( ...
        data_folder_override, ...
        deformation_comparison_plot_folder, ...
        save_deformation_comparison_plots);
end

if make_pole_comparison_plots
    plot_pole_comparisons( ...
        data_folder_override, ...
        pole_comparison_plot_folder, ...
        save_pole_comparison_plots, ...
        font_name, ...
        tick_font_size, ...
        label_font_size, ...
        title_font_size, ...
        legend_font_size, ...
        axes_line_width, ...
        line_width, ...
        pole_marker_size, ...
        nominal_pole_marker_size, ...
        deformed_pole_color, ...
        nominal_pole_color, ...
        grid_state);
end

if make_config_plots && any([ ...
        make_shape_panel, ...
        make_pole_panel, ...
        make_accuracy_text_panel])

    true_ids = unique(trials.true_idx, 'stable');

    if ~isempty(configs_to_plot)
        true_ids = true_ids(ismember(true_ids, configs_to_plot));
    end

    if save_config_plots && ~exist(config_plot_folder, 'dir')
        mkdir(config_plot_folder)
    end

    npanels = nnz([ ...
        make_shape_panel, ...
        make_pole_panel, ...
        make_accuracy_text_panel]);

    fprintf('\nPlotting %d configuration diagnostic figure(s).\n', ...
        numel(true_ids))

    for ii = 1:numel(true_ids)
        true_id = true_ids(ii);
        config_name = eval_names(true_id);

        file_path = resolve_data_file( ...
            eval_objects, ...
            true_id, ...
            data_folder_override, ...
            {'eval_file_path', 'file_path'}, ...
            {'eval_file_name', 'file_name'});

        D = load_diagnostics(file_path);
        [nominal_D, nominal_name] = load_nominal_reference( ...
            config_name, ...
            library_objects, ...
            library_names, ...
            data_folder_override);

        if strcmp(config_figure_mode, 'separate')
            fig_id = fig_config_start + ii - 1;
        elseif strcmp(config_figure_mode, 'reuse')
            fig_id = fig_config_start;
        else
            error(['config_figure_mode must be ', ...
                   '''separate'' or ''reuse''.'])
        end

        fh = newfig( ...
            fig_id, char(config_name), ...
            dock_figures, config_figure_position);

        tiledlayout( ...
            1, npanels, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact')

        if make_shape_panel
            nexttile
            hold on

            if isempty(D.curveX) || isempty(D.curveY)
                text(0.5, 0.5, 'curveX/curveY missing', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', text_font_size)
                axis off
            else
                if ~isempty(nominal_D.curveX)
                    plot( ...
                        nominal_D.curveX, nominal_D.curveY, '-', ...
                        'Color', nominal_shape_color, ...
                        'LineWidth', line_width, ...
                        'DisplayName', 'nominal')
                end

                plot( ...
                    D.curveX, D.curveY, '.', ...
                    'Color', deformed_shape_color, ...
                    'MarkerSize', shape_marker_size, ...
                    'DisplayName', 'configuration')

                axis equal
                style_axes( ...
                    gca, font_name, tick_font_size, ...
                    axes_line_width, grid_state)
                xlabel('x', 'FontSize', label_font_size)
                ylabel('y', 'FontSize', label_font_size)
                legend( ...
                    'Location', 'best', ...
                    'Interpreter', 'none', ...
                    'FontSize', legend_font_size)
            end

            hold off
            title( ...
                'nominal versus deformed shape', ...
                'Interpreter', 'none', ...
                'FontSize', title_font_size)
        end

        if make_pole_panel
            nexttile

            if isempty(D.poles) && isempty(nominal_D.poles)
                text(0.5, 0.5, 'poles missing', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', text_font_size)
                axis off
            else
                hold on

                if ~isempty(D.poles)
                    plot( ...
                        real(D.poles), imag(D.poles), 'o', ...
                        'Color', deformed_pole_color, ...
                        'MarkerSize', pole_marker_size, ...
                        'LineWidth', line_width, ...
                        'DisplayName', 'physical-geometry poles')
                end

                if ~isempty(nominal_D.poles)
                    plot( ...
                        real(nominal_D.poles), ...
                        imag(nominal_D.poles), ...
                        'x', ...
                        'Color', nominal_pole_color, ...
                        'MarkerSize', nominal_pole_marker_size, ...
                        'LineWidth', line_width, ...
                        'DisplayName', 'nominal library poles')
                end

                style_axes( ...
                    gca, font_name, tick_font_size, ...
                    axes_line_width, grid_state)
                xlabel('Re pole', 'FontSize', label_font_size)
                ylabel('Im pole', 'FontSize', label_font_size)

                xl = xlim;
                yl = ylim;
                plot(xl, [0, 0], 'k:', 'HandleVisibility', 'off')
                plot([0, 0], yl, 'k:', 'HandleVisibility', 'off')
                xlim(xl)
                ylim(yl)

                legend( ...
                    'Location', 'best', ...
                    'Interpreter', 'none', ...
                    'FontSize', legend_font_size)

                if strlength(nominal_name) > 0
                    subtitle( ...
                        "reference: " + nominal_name, ...
                        'Interpreter', 'none', ...
                        'FontSize', subtitle_font_size)
                end

                hold off
            end

            title( ...
                'physical versus nominal poles', ...
                'Interpreter', 'none', ...
                'FontSize', title_font_size)
        end

        if make_accuracy_text_panel
            nexttile

            mask_config = trials.true_idx == true_id;

            if isempty(sigmas_for_config_accuracy)
                mask_sigma = true(num_trials, 1);
                sigma_text = "all sigmas";
            else
                mask_sigma = ismember( ...
                    trials.sigma, sigmas_for_config_accuracy);
                sigma_text = "selected sigmas";
            end

            mask = mask_config & mask_sigma;

            if ~any(mask)
                text(0.5, 0.5, 'no matching trials', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', text_font_size)
            else
                metadata = deformation_metadata_from_name(config_name);
                ntrials = nnz(mask);
                ncorrect = nnz(trials.correct(mask));
                receiver_ids = unique( ...
                    trials.receiver_idx(mask), 'stable');
                sigmas_used = unique(trials.sigma(mask), 'stable');

                lines = strings(0, 1);
                lines(end + 1) = sprintf( ...
                    'accuracy: %.2f%%', 100*ncorrect/ntrials);
                lines(end + 1) = sprintf( ...
                    'correct: %d / %d', ncorrect, ntrials);
                lines(end + 1) = sprintf( ...
                    'receiver signals: %d', numel(receiver_ids));
                lines(end + 1) = "sigmas: " + sigma_text;
                lines(end + 1) = "";
                lines(end + 1) = "type: " + metadata.config_kind;
                lines(end + 1) = sprintf( ...
                    'modes: %d', metadata.num_modes);
                lines(end + 1) = sprintf( ...
                    'amplitude: %.4g', metadata.amplitude);
                lines(end + 1) = sprintf( ...
                    'realization: %d', metadata.realization_id);
                lines(end + 1) = "";

                for isig = 1:numel(sigmas_used)
                    sigma = sigmas_used(isig);
                    smask = mask & trials.sigma == sigma;

                    if sigma == 0
                        sigma_label = sprintf('%s = 0', sigma_char);
                    else
                        sigma_label = sprintf( ...
                            '%s = %.1e', sigma_char, sigma);
                    end

                    lines(end + 1) = sprintf( ...
                        '%s: %.2f%%  (%d/%d)', ...
                        sigma_label, ...
                        100*mean(trials.correct(smask)), ...
                        nnz(trials.correct(smask)), ...
                        nnz(smask));
                end

                text(0.05, 0.95, lines, ...
                    'Units', 'normalized', ...
                    'VerticalAlignment', 'top', ...
                    'FontName', 'Consolas', ...
                    'FontSize', text_font_size, ...
                    'Interpreter', interpreter)
            end

            axis off
            title( ...
                'classification accuracy', ...
                'Interpreter', 'none', ...
                'FontSize', title_font_size)
        end

        sgtitle( ...
            config_name, ...
            'Interpreter', 'none', ...
            'FontSize', title_font_size)
        drawnow

        if save_config_plots
            fname = matlab.lang.makeValidName(char(config_name));
            saveas(fh, fullfile(config_plot_folder, [fname, '.png']))
        end
    end
end

function T = class_accuracy_by_sigma( ...
    trials, object_types, config_kinds, library_names)
%CLASS_ACCURACY_BY_SIGMA  Accuracy for each true class and noise level.

    mask = config_kinds == "deformed";
    classes = unique(object_types(mask), 'stable');
    sigmas = unique(trials.sigma, 'stable');

    object_type = strings(0, 1);
    sigma = zeros(0, 1);
    accuracy = zeros(0, 1);
    ntrials = zeros(0, 1);

    for iclass = 1:numel(classes)
        for isig = 1:numel(sigmas)
            current = mask & ...
                object_types == classes(iclass) & ...
                trials.sigma == sigmas(isig);

            object_type(end + 1, 1) = classes(iclass); %#ok<AGROW>
            sigma(end + 1, 1) = sigmas(isig); %#ok<AGROW>
            accuracy(end + 1, 1) = mean(trials.correct(current)); %#ok<AGROW>
            ntrials(end + 1, 1) = nnz(current); %#ok<AGROW>
        end
    end

    T = table(object_type, sigma, accuracy, ntrials);

    % Confirm that the prediction indices refer to the nominal library.
    if any(trials.pred_idx < 1 | trials.pred_idx > numel(library_names))
        warning('Some prediction indices do not index the nominal library.')
    end
end

function plot_paper_geometry_examples( ...
    input_folder, requested_realization, fig_id, dock_figures, ...
    font_name, tick_font_size, label_font_size, title_font_size, ...
    legend_font_size, axes_line_width, line_width, save_plot, output_folder)
%PLOT_PAPER_GEOMETRY_EXAMPLES
%   Show nominal targets and all tested amplitudes at the largest mode count.

    records = load_deformation_records(input_folder);
    deformed = ~[records.is_nominal];
    max_modes = max([records(deformed).num_modes]);
    realization = select_example_realization( ...
        records, max_modes, requested_realization);

    object_tags = unique(string({records.object_tag}), 'stable');
    amplitudes = unique([records(deformed & ...
        [records.num_modes] == max_modes).amplitude]);
    amplitudes = sort(amplitudes);

    [x_limits, y_limits] = geometry_example_limits( ...
        records, object_tags, amplitudes, max_modes, realization);

    fig = newfig(fig_id, 'Paper: deformation examples', ...
        dock_figures, [100, 100, 850, 460]);
    layout = tiledlayout(fig, numel(object_tags), numel(amplitudes) + 1, ...
        'TileSpacing', 'none', 'Padding', 'compact');

    nominal_color = [0.15, 0.15, 0.15];
    deformed_color = [0.00, 0.45, 0.74];

    for itag = 1:numel(object_tags)
        tag = object_tags(itag);
        nominal_idx = find(string({records.object_tag}) == tag & ...
            [records.is_nominal], 1);
        nominal = records(nominal_idx);
        [xn, yn] = order_open_boundary(nominal.curveX, nominal.curveY);

        ax = nexttile(layout);
        plot(ax, xn, yn, '-', 'Color', nominal_color, ...
            'LineWidth', line_width)
        format_geometry_axis( ...
            ax, x_limits, y_limits, ...
            itag == numel(object_tags), true, ...
            font_name, tick_font_size, label_font_size, axes_line_width)
        title(ax, 'nominal', 'FontSize', title_font_size)
        text(ax, -0.42, 0.5, object_display_name(tag), ...
            'Units', 'normalized', ...
            'Rotation', 90, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontName', font_name, ...
            'FontSize', label_font_size, ...
            'FontWeight', 'bold', ...
            'Clipping', 'off')

        for iamp = 1:numel(amplitudes)
            idx = find(string({records.object_tag}) == tag & ...
                ~[records.is_nominal] & ...
                [records.num_modes] == max_modes & ...
                [records.amplitude] == amplitudes(iamp) & ...
                [records.realization_id] == realization, 1);

            if isempty(idx)
                error(['No geometry found for object %s, M = %d, ', ...
                    'amplitude = %g, realization = %d.'], ...
                    tag, max_modes, amplitudes(iamp), realization)
            end

            [xd, yd] = order_open_boundary( ...
                records(idx).curveX, records(idx).curveY);
            ax = nexttile(layout);
            hold(ax, 'on')
            plot(ax, xn, yn, '--', 'Color', nominal_color, ...
                'LineWidth', line_width, 'DisplayName', 'nominal')
            plot(ax, xd, yd, '-', 'Color', deformed_color, ...
                'LineWidth', line_width, 'DisplayName', 'deformed')
            hold(ax, 'off')
            format_geometry_axis( ...
                ax, x_limits, y_limits, ...
                itag == numel(object_tags), false, ...
                font_name, tick_font_size, label_font_size, axes_line_width)
            title(ax, sprintf('\\epsilon = %.3g', amplitudes(iamp)), ...
                'FontSize', title_font_size)

            if itag == 1 && iamp == numel(amplitudes)
                legend(ax, 'Location', 'best', ...
                    'FontSize', legend_font_size)
            end
        end
    end

    % title(layout, sprintf( ...
    %     'Representative deformations: M = %d, realization = %d', ...
    %     max_modes, realization), 'FontSize', title_font_size)

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'deformation_geometry_examples.pdf'), 'ContentType', 'vector')
    end
end

function plot_paper_accuracy_by_sigma( ...
    trials, config_kinds, num_modes, amplitudes, fig_id, ...
    dock_figures, font_name, tick_font_size, label_font_size, ...
    title_font_size, legend_font_size, axes_line_width, line_width, ...
    marker_size, save_plot, output_folder)
%PLOT_PAPER_ACCURACY_BY_SIGMA
%   Pooled classification accuracy for each fixed Fourier mode count,
%   deformation amplitude, and noise level.

    deformed_mask = config_kinds == "deformed";
    nominal_mask = config_kinds == "nominal";
    modes = sort(unique(num_modes(deformed_mask)));
    sigmas = sort(unique(trials.sigma));
    [xplot, ~, xtick_values, xtick_labels] = ...
        sigma_plot_data(sigmas, zeros(size(sigmas)));

    fig = newfig(fig_id, 'Paper: deformation accuracy', ...
        dock_figures, [100, 100, 1050, 390]);
    layout = tiledlayout(fig, 1, numel(modes), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    all_amplitudes = sort(unique(amplitudes(deformed_mask)));
    colors = lines(numel(all_amplitudes));

    for imode = 1:numel(modes)
        ax = nexttile(layout);
        hold(ax, 'on')

        baseline = nan(numel(sigmas), 1);
        for isig = 1:numel(sigmas)
            mask = nominal_mask & trials.sigma == sigmas(isig);
            baseline(isig) = mean(trials.correct(mask));
        end

        plot(ax, xplot, baseline, 'k--', ...
            'LineWidth', line_width, 'DisplayName', 'nominal baseline')

        mode_mask = deformed_mask & num_modes == modes(imode);
        mode_amplitudes = sort(unique(amplitudes(mode_mask)));

        for iamp = 1:numel(mode_amplitudes)
            amp = mode_amplitudes(iamp);
            color_idx = find(all_amplitudes == amp, 1);
            accuracy = nan(numel(sigmas), 1);

            for isig = 1:numel(sigmas)
                mask = mode_mask & amplitudes == amp & ...
                    trials.sigma == sigmas(isig);
                accuracy(isig) = mean(trials.correct(mask));
            end

            color = colors(color_idx, :);
            plot(ax, xplot, accuracy, '-o', 'Color', color, ...
                'LineWidth', line_width, 'MarkerSize', marker_size, ...
                'MarkerFaceColor', 'w', ...
                'DisplayName', sprintf('\\epsilon = %.3g', amp))
        end

        set(ax, 'XScale', 'log')
        xlim(ax, [min(xplot), max(xplot)])
        ylim(ax, [0, 1])
        xticks(ax, xtick_values)
        xticklabels(ax, cellstr(xtick_labels))
        style_axes(ax, font_name, tick_font_size, ...
            axes_line_width, 'on')
        xlabel(ax, 'noise level $\sigma$', ...
            'Interpreter', 'latex', 'FontSize', label_font_size)
        title(ax, sprintf('M = %d', modes(imode)), ...
            'FontSize', title_font_size)

        if imode == 1
            ylabel(ax, 'classification accuracy', ...
                'FontSize', label_font_size)
        end

        if imode == numel(modes)
            legend(ax, 'Location', 'best', ...
                'FontSize', legend_font_size)
        end
        hold(ax, 'off')
    end

    % title(layout, 'Deformed-target classification accuracy', ...
    %     'FontSize', title_font_size)

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'deformation_accuracy_by_noise.pdf'), 'ContentType', 'vector')
    end
end

function T = build_geometry_mismatch_table( ...
    eval_objects, eval_names, library_objects, library_names, ...
    data_folder_override)
%BUILD_GEOMETRY_MISMATCH_TABLE
%   Compute a normalized discrete Hausdorff distance for every deformation.

    neval = height(eval_objects);
    eval_object_type = strings(neval, 1);
    eval_num_modes = zeros(neval, 1);
    eval_amplitudes = zeros(neval, 1);
    eval_realization_ids = zeros(neval, 1);
    is_deformed = false(neval, 1);

    for ii = 1:neval
        metadata = deformation_metadata_from_name(eval_names(ii));
        eval_object_type(ii) = metadata.object_type;
        eval_num_modes(ii) = metadata.num_modes;
        eval_amplitudes(ii) = metadata.amplitude;
        eval_realization_ids(ii) = metadata.realization_id;
        is_deformed(ii) = metadata.config_kind == "deformed";
    end

    ids = find(is_deformed);
    normalized_hausdorff = nan(numel(ids), 1);

    for ii = 1:numel(ids)
        idx = ids(ii);
        file_path = resolve_data_file(eval_objects, idx, ...
            data_folder_override, {'eval_file_path', 'file_path'}, ...
            {'eval_file_name', 'file_name'});
        D = load_diagnostics(file_path);
        [nominal_D, ~] = load_nominal_reference(eval_names(idx), ...
            library_objects, library_names, data_folder_override);

        normalized_hausdorff(ii) = normalized_discrete_hausdorff( ...
            D.curveX, D.curveY, nominal_D.curveX, nominal_D.curveY);
    end

    true_idx = ids(:);
    true_name = eval_names(ids);
    object_type = eval_object_type(ids);
    num_modes_out = eval_num_modes(ids);
    amplitude = eval_amplitudes(ids);
    realization_id = eval_realization_ids(ids);
    T = table(true_idx, true_name, object_type, num_modes_out, ...
        amplitude, realization_id, normalized_hausdorff, ...
        'VariableNames', {'true_idx', 'true_name', 'object_type', ...
        'num_modes', 'amplitude', 'realization_id', ...
        'normalized_hausdorff'});
end

function plot_paper_accuracy_by_mismatch( ...
    trials, mismatch, object_types, config_kinds, requested_sigma, ...
    fig_id, dock_figures, ...
    font_name, tick_font_size, label_font_size, title_font_size, ...
    legend_font_size, axes_line_width, marker_size, save_plot, output_folder)
%PLOT_PAPER_ACCURACY_BY_MISMATCH
%   Isolate geometric model mismatch at zero or the lowest available noise.

    sigmas = sort(unique(trials.sigma));
    if isempty(requested_sigma)
        if any(sigmas == 0)
            sigma = 0;
        else
            sigma = sigmas(1);
        end
    else
        sigma = requested_sigma;
        if ~any(sigmas == sigma)
            error('mismatch_sigma = %g is not present in the results.', sigma)
        end
    end

    classes = unique(mismatch.object_type, 'stable');
    modes = sort(unique(mismatch.num_modes));
    colors = lines(numel(modes));
    markers = {'o', 's', '^', 'd', 'v', '>'};

    fig = newfig(fig_id, 'Paper: accuracy by geometry mismatch', ...
        dock_figures, [100, 100, 760, 370]);
    layout = tiledlayout(fig, 1, numel(classes), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for iclass = 1:numel(classes)
        ax = nexttile(layout);
        hold(ax, 'on')
        class_rows = mismatch.object_type == classes(iclass);

        for imode = 1:numel(modes)
            rows = find(class_rows & mismatch.num_modes == modes(imode));
            accuracy = nan(numel(rows), 1);

            for ii = 1:numel(rows)
                mask = trials.true_idx == mismatch.true_idx(rows(ii)) & ...
                    trials.sigma == sigma;
                accuracy(ii) = mean(trials.correct(mask));
            end

            scatter(ax, mismatch.normalized_hausdorff(rows), accuracy, ...
                marker_size^2, colors(imode, :), markers{imode}, ...
                'filled', 'MarkerFaceAlpha', 0.65, ...
                'DisplayName', sprintf('M = %d', modes(imode)))
        end

        nominal_mask = config_kinds == "nominal" & ...
            object_types == classes(iclass) & trials.sigma == sigma;
        nominal_accuracy = mean(trials.correct(nominal_mask));
        yline(ax, nominal_accuracy, 'k--', 'LineWidth', 1.25, ...
            'DisplayName', 'nominal baseline')

        ylim(ax, [0, 1])
        style_axes(ax, font_name, tick_font_size, ...
            axes_line_width, 'on')
        xlabel(ax, 'normalized discrete Hausdorff distance', ...
            'FontSize', label_font_size)
        title(ax, object_display_name(classes(iclass)), ...
            'FontSize', title_font_size)

        if iclass == 1
            ylabel(ax, 'classification accuracy', ...
                'FontSize', label_font_size)
        else
            legend(ax, 'Location', 'best', ...
                'FontSize', legend_font_size)
        end
        hold(ax, 'off')
    end

    title(layout, sprintf('Geometric model mismatch at \\sigma = %g', sigma), ...
        'FontSize', title_font_size)

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'deformation_accuracy_by_geometry_mismatch.pdf'), ...
            'ContentType', 'vector')
    end
end

function plot_paper_pole_examples( ...
    input_folder, requested_realization, fig_id, dock_figures, ...
    font_name, tick_font_size, label_font_size, title_font_size, ...
    legend_font_size, axes_line_width, line_width, pole_marker_size, ...
    nominal_pole_marker_size, save_plot, output_folder)
%PLOT_PAPER_POLE_EXAMPLES
%   Compare nominal poles with the mildest and strongest tested deformation.

    records = load_deformation_records(input_folder);
    deformed = ~[records.is_nominal];
    max_modes = max([records(deformed).num_modes]);
    realization = select_example_realization( ...
        records, max_modes, requested_realization);
    object_tags = unique(string({records.object_tag}), 'stable');
    amplitudes = sort(unique([records(deformed & ...
        [records.num_modes] == max_modes).amplitude]));
    amplitudes = unique([amplitudes(1), amplitudes(end)]);

    fig = newfig(fig_id, 'Paper: representative pole changes', ...
        dock_figures, [100, 100, 760, 600]);
    layout = tiledlayout(fig, numel(object_tags), numel(amplitudes), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for itag = 1:numel(object_tags)
        tag = object_tags(itag);
        nominal_idx = find(string({records.object_tag}) == tag & ...
            [records.is_nominal], 1);

        for iamp = 1:numel(amplitudes)
            idx = find(string({records.object_tag}) == tag & ...
                ~[records.is_nominal] & ...
                [records.num_modes] == max_modes & ...
                [records.amplitude] == amplitudes(iamp) & ...
                [records.realization_id] == realization, 1);
            ax = nexttile(layout);
            hold(ax, 'on')
            plot(ax, real(records(nominal_idx).poles), ...
                imag(records(nominal_idx).poles), 'x', ...
                'Color', [0.85, 0.325, 0.098], ...
                'MarkerSize', nominal_pole_marker_size, ...
                'LineWidth', line_width, 'DisplayName', 'nominal')
            plot(ax, real(records(idx).poles), imag(records(idx).poles), 'o', ...
                'Color', [0.00, 0.45, 0.74], ...
                'MarkerSize', pole_marker_size, 'LineWidth', line_width, ...
                'DisplayName', 'deformed')
            style_axes(ax, font_name, tick_font_size, ...
                axes_line_width, 'on')
            xlabel(ax, 'Re pole', 'FontSize', label_font_size)
            ylabel(ax, 'Im pole', 'FontSize', label_font_size)
            title(ax, sprintf('%s, \\epsilon = %.3g', ...
                object_display_name(tag), amplitudes(iamp)), ...
                'FontSize', title_font_size)

            if itag == 1 && iamp == numel(amplitudes)
                legend(ax, 'Location', 'best', ...
                    'FontSize', legend_font_size)
            end
            hold(ax, 'off')
        end
    end

    title(layout, sprintf('Pole perturbations: M = %d, realization = %d', ...
        max_modes, realization), 'FontSize', title_font_size)

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'deformation_representative_poles.pdf'), 'ContentType', 'vector')
    end
end

function records = load_deformation_records(input_folder)
%LOAD_DEFORMATION_RECORDS  Load geometry, pole, and deformation metadata.

    files = dir(fullfile(input_folder, '*.mat'));
    records = struct('name', {}, 'object_tag', {}, 'is_nominal', {}, ...
        'num_modes', {}, 'amplitude', {}, 'realization_id', {}, ...
        'curveX', {}, 'curveY', {}, 'poles', {});

    for ifile = 1:numel(files)
        file_path = fullfile(files(ifile).folder, files(ifile).name);
        S = load(file_path, 'config_name', 'config_info', ...
            'curveX', 'curveY', 'pols');
        required = {'config_name', 'config_info', 'curveX', 'curveY', 'pols'};

        if ~all(isfield(S, required))
            continue
        end

        info = S.config_info;
        required_info = {'object_tag', 'is_nominal', 'num_modes', ...
            'deformation_amplitude', 'realization_id'};
        if ~all(isfield(info, required_info))
            continue
        end

        record.name = string(S.config_name);
        record.object_tag = string(info.object_tag);
        record.is_nominal = logical(info.is_nominal);
        record.num_modes = info.num_modes;
        record.amplitude = info.deformation_amplitude;
        record.realization_id = info.realization_id;
        record.curveX = S.curveX(:);
        record.curveY = S.curveY(:);
        record.poles = S.pols(:);
        records(end + 1) = record; %#ok<AGROW>
    end

    if isempty(records)
        error('No complete deformation records were found in %s.', input_folder)
    end
end

function realization = select_example_realization( ...
    records, num_modes, requested_realization)
%SELECT_EXAMPLE_REALIZATION  Find a realization present for every example.

    deformed = ~[records.is_nominal] & [records.num_modes] == num_modes;
    tags = unique(string({records(deformed).object_tag}), 'stable');
    amplitudes = unique([records(deformed).amplitude]);
    candidates = unique([records(deformed).realization_id]);

    valid = false(size(candidates));
    for ii = 1:numel(candidates)
        valid(ii) = true;
        for itag = 1:numel(tags)
            for iamp = 1:numel(amplitudes)
                present = any(string({records.object_tag}) == tags(itag) & ...
                    ~[records.is_nominal] & ...
                    [records.num_modes] == num_modes & ...
                    [records.amplitude] == amplitudes(iamp) & ...
                    [records.realization_id] == candidates(ii));
                valid(ii) = valid(ii) && present;
            end
        end
    end
    candidates = sort(candidates(valid));

    if isempty(candidates)
        error('No realization is common to all requested paper examples.')
    end

    if isempty(requested_realization)
        realization = candidates(1);
    elseif ismember(requested_realization, candidates)
        realization = requested_realization;
    else
        error('Requested paper realization %d is not common to all examples.', ...
            requested_realization)
    end
end

function [x_limits, y_limits] = geometry_example_limits( ...
    records, object_tags, amplitudes, num_modes, realization)
%GEOMETRY_EXAMPLE_LIMITS  Common square limits for all example panels.

    x_all = zeros(0, 1);
    y_all = zeros(0, 1);

    for itag = 1:numel(object_tags)
        tag = object_tags(itag);
        nominal_idx = find(string({records.object_tag}) == tag & ...
            [records.is_nominal], 1);
        x_all = [x_all; records(nominal_idx).curveX(:)]; %#ok<AGROW>
        y_all = [y_all; records(nominal_idx).curveY(:)]; %#ok<AGROW>

        for iamp = 1:numel(amplitudes)
            idx = find(string({records.object_tag}) == tag & ...
                ~[records.is_nominal] & ...
                [records.num_modes] == num_modes & ...
                [records.amplitude] == amplitudes(iamp) & ...
                [records.realization_id] == realization, 1);
            x_all = [x_all; records(idx).curveX(:)]; %#ok<AGROW>
            y_all = [y_all; records(idx).curveY(:)]; %#ok<AGROW>
        end
    end

    valid = isfinite(x_all) & isfinite(y_all);
    x_all = x_all(valid);
    y_all = y_all(valid);

    x_center = 0.5*(min(x_all) + max(x_all));
    y_center = 0.5*(min(y_all) + max(y_all));
    half_width = 0.525*max( ...
        max(x_all) - min(x_all), ...
        max(y_all) - min(y_all));

    x_limits = x_center + [-half_width, half_width];
    y_limits = y_center + [-half_width, half_width];
end

function format_geometry_axis( ...
    ax, x_limits, y_limits, show_x_axis, show_y_axis, ...
    font_name, tick_font_size, label_font_size, axes_line_width)
%FORMAT_GEOMETRY_AXIS  Square, boxed geometry panel at a common scale.

    axis(ax, 'equal')
    pbaspect(ax, [1, 1, 1])
    xlim(ax, x_limits)
    ylim(ax, y_limits)
    box(ax, 'on')
    set(ax, ...
        'FontName', font_name, ...
        'FontSize', tick_font_size, ...
        'LineWidth', axes_line_width)

    if show_x_axis
        xlabel(ax, 'x', 'FontSize', label_font_size)
    else
        xticks(ax, [])
    end

    if show_y_axis
        ylabel(ax, 'y', 'FontSize', label_font_size)
    else
        yticks(ax, [])
    end
end

function value = normalized_discrete_hausdorff(x, y, xn, yn)
%NORMALIZED_DISCRETE_HAUSDORFF
%   Symmetric point-cloud Hausdorff distance divided by the diagonal of the
%   nominal geometry's bounding box. Both curves are first resampled at
%   uniform arc-length points to reduce sensitivity to their discretizations.

    A = resample_open_curve(x, y, 1200);
    B = resample_open_curve(xn, yn, 1200);

    scale = hypot(max(B(:, 1)) - min(B(:, 1)), ...
        max(B(:, 2)) - min(B(:, 2)));
    if scale == 0
        error('Nominal geometry has zero bounding-box diagonal.')
    end

    value = max(directed_cloud_distance(A, B), ...
        directed_cloud_distance(B, A))/scale;
end

function P = resample_open_curve(x, y, npoints)
%RESAMPLE_OPEN_CURVE  Uniform arc-length samples of an ordered open curve.

    [x, y] = order_open_boundary(x, y);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    s = [0; cumsum(hypot(diff(x), diff(y)))];
    [s, unique_idx] = unique(s, 'stable');
    x = x(unique_idx);
    y = y(unique_idx);

    if numel(s) < 2 || s(end) == 0
        error('Cannot resample a degenerate boundary.')
    end

    sq = linspace(0, s(end), npoints).';
    P = [interp1(s, x, sq, 'linear'), interp1(s, y, sq, 'linear')];
end

function value = directed_cloud_distance(A, B)
%DIRECTED_CLOUD_DISTANCE  Maximum nearest-neighbor distance from A to B.

    chunk_size = 200;
    nearest = inf(size(A, 1), 1);

    for first = 1:chunk_size:size(A, 1)
        last = min(first + chunk_size - 1, size(A, 1));
        dx = A(first:last, 1) - B(:, 1).';
        dy = A(first:last, 2) - B(:, 2).';
        nearest(first:last) = sqrt(min(dx.^2 + dy.^2, [], 2));
    end

    value = max(nearest);
end

function plot_pole_comparisons( ...
    input_folder, output_folder, save_plots, ...
    font_name, tick_font_size, label_font_size, ...
    title_font_size, legend_font_size, axes_line_width, ...
    line_width, pole_marker_size, nominal_pole_marker_size, ...
    deformed_pole_color, nominal_pole_color, grid_state)
%PLOT_POLE_COMPARISONS
%   Overlay each deformed object's poles with the poles of its matching
%   unperturbed object. One figure is created per object class, with one
%   tile per deformation configuration.

    files = dir(fullfile(input_folder, '*.mat'));

    if isempty(files)
        warning('No .mat files found for pole plotting.');
        return
    end

    records = struct( ...
        'name', {}, ...
        'object_tag', {}, ...
        'is_nominal', {}, ...
        'num_modes', {}, ...
        'amplitude', {}, ...
        'realization_id', {}, ...
        'poles', {});

    for ifile = 1:numel(files)
        file_path = fullfile(files(ifile).folder, files(ifile).name);
        loaded = load(file_path, 'config_name', 'config_info', 'pols');

        if ~all(isfield(loaded, {'config_name', 'config_info', 'pols'}))
            warning( ...
                'Skipping "%s" because pole metadata is missing.', ...
                files(ifile).name);
            continue
        end

        info = loaded.config_info;

        required_info = { ...
            'object_tag', 'is_nominal', 'num_modes', ...
            'deformation_amplitude', 'realization_id'};

        if ~all(isfield(info, required_info))
            warning( ...
                'Skipping "%s" because config_info is incomplete.', ...
                files(ifile).name);
            continue
        end

        record.name = string(loaded.config_name);
        record.object_tag = string(info.object_tag);
        record.is_nominal = logical(info.is_nominal);
        record.num_modes = info.num_modes;
        record.amplitude = info.deformation_amplitude;
        record.realization_id = info.realization_id;
        record.poles = loaded.pols(:);

        records(end + 1) = record; %#ok<AGROW>
    end

    if isempty(records)
        warning('No valid pole records were loaded.');
        return
    end

    if save_plots && ~exist(output_folder, 'dir')
        mkdir(output_folder)
    end

    object_tags = unique(string({records.object_tag}), 'stable');

    for itag = 1:numel(object_tags)
        object_tag = object_tags(itag);
        same_object = string({records.object_tag}) == object_tag;
        nominal_indices = find(same_object & [records.is_nominal]);
        deformed_indices = find(same_object & ~[records.is_nominal]);

        if numel(nominal_indices) ~= 1
            warning( ...
                'Expected one nominal pole set for "%s"; found %d.', ...
                object_tag, numel(nominal_indices));
            continue
        end

        if isempty(deformed_indices)
            warning('No deformed pole sets found for "%s".', object_tag);
            continue
        end

        sort_data = [ ...
            [records(deformed_indices).num_modes].', ...
            -[records(deformed_indices).amplitude].', ...
            [records(deformed_indices).realization_id].'];
        [~, order] = sortrows(sort_data, [1, 2, 3]);
        deformed_indices = deformed_indices(order);

        nplots = numel(deformed_indices);
        ncols = min(3, nplots);
        nrows = ceil(nplots/ncols);

        fig = figure( ...
            'Color', 'w', ...
            'Name', sprintf('%s pole comparison', object_tag));
        layout = tiledlayout( ...
            fig, nrows, ncols, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        nominal = records(nominal_indices);

        for iplot = 1:nplots
            deformed = records(deformed_indices(iplot));
            ax = nexttile(layout);
            hold(ax, 'on')

            plot( ...
                ax, real(nominal.poles), imag(nominal.poles), 'x', ...
                'Color', nominal_pole_color, ...
                'MarkerSize', nominal_pole_marker_size, ...
                'LineWidth', line_width, ...
                'DisplayName', 'unperturbed poles');

            plot( ...
                ax, real(deformed.poles), imag(deformed.poles), 'o', ...
                'Color', deformed_pole_color, ...
                'MarkerSize', pole_marker_size, ...
                'LineWidth', line_width, ...
                'MarkerFaceColor', 'none', ...
                'DisplayName', 'deformed poles');

            style_axes( ...
                ax, font_name, tick_font_size, ...
                axes_line_width, grid_state)
            box(ax, 'on')
            xlabel(ax, 'Re pole', 'FontSize', label_font_size)
            ylabel(ax, 'Im pole', 'FontSize', label_font_size)

            xl = xlim(ax);
            yl = ylim(ax);
            plot(ax, xl, [0, 0], 'k:', 'HandleVisibility', 'off')
            plot(ax, [0, 0], yl, 'k:', 'HandleVisibility', 'off')
            xlim(ax, xl)
            ylim(ax, yl)

            title( ...
                ax, sprintf('m = %d, amplitude = %.4g, realization = %d', ...
                deformed.num_modes, deformed.amplitude, ...
                deformed.realization_id), ...
                'Interpreter', 'none', ...
                'FontSize', title_font_size)

            if iplot == 1
                legend( ...
                    ax, 'Location', 'best', ...
                    'Interpreter', 'none', ...
                    'FontSize', legend_font_size)
            end

            hold(ax, 'off')
        end

        sgtitle( ...
            layout, sprintf('%s: deformed versus unperturbed poles', ...
            object_display_name(object_tag)), ...
            'Interpreter', 'none', ...
            'FontSize', title_font_size)

        if save_plots
            file_stem = matlab.lang.makeValidName( ...
                sprintf('%s_pole_comparison', object_tag));
            exportgraphics( ...
                fig, fullfile(output_folder, [file_stem, '.png']), ...
                'Resolution', 300)
        end
    end
end

function plot_deformation_comparisons( ...
    input_folder, output_folder, save_plots)
%PLOT_DEFORMATION_COMPARISONS
%   Compare every deformed boundary against its nominal boundary.
%
%   One figure is created for each object class. Each tile contains one
%   nominal/deformed pair.

    files = dir(fullfile(input_folder, '*.mat'));

    if isempty(files)
        warning('No .mat files found for geometry plotting.');
        return
    end

    records = struct( ...
        'name', {}, ...
        'object_tag', {}, ...
        'is_nominal', {}, ...
        'num_modes', {}, ...
        'amplitude', {}, ...
        'realization_id', {}, ...
        'curveX', {}, ...
        'curveY', {});

    for ifile = 1:numel(files)
        file_path = fullfile(files(ifile).folder, files(ifile).name);

        loaded = load( ...
            file_path, ...
            'config_name', ...
            'config_info', ...
            'curveX', ...
            'curveY');

        required_fields = { ...
            'config_name', ...
            'config_info', ...
            'curveX', ...
            'curveY'};

        if ~all(isfield(loaded, required_fields))
            warning( ...
                'Skipping "%s" because geometry information is missing.', ...
                files(ifile).name);
            continue
        end

        info = loaded.config_info;

        record.name = string(loaded.config_name);
        record.object_tag = string(info.object_tag);
        record.is_nominal = logical(info.is_nominal);
        record.num_modes = info.num_modes;
        record.amplitude = info.deformation_amplitude;
        record.realization_id = info.realization_id;
        record.curveX = loaded.curveX(:);
        record.curveY = loaded.curveY(:);

        records(end + 1) = record; %#ok<AGROW>
    end

    if isempty(records)
        warning('No valid geometry records were loaded.');
        return
    end

    object_tags = unique( ...
        string({records.object_tag}), ...
        'stable');

    nominal_color = [0.15, 0.15, 0.15];
    deformed_color = [0.00, 0.45, 0.74];

    for itag = 1:numel(object_tags)
        object_tag = object_tags(itag);

        same_object = string({records.object_tag}) == object_tag;
        nominal_mask = same_object & [records.is_nominal];
        deformed_mask = same_object & ~[records.is_nominal];

        nominal_index = find(nominal_mask);

        if numel(nominal_index) ~= 1
            warning( ...
                'Expected one nominal geometry for "%s"; found %d.', ...
                object_tag, numel(nominal_index));
            continue
        end

        deformed_indices = find(deformed_mask);

        if isempty(deformed_indices)
            warning('No deformations found for "%s".', object_tag);
            continue
        end

        % Sort first by number of modes and then by decreasing amplitude.
        sort_data = [ ...
            [records(deformed_indices).num_modes].', ...
            -[records(deformed_indices).amplitude].', ...
            [records(deformed_indices).realization_id].'];

        [~, order] = sortrows(sort_data, [1, 2, 3]);
        deformed_indices = deformed_indices(order);

        nplots = numel(deformed_indices);
        ncols = min(3, nplots);
        nrows = ceil(nplots/ncols);

        fig = figure( ...
            'Color', 'w', ...
            'Name', sprintf('%s deformation comparison', object_tag));

        layout = tiledlayout( ...
            fig, ...
            nrows, ...
            ncols, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        nominal = records(nominal_index);

        for iplot = 1:nplots
            deformed = records(deformed_indices(iplot));

            [nominal_x, nominal_y] = order_open_boundary( ...
                nominal.curveX, nominal.curveY);
            [deformed_x, deformed_y] = order_open_boundary( ...
                deformed.curveX, deformed.curveY);

            ax = nexttile(layout);
            hold(ax, 'on');

            plot( ...
                ax, ...
                nominal_x, ...
                nominal_y, ...
                '-', ...
                'Color', nominal_color, ...
                'LineWidth', 1.5, ...
                'DisplayName', 'nominal');

            plot( ...
                ax, ...
                deformed_x, ...
                deformed_y, ...
                '-', ...
                'Color', deformed_color, ...
                'LineWidth', 1.2, ...
                'DisplayName', 'deformed');

            axis(ax, 'equal');
            grid(ax, 'on');
            box(ax, 'on');

            title( ...
                ax, ...
                sprintf( ...
                    'm = %d, \\epsilon = %.3g, r = %d', ...
                    deformed.num_modes, ...
                    deformed.amplitude, ...
                    deformed.realization_id));

            if iplot == 1
                legend(ax, 'Location', 'best');
            end
        end

        title(layout, sprintf('%s: nominal versus deformed', object_tag));

        if save_plots
            if ~exist(output_folder, 'dir')
                mkdir(output_folder);
            end

            figure_path = fullfile( ...
                output_folder, ...
                sprintf('%s_geometry_comparison.png', object_tag));

            exportgraphics( ...
                fig, ...
                figure_path, ...
                'Resolution', 300);

            fprintf('Saved geometry comparison to:\n  %s\n', ...
                figure_path);
        end
    end
end

function [x_ordered, y_ordered] = order_open_boundary(x, y)
%ORDER_OPEN_BOUNDARY  Order star-shaped cavity nodes without closing its gap.
%
%   The saved global arrays concatenate nodes in patch-index order, which is
%   not necessarily geometric order. Sorting by polar angle reconstructs the
%   boundary order. The sequence is rotated after the largest angular gap so
%   MATLAB does not draw a line across the cavity opening.

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        x_ordered = x;
        y_ordered = y;
        return
    end

    theta = atan2(y, x);
    [theta, sort_idx] = sort(theta);

    x = x(sort_idx);
    y = y(sort_idx);

    angular_gaps = [ ...
        diff(theta); ...
        theta(1) + 2*pi - theta(end)];

    [~, gap_idx] = max(angular_gaps);
    start_idx = mod(gap_idx, numel(theta)) + 1;
    order = [start_idx:numel(theta), 1:start_idx - 1];

    x_ordered = x(order);
    y_ordered = y(order);
end


function [snr_table, axis_label] = extract_snr_by_sigma(results)
%EXTRACT_SNR_BY_SIGMA  Normalize the saved representative-SNR mapping.
%
%   The updated driver saves results.representative_snr. This function also
%   accepts several older result-field names and common SNR column names.

    result_fields = { ...
        'representative_snr', ...
        'snr_by_sigma', ...
        'snr_table', ...
        'snr'};

    raw = [];

    for ii = 1:numel(result_fields)
        if isfield(results, result_fields{ii})
            raw = results.(result_fields{ii});
            break
        end
    end

    if isempty(raw)
        error([ ...
            'No representative SNR mapping was found in the result file.\n', ...
            'Rerun the updated identification driver so it saves ', ...
            'results.representative_snr.'])
    end

    if istable(raw)
        T = raw;
    elseif isstruct(raw)
        T = struct2table(raw);
    else
        error('The saved representative SNR data must be a table or struct.')
    end

    vars = T.Properties.VariableNames;
    sigma_name = first_matching_variable(vars, { ...
        'sigma', 'noise_sigma', 'noise_level'});

    snr_candidates = { ...
        'snr_late_db', ...
        'late_time_snr_db', ...
        'snr_late_time_db', ...
        'late_snr_db', ...
        'snr_full_db', ...
        'snr_db', ...
        'snr_late', ...
        'late_time_snr', ...
        'snr_late_time', ...
        'late_snr', ...
        'snr_full', ...
        'snr'};
    snr_name = first_matching_variable(vars, snr_candidates);

    if strlength(sigma_name) == 0 || strlength(snr_name) == 0
        error([ ...
            'Could not identify sigma and SNR columns in the saved ', ...
            'representative-SNR data. Available columns: %s'], ...
            strjoin(vars, ', '))
    end

    sigma = T.(sigma_name);
    snr = T.(snr_name);

    sigma = sigma(:);
    snr = snr(:);

    if numel(sigma) ~= numel(snr)
        error('The representative sigma and SNR arrays have different sizes.')
    end

    [sigma, order] = sort(sigma);
    snr = snr(order);

    if numel(unique(sigma)) ~= numel(sigma)
        error('The representative-SNR table contains duplicate sigma values.')
    end

    snr_table = table(sigma, snr);

    if contains(lower(snr_name), 'db')
        axis_label = 'representative late-time SNR (dB)';
    else
        axis_label = 'representative late-time SNR';
    end
end

function name = first_matching_variable(vars, candidates)
%FIRST_MATCHING_VARIABLE  Return the first matching table variable name.

    name = "";

    for ii = 1:numel(candidates)
        idx = find(strcmpi(vars, candidates{ii}), 1);

        if ~isempty(idx)
            name = string(vars{idx});
            return
        end
    end
end

function plot_accuracy_distribution_by_snr( ...
    trials, config_kinds, num_modes, amplitudes, realization_ids, ...
    snr_by_sigma, snr_axis_label, fig_id, dock_figures, fig_position, ...
    font_name, tick_font_size, label_font_size, title_font_size, ...
    axes_line_width, line_width, marker_size, grid_state)
%PLOT_ACCURACY_DISTRIBUTION_BY_SNR
%   Plot deformation accuracy versus SNR in one panel per Fourier mode count.
%
%   Within each panel, one curve is shown for each deformation amplitude.
%   Accuracy is aggregated over all realizations and over both object classes.
%   This avoids hiding systematic class collapse while keeping the figure
%   readable.

    deformed_mask = config_kinds == "deformed";
    mode_values = unique(num_modes(deformed_mask));
    mode_values = mode_values(isfinite(mode_values));
    mode_values = sort(mode_values(:));

    if isempty(mode_values)
        warning('No deformed configurations were found for the SNR plot.')
        return
    end

    sigma_values = unique(trials.sigma, 'stable');
    nsigma = numel(sigma_values);
    snr_values = nan(nsigma, 1);

    for isig = 1:nsigma
        sigma = sigma_values(isig);
        snr_idx = find(snr_by_sigma.sigma == sigma, 1);

        if isempty(snr_idx)
            error('No representative SNR was saved for sigma = %g.', sigma)
        end

        snr_values(isig) = snr_by_sigma.snr(snr_idx);
    end

    finite_snr = isfinite(snr_values);
    omitted_infinite = any(~finite_snr);
    sigma_values = sigma_values(finite_snr);
    snr_values = snr_values(finite_snr);

    [snr_values, snr_order] = sort(snr_values);
    sigma_values = sigma_values(snr_order);

    newfig( ...
        fig_id, 'Accuracy by SNR and Fourier modes', ...
        dock_figures, fig_position);

    ncols = min(3, numel(mode_values));
    nrows = ceil(numel(mode_values)/ncols);
    tl = tiledlayout(nrows, ncols, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for imode = 1:numel(mode_values)
        mode_value = mode_values(imode);
        ax = nexttile(tl);
        hold(ax, 'on')

        mode_mask = deformed_mask & num_modes == mode_value;
        amplitude_values = unique(amplitudes(mode_mask));
        amplitude_values = amplitude_values(isfinite(amplitude_values));
        amplitude_values = sort(amplitude_values(:));

        for iamp = 1:numel(amplitude_values)
            amplitude_value = amplitude_values(iamp);
            accuracy = nan(numel(sigma_values), 1);

            for isig = 1:numel(sigma_values)
                mask = mode_mask & ...
                    amplitudes == amplitude_value & ...
                    trials.sigma == sigma_values(isig);

                if any(mask)
                    % Aggregates over every realization and both object classes.
                    accuracy(isig) = mean(trials.correct(mask));
                end
            end

            plot( ...
                ax, snr_values, accuracy, '-o', ...
                'LineWidth', line_width, ...
                'MarkerSize', marker_size, ...
                'DisplayName', sprintf('amplitude = %.3g', amplitude_value))
        end

        hold(ax, 'off')
        ylim(ax, [0, 1])
        style_axes( ...
            ax, font_name, tick_font_size, ...
            axes_line_width, grid_state)
        title( ...
            ax, sprintf('%d Fourier modes', mode_value), ...
            'FontSize', title_font_size, ...
            'Interpreter', 'none')

        if imode > (nrows - 1)*ncols
            xlabel(ax, snr_axis_label, 'FontSize', label_font_size)
        end

        if mod(imode - 1, ncols) == 0
            ylabel(ax, 'classification accuracy', ...
                'FontSize', label_font_size)
        end

        legend(ax, 'Location', 'best', 'Interpreter', 'none')
    end

    title( ...
        tl, 'Classification accuracy by deformation amplitude and Fourier modes', ...
        'FontSize', title_font_size, ...
        'Interpreter', 'none')
    drawnow

    if omitted_infinite
        fprintf([ ...
            '\nThe zero-noise point was omitted from the SNR plot ', ...
            'because its representative SNR is infinite.\n'])
    end
end

function plot_accuracy_by_sigma( ...
    acc, fig_id, dock_figures, fig_position, ...
    font_name, axes_line_width)
%PLOT_ACCURACY_BY_SIGMA  Plot accuracy across every evaluation trial.

    [xplot, yplot, tick_values, tick_labels] = ...
        sigma_plot_data(acc.sigma, acc.accuracy);

    newfig( ...
        fig_id, 'Accuracy by sigma', ...
        dock_figures, fig_position);

    semilogx( ...
        xplot, yplot, '-o', ...
        'LineWidth', 2, ...
        'MarkerSize', 7)

    grid on
    xlabel('noise level')
    ylabel('classification accuracy')
    title('Overall classification accuracy by noise level')

    xlim([min(xplot), max(xplot)])
    ylim([0, 1])
    pbaspect([1, 1, 1])
    xticks(tick_values)
    xticklabels(cellstr(tick_labels))

    set(gca, ...
        'FontSize', 15, ...
        'FontName', font_name, ...
        'LineWidth', axes_line_width)
    drawnow
end

function [xplot, yplot, tick_values, tick_labels] = ...
    sigma_plot_data(sigmas, accuracy)
%SIGMA_PLOT_DATA  Place sigma = 0 one log interval before the first sigma.

    sigmas = sigmas(:);
    accuracy = accuracy(:);

    zero_mask = sigmas == 0;
    pos_mask = sigmas > 0;

    sigmas_pos = sigmas(pos_mask);
    accuracy_pos = accuracy(pos_mask);
    [sigmas_pos, idx] = sort(sigmas_pos);
    accuracy_pos = accuracy_pos(idx);

    if any(zero_mask)
        if numel(sigmas_pos) < 2
            error('Need at least two positive sigma values to place sigma = 0.')
        end

        sigma_zero_plot = sigmas_pos(1)^2/sigmas_pos(2);
        xplot = [sigma_zero_plot; sigmas_pos];
        yplot = [accuracy(find(zero_mask, 1)); accuracy_pos];
        tick_values = xplot;
        tick_labels = ["0"; compose('%.0e', sigmas_pos)];
    else
        xplot = sigmas_pos;
        yplot = accuracy_pos;
        tick_values = xplot;
        tick_labels = compose('%.0e', sigmas_pos);
    end
end

function T = deformation_accuracy_table( ...
    object_types, num_modes, amplitudes, correct)
%DEFORMATION_ACCURACY_TABLE
%   Compute accuracy by object class, mode count, and amplitude.

    combinations = unique( ...
        table(object_types, num_modes, amplitudes), ...
        'rows', ...
        'stable');

    accuracy = zeros(height(combinations), 1);
    ntrials = zeros(height(combinations), 1);

    for ii = 1:height(combinations)
        mask = object_types == combinations.object_types(ii) & ...
            num_modes == combinations.num_modes(ii) & ...
            amplitudes == combinations.amplitudes(ii);

        accuracy(ii) = mean(correct(mask));
        ntrials(ii) = nnz(mask);
    end

    T = combinations;
    T.Properties.VariableNames = { ...
        'object_type', 'num_modes', 'amplitude'};
    T.accuracy = accuracy;
    T.ntrials = ntrials;
    T = sortrows(T, {'object_type', 'num_modes', 'amplitude'});
end

function metadata = deformation_metadata_from_name(name)
%DEFORMATION_METADATA_FROM_NAME  Parse nominal and deformation names.

    name = char(string(name));

    nominal_tokens = regexp( ...
        name, ...
        '^(.*)_nominal$', ...
        'tokens', ...
        'once');

    if ~isempty(nominal_tokens)
        metadata.object_type = string(nominal_tokens{1});
        metadata.config_kind = "nominal";
        metadata.num_modes = 0;
        metadata.realization_id = 0;
        metadata.amplitude = 0;
        return
    end

    deformation_tokens = regexp( ...
        name, ...
        '^(.*)_m(\d+)_r(\d+)_e(\d+)p(\d+)$', ...
        'tokens', ...
        'once');

    if isempty(deformation_tokens)
        error('Could not parse deformation configuration name "%s".', name)
    end

    whole_part = deformation_tokens{4};
    fractional_part = deformation_tokens{5};

    metadata.object_type = string(deformation_tokens{1});
    metadata.config_kind = "deformed";
    metadata.num_modes = str2double(deformation_tokens{2});
    metadata.realization_id = str2double(deformation_tokens{3});
    metadata.amplitude = str2double( ...
        [whole_part, '.', fractional_part]);
end

function label = object_display_name(object_type)
%OBJECT_DISPLAY_NAME  Convert generated object tags to display names.

    if object_type == "circ_o035"
        label = 'circular cavity';
    elseif object_type == "ell_o045"
        label = 'elliptic cavity';
    else
        label = char(object_type);
    end
end

function [D, nominal_name] = load_nominal_reference( ...
    config_name, library_objects, library_names, data_folder_override)
%LOAD_NOMINAL_REFERENCE  Load the matching nominal geometry and poles.

    target = deformation_metadata_from_name(config_name);

    for ii = 1:height(library_objects)
        candidate = deformation_metadata_from_name(library_names(ii));

        if candidate.config_kind ~= "nominal"
            continue
        end

        if candidate.object_type ~= target.object_type
            continue
        end

        file_path = resolve_data_file( ...
            library_objects, ...
            ii, ...
            data_folder_override, ...
            {'file_path'}, ...
            {'file_name'});

        D = load_diagnostics(file_path);
        nominal_name = library_names(ii);
        return
    end

    D = empty_diagnostics();
    nominal_name = "";
end

function fh = newfig(fig_id, fig_name, dock_figures, fig_position)
%NEWFIG  Create or clear one controlled figure.

    fh = figure(fig_id);
    clf(fh)
    set(fh, 'Name', fig_name, 'NumberTitle', 'on')

    if dock_figures
        set(fh, 'WindowStyle', 'docked')
    else
        set(fh, ...
            'WindowStyle', 'normal', ...
            'Position', fig_position)
    end
end

function style_axes( ...
    ax, font_name, tick_font_size, axes_line_width, grid_state)
%STYLE_AXES  Apply common axes styling.

    set(ax, ...
        'FontName', font_name, ...
        'FontSize', tick_font_size, ...
        'LineWidth', axes_line_width)
    grid(ax, grid_state)
end

function T = grouped_accuracy(keys, correct, key_name)
%GROUPED_ACCURACY  Compute accuracy grouped by a vector of keys.

    vals = unique(keys, 'stable');
    accuracy = zeros(numel(vals), 1);
    ntrials = zeros(numel(vals), 1);

    for ii = 1:numel(vals)
        if isnumeric(keys) && isnan(vals(ii))
            mask = isnan(keys);
        else
            mask = keys == vals(ii);
        end

        accuracy(ii) = mean(correct(mask));
        ntrials(ii) = nnz(mask);
    end

    T = table( ...
        vals(:), accuracy, ntrials, ...
        'VariableNames', {key_name, 'accuracy', 'ntrials'});
end

function counts = prediction_counts(trials, library_names)
%PREDICTION_COUNTS  Count predictions by nominal library object.

    pred_ids = unique(trials.pred_idx, 'stable');
    pred_name = library_names(pred_ids);
    count = zeros(numel(pred_ids), 1);

    for ii = 1:numel(pred_ids)
        count(ii) = nnz(trials.pred_idx == pred_ids(ii));
    end

    fraction = count/height(trials);
    counts = table( ...
        pred_ids, pred_name, count, fraction, ...
        'VariableNames', { ...
            'pred_idx', 'pred_name', 'count', 'fraction'});
end

function conf = confusion_counts(trials, eval_names, library_names)
%CONFUSION_COUNTS  Count true geometry versus predicted nominal class.

    true_ids = unique(trials.true_idx, 'stable');
    pred_ids = unique(trials.pred_idx, 'stable');
    counts = zeros(numel(true_ids), numel(pred_ids));

    for ii = 1:numel(true_ids)
        for jj = 1:numel(pred_ids)
            counts(ii, jj) = nnz( ...
                trials.true_idx == true_ids(ii) & ...
                trials.pred_idx == pred_ids(jj));
        end
    end

    pred_vars = matlab.lang.makeValidName( ...
        cellstr(library_names(pred_ids)));
    pred_vars = matlab.lang.makeUniqueStrings(pred_vars);

    conf = array2table( ...
        counts, ...
        'VariableNames', pred_vars, ...
        'RowNames', cellstr(eval_names(true_ids)));
end

function D = load_diagnostics(file_path)
%LOAD_DIAGNOSTICS  Load plotting data from one original .mat file.

    if ~isfile(file_path)
        error('Could not find original data file:\n  %s', file_path);
    end

    S = load( ...
        file_path, ...
        'curveX', ...
        'curveY', ...
        'pols', ...
        'receiver_angles');

    D.curveX = field_or_empty(S, 'curveX');
    D.curveY = field_or_empty(S, 'curveY');
    D.poles = field_or_empty(S, 'pols');
    D.receiver_angles = field_or_empty(S, 'receiver_angles');
    D.poles = D.poles(:);
    D.receiver_angles = D.receiver_angles(:);
end

function file_path = resolve_data_file( ...
    T, idx, folder_override, path_cols, name_cols)
%RESOLVE_DATA_FILE  Use saved path first, then folder_override/file_name.

    file_path = char(table_string(T, path_cols, idx));

    if isfile(file_path)
        return
    end

    file_name = table_string(T, name_cols, idx);
    candidate = fullfile(char(folder_override), char(file_name));

    if strlength(string(folder_override)) > 0 && isfile(candidate)
        file_path = candidate;
        return
    end

    error([ ...
        'Could not resolve original data file.\n', ...
        'Saved path:\n  %s\n', ...
        'Fallback path:\n  %s'], ...
        file_path, candidate);
end

function x = field_or_empty(S, field_name)
%FIELD_OR_EMPTY  Return a struct field or an empty array.

    if isfield(S, field_name)
        x = S.(field_name);
    else
        x = [];
    end
end

function D = empty_diagnostics()
%EMPTY_DIAGNOSTICS  Empty diagnostic data.

    D.curveX = [];
    D.curveY = [];
    D.poles = [];
    D.receiver_angles = [];
end

function vals = table_strings(T, candidates)
%TABLE_STRINGS  Return one string column from a table.

    vals = strings(height(T), 1);

    for ii = 1:height(T)
        vals(ii) = table_string(T, candidates, ii);
    end
end

function val = table_string(T, candidates, idx)
%TABLE_STRING  Return table value from the first available column.

    vars = T.Properties.VariableNames;

    for jj = 1:numel(candidates)
        if any(strcmp(vars, candidates{jj}))
            val = string(T.(candidates{jj})(idx));
            return
        end
    end

    val = "";
end
