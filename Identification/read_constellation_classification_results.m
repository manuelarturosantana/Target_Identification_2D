% read_constellation_classification_results.m
%
% Read and plot same-object constellation classification results produced by
% run_identification_single_poles_constellations.m.

clear
clc
close all

% ---- user-defined parameters ----

results_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'single_poles_constellation_results.mat');

data_folder_override = fullfile( ...
    'local', ...
    'data_same_object_constellations_2d_fast');

% ---- paper plots ----

make_paper_geometry_plot = true;
make_paper_accuracy_noise_plot = true;
make_paper_accuracy_count_plot = true;
make_paper_class_accuracy_plot = true;
make_paper_pole_plot = false;

save_paper_plots = true;
paper_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'constellation_paper_figures');

% Geometry and pole examples use the largest available object count.
% Empty selects sigma = 0 when it is available.
summary_sigma = [];

% ---- supplementary plots ----

make_supplementary_all_geometry_plot = false;
make_supplementary_all_pole_plot = false;
save_supplementary_plots = false;

supplementary_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'constellation_supplementary_figures');

dock_figures = true;

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
pole_marker_size = 9;
nominal_pole_marker_size = 11;

nominal_shape_color = [0.15, 0.15, 0.15];
constellation_shape_color = [0.00, 0.45, 0.74];
nominal_pole_color = [0.85, 0.325, 0.098];
constellation_pole_color = [0.00, 0.45, 0.74];

fig_paper_geometry = 1;
fig_paper_accuracy_noise = 2;
fig_paper_accuracy_count = 3;
fig_paper_class_accuracy = 4;
fig_paper_poles = 5;
fig_supplementary_geometry = 20;
fig_supplementary_poles = 30;

% ---------------------------------

if dock_figures
    set(groot, 'DefaultFigureWindowStyle', 'docked')
else
    set(groot, 'DefaultFigureWindowStyle', 'normal')
end

S = load(results_path, 'results');
results = S.results;
trials = results.trials;

if ~isfield(results, 'eval_metadata') || ...
        ~isfield(results, 'library_metadata')
    error([ ...
        'Constellation metadata is missing. Rerun identification with ', ...
        'run_identification_single_poles_constellations.m.'])
end

metadata = results.eval_metadata;
library_metadata = results.library_metadata;
validate_metadata(metadata, trials);

trial_family = metadata.orientation_family(trials.true_idx);
trial_count = metadata.num_objects(trials.true_idx);
trial_type = metadata.object_type(trials.true_idx);

families = unique(metadata.orientation_family, 'stable');
families = families(families ~= "single");
object_types = unique(metadata.object_type, 'stable');
object_counts = sort(unique(metadata.num_objects));
constellation_counts = object_counts(object_counts > 1);

summary_sigma = choose_summary_sigma(trials.sigma, summary_sigma);

fprintf('\nLoaded results from:\n  %s\n', results_path)
fprintf('Overall accuracy: %.4f\n', mean(trials.correct))
fprintf('Summary noise level: %.4g\n', summary_sigma)

fprintf('\nAccuracy by noise level:\n')
disp(grouped_accuracy(trials.sigma, trials.correct, 'sigma'))

fprintf('\nAccuracy by orientation family at the summary noise level:\n')
mask_summary = trials.sigma == summary_sigma;
disp(grouped_accuracy( ...
    trial_family(mask_summary), ...
    trials.correct(mask_summary), ...
    'orientation_family'))

fprintf('\nAccuracy by object count at the summary noise level:\n')
disp(grouped_accuracy( ...
    trial_count(mask_summary), ...
    trials.correct(mask_summary), ...
    'num_objects'))

fprintf('\nAccuracy by object class at the summary noise level:\n')
disp(grouped_accuracy( ...
    trial_type(mask_summary), ...
    trials.correct(mask_summary), ...
    'object_type'))

need_records = any([ ...
    make_paper_geometry_plot, ...
    make_paper_pole_plot, ...
    make_supplementary_all_geometry_plot, ...
    make_supplementary_all_pole_plot]);

records = struct([]);
if need_records
    records = load_constellation_records(data_folder_override);
end

if save_paper_plots && ~exist(paper_plot_folder, 'dir')
    mkdir(paper_plot_folder)
end

% ---- paper plots ----

if make_paper_geometry_plot
    plot_paper_geometry_examples( ...
        records, metadata, families, object_types, ...
        max(constellation_counts), fig_paper_geometry, dock_figures, ...
        font_name, tick_font_size, label_font_size, title_font_size, ...
        axes_line_width, line_width, constellation_shape_color, ...
        save_paper_plots, paper_plot_folder)
end

if make_paper_accuracy_noise_plot
    plot_accuracy_by_noise( ...
        trials, trial_family, trial_count, families, ...
        constellation_counts, fig_paper_accuracy_noise, dock_figures, ...
        font_name, tick_font_size, label_font_size, title_font_size, ...
        legend_font_size, axes_line_width, line_width, marker_size, ...
        save_paper_plots, paper_plot_folder)
end

if make_paper_accuracy_count_plot
    plot_accuracy_by_count( ...
        trials, trial_family, trial_count, families, object_counts, ...
        summary_sigma, fig_paper_accuracy_count, dock_figures, ...
        font_name, tick_font_size, label_font_size, legend_font_size, ...
        axes_line_width, line_width, marker_size, ...
        save_paper_plots, paper_plot_folder)
end

if make_paper_class_accuracy_plot
    plot_class_accuracy_by_count( ...
        trials, trial_family, trial_count, trial_type, ...
        families, object_types, object_counts, summary_sigma, ...
        fig_paper_class_accuracy, dock_figures, ...
        font_name, tick_font_size, label_font_size, title_font_size, ...
        legend_font_size, axes_line_width, line_width, marker_size, ...
        save_paper_plots, paper_plot_folder)
end

if make_paper_pole_plot
    plot_paper_pole_examples( ...
        records, metadata, library_metadata, families, object_types, ...
        max(constellation_counts), fig_paper_poles, dock_figures, ...
        font_name, tick_font_size, label_font_size, title_font_size, ...
        legend_font_size, axes_line_width, line_width, ...
        pole_marker_size, nominal_pole_marker_size, ...
        nominal_pole_color, constellation_pole_color, ...
        save_paper_plots, paper_plot_folder)
end

% ---- supplementary plots ----

if save_supplementary_plots && ~exist(supplementary_plot_folder, 'dir')
    mkdir(supplementary_plot_folder)
end

if make_supplementary_all_geometry_plot
    plot_all_geometry_examples( ...
        records, metadata, families, object_types, constellation_counts, ...
        fig_supplementary_geometry, dock_figures, font_name, ...
        tick_font_size, label_font_size, title_font_size, ...
        axes_line_width, line_width, constellation_shape_color, ...
        save_supplementary_plots, supplementary_plot_folder)
end

if make_supplementary_all_pole_plot
    plot_all_pole_examples( ...
        records, metadata, families, object_types, constellation_counts, ...
        fig_supplementary_poles, dock_figures, font_name, ...
        tick_font_size, label_font_size, title_font_size, ...
        axes_line_width, line_width, pole_marker_size, ...
        nominal_pole_marker_size, nominal_pole_color, ...
        constellation_pole_color, save_supplementary_plots, ...
        supplementary_plot_folder)
end

function validate_metadata(metadata, trials)
%VALIDATE_METADATA  Confirm that true_idx indexes the saved metadata.

    required = { ...
        'object_idx', ...
        'object_name', ...
        'orientation_family', ...
        'num_objects', ...
        'object_type'};

    if ~all(ismember(required, metadata.Properties.VariableNames))
        error('results.eval_metadata is missing required columns.')
    end

    if any(trials.true_idx < 1 | trials.true_idx > height(metadata))
        error('Some trials.true_idx values do not index eval_metadata.')
    end
end

function sigma = choose_summary_sigma(available_sigmas, requested_sigma)
%CHOOSE_SUMMARY_SIGMA  Use zero noise when available, otherwise the minimum.

    available_sigmas = unique(available_sigmas);

    if isempty(requested_sigma)
        if any(available_sigmas == 0)
            sigma = 0;
        else
            sigma = min(available_sigmas);
        end
    elseif any(available_sigmas == requested_sigma)
        sigma = requested_sigma;
    else
        error('summary_sigma = %g is not present in the results.', ...
            requested_sigma)
    end
end

function plot_paper_geometry_examples( ...
    records, metadata, families, object_types, num_objects, ...
    fig_id, dock_figures, font_name, tick_font_size, ...
    label_font_size, title_font_size, axes_line_width, line_width, ...
    shape_color, save_plot, output_folder)
%PLOT_PAPER_GEOMETRY_EXAMPLES  One maximum-count geometry per family/class.

    selected = select_record_grid( ...
        records, metadata, families, object_types, num_objects);
    [x_limits, y_limits] = common_geometry_limits(selected);

    fig = newfig(fig_id, 'Paper: constellation geometries', ...
        dock_figures, [100, 100, 1100, 550]);
    layout = tiledlayout(fig, numel(object_types), numel(families), ...
        'TileSpacing', 'none', 'Padding', 'compact');

    for itype = 1:numel(object_types)
        for ifamily = 1:numel(families)
            record = selected(itype, ifamily);
            ax = nexttile(layout);
            plot_constellation_geometry(ax, record, shape_color, line_width)
            format_geometry_axis( ...
                ax, x_limits, y_limits, ...
                itype == numel(object_types), ifamily == 1, ...
                font_name, tick_font_size, label_font_size, axes_line_width)

            if itype == 1
                title(ax, family_display_name(families(ifamily)), ...
                    'FontSize', title_font_size)
            end

            if ifamily == 1
                text(ax, -0.42, 0.5, object_display_name(object_types(itype)), ...
                    'Units', 'normalized', ...
                    'Rotation', 90, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontName', font_name, ...
                    'FontSize', label_font_size, ...
                    'FontWeight', 'bold', ...
                    'Clipping', 'off')
            end
        end
    end

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'constellation_geometry_examples.pdf'), 'ContentType', 'vector')
    end
end

function plot_accuracy_by_noise( ...
    trials, trial_family, trial_count, families, object_counts, ...
    fig_id, dock_figures, font_name, tick_font_size, label_font_size, ...
    title_font_size, legend_font_size, axes_line_width, line_width, ...
    marker_size, save_plot, output_folder)
%PLOT_ACCURACY_BY_NOISE  Pooled accuracy by family and object count.

    sigmas = sort(unique(trials.sigma));
    [xplot, xtick_values, xtick_labels] = sigma_plot_data(sigmas);
    colors = lines(numel(object_counts));

    fig = newfig(fig_id, 'Paper: constellation accuracy by noise', ...
        dock_figures, [100, 100, 1200, 360]);
    layout = tiledlayout(fig, 1, numel(families), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for ifamily = 1:numel(families)
        ax = nexttile(layout);
        hold(ax, 'on')

        baseline = accuracy_for_sigmas( ...
            trials, trial_family == "single", sigmas);
        plot(ax, xplot, baseline, 'k--', ...
            'LineWidth', line_width, 'DisplayName', 'N = 1')

        for icount = 1:numel(object_counts)
            count = object_counts(icount);
            mask = trial_family == families(ifamily) & ...
                trial_count == count;
            accuracy = accuracy_for_sigmas(trials, mask, sigmas);
            plot(ax, xplot, accuracy, '-o', ...
                'Color', colors(icount, :), ...
                'LineWidth', line_width, ...
                'MarkerSize', marker_size, ...
                'MarkerFaceColor', 'w', ...
                'DisplayName', sprintf('N = %d', count))
        end

        set(ax, 'XScale', 'log')
        xlim(ax, [min(xplot), max(xplot)])
        ylim(ax, [0, 1])
        xticks(ax, xtick_values)
        xticklabels(ax, cellstr(xtick_labels))
        style_axes(ax, font_name, tick_font_size, axes_line_width)
        xlabel(ax, 'noise level $\sigma$', ...
            'Interpreter', 'latex', 'FontSize', label_font_size)
        title(ax, family_display_name(families(ifamily)), ...
            'FontSize', title_font_size)

        if ifamily == 1
            ylabel(ax, 'classification accuracy', ...
                'FontSize', label_font_size)
        end

        if ifamily == numel(families)
            legend(ax, 'Location', 'best', 'FontSize', legend_font_size)
        end
        hold(ax, 'off')
    end

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'constellation_accuracy_by_noise.pdf'), 'ContentType', 'vector')
    end
end

function plot_accuracy_by_count( ...
    trials, trial_family, trial_count, families, object_counts, sigma, ...
    fig_id, dock_figures, font_name, tick_font_size, label_font_size, ...
    legend_font_size, axes_line_width, line_width, marker_size, ...
    save_plot, output_folder)
%PLOT_ACCURACY_BY_COUNT  Overall low-noise accuracy versus object count.

    fig = newfig(fig_id, 'Paper: accuracy by object count', ...
        dock_figures, [100, 100, 620, 500]);
    ax = axes(fig);
    hold(ax, 'on')
    colors = lines(numel(families));

    for ifamily = 1:numel(families)
        accuracy = nan(numel(object_counts), 1);
        for ii = 1:numel(object_counts)
            if object_counts(ii) == 1
                mask = trial_family == "single";
            else
                mask = trial_family == families(ifamily) & ...
                    trial_count == object_counts(ii);
            end
            mask = mask & trials.sigma == sigma;
            accuracy(ii) = mean(trials.correct(mask));
        end

        plot(ax, object_counts, accuracy, '-o', ...
            'Color', colors(ifamily, :), ...
            'LineWidth', line_width, ...
            'MarkerSize', marker_size, ...
            'MarkerFaceColor', 'w', ...
            'DisplayName', family_display_name(families(ifamily)))
    end

    xlim(ax, [min(object_counts), max(object_counts)])
    ylim(ax, [0, 1])
    xticks(ax, object_counts)
    style_axes(ax, font_name, tick_font_size, axes_line_width)
    xlabel(ax, 'number of objects', 'FontSize', label_font_size)
    ylabel(ax, 'classification accuracy', 'FontSize', label_font_size)
    legend(ax, 'Location', 'best', 'FontSize', legend_font_size)
    hold(ax, 'off')

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'constellation_accuracy_by_object_count.pdf'), ...
            'ContentType', 'vector')
    end
end

function plot_class_accuracy_by_count( ...
    trials, trial_family, trial_count, trial_type, families, ...
    object_types, object_counts, sigma, fig_id, dock_figures, ...
    font_name, tick_font_size, label_font_size, title_font_size, ...
    legend_font_size, axes_line_width, line_width, marker_size, ...
    save_plot, output_folder)
%PLOT_CLASS_ACCURACY_BY_COUNT  Low-noise accuracy separated by true class.

    fig = newfig(fig_id, 'Paper: class accuracy by object count', ...
        dock_figures, [100, 100, 900, 430]);
    layout = tiledlayout(fig, 1, numel(object_types), ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    colors = lines(numel(families));

    for itype = 1:numel(object_types)
        ax = nexttile(layout);
        hold(ax, 'on')

        for ifamily = 1:numel(families)
            accuracy = nan(numel(object_counts), 1);
            for ii = 1:numel(object_counts)
                if object_counts(ii) == 1
                    mask = trial_family == "single";
                else
                    mask = trial_family == families(ifamily) & ...
                        trial_count == object_counts(ii);
                end
                mask = mask & ...
                    trial_type == object_types(itype) & ...
                    trials.sigma == sigma;
                accuracy(ii) = mean(trials.correct(mask));
            end

            plot(ax, object_counts, accuracy, '-o', ...
                'Color', colors(ifamily, :), ...
                'LineWidth', line_width, ...
                'MarkerSize', marker_size, ...
                'MarkerFaceColor', 'w', ...
                'DisplayName', family_display_name(families(ifamily)))
        end

        xlim(ax, [min(object_counts), max(object_counts)])
        ylim(ax, [0, 1])
        xticks(ax, object_counts)
        style_axes(ax, font_name, tick_font_size, axes_line_width)
        xlabel(ax, 'number of objects', 'FontSize', label_font_size)
        title(ax, object_display_name(object_types(itype)), ...
            'FontSize', title_font_size)

        if itype == 1
            ylabel(ax, 'classification accuracy', ...
                'FontSize', label_font_size)
        else
            legend(ax, 'Location', 'best', 'FontSize', legend_font_size)
        end
        hold(ax, 'off')
    end

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'constellation_class_accuracy_by_object_count.pdf'), ...
            'ContentType', 'vector')
    end
end

function plot_paper_pole_examples( ...
    records, metadata, library_metadata, families, object_types, ...
    num_objects, fig_id, dock_figures, font_name, tick_font_size, ...
    label_font_size, title_font_size, legend_font_size, ...
    axes_line_width, line_width, pole_marker_size, ...
    nominal_pole_marker_size, nominal_color, constellation_color, ...
    save_plot, output_folder)
%PLOT_PAPER_POLE_EXAMPLES  Single-object and maximum-count pole overlays.

    fig = newfig(fig_id, 'Paper: constellation poles', ...
        dock_figures, [100, 100, 1100, 560]);
    layout = tiledlayout(fig, numel(object_types), numel(families), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for itype = 1:numel(object_types)
        nominal_name = library_metadata.object_name( ...
            library_metadata.object_type == object_types(itype));
        nominal = find_record(records, nominal_name);

        for ifamily = 1:numel(families)
            row = metadata.object_type == object_types(itype) & ...
                metadata.orientation_family == families(ifamily) & ...
                metadata.num_objects == num_objects;
            constellation = find_record(records, metadata.object_name(row));
            ax = nexttile(layout);
            hold(ax, 'on')
            plot(ax, real(nominal.poles), imag(nominal.poles), 'x', ...
                'Color', nominal_color, ...
                'MarkerSize', nominal_pole_marker_size, ...
                'LineWidth', line_width, 'DisplayName', 'single object')
            plot(ax, real(constellation.poles), ...
                imag(constellation.poles), 'o', ...
                'Color', constellation_color, ...
                'MarkerSize', pole_marker_size, ...
                'LineWidth', line_width, 'DisplayName', 'constellation')
            style_axes(ax, font_name, tick_font_size, axes_line_width)
            xlabel(ax, 'Re pole', 'FontSize', label_font_size)
            ylabel(ax, 'Im pole', 'FontSize', label_font_size)

            if itype == 1
                title(ax, family_display_name(families(ifamily)), ...
                    'FontSize', title_font_size)
            end
            if itype == 1 && ifamily == numel(families)
                legend(ax, 'Location', 'best', ...
                    'FontSize', legend_font_size)
            end
            hold(ax, 'off')
        end
    end

    if save_plot
        exportgraphics(fig, fullfile(output_folder, ...
            'constellation_pole_examples.pdf'), 'ContentType', 'vector')
    end
end

function plot_all_geometry_examples( ...
    records, metadata, families, object_types, object_counts, ...
    fig_start, dock_figures, font_name, tick_font_size, ...
    label_font_size, title_font_size, axes_line_width, line_width, ...
    shape_color, save_plot, output_folder)
%PLOT_ALL_GEOMETRY_EXAMPLES  Supplementary geometry grids by object class.

    for itype = 1:numel(object_types)
        fig = newfig(fig_start + itype - 1, ...
            'Supplementary: all constellation geometries', ...
            dock_figures, [100, 100, 1100, 800]);
        layout = tiledlayout(fig, numel(object_counts), numel(families), ...
            'TileSpacing', 'compact', 'Padding', 'compact');

        for icount = 1:numel(object_counts)
            for ifamily = 1:numel(families)
                row = metadata.object_type == object_types(itype) & ...
                    metadata.orientation_family == families(ifamily) & ...
                    metadata.num_objects == object_counts(icount);
                record = find_record(records, metadata.object_name(row));
                ax = nexttile(layout);
                plot_constellation_geometry(ax, record, shape_color, line_width)
                axis(ax, 'equal')
                box(ax, 'on')
                style_axes(ax, font_name, tick_font_size, axes_line_width)
                title(ax, sprintf('%s, N = %d', ...
                    family_display_name(families(ifamily)), ...
                    object_counts(icount)), 'FontSize', title_font_size)
                xlabel(ax, 'x', 'FontSize', label_font_size)
                ylabel(ax, 'y', 'FontSize', label_font_size)
            end
        end

        if save_plot
            exportgraphics(fig, fullfile(output_folder, sprintf( ...
                'all_%s_constellation_geometries.pdf', ...
                object_types(itype))), 'ContentType', 'vector')
        end
    end
end

function plot_all_pole_examples( ...
    records, metadata, families, object_types, object_counts, ...
    fig_start, dock_figures, font_name, tick_font_size, ...
    label_font_size, title_font_size, axes_line_width, line_width, ...
    pole_marker_size, nominal_pole_marker_size, nominal_color, ...
    constellation_color, save_plot, output_folder)
%PLOT_ALL_POLE_EXAMPLES  Supplementary pole grids by object class.

    for itype = 1:numel(object_types)
        nominal_row = metadata.object_type == object_types(itype) & ...
            metadata.orientation_family == "single";
        nominal = find_record(records, metadata.object_name(nominal_row));
        fig = newfig(fig_start + itype - 1, ...
            'Supplementary: all constellation poles', ...
            dock_figures, [100, 100, 1100, 800]);
        layout = tiledlayout(fig, numel(object_counts), numel(families), ...
            'TileSpacing', 'compact', 'Padding', 'compact');

        for icount = 1:numel(object_counts)
            for ifamily = 1:numel(families)
                row = metadata.object_type == object_types(itype) & ...
                    metadata.orientation_family == families(ifamily) & ...
                    metadata.num_objects == object_counts(icount);
                constellation = find_record(records, metadata.object_name(row));
                ax = nexttile(layout);
                hold(ax, 'on')
                plot(ax, real(nominal.poles), imag(nominal.poles), 'x', ...
                    'Color', nominal_color, ...
                    'MarkerSize', nominal_pole_marker_size, ...
                    'LineWidth', line_width)
                plot(ax, real(constellation.poles), ...
                    imag(constellation.poles), 'o', ...
                    'Color', constellation_color, ...
                    'MarkerSize', pole_marker_size, ...
                    'LineWidth', line_width)
                style_axes(ax, font_name, tick_font_size, axes_line_width)
                title(ax, sprintf('%s, N = %d', ...
                    family_display_name(families(ifamily)), ...
                    object_counts(icount)), 'FontSize', title_font_size)
                xlabel(ax, 'Re pole', 'FontSize', label_font_size)
                ylabel(ax, 'Im pole', 'FontSize', label_font_size)
                hold(ax, 'off')
            end
        end

        if save_plot
            exportgraphics(fig, fullfile(output_folder, sprintf( ...
                'all_%s_constellation_poles.pdf', object_types(itype))), ...
                'ContentType', 'vector')
        end
    end
end

function selected = select_record_grid( ...
    records, metadata, families, object_types, num_objects)
%SELECT_RECORD_GRID  Records for every object-type/family combination.

    selected = repmat(records(1), numel(object_types), numel(families));

    for itype = 1:numel(object_types)
        for ifamily = 1:numel(families)
            row = metadata.object_type == object_types(itype) & ...
                metadata.orientation_family == families(ifamily) & ...
                metadata.num_objects == num_objects;

            if nnz(row) ~= 1
                error(['Expected one configuration for type %s, family %s, ', ...
                    'and N = %d; found %d.'], object_types(itype), ...
                    families(ifamily), num_objects, nnz(row))
            end
            selected(itype, ifamily) = ...
                find_record(records, metadata.object_name(row));
        end
    end
end

function records = load_constellation_records(input_folder)
%LOAD_CONSTELLATION_RECORDS  Load saved geometry, poles, and configuration.

    files = dir(fullfile(input_folder, '*.mat'));
    records = struct('name', {}, 'config_info', {}, ...
        'curveX', {}, 'curveY', {}, 'poles', {});

    for ifile = 1:numel(files)
        file_path = fullfile(files(ifile).folder, files(ifile).name);
        S = load(file_path, ...
            'config_name', 'config_info', 'curveX', 'curveY', 'pols');
        required = {'config_name', 'config_info', 'curveX', 'curveY', 'pols'};

        if ~all(isfield(S, required))
            warning('Skipping incomplete data file "%s".', files(ifile).name)
            continue
        end

        record.name = string(S.config_name);
        record.config_info = S.config_info;
        record.curveX = S.curveX(:);
        record.curveY = S.curveY(:);
        record.poles = S.pols(:);
        records(end + 1) = record; %#ok<AGROW>
    end

    if isempty(records)
        error('No complete constellation data found in %s.', input_folder)
    end
end

function record = find_record(records, name)
%FIND_RECORD  Return the unique record with the requested configuration name.

    name = string(name);
    if numel(name) ~= 1
        error('Expected one configuration name; received %d.', numel(name))
    end

    idx = find(string({records.name}) == name);
    if numel(idx) ~= 1
        error('Expected one data record for "%s"; found %d.', name, numel(idx))
    end
    record = records(idx);
end

function plot_constellation_geometry(ax, record, color, line_width)
%PLOT_CONSTELLATION_GEOMETRY  Plot each disconnected cavity independently.

    centers = record.config_info.centers;
    points = [record.curveX(:), record.curveY(:)];
    valid = all(isfinite(points), 2);
    points = points(valid, :);

    distance_squared = zeros(size(points, 1), size(centers, 1));
    for ii = 1:size(centers, 1)
        distance_squared(:, ii) = ...
            (points(:, 1) - centers(ii, 1)).^2 + ...
            (points(:, 2) - centers(ii, 2)).^2;
    end
    [~, assignments] = min(distance_squared, [], 2);

    hold(ax, 'on')
    for ii = 1:size(centers, 1)
        component = points(assignments == ii, :);
        component = order_open_component(component, centers(ii, :));
        plot(ax, component(:, 1), component(:, 2), '-', ...
            'Color', color, 'LineWidth', line_width)
    end
    hold(ax, 'off')
end

function points = order_open_component(points, center)
%ORDER_OPEN_COMPONENT  Sort one cavity and break at its angular opening.

    theta = atan2(points(:, 2) - center(2), points(:, 1) - center(1));
    [theta, order] = sort(theta);
    points = points(order, :);

    gaps = [diff(theta); theta(1) + 2*pi - theta(end)];
    [~, gap_idx] = max(gaps);
    start_idx = mod(gap_idx, numel(theta)) + 1;
    order = [start_idx:numel(theta), 1:start_idx - 1];
    points = points(order, :);
end

function [x_limits, y_limits] = common_geometry_limits(records)
%COMMON_GEOMETRY_LIMITS  Common square limits for a record array.

    x_all = vertcat(records.curveX);
    y_all = vertcat(records.curveY);
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
%FORMAT_GEOMETRY_AXIS  Square boxed axes with exterior labels only.

    axis(ax, 'equal')
    pbaspect(ax, [1, 1, 1])
    xlim(ax, x_limits)
    ylim(ax, y_limits)
    box(ax, 'on')
    set(ax, 'FontName', font_name, 'FontSize', tick_font_size, ...
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

function accuracy = accuracy_for_sigmas(trials, base_mask, sigmas)
%ACCURACY_FOR_SIGMAS  Pooled accuracy for a fixed configuration grouping.

    accuracy = nan(numel(sigmas), 1);
    for isig = 1:numel(sigmas)
        mask = base_mask & trials.sigma == sigmas(isig);
        accuracy(isig) = mean(trials.correct(mask));
    end
end

function [xplot, tick_values, tick_labels] = sigma_plot_data(sigmas)
%SIGMA_PLOT_DATA  Place zero one log interval before the first positive sigma.

    sigmas = sigmas(:);
    positive = sort(sigmas(sigmas > 0));

    if any(sigmas == 0)
        if numel(positive) < 2
            error('At least two positive sigma values are required.')
        end
        zero_position = positive(1)^2/positive(2);
        xplot = [zero_position; positive];
        tick_labels = ["0"; compose('%.0e', positive)];
    else
        xplot = positive;
        tick_labels = compose('%.0e', positive);
    end
    tick_values = xplot;
end

function T = grouped_accuracy(keys, correct, key_name)
%GROUPED_ACCURACY  Accuracy grouped by a numeric or string key.

    values = unique(keys, 'stable');
    accuracy = zeros(numel(values), 1);
    ntrials = zeros(numel(values), 1);

    for ii = 1:numel(values)
        mask = keys == values(ii);
        accuracy(ii) = mean(correct(mask));
        ntrials(ii) = nnz(mask);
    end

    T = table(values(:), accuracy, ntrials, ...
        'VariableNames', {key_name, 'accuracy', 'ntrials'});
end

function label = family_display_name(family)
%FAMILY_DISPLAY_NAME  Readable orientation-family name.

    switch string(family)
        case "same_direction"
            label = 'same direction';
        case "facing_away_origin"
            label = 'facing away';
        case "random_orientation"
            label = 'random orientation';
        case "facing_origin"
            label = 'facing origin';
        case "single"
            label = 'single object';
        otherwise
            label = char(string(family));
    end
end

function label = object_display_name(object_type)
%OBJECT_DISPLAY_NAME  Readable object-class name.

    switch string(object_type)
        case "circ_o050"
            label = 'circular cavities';
        case "ell_o075"
            label = 'elliptic cavities';
        otherwise
            label = char(string(object_type));
    end
end

function fig = newfig(fig_id, name, dock_figures, position)
%NEWFIG  Create or clear one controlled figure.

    fig = figure(fig_id);
    clf(fig)
    set(fig, 'Name', name, 'NumberTitle', 'on', 'Color', 'w')

    if dock_figures
        set(fig, 'WindowStyle', 'docked')
    else
        set(fig, 'WindowStyle', 'normal', 'Position', position)
    end
end

function style_axes(ax, font_name, tick_font_size, axes_line_width)
%STYLE_AXES  Common non-geometry axes formatting.

    set(ax, 'FontName', font_name, 'FontSize', tick_font_size, ...
        'LineWidth', axes_line_width)
    grid(ax, 'on')
    box(ax, 'on')
end
