% read_constellation_classification_results_script_v6.m
%
% Read and summarize constellation classification results.
%
% This is a script. Edit the parameters below and run it directly.

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

% Empty means use all sigmas in the confusion table.
sigma_to_report = [];

make_summary_plots = true;
print_latex_config_type_table = true;
make_config_plots = false;

make_shape_panel = true;
make_pole_panel = true;
make_accuracy_text_panel = true;

% Empty means plot all evaluation configurations.
% These are true_idx values.
configs_to_plot = [];

% Empty means use all sigmas for per-configuration accuracy.
sigmas_for_config_accuracy = [];

% Dock every figure into the MATLAB desktop.
dock_figures = true;

% 'separate' gives one docked tab per configuration: 10, 11, 12, ...
% 'reuse' overwrites figure 10 at each configuration.
config_figure_mode = 'separate';

save_config_plots = false;
config_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'single_poles_config_diagnostics');

% ---- plotting parameters ----

font_name = 'Times New Roman';
tick_font_size = 10;
label_font_size = 12;
title_font_size = 12;
subtitle_font_size = 10;
legend_font_size = 9;
text_font_size = 10;

line_width = 1.25;
axes_line_width = 0.75;
marker_size = 6;
shape_marker_size = 5;
pole_marker_size = 6;
single_pole_marker_size = 8;

config_pole_marker = 'o';
single_pole_marker = 'x';
config_pole_color = [0.0000, 0.4470, 0.7410];
single_pole_color = [0.8500, 0.3250, 0.0980];

sigma_char = char(963);              % Unicode sigma; no TeX interpreter needed
sigma_xlabel = sigma_char;
accuracy_ylabel = 'average classification accuracy';
num_objects_xlabel = 'number of objects';
config_type_xlabel = 'configuration type';
shape_xlabel = 'x';
shape_ylabel = 'y';
pole_xlabel = 'Re pole';
pole_ylabel = 'Im pole';

summary_figure_position = [100, 100, 520, 400];
config_figure_position = [100, 100, 1150, 360];

fig_accuracy_sigma = 1;
fig_accuracy_nobj = 2;
fig_accuracy_config_type = 3;
fig_config_start = 10;

interpreter = 'none';
grid_state = 'on';

config_type_table_caption = ...
    'Average classification accuracy by configuration type.';
config_type_table_label = 'tab:config-type-accuracy';

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

fprintf('\nLoaded results from:\n  %s\n', results_path);

if isfield(results, 'method_name')
    fprintf('\nMethod: %s\n', results.method_name);
end

if isfield(results, 'method_description')
    fprintf('%s\n', results.method_description);
end

fprintf('\nOverall accuracy: %.4f\n', mean(trials.correct));

acc_by_sigma = grouped_accuracy(trials.sigma, trials.correct, 'sigma');

acc_by_config = grouped_accuracy(trials.true_idx, trials.correct, 'true_idx');
acc_by_config.true_name = eval_names(acc_by_config.true_idx);
acc_by_config = movevars(acc_by_config, 'true_name', 'After', 'true_idx');

true_types = strings(height(trials), 1);
num_objects = zeros(height(trials), 1);
config_types = strings(height(trials), 1);

for ii = 1:height(trials)
    true_types(ii) = object_type_from_name(true_names(ii));
    num_objects(ii) = num_objects_from_name(true_names(ii));
    config_types(ii) = config_type_from_name(true_names(ii));
end

acc_by_type = grouped_accuracy(true_types, trials.correct, 'object_type');
acc_by_nobj = grouped_accuracy(num_objects, trials.correct, 'num_objects');
acc_by_config_type = grouped_accuracy(config_types, trials.correct, 'config_type');
pred_counts = prediction_counts(trials, library_names);

fprintf('\nAccuracy by sigma:\n')
disp(acc_by_sigma)

fprintf('\nAccuracy by true configuration:\n')
disp(acc_by_config)

fprintf('\nAccuracy by object type:\n')
disp(acc_by_type)

fprintf('\nAccuracy by number of objects:\n')
disp(acc_by_nobj)

fprintf('\nAccuracy by configuration type:\n')
disp(acc_by_config_type)

fprintf(['\nThe configuration-type plot/table shows average classification accuracy, ', ...
         'computed as mean(trials.correct) within each true configuration type.\n'])

if print_latex_config_type_table
    fprintf('\nLaTeX table: average classification accuracy by configuration type\n');
    fprintf('\\begin{tabular}{lc}\n');
    fprintf('\\hline\n');
    fprintf('Configuration type & Accuracy (\\%%) \\\\\n');
    fprintf('\\hline\n');

    for i = 1:height(acc_by_config_type)
        config_label = char(acc_by_config_type.config_type(i));

        % display formatting only
        config_label = strrep(config_label, '_', ' ');
        config_label = lower(config_label);

        if ~isempty(config_label)
            config_label(1) = upper(config_label(1));
        end

        fprintf('%s & %.2f \\\\\n', ...
            config_label, ...
            100 * acc_by_config_type.accuracy(i));
    end

    fprintf('\\hline\n');
    fprintf('\\end{tabular}\n');
end

fprintf('\nPrediction counts:\n')
disp(pred_counts)

if isempty(sigma_to_report)
    mask_report = true(height(trials), 1);
else
    mask_report = trials.sigma == sigma_to_report;

    if ~any(mask_report)
        error('Requested sigma_to_report = %g was not found.', sigma_to_report);
    end
end

fprintf('\nConfusion table')
if ~isempty(sigma_to_report)
    fprintf(' for sigma = %g', sigma_to_report)
end
fprintf(':\n')

disp(confusion_counts(trials(mask_report, :), eval_names, library_names))

if make_summary_plots
    sigmas = acc_by_sigma.sigma(:);
    accuracy = acc_by_sigma.accuracy(:);

    zero_mask = sigmas == 0;
    pos_mask = sigmas > 0;
    sigmas_pos = sigmas(pos_mask);
    accuracy_pos = accuracy(pos_mask);
    [sigmas_pos, idx] = sort(sigmas_pos);
    accuracy_pos = accuracy_pos(idx);

    if any(zero_mask)
        if numel(sigmas_pos) < 2
            error('Need at least two positive sigma values to place sigma = 0.');
        end

        sigma_zero_plot = sigmas_pos(1)^2 / sigmas_pos(2);
        xplot = [sigma_zero_plot; sigmas_pos];
        yplot = [accuracy(find(zero_mask, 1, 'first')); accuracy_pos];
        xticks_vals = xplot;
        xticks_labels = ["0"; string(sigmas_pos)];
    else
        xplot = sigmas_pos;
        yplot = accuracy_pos;
        xticks_vals = xplot;
        xticks_labels = string(sigmas_pos);
    end

    newfig(fig_accuracy_sigma, 'Accuracy by sigma', dock_figures, summary_figure_position);
    
    semilogx(xplot, yplot, '-o', ...
        'LineWidth', 2, ...
        'MarkerSize', 7)
    
    grid on
    xlabel("Noise Level")
    ylabel("Overall Accuracy")
    title("Classification accuracy by noise level")
    
    xlim([min(xplot), max(xplot)])
    ylim([0, 1])
    pbaspect([1 1 1])
    
    xticks(xticks_vals)
    xticklabels(cellstr(xticks_labels))
    
    ax = gca;
    set(ax, ...
        'FontSize', 15, ...
        'FontName', font_name, ...
        'LineWidth', axes_line_width)
    
    drawnow

    newfig(fig_accuracy_nobj, 'Accuracy by number of objects', dock_figures, summary_figure_position);
    plot(acc_by_nobj.num_objects, acc_by_nobj.accuracy, '-o', ...
        'LineWidth', line_width, ...
        'MarkerSize', marker_size)
    style_axes(gca, font_name, tick_font_size, axes_line_width, grid_state)
    xlabel(num_objects_xlabel, 'FontSize', label_font_size, 'Interpreter', interpreter)
    ylabel(accuracy_ylabel, 'FontSize', label_font_size, 'Interpreter', interpreter)
    title('Classification accuracy by constellation size', 'FontSize', title_font_size, 'Interpreter', 'none')
    drawnow

    newfig(fig_accuracy_config_type, 'Accuracy by configuration type', dock_figures, summary_figure_position);
    cats = strings(height(acc_by_config_type), 1);

    for i = 1:height(acc_by_config_type)
        cats(i) = acc_by_config_type.config_type(i);
        cats(i) = replace(cats(i), '_', ' ');
        cats(i) = lower(cats(i));

        if strlength(cats(i)) > 0
            tmp = char(cats(i));
            tmp(1) = upper(tmp(1));
            cats(i) = string(tmp);
        end
    end

    cats = cellstr(cats);
    bar(categorical(cats, cats, 'Ordinal', true), acc_by_config_type.accuracy)
    ylim([0, 1])
    style_axes(gca, font_name, tick_font_size, axes_line_width, grid_state)
    xlabel(config_type_xlabel, 'FontSize', label_font_size, 'Interpreter', interpreter)
    ylabel(accuracy_ylabel, 'FontSize', label_font_size, 'Interpreter', interpreter)
    title('Average classification accuracy by configuration type', 'FontSize', title_font_size, 'Interpreter', 'none')
    drawnow
end

if make_config_plots && any([make_shape_panel, make_pole_panel, make_accuracy_text_panel])
    true_ids = unique(trials.true_idx, 'stable');

    if ~isempty(configs_to_plot)
        true_ids = true_ids(ismember(true_ids, configs_to_plot));
    end

    if save_config_plots && ~exist(config_plot_folder, 'dir')
        mkdir(config_plot_folder)
    end

    npanels = nnz([make_shape_panel, make_pole_panel, make_accuracy_text_panel]);
    fprintf('\nPlotting %d configuration diagnostic figure(s).\n', numel(true_ids))

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
        [single_D, single_name] = load_single_reference( ...
            config_name, ...
            library_objects, ...
            library_names, ...
            data_folder_override);

        if strcmp(config_figure_mode, 'separate')
            fig_id = fig_config_start + ii - 1;
        elseif strcmp(config_figure_mode, 'reuse')
            fig_id = fig_config_start;
        else
            error('config_figure_mode must be ''separate'' or ''reuse''.')
        end

        fh = newfig(fig_id, char(config_name), dock_figures, config_figure_position);
        tiledlayout(1, npanels, 'TileSpacing', 'compact', 'Padding', 'compact')

        if make_shape_panel
            nexttile

            if isempty(D.curveX) || isempty(D.curveY)
                text(0.5, 0.5, 'curveX/curveY missing', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', text_font_size)
                axis off
            else
                plot(D.curveX, D.curveY, '.', 'MarkerSize', shape_marker_size)
                axis equal
                style_axes(gca, font_name, tick_font_size, axes_line_width, grid_state)
                xlabel(shape_xlabel, 'FontSize', label_font_size, 'Interpreter', interpreter)
                ylabel(shape_ylabel, 'FontSize', label_font_size, 'Interpreter', interpreter)
            end

            title('shape', 'Interpreter', 'none', 'FontSize', title_font_size)
        end

        if make_pole_panel
            nexttile

            if isempty(D.poles) && isempty(single_D.poles)
                text(0.5, 0.5, 'poles missing', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', text_font_size)
                axis off
            else
                hold on

                if ~isempty(D.poles)
                    plot(real(D.poles), imag(D.poles), config_pole_marker, ...
                        'Color', config_pole_color, ...
                        'MarkerSize', pole_marker_size, ...
                        'LineWidth', line_width, ...
                        'DisplayName', 'configuration poles')
                end

                if ~isempty(single_D.poles)
                    plot(real(single_D.poles), imag(single_D.poles), single_pole_marker, ...
                        'Color', single_pole_color, ...
                        'MarkerSize', single_pole_marker_size, ...
                        'LineWidth', line_width, ...
                        'DisplayName', 'single-object poles')
                end

                style_axes(gca, font_name, tick_font_size, axes_line_width, grid_state)
                xlabel(pole_xlabel, 'FontSize', label_font_size, 'Interpreter', interpreter)
                ylabel(pole_ylabel, 'FontSize', label_font_size, 'Interpreter', interpreter)

                xl = xlim;
                yl = ylim;
                plot(xl, [0 0], 'k:', 'HandleVisibility', 'off')
                plot([0 0], yl, 'k:', 'HandleVisibility', 'off')
                xlim(xl)
                ylim(yl)

                legend('Location', 'best', 'Interpreter', 'none', 'FontSize', legend_font_size)

                if strlength(single_name) > 0
                    subtitle("reference: " + single_name, 'Interpreter', 'none', 'FontSize', subtitle_font_size)
                else
                    subtitle('no single-object reference found', 'Interpreter', 'none', 'FontSize', subtitle_font_size)
                end

                hold off
            end

            title('poles', 'Interpreter', 'none', 'FontSize', title_font_size)
        end

        if make_accuracy_text_panel
            nexttile

            mask_config = trials.true_idx == true_id;

            if isempty(sigmas_for_config_accuracy)
                mask_sigma = true(height(trials), 1);
                sigma_text = "all sigmas";
            else
                mask_sigma = false(height(trials), 1);

                for isig = 1:numel(sigmas_for_config_accuracy)
                    mask_sigma = mask_sigma | trials.sigma == sigmas_for_config_accuracy(isig);
                end

                sigma_text = "selected sigmas";
            end

            mask = mask_config & mask_sigma;

            if ~any(mask)
                text(0.5, 0.5, 'no matching trials', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', text_font_size)
            else
                ntrials = nnz(mask);
                ncorrect = nnz(trials.correct(mask));
                receiver_ids = unique(trials.receiver_idx(mask), 'stable');
                sigmas_used = unique(trials.sigma(mask), 'stable');

                lines = strings(0, 1);
                lines(end + 1) = sprintf('accuracy: %.2f%%', 100 * ncorrect / ntrials);
                lines(end + 1) = sprintf('correct: %d / %d', ncorrect, ntrials);
                lines(end + 1) = sprintf('receiver signals: %d', numel(receiver_ids));
                lines(end + 1) = sprintf('noise levels: %d', numel(sigmas_used));
                lines(end + 1) = "";
                lines(end + 1) = "sigmas: " + sigma_text;
                lines(end + 1) = "";

                for isig = 1:numel(sigmas_used)
                    sigma = sigmas_used(isig);
                    smask = mask & trials.sigma == sigma;

                    if sigma == 0
                        sigma_label = sprintf('%s = 0', sigma_char);
                    else
                        sigma_label = sprintf('%s = %.1e', sigma_char, sigma);
                    end

                    lines(end + 1) = sprintf('%s: %.2f%%  (%d/%d)', ...
                        sigma_label, ...
                        100 * mean(trials.correct(smask)), ...
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
            title('classification accuracy', 'Interpreter', 'none', 'FontSize', title_font_size)
        end

        sgtitle(config_name, 'Interpreter', 'none', 'FontSize', title_font_size)
        drawnow

        if save_config_plots
            fname = matlab.lang.makeValidName(char(config_name));
            saveas(fh, fullfile(config_plot_folder, [fname, '.png']))
        end
    end
end

function fh = newfig(fig_id, fig_name, dock_figures, fig_position)
%NEWFIG  Create or clear one controlled figure.

    fh = figure(fig_id);
    clf(fh)
    set(fh, 'Name', fig_name, 'NumberTitle', 'on')

    if dock_figures
        set(fh, 'WindowStyle', 'docked')
    else
        set(fh, 'WindowStyle', 'normal', 'Position', fig_position)
    end
end

function style_axes(ax, font_name, tick_font_size, axes_line_width, grid_state)
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

    T = table(vals(:), accuracy, ntrials, ...
        'VariableNames', {key_name, 'accuracy', 'ntrials'});
end

function counts = prediction_counts(trials, library_names)
%PREDICTION_COUNTS  Count predictions by library object.

    pred_ids = unique(trials.pred_idx, 'stable');
    pred_name = library_names(pred_ids);
    count = zeros(numel(pred_ids), 1);

    for ii = 1:numel(pred_ids)
        count(ii) = nnz(trials.pred_idx == pred_ids(ii));
    end

    fraction = count / height(trials);
    counts = table(pred_ids, pred_name, count, fraction, ...
        'VariableNames', {'pred_idx', 'pred_name', 'count', 'fraction'});
end

function conf = confusion_counts(trials, eval_names, library_names)
%CONFUSION_COUNTS  Count true configuration versus predicted library class.

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

    pred_vars = matlab.lang.makeValidName(cellstr(library_names(pred_ids)));
    pred_vars = matlab.lang.makeUniqueStrings(pred_vars);
    conf = array2table(counts, ...
        'VariableNames', pred_vars, ...
        'RowNames', cellstr(eval_names(true_ids)));
end

function D = load_diagnostics(file_path)
%LOAD_DIAGNOSTICS  Load plotting data from one original .mat file.

    if ~isfile(file_path)
        error('Could not find original data file:\n  %s', file_path);
    end

    S = load(file_path, 'curveX', 'curveY', 'pols', 'receiver_angles');
    D.curveX = field_or_empty(S, 'curveX');
    D.curveY = field_or_empty(S, 'curveY');
    D.poles = field_or_empty(S, 'pols');
    D.receiver_angles = field_or_empty(S, 'receiver_angles');
    D.poles = D.poles(:);
    D.receiver_angles = D.receiver_angles(:);
end

function [D, single_name] = load_single_reference( ...
    config_name, ...
    library_objects, ...
    library_names, ...
    data_folder_override)
%LOAD_SINGLE_REFERENCE  Load matching single-object pole data.

    target_type = object_type_from_name(config_name);

    if target_type == "unknown" || target_type == "mixed"
        D = empty_diagnostics();
        single_name = "";
        return
    end

    for ii = 1:height(library_objects)
        lib_name = library_names(ii);

        if ~is_single_object_name(lib_name)
            continue
        end

        if object_type_from_name(lib_name) ~= target_type
            continue
        end

        file_path = resolve_data_file( ...
            library_objects, ...
            ii, ...
            data_folder_override, ...
            {'file_path'}, ...
            {'file_name'});

        D = load_diagnostics(file_path);
        single_name = lib_name;
        return
    end

    D = empty_diagnostics();
    single_name = "";
end

function file_path = resolve_data_file(T, idx, folder_override, path_cols, name_cols)
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

    error(['Could not resolve original data file.\n', ...
           'Saved path:\n  %s\n', ...
           'Fallback path:\n  %s'], ...
           file_path, candidate);
end

function x = field_or_empty(S, field_name)
%FIELD_OR_EMPTY  Return struct field or empty array.

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
%TABLE_STRING  Return table value from first available column.

    vars = T.Properties.VariableNames;

    for jj = 1:numel(candidates)
        if any(strcmp(vars, candidates{jj}))
            val = string(T.(candidates{jj})(idx));
            return
        end
    end

    error('Could not find any expected table column.');
end

function tf = is_single_object_name(name)
%IS_SINGLE_OBJECT_NAME  True for single-object library names.

    name = string(name);
    tf = startsWith(name, "single_") || num_objects_from_name(name) == 1;
end

function type = object_type_from_name(name)
%OBJECT_TYPE_FROM_NAME  Extract object type from generated names.

    tags = object_tags_from_name(name);
    tags = unique(tags, 'stable');

    if isempty(tags)
        type = "unknown";
    elseif numel(tags) == 1
        type = tags(1);
    else
        type = "mixed";
    end
end

function nobj = num_objects_from_name(name)
%NUM_OBJECTS_FROM_NAME  Extract N from generated constellation names.

    name = string(name);
    tokens = regexp(char(name), '_N(\d+)_', 'tokens', 'once');

    if ~isempty(tokens)
        nobj = str2double(tokens{1});
    elseif startsWith(name, "single_")
        nobj = 1;
    else
        tags = object_tags_from_name(name);
        nobj = numel(tags);

        if nobj == 0
            nobj = NaN;
        end
    end
end

function config_type = config_type_from_name(name)
%CONFIG_TYPE_FROM_NAME  Extract configuration type from name.

    name = string(name);

    if startsWith(name, "single_")
        config_type = "single";
    elseif startsWith(name, "same_direction")
        config_type = "same_direction";
    elseif startsWith(name, "facing_away_origin")
        config_type = "facing_away_origin";
    elseif startsWith(name, "random_orientation")
        config_type = "random_orientation";
    elseif startsWith(name, "facing_origin")
        config_type = "facing_origin";
    elseif startsWith(name, "same_dir")
        config_type = "same_dir";
    elseif startsWith(name, "facing")
        config_type = "facing";
    else
        config_type = "unknown";
    end
end


function tags = object_tags_from_name(name)
%OBJECT_TAGS_FROM_NAME  Extract known object tags from a name.

    hits = regexp(char(string(name)), '(circ_o050|ell_o075)', 'match');
    tags = string(hits(:));
end
