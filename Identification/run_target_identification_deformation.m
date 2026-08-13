% run_identification_nominal_poles_deformations.m
%
% Run GLRT target identification for deformed cavities using only the poles
% of the corresponding nominal cavities.

clear
clc

% ---- user-defined parameters ----

% Local folder containing downloaded .mat files.
input_folder = fullfile( ...
    'local', ...
    'data_deformation_experiment_2d');

output_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'nominal_poles_deformation_results.mat');

geometry_plot_folder = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'deformation_geometry_comparisons');

% classification params
late_time = 80.0;
sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1];

% SNR based on a representative nominal configuration.
snr_object_name = 'circ_o035_nominal';
snr_receiver_angle = 3*pi/2;

save_geometry_plots = true;

% ---------------------------------

params = build_id_params( ...
    'input_folder', input_folder, ...
    'output_path', output_path, ...
    'late_time', late_time, ...
    'noise_sigmas', sigmas, ...
    'snr_object_name', snr_object_name, ...
    'snr_receiver_angle', snr_receiver_angle, ...
    'show_progress', true);

% Load nominal and deformed configurations.
data_all = load_pole_data(params);

% Build the pole library using only nominal configurations.
nominal_library_mask = arrayfun( ...
    @is_nominal_configuration, ...
    data_all);

data_nominal_library = data_all(nominal_library_mask);

if isempty(data_nominal_library)
    error('No nominal configurations found in input folder.');
end

fprintf('\nNominal pole library:\n');
for ii = 1:numel(data_nominal_library)
    fprintf('  %s\n', data_nominal_library(ii).name);
end

fprintf('\nEvaluation set:\n');
for ii = 1:numel(data_all)
    fprintf('  %s\n', data_all(ii).name);
end

% Build the GLRT model using only nominal poles.
model_nominal_poles = build_glrt_model( ...
    data_nominal_library, ...
    params);

% Map each deformation to its corresponding nominal object class.
true_to_library_map = build_nominal_poles_truth_map( ...
    data_all, ...
    data_nominal_library);

model_nominal_poles = attach_glrt_classification_map( ...
    model_nominal_poles, ...
    data_all, ...
    data_nominal_library, ...
    true_to_library_map);

% Classify every nominal and deformed signal using nominal poles.
results = classify_glrt_data( ...
    data_all, ...
    model_nominal_poles, ...
    params);

results.method_name = "nominal_poles_deformations";
results.method_description = ...
    "Evaluation uses nominal and deformed cavities. " + ...
    "The GLRT library uses only nominal-cavity poles.";

results.params = params;
results.nominal_library_names = string( ...
    {data_nominal_library.name}).';

results.accuracy_by_sigma = accuracy_by_sigma_table( ...
    results.trials);

% Save the representative SNR corresponding to each noise sigma.
% The reader uses this mapping for accuracy-versus-SNR plots.
results.representative_snr = compute_representative_snr( ...
    data_all, ...
    model_nominal_poles, ...
    params);

% Save the classification results.
output_folder = fileparts(params.output_path);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

save(params.output_path, 'results');

fprintf('\nSaved classification results to:\n  %s\n', ...
    params.output_path);

% Compare the nominal and deformed boundaries.
plot_deformation_comparisons( ...
    input_folder, ...
    geometry_plot_folder, ...
    save_geometry_plots);

function tf = is_nominal_configuration(d)
%IS_NOMINAL_CONFIGURATION  True for a nominal deformation experiment file.

    name = string(d.name);
    tf = endsWith(name, "_nominal");
end

function map = build_nominal_poles_truth_map(eval_data, library_data)
%BUILD_NOMINAL_POLES_TRUTH_MAP
%   Map each evaluation geometry to its matching nominal library class.
%
%   For example:
%       circ_o035_m02_r001_e0p0200 -> circ_o035_nominal
%       ell_o045_m06_r001_e0p0050  -> ell_o045_nominal

    neval = numel(eval_data);
    nlib = numel(library_data);

    map = false(neval, nlib);

    library_keys = strings(nlib, 1);

    for j = 1:nlib
        library_keys(j) = canonical_object_key( ...
            library_data(j).name);
    end

    if numel(unique(library_keys)) ~= nlib
        error('The nominal pole library contains duplicate object classes.');
    end

    for i = 1:neval
        true_key = canonical_object_key(eval_data(i).name);
        map(i, :) = true_key == library_keys;

        if ~any(map(i, :))
            error( ...
                'No nominal pole-library match found for "%s".', ...
                eval_data(i).name);
        end

        if nnz(map(i, :)) > 1
            error( ...
                'Multiple nominal pole-library matches found for "%s".', ...
                eval_data(i).name);
        end
    end
end

function key = canonical_object_key(name)
%CANONICAL_OBJECT_KEY
%   Extract the underlying object tag from nominal or deformed names.
%
%   Examples:
%       circ_o035_nominal             -> circ_o035
%       circ_o035_m02_r001_e0p0200    -> circ_o035

    name = char(string(name));

    nominal_tokens = regexp( ...
        name, ...
        '^(.*)_nominal$', ...
        'tokens', ...
        'once');

    if ~isempty(nominal_tokens)
        key = string(nominal_tokens{1});
        return
    end

    deformed_tokens = regexp( ...
        name, ...
        '^(.*)_m\d+_r\d+_e[^_]+$', ...
        'tokens', ...
        'once');

    if ~isempty(deformed_tokens)
        key = string(deformed_tokens{1});
        return
    end

    error('Could not extract object tag from name "%s".', name);
end

function acc = accuracy_by_sigma_table(trials)
%ACCURACY_BY_SIGMA_TABLE  Compute accuracy grouped by noise sigma.

    sigmas = unique(trials.sigma, 'stable');
    accuracy = zeros(numel(sigmas), 1);
    ntrials = zeros(numel(sigmas), 1);

    for isig = 1:numel(sigmas)
        mask = trials.sigma == sigmas(isig);

        accuracy(isig) = mean(trials.correct(mask));
        ntrials(isig) = nnz(mask);
    end

    acc = table( ...
        sigmas, ...
        accuracy, ...
        ntrials, ...
        'VariableNames', { ...
            'sigma', ...
            'accuracy', ...
            'ntrials'});
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
