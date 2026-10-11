function partial_results = run_plane_identification_partial_circumference()
%RUN_PLANE_IDENTIFICATION_PARTIAL_CIRCUMFERENCE  Identify planes on angle subsets.
%
% For a retained circumference percentage X, observations are restricted
% to the centered angular interval
%
%                     [-pi*X/100, pi*X/100].
%
% Thus X = 25 keeps [-pi/4, pi/4], and X = 50 keeps
% [-pi/2, pi/2]. Excluded classes are removed from both the evaluation data
% and the pole library; only observations in the interval are classified.

clc

cfg.input_folder = 'local/plane_data_360';
cfg.output_path = fullfile('Identification', 'classification_results', ...
    'plane_partial_circumference_results.mat');
cfg.percentages = [25, 50, 75, 100];
cfg.late_time = 40.0;
cfg.noise_sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
cfg.excluded_object_names = "twinjet1";
cfg.random_seed = 2026;
cfg.show_progress = true;

validate_percentages(cfg.percentages);

params = build_id_params( ...
    'input_folder', cfg.input_folder, ...
    'output_path', cfg.output_path, ...
    'late_time', cfg.late_time, ...
    'noise_sigmas', cfg.noise_sigmas, ...
    'random_seed', cfg.random_seed, ...
    'snr_report_type', "averaged", ...
    'show_progress', cfg.show_progress);

library_data = load_pole_data(params);
[library_data, excluded_data] = exclude_objects( ...
    library_data, cfg.excluded_object_names);
model = build_glrt_model(library_data, params);
percentages = unique(cfg.percentages(:).', 'stable');

partial_results = struct();
partial_results.configuration = cfg;
partial_results.percentages = percentages;
partial_results.excluded_objects = string({excluded_data.name}).';
partial_results.runs = struct([]);

for ipercentage = 1:numel(percentages)
    percentage = percentages(ipercentage);
    half_width = pi * percentage / 100;
    [evaluation_data, receiver_indices] = select_observations( ...
        library_data, half_width);

    fprintf('\nRetaining %g%% of the circumference: [%g, %g] rad\n', ...
        percentage, -half_width, half_width);

    results = classify_glrt_data(evaluation_data, model, params);
    results.observation_percentage = percentage;
    results.observation_angle_interval = [-half_width, half_width];
    results.original_receiver_indices = receiver_indices;

    partial_results.runs(ipercentage).percentage = percentage;
    partial_results.runs(ipercentage).angle_interval = ...
        [-half_width, half_width];
    partial_results.runs(ipercentage).receiver_count_by_object = ...
        cellfun(@numel, receiver_indices).';
    partial_results.runs(ipercentage).results = results;
end

output_folder = fileparts(cfg.output_path);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

save(cfg.output_path, 'partial_results', '-v7.3');
fprintf('\nSaved partial-circumference results to:\n  %s\n', ...
    cfg.output_path);
end


function [included, excluded] = exclude_objects(data, excluded_names)
%EXCLUDE_OBJECTS  Remove named classes from evaluation and the GLRT library.

% Canonical matching lets names such as "twinjet1" match the dataset name
% "twinjet_plane1" and its corresponding file stem.

excluded_names = string(excluded_names(:));
excluded_keys = strings(size(excluded_names));

for iname = 1:numel(excluded_names)
    excluded_keys(iname) = canonical_object_name(excluded_names(iname));
end

exclude = false(size(data));
available_keys = strings(numel(data), 2);

for iobject = 1:numel(data)
    [~, file_stem] = fileparts(data(iobject).file_name);
    available_keys(iobject, :) = [ ...
        canonical_object_name(data(iobject).name), ...
        canonical_object_name(file_stem)];
    exclude(iobject) = any(ismember(available_keys(iobject, :), excluded_keys));
end

for iname = 1:numel(excluded_keys)
    if ~any(available_keys == excluded_keys(iname), 'all')
        error('run_plane_identification_partial_circumference:MissingExclusion', ...
            'Could not find excluded object "%s" in the input data.', ...
            excluded_names(iname));
    end
end

excluded = data(exclude);
included = data(~exclude);

if isempty(included)
    error('run_plane_identification_partial_circumference:EmptyLibrary', ...
        'Object exclusions removed every classification class.');
end

fprintf('Excluded classification class%s: %s\n', ...
    plural_suffix(numel(excluded)), strjoin(string({excluded.name}), ', '));
end


function key = canonical_object_name(name)
%CANONICAL_OBJECT_NAME  Normalize plane object aliases for matching.

key = lower(string(name));
key = erase(key, "plane");
key = regexprep(key, '[^a-z0-9]', '');
end


function suffix = plural_suffix(count)
%PLURAL_SUFFIX  Return an s when a printed noun is plural.

if count == 1
    suffix = '';
else
    suffix = 'es';
end
end


function validate_percentages(percentages)
%VALIDATE_PERCENTAGES  Require unique, finite percentages in (0, 100].

if isempty(percentages) || ~isnumeric(percentages) ...
        || any(~isfinite(percentages)) ...
        || any(percentages <= 0 | percentages > 100)
    error('run_plane_identification_partial_circumference:BadPercentage', ...
        'cfg.percentages must contain finite values in (0, 100].');
end
end


function [selected, receiver_indices] = select_observations(data, half_width)
%SELECT_OBSERVATIONS  Keep receiver angles in the centered periodic window.

selected = data;
receiver_indices = cell(numel(data), 1);
angle_tolerance = 100 * eps(pi);

for iobject = 1:numel(data)
    if isempty(data(iobject).receiver_angles)
        error('select_observations:MissingReceiverAngles', ...
            '%s has no receiver_angles field.', data(iobject).file_path);
    end

    wrapped_angles = wrap_angle(data(iobject).receiver_angles);
    keep = abs(wrapped_angles) <= half_width + angle_tolerance;

    if ~any(keep)
        error('select_observations:NoObservations', ...
            'The selected interval contains no observations for %s.', ...
            data(iobject).name);
    end

    receiver_indices{iobject} = find(keep);
    selected(iobject).signals = data(iobject).signals(keep, :);
    selected(iobject).receiver_angles = ...
        data(iobject).receiver_angles(keep);
end
end


function angle = wrap_angle(angle)
%WRAP_ANGLE  Map angles to [-pi, pi).

angle = mod(angle + pi, 2 * pi) - pi;
end
