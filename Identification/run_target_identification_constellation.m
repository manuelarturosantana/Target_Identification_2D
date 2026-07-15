% run_identification_single_poles_constellations.m
%
% Run GLRT target identification for same-object constellations using only
% the poles of the single-object configurations.

clear
clc

% ---- user-defined parameters ----

% Local folder containing downloaded .mat files.
input_folder = fullfile( ...
    'local', ...
    'data_same_object_constellations_2d_fast');

output_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'single_poles_constellation_results.mat');

% classification params
late_time = 80.0;
sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-1, 1e-1, 5e-1];

% SNR based on representative single-object library element.
% For the current generated names, use e.g.
%   single_N1_circ_o050
%   single_N1_ell_o075
snr_object_name = 'single_N1_circ_o050';
snr_receiver_angle = 3*pi/2;

% ---------------------------------

params = build_id_params( ...
    'input_folder', input_folder, ...
    'output_path', output_path, ...
    'late_time', late_time, ...
    'noise_sigmas', sigmas, ...
    'snr_object_name', snr_object_name, ...
    'snr_receiver_angle', snr_receiver_angle, ...
    'show_progress', true);

% Load all configurations: singles and constellations.
data_all = load_pole_data(params);

% Build pole library using only single-object configurations.
single_library_mask = arrayfun(@is_single_object_config, data_all);
data_single_library = data_all(single_library_mask);

if isempty(data_single_library)
    error('No single-object configurations found in input folder.');
end

fprintf('\nSingle-object pole library:\n');
for ii = 1:numel(data_single_library)
    fprintf('  %s\n', data_single_library(ii).name);
end

fprintf('\nEvaluation set:\n');
for ii = 1:numel(data_all)
    fprintf('  %s\n', data_all(ii).name);
end

% Build GLRT model from single-object poles only.
model_single_poles = build_glrt_model(data_single_library, params);

% Map each constellation to the correct single-object class.
true_to_library_map = build_single_poles_truth_map( ...
    data_all, ...
    data_single_library);

model_single_poles = attach_glrt_classification_map( ...
    model_single_poles, ...
    data_all, ...
    data_single_library, ...
    true_to_library_map);

% Classify every signal using only the single-object pole library.
results = classify_glrt_data( ...
    data_all, ...
    model_single_poles, ...
    params);

results.method_name = "single_poles_constellations";
results.method_description = ...
    "Evaluation uses all same-object constellations. Library uses only single-object poles.";

results.params = params;
results.accuracy_by_sigma = accuracy_by_sigma_table(results.trials);

% Save.
output_folder = fileparts(params.output_path);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

save(params.output_path, 'results');

fprintf('\nSaved classification results to:\n  %s\n', params.output_path);

function tf = is_single_object_config(d)
%IS_SINGLE_OBJECT_CONFIG  True for generated single-object configurations.

    name = string(d.name);

    if startsWith(name, "single_")
        tf = true;
        return
    end

    nobj = get_num_objects_from_name(name);
    tf = nobj == 1;
end

function map = build_single_poles_truth_map(eval_data, library_data)
%BUILD_SINGLE_POLES_TRUTH_MAP
%   Map each evaluation object to the matching single-object library class.
%
%   For same-object constellations:
%       same_direction_N4_circ_o050 -> single_N1_circ_o050
%       facing_origin_N3_ell_o075   -> single_N1_ell_o075

    neval = numel(eval_data);
    nlib = numel(library_data);

    map = false(neval, nlib);

    for i = 1:neval
        true_key = canonical_single_object_key(eval_data(i).name);

        for j = 1:nlib
            pred_key = canonical_single_object_key(library_data(j).name);
            map(i, j) = true_key == pred_key;
        end

        if ~any(map(i, :))
            error('No valid single-poles library match found for "%s".', ...
                eval_data(i).name);
        end
    end
end

function key = canonical_single_object_key(name)
%CANONICAL_SINGLE_OBJECT_KEY
%   Canonical object type used for single-pole classification.

    tags = get_object_tags_from_name(name);
    tags = unique(tags, 'stable');

    if numel(tags) ~= 1
        error(['Single-poles classification expects same-object ', ...
               'configurations only. Got "%s".'], string(name));
    end

    key = tags(1);
end

function tags = get_object_tags_from_name(name)
%GET_OBJECT_TAGS_FROM_NAME  Extract object tags from generated names.

    name = char(string(name));

    expr = '(circ_o050|ell_o075)';
    hits = regexp(name, expr, 'match');

    tags = string(hits(:));

    if isempty(tags)
        error('Could not extract object tag from name "%s".', name);
    end
end

function nobj = get_num_objects_from_name(name)
%GET_NUM_OBJECTS_FROM_NAME  Extract N from names like same_direction_N3_circ_o050.

    name = char(string(name));

    tokens = regexp(name, '_N(\d+)_', 'tokens', 'once');

    if isempty(tokens)
        nobj = NaN;
    else
        nobj = str2double(tokens{1});
    end
end

function acc = accuracy_by_sigma_table(trials)
%ACCURACY_BY_SIGMA_TABLE  Compute accuracy grouped by sigma.

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