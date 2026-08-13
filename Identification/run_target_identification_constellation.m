% run_identification_single_poles_constellations.m
%
% Run GLRT target identification for same-object constellations using only
% the poles of the two single-object configurations.

clear
clc

% ---- user-defined parameters ----

input_folder = fullfile( ...
    'local', ...
    'data_same_object_constellations_2d_fast');

output_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'single_poles_constellation_results.mat');

late_time = 80.0;
noise_sigmas = [ ...
    0, ...
    1e-6, ...
    1e-5, ...
    1e-4, ...
    1e-3, ...
    1e-2, ...
    1e-1, ...
    5e-1];

% This is a reference SNR, not the SNR of every constellation signal.
snr_object_name = 'single_N1_circ_o050';
snr_receiver_angle = 3*pi/2;

% ---------------------------------

params = build_id_params( ...
    'input_folder', input_folder, ...
    'output_path', output_path, ...
    'late_time', late_time, ...
    'noise_sigmas', noise_sigmas, ...
    'snr_object_name', snr_object_name, ...
    'snr_receiver_angle', snr_receiver_angle, ...
    'show_progress', true);

% Load every single-object and constellation configuration.
data_all = load_pole_data(params);
eval_metadata = build_constellation_metadata(data_all);

% The pole library contains only the two single-object configurations.
single_library_mask = ...
    eval_metadata.orientation_family == "single" & ...
    eval_metadata.num_objects == 1;

data_single_library = data_all(single_library_mask);
library_metadata = eval_metadata(single_library_mask, :);
library_metadata.object_idx = (1:height(library_metadata)).';

if isempty(data_single_library)
    error('No single-object configurations found in input folder.')
end

if numel(unique(library_metadata.object_type)) ~= ...
        numel(data_single_library)
    error('The single-object library contains duplicate object classes.')
end

fprintf('\nSingle-object pole library:\n')
for ii = 1:numel(data_single_library)
    fprintf('  %s\n', data_single_library(ii).name)
end

fprintf('\nEvaluation set:\n')
for ii = 1:numel(data_all)
    fprintf('  %s\n', data_all(ii).name)
end

model_single_poles = build_glrt_model(data_single_library, params);

true_to_library_map = build_single_poles_truth_map( ...
    eval_metadata, ...
    library_metadata);

model_single_poles = attach_glrt_classification_map( ...
    model_single_poles, ...
    data_all, ...
    data_single_library, ...
    true_to_library_map);

results = classify_glrt_data( ...
    data_all, ...
    model_single_poles, ...
    params);

results.method_name = "single_poles_constellations";
results.method_description = [ ...
    "Evaluation uses every same-object constellation. " + ...
    "The GLRT library uses only single-object poles."];

results.params = params;
results.eval_metadata = eval_metadata;
results.library_metadata = library_metadata;
results.true_to_library_map = true_to_library_map;
results.accuracy_by_sigma = accuracy_by_sigma_table(results.trials);

output_folder = fileparts(params.output_path);
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

save(params.output_path, 'results')

fprintf('\nSaved classification results to:\n  %s\n', params.output_path)

function metadata = build_constellation_metadata(data)
%BUILD_CONSTELLATION_METADATA  Parse the generated configuration names.

    nconfigs = numel(data);

    object_idx = (1:nconfigs).';
    object_name = strings(nconfigs, 1);
    orientation_family = strings(nconfigs, 1);
    num_objects = zeros(nconfigs, 1);
    object_type = strings(nconfigs, 1);

    expression = [ ...
        '^(single|same_direction|facing_away_origin|', ...
        'random_orientation|facing_origin)_N(\d+)_(.+)$'];

    for ii = 1:nconfigs
        object_name(ii) = string(data(ii).name);
        tokens = regexp( ...
            char(object_name(ii)), ...
            expression, ...
            'tokens', ...
            'once');

        if isempty(tokens)
            error('Could not parse constellation name "%s".', ...
                object_name(ii))
        end

        orientation_family(ii) = string(tokens{1});
        num_objects(ii) = str2double(tokens{2});
        object_type(ii) = string(tokens{3});
    end

    metadata = table( ...
        object_idx, ...
        object_name, ...
        orientation_family, ...
        num_objects, ...
        object_type);
end

function map = build_single_poles_truth_map( ...
    eval_metadata, library_metadata)
%BUILD_SINGLE_POLES_TRUTH_MAP
%   Map each constellation to the single-object library entry having the
%   same underlying object type.

    neval = height(eval_metadata);
    nlib = height(library_metadata);
    map = false(neval, nlib);

    for ii = 1:neval
        for jj = 1:nlib
            map(ii, jj) = ...
                eval_metadata.object_type(ii) == ...
                library_metadata.object_type(jj);
        end

        if nnz(map(ii, :)) ~= 1
            error([ ...
                'Expected exactly one library match for "%s"; ', ...
                'found %d.'], ...
                eval_metadata.object_name(ii), ...
                nnz(map(ii, :)))
        end
    end
end

function acc = accuracy_by_sigma_table(trials)
%ACCURACY_BY_SIGMA_TABLE  Compute accuracy grouped by noise level.

    sigma = unique(trials.sigma, 'stable');
    accuracy = zeros(numel(sigma), 1);
    ntrials = zeros(numel(sigma), 1);

    for isig = 1:numel(sigma)
        mask = trials.sigma == sigma(isig);
        accuracy(isig) = mean(trials.correct(mask));
        ntrials(isig) = nnz(mask);
    end

    acc = table(sigma, accuracy, ntrials);
end
