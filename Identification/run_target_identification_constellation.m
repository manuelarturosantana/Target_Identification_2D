% run_target_identification_constellation.m
%
% Run GLRT target identification for cavity constellations using four methods:
%
%   Method 1:
%       Full objects and full poles.
%       Naive GLRT classification.
%       Permutations are treated as different classes.
%
%   Method 2:
%       Full objects and full poles.
%       Permutations are treated as the same class.
%       Offset and orientation are still part of the class.
%
%   Method 3:
%       Full objects and full poles.
%       Permutations, offset, and orientation are ignored.
%       Only the object multiset matters.
%
%   Method 4:
%       Single-object pole classification.
%       Library contains only single-object poles.
%       Evaluation contains only singles and same-object constellations.

clear
clc

% ---- user-defined parameters ----
% paths
input_folder = '/scratch/vhojas/data_same_object_constellations_2d_fast';

output_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'constellation_four_methods_results.mat');

% classification params
late_time = 80.0;
sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];

% SNR based on representative object
% If snr_object_name is empty, the first object is used.
% If snr_receiver_angle is empty, the first receiver is used.
snr_object_name = 'single_circ_o050';
snr_receiver_angle = 3*pi/2;

% build params struct
params = build_id_params( ...
    'input_folder', input_folder, ...
    'output_path', output_path, ...
    'late_time', late_time, ...
    'noise_sigmas', sigmas, ...
    'snr_object_name', snr_object_name, ...
    'snr_receiver_angle', snr_receiver_angle, ...
    'show_progress', true);

% load objects and spectral data
data_all = load_pole_data(params);

% ============================================================
% METHOD 1: FULL OBJECTS, NAIVE CLASSIFICATION
% ============================================================

fprintf('\n============================================================\n');
fprintf('Method 1: full objects, naive GLRT classification\n');
fprintf('============================================================\n');

model_full_naive = build_glrt_model(data_all, params);

% No map is attached here. This preserves the original classifier behavior:
% correct means pred_idx == true_idx.
results_full_naive = classify_glrt_data( ...
    data_all, ...
    model_full_naive, ...
    params);

results_full_naive.method_name = "full_objects_naive";
results_full_naive.method_description = ...
    "Full object library. Permutations, offset, and orientation are all distinct.";

results_full_naive.accuracy_by_sigma = ...
    accuracy_by_sigma_table(results_full_naive.trials);

% ============================================================
% METHOD 2: REMOVE PERMUTATIONS ONLY
% ============================================================

fprintf('\n============================================================\n');
fprintf('Method 2: full objects, permutations treated as equivalent\n');
fprintf('============================================================\n');

model_remove_permutations = build_glrt_model(data_all, params);

% This map ignores object ordering only.
% It still keeps offset and orientation/family.
map_remove_permutations = build_permutation_only_map(data_all, data_all);

model_remove_permutations = attach_glrt_classification_map( ...
    model_remove_permutations, ...
    data_all, ...
    data_all, ...
    map_remove_permutations);

results_remove_permutations = classify_glrt_data( ...
    data_all, ...
    model_remove_permutations, ...
    params);

results_remove_permutations.method_name = "remove_permutations";
results_remove_permutations.method_description = ...
    "Full object library. Object-order permutations are equivalent. Offset and orientation remain distinct.";

results_remove_permutations.accuracy_by_sigma = ...
    accuracy_by_sigma_table(results_remove_permutations.trials);

% ============================================================
% METHOD 3: REMOVE PERMUTATIONS, OFFSET, AND ORIENTATION
% ============================================================

fprintf('\n============================================================\n');
fprintf('Method 3: full objects, only object multiset matters\n');
fprintf('============================================================\n');

model_object_multiset = build_glrt_model(data_all, params);

% This map ignores:
%   - object order,
%   - aligned vs offset,
%   - same-direction vs facing.
%
% It keeps only:
%   - single vs pair,
%   - circle/ellipse object multiset.
map_object_multiset = build_object_multiset_map(data_all, data_all);

model_object_multiset = attach_glrt_classification_map( ...
    model_object_multiset, ...
    data_all, ...
    data_all, ...
    map_object_multiset);

results_object_multiset = classify_glrt_data( ...
    data_all, ...
    model_object_multiset, ...
    params);

results_object_multiset.method_name = "object_multiset";
results_object_multiset.method_description = ...
    "Full object library. Permutations, offset, and orientation are ignored.";

results_object_multiset.accuracy_by_sigma = ...
    accuracy_by_sigma_table(results_object_multiset.trials);

% ============================================================
% METHOD 4: SINGLE-OBJECT POLE CLASSIFICATION
% ============================================================

fprintf('\n============================================================\n');
fprintf('Method 4: single-object pole classification\n');
fprintf('============================================================\n');

single_library_mask = arrayfun(@is_single_object_config, data_all);
single_eval_mask = arrayfun(@is_single_poles_eval_config, data_all);

data_single_library = data_all(single_library_mask);
data_single_eval = data_all(single_eval_mask);

fprintf('\nSingle-object pole library:\n');
for ii = 1:numel(data_single_library)
    fprintf('  %s\n', data_single_library(ii).name);
end

fprintf('\nSingle-poles evaluation set:\n');
for ii = 1:numel(data_single_eval)
    fprintf('  %s\n', data_single_eval(ii).name);
end

model_single_poles = build_glrt_model(data_single_library, params);

map_single_poles = build_single_poles_truth_map( ...
    data_single_eval, ...
    data_single_library);

model_single_poles = attach_glrt_classification_map( ...
    model_single_poles, ...
    data_single_eval, ...
    data_single_library, ...
    map_single_poles);

results_single_poles = classify_glrt_data( ...
    data_single_eval, ...
    model_single_poles, ...
    params);

results_single_poles.method_name = "single_poles";
results_single_poles.method_description = ...
    "Single-object pole library. Evaluation excludes mixed-object constellations.";

results_single_poles.accuracy_by_sigma = ...
    accuracy_by_sigma_table(results_single_poles.trials);

% ============================================================
% SAVE
% ============================================================

shape_info = collect_shape_info(data_all);

results = struct();
results.params = params;

results.method_full_naive = results_full_naive;
results.method_remove_permutations = results_remove_permutations;
results.method_object_multiset = results_object_multiset;
results.method_single_poles = results_single_poles;

results.shape_info = shape_info;

output_folder = fileparts(params.output_path);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

save(params.output_path, 'results');

fprintf('\nSaved classification results to:\n  %s\n', params.output_path);

% ============================================================
% NEW / CHANGED LOCAL HELPERS
% ============================================================

function map = build_permutation_only_map(eval_data, library_data)
%BUILD_PERMUTATION_ONLY_MAP
%   Truth map for method 2.
%
%   Ignores object order only.
%   Keeps:
%       - single vs pair,
%       - family/orientation,
%       - aligned vs offset.

    neval = numel(eval_data);
    nlib = numel(library_data);

    map = false(neval, nlib);

    for i = 1:neval
        true_key = canonical_permutation_only_key(eval_data(i).name);

        for j = 1:nlib
            pred_key = canonical_permutation_only_key(library_data(j).name);
            map(i, j) = true_key == pred_key;
        end

        if ~any(map(i, :))
            error('No valid permutation-equivalent match found for "%s".', ...
                eval_data(i).name);
        end
    end
end

function map = build_object_multiset_map(eval_data, library_data)
%BUILD_OBJECT_MULTISET_MAP
%   Truth map for method 3.
%
%   Ignores:
%       - object order,
%       - aligned vs offset,
%       - same-direction vs facing.
%
%   Keeps only the object multiset.

    neval = numel(eval_data);
    nlib = numel(library_data);

    map = false(neval, nlib);

    for i = 1:neval
        true_key = canonical_object_multiset_key(eval_data(i).name);

        for j = 1:nlib
            pred_key = canonical_object_multiset_key(library_data(j).name);
            map(i, j) = true_key == pred_key;
        end

        if ~any(map(i, :))
            error('No valid object-multiset match found for "%s".', ...
                eval_data(i).name);
        end
    end
end

function key = canonical_permutation_only_key(name)
%CANONICAL_PERMUTATION_ONLY_KEY
%   Canonical key for method 2.
%
%   Object order is ignored.
%   Offset and orientation/family are retained.

    name = string(name);

    tags = get_object_tags_from_name(name);
    tags = sort(tags);

    if startsWith(name, "single_")
        key = "single_" + strjoin(tags, "_");
        return
    end

    family = extract_family_key(name);
    placement = extract_placement_key(name);

    key = family + "_" + placement + "_" + strjoin(tags, "_");
end

function key = canonical_object_multiset_key(name)
%CANONICAL_OBJECT_MULTISET_KEY
%   Canonical key for method 3.
%
%   Object order, offset, and orientation/family are ignored.

    name = string(name);

    tags = get_object_tags_from_name(name);
    tags = sort(tags);

    if numel(tags) == 1
        key = "single_" + tags(1);
    else
        key = "pair_" + strjoin(tags, "_");
    end
end

function tags = get_object_tags_from_name(name)
%GET_OBJECT_TAGS_FROM_NAME
%   Extract object tags from generated configuration names.

    name = char(string(name));

    expr = '(circ_o050|ell_o075)';
    hits = regexp(name, expr, 'match');

    tags = string(hits(:));

    if isempty(tags)
        error('Could not extract object tag from name "%s".', name);
    end
end

function family = extract_family_key(name)
%EXTRACT_FAMILY_KEY
%   Extract orientation/family from configuration name.

    name = string(name);

    if startsWith(name, "same_dir_same_obj")
        family = "same_dir_same_obj";
    elseif startsWith(name, "same_dir_diff_obj")
        family = "same_dir_diff_obj";
    elseif startsWith(name, "facing_same_obj")
        family = "facing_same_obj";
    elseif startsWith(name, "facing_diff_obj")
        family = "facing_diff_obj";
    else
        family = "unknown_family";
    end
end

function placement = extract_placement_key(name)
%EXTRACT_PLACEMENT_KEY
%   Extract aligned/offset placement from configuration name.

    name = string(name);

    if contains(name, "_aligned_")
        placement = "aligned";
    elseif contains(name, "_offset_")
        placement = "offset";
    elseif startsWith(name, "single_")
        placement = "single";
    else
        placement = "unknown_placement";
    end
end

% ============================================================
% LOCAL HELPERS
% ============================================================

function map = build_constellation_permutation_map(eval_data, library_data)
%BUILD_CONSTELLATION_PERMUTATION_MAP
%   Build truth map for whole-constellation classification.
%
%   Same object-order permutations are counted as equivalent.
%
%   Example:
%       same_dir_diff_obj_aligned_circ_o050_ell_o075
%       same_dir_diff_obj_aligned_ell_o075_circ_o050
%
%   get the same canonical key.

    neval = numel(eval_data);
    nlib = numel(library_data);

    map = false(neval, nlib);

    for i = 1:neval
        true_key = canonical_constellation_key(eval_data(i).name);

        for j = 1:nlib
            pred_key = canonical_constellation_key(library_data(j).name);
            map(i, j) = true_key == pred_key;
        end

        if ~any(map(i, :))
            error('No valid constellation library match found for "%s".', ...
                eval_data(i).name);
        end
    end
end

function map = build_single_poles_truth_map(eval_data, library_data)
%BUILD_SINGLE_POLES_TRUTH_MAP
%   Build truth map for single-poles classification.
%
%   Evaluation objects may be same-object constellations, but predictions
%   are only allowed to be single-object pole-library entries.

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

function tf = is_single_object_config(d)
%IS_SINGLE_OBJECT_CONFIG  True for single-object configurations.

    name = string(d.name);
    tf = startsWith(name, "single_");
end

function tf = is_single_poles_eval_config(d)
%IS_SINGLE_POLES_EVAL_CONFIG
%   True for objects evaluated in the single-poles method.
%
%   Keeps:
%       single_circ_o050
%       single_ell_o075
%       circ+circ constellations
%       ell+ell constellations
%
%   Excludes:
%       circ+ell constellations
%       ell+circ constellations

    name = string(d.name);

    if startsWith(name, "single_")
        tf = true;
        return
    end

    tags = get_object_tags_from_name(name);

    if numel(tags) ~= 2
        tf = false;
        return
    end

    tf = tags(1) == tags(2);
end

function key = canonical_constellation_key(name)
%CANONICAL_CONSTELLATION_KEY
%   Canonical label for whole-constellation classification.
%
%   Object order is ignored, but family and placement are retained.

    name = string(name);

    tags = get_object_tags_from_name(name);
    tags = sort(tags);

    if startsWith(name, "single_")
        key = "single_" + strjoin(tags, "_");
        return
    end

    family = extract_family_key(name);
    placement = extract_placement_key(name);

    key = family + "_" + placement + "_" + strjoin(tags, "_");
end

function key = canonical_single_object_key(name)
%CANONICAL_SINGLE_OBJECT_KEY
%   Canonical label for single-poles classification.
%
%   Same-object constellations are mapped to their single-object type.
%
%   Examples:
%       single_circ_o050
%           -> circ_o050
%
%       same_dir_same_obj_aligned_circ_o050_circ_o050
%           -> circ_o050

    name = string(name);

    tags = get_object_tags_from_name(name);
    tags = unique(tags, 'stable');

    if numel(tags) ~= 1
        error(['Single-poles method received mixed-object label "%s". ', ...
               'Mixed constellations should have been excluded.'], name);
    end

    key = tags(1);
end

function acc = accuracy_by_sigma_table(trials)
%ACCURACY_BY_SIGMA_TABLE  Compute accuracy grouped by sigma.

    sigmas = unique(trials.sigma, 'stable');
    acc = zeros(numel(sigmas), 1);

    for isig = 1:numel(sigmas)
        mask = trials.sigma == sigmas(isig);
        acc(isig) = mean(trials.correct(mask));
    end
end

function shape_info = collect_shape_info(data)
%COLLECT_SHAPE_INFO
%   Collect curve data for plotting.
%
%   This works even if load_pole_data does not load curveX/curveY, because
%   it falls back to loading them directly from data(ii).file_path.

    shape_info = struct([]);

    for ii = 1:numel(data)
        shape_info(ii).name = string(data(ii).name);
        shape_info(ii).file_path = string(data(ii).file_path);

        curveX = [];
        curveY = [];

        if isfield(data, 'curveX') && ~isempty(data(ii).curveX)
            curveX = data(ii).curveX;
        end

        if isfield(data, 'curveY') && ~isempty(data(ii).curveY)
            curveY = data(ii).curveY;
        end

        if (isempty(curveX) || isempty(curveY)) && isfile(data(ii).file_path)
            S = load(data(ii).file_path, 'curveX', 'curveY');

            if isfield(S, 'curveX')
                curveX = S.curveX;
            end

            if isfield(S, 'curveY')
                curveY = S.curveY;
            end
        end

        shape_info(ii).curveX = curveX;
        shape_info(ii).curveY = curveY;
    end
end