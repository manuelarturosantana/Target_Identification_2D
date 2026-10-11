function run_identification_01_hard_visibility_filter()
%RUN_IDENTIFICATION_01_HARD_VISIBILITY_FILTER  Hard angle-dependent high-SNR candidate filter.
%
% This file is self-contained except for the input MAT files. Each MAT file
% must contain config_name, pols, r_res_all, receiver_angles, ts, and uff_all.

clc

cfg.library_folder = 'local/data_deformation_experiment_2d';
cfg.evaluation_folder = cfg.library_folder;
cfg.library_pattern = '*.mat';
cfg.evaluation_pattern = '*.mat';
cfg.truth_map_path = '';
cfg.output_folder = fullfile('Identification', 'classification_results', 'angle_upgrades');
cfg.late_time = 80.0;
cfg.noise_sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
cfg.svd_tol = 1e-10;
cfg.ridge = 1e-10;
cfg.random_seed = 1;
cfg.show_progress = true;
cfg.visibility_fraction = 0.10;

method_id = 1;
method_name = "hard_visibility_filter";

[evaluation_data, library_data, truth_map] = load_experiment(cfg);
models = build_pole_models(library_data, cfg);
method_data = prepare_method(library_data, models, cfg);
results = classify_method(evaluation_data, library_data, models, ...
    method_data, truth_map, cfg);

results.method_id = method_id;
results.method_name = method_name;
results.configuration = cfg;

if ~exist(cfg.output_folder, 'dir')
    mkdir(cfg.output_folder);
end

output_path = fullfile(cfg.output_folder, ...
    'identification_01_hard_visibility_filter.mat');
save(output_path, 'results');
fprintf('Saved method 01 results to:\n  %s\n', output_path);
end


function [evaluation_data, library_data, truth_map] = load_experiment(cfg)
%LOAD_EXPERIMENT  Load library and evaluation objects and construct truth map.

library_data = load_object_folder(cfg.library_folder, cfg.library_pattern);

if strcmp(cfg.library_folder, cfg.evaluation_folder) ...
        && strcmp(cfg.library_pattern, cfg.evaluation_pattern)
    evaluation_data = library_data;
else
    evaluation_data = load_object_folder(cfg.evaluation_folder, cfg.evaluation_pattern);
end

nlib = numel(library_data);
neval = numel(evaluation_data);

if strlength(string(cfg.truth_map_path)) > 0
    map_file = load(cfg.truth_map_path);

    if ~isfield(map_file, 'true_to_library_map')
        error('load_experiment:MissingTruthMap', ...
            'Truth-map file must contain true_to_library_map.');
    end

    truth_map = logical(map_file.true_to_library_map);
else
    truth_map = false(neval, nlib);

    for i = 1:neval
        for j = 1:nlib
            truth_map(i, j) = evaluation_data(i).name == library_data(j).name;
        end
    end

    if any(~any(truth_map, 2))
        same_files = neval == nlib && all([evaluation_data.file_name] == [library_data.file_name]);

        if same_files
            truth_map = eye(neval, nlib) > 0;
        else
            error('load_experiment:UnmappedEvaluationObject', ...
                ['At least one evaluation config_name has no library match. ', ...
                 'Set cfg.truth_map_path to a MAT file containing ', ...
                 'true_to_library_map with size neval-by-nlib.']);
        end
    end
end

if ~isequal(size(truth_map), [neval, nlib])
    error('load_experiment:BadTruthMapSize', ...
        'true_to_library_map must have size %d-by-%d.', neval, nlib);
end
end

function data = load_object_folder(folder, pattern)
%LOAD_OBJECT_FOLDER  Load the required fields from every matching MAT file.

files = dir(fullfile(folder, pattern));

if isempty(files)
    error('load_object_folder:NoFiles', ...
        'No files matching %s were found in %s.', pattern, folder);
end

[~, order] = sort(string({files.name}));
files = files(order);
data = struct([]);

required = {'config_name', 'pols', 'r_res_all', ...
    'receiver_angles', 'ts', 'uff_all'};

for j = 1:numel(files)
    file_path = fullfile(files(j).folder, files(j).name);
    raw = load(file_path);

    for k = 1:numel(required)
        if ~isfield(raw, required{k})
            error('load_object_folder:MissingField', ...
                '%s does not contain required field %s.', ...
                file_path, required{k});
        end
    end

    signals = raw.uff_all;
    angles = raw.receiver_angles(:);
    poles = raw.pols(:);
    residues = raw.r_res_all;
    ts = raw.ts(:);

    if size(signals, 1) ~= numel(angles)
        error('load_object_folder:AngleCountMismatch', ...
            '%s has %d signals but %d receiver angles.', ...
            file_path, size(signals, 1), numel(angles));
    end

    if size(signals, 2) ~= numel(ts)
        error('load_object_folder:TimeCountMismatch', ...
            '%s has %d time samples in uff_all but %d entries in ts.', ...
            file_path, size(signals, 2), numel(ts));
    end

    if ~isequal(size(residues), [numel(angles), numel(poles)])
        error('load_object_folder:ResidueSizeMismatch', ...
            '%s must have size(r_res_all) = nreceiver-by-npole.', file_path);
    end

    data(j).name = string(raw.config_name);
    data(j).file_name = string(files(j).name);
    data(j).file_path = string(file_path);
    data(j).poles = poles;
    data(j).residues = residues;
    data(j).receiver_angles = angles;
    data(j).ts = ts;
    data(j).signals = signals;
end
end

function models = build_pole_models(library_data, cfg)
%BUILD_POLE_MODELS  Build real cosine-sine pole bases for every library object.

models = struct([]);

for j = 1:numel(library_data)
    ts = library_data(j).ts;
    late_idx = find(ts >= cfg.late_time, 1, 'first');

    if isempty(late_idx)
        error('build_pole_models:InvalidLateTime', ...
            'cfg.late_time exceeds the final sample for %s.', library_data(j).name);
    end

    tt = ts(late_idx:end);
    Z = real_pole_basis(tt, library_data(j).poles);
    [U, S, ~] = svd(Z, 'econ');
    singular_values = diag(S);

    if isempty(singular_values)
        keep = false(size(singular_values));
    else
        keep = singular_values > cfg.svd_tol * max(singular_values);
    end

    models(j).late_idx = late_idx;
    models(j).time_window = tt;
    models(j).Z = Z;
    models(j).U = U(:, keep);
    models(j).rank = nnz(keep);
    models(j).fit = (Z' * Z + cfg.ridge * eye(size(Z, 2))) \ Z';
end
end

function Z = real_pole_basis(times, poles)
%REAL_POLE_BASIS  Build the real SEM basis for complex-frequency poles.

poles = poles(:).';
decay = imag(poles);
oscillation = real(poles);
envelope = exp(times(:) * decay);
Z = [envelope .* cos(times(:) * oscillation), ...
     envelope .* sin(times(:) * oscillation)];
end

function q = residue_signal(times, poles, residues)
%RESIDUE_SIGNAL  Synthesize a real late-time signal from poles and residues.

Z = real_pole_basis(times, poles);
coefficients = [real(residues(:)); imag(residues(:))];
q = Z * coefficients;
end

function score = projection_score(f_full, model)
%PROJECTION_SCORE  Return the normalized pole-subspace projection energy.

f = real(f_full(model.late_idx:end));
score = norm(model.U' * f, 2)^2 / max(norm(f, 2)^2, eps);
end

function idx = nearest_angle_index(angles, angle)
%NEAREST_ANGLE_INDEX  Find the nearest direction on the periodic circle.

distance = abs(wrap_angle(angles(:) - angle));
[~, idx] = min(distance);
end

function angle = wrap_angle(angle)
%WRAP_ANGLE  Map angles to [-pi, pi).

angle = mod(angle + pi, 2 * pi) - pi;
end

function energy = clean_late_energy(object, late_idx)
%CLEAN_LATE_ENERGY  Compute clean late-time energy at every receiver angle.

signals = real(object.signals(:, late_idx:end));
energy = sum(signals.^2, 2);
end

function [magnitude, coefficients] = fitted_mode_magnitudes(f_full, model)
%FITTED_MODE_MAGNITUDES  Fit real modal coefficients and return pair magnitudes.

f = real(f_full(model.late_idx:end));
coefficients = model.fit * f;
npoles = numel(coefficients) / 2;
magnitude = hypot(coefficients(1:npoles), coefficients(npoles + 1:end));
end

function vector = unit_vector(vector)
%UNIT_VECTOR  Normalize a vector without producing NaNs at zero norm.

vector = vector(:);
vector = vector / max(norm(vector, 2), eps);
end

function value = safe_log(value)
%SAFE_LOG  Compute a finite logarithm for nonnegative numerical quantities.

value = log(max(real(value), realmin('double')));
end


function results = classify_method(evaluation_data, library_data, models, ...
        method_data, truth_map, cfg)
%CLASSIFY_METHOD  Classify every single-look evaluation signal.

rng(cfg.random_seed);
neval = numel(evaluation_data);
nlib = numel(library_data);
sigmas = cfg.noise_sigmas(:);
ntrials_per_sigma = sum(arrayfun(@(x) size(x.signals, 1), evaluation_data));
ntrials = numel(sigmas) * ntrials_per_sigma;

sigma_col = zeros(ntrials, 1);
true_idx = zeros(ntrials, 1);
pred_idx = zeros(ntrials, 1);
receiver_idx = zeros(ntrials, 1);
receiver_angle = zeros(ntrials, 1);
score_max = zeros(ntrials, 1);
score_second = zeros(ntrials, 1);
correct = false(ntrials, 1);
scores_all = zeros(ntrials, nlib);

row = 0;
start_time = tic;
last_report_time = tic;

for isig = 1:numel(sigmas)
    sigma = sigmas(isig);

    for itrue = 1:neval
        for irec = 1:size(evaluation_data(itrue).signals, 1)
            row = row + 1;
            clean = real(evaluation_data(itrue).signals(irec, :)).';
            f_full = clean + sigma * randn(size(clean));
            angle = evaluation_data(itrue).receiver_angles(irec);
            scores = score_candidates(f_full, angle, sigma, library_data, ...
                models, method_data, cfg);

            [ordered, order] = sort(scores, 'descend');
            pred = order(1);

            sigma_col(row) = sigma;
            true_idx(row) = itrue;
            pred_idx(row) = pred;
            receiver_idx(row) = irec;
            receiver_angle(row) = angle;
            score_max(row) = ordered(1);

            if numel(ordered) >= 2
                score_second(row) = ordered(2);
            else
                score_second(row) = NaN;
            end

            correct(row) = truth_map(itrue, pred);
            scores_all(row, :) = scores;

            if cfg.show_progress && ...
                    (toc(last_report_time) > 0.5 || row == ntrials)
                fprintf('\r  %6d / %6d  (%6.2f%%)  elapsed: %7.1fs', ...
                    row, ntrials, 100 * row / ntrials, toc(start_time));
                last_report_time = tic;
            end
        end
    end
end

if cfg.show_progress
    fprintf('\n');
end

trials = table(sigma_col, true_idx, pred_idx, receiver_idx, ...
    receiver_angle, score_max, score_second, correct, ...
    'VariableNames', {'sigma', 'true_idx', 'pred_idx', 'receiver_idx', ...
    'receiver_angle', 'score_max', 'score_second', 'correct'});

results = package_results(trials, scores_all, evaluation_data, library_data);
end


function method_data = prepare_method(library_data, models, cfg)
%PREPARE_METHOD  Mark high-visibility aspects for every candidate.

for j = 1:numel(library_data)
    energy = clean_late_energy(library_data(j), models(j).late_idx);
    method_data(j).energy = energy;
    method_data(j).high = energy >= cfg.visibility_fraction * max(energy);
end
end

function scores = score_candidates(f_full, angle, ~, library_data, models, ...
        method_data, ~)
%SCORE_CANDIDATES  Discard candidates expected to have low directional SNR.

scores = zeros(1, numel(library_data));
high = false(1, numel(library_data));

for j = 1:numel(library_data)
    scores(j) = projection_score(f_full, models(j));
    idx = nearest_angle_index(library_data(j).receiver_angles, angle);
    high(j) = method_data(j).high(idx);
end

if any(high)
    scores(~high) = -Inf;
end
end


function results = package_results(trials, scores, evaluation_data, library_data)
%PACKAGE_RESULTS  Build a common result structure for every method.

nlib = numel(library_data);
object_idx = (1:nlib).';
object_name = strings(nlib, 1);
file_name = strings(nlib, 1);
file_path = strings(nlib, 1);

for j = 1:nlib
    object_name(j) = library_data(j).name;
    file_name(j) = library_data(j).file_name;
    file_path(j) = library_data(j).file_path;
end

objects = table(object_idx, object_name, file_name, file_path);

neval = numel(evaluation_data);
eval_idx = (1:neval).';
eval_object_name = strings(neval, 1);
eval_file_name = strings(neval, 1);
eval_file_path = strings(neval, 1);

for j = 1:neval
    eval_object_name(j) = evaluation_data(j).name;
    eval_file_name(j) = evaluation_data(j).file_name;
    eval_file_path(j) = evaluation_data(j).file_path;
end

eval_objects = table(eval_idx, eval_object_name, eval_file_name, eval_file_path);

results.trials = trials;
results.objects = objects;
results.eval_objects = eval_objects;
results.scores = scores;
results.accuracy = mean(trials.correct);
end
