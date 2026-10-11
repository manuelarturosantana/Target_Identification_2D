function run_identification_02_residue_comparison()
%RUN_IDENTIFICATION_02_RESIDUE_COMPARISON  Directional fitted-residue comparison.
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
cfg.ridge = 1e-10;
cfg.random_seed = 1;
cfg.show_progress = true;
cfg.residue_floor = 1e-12;

method_id = 2;
method_name = "fitted_residue_comparison";

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
    'identification_02_residue_comparison.mat');
save(output_path, 'results');
fprintf('Saved method 02 results to:\n  %s\n', output_path);
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
%BUILD_POLE_MODELS  Build complex pole bases for every library object.

models = struct([]);

for j = 1:numel(library_data)
    ts = library_data(j).ts;
    late_idx = find(ts >= cfg.late_time, 1, 'first');

    if isempty(late_idx)
        error('build_pole_models:InvalidLateTime', ...
            'cfg.late_time exceeds the final sample for %s.', library_data(j).name);
    end

    tt = ts(late_idx:end);
    Z = pole_basis(tt, library_data(j).poles);
    models(j).late_idx = late_idx;
    models(j).time_window = tt;
    models(j).Z = Z;
    models(j).fit = (Z' * Z + cfg.ridge * eye(size(Z, 2))) \ Z';
end
end

function Z = pole_basis(times, poles)
%POLE_BASIS  Build exp(-1i * phi_k * t) for every pole phi_k.

Z = exp(-1i * times(:) * poles(:).');
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

function residues = fit_residues(f_full, model)
%FIT_RESIDUES  Fit r_k in sum_k r_k exp(-1i * phi_k * t).

f = f_full(model.late_idx:end);
residues = model.fit * f;
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
            clean = evaluation_data(itrue).signals(irec, :).';
            noise = (sigma / sqrt(2)) * ...
                (randn(size(clean)) + 1i * randn(size(clean)));
            f_full = clean + noise;
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
%PREPARE_METHOD  Cache directional library residue vectors.

method_data = repmat(struct( ...
    'residues', [], 'norm_squared', []), numel(library_data), 1);

for j = 1:numel(library_data)
    nreceiver = size(library_data(j).signals, 1);
    npoles = numel(library_data(j).poles);
    fitted_library_residues = zeros(nreceiver, npoles);

    for ireceiver = 1:nreceiver
        fitted_library_residues(ireceiver, :) = fit_residues( ...
            library_data(j).signals(ireceiver, :).', models(j)).';
    end

    % r_res_all contains the spatial scattering residue before the
    % incident-pulse and time-shift factors are applied. Fitting the clean
    % library waveform here puts both sides of the comparison in the same
    % exp(-1i * phi_k * t) coefficient convention, including for older MAT
    % files that do not store their pulse metadata.
    method_data(j).residues = fitted_library_residues;
    method_data(j).norm_squared = ...
        sum(abs(fitted_library_residues).^2, 2) + cfg.residue_floor;
end
end

function scores = score_candidates(f_full, angle, ~, library_data, models, ...
        method_data, cfg)
%SCORE_CANDIDATES  Compare fitted residues with the directional library residues.

scores = zeros(1, numel(library_data));

for j = 1:numel(library_data)
    idx = nearest_angle_index(library_data(j).receiver_angles, angle);
    fitted = fit_residues(f_full, models(j));
    expected = method_data(j).residues(idx, :).';

    % Remove an arbitrary global complex scale, then decide solely from
    % the relative fitted-versus-library residue mismatch.
    scale = (expected' * fitted) / method_data(j).norm_squared(idx);
    mismatch = norm(fitted - scale * expected, 2)^2 ...
        / max(norm(fitted, 2)^2, cfg.residue_floor);
    scores(j) = -real(mismatch);
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
