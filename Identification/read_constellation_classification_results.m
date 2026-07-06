% read_constellation_classification_results.m
%
% Read and summarize constellation classification results.

clear
clc

% ---- user-defined parameters ----

results_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'single_poles_constellation_results.mat');

% Empty means use all sigmas.
% Example: sigma_to_plot = 1e-4;
sigma_to_plot = [];

make_plots = true;

% ---------------------------------

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

fprintf('\nLoaded results from:\n  %s\n', results_path);

if isfield(results, 'method_name')
    fprintf('\nMethod: %s\n', results.method_name);
end

if isfield(results, 'method_description')
    fprintf('%s\n', results.method_description);
end

fprintf('\nOverall accuracy: %.4f\n', mean(trials.correct));

fprintf('\nAccuracy by sigma:\n');
acc_by_sigma = accuracy_by_sigma_table(trials);
disp(acc_by_sigma);

fprintf('\nAccuracy by true configuration:\n');
acc_by_config = accuracy_by_true_config_table(trials, eval_objects);
disp(acc_by_config);

fprintf('\nAccuracy by object type:\n');
acc_by_type = accuracy_by_object_type_table(trials, eval_objects);
disp(acc_by_type);

fprintf('\nAccuracy by number of objects:\n');
acc_by_nobj = accuracy_by_num_objects_table(trials, eval_objects);
disp(acc_by_nobj);

fprintf('\nPrediction counts:\n');
pred_counts = prediction_count_table(trials, library_objects);
disp(pred_counts);

if ~isempty(sigma_to_plot)
    mask = trials.sigma == sigma_to_plot;

    if ~any(mask)
        error('Requested sigma_to_plot = %g was not found.', sigma_to_plot);
    end
else
    mask = true(height(trials), 1);
end

fprintf('\nConfusion table');
if ~isempty(sigma_to_plot)
    fprintf(' for sigma = %g', sigma_to_plot);
end
fprintf(':\n');

conf = confusion_table(trials(mask, :), eval_objects, library_objects);
disp(conf);

if make_plots
    plot_accuracy_by_sigma(acc_by_sigma);
    plot_accuracy_by_num_objects(acc_by_nobj);
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

function acc = accuracy_by_true_config_table(trials, eval_objects)
%ACCURACY_BY_TRUE_CONFIG_TABLE  Compute accuracy grouped by true config.

    true_ids = unique(trials.true_idx, 'stable');

    true_name = strings(numel(true_ids), 1);
    accuracy = zeros(numel(true_ids), 1);
    ntrials = zeros(numel(true_ids), 1);

    for ii = 1:numel(true_ids)
        id = true_ids(ii);
        mask = trials.true_idx == id;

        true_name(ii) = get_eval_name(eval_objects, id);
        accuracy(ii) = mean(trials.correct(mask));
        ntrials(ii) = nnz(mask);
    end

    acc = table( ...
        true_ids, ...
        true_name, ...
        accuracy, ...
        ntrials, ...
        'VariableNames', { ...
            'true_idx', ...
            'true_name', ...
            'accuracy', ...
            'ntrials'});
end

function acc = accuracy_by_object_type_table(trials, eval_objects)
%ACCURACY_BY_OBJECT_TYPE_TABLE  Compute accuracy grouped by circ/ell type.

    true_type = strings(height(trials), 1);

    for ii = 1:height(trials)
        name = get_eval_name(eval_objects, trials.true_idx(ii));
        true_type(ii) = get_object_type_from_name(name);
    end

    types = unique(true_type, 'stable');

    accuracy = zeros(numel(types), 1);
    ntrials = zeros(numel(types), 1);

    for itype = 1:numel(types)
        mask = true_type == types(itype);
        accuracy(itype) = mean(trials.correct(mask));
        ntrials(itype) = nnz(mask);
    end

    acc = table( ...
        types, ...
        accuracy, ...
        ntrials, ...
        'VariableNames', { ...
            'object_type', ...
            'accuracy', ...
            'ntrials'});
end

function acc = accuracy_by_num_objects_table(trials, eval_objects)
%ACCURACY_BY_NUM_OBJECTS_TABLE  Compute accuracy grouped by constellation size.

    num_objects = zeros(height(trials), 1);

    for ii = 1:height(trials)
        name = get_eval_name(eval_objects, trials.true_idx(ii));
        num_objects(ii) = get_num_objects_from_name(name);
    end

    nvals = unique(num_objects, 'stable');

    accuracy = zeros(numel(nvals), 1);
    ntrials = zeros(numel(nvals), 1);

    for ii = 1:numel(nvals)
        mask = num_objects == nvals(ii);
        accuracy(ii) = mean(trials.correct(mask));
        ntrials(ii) = nnz(mask);
    end

    acc = table( ...
        nvals, ...
        accuracy, ...
        ntrials, ...
        'VariableNames', { ...
            'num_objects', ...
            'accuracy', ...
            'ntrials'});
end

function counts = prediction_count_table(trials, library_objects)
%PREDICTION_COUNT_TABLE  Count predictions by library object.

    pred_ids = unique(trials.pred_idx, 'stable');

    pred_name = strings(numel(pred_ids), 1);
    count = zeros(numel(pred_ids), 1);
    fraction = zeros(numel(pred_ids), 1);

    for ii = 1:numel(pred_ids)
        id = pred_ids(ii);
        mask = trials.pred_idx == id;

        pred_name(ii) = get_library_name(library_objects, id);
        count(ii) = nnz(mask);
        fraction(ii) = count(ii) / height(trials);
    end

    counts = table( ...
        pred_ids, ...
        pred_name, ...
        count, ...
        fraction, ...
        'VariableNames', { ...
            'pred_idx', ...
            'pred_name', ...
            'count', ...
            'fraction'});
end

function conf = confusion_table(trials, eval_objects, library_objects)
%CONFUSION_TABLE  Count true configuration versus predicted library class.

    true_ids = unique(trials.true_idx, 'stable');
    pred_ids = unique(trials.pred_idx, 'stable');

    true_names = strings(numel(true_ids), 1);
    pred_names = strings(numel(pred_ids), 1);

    for ii = 1:numel(true_ids)
        true_names(ii) = get_eval_name(eval_objects, true_ids(ii));
    end

    for jj = 1:numel(pred_ids)
        pred_names(jj) = get_library_name(library_objects, pred_ids(jj));
    end

    counts = zeros(numel(true_ids), numel(pred_ids));

    for ii = 1:numel(true_ids)
        for jj = 1:numel(pred_ids)
            counts(ii, jj) = nnz( ...
                trials.true_idx == true_ids(ii) & ...
                trials.pred_idx == pred_ids(jj));
        end
    end

    conf = array2table(counts, ...
        'VariableNames', matlab.lang.makeValidName(pred_names), ...
        'RowNames', cellstr(true_names));
end

function plot_accuracy_by_sigma(acc_by_sigma)
%PLOT_ACCURACY_BY_SIGMA  Plot accuracy versus sigma.
%
%   sigma = 0 is plotted manually one log-step to the left of the smallest
%   positive sigma. Its y-value is the actual accuracy for sigma = 0.

    sigmas = acc_by_sigma.sigma(:);
    accuracy = acc_by_sigma.accuracy(:);

    zero_mask = sigmas == 0;
    pos_mask = sigmas > 0;

    sigmas_pos = sigmas(pos_mask);
    accuracy_pos = accuracy(pos_mask);

    [sigmas_pos, idx] = sort(sigmas_pos);
    accuracy_pos = accuracy_pos(idx);

    if any(zero_mask)
        accuracy_zero = accuracy(find(zero_mask, 1, 'first'));

        if numel(sigmas_pos) < 2
            error('Need at least two positive sigma values to place sigma = 0.');
        end

        % Place sigma = 0 one equal log-step to the left.
        sigma_zero_plot = sigmas_pos(1)^2 / sigmas_pos(2);

        xplot = [sigma_zero_plot; sigmas_pos];
        yplot = [accuracy_zero; accuracy_pos];

        xticks_vals = xplot;
        xticks_labels = ["0"; string(sigmas_pos)];
    else
        xplot = sigmas_pos;
        yplot = accuracy_pos;

        xticks_vals = xplot;
        xticks_labels = string(sigmas_pos);
    end

    figure
    semilogx(xplot, yplot, '-o')
    grid on
    xlabel('\sigma')
    ylabel('accuracy')
    title('Classification accuracy by noise level')

    xticks(xticks_vals)
    xticklabels(xticks_labels)
end

function plot_accuracy_by_num_objects(acc_by_nobj)
%PLOT_ACCURACY_BY_NUM_OBJECTS  Plot accuracy versus constellation size.

    figure
    plot(acc_by_nobj.num_objects, acc_by_nobj.accuracy, '-o')
    grid on
    xlabel('number of objects')
    ylabel('accuracy')
    title('Classification accuracy by constellation size')
end

function name = get_eval_name(eval_objects, idx)
%GET_EVAL_NAME  Return evaluation object name for either result format.

    if any(strcmp(eval_objects.Properties.VariableNames, 'eval_object_name'))
        name = string(eval_objects.eval_object_name(idx));
    elseif any(strcmp(eval_objects.Properties.VariableNames, 'object_name'))
        name = string(eval_objects.object_name(idx));
    else
        error('Could not find evaluation object name column.');
    end
end

function name = get_library_name(library_objects, idx)
%GET_LIBRARY_NAME  Return library object name.

    if any(strcmp(library_objects.Properties.VariableNames, 'object_name'))
        name = string(library_objects.object_name(idx));
    else
        error('Could not find library object name column.');
    end
end

function type = get_object_type_from_name(name)
%GET_OBJECT_TYPE_FROM_NAME  Extract circ/ell type from generated names.

    name = string(name);

    if contains(name, "circ_o050")
        type = "circ_o050";
    elseif contains(name, "ell_o075")
        type = "ell_o075";
    else
        type = "unknown";
    end
end

function nobj = get_num_objects_from_name(name)
%GET_NUM_OBJECTS_FROM_NAME  Extract N from generated constellation names.

    name = char(string(name));

    tokens = regexp(name, '_N(\d+)_', 'tokens', 'once');

    if isempty(tokens)
        nobj = NaN;
    else
        nobj = str2double(tokens{1});
    end
end