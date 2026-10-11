function comparison = read_identification_results()
%READ_IDENTIFICATION_RESULTS  Load and compare all numbered method results.

% Plot selection
% Leave empty to plot every method.
% Example: plot_method_ids = [0, 2];
plot_method_ids = [0, 1, 2];

results_folder = fullfile( ...
    'Identification', 'classification_results', 'angle_upgrades');
output_path = fullfile( ...
    results_folder, 'identification_comparison.mat');
figure_path = fullfile( ...
    results_folder, 'identification_comparison.pdf');

files = dir(fullfile(results_folder, 'identification_*.mat'));
allowed_files = [ ...
    "identification_00_baseline.mat", ...
    "identification_01_hard_visibility_filter.mat", ...
    "identification_02_residue_comparison.mat"];
valid_name = ismember(string({files.name}), allowed_files);
files = files(valid_name);

if isempty(files)
    error('read_identification_results:NoResults', ...
        'No numbered result files were found in %s.', results_folder);
end

records = struct([]);
all_sigmas = [];

for i = 1:numel(files)
    loaded = load(fullfile(files(i).folder, files(i).name));

    if ~isfield(loaded, 'results')
        warning('Skipping %s because it has no results variable.', ...
            files(i).name);
        continue
    end

    result = loaded.results;

    if ~ismember(result.method_id, [0, 1, 2])
        continue
    end

    records(end + 1).method_id = result.method_id;
    records(end).method_name = string(result.method_name);
    records(end).file_name = string(files(i).name);
    records(end).accuracy = mean(result.trials.correct);
    records(end).trials = result.trials;

    all_sigmas = union(all_sigmas, result.trials.sigma);
end

if isempty(records)
    error('read_identification_results:NoValidResults', ...
        'No valid result structures were found.');
end

[~, order] = sort([records.method_id]);
records = records(order);

available_method_ids = [records.method_id];

if isempty(plot_method_ids)
    plot_indices = 1:numel(records);
else
    plot_method_ids = unique(plot_method_ids(:).', 'stable');
    [found, plot_indices] = ismember( ...
        plot_method_ids, available_method_ids);

    if any(~found)
        warning('read_identification_results:UnknownPlotMethod', ...
            'Ignoring unavailable method IDs: %s.', ...
            mat2str(plot_method_ids(~found)));
    end

    plot_indices = plot_indices(found);

    if isempty(plot_indices)
        error('read_identification_results:NoPlotMethods', ...
            'None of the selected method IDs have result files.');
    end
end

method_id = [records.method_id].';
method_name = string({records.method_name}).';
file_name = string({records.file_name}).';
overall_accuracy = [records.accuracy].';

summary = table( ...
    method_id, method_name, file_name, overall_accuracy);

all_sigmas = all_sigmas(:);
accuracy = NaN(numel(all_sigmas), numel(records));
trial_count = zeros(numel(all_sigmas), numel(records));

for j = 1:numel(records)
    trials = records(j).trials;

    for i = 1:numel(all_sigmas)
        rows = trials.sigma == all_sigmas(i);
        trial_count(i, j) = nnz(rows);

        if any(rows)
            accuracy(i, j) = mean(trials.correct(rows));
        end
    end
end

accuracy_by_sigma = array2table(accuracy, ...
    'VariableNames', matlab.lang.makeValidName(method_name));
accuracy_by_sigma = addvars( ...
    accuracy_by_sigma, all_sigmas, ...
    'Before', 1, 'NewVariableNames', 'sigma');

count_by_sigma = array2table(trial_count, ...
    'VariableNames', matlab.lang.makeValidName(method_name));
count_by_sigma = addvars( ...
    count_by_sigma, all_sigmas, ...
    'Before', 1, 'NewVariableNames', 'sigma');

comparison.summary = summary;
comparison.accuracy_by_sigma = accuracy_by_sigma;
comparison.count_by_sigma = count_by_sigma;
comparison.plotted_method_ids = ...
    [records(plot_indices).method_id].';

disp(summary)
disp(accuracy_by_sigma)

% Use 1e-7 only as the plotting position for sigma = 0.
zero_plot_sigma = 1e-7;
plot_sigmas = all_sigmas;
plot_sigmas(all_sigmas == 0) = zero_plot_sigma;

figure('Color', 'w');
hold on

for j = plot_indices
    semilogx(plot_sigmas, accuracy(:, j), '-o', ...
        'DisplayName', sprintf('%02d %s', ...
        records(j).method_id, records(j).method_name));
end

xlabel('Noise standard deviation')
ylabel('Classification accuracy')
ylim([0, 1])

ax = gca;
set(ax, 'XScale', 'log');

% Include the artificial position and label it as sigma = 0.
tick_values = unique([zero_plot_sigma, get(ax, 'XTick')]);
set(ax, 'XTick', tick_values);

tick_labels = string(get(ax, 'XTickLabel'));
tick_labels(tick_values == zero_plot_sigma) = "0";
set(ax, 'XTickLabel', tick_labels);

grid on
legend('Location', 'eastoutside', 'Interpreter', 'none')
title('Comparison of selected identification methods')

exportgraphics(gcf, figure_path, 'ContentType', 'vector');

save(output_path, 'comparison');

fprintf('Saved comparison to:\n  %s\n', output_path);
fprintf('Saved comparison figure to:\n  %s\n', figure_path);
end
