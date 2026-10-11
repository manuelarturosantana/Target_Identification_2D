function comparison = read_plane_identification_partial_circumference()
%READ_PLANE_IDENTIFICATION_PARTIAL_CIRCUMFERENCE  Plot accuracy for each X%.

results_path = fullfile('Identification', 'classification_results', ...
    'plane_partial_circumference_results.mat');
figure_path = fullfile('Identification', 'classification_results', ...
    'plane_partial_circumference_comparison.pdf');
comparison_path = fullfile('Identification', 'classification_results', ...
    'plane_partial_circumference_comparison.mat');

loaded = load(results_path, 'partial_results');

if ~isfield(loaded, 'partial_results') ...
        || ~isfield(loaded.partial_results, 'runs') ...
        || isempty(loaded.partial_results.runs)
    error('read_plane_identification_partial_circumference:BadResults', ...
        '%s does not contain any partial-circumference runs.', results_path);
end

runs = loaded.partial_results.runs;
percentages = [runs.percentage].';
[percentages, order] = sort(percentages);
runs = runs(order);
all_sigmas = [];

for irun = 1:numel(runs)
    all_sigmas = union(all_sigmas, runs(irun).results.trials.sigma);
end

all_sigmas = all_sigmas(:);
accuracy = NaN(numel(all_sigmas), numel(runs));
trial_count = zeros(numel(all_sigmas), numel(runs));
snr_full_by_run = NaN(numel(all_sigmas), numel(runs));
snr_late_by_run = NaN(numel(all_sigmas), numel(runs));
zero_angle_by_run = NaN(numel(runs), 1);

for irun = 1:numel(runs)
    trials = runs(irun).results.trials;
    [run_sigmas, run_snr_full, run_snr_late, zero_angle_by_run(irun)] = ...
        representative_zero_angle_snr(trials);

    for isigma = 1:numel(all_sigmas)
        rows = trials.sigma == all_sigmas(isigma);
        trial_count(isigma, irun) = nnz(rows);

        if any(rows)
            accuracy(isigma, irun) = mean(trials.correct(rows));
        end


        snr_row = find(run_sigmas == all_sigmas(isigma), 1);

        if ~isempty(snr_row)
            snr_full_by_run(isigma, irun) = run_snr_full(snr_row);
            snr_late_by_run(isigma, irun) = run_snr_late(snr_row);
        end
    end
end


assert_consistent_snr(snr_full_by_run, 'full-time');
assert_consistent_snr(snr_late_by_run, 'late-time');
representative_snr_full = snr_full_by_run(:, 1);
representative_snr_late = snr_late_by_run(:, 1);

column_names = strings(1, numel(percentages));

for irun = 1:numel(percentages)
    column_names(irun) = sprintf('X%g_percent', percentages(irun));
end

column_names = matlab.lang.makeValidName(column_names);
accuracy_by_sigma = array2table(accuracy, ...
    'VariableNames', cellstr(column_names));
accuracy_by_sigma = addvars(accuracy_by_sigma, all_sigmas, ...
    'Before', 1, 'NewVariableNames', 'sigma');

count_by_sigma = array2table(trial_count, ...
    'VariableNames', cellstr(column_names));
count_by_sigma = addvars(count_by_sigma, all_sigmas, ...
    'Before', 1, 'NewVariableNames', 'sigma');

angle_interval = vertcat(runs.angle_interval);
receiver_count_per_object = vertcat(runs.receiver_count_by_object);
summary = table(percentages, angle_interval(:, 1), angle_interval(:, 2), ...
    receiver_count_per_object, ...
    'VariableNames', {'percentage', 'angle_min', 'angle_max', ...
    'receiver_count_by_object'});

comparison = struct();
comparison.summary = summary;
comparison.accuracy_by_sigma = accuracy_by_sigma;
comparison.count_by_sigma = count_by_sigma;
comparison.representative_snr_by_sigma = table( ...
    all_sigmas, representative_snr_full, representative_snr_late, ...
    'VariableNames', {'sigma', 'snr_full_db', 'snr_late_db'});
comparison.representative_receiver_angle = zero_angle_by_run(1);

disp(summary)
disp(accuracy_by_sigma)

positive_sigmas = all_sigmas(all_sigmas > 0);

if isempty(positive_sigmas)
    zero_plot_sigma = 1e-7;
else
    zero_plot_sigma = min(positive_sigmas) / 10;
end

plot_sigmas = all_sigmas;
plot_sigmas(all_sigmas == 0) = zero_plot_sigma;

figure('Color', 'w');
ax_sigma = axes('Position', [0.11, 0.32, 0.66, 0.58]);
hold(ax_sigma, 'on')

for irun = 1:numel(runs)
    semilogx(ax_sigma, plot_sigmas, accuracy(:, irun), '-o', ...
        'DisplayName', sprintf('%g%% circumference', percentages(irun)));
end

xlabel(ax_sigma, 'Noise standard deviation, \sigma')
ylabel(ax_sigma, 'Classification accuracy')
ylim(ax_sigma, [0, 1])
grid(ax_sigma, 'on')
legend(ax_sigma, 'Location', 'eastoutside')
title(ax_sigma, 'Plane identification by retained observation circumference')

% The outside legend changes the main axes width lazily. Resolve that
% layout before copying the main x-axis geometry to the two SNR axes.
drawnow

ax_sigma.XScale = 'log';
ax_sigma.XTick = plot_sigmas;
ax_sigma.XTickLabel = sigma_tick_labels(all_sigmas);

ax_full = bottom_snr_axis(ax_sigma, plot_sigmas, ...
    representative_snr_full, ...
    'Representative full-time SNR at 0^\circ [dB]', 0.13);
ax_late = bottom_snr_axis(ax_sigma, plot_sigmas, ...
    representative_snr_late, ...
    'Representative late-time SNR at 0^\circ [dB]', 0.26);
linkaxes([ax_sigma, ax_full, ax_late], 'x')

% Keep all three scales exactly aligned even if labels affect layout.
drawnow
main_position = ax_sigma.Position;
ax_full.Position([1, 3]) = main_position([1, 3]);
ax_late.Position([1, 3]) = main_position([1, 3]);

exportgraphics(gcf, figure_path, 'ContentType', 'vector');
save(comparison_path, 'comparison');

fprintf('Saved comparison to:\n  %s\n', comparison_path);
fprintf('Saved comparison figure to:\n  %s\n', figure_path);
end


function [sigmas, snr_full, snr_late, selected_angle] = ...
        representative_zero_angle_snr(trials)
%REPRESENTATIVE_ZERO_ANGLE_SNR  Average class SNRs at the receiver nearest 0.

required = {'sigma', 'receiver_angle', 'snr_full_db', 'snr_late_db'};

if ~all(ismember(required, trials.Properties.VariableNames))
    error('read_plane_identification_partial_circumference:MissingTrialSNR', ...
        ['Trial data is missing receiver angle or SNR columns. Re-run ', ...
         'run_plane_identification_partial_circumference first.']);
end

angle_distance = abs(angle(exp(1i * trials.receiver_angle)));
minimum_distance = min(angle_distance);
angle_tolerance = max(100 * eps(pi), 100 * eps(minimum_distance));
angle_rows = abs(angle_distance - minimum_distance) <= angle_tolerance;
selected_angles = trials.receiver_angle(angle_rows);
selected_angle = angle(exp(1i * selected_angles(1)));
sigmas = unique(trials.sigma);
snr_full = NaN(size(sigmas));
snr_late = NaN(size(sigmas));

for isigma = 1:numel(sigmas)
    rows = angle_rows & trials.sigma == sigmas(isigma);
    snr_full(isigma) = mean(trials.snr_full_db(rows));
    snr_late(isigma) = mean(trials.snr_late_db(rows));
end
end


function assert_consistent_snr(snr_by_run, description)
%ASSERT_CONSISTENT_SNR  The 0-degree signal must be common to every aperture.

reference = snr_by_run(:, 1);
matches = (snr_by_run == reference) | ...
    (isinf(snr_by_run) & isinf(reference));

if ~all(matches, 'all')
    error('read_plane_identification_partial_circumference:InconsistentSNR', ...
        'Representative %s SNR differs between aperture runs.', description);
end
end


function labels = sigma_tick_labels(sigmas)
%SIGMA_TICK_LABELS  Format the fixed noise levels, including sigma = 0.

labels = strings(size(sigmas));

for isigma = 1:numel(sigmas)
    if sigmas(isigma) == 0
        labels(isigma) = "0";
    else
        labels(isigma) = sprintf('%.0e', sigmas(isigma));
    end
end
end


function ax = bottom_snr_axis(base_axis, x_ticks, snr_values, label, offset)
%BOTTOM_SNR_AXIS  Add a negligible-height x-axis below the main plot.

position = base_axis.Position;
position = [position(1), position(2) - offset, position(3), 1e-12];
ax = axes( ...
    'Position', position, ...
    'Color', 'none', ...
    'XAxisLocation', 'bottom', ...
    'YAxisLocation', 'right', ...
    'YTick', [], ...
    'XScale', 'log', ...
    'XLim', base_axis.XLim, ...
    'XTick', x_ticks, ...
    'XTickLabel', snr_tick_labels(snr_values), ...
    'Box', 'off');
xlabel(ax, label)
end


function labels = snr_tick_labels(snr_values)
%SNR_TICK_LABELS  Format SNR values for the two upper axes.

labels = strings(size(snr_values));

for ivalue = 1:numel(snr_values)
    if isinf(snr_values(ivalue))
        labels(ivalue) = "Inf";
    else
        labels(ivalue) = sprintf('%.1f', snr_values(ivalue));
    end
end
end
