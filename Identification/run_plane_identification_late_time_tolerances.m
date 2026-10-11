function comparison = run_plane_identification_late_time_tolerances()
%RUN_PLANE_IDENTIFICATION_LATE_TIME_TOLERANCES  Stress-test pulse cutoff.
%
%   Compute the exact time at which the Gaussian incident-field envelope
%   falls below each requested absolute tolerance, then run the same plane
%   GLRT at every resulting global cutoff. Each run builds one SVD basis per
%   library object and reuses it for all receiver signals.

cfg = struct();
cfg.input_folder = fullfile('local', 'plane_data_360');
cfg.output_path = fullfile('Identification', 'classification_results', ...
    'plane_late_time_tolerance_comparison.mat');
cfg.figure_path = fullfile('Identification', 'classification_results', ...
    'plane_late_time_tolerance_comparison.pdf');
cfg.noise_sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
cfg.svd_tol = 1e-12;
cfg.random_seed = 2026;
cfg.show_progress = true;

% These reproduce the incident Gaussian used for the plane data.
cfg.gaussian_center_time = 30.0;
cfg.gaussian_frequency_half_width = 5.0;
cfg.gaussian_frequency_tolerance = 1e-10;
cfg.turnoff_tolerances = [1e-16, 1e-12, 1e-8];

output_folder = fileparts(cfg.output_path);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% gauss_picker chooses s so exp(-dist^2/s) equals the frequency tolerance.
gaussian_parameter = -cfg.gaussian_frequency_half_width^2 / ...
    log(cfg.gaussian_frequency_tolerance);
peak_envelope = sqrt(gaussian_parameter) / (2 * sqrt(pi));
turnoff_delay = gaussian_turnoff_delay( ...
    gaussian_parameter, cfg.turnoff_tolerances, peak_envelope);
late_time = cfg.gaussian_center_time + turnoff_delay;
turnoff_table = table(cfg.turnoff_tolerances(:), turnoff_delay(:), ...
    late_time(:), 'VariableNames', ...
    {'tolerance', 'turnoff_delay', 'late_time'});
disp(turnoff_table);

base_params = build_id_params( ...
    'input_folder', cfg.input_folder, ...
    'output_path', cfg.output_path, ...
    'late_time', late_time(1), ...
    'noise_sigmas', cfg.noise_sigmas, ...
    'svd_tol', cfg.svd_tol, ...
    'random_seed', cfg.random_seed, ...
    'noise_on_full_time_grid', true, ...
    'snr_report_type', "averaged", ...
    'show_progress', cfg.show_progress);

fprintf('Loading plane pole and time-domain data...\n');
data = load_pole_data(base_params);
runs = struct([]);

for k = 1:numel(late_time)
    params = base_params;
    params.late_time = late_time(k);
    fprintf(['\nTolerance %.1e: turnoff delay %.12g, ', ...
        'late time %.12g\n'], cfg.turnoff_tolerances(k), ...
        turnoff_delay(k), late_time(k));

    % Build once for this cutoff, then reuse all library SVDs during the run.
    model = build_glrt_model(data, params);
    results = classify_glrt_data(data, model, params);
    results.configuration = params;

    runs(k).tolerance = cfg.turnoff_tolerances(k);
    runs(k).turnoff_delay = turnoff_delay(k);
    runs(k).late_time = late_time(k);
    runs(k).results = results;
end

summary = build_tolerance_summary(runs);
disp(summary);

comparison = struct();
comparison.configuration = cfg;
comparison.gaussian_parameter = gaussian_parameter;
comparison.peak_envelope = peak_envelope;
comparison.turnoff_table = turnoff_table;
comparison.runs = runs;
comparison.summary = summary;

save(cfg.output_path, 'comparison', '-v7.3');
make_tolerance_figure(comparison, cfg.figure_path);
fprintf('Saved tolerance comparison to:\n  %s\n', cfg.output_path);
fprintf('Saved tolerance figure to:\n  %s\n', cfg.figure_path);
end


function delay = gaussian_turnoff_delay(parameter, tolerance, peak_envelope)
%GAUSSIAN_TURNOFF_DELAY  Solve |a(t)| = tolerance after the pulse center.
%
%   |a(t)| = peak_envelope * exp(-parameter * t^2 / 4).

if any(tolerance <= 0 | tolerance >= peak_envelope)
    error('gaussian_turnoff_delay:BadTolerance', ...
        'Each tolerance must lie strictly between zero and the peak envelope.');
end

delay = sqrt(-(4 / parameter) * log(tolerance / peak_envelope));
end


function summary = build_tolerance_summary(runs)
%BUILD_TOLERANCE_SUMMARY  Accuracy for every tolerance and noise level.

sigma = unique(runs(1).results.trials.sigma, 'sorted');
nrun = numel(runs);
accuracy = nan(numel(sigma), nrun);

for irun = 1:nrun
    trials = runs(irun).results.trials;

    for isig = 1:numel(sigma)
        mask = trials.sigma == sigma(isig);
        accuracy(isig, irun) = mean(trials.correct(mask));
    end
end

summary = table(sigma, accuracy);
summary.Properties.UserData.tolerances = [runs.tolerance];
summary.Properties.UserData.late_times = [runs.late_time];
end


function make_tolerance_figure(comparison, figure_path)
%MAKE_TOLERANCE_FIGURE  Plot classification sensitivity to pulse tolerance.

summary = comparison.summary;
runs = comparison.runs;
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 450]);
ax = axes(fig);
x = 1:height(summary);
hold(ax, 'on');

for k = 1:numel(runs)
    label = sprintf('tol = %.0e, t_L = %.3f', ...
        runs(k).tolerance, runs(k).late_time);
    plot(ax, x, 100 * summary.accuracy(:, k), '-o', ...
        'LineWidth', 1.5, 'DisplayName', label);
end

hold(ax, 'off');
grid(ax, 'on');
xticks(ax, x);
xticklabels(ax, compose('%g', summary.sigma));
xlabel(ax, 'Noise standard deviation');
ylabel(ax, 'Classification accuracy [%]');
ylim(ax, [0, 100]);
title(ax, 'Sensitivity to Gaussian turn-off tolerance');
legend(ax, 'Location', 'best');
exportgraphics(fig, figure_path, 'ContentType', 'vector');
close(fig);
end
