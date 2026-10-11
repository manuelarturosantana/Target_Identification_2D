function comparison = run_plane_identification_idk()
%RUN_PLANE_IDENTIFICATION_IDK  Compare ordinary and selective plane GLRT.
%
%   Run plane identification once with the fixed late time t = 40. Then
%   apply a matched-subspace significance gate to the saved projection
%   scores. The selective method says IDK when the winning pole subspace
%   does not contain statistically significant energy above noise. Three
%   family-wise alpha values are evaluated without rerunning the GLRT.

cfg = struct();
cfg.input_folder = fullfile('local', 'plane_data_360');
cfg.output_path = fullfile('Identification', 'classification_results', ...
    'plane_idk_comparison.mat');
cfg.figure_path = fullfile('Identification', 'classification_results', ...
    'plane_idk_comparison.pdf');
cfg.late_time = 40.0;
cfg.noise_sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
cfg.svd_tol = 1e-12;
cfg.random_seed = 2026;
cfg.show_progress = true;
% Span the conservative regime through the no-rejection endpoint. Small
% alpha values behave similarly for these data, so include larger values to
% make the IDK/accuracy tradeoff visible. Alpha = 0.01 is the original value.
cfg.reject_familywise_alphas = [1e-6, 0.01, 0.1, 0.25, 0.5, 0.9, 1.0];
cfg.reject_alpha_grid = unique([logspace(-6, 0, 25).'; ...
    cfg.reject_familywise_alphas(:)]);

output_folder = fileparts(cfg.output_path);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

params = build_id_params( ...
    'input_folder', cfg.input_folder, ...
    'output_path', cfg.output_path, ...
    'late_time', cfg.late_time, ...
    'noise_sigmas', cfg.noise_sigmas, ...
    'svd_tol', cfg.svd_tol, ...
    'random_seed', cfg.random_seed, ...
    'noise_on_full_time_grid', true, ...
    'snr_report_type', "averaged", ...
    'show_progress', cfg.show_progress);

fprintf('Loading plane pole and time-domain data...\n');
data = load_pole_data(params);

% One basis per library object, reused for every signal and noise level.
model = build_glrt_model(data, params);
fprintf('Running fixed-cutoff plane identification at t = %.6g...\n', ...
    cfg.late_time);
ordinary = classify_glrt_data(data, model, params);
ordinary.configuration = params;

reject_options = struct();
reject_options.alpha_grid = cfg.reject_alpha_grid;
reject_options.familywise_alpha = cfg.reject_familywise_alphas(1);
selective_template = apply_glrt_reject_option( ...
    ordinary, model, reject_options);
selective = repmat(selective_template, size(cfg.reject_familywise_alphas));

for k = 1:numel(cfg.reject_familywise_alphas)
    selective(k) = apply_alpha_threshold( ...
        selective_template, cfg.reject_familywise_alphas(k));
end

summary = build_idk_summary(ordinary, selective);
disp(summary);

for k = 1:numel(selective)
    fprintf(['Selective GLRT at family-wise alpha = %.3g:\n', ...
        '  IDK rate:             %.2f%%\n', ...
        '  coverage:             %.2f%%\n', ...
        '  accuracy when known:  %.2f%%\n'], ...
        selective(k).familywise_alpha, 100 * selective(k).idk_rate, ...
        100 * selective(k).coverage, ...
        100 * selective(k).accuracy_when_known);
end

comparison = struct();
comparison.configuration = cfg;
comparison.ordinary = ordinary;
comparison.selective = selective;
comparison.summary = summary;

save(cfg.output_path, 'comparison', '-v7.3');
make_idk_figure(comparison, cfg.figure_path);
fprintf('Saved IDK comparison to:\n  %s\n', cfg.output_path);
fprintf('Saved IDK figure to:\n  %s\n', cfg.figure_path);
end


function selective = apply_alpha_threshold(selective, familywise_alpha)
%APPLY_ALPHA_THRESHOLD  Reuse detection p-values at a new alpha threshold.

trials = selective.trials;
knows = trials.adjusted_p_value <= familywise_alpha;
trials.knows = knows;
trials.idk = ~knows;

selective.trials = trials;
selective.familywise_alpha = familywise_alpha;
selective.idk_rate = mean(trials.idk);
selective.coverage = mean(knows);

if any(knows)
    selective.accuracy_when_known = mean(trials.correct(knows));
else
    selective.accuracy_when_known = NaN;
end

selective.summary = summarize_selective_trials(trials);
end


function summary = summarize_selective_trials(trials)
%SUMMARIZE_SELECTIVE_TRIALS  Selective metrics at each noise level.

sigma = unique(trials.sigma, 'sorted');
ntrials = zeros(size(sigma));
nknown = zeros(size(sigma));
nidk = zeros(size(sigma));
coverage = zeros(size(sigma));
idk_rate = zeros(size(sigma));
accuracy_when_known = nan(size(sigma));

for k = 1:numel(sigma)
    mask = trials.sigma == sigma(k);
    known_mask = mask & trials.knows;
    ntrials(k) = nnz(mask);
    nknown(k) = nnz(known_mask);
    nidk(k) = ntrials(k) - nknown(k);
    coverage(k) = nknown(k) / ntrials(k);
    idk_rate(k) = nidk(k) / ntrials(k);

    if nknown(k) > 0
        accuracy_when_known(k) = mean(trials.correct(known_mask));
    end
end

summary = table(sigma, ntrials, nknown, nidk, coverage, idk_rate, ...
    accuracy_when_known);
end


function summary = build_idk_summary(ordinary, selective)
%BUILD_IDK_SUMMARY  Metrics for every noise level and alpha value.

sigma = unique(ordinary.trials.sigma, 'sorted');
ordinary_accuracy = nan(size(sigma));

for k = 1:numel(sigma)
    mask = ordinary.trials.sigma == sigma(k);
    ordinary_accuracy(k) = mean(ordinary.trials.correct(mask));
end

nalpha = numel(selective);
nsigma = numel(sigma);
familywise_alpha = zeros(nalpha * nsigma, 1);
summary = table();

for k = 1:nalpha
    selective_summary = selective(k).summary;

    if ~isequal(sigma, selective_summary.sigma)
        error('build_idk_summary:SigmaMismatch', ...
            'Ordinary and selective summaries use different noise levels.');
    end

    rows = (k - 1) * nsigma + (1:nsigma);
    familywise_alpha(rows) = selective(k).familywise_alpha;
    block = table(familywise_alpha(rows), sigma, ...
        selective_summary.ntrials, ordinary_accuracy, ...
        selective_summary.nknown, selective_summary.nidk, ...
        selective_summary.coverage, selective_summary.idk_rate, ...
        selective_summary.accuracy_when_known, ...
        'VariableNames', {'familywise_alpha', 'sigma', 'ntrials', ...
        'ordinary_accuracy', 'nknown', 'nidk', 'coverage', 'idk_rate', ...
        'accuracy_when_known'});
    summary = [summary; block]; %#ok<AGROW>
end
end


function make_idk_figure(comparison, figure_path)
%MAKE_IDK_FIGURE  Plot conditional accuracy and abstention by alpha.

summary = comparison.summary;
selective = comparison.selective;
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1100, 460]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
sigma = unique(summary.sigma, 'sorted');
x = 1:numel(sigma);
colors = lines(numel(selective));

ax_accuracy = nexttile(layout);
first_alpha = summary.familywise_alpha == selective(1).familywise_alpha;
plot(ax_accuracy, x, 100 * summary.ordinary_accuracy(first_alpha), '-o', ...
    'LineWidth', 1.5, 'DisplayName', 'Always decide');
hold(ax_accuracy, 'on');
for k = 1:numel(selective)
    mask = summary.familywise_alpha == selective(k).familywise_alpha;
    plot(ax_accuracy, x, 100 * summary.accuracy_when_known(mask), '-s', ...
        'Color', colors(k, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Known, \\alpha = %g', ...
        selective(k).familywise_alpha));
end
hold(ax_accuracy, 'off');
grid(ax_accuracy, 'on');
xticks(ax_accuracy, x);
xticklabels(ax_accuracy, compose('%g', sigma));
xlabel(ax_accuracy, 'Noise standard deviation');
ylabel(ax_accuracy, 'Accuracy [%]');
ylim(ax_accuracy, [0, 100]);
legend(ax_accuracy, 'Location', 'best');
title(ax_accuracy, ...
    'Overall accuracy vs. accuracy conditional on ``known''');

ax_idk = nexttile(layout);
hold(ax_idk, 'on');
for k = 1:numel(selective)
    mask = summary.familywise_alpha == selective(k).familywise_alpha;
    plot(ax_idk, x, 100 * summary.idk_rate(mask), '-o', ...
        'Color', colors(k, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\alpha = %g', ...
        selective(k).familywise_alpha));
end
hold(ax_idk, 'off');
grid(ax_idk, 'on');
xticks(ax_idk, x);
xticklabels(ax_idk, compose('%g', sigma));
xlabel(ax_idk, 'Noise standard deviation');
ylabel(ax_idk, 'IDK rate [%]');
ylim(ax_idk, [0, 100]);
legend(ax_idk, 'Location', 'best');
title(ax_idk, 'Family-wise rejection thresholds');

exportgraphics(fig, figure_path, 'ContentType', 'vector');
close(fig);
end
