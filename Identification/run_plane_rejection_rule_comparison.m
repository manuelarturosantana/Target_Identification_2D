function comparison = run_plane_rejection_rule_comparison()
%RUN_PLANE_REJECTION_RULE_COMPARISON  Compare angle and noise IDK rules.
%
% Four rules classify the same noisy observations:
%   1. no rejection;
%   2. a candidate-orientation gate;
%   3. a Bonferroni-corrected noise-significance gate; and
%   4. the angle gate followed by the noise-significance gate, with the
%      multiplicity correction based only on angle-eligible candidates.
%
% The classifier uses real(f) after injecting circular complex noise as
%
%        noise = (sigma/sqrt(2)) * (randn + 1i*randn).
%
% Consequently, real(noise) ~ N(0, sigma^2/2). For an orthonormal real
% pole-subspace basis U of rank r, under a noise-only observation,
%
%        ||U' * noise||^2 / (sigma^2/2) ~ chi2(r).
%
% The angle decision uses only receiver angle and candidate boresight. The
% true class is never used to decide eligibility or rejection.

cfg = struct();
cfg.input_folder = fullfile('local', 'plane_data_360');
cfg.output_folder = fullfile('Identification', 'classification_results', ...
    'plane_rejection_comparison');
cfg.output_path = fullfile(cfg.output_folder, ...
    'plane_rejection_rule_comparison.mat');
cfg.report_path = fullfile(cfg.output_folder, ...
    'plane_rejection_rule_comparison.txt');
cfg.figure_path = fullfile(cfg.output_folder, ...
    'accuracy_vs_noise_level.pdf');
cfg.late_time = 40.0;
cfg.noise_sigmas = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
cfg.svd_tol = 1e-12;
cfg.random_seed = 2026;
cfg.show_progress = true;
cfg.angle_cutoff_degrees = 90;
cfg.familywise_alpha = 0.01;
% All plane library geometries are registered to the same zero-degree
% boresight. A scalar is expanded to one value per candidate class.
cfg.candidate_boresight_degrees = 0;
cfg.snr_bin_edges_db = [-Inf, -100, -60, -30, -10, 10, 30, 50, Inf];

if ~exist(cfg.output_folder, 'dir')
    mkdir(cfg.output_folder);
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
model = build_glrt_model(data, params);

if numel(data) ~= numel(model)
    error('run_plane_rejection_rule_comparison:LibraryMismatch', ...
        'This experiment requires one library class per evaluation geometry.');
end

boresight = expand_candidate_boresights( ...
    cfg.candidate_boresight_degrees, numel(model));

fprintf('Classifying shared signal-plus-noise observations...\n');
ordinary = classify_glrt_data(data, model, params);
signal_variants = evaluate_signal_variants(ordinary, model, boresight, cfg);

summaries = build_signal_summaries(ordinary, signal_variants, cfg);
disagreement = build_disagreement_summaries(ordinary, signal_variants);

comparison = struct();
comparison.configuration = cfg;
comparison.noise_distribution = struct( ...
    'sample_variance', 'sigma^2 / 2', ...
    'statistic', 'norm(U''*f)^2 / (sigma^2/2)', ...
    'null_distribution', 'chi-square with retained subspace rank df', ...
    'multiplicity_correction', ...
    'Bonferroni over the number of eligible candidate classes tested');
comparison.objects = ordinary.objects;
comparison.signal_trials = ordinary.trials;
comparison.signal_scores = ordinary.scores;
comparison.signal_variants = signal_variants;
comparison.summaries = summaries;
comparison.disagreement = disagreement;

save(cfg.output_path, 'comparison', '-v7.3');
make_accuracy_and_rejection_figure(summaries.by_sigma, cfg.figure_path);
write_report(cfg.report_path, summaries.overall, summaries.by_sigma, ...
    disagreement.overall, cfg);

disp(summaries.overall);
fprintf('Saved comparison to:\n  %s\n', cfg.output_path);
fprintf('Saved text report to:\n  %s\n', cfg.report_path);
fprintf('Saved figure to:\n  %s\n', cfg.figure_path);
end


function boresight = expand_candidate_boresights(value_degrees, nclass)
%EXPAND_CANDIDATE_BORESIGHTS  Return one known boresight per library class.

validateattributes(value_degrees, {'numeric'}, {'vector', 'finite'});

if isscalar(value_degrees)
    value_degrees = repmat(value_degrees, 1, nclass);
elseif numel(value_degrees) ~= nclass
    error('expand_candidate_boresights:BadSize', ...
        'Provide either one boresight or one per candidate class.');
end

boresight = deg2rad(value_degrees(:).');
end


function variants = evaluate_signal_variants(results, model, boresight, cfg)
%EVALUATE_SIGNAL_VARIANTS  Apply all four decisions to shared score rows.

trials = results.trials;
scores = results.scores;
ntrial = height(trials);
nclass = size(scores, 2);
ranks = [model.retained_rank];
all_eligible = true(ntrial, nclass);
angle_eligible = candidate_angle_eligibility( ...
    trials.receiver_angle, boresight, cfg.angle_cutoff_degrees);

[plain_pred, plain_has_candidate] = select_winner(scores, all_eligible);
[angle_pred, angle_has_candidate] = select_winner(scores, angle_eligible);

[plain_p, plain_adjusted_p] = corrected_winner_p_values( ...
    scores, plain_pred, trials.sigma, ranks, sum(all_eligible, 2), ...
    plain_has_candidate);
[angle_p, angle_adjusted_p] = corrected_winner_p_values( ...
    scores, angle_pred, trials.sigma, ranks, sum(angle_eligible, 2), ...
    angle_has_candidate);

variants = repmat(empty_signal_variant(ntrial), 1, 4);
variants(1) = make_signal_variant("no_rejection", "No rejection", ...
    plain_pred, plain_has_candidate, trials.true_idx, nan(ntrial, 1), ...
    nan(ntrial, 1), sum(all_eligible, 2));
variants(2) = make_signal_variant("angle_only", "Angle only", ...
    angle_pred, angle_has_candidate, trials.true_idx, nan(ntrial, 1), ...
    nan(ntrial, 1), sum(angle_eligible, 2));
variants(3) = make_signal_variant("noise_only", "Noise significance only", ...
    plain_pred, plain_has_candidate & plain_adjusted_p <= cfg.familywise_alpha, ...
    trials.true_idx, plain_p, plain_adjusted_p, sum(all_eligible, 2));
variants(4) = make_signal_variant("combined", "Angle + noise significance", ...
    angle_pred, angle_has_candidate & angle_adjusted_p <= cfg.familywise_alpha, ...
    trials.true_idx, angle_p, angle_adjusted_p, sum(angle_eligible, 2));
end


function variant = empty_signal_variant(ntrial)
%EMPTY_SIGNAL_VARIANT  Preallocate a consistently shaped variant struct.

variant = struct( ...
    'name', "", ...
    'label', "", ...
    'pred_idx', zeros(ntrial, 1), ...
    'committed', false(ntrial, 1), ...
    'correct', false(ntrial, 1), ...
    'raw_p_value', nan(ntrial, 1), ...
    'adjusted_p_value', nan(ntrial, 1), ...
    'eligible_count', zeros(ntrial, 1));
end


function variant = make_signal_variant(name, label, pred_idx, committed, ...
        true_idx, raw_p_value, adjusted_p_value, eligible_count)
%MAKE_SIGNAL_VARIANT  Package one signal decision rule.

variant = struct();
variant.name = name;
variant.label = label;
variant.pred_idx = pred_idx;
variant.committed = committed;
variant.correct = committed & pred_idx == true_idx;
variant.raw_p_value = raw_p_value;
variant.adjusted_p_value = adjusted_p_value;
variant.eligible_count = eligible_count;
end


function eligible = candidate_angle_eligibility(receiver_angle, ...
        candidate_boresight, cutoff_degrees)
%CANDIDATE_ANGLE_ELIGIBILITY  Gate candidates using prediction-time angles.

cutoff = deg2rad(cutoff_degrees);
offset = wrap_angle(receiver_angle(:) - candidate_boresight);
eligible = abs(offset) <= cutoff + 100 * eps(pi);
end


function [pred_idx, has_candidate] = select_winner(scores, eligible)
%SELECT_WINNER  Select the highest score after candidate eligibility.

masked_scores = scores;
masked_scores(~eligible) = -Inf;
[~, pred_idx] = max(masked_scores, [], 2);
has_candidate = any(eligible, 2);
pred_idx(~has_candidate) = 0;
end


function [raw_p, adjusted_p] = corrected_winner_p_values( ...
        scores, pred_idx, sigma, ranks, eligible_count, has_candidate)
%CORRECTED_WINNER_P_VALUES  Test winning projected energy under noise-only.

ntrial = size(scores, 1);
raw_p = ones(ntrial, 1);
adjusted_p = ones(ntrial, 1);
active = find(has_candidate);

if isempty(active)
    return
end

linear_idx = sub2ind(size(scores), active, pred_idx(active));
winning_score = scores(linear_idx);
winning_rank = ranks(pred_idx(active)).';
active_sigma = sigma(active);
active_p = ones(size(active));

noisy = active_sigma > 0 & winning_rank > 0;
statistic = nan(size(active));
statistic(noisy) = winning_score(noisy).^2 ./ ...
    (active_sigma(noisy).^2 / 2);

% Group by rank so gammainc receives a scalar shape parameter. This is much
% faster than passing a different shape parameter for every observation.
unique_ranks = unique(winning_rank(noisy));
for irank = 1:numel(unique_ranks)
    rank_value = unique_ranks(irank);
    mask = noisy & winning_rank == rank_value;
    active_p(mask) = gammainc(statistic(mask) / 2, ...
        rank_value / 2, 'upper');
end

noiseless = active_sigma == 0 & winning_rank > 0;
nonzero_floor = sqrt(eps(class(winning_score)));
active_p(noiseless & winning_score > nonzero_floor) = 0;

raw_p(active) = active_p;
adjusted_p(active) = min(1, eligible_count(active) .* active_p);
end


function summaries = build_signal_summaries(results, variants, cfg)
%BUILD_SIGNAL_SUMMARIES  Aggregate selective metrics across all dimensions.

trials = results.trials;
angle_degrees = round(rad2deg(wrap_angle(trials.receiver_angle)), 10);
abs_angle_degrees = round(abs(angle_degrees), 10);
geometry = results.objects.object_name(trials.true_idx);

summaries = struct();
summaries.overall = metric_table(variants, ones(height(trials), 1));

[group, sigma] = findgroups(trials.sigma);
summaries.by_sigma = add_group_column(metric_table(variants, group), ...
    'sigma', sigma);

[group, geometry_name] = findgroups(geometry);
summaries.by_geometry = add_group_column(metric_table(variants, group), ...
    'geometry', geometry_name);

[group, receiver_angle_degrees] = findgroups(angle_degrees);
summaries.by_angle = add_group_column(metric_table(variants, group), ...
    'receiver_angle_degrees', receiver_angle_degrees);

[group, absolute_angle_degrees] = findgroups(abs_angle_degrees);
summaries.by_abs_angle = add_group_column(metric_table(variants, group), ...
    'absolute_angle_degrees', absolute_angle_degrees);

[group, sigma, geometry_name] = findgroups(trials.sigma, geometry);
table_by_group = metric_table(variants, group);
table_by_group = add_group_column(table_by_group, 'geometry', geometry_name);
summaries.by_sigma_geometry = add_group_column( ...
    table_by_group, 'sigma', sigma);

[group, sigma, receiver_angle_degrees] = findgroups( ...
    trials.sigma, angle_degrees);
table_by_group = metric_table(variants, group);
table_by_group = add_group_column(table_by_group, ...
    'receiver_angle_degrees', receiver_angle_degrees);
summaries.by_sigma_angle = add_group_column( ...
    table_by_group, 'sigma', sigma);

finite_snr = isfinite(trials.snr_late_db);
snr_group = discretize(trials.snr_late_db(finite_snr), ...
    cfg.snr_bin_edges_db);
snr_variants = subset_signal_variants(variants, finite_snr);
present_bins = unique(snr_group, 'sorted');
compact_group = zeros(size(snr_group));
snr_center_db = zeros(numel(present_bins), 1);
snr_lower_db = zeros(numel(present_bins), 1);
snr_upper_db = zeros(numel(present_bins), 1);

for k = 1:numel(present_bins)
    compact_group(snr_group == present_bins(k)) = k;
    values = trials.snr_late_db(finite_snr);
    snr_center_db(k) = median(values(snr_group == present_bins(k)));
    snr_lower_db(k) = cfg.snr_bin_edges_db(present_bins(k));
    snr_upper_db(k) = cfg.snr_bin_edges_db(present_bins(k) + 1);
end

table_by_group = metric_table(snr_variants, compact_group);
table_by_group = add_group_column(table_by_group, ...
    'snr_upper_db', snr_upper_db);
table_by_group = add_group_column(table_by_group, ...
    'snr_lower_db', snr_lower_db);
summaries.by_snr = add_group_column(table_by_group, ...
    'snr_center_db', snr_center_db);
end


function subset = subset_signal_variants(variants, mask)
%SUBSET_SIGNAL_VARIANTS  Apply a common row subset to every variant.

subset = variants;
fields = {'pred_idx', 'committed', 'correct', 'raw_p_value', ...
    'adjusted_p_value', 'eligible_count'};

for ivariant = 1:numel(variants)
    for ifield = 1:numel(fields)
        field = fields{ifield};
        subset(ivariant).(field) = variants(ivariant).(field)(mask);
    end
end
end


function summary = metric_table(variants, group)
%METRIC_TABLE  Compute coverage and accuracy metrics for numbered groups.

group = group(:);
ngroup = max(group);
nvariant = numel(variants);
variant = strings(ngroup * nvariant, 1);
label = strings(ngroup * nvariant, 1);
group_index = zeros(ngroup * nvariant, 1);
ntrials = zeros(ngroup * nvariant, 1);
ncommitted = zeros(ngroup * nvariant, 1);
ncorrect_committed = zeros(ngroup * nvariant, 1);
coverage = zeros(ngroup * nvariant, 1);
accuracy_when_committed = nan(ngroup * nvariant, 1);
overall_correct_rate = zeros(ngroup * nvariant, 1);
idk_rate = zeros(ngroup * nvariant, 1);

for ivariant = 1:nvariant
    rows = (ivariant - 1) * ngroup + (1:ngroup);
    committed = variants(ivariant).committed;
    correct_committed = variants(ivariant).correct & committed;
    group_index(rows) = (1:ngroup).';
    variant(rows) = variants(ivariant).name;
    label(rows) = variants(ivariant).label;
    ntrials(rows) = accumarray(group, 1, [ngroup, 1]);
    ncommitted(rows) = accumarray(group, double(committed), [ngroup, 1]);
    ncorrect_committed(rows) = accumarray( ...
        group, double(correct_committed), [ngroup, 1]);
    coverage(rows) = ncommitted(rows) ./ ntrials(rows);
    known = ncommitted(rows) > 0;
    accuracy_when_committed(rows(known)) = ...
        ncorrect_committed(rows(known)) ./ ncommitted(rows(known));
    overall_correct_rate(rows) = ncorrect_committed(rows) ./ ntrials(rows);
    idk_rate(rows) = 1 - coverage(rows);
end

summary = table(group_index, variant, label, ntrials, ncommitted, ...
    ncorrect_committed, coverage, accuracy_when_committed, ...
    overall_correct_rate, idk_rate);
end


function output = add_group_column(input, name, group_values)
%ADD_GROUP_COLUMN  Expand group labels across variant blocks.

nvariant = numel(unique(input.variant, 'stable'));
expanded = repmat(group_values(:), nvariant, 1);
input.(name) = expanded;
output = movevars(input, name, 'Before', 'variant');
end


function summaries = build_disagreement_summaries(results, variants)
%BUILD_DISAGREEMENT_SUMMARIES  Compare accept/reject outcomes of both gates.

angle_accepts = variants(2).committed;
noise_accepts = variants(3).committed;
category = strings(height(results.trials), 1);
category(angle_accepts & noise_accepts) = "both_accept";
category(angle_accepts & ~noise_accepts) = "angle_only_accepts";
category(~angle_accepts & noise_accepts) = "noise_only_accepts";
category(~angle_accepts & ~noise_accepts) = "neither_accepts";

names = ["both_accept"; "angle_only_accepts"; ...
    "noise_only_accepts"; "neither_accepts"];
count = zeros(numel(names), 1);
for k = 1:numel(names)
    count(k) = nnz(category == names(k));
end
rate = count / numel(category);

summaries = struct();
summaries.category_by_trial = category;
summaries.overall = table(names, count, rate, ...
    'VariableNames', {'category', 'count', 'rate'});

angle_degrees = round(abs(rad2deg(wrap_angle( ...
    results.trials.receiver_angle))), 10);
[group, absolute_angle_degrees] = findgroups(angle_degrees);
ngroup = numel(absolute_angle_degrees);
group_count = accumarray(group, 1, [ngroup, 1]);
both_accept = accumarray(group, double(category == names(1)), [ngroup, 1]);
angle_only_accepts = accumarray( ...
    group, double(category == names(2)), [ngroup, 1]);
noise_only_accepts = accumarray( ...
    group, double(category == names(3)), [ngroup, 1]);
neither_accepts = accumarray( ...
    group, double(category == names(4)), [ngroup, 1]);

summaries.by_abs_angle = table(absolute_angle_degrees, group_count, ...
    both_accept ./ group_count, angle_only_accepts ./ group_count, ...
    noise_only_accepts ./ group_count, neither_accepts ./ group_count, ...
    'VariableNames', {'absolute_angle_degrees', 'ntrials', ...
    'both_accept', 'angle_only_accepts', 'noise_only_accepts', ...
    'neither_accepts'});
end


function make_accuracy_and_rejection_figure(summary, figure_path)
%MAKE_ACCURACY_AND_REJECTION_FIGURE  Plot accuracy and IDK on shared sigma.

sigmas = unique(summary.sigma, 'sorted');
positive_sigmas = sigmas(sigmas > 0);

if isempty(positive_sigmas)
    zero_plot_sigma = 1e-7;
else
    zero_plot_sigma = min(positive_sigmas) / 10;
end

plot_sigmas = sigmas;
plot_sigmas(sigmas == 0) = zero_plot_sigma;
variant_names = unique(summary.variant, 'stable');
colors = lines(numel(variant_names));

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 850, 820]);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
ax_accuracy = nexttile(layout);
ax_idk = nexttile(layout);

plot_sigma_metric(ax_accuracy, summary, variant_names, colors, ...
    plot_sigmas, 'accuracy_when_committed');
plot_sigma_metric(ax_idk, summary, variant_names, colors, ...
    plot_sigmas, 'idk_rate');

ylabel(ax_accuracy, 'Accuracy among committed predictions');
ylabel(ax_idk, 'Rejected / IDK fraction');
xlabel(ax_idk, 'Noise level (\sigma)');
set(ax_accuracy, 'XTickLabel', []);
legend(ax_accuracy, 'Location', 'southwest');
linkaxes([ax_accuracy, ax_idk], 'x');

for ax = [ax_accuracy, ax_idk]
    ylim(ax, [0, 1]);
    grid(ax, 'on');
    ax.XScale = 'log';
    ax.XTick = plot_sigmas;
    ax.FontSize = 14;
end
ax_idk.XTickLabel = sigma_tick_labels(sigmas);

exportgraphics(fig, figure_path, 'ContentType', 'vector');
close(fig);
end


function plot_sigma_metric(ax, summary, variant_names, colors, ...
        plot_sigmas, metric_name)
%PLOT_SIGMA_METRIC  Draw one metric for all four methods with -o styling.

hold(ax, 'on');
for ivariant = 1:numel(variant_names)
    mask = summary.variant == variant_names(ivariant);
    variant_summary = sortrows(summary(mask, :), 'sigma');
    semilogx(ax, plot_sigmas, variant_summary.(metric_name), '-o', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'Color', colors(ivariant, :), ...
        'DisplayName', variant_summary.label(1));
end
hold(ax, 'off');
end


function labels = sigma_tick_labels(sigmas)
%SIGMA_TICK_LABELS  Format fixed logarithmic ticks, including sigma = 0.

labels = strings(size(sigmas));
for k = 1:numel(sigmas)
    if sigmas(k) == 0
        labels(k) = "0";
    else
        labels(k) = sprintf('%.0e', sigmas(k));
    end
end
end


function write_report(report_path, overall, by_sigma, disagreement, cfg)
%WRITE_REPORT  Save a compact, reproducible text summary.

file_id = fopen(report_path, 'w');
if file_id < 0
    error('write_report:OpenFailed', 'Could not open %s.', report_path);
end
cleanup = onCleanup(@() fclose(file_id));

fprintf(file_id, 'Plane late-time pole-subspace rejection comparison\n');
fprintf(file_id, 'Late-time cutoff: %.6g\n', cfg.late_time);
fprintf(file_id, 'Angle cutoff: %.6g degrees\n', ...
    cfg.angle_cutoff_degrees);
fprintf(file_id, 'Family-wise alpha: %.6g\n\n', cfg.familywise_alpha);
fprintf(file_id, ['Noise derivation: after circular complex noise with ', ...
    'scale sigma/sqrt(2) is reduced to its real part, each sample has ', ...
    'variance sigma^2/2. Thus norm(U''*noise)^2/(sigma^2/2) is ', ...
    'chi-square with rank(U) degrees of freedom. The winning-class ', ...
    'p-value is Bonferroni corrected by the number of eligible classes ', ...
    'tested.\n\n']);

fprintf(file_id, ['Variant\tCoverage(%%)\tAccuracy committed(%%)\t', ...
    'Overall correct(%%)\tIDK(%%)\n']);
for k = 1:height(overall)
    fprintf(file_id, '%s\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        overall.label(k), 100 * overall.coverage(k), ...
        100 * overall.accuracy_when_committed(k), ...
        100 * overall.overall_correct_rate(k), ...
        100 * overall.idk_rate(k));
end

fprintf(file_id, '\nMetrics by noise level\n');
fprintf(file_id, ['Sigma\tVariant\tCoverage(%%)\t', ...
    'Accuracy committed(%%)\tOverall correct(%%)\tIDK(%%)\n']);
for k = 1:height(by_sigma)
    fprintf(file_id, '%.6g\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        by_sigma.sigma(k), by_sigma.label(k), ...
        100 * by_sigma.coverage(k), ...
        100 * by_sigma.accuracy_when_committed(k), ...
        100 * by_sigma.overall_correct_rate(k), ...
        100 * by_sigma.idk_rate(k));
end

fprintf(file_id, '\nAngle/noise gate agreement on signal observations\n');
fprintf(file_id, 'Category\tRate(%%)\tCount\n');
for k = 1:height(disagreement)
    fprintf(file_id, '%s\t%.4f\t%d\n', disagreement.category(k), ...
        100 * disagreement.rate(k), disagreement.count(k));
end

clear cleanup
end


function value = wrap_angle(value)
%WRAP_ANGLE  Map angles to [-pi, pi).

value = mod(value + pi, 2 * pi) - pi;
end
