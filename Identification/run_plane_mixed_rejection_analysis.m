function analysis = run_plane_mixed_rejection_analysis()
%RUN_PLANE_MIXED_REJECTION_ANALYSIS  Combine angle and noise-significance IDK.
%
% Re-scores the original plane experiment stored in
% plane_results_averaged.mat with a two-part reject option:
%
%   1. angle gate:      commit only if |receiver angle| <= b, where 0 deg
%                       is nose-on and 180 deg is tail-on (monostatic).
%   2. noise gate:      commit only if the winning subspace energy is
%                       significant under the noise-only null,
%
%          T = ||U' f||^2 / (sigma^2/2) ~ chi2(rank(U)),
%
%                       after a Bonferroni correction over the library,
%                       at family-wise level alpha.
%
% b = 180 disables the angle gate and alpha = 1 disables the noise gate.
% Neither gate changes the GLRT argmax, so every rejected sample can be
% labelled by what the classifier *would* have answered.
%
% Rejected samples are split by cause:
%   angle_only  angle gate fails, noise gate passes
%   noise_only  angle gate passes, noise gate fails
%   both        both gates fail
% and by outcome: lost (would have been correct) or avoided (would have
% been wrong).
%
% Library ranks are not stored in the results file. They are estimated
% from the saved scores themselves, using that every sigma shares the same
% (object, receiver) row order:
%
%   E[score_sigma^2 - score_0^2] = rank * sigma^2 / 2.
%
% This matters because twinjet_plane1 was regenerated after the local copy
% of plane_data_360 was made; its library (poles) is not available locally.

cfg = struct();
cfg.results_path = fullfile('Identification', 'classification_results', ...
    'plane_results_averaged.mat');
cfg.output_folder = fullfile('Identification', 'classification_results', ...
    'plane_mixed_rejection');
cfg.output_path = fullfile(cfg.output_folder, ...
    'plane_mixed_rejection_analysis.mat');
cfg.csv_path = fullfile(cfg.output_folder, ...
    'plane_mixed_rejection_by_sigma.csv');
cfg.report_path = fullfile(cfg.output_folder, ...
    'plane_mixed_rejection_report.txt');
cfg.angle_bounds_degrees = [15, 30, 60, 90, 120, 150, 180];
cfg.familywise_alphas = [1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 2e-2, 5e-2, 1e-1, 1];
cfg.rank_estimation_sigma = 1e-1;
cfg.abs_angle_bin_edges_degrees = 0:15:180;

if ~exist(cfg.output_folder, 'dir')
    mkdir(cfg.output_folder);
end

loaded = load(cfg.results_path);
results = loaded.results;
trials = results.trials;
nlib = size(results.scores, 2);

[ranks, rank_estimate, rank_se] = estimate_library_ranks(results, ...
    cfg.rank_estimation_sigma);

angle_degrees = rad2deg(results.snr.receiver_angle(trials.receiver_idx));
angle_degrees = mod(angle_degrees(:) + 180, 360) - 180;
abs_angle = abs(angle_degrees);

adjusted_p = bonferroni_winner_p_values(results.scores, trials.pred_idx, ...
    trials.sigma, ranks, nlib);

sigmas = unique(trials.sigma, 'sorted');
[by_sigma, pooled, pooled_noisy] = tabulate_grid(trials, abs_angle, ...
    adjusted_p, sigmas, cfg);

by_angle_bin = tabulate_angle_bins(trials, abs_angle, adjusted_p, ...
    sigmas, cfg);

statistical_by_receiver = tabulate_statistical_by_receiver(trials, ...
    angle_degrees, adjusted_p, sigmas, cfg.familywise_alphas);

analysis = struct();
analysis.configuration = cfg;
analysis.library_ranks = table(results.objects.object_name, ...
    rank_estimate(:), rank_se(:), ranks(:), 'VariableNames', ...
    {'object_name', 'rank_estimate', 'rank_standard_error', 'rank'});
analysis.trial_angle_degrees = angle_degrees;
analysis.trial_abs_angle_degrees = abs_angle;
analysis.trial_adjusted_p_value = adjusted_p;
analysis.by_sigma = by_sigma;
analysis.pooled = pooled;
analysis.pooled_noisy = pooled_noisy;
analysis.by_angle_bin = by_angle_bin;
analysis.statistical_by_receiver = statistical_by_receiver;

save(cfg.output_path, 'analysis');
writetable(by_sigma, cfg.csv_path);
write_report(cfg.report_path, analysis);

fprintf('Saved analysis to:\n  %s\n  %s\n  %s\n', cfg.output_path, ...
    cfg.csv_path, cfg.report_path);
end


function [ranks, estimate, standard_error] = estimate_library_ranks( ...
        results, sigma)
%ESTIMATE_LIBRARY_RANKS  Recover rank(U) from paired noisy/clean scores.

trials = results.trials;
clean = results.scores(trials.sigma == 0, :);
noisy = results.scores(trials.sigma == sigma, :);

if isempty(clean) || ~isequal(size(clean), size(noisy)) || ...
        ~isequal(trials.receiver_idx(trials.sigma == 0), ...
        trials.receiver_idx(trials.sigma == sigma))
    error('estimate_library_ranks:RowMismatch', ...
        'Clean and noisy rows must share the same object/receiver order.');
end

excess = (noisy.^2 - clean.^2) / (sigma^2 / 2);
estimate = mean(excess, 1);
standard_error = std(excess, 0, 1) / sqrt(size(excess, 1));
ranks = round(estimate);

if any(abs(estimate - ranks) > 4 * standard_error + 0.25)
    warning('estimate_library_ranks:NotInteger', ...
        'Some rank estimates are far from an integer: %s', ...
        mat2str(estimate, 4));
end
end


function adjusted_p = bonferroni_winner_p_values(scores, pred_idx, ...
        sigma, ranks, nlib)
%BONFERRONI_WINNER_P_VALUES  Noise-only p-value of the winning subspace.

ntrial = size(scores, 1);
winning_score = scores(sub2ind(size(scores), (1:ntrial).', pred_idx));
winning_rank = ranks(pred_idx);
winning_rank = winning_rank(:);
p = ones(ntrial, 1);

noisy = sigma > 0;
statistic = winning_score.^2 ./ (sigma.^2 / 2);

for rank_value = unique(winning_rank(noisy)).'
    mask = noisy & winning_rank == rank_value;
    p(mask) = gammainc(statistic(mask) / 2, rank_value / 2, 'upper');
end

% With no noise the null variance is zero, so any nonzero response is
% significant.
p(~noisy & winning_score > sqrt(eps)) = 0;
adjusted_p = min(1, nlib * p);
end


function [by_sigma, pooled, pooled_noisy] = tabulate_grid(trials, ...
        abs_angle, adjusted_p, sigmas, cfg)
%TABULATE_GRID  Count every accept/reject outcome on the (sigma,b,alpha) grid.

bounds = cfg.angle_bounds_degrees;
alphas = cfg.familywise_alphas;
rows = cell(numel(sigmas), numel(bounds), numel(alphas));
pooled_rows = cell(numel(bounds), numel(alphas));
noisy_rows = cell(numel(bounds), numel(alphas));

for ib = 1:numel(bounds)
    angle_ok = abs_angle <= bounds(ib) + 1e-9;

    for ia = 1:numel(alphas)
        noise_ok = adjusted_p <= alphas(ia);

        for is = 1:numel(sigmas)
            mask = trials.sigma == sigmas(is);
            rows{is, ib, ia} = count_outcomes(mask, angle_ok, noise_ok, ...
                trials.correct, sigmas(is), bounds(ib), alphas(ia));
        end

        pooled_rows{ib, ia} = count_outcomes(true(height(trials), 1), ...
            angle_ok, noise_ok, trials.correct, NaN, bounds(ib), ...
            alphas(ia));
        noisy_rows{ib, ia} = count_outcomes(trials.sigma > 0, ...
            angle_ok, noise_ok, trials.correct, NaN, bounds(ib), ...
            alphas(ia));
    end
end

by_sigma = sortrows(vertcat(rows{:}), {'sigma', 'angle_bound_deg', ...
    'familywise_alpha'});
pooled = sortrows(vertcat(pooled_rows{:}), {'angle_bound_deg', ...
    'familywise_alpha'});
pooled_noisy = sortrows(vertcat(noisy_rows{:}), {'angle_bound_deg', ...
    'familywise_alpha'});
end


function row = count_outcomes(mask, angle_ok, noise_ok, correct, sigma, ...
        bound, alpha)
%COUNT_OUTCOMES  One summary row for a subset of trials and one rule.

committed = mask & angle_ok & noise_ok;
angle_only = mask & ~angle_ok & noise_ok;
noise_only = mask & angle_ok & ~noise_ok;
both = mask & ~angle_ok & ~noise_ok;

ntrials = nnz(mask);
ncommitted = nnz(committed);
nrejected = ntrials - ncommitted;
ncorrect_committed = nnz(committed & correct);

row = table(sigma, bound, alpha, ntrials, ncommitted, nrejected, ...
    ncorrect_committed, ncommitted - ncorrect_committed, ...
    nnz(mask & correct), ...
    nnz(angle_only), nnz(noise_only), nnz(both), ...
    nnz(angle_only & correct), nnz(noise_only & correct), ...
    nnz(both & correct), ...
    nnz(angle_only & ~correct), nnz(noise_only & ~correct), ...
    nnz(both & ~correct), ...
    'VariableNames', {'sigma', 'angle_bound_deg', 'familywise_alpha', ...
    'ntrials', 'ncommitted', 'nrejected', 'ncorrect_committed', ...
    'nwrong_committed', 'ncorrect_without_rejection', ...
    'nrej_angle_only', 'nrej_noise_only', 'nrej_both', ...
    'lost_angle_only', 'lost_noise_only', 'lost_both', ...
    'avoided_angle_only', 'avoided_noise_only', 'avoided_both'});

row.coverage = ncommitted / ntrials;
row.accuracy_when_committed = ncorrect_committed / max(ncommitted, 1);
row.accuracy_when_committed(ncommitted == 0) = NaN;
row.overall_correct_rate = ncorrect_committed / ntrials;
row.nlost = row.lost_angle_only + row.lost_noise_only + row.lost_both;
row.navoided = row.avoided_angle_only + row.avoided_noise_only + ...
    row.avoided_both;
% Fraction of rejections that removed an answer that would have been wrong.
row.rejection_precision = row.navoided / max(nrejected, 1);
row.rejection_precision(nrejected == 0) = NaN;
end


function summary = tabulate_angle_bins(trials, abs_angle, adjusted_p, ...
        sigmas, cfg)
%TABULATE_ANGLE_BINS  Accuracy and noise-gate pass rate per viewing sector.

edges = cfg.abs_angle_bin_edges_degrees;
bin = discretize(abs_angle, edges, 'IncludedEdge', 'right');
bin(abs_angle == 0) = 1;
alphas = cfg.familywise_alphas;
nbin = numel(edges) - 1;
parts = cell(numel(sigmas), 1);

for is = 1:numel(sigmas)
    sigma = repmat(sigmas(is), nbin, 1);
    lower_deg = edges(1:end-1).';
    upper_deg = edges(2:end).';
    ntrials = zeros(nbin, 1);
    accuracy = zeros(nbin, 1);
    pass_rate = zeros(nbin, numel(alphas));

    for k = 1:nbin
        mask = trials.sigma == sigmas(is) & bin == k;
        ntrials(k) = nnz(mask);
        accuracy(k) = mean(trials.correct(mask));

        for ia = 1:numel(alphas)
            pass_rate(k, ia) = mean(adjusted_p(mask) <= alphas(ia));
        end
    end

    parts{is} = table(sigma, lower_deg, upper_deg, ntrials, accuracy, ...
        pass_rate);
end

summary = vertcat(parts{:});
end


function summary = tabulate_statistical_by_receiver(trials, ...
        angle_degrees, adjusted_p, sigmas, alphas)
%TABULATE_STATISTICAL_BY_RECEIVER  Noise-gate outcomes per receiver angle.
%
% One row per (sigma, alpha, receiver), counting over objects the samples
% the noise gate alone keeps or rejects, split by whether the GLRT answer
% is correct. Angles are signed, in (-180, 180], 0 deg = nose-on.

receivers = unique(trials.receiver_idx, 'sorted');
nrec = numel(receivers);
[~, receiver_row] = ismember(trials.receiver_idx, receivers);
receiver_angle = accumarray(receiver_row, angle_degrees, [nrec, 1], ...
    @(a) a(1));
parts = cell(numel(sigmas), numel(alphas));

for is = 1:numel(sigmas)
    in_sigma = trials.sigma == sigmas(is);
    count = @(mask) accumarray(receiver_row(mask), 1, [nrec, 1]);

    for ia = 1:numel(alphas)
        kept = adjusted_p <= alphas(ia);
        parts{is, ia} = table( ...
            repmat(sigmas(is), nrec, 1), repmat(alphas(ia), nrec, 1), ...
            receivers, receiver_angle, count(in_sigma), ...
            count(in_sigma & kept & trials.correct), ...
            count(in_sigma & kept & ~trials.correct), ...
            count(in_sigma & ~kept & trials.correct), ...
            count(in_sigma & ~kept & ~trials.correct), ...
            'VariableNames', {'sigma', 'familywise_alpha', ...
            'receiver_idx', 'angle_deg', 'ntrials', 'kept_correct', ...
            'kept_wrong', 'rejected_correct', 'rejected_wrong'});
    end
end

summary = vertcat(parts{:});
end


function write_report(report_path, analysis)
%WRITE_REPORT  Human-readable account of the mixed reject option.

cfg = analysis.configuration;
file_id = fopen(report_path, 'w');

if file_id < 0
    error('write_report:OpenFailed', 'Could not open %s.', report_path);
end

cleanup = onCleanup(@() fclose(file_id));

fprintf(file_id, 'Mixed angle + noise-significance reject option\n');
fprintf(file_id, 'Source: %s\n\n', cfg.results_path);
fprintf(file_id, 'Estimated library ranks (sigma = %g):\n', ...
    cfg.rank_estimation_sigma);
ranks = analysis.library_ranks;

for k = 1:height(ranks)
    fprintf(file_id, '  %-16s %7.2f +- %4.2f -> %d\n', ...
        ranks.object_name(k), ranks.rank_estimate(k), ...
        ranks.rank_standard_error(k), ranks.rank(k));
end

write_block(file_id, 'Pooled over all sigma', analysis.pooled);
write_block(file_id, 'Pooled over sigma > 0', analysis.pooled_noisy);

sigmas = unique(analysis.by_sigma.sigma, 'sorted');

for is = 1:numel(sigmas)
    rows = analysis.by_sigma(analysis.by_sigma.sigma == sigmas(is), :);
    write_block(file_id, sprintf('sigma = %g', sigmas(is)), rows);
end

clear cleanup
end


function write_block(file_id, title, rows)
%WRITE_BLOCK  One grid of rules as a tab-separated table.

fprintf(file_id, '\n=== %s ===\n', title);
fprintf(file_id, ['bound\talpha\tcover%%\tacc|commit%%\tcorrect/all%%\t', ...
    'rejected\tlost(angle,noise,both)\tavoided(angle,noise,both)\t', ...
    'rej.precision%%\n']);

for k = 1:height(rows)
    r = rows(k, :);
    fprintf(file_id, ['%g\t%g\t%.1f\t%.1f\t%.1f\t%d\t%d (%d,%d,%d)\t', ...
        '%d (%d,%d,%d)\t%.1f\n'], r.angle_bound_deg, ...
        r.familywise_alpha, 100 * r.coverage, ...
        100 * r.accuracy_when_committed, 100 * r.overall_correct_rate, ...
        r.nrejected, r.nlost, r.lost_angle_only, r.lost_noise_only, ...
        r.lost_both, r.navoided, r.avoided_angle_only, ...
        r.avoided_noise_only, r.avoided_both, ...
        100 * r.rejection_precision);
end
end
