function selective = apply_glrt_reject_option(results, model, options)
%APPLY_GLRT_REJECT_OPTION  Abstain when no library subspace is detected.
%
%   selective = APPLY_GLRT_REJECT_OPTION(results, model, options) applies a
%   matched-subspace detection gate to an existing GLRT classification run.
%   It does not repeat any signal projections or SVDs.
%
%   Under the noise-only null and the real-noise convention used by
%   classify_glrt_data,
%
%       norm(U' * noise)^2 / (sigma^2 / 2) ~ chi2(rank(U)).
%
%   A prediction is accepted only when its winning subspace passes this
%   test after a Bonferroni correction across the library. Otherwise the
%   decision is IDK. The correction controls the family-wise false trapping
%   declaration probability at options.familywise_alpha.
%
%   References:
%     Scharf & Friedlander, "Matched subspace detectors," IEEE TSP, 1994.
%     DOI: 10.1109/78.301849
%     Chow, "On optimum recognition error and reject tradeoff," IEEE TIT,
%     1970. DOI: 10.1109/TIT.1970.1054406

    arguments
        results (1,1) struct
        model (1,:) struct
        options (1,1) struct
    end

    if ~isfield(options, 'familywise_alpha')
        error('apply_glrt_reject_option:MissingAlpha', ...
            'options.familywise_alpha is required.');
    end

    if ~isfield(options, 'alpha_grid')
        options.alpha_grid = logspace(-6, 0, 25).';
    end

    validateattributes(options.familywise_alpha, {'numeric'}, ...
        {'scalar', '>', 0, '<=', 1});

    trials = results.trials;
    ntrial = height(trials);
    nlib = numel(model);

    if size(results.scores, 1) ~= ntrial || ...
            size(results.scores, 2) ~= nlib
        error('apply_glrt_reject_option:ScoreSizeMismatch', ...
            'results.scores must have one row per trial and one column per model.');
    end

    ranks = [model.retained_rank];
    winning_linear_idx = sub2ind(size(results.scores), ...
        (1:ntrial).', trials.pred_idx);
    winning_score = results.scores(winning_linear_idx);
    winning_rank = ranks(trials.pred_idx).';
    sigma = trials.sigma;

    detection_statistic = nan(ntrial, 1);
    detection_p_value = ones(ntrial, 1);
    noisy = sigma > 0 & winning_rank > 0;
    detection_statistic(noisy) = winning_score(noisy).^2 ./ ...
        (sigma(noisy).^2 / 2);

    % gammainc(x/2, nu/2, 'upper') is the chi-square upper-tail
    % probability and does not require Statistics and Machine Learning Toolbox.
    detection_p_value(noisy) = gammainc( ...
        detection_statistic(noisy) / 2, winning_rank(noisy) / 2, 'upper');

    % With no injected noise, accept every nonzero subspace response. The
    % statistical threshold is otherwise undefined because its null variance
    % is zero.
    noiseless = sigma == 0 & winning_rank > 0;
    nonzero_floor = sqrt(eps(class(winning_score)));
    detection_statistic(noiseless & winning_score > nonzero_floor) = Inf;
    detection_p_value(noiseless & winning_score > nonzero_floor) = 0;

    adjusted_p_value = min(1, nlib * detection_p_value);
    knows = adjusted_p_value <= options.familywise_alpha;
    idk = ~knows;

    trials.winning_rank = winning_rank;
    trials.detection_statistic = detection_statistic;
    trials.detection_p_value = detection_p_value;
    trials.adjusted_p_value = adjusted_p_value;
    trials.knows = knows;
    trials.idk = idk;

    selective = struct();
    selective.trials = trials;
    selective.familywise_alpha = options.familywise_alpha;
    selective.library_size = nlib;
    selective.idk_rate = mean(idk);
    selective.coverage = mean(knows);

    if any(knows)
        selective.accuracy_when_known = mean(trials.correct(knows));
    else
        selective.accuracy_when_known = NaN;
    end

    selective.summary = summarize_by_sigma(trials);
    selective.risk_coverage = build_risk_coverage( ...
        trials, adjusted_p_value, options.alpha_grid);
    selective.reference_dois = [ ...
        "10.1109/78.301849"; ...
        "10.1109/TIT.1970.1054406"];
end


function summary = summarize_by_sigma(trials)
%SUMMARIZE_BY_SIGMA  Coverage and conditional accuracy at each noise level.

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


function curve = build_risk_coverage(trials, adjusted_p_value, alpha_grid)
%BUILD_RISK_COVERAGE  Evaluate the reject tradeoff without rerunning GLRT.

alpha = unique(alpha_grid(:), 'sorted');

if any(alpha <= 0 | alpha > 1)
    error('apply_glrt_reject_option:BadAlphaGrid', ...
        'Every alpha_grid entry must lie in (0, 1].');
end

coverage = zeros(size(alpha));
accuracy_when_known = nan(size(alpha));
risk_when_known = nan(size(alpha));

for k = 1:numel(alpha)
    knows = adjusted_p_value <= alpha(k);
    coverage(k) = mean(knows);

    if any(knows)
        accuracy_when_known(k) = mean(trials.correct(knows));
        risk_when_known(k) = 1 - accuracy_when_known(k);
    end
end

curve = table(alpha, coverage, accuracy_when_known, risk_when_known);
end
