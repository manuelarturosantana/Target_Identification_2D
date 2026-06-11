function results = classify_glrt_data(data, model, params)
%CLASSIFY_GLRT_DATA  Classify all receiver signals using a GLRT model.
%
%   results = CLASSIFY_GLRT_DATA(data, model, params) classifies every
%   receiver signal in data against every GLRT basis in model.
%s
%   The output intentionally stores only compact classification information:
%   true labels, predicted labels, receiver indices, GLRT scores, object
%   names/file paths, and representative SNR values.

    rng(params.random_seed);

    nobj = numel(data);
    sigmas = params.noise_sigmas(:);
    nsigma = numel(sigmas);

    late_idx = model(1).late_idx;

    % Count total trials.
    ntrials_per_sigma = 0;

    for j = 1:nobj
        ntrials_per_sigma = ntrials_per_sigma + size(data(j).signals, 1);
    end

    ntrials = nsigma * ntrials_per_sigma;

    sigma_col = zeros(ntrials, 1);
    true_idx = zeros(ntrials, 1);
    pred_idx = zeros(ntrials, 1);
    receiver_idx = zeros(ntrials, 1);
    score_max = zeros(ntrials, 1);
    correct = false(ntrials, 1);
    scores_all = zeros(ntrials, nobj);

    row = 0;
    start_time = tic;
    last_report_time = tic;

    if params.show_progress
        fprintf('Classifying %d signals...\n', ntrials);
    end

    for isig = 1:nsigma
        sigma = sigmas(isig);

        for itrue = 1:nobj
            signals = data(itrue).signals;
            nrec = size(signals, 1);

            for irec = 1:nrec
                row = row + 1;

                f = signals(irec, late_idx:end).';

                if sigma > 0
                    noise = (sigma / sqrt(2)) * ...
                        (randn(size(f)) + 1i * randn(size(f)));
                    f = f + noise;
                end

                scores = zeros(1, nobj);

                for ipred = 1:nobj
                    Ur = model(ipred).Ur;

                    if size(Ur, 1) ~= numel(f)
                        error('classify_glrt_data:DimensionMismatch', ...
                              ['Signal length does not match GLRT basis ' ...
                               'for model entry %d.'], ipred);
                    end

                    scores(ipred) = norm(Ur' * f, 2);
                end

                [score_max(row), pred_idx(row)] = max(scores);

                sigma_col(row) = sigma;
                true_idx(row) = itrue;
                receiver_idx(row) = irec;
                correct(row) = pred_idx(row) == itrue;
                scores_all(row, :) = scores;

                % Progress report.
                if params.show_progress && ...
                        (toc(last_report_time) > 0.5 || row == ntrials)

                    elapsed = toc(start_time);
                    frac = row / ntrials;

                    if frac > 0
                        eta = elapsed * (1 - frac) / frac;
                    else
                        eta = NaN;
                    end

                    fprintf( ...
                        ['\r  %6d / %6d  (%6.2f%%)  ' ...
                         'elapsed: %7.1fs  eta: %7.1fs'], ...
                        row, ntrials, 100 * frac, elapsed, eta);

                    last_report_time = tic;
                end
            end
        end
    end

    if params.show_progress
        fprintf('\n');
    end

    trials = table( ...
        sigma_col, ...
        true_idx, ...
        pred_idx, ...
        receiver_idx, ...
        score_max, ...
        correct, ...
        'VariableNames', { ...
            'sigma', ...
            'true_idx', ...
            'pred_idx', ...
            'receiver_idx', ...
            'score_max', ...
            'correct'});

    object_idx = (1:nobj).';
    object_name = strings(nobj, 1);
    file_name = strings(nobj, 1);
    file_path = strings(nobj, 1);

    for j = 1:nobj
        object_name(j) = data(j).name;
        file_name(j) = data(j).file_name;
        file_path(j) = data(j).file_path;
    end

    objects = table( ...
        object_idx, ...
        object_name, ...
        file_name, ...
        file_path, ...
        'VariableNames', { ...
            'object_idx', ...
            'object_name', ...
            'file_name', ...
            'file_path'});

    results = struct();
    results.trials = trials;
    results.objects = objects;
    results.scores = scores_all;
    if strcmp(params.snr_report_type, "representative_object")
        results.snr = compute_representative_snr(data, model, params);
    elseif strcmp(params.snr_report_type, "averaged")
        results.snr = compute_average_snr(data, model, params);
    else
        error(params.snr_report_type + " is not an implemented report type")
    end
    results.accuracy = mean(trials.correct);

end