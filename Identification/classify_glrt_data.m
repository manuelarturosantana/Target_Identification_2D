function results = classify_glrt_data(data, model, params)
%CLASSIFY_GLRT_DATA  Classify all receiver signals using a GLRT model.
%
%   results = CLASSIFY_GLRT_DATA(data, model, params) classifies every
%   receiver signal in data against every GLRT basis in model.
%
%   Default behavior:
%       If no classification map is attached to model, this function
%       requires numel(data) == numel(model) and behaves as before.
%
%   Flexible behavior:
%       If model(1).classification.true_to_library_map exists, then data
%       may have a different length from model. The map determines which
%       model/library predictions count as correct.

    rng(params.random_seed);

    neval = numel(data);
    nlib = numel(model);

    has_classification_map = ...
        isfield(model, 'classification') && ...
        ~isempty(model(1).classification) && ...
        isfield(model(1).classification, 'true_to_library_map');

    if has_classification_map
        true_to_library_map = logical(model(1).classification.true_to_library_map);

        if ~isequal(size(true_to_library_map), [neval, nlib])
            error('classify_glrt_data:BadTruthMapSize', ...
                  ['model(1).classification.true_to_library_map must have ', ...
                   'size numel(data)-by-numel(model), i.e. %d-by-%d.'], ...
                  neval, nlib);
        end

        if isfield(model(1).classification, 'library_data')
            library_data = model(1).classification.library_data;
        else
            error('classify_glrt_data:MissingLibraryData', ...
                  ['Flexible classification requires ', ...
                   'model(1).classification.library_data.']);
        end

        if numel(library_data) ~= nlib
            error('classify_glrt_data:LibraryDataMismatch', ...
                  'numel(model(1).classification.library_data) must equal numel(model).');
        end
    else
        if neval ~= nlib
            error('classify_glrt_data:LengthMismatch', ...
                  ['numel(data) and numel(model) differ. ', ...
                   'Attach model(1).classification.true_to_library_map first.']);
        end

        true_to_library_map = eye(neval, nlib) > 0;
        library_data = data;
    end

    sigmas = params.noise_sigmas(:);
    nsigma = numel(sigmas);

    late_idx = model(1).late_idx;

    % Count total trials.
    ntrials_per_sigma = 0;

    for j = 1:neval
        ntrials_per_sigma = ntrials_per_sigma + size(data(j).signals, 1);
    end

    ntrials = nsigma * ntrials_per_sigma;

    sigma_col = zeros(ntrials, 1);
    true_idx = zeros(ntrials, 1);      % index into data
    pred_idx = zeros(ntrials, 1);      % index into model/library_data
    receiver_idx = zeros(ntrials, 1);
    score_max = zeros(ntrials, 1);
    correct = false(ntrials, 1);
    scores_all = zeros(ntrials, nlib);

    row = 0;
    start_time = tic;
    last_report_time = tic;

    if params.show_progress
        fprintf('Classifying %d signals...\n', ntrials);
    end

    for isig = 1:nsigma
        sigma = sigmas(isig);

        for itrue = 1:neval
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

                scores = zeros(1, nlib);

                for ipred = 1:nlib
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
                correct(row) = true_to_library_map(itrue, pred_idx(row));
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

    % Preserve the original trials table fields.
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

    % In the old path, this is exactly the old object table.
    % In the flexible path, this is the prediction library.
    object_idx = (1:nlib).';
    object_name = strings(nlib, 1);
    file_name = strings(nlib, 1);
    file_path = strings(nlib, 1);

    for j = 1:nlib
        object_name(j) = library_data(j).name;
        file_name(j) = library_data(j).file_name;
        file_path(j) = library_data(j).file_path;
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
        results.snr = compute_representative_snr(library_data, model, params);
    elseif strcmp(params.snr_report_type, "averaged")
        results.snr = compute_average_snr(library_data, model, params);
    else
        error(params.snr_report_type + " is not an implemented report type")
    end

    results.accuracy = mean(trials.correct);

    % Only add extra bookkeeping in flexible mode.
    % This avoids changing the old result structure.
    if has_classification_map
        eval_idx = (1:neval).';
        eval_object_name = strings(neval, 1);
        eval_file_name = strings(neval, 1);
        eval_file_path = strings(neval, 1);

        for j = 1:neval
            eval_object_name(j) = data(j).name;
            eval_file_name(j) = data(j).file_name;
            eval_file_path(j) = data(j).file_path;
        end

        eval_objects = table( ...
            eval_idx, ...
            eval_object_name, ...
            eval_file_name, ...
            eval_file_path, ...
            'VariableNames', { ...
                'eval_idx', ...
                'eval_object_name', ...
                'eval_file_name', ...
                'eval_file_path'});

        results.eval_objects = eval_objects;
        results.library_objects = objects;
        results.true_to_library_map = true_to_library_map;
    end

end