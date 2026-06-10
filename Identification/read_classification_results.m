% analyze_classification_results.m
%
% Analyze GLRT target identification results.

clear
clc

% ---- user-defined parameters ----
results_path = fullfile( ...
    'Identification', ...
    'classification_results', ...
    'constellation_results.mat');

% Use 'late' or 'full' for the SNR axis.
snr_axis = 'late';

% Pole plots for misclassifications.
plot_misclassified_poles = true;
max_misclassification_plots = 10;

% Leave empty to inspect all sigmas.
% Example: sigma_to_inspect = 1e-3;
sigma_to_inspect = [];

% load classification results
S = load(results_path, 'results');
results = S.results;

trials = results.trials;
objects = results.objects;
snr = results.snr;

% compute accuracy by sigma
sigmas = unique(trials.sigma);
nsigma = numel(sigmas);

accuracy = zeros(nsigma, 1);
ntrials = zeros(nsigma, 1);
nwrong = zeros(nsigma, 1);

for k = 1:nsigma
    idx = trials.sigma == sigmas(k);

    ntrials(k) = nnz(idx);
    nwrong(k) = nnz(~trials.correct(idx));
    accuracy(k) = mean(trials.correct(idx));
end

summary = table( ...
    sigmas, ...
    accuracy, ...
    ntrials, ...
    nwrong, ...
    'VariableNames', { ...
        'sigma', ...
        'accuracy', ...
        'ntrials', ...
        'nwrong'});

disp(summary)

% attach representative SNR values
snr_full_db = nan(nsigma, 1);
snr_late_db = nan(nsigma, 1);

for k = 1:nsigma
    idx = snr.sigma == sigmas(k);

    if any(idx)
        snr_full_db(k) = snr.snr_full_db(find(idx, 1));
        snr_late_db(k) = snr.snr_late_db(find(idx, 1));
    end
end

% plot accuracy with sigma bottom axis and SNR top axis
positive_sigmas = sigmas(sigmas > 0);

if isempty(positive_sigmas)
    warning('No positive sigma values available for log-scale plotting.');
else
    min_positive_sigma = min(positive_sigmas);
    sigma_zero_plot = min_positive_sigma / 10;

    sigmas_plot = sigmas;
    sigmas_plot(sigmas_plot == 0) = sigma_zero_plot;

    [sigmas_plot, order] = sort(sigmas_plot);
    accuracy_plot = accuracy(order);

    if strcmpi(snr_axis, 'late')
        snr_plot = snr_late_db(order);
        snr_label = 'Late-time SNR [dB]';
    elseif strcmpi(snr_axis, 'full')
        snr_plot = snr_full_db(order);
        snr_label = 'Full-time SNR [dB]';
    else
        error('Unknown snr_axis "%s". Use "late" or "full".', snr_axis);
    end

    sigma_labels = strings(size(sigmas_plot));

    for k = 1:numel(sigmas_plot)
        if sigmas_plot(k) == sigma_zero_plot
            sigma_labels(k) = "0";
        else
            sigma_labels(k) = sprintf('%.0e', sigmas_plot(k));
        end
    end

    snr_labels = strings(size(snr_plot));

    for k = 1:numel(snr_plot)
        if isinf(snr_plot(k))
            snr_labels(k) = "Inf";
        else
            snr_labels(k) = sprintf('%.1f', snr_plot(k));
        end
    end

    figure

    ax1 = axes;
    semilogx(ax1, sigmas_plot, accuracy_plot, '-o', 'LineWidth', 1.5)
    grid(ax1, 'on')
    xlabel(ax1, '\sigma')
    ylabel(ax1, 'Classification accuracy')
    title(ax1, 'GLRT classification accuracy')
    ylim(ax1, [0, 1])

    ax1.XTick = sigmas_plot;
    ax1.XTickLabel = sigma_labels;

    ax2 = axes( ...
        'Position', ax1.Position, ...
        'Color', 'none', ...
        'XAxisLocation', 'top', ...
        'YAxisLocation', 'right', ...
        'YTick', [], ...
        'XScale', 'log', ...
        'XLim', ax1.XLim);

    ax2.XTick = sigmas_plot;
    ax2.XTickLabel = snr_labels;
    xlabel(ax2, snr_label)

    linkaxes([ax1, ax2], 'x')
end

% optionally plot poles of misclassified true/predicted pairs
if plot_misclassified_poles
    wrong = trials(~trials.correct, :);

    if ~isempty(sigma_to_inspect)
        wrong = wrong(wrong.sigma == sigma_to_inspect, :);
    end

    if isempty(wrong)
        fprintf('No misclassifications found.\n');
    else
        pairs = unique([wrong.true_idx, wrong.pred_idx], 'rows');
        nplots = min(size(pairs, 1), max_misclassification_plots);

        fprintf('Plotting %d misclassified true/predicted pole pairs.\n', ...
            nplots);

        for k = 1:nplots
            true_idx = pairs(k, 1);
            pred_idx = pairs(k, 2);

            true_name = objects.object_name(true_idx);
            pred_name = objects.object_name(pred_idx);

            true_path = objects.file_path(true_idx);
            pred_path = objects.file_path(pred_idx);

            T = load(true_path, 'pols');
            P = load(pred_path, 'pols');

            figure
            plot(real(T.pols), imag(T.pols), 'o', 'DisplayName', 'true')
            hold on
            plot(real(P.pols), imag(P.pols), 'x', 'DisplayName', 'predicted')
            grid on
            xlabel('Re pole')
            ylabel('Im pole')
            legend('Location', 'best')

            title({ ...
                'Misclassified pole comparison', ...
                ['true: ', char(true_name)], ...
                ['predicted: ', char(pred_name)]});
        end
    end
end