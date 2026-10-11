% Experiment figures for the plane reject options (not paper figures).
%
% Reads the tables written by run_plane_mixed_rejection_analysis.m (run
% automatically if missing) and compares three criteria:
%   angle        commit only if |receiver angle| <= b (0 deg = nose-on)
%   statistical  commit only if the winning subspace energy is significant
%                under the noise-only null at family-wise level alpha
%   combined     both gates
%
% Figures 301-303: accuracy among kept samples vs sigma, formatted like
% Paper_Figure_Generator/all_plane_results_with_snr_axes.m.
% Figure 304: kept correct / kept wrong / rejected but would have been
% correct / rejected and would have been wrong, per sigma, one panel per
% criterion.
% Figure 305: the same four outcomes under the statistical criterion vs
% signed receiver angle (0 deg = nose-on), one panel per sigma.
clear;
script_directory = fileparts(mfilename('fullpath'));
project_root = fileparts(script_directory);
results_file = fullfile(script_directory, 'classification_results', ...
    'plane_results_averaged.mat');
output_folder = fullfile(script_directory, 'classification_results', ...
    'plane_mixed_rejection');
analysis_file = fullfile(output_folder, ...
    'plane_mixed_rejection_analysis.mat');

% Statistical gate level used by the statistical and combined criteria.
familywise_alpha = 1e-2;
% Angle windows drawn in the angle and combined accuracy figures.
angle_bounds = [15, 30, 60, 90];
% Alphas drawn in the statistical accuracy figure.
statistical_alphas = [1e-3, 1e-2, 2e-2, 5e-2];
% Angle window used by the angle and combined outcome panels.
outcome_angle_bound = 90;
% Width of the receiver-angle sectors in the angle outcome figure (about
% 10 samples per degree per sigma).
angle_bin_width = 10;
save_figures = true;

if isfile(analysis_file)
    loaded = load(analysis_file);
end

% Rerun if missing or written before the per-receiver table was added.
if ~isfile(analysis_file) || ...
        ~isfield(loaded.analysis, 'statistical_by_receiver')
    current_folder = cd(project_root);
    restore_folder = onCleanup(@() cd(current_folder));
    run_plane_mixed_rejection_analysis();
    clear restore_folder
    loaded = load(analysis_file);
end

by_sigma = loaded.analysis.by_sigma;
statistical_by_receiver = loaded.analysis.statistical_by_receiver;
loaded = load(results_file);
results = loaded.results;
sigmas = unique(by_sigma.sigma, 'sorted');

%% Accuracy among kept samples for each criterion.
no_rejection = select_rule(by_sigma, 180, 1).accuracy_when_committed;
angle_labels = [compose('[-%d^\\circ, %d^\\circ]', ...
    angle_bounds.', angle_bounds.').', {'No rejection'}];
alpha_labels = [compose('\\alpha = %g', statistical_alphas.').', ...
    {'No rejection'}];

angle_accuracy = rule_accuracy(by_sigma, angle_bounds, 1);
statistical_accuracy = zeros(numel(sigmas), numel(statistical_alphas));

for ia = 1:numel(statistical_alphas)
    statistical_accuracy(:, ia) = select_rule(by_sigma, 180, ...
        statistical_alphas(ia)).accuracy_when_committed;
end

combined_accuracy = [rule_accuracy(by_sigma, angle_bounds, ...
    familywise_alpha), ...
    select_rule(by_sigma, 180, familywise_alpha).accuracy_when_committed];
combined_labels = [angle_labels(1:end-1), ...
    {'[-180^\circ, 180^\circ]'}, {'No rejection'}];

snr = representative_snr(results, sigmas);

figures = gobjects(1, 5);
figures(1) = plot_accuracy(301, 'Angle-based rejection', sigmas, ...
    [angle_accuracy, no_rejection], angle_labels, snr, ...
    'angle gate only');
figures(2) = plot_accuracy(302, 'Statistical rejection', sigmas, ...
    [statistical_accuracy, no_rejection], alpha_labels, snr, ...
    'statistical gate only');
figures(3) = plot_accuracy(303, 'Combined rejection', sigmas, ...
    [combined_accuracy, no_rejection], combined_labels, snr, ...
    sprintf('angle gate + statistical gate, \\alpha = %g', ...
    familywise_alpha));

%% Outcome of every sample under each criterion.
criteria = struct( ...
    'title', { ...
        sprintf('Angle: [-%d^\\circ, %d^\\circ]', ...
            outcome_angle_bound, outcome_angle_bound), ...
        sprintf('Statistical: \\alpha = %g', familywise_alpha), ...
        sprintf('Combined: [-%d^\\circ, %d^\\circ], \\alpha = %g', ...
            outcome_angle_bound, outcome_angle_bound, familywise_alpha)}, ...
    'rows', { ...
        select_rule(by_sigma, outcome_angle_bound, 1), ...
        select_rule(by_sigma, 180, familywise_alpha), ...
        select_rule(by_sigma, outcome_angle_bound, familywise_alpha)});
figures(4) = plot_outcomes(304, criteria, sigmas);

%% Statistical outcomes vs receiver angle.
figures(5) = plot_angle_outcomes(305, statistical_by_receiver, sigmas, ...
    familywise_alpha, angle_bin_width);

for k = 1:numel(criteria)
    fprintf('\n%s\n', strrep(strrep(criteria(k).title, '^\circ', ''), ...
        '\alpha', 'alpha'));
    disp(criteria(k).rows(:, {'sigma', 'ncorrect_committed', ...
        'nwrong_committed', 'nlost', 'navoided', 'coverage', ...
        'accuracy_when_committed'}));
end

if save_figures
    names = {'accuracy_angle', 'accuracy_statistical', ...
        'accuracy_combined', 'outcomes', 'angle_outcomes_statistical'};
    tag = sprintf('_alpha_%g_bound_%d', familywise_alpha, ...
        outcome_angle_bound);

    for k = 1:numel(figures)
        exportgraphics(figures(k), fullfile(output_folder, ...
            [names{k}, tag, '.pdf']), 'ContentType', 'vector');
    end
end


function rows = select_rule(by_sigma, bound, alpha)
%SELECT_RULE  Rows of one (angle bound, alpha) rule, sorted by sigma.

mask = by_sigma.angle_bound_deg == bound & ...
    abs(by_sigma.familywise_alpha - alpha) < 1e-9 * alpha;
rows = sortrows(by_sigma(mask, :), 'sigma');

if isempty(rows)
    error('select_rule:Missing', ...
        'No rows for bound %g, alpha %g. Add it to the analysis grid.', ...
        bound, alpha);
end
end


function accuracy = rule_accuracy(by_sigma, bounds, alpha)
%RULE_ACCURACY  Accuracy among kept samples, one column per angle bound.

accuracy = [];

for ib = 1:numel(bounds)
    rows = select_rule(by_sigma, bounds(ib), alpha);
    accuracy = [accuracy, rows.accuracy_when_committed]; %#ok<AGROW>
end
end


function snr = representative_snr(results, sigmas)
%REPRESENTATIVE_SNR  Mean SNR over planes at the receiver nearest 0 deg.

receiver_distance = abs(angle(exp(1i * results.snr.receiver_angle(:))));
[~, zero_receiver_idx] = min(receiver_distance);
[~, rows] = ismember(sigmas, results.snr.sigma(:));
snr.full = results.snr.snr_full_mean_db(rows, zero_receiver_idx);
snr.late = results.snr.snr_late_mean_db(rows, zero_receiver_idx);
end


function fig = plot_accuracy(number, name, sigmas, accuracy, labels, ...
        snr, note)
%PLOT_ACCURACY  Accuracy vs sigma in the all_plane_results style.

plot_styles = {'-o', '-s', '-^', '-d', '-v', '-p'};
plot_sigmas = display_sigmas(sigmas);

fig = figure(number); clf;
fig.Name = name;
fig.Color = 'w';
ax_sigma = axes(fig);
hold(ax_sigma, 'on');

for k = 1:size(accuracy, 2) - 1
    semilogx(ax_sigma, plot_sigmas, accuracy(:, k), plot_styles{k}, ...
        'LineWidth', 2.5);
end

% The last column is always the no-rejection reference.
semilogx(ax_sigma, plot_sigmas, accuracy(:, end), '--k*', ...
    'LineWidth', 2.5);
hold(ax_sigma, 'off');

ax_sigma.XScale = 'log';
grid(ax_sigma, 'on');
xlabel(ax_sigma, "Noise Level (\sigma)", "FontSize", 14);
ylabel(ax_sigma, "Correct Classification Rate (kept samples)", ...
    "FontSize", 14);
legend(ax_sigma, labels, 'Location', 'SouthWest', 'FontSize', 12);
ax_sigma.FontSize = 14;
ax_sigma.YLim = [0, 1.02];
ax_sigma.XTick = plot_sigmas;
ax_sigma.XTickLabel = sigma_tick_labels(sigmas);
x_padding_factor = 1.15;
ax_sigma.XLim = [min(plot_sigmas) / x_padding_factor, ...
    max(plot_sigmas) * x_padding_factor];
text(ax_sigma, 0.98, 0.04, note, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 12);

main_position = ax_sigma.Position;
main_position(4) = main_position(4) - 0.22;
ax_sigma.Position = main_position;

ax_full = top_snr_axis(ax_sigma, plot_sigmas, snr.full, ...
    'Representative full-time SNR at 0^\circ [dB]', 0.03);
ax_late = top_snr_axis(ax_sigma, plot_sigmas, snr.late, ...
    'Representative late-time SNR at 0^\circ [dB]', 0.16);
linkaxes([ax_sigma, ax_full, ax_late], 'x');

drawnow;
main_position = ax_sigma.Position;
ax_full.Position([1, 3]) = main_position([1, 3]);
ax_late.Position([1, 3]) = main_position([1, 3]);
end


function fig = plot_outcomes(number, criteria, sigmas)
%PLOT_OUTCOMES  Stacked outcome fractions per sigma, one panel per rule.

outcome_labels = {'kept, correct', 'kept, wrong', ...
    'rejected, would have been correct', ...
    'rejected, would have been wrong'};
outcome_colors = [0.20, 0.55, 0.30; 0.80, 0.25, 0.20; ...
    0.60, 0.80, 0.60; 0.95, 0.70, 0.65];

fig = figure(number); clf;
fig.Name = 'Outcome accounting';
fig.Color = 'w';
fig.Position(3:4) = [1400, 480];
layout = tiledlayout(fig, 1, numel(criteria), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for k = 1:numel(criteria)
    rows = criteria(k).rows;
    outcome = [rows.ncorrect_committed, rows.nwrong_committed, ...
        rows.nlost, rows.navoided] ./ rows.ntrials;
    ax = nexttile(layout);
    bars = bar(ax, 1:numel(sigmas), outcome, 'stacked', 'BarWidth', 0.7);

    for ib = 1:numel(bars)
        bars(ib).FaceColor = outcome_colors(ib, :);
        bars(ib).EdgeColor = 'w';
    end

    ax.XTick = 1:numel(sigmas);
    ax.XTickLabel = sigma_tick_labels(sigmas);
    ax.YLim = [0, 1];
    ax.FontSize = 13;
    grid(ax, 'on');
    ax.XGrid = 'off';
    title(ax, criteria(k).title, 'FontWeight', 'normal');
    xlabel(ax, 'Noise Level (\sigma)');

    if k == 1
        ylabel(ax, 'Fraction of samples');
    else
        ax.YTickLabel = [];
    end
end

lgd = legend(bars, outcome_labels, 'Orientation', 'horizontal', ...
    'FontSize', 12);
lgd.Layout.Tile = 'south';
end


function fig = plot_angle_outcomes(number, by_receiver, sigmas, alpha, ...
        bin_width)
%PLOT_ANGLE_OUTCOMES  Statistical-gate outcome fractions vs receiver angle.

outcome_names = {'rejected_correct', 'rejected_wrong', 'kept_correct', ...
    'kept_wrong'};
outcome_labels = {'rejected, would have been correct', ...
    'rejected, would have been wrong', 'kept, correct', 'kept, wrong'};
outcome_styles = {'--', '--', '-', '-'};
outcome_colors = [0.20, 0.55, 0.30; 0.80, 0.25, 0.20; ...
    0.20, 0.55, 0.30; 0.80, 0.25, 0.20];
edges = -180:bin_width:180;
centers = edges(1:end-1) + bin_width / 2;

fig = figure(number); clf;
fig.Name = 'Statistical outcomes vs angle';
fig.Color = 'w';
fig.Position(3:4) = [1500, 700];
layout = tiledlayout(fig, 2, ceil(numel(sigmas) / 2), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf(['Statistical rejection, \\alpha = %g: ', ...
    'outcome fractions per %d^\\circ sector'], alpha, bin_width), ...
    'FontSize', 14);

for is = 1:numel(sigmas)
    rows = by_receiver(by_receiver.sigma == sigmas(is) & ...
        abs(by_receiver.familywise_alpha - alpha) < 1e-9 * alpha, :);

    if isempty(rows)
        error('plot_angle_outcomes:Missing', ...
            'No rows for alpha %g. Add it to the analysis grid.', alpha);
    end

    bin = discretize(rows.angle_deg, edges, 'IncludedEdge', 'right');
    bin(rows.angle_deg == -180) = 1;
    ntrials = accumarray(bin, rows.ntrials, [numel(centers), 1]);

    ax = nexttile(layout);
    hold(ax, 'on');
    lines = gobjects(1, numel(outcome_names));

    for k = 1:numel(outcome_names)
        fraction = accumarray(bin, rows.(outcome_names{k}), ...
            [numel(centers), 1]) ./ ntrials;
        lines(k) = plot(ax, centers, fraction, outcome_styles{k}, ...
            'Color', outcome_colors(k, :), 'LineWidth', 2);
    end

    hold(ax, 'off');
    ax.XLim = [-180, 180];
    ax.XTick = -180:90:180;
    ax.YLim = [0, 1];
    ax.FontSize = 12;
    grid(ax, 'on');
    title(ax, sprintf('\\sigma = %s', sigma_tick_labels(sigmas(is))), ...
        'FontWeight', 'normal');
    xlabel(ax, 'Receiver angle [deg] (0 = nose-on)');
    ylabel(ax, 'Fraction of samples');
end

lgd = legend(lines, outcome_labels, 'Orientation', 'horizontal', ...
    'FontSize', 12);
lgd.Layout.Tile = 'south';
end


function plot_sigmas = display_sigmas(sigmas)
%DISPLAY_SIGMAS  Give sigma = 0 a display-only spot on the log axis.

plot_sigmas = sigmas;
plot_sigmas(sigmas == 0) = min(sigmas(sigmas > 0)) / 10;
end


function labels = sigma_tick_labels(sigma_values)
%SIGMA_TICK_LABELS Format fixed logarithmic ticks, including sigma = 0.

labels = strings(size(sigma_values));

for ivalue = 1:numel(sigma_values)
    if sigma_values(ivalue) == 0
        labels(ivalue) = "0";
    else
        labels(ivalue) = sprintf('%.0e', sigma_values(ivalue));
    end
end
end


function labels = snr_tick_labels(snr_values)
%SNR_TICK_LABELS Format SNR values for the two upper axes.

labels = strings(size(snr_values));

for ivalue = 1:numel(snr_values)
    if isinf(snr_values(ivalue))
        labels(ivalue) = "Inf";
    else
        labels(ivalue) = sprintf('%.1f', snr_values(ivalue));
    end
end
end


function ax = top_snr_axis(base_axis, x_ticks, snr_values, label, offset)
%TOP_SNR_AXIS Add a negligible-height x-axis above the main plot.

position = base_axis.Position;
position = [position(1), position(2) + position(4) + offset, ...
    position(3), 1e-12];
ax = axes( ...
    'Parent', base_axis.Parent, ...
    'Position', position, ...
    'Color', 'none', ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'YColor', 'none', ...
    'YTick', [], ...
    'XScale', 'log', ...
    'XLim', base_axis.XLim, ...
    'XTick', x_ticks, ...
    'XTickLabel', snr_tick_labels(snr_values), ...
    'FontSize', base_axis.FontSize, ...
    'Box', 'off');
xlabel(ax, label, 'FontSize', 14);
end
