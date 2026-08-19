% read_averaged_classification_results.m
%

% raw data lives in /scratch/msantana/plane_data_360
clear;
results_file = 'Identification/classification_results/plane_results_averaged.mat';

% Define angle range half-widths in degrees
angle_bounds = [15, 30, 60, 90];
num_bounds   = length(angle_bounds);

tiledlayout(2,2);

loop_cnt = 0; % Track loop iteration for transposing the layout
obs_indicies = [6, 9]; % Failure example
ax_top = gobjects(1, 2); % Store axes handles for top row to sync limits

for obstacle_idx = obs_indicies
    loop_cnt = loop_cnt + 1;
    
    if obstacle_idx < 6
        plane_name = "Quadjet " + num2str(obstacle_idx);
    else
        plane_name = "Twinjet " + num2str(obstacle_idx - 5);
    end

    snr_mode = 'late'; % 'early' or 'late'
    show_legend = true;
    results = load(results_file);
    results = results.results;
    
    if strcmp(snr_mode,'late')
        snr_data = results.snr.snr_late_db;
    else
        snr_data = results.snr.snr_full_db;
    end
    
    trials = results.trials;
    sigmas = unique(trials.sigma);
    nsigma = numel(sigmas);
    
    accuracy     = zeros(nsigma, 1);
    ntrials      = zeros(nsigma, 1);
    nwrong       = zeros(nsigma, 1);
    correct_inds = zeros(nsigma-1, 359);
    
    % Store accuracy across angle ranges:
    % Columns: [-15, 15], [-30, 30], [-60, 60], [-90, 90], All angles
    accuracy_ranges = zeros(nsigma, num_bounds + 1);

    % Convert receiver angles to degrees and wrap to [-180, 180] range
    rec_angles_deg = rad2deg(results.snr.receiver_angle);
    angles_wrapped = mod(rec_angles_deg + 180, 360) - 180;

    %% Data Extraction Loop
    inds = trials.true_idx == obstacle_idx;
    for k = 1:nsigma
        idx = (trials.sigma == sigmas(k)) & inds;
        if k ~= 1
            correct_inds(k-1, :) = trials.correct(idx);
        end
        ntrials(k)  = nnz(idx);
        nwrong(k)   = nnz(~trials.correct(idx));
        accuracy(k) = mean(trials.correct(idx));

        correct_k = trials.correct(idx);

        % Calculate accuracy for each angle range (-bound, bound)
        for r = 1:num_bounds
            mask = abs(angles_wrapped) <= angle_bounds(r);
            accuracy_ranges(k, r) = mean(correct_k(mask));
        end
        % Full sweep accuracy
        accuracy_ranges(k, num_bounds + 1) = accuracy(k);
    end

    %% Subplot 1 & 2 (Top Row): SNR vs Receiver Angle
    ax_top(loop_cnt) = nexttile(loop_cnt); % Transposed: Tile 1 or Tile 2

    % Predefine explicit colors for consistent matching across plots
    num_snr_curves = length(results.snr.sigma) - 1;
    snr_colors = lines(num_snr_curves);

    % Sort by wrapped angle so the line plots cleanly over [-180, 180]
    [angles_sorted, sort_idx] = sort(angles_wrapped);

    legend_data = {};
    hold on
    for ii = 1:num_snr_curves
        % Plot line using explicit color array
        plot(angles_sorted, snr_data(ii+1, sort_idx, obstacle_idx), ...
             'Color', snr_colors(ii, :), 'LineWidth', 1.5);
        
        legend_data{end+1} = "$\sigma$ = " + num2str(results.snr.sigma(ii+1));
        
        inds_missed = ~correct_inds(ii, sort_idx);
        plot(angles_sorted(inds_missed), snr_data(ii+1, sort_idx(inds_missed), obstacle_idx), 'ko');
        
        if sum(~correct_inds(ii,:) ~= 0)
            legend_data{end+1} = "";
        end
    end
    hold off

    if obstacle_idx == obs_indicies(2)
        ldg = legend(legend_data, 'Location', 'SouthEast', 'interpreter', 'latex');
        ldg.FontSize = 10;
    end

    xlabel("Scatterer Orientation in Degrees")
    
    % Apply Y-axis label only on Tile 1; remove ticks on Tile 2
    if loop_cnt == 1
        ylabel("Late Time SNR")
    else
        yticklabels([])
    end

    xlim([-180, 180])
    title(plane_name)
    ax = gca; ax.FontSize = 15;
    grid on

    %% Subplot 3 & 4 (Bottom Row): Accuracy vs Noise Level by Angle Range
    nexttile(loop_cnt + 2) % Transposed: Tile 3 or Tile 4

    plot_styles = {'-o', '-s', '-^', '-d', '--k*'};
    range_labels = { ...
        '[-15^\circ, 15^\circ]', ...
        '[-30^\circ, 30^\circ]', ...
        '[-60^\circ, 60^\circ]', ...
        '[-90^\circ, 90^\circ]', ...
        '[-180^\circ, 180^\circ]'};

    % Single semilogx call to preserve log scaling robustly
    semilogx( ...
        sigmas, accuracy_ranges(:, 1), plot_styles{1}, ...
        sigmas, accuracy_ranges(:, 2), plot_styles{2}, ...
        sigmas, accuracy_ranges(:, 3), plot_styles{3}, ...
        sigmas, accuracy_ranges(:, 4), plot_styles{4}, ...
        sigmas, accuracy_ranges(:, 5), plot_styles{5}, ...
        'LineWidth', 2);

    grid on
    xlabel("Noise Level")
    ylabel("Accuracy")
    ax = gca; ax.FontSize = 15;

    if obstacle_idx == obs_indicies(2)
        legend(range_labels, 'Location', 'best')
    end
end

% Link Y-axis limits across top row
linkaxes(ax_top, 'y');