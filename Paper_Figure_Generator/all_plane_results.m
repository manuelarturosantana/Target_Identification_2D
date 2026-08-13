% evaluate_all_planes_accuracy.m
clear;
results_file = 'Identification/classification_results/plane_results_averaged.mat';

% Define return angle half-widths in degrees to test
angle_bounds = [15, 30, 60, 90]; 
num_bounds   = length(angle_bounds);

results = load(results_file);
results = results.results;
trials  = results.trials;

sigmas = unique(trials.sigma);
nsigma = numel(sigmas);

all_obstacles = unique(trials.true_idx);
num_obstacles = numel(all_obstacles);

% Convert receiver angles to degrees and wrap to [-180, 180] range
rec_angles_deg = rad2deg(results.snr.receiver_angle);
angles_wrapped = mod(rec_angles_deg + 180, 360) - 180;

%% 1. Preallocate Accuracy Data Matrices
% Overall accuracy pooled over all obstacles: (nsigma x num_bounds+1)
overall_accuracy = zeros(nsigma, num_bounds + 1);

% Individual plane accuracy: (nsigma x num_bounds+1 x num_obstacles)
per_plane_accuracy = zeros(nsigma, num_bounds + 1, num_obstacles);

%% 2. Calculate Accuracy across All Noise Levels, Angle Ranges, and Planes
for k = 1:nsigma
    % A. Calculate bounded angle accuracies
    for r = 1:num_bounds
        mask_r = abs(angles_wrapped) <= angle_bounds(r);
        pooled_correct_r = [];
        
        for obs_i = 1:num_obstacles
            obs_idx = all_obstacles(obs_i);
            idx = (trials.sigma == sigmas(k)) & (trials.true_idx == obs_idx);
            correct_k = trials.correct(idx);
            
            % Save per-plane accuracy for range r
            per_plane_accuracy(k, r, obs_i) = mean(correct_k(mask_r));
            pooled_correct_r = [pooled_correct_r; correct_k(mask_r)];
        end
        
        % Pool all planes together for range r
        overall_accuracy(k, r) = mean(pooled_correct_r);
    end
    
    % B. Calculate full sweep (All angles) accuracy
    pooled_correct_all = [];
    for obs_i = 1:num_obstacles
        obs_idx = all_obstacles(obs_i);
        idx = (trials.sigma == sigmas(k)) & (trials.true_idx == obs_idx);
        correct_k = trials.correct(idx);
        
        per_plane_accuracy(k, num_bounds + 1, obs_i) = mean(correct_k);
        pooled_correct_all = [pooled_correct_all; correct_k];
    end
    overall_accuracy(k, num_bounds + 1) = mean(pooled_correct_all);
end

%% Plot Setup Configurations
plot_styles = {'-o', '-s', '-^', '-d', '--k*'};
range_labels = { ...
    '[-15^\circ, 15^\circ]', ...
    '[-30^\circ, 30^\circ]', ...
    '[-60^\circ, 60^\circ]', ...
    '[-90^\circ, 90^\circ]', ...
    '[-180^\circ, 180^\circ]'};

%% Figure 1: Aggregate Performance across ALL Planes
fig1 = figure(101); clf;
fig1.Name = 'Overall Classification Accuracy across All Planes';

semilogx( ...
    sigmas, overall_accuracy(:, 1), plot_styles{1}, ...
    sigmas, overall_accuracy(:, 2), plot_styles{2}, ...
    sigmas, overall_accuracy(:, 3), plot_styles{3}, ...
    sigmas, overall_accuracy(:, 4), plot_styles{4}, ...
    sigmas, overall_accuracy(:, 5), plot_styles{5}, ...
    'LineWidth', 2.5);

grid on;
xlabel("Noise Level (\sigma)", "FontSize", 14);
ylabel("Overall Correct Classification Rate", "FontSize", 14);
% title("Method Effectiveness Pooled Across All 10 Aircraft Geometries", "FontSize", 16);
legend(range_labels, 'Location', 'SouthWest', 'FontSize', 12);
ax = gca; ax.FontSize = 14;

%% Figure 2: Plane-by-Plane Breakdown Grid (2x5 Layout)
fig2 = figure(102); clf;
fig2.Name = 'Per-Plane Accuracy Breakdown';
tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

for obs_i = 1:num_obstacles
    obs_idx = all_obstacles(obs_i);
    if obs_idx <= 5
        plane_title = sprintf("Quadjet %d", obs_idx);
    else
        plane_title = sprintf("Twinjet %d", obs_idx - 5);
    end
    
    nexttile(obs_i);
    
    semilogx( ...
        sigmas, per_plane_accuracy(:, 1, obs_i), plot_styles{1}, ...
        sigmas, per_plane_accuracy(:, 2, obs_i), plot_styles{2}, ...
        sigmas, per_plane_accuracy(:, 3, obs_i), plot_styles{3}, ...
        sigmas, per_plane_accuracy(:, 4, obs_i), plot_styles{4}, ...
        sigmas, per_plane_accuracy(:, 5, obs_i), plot_styles{5}, ...
        'LineWidth', 1.5);
        
    grid on;
    title(plane_title, 'FontSize', 12);
    xlim([min(sigmas), max(sigmas)]);
    ylim([0, 1.05]);
    
    % Label outer edges only for clean presentation
    if obs_i == 1 || obs_i == 6
        ylabel("Accuracy", "FontSize", 12);
    else
        yticklabels([]);
    end
    
    if obs_i >= 6
        xlabel("Noise Level (\sigma)", "FontSize", 11);
    else
        xticklabels([]);
    end
    
    % Add legend to the last tile
    if obs_i == num_obstacles
        legend(range_labels, 'Location', 'SouthWest', 'FontSize', 9);
    end
    
    ax = gca; ax.FontSize = 11;
end

sgtitle("Classification Accuracy by Angle Range for Each Plane Geometry", "FontSize", 16);