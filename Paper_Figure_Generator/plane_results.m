% read_averaged_classification_results.m
%
clear;
results_file   = 'Identification/classification_results/plane_results_averaged.mat';
tiledlayout(2,2);

loop_cnt = 0; % Track loop iteration for transposing the layout
% obs_indicies = [3,9]; 
obs_indicies = [1,6]; % Failure example
for obstacle_idx = obs_indicies
loop_cnt = loop_cnt + 1;
obstacle_idx
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
accuracy = zeros(nsigma, 1);
ntrials  = zeros(nsigma, 1);
nwrong   = zeros(nsigma, 1);
correct_inds = zeros(nsigma-1, 359);
%%
inds = trials.true_idx == obstacle_idx;
for k = 1:nsigma
    idx = (trials.sigma == sigmas(k)) & inds;
    if k~=1
        correct_inds(k-1, :) = trials.correct(idx);
    end
    ntrials(k) = nnz(idx);
    nwrong(k) = nnz(~trials.correct(idx));
    accuracy(k) = mean(trials.correct(idx));
end
%% Plots

grid on

nexttile(loop_cnt) % Transposed: Forces the first plot into Tile 1 (Iter 1) or Tile 2 (Iter 2)

legend_data = {};
hold on
for ii = 1:length(results.snr.sigma) - 1
    plot(rad2deg(results.snr.receiver_angle), snr_data(ii+1,:,obstacle_idx))
    legend_data{end+1} = "$\sigma$ = " + num2str(results.snr.sigma(ii+1));
    inds_missed = ~correct_inds(ii,:);
    plot(rad2deg(results.snr.receiver_angle(inds_missed)), snr_data(ii+1,inds_missed,obstacle_idx),'ko')
    if sum(~correct_inds(ii,:) ~= 0)
        legend_data{end+1} = "";
    end
end
hold off
if obstacle_idx == obs_indicies(2)
    ldg = legend(legend_data,'Location','SouthEast','interpreter','latex');
end

xlabel("Scatterer Orientation in Degrees")
ylabel("Late Time SNR")
xlim([0,360])
plane_name
title(plane_name)
ax = gca; ax.set("Fontsize",15)

if obstacle_idx == obs_indicies(2)
    ldg.FontSize = 10;
end

grid on

nexttile(loop_cnt + 2) % Transposed: Forces the second plot into Tile 3 (Iter 1) or Tile 4 (Iter 2)

semilogx(sigmas, accuracy,'-o','linewidth',2)
grid on
xlabel("Noise Level")
ylabel("Overall Accuracy")
ax = gca; ax.set("FontSize",15)
end


