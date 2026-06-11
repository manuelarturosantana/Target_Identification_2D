% read_averaged_classification_results.m


clear;

results_file   = 'Identification/classification_results/plane_results_averaged.mat';
objects_folder = '/scratch/msantana/plane_data_360';
obstacle_idx = 4;
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

figure(1)
clf
tiledlayout(1,2);
nexttile

hold on
for ii = 1:length(results.snr.sigma) - 1
    plot(rad2deg(results.snr.receiver_angle), snr_data(ii+1,:,obstacle_idx))
    inds_missed = ~correct_inds(ii,:);
    plot(rad2deg(results.snr.receiver_angle(inds_missed)), snr_data(ii+1,inds_missed,obstacle_idx),'ko')
end
hold off

legend_data = {};
for ii = 1:length(results.snr.sigma) - 1
    legend_data{end+1} = "$\sigma$ = " + num2str(results.snr.sigma(ii+1));
    legend_data{end+1} = "";
end
legend(legend_data,'interpreter','latex')
xlabel("Reciever Angle")
ylabel("Late Time SNR")
grid on

nexttile



semilogx(sigmas, accuracy,'-o')
xlabel("Noise Level")
ylabel("Accuracy")
ax = gca; ax.set("FontSize",15)

