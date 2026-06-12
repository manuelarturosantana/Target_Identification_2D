% read_averaged_classification_results.m


clear;

results_file   = 'Identification/classification_results/plane_results_averaged.mat';
objects_folder = '/scratch/msantana/plane_data_360';

% for obstacle_idx = 1:10
for obstacle_idx = 10

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

fig = figure(1);
clf
tiledlayout(1,3);
nexttile
curve = build_plane(obstacle_idx);

opts = struct(); opts.show_nodes= false; opts.plot_color='black'; opts.line_width=2;
opts.no_gaps = true;
plot_curve(curve, opts);
grid on

nexttile

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

legend(legend_data,'Location','SouthEast','interpreter','latex')
xlabel("Reciever Angle")
ylabel("Late Time SNR")
ax = gca; ax.set("Fontsize",15)
grid on

nexttile



semilogx(sigmas, accuracy,'-o','linewidth',2)
grid on
xlabel("Noise Level")
ylabel("Correct Classification Rate over All Angles")
ax = gca; ax.set("FontSize",15)

%%
if obstacle_idx <= 5
    fig_name = sprintf("twinjet_%d", obstacle_idx);
else
    fig_name = sprintf("quadjet_%d", obstacle_idx - 5);
end

savefig("Plane_Classification_Figures/" + fig_name + ".fig");
exportgraphics(fig, "Plane_Classification_Figures/" + fig_name + ".png")
end
%%
function Curve = build_plane(idx)
    params = build_params( ...
    'k', 50, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 160, ...
    'n', 20, ...
    'p', 3, ...
    'p_edge', 2);

  

    if idx <= 5
        preset_name = sprintf("plane%d", idx);
        gparams     = build_twinjet_params('preset', preset_name);
        geometry    = twinjet_plane_geometry(gparams);
    else
        preset_name = sprintf("plane%d", idx - 5);
        gparams     = build_quadjet_params('preset', preset_name);
        geometry    = quadjet_plane_geometry(gparams);
    end

    Curve = build_curve(geometry, params);
    Curve = refine_curve_wavenumber(Curve, 2*pi/params.k, 1);

end

