% read_averaged_classification_results.m

%%
clear;

results_file = 'Identification/classification_results/plane_results_averaged.mat';

for obstacle_idx = 1:10
% for obstacle_idx = 1

snr_mode = 'late'; % 'early' or 'late'
show_legend = true;

% Define angle ranges in degrees centered around 0/360: (-15,15), (-30,30), (-60,60), (-90,90)
range_limits = [15, 30, 60, 90]; 
num_ranges = length(range_limits);

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

accuracy_all    = zeros(nsigma, 1);
accuracy_ranges = zeros(nsigma, num_ranges);
ntrials         = zeros(nsigma, 1);
nwrong          = zeros(nsigma, 1);
correct_inds    = zeros(nsigma-1, 359);

%%
inds = trials.true_idx == obstacle_idx;
rec_deg = rad2deg(results.snr.receiver_angle(:)); % Receiver angles in degrees (0 to 360)

for k = 1:nsigma
    idx = (trials.sigma == sigmas(k)) & inds;
    if k ~= 1
        correct_inds(k-1, :) = trials.correct(idx);
    end
    ntrials(k) = nnz(idx);
    nwrong(k) = nnz(~trials.correct(idx));

    correct_k = trials.correct(idx);
    correct_k = correct_k(:); % ensure column vector
    
    % Classification accuracy over all angles
    accuracy_all(k) = mean(correct_k);

    % Classification accuracy over specified angular ranges
    for r = 1:num_ranges
        lim = range_limits(r);
        % Mask for angles in [0, lim] U [360 - lim, 360]
        range_mask = (rec_deg <= lim) | (rec_deg >= 360 - lim);
        accuracy_ranges(k, r) = mean(correct_k(range_mask));
    end
end


%% Plots

fig = figure(obstacle_idx);
clf
tiledlayout(1,3);
nexttile
% curve = build_plane(obstacle_idx);
% 
% opts = struct(); opts.show_nodes = false; opts.plot_color = 'black'; opts.line_width = 2;
% opts.no_gaps = true;
% plot_curve2(curve, opts);
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
xlabel("Receiver Angle")
ylabel("Late Time SNR")
xlim([0,360])
ax = gca; ax.set("FontSize",15)
grid on

nexttile

% Plot accuracy curves across noise levels for each angle range
semilogx(sigmas, accuracy_all, '-o', 'LineWidth', 2)
hold on
markers = {'-s', '-^', '-d', '-v'};
for r = 1:num_ranges
    semilogx(sigmas, accuracy_ranges(:, r), markers{r}, 'LineWidth', 2)
end
hold off
grid on
xlabel("Noise Level")
ylabel("Correct Classification Rate")
legend({'All angles', ...
        '(-15^\circ, 15^\circ)', ...
        '(-30^\circ, 30^\circ)', ...
        '(-60^\circ, 60^\circ)', ...
        '(-90^\circ, 90^\circ)'}, ...
       'Location', 'best', 'Interpreter', 'latex')
ax = gca; ax.set("FontSize",15)

%%
if obstacle_idx <= 5
    fig_name = sprintf("quadjet_%d", obstacle_idx);
else
    fig_name = sprintf("twinjet_%d", obstacle_idx - 5);
end

savefig("Plane_Classification_Figures/" + fig_name + ".fig");
exportgraphics(fig, "Plane_Classification_Figures/" + fig_name + ".png")
end

return
%% Bonus to figure out why the quadjet is miss classifying
inds1 = trials.true_idx == 6;
inds2 = trials.sigma == 1e-6;
trials.pred_idx(inds1 & inds2)
sum(trials.correct(inds1 & inds2)) / 360

%%
inds3 = trials.pred_idx == 1;
sum(trials.pred_idx(inds1 & inds2 & inds3))
360 - sum(trials.correct(inds1 & inds2))

%%
load('/scratch/msantana/plane_data_360/quadjet_plane3.mat')
figure(1)
clf
plot(pols,'*')
hold on
load('/scratch/msantana/plane_data_360/twinjet_plane1.mat')
plot(pols,'^')
hold off

%% No noise classification error
indx = trials.sigma == 0;
percent_correct = sum(trials.correct(indx)) / length(trials.correct(indx))

%% Who mis classifies at 0?
inds1 = trials.sigma == 0;
inds2 = trials.correct == false;
trials.true_idx(inds1 & inds2)
trials.pred_idx(inds1 & inds2)
mod(find(inds1 & inds2),360)

%%
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
        gparams     = build_quadjet_params('preset', preset_name);
        geometry    = quadjet_plane_geometry(gparams);
    else
        preset_name = sprintf("plane%d", idx - 5);
        gparams     = build_twinjet_params('preset', preset_name);
        geometry    = twinjet_plane_geometry(gparams);
   
    end

    Curve = build_curve(geometry, params);
    Curve = refine_curve_wavenumber(Curve, 2*pi/params.k, 1);

end

function plot_curve2(Curve, opts)
if nargin < 2
    opts = struct();
end

if ~isfield(opts, "nplot"), opts.nplot = 200; end
if ~isfield(opts, "show_nodes"), opts.show_nodes = true; end
if ~isfield(opts, "show_labels"), opts.show_labels = false; end
if ~isfield(opts, "node_size"), opts.node_size = 18; end
if ~isfield(opts, "line_width"), opts.line_width = 1.25; end
if ~isfield(opts, "no_gaps"), opts.no_gaps = false; end

patches = Curve.patches;

hold_state = ishold;
hold on;

for j = 1:numel(patches)
    P = patches{j};

    tmin = min(P.tj);
    tmax = max(P.tj);
    if (opts.no_gaps)
        t = linspace(-1, 1, opts.nplot);
        if (opts.show_nodes)
            warning("plot_curve: Show nodes and no_gaps both set to true in opts. Setting no_gaps=true does not use discretization nodes.")
        end
    else
        t = linspace(tmin, tmax, opts.nplot);
    end
    

    x = P.x(t);
    y = P.y(t);

    if ~isfield(opts, "plot_color")
        plot(x, y, "-", "LineWidth", opts.line_width);
    else
        plot(x, y, "-", "LineWidth", opts.line_width,'Color',opts.plot_color);
    end

    if opts.show_nodes
        plot(P.xj, P.yj, ".", "MarkerSize", opts.node_size);
    end

    if opts.show_labels
        tm = 0.5 * (tmin + tmax);
        xm = P.x(tm);
        ym = P.y(tm);
        text(xm, ym, sprintf("%d", j), ...
            "VerticalAlignment", "bottom", ...
            "HorizontalAlignment", "center");
    end
end

axis equal;
grid on;

if ~hold_state
    hold off;
end

end