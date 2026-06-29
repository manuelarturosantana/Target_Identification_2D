% generates data for 5 twinjet and 5 quadjet geometries
% The first interesting thing is that changing n_polar and n messed up the
% geometry.
% 
% 
% load("data_gen_planes_check/quadjet_plane_1.mat")
% pols2use = pols;
% load("./data_gen_planes_check/twinjet_plane4.mat")
%%
% figure(1)
% clf
% plot(ts,uff_all(1,:))


%%
% paths
% addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
% addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
% addpath(genpath("/home/vhojas/chebfun"))

% ------------------------------------------------------------
% Output directory
% ------------------------------------------------------------
save_data_dir = "./data_gen_planes_check";
if ~exist(save_data_dir, "dir")
    mkdir(save_data_dir);
end

% ------------------------------------------------------------
% Experiment params
% ------------------------------------------------------------
w_0              = 50;            % center frequency of the Gaussian
freq_gauss_width = 5;             % Gaussian frequency width
gauss_picker_tol = 1e-10;         % Gaussian truncation tolerance
numw             = 400;           % number of points for frequency integrals
t_end            = 400;           % final time
t_offset         = 30;            % time offset
num_receivers_per_obs = 1;
are_recivers_random = true;



%%
% receiver angles
angle1 = -5*pi/180;
angle2 =  5*pi/180;

wimag  = -0.3;

% time-domain and AAA derived params
[sigmas, wlims] = gauss_picker(w_0, freq_gauss_width, ...
    gauss_picker_tol);
ws = linspace(wlims(1), wlims(2), 1000);
psp = problem_data( ...
    'wlims', wlims, ...
    'sigmas', sigmas, ...
    'w_0', w_0, ...
    'xs', 0, ...
    'ys', 0, ...
    't_0', t_offset);

% time interval and number of time samples
lambda = (2 * pi) / wlims(2);
tlims = [0, t_end];
numt = ceil((tlims(2) / lambda) * 15);

%% 
vals = [];
for ii = 1:length(ws)
    vals(ii) = psp.B_w(ws(ii));
end
% 
figure(1)
clf
semilogy(ws,abs(vals))
title("Frequency Domain Gaussian Decay")
xlabel("ws")
ylabel("Magnitude of Gaussian")


ts = linspace(1,t_end,30000);
vals = psp.a_t(ts);
figure(2)
clf
% semilogy(ts,abs(vals))
plot(ts,real(vals))
title("Time Domain Gaussian Decay")
xlabel("ws")
ylabel("Magnitude of Gaussian")

return
%%
% ------------------------------------------------------------
% Solver params
% ------------------------------------------------------------
params = build_params( ...
    'k', w_0 + freq_gauss_width + 30, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

% params = build_params( ...
%     'k', w_0 + freq_gauss_width, ...
%     'formulation', 'sl', ...
%     'solver', 'direct', ...
%     'delta', 0.05, ...
%     'npolar', 160, ...
%     'n', 24, ...
%     'p', 6, ...
%     'p_edge', 2);

% ------------------------------------------------------------
% Build all geometries
% ------------------------------------------------------------
num_presets = 5;

all_geometries = cell(1, 1);
all_geometry_names = cell(1, 1);

idx = 1;

% twinjet presets
% preset_name = sprintf("plane%d",4);
% gparams = build_twinjet_params('preset', preset_name);
% all_geometries{1} = twinjet_plane_geometry(gparams);
% all_geometry_names{1} = sprintf("twinjet_%s", preset_name);

for i = 1:num_presets
    preset_name = sprintf("plane%d", i);
    gparams = build_twinjet_params('preset', preset_name);

    all_geometries{idx} = twinjet_plane_geometry(gparams);
    all_geometry_names{idx} = sprintf("twinjet_%s", preset_name);

    idx = idx + 1;
end

%quadjet presets
for i = 1:num_presets
    preset_name = sprintf("plane%d", i);
    gparams = build_quadjet_params('preset', preset_name);

    all_geometries{idx} = quadjet_plane_geometry(gparams);
    all_geometry_names{idx} = sprintf("quadjet_%s", preset_name);

    idx = idx + 1;
end

% ------------------------------------------------------------
% Run experiment
% ------------------------------------------------------------
start = tic;

% for ii = 1:length(all_geometries)
for ii = 1:1
    geom = all_geometries{ii};

    fprintf("\n==================================================\n");
    fprintf("Running object %d / %d: %s\n", ii, ...
        length(all_geometries), all_geometry_names{ii});
    fprintf("==================================================\n");

    lp = RP2LP(params, geom);

    % total number of points in the curve
    N = lp.RP_Curve.global.N;

    % rule of thumb for computing complex resonances
    w_len = wlims(2) - wlims(1);
    num_sub_ints = ceil(w_len / 2.0);
    jump = w_len / num_sub_ints;

    pols_wp = [];
    pstart = tic;

    for jj = 1:num_sub_ints
        wlims_loc = [wlims(1) + (jj - 1) * jump, ...
                     wlims(1) + jj * jump];

        %dummy problem data for pole computation
        psp = problem_data( ...
            'wlims', wlims_loc, ...
            'w_0', (wlims_loc(2) + wlims_loc(1)) / 2.0, ...
            'numw', numw, ...
            'xs', 0, ...
            'ys', 0, ...
            'wimag', wimag, ...
            'N', N);

        pols_wp_loc = comp_poles(psp, lp);
        pols_wp_loc = sort(pols_wp_loc, ...
            "comparisonMethod", "real");
        pols_wp = [pols_wp; pols_wp_loc(:)];
    end
    % 
    pol_comp_time = toc(pstart)
    % pols_wp = pols2use;

    % dummy problem for LU decompositions
    psp = problem_data( ...
        'wlims', wlims, ...
        'w_0', w_0, ...
        'numw', numw, ...
        'xs', 0, ...
        'ys', 0, ...
        'wimag', wimag, ...
        'N', N);

    tic
    [L_rl, U_rl] = comp_bie_LU(psp.ws, lp);
    LU_time_real_line = toc

    tic
    [L_pole, U_pole] = comp_pole_LU(psp, lp, pols_wp);
    LU_time_pols = toc

    if are_recivers_random
        [xs, ys] = ff_points_random(angle1, angle2, ...
            num_receivers_per_obs);
    else
        [xs, ys] = ff_points_equally_spaced(angle1, angle2, ...
            num_receivers_per_obs);
    end

    uff_all = [];
    r_res_all = [];
%%
    for rind = 1:num_receivers_per_obs
        x = xs(rind);
        y = ys(rind);

        psp = problem_data( ...
            'wlims', wlims, ...
            'N', N, ...
            'Tlims', tlims, ...
            'kappa', [-x, -y], ...
            'numt', numt, ...
            'numw', numw, ...
            'w_0', w_0, ...
            't_0', t_offset, ...
            'sigmas', sigmas, ...
            'wimag', wimag, ...
            'is_far_field', true, ... // to make false some small things in the solver need to be changed. See RP2LP.m
            'xs', x, ...
            'ys', y);
% This to check if we get the early bumps in the near field.
         % psp = problem_data( ...
         %    'wlims', wlims, ...
         %    'N', N, ...
         %    'Tlims', tlims, ...
         %    'kappa', [-x, -y], ...
         %    'numt', numt, ...
         %    'numw', numw, ...
         %    'w_0', w_0, ...
         %    't_0', t_offset, ...
         %    'sigmas', sigmas, ...
         %    'wimag', wimag, ...
         %    'is_far_field', false, ... // to make false some small things in the solver need to be changed. See RP2LP.m
         %    'xs', x * 4, ...
         %    'ys', y * 4, ...
         %    'is_rp_curve', true, ...
         %    'is_open_curve', false);

        tic
        % [uff, ~, r_res] = cts_wpols_LU(psp, lp, pols_wp, ...
        %     L_rl, U_rl, L_pole, U_pole);
        [usol, scat_sol, r_res, ~, smoothint, f_sols, zeroint] = cts_wpols_LU(psp, lp, pols_wp, ...
            L_rl, U_rl, L_pole, U_pole);
        single_reciever_comp_time = toc

        % uff = uff(:).';
        % r_res = r_res(:).';
        % 
        % uff_all = [uff_all; uff];
        r_res_all = [r_res_all; r_res];
    end

    ts = psp.ts;
    curveX = lp.curve.X;
    curveY = lp.curve.Y;
    pols = pols_wp;
    % 
    % save(fullfile(save_data_dir, all_geometry_names{ii} + ".mat"), ...
    %     "r_res_all", "uff_all", "pols", "ts", "xs", "ys", ...
    %     "curveX", "curveY");
end

total_time = toc(start)


%%

figure(3)
clf
% plot(ts,uff_all(1,:))
plot(psp.ts,usol)

figure(4)
clf
plot(psp.ws, f_sols,'-o')

sub_sol = f_sols;
for ii = 1:length(pols)
    for jj = 1:length(psp.ws)
        sub_sol(jj) = sub_sol(jj) -  (r_res(ii) / (psp.ws(jj) - pols(ii)));
    end
end

figure(5)
clf
plot(psp.ws, sub_sol,'-o')
return
% figure(5)
% clf
% plot(psp.ws, smoothint)

%%
for ii = 1:length(pols)
    S = svd(lp.bie_mat(pols(ii)));
    S(end)
end
