% generates data for 5 twinjet
% Now we use a refined curve to see if it makes the data better.

% paths
addpath(genpath('/home/msantana/Progamming/rp2d_matlab/src'))
addpath(genpath('/home/msantana/Progamming/Target_Identification_2D'))
addpath(genpath('/home/msantana/Progamming/MATLAB_PACKAGES'))

% ------------------------------------------------------------
% Output directory
% ------------------------------------------------------------
save_data_dir = "./plane_data";
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
num_receivers_per_obs = 359;
are_recivers_random = false;

% receiver angles
angle1 = 0;
angle2 = 2*pi - 0.0174533; % This ensures that we don't double count the 360 degrees

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

% ------------------------------------------------------------
% Solver params
% ------------------------------------------------------------
params = build_params( ...
    'k', w_0 + freq_gauss_width, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 160, ...
    'n', 20, ...
    'p', 3, ...
    'p_edge', 2);

% ------------------------------------------------------------
% Build all geometries
% ------------------------------------------------------------
num_presets = 5;

all_geometries = cell(1, 2 * num_presets);
all_geometry_names = cell(1, 2 * num_presets);

idx = 1;

% twinjet presets
for i = 1:num_presets
    preset_name = sprintf("plane%d", i);
    gparams = build_twinjet_params('preset', preset_name);

    all_geometries{idx} = twinjet_plane_geometry(gparams);
    all_geometry_names{idx} = sprintf("twinjet_%s", preset_name);

    idx = idx + 1;
end

%quadjet presets
% for i = 1:num_presets
%     preset_name = sprintf("plane%d", i);
%     gparams = build_quadjet_params('preset', preset_name);
% 
%     all_geometries{idx} = quadjet_plane_geometry(gparams);
%     all_geometry_names{idx} = sprintf("quadjet_%s", preset_name);
% 
%     idx = idx + 1;
% end

% ------------------------------------------------------------
% Run experiment
% ------------------------------------------------------------
start = tic;

geom = all_geometries{1};

fprintf("\n==================================================\n");
fprintf("Running object %d / %d: %s\n", 1, ...
    length(all_geometries), all_geometry_names{1});
fprintf("==================================================\n");

% Build a refined curve.
C = build_curve(geom, params);
C = refine_curve_wavenumber(C, 2*pi/params.k, 1);

for jj = 1:2
    C = bisect_all_patches(C);
    C = build_global_indexing(C);
    C = get_close_points(C);
end

lp = RP2LP(params, geom);
lp = lp.update_RP_Curve(C);

% total number of points in the curve
N = lp.RP_Curve.global.N

% rule of thumb for computing complex resonances
w_len = wlims(2) - wlims(1);
num_sub_ints = ceil(w_len / 2.0);
jump = w_len / num_sub_ints;

pols_wp = [];
pstart = tic;

load("/scratch/msantana/plane_data_360/twinjet_plane1.mat");
pols_wp = pols;
pol_comp_time = toc(pstart)

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
f_sols_all = [];

rind = 15;
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
    'is_far_field', true, ...
    'xs', x, ...
    'ys', y);

tic
[uff, scat_sol, r_res, ~, smoothint, f_sols, zeroint] = cts_wpols_LU(psp, lp, pols_wp, ...
    L_rl, U_rl, L_pole, U_pole);
single_reciever_comp_time = toc

uff = uff(:).';
r_res = r_res(:).';
f_sols = f_sols(:).';
%%
Bkslow = comp_bkslow(psp);

% F3 compute the densitites
densities = comp_dens(psp,lp,L_rl, U_rl);

% F4-5 Compute The frequency domain solutions.
[f_sols_b] = comp_freq_sols(psp,lp, densities, Bkslow);

% T0-T5 Compute Time Solution
% We pass in empty for the last two so the parfor loop works.
[usol_b, scat_sol_b] = comp_time_sol(psp,lp,f_sols);
return
%% Todo: Check if the pole subtraction added a spike. Check if SEM matches the origonal signal.

figure(1)
clf
title("Frequency domain without subtraction")
plot(psp.ws, real(f_sols_b))

figure(2)
clf
title("Frequency domain with subtraction")
plot(psp.ws, real(f_sols))

%%
figure(3)
clf
title("Time domain without subtraction")
plot(psp.ts,real(usol_b))

figure(4)
clf
title("Time domain with subtraction")
plot(psp.ts,uff)





