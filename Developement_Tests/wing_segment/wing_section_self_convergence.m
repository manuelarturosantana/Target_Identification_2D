% test_plane_scatter_self_convergence_field.m
clear; clc;
addpath(genpath('/home/msantana/Progamming/rp2d_matlab/src'))
addpath(genpath('/home/msantana/Progamming/Target_Identification_2D'))
addpath(genpath('/home/msantana/Progamming/MATLAB_PACKAGES'))

rng(1);

fprintf("==============================================\n");
fprintf(" Self-convergence test (field at target pts)\n");
fprintf("==============================================\n\n");

% --------------------------------------------------
% parameters
% --------------------------------------------------
params = build_params();
params.n = 12;
params.p_edge = 2;
params.p = 3;
params.delta = 0.05;
% params.k = 2*pi;
params.k = 50;
params.npolar = 160;
params.solver = 'direct';

xhat = [1/sqrt(2); 1/sqrt(2)];

% --------------------------------------------------
% refinement levels
% --------------------------------------------------
NREF = 4;   % comparisons; finest is reference

% --------------------------------------------------
% base geometry
% --------------------------------------------------
gparams = build_twinjet_params();
geom = wing_section_geometry(gparams);

C0 = build_curve(geom, params);
C0 = refine_curve_wavenumber(C0, 2*pi/params.k, 2);

Curves = cell(NREF + 1, 1);
Phis   = cell(NREF + 1, 1);

Curves{1} = C0;

fprintf("Building nested curves + solving...\n");
for lev = 1:(NREF + 1)
    C = Curves{lev};

    % solve on this level
    phi = solve_planewave_scattering(C, xhat, params);

    Curves{lev} = C;
    Phis{lev} = phi;

    fprintf("  Level %d: patches = %4d, N = %6d\n", ...
        lev - 1, numel(C.patches), C.global.N);

    % refine for next level (ensure bookkeeping after bisection)
    if lev <= NREF
        Cn = bisect_all_patches(C);
        Cn = build_global_indexing(Cn);
        Cn = get_close_points(Cn);
        Curves{lev + 1} = Cn;
    end
end
%%
% --------------------------------------------------
% choose fixed exterior target points
% constraints:
%   1) r = sqrt(x^2+y^2) >= 1
%   2) distance from boundary nodes >= c*delta (avoid near-singular eval)
% --------------------------------------------------
nt_want = 12;
L = 1.5;
ncand = 600;

Cref = Curves{NREF + 1};
phiref = Phis{NREF + 1};

[xb, yb] = get_global_xy(Cref);

xc = -L + (2*L) * rand(ncand, 1);
yc = -L + (2*L) * rand(ncand, 1);

rc = sqrt(xc.^2 + yc.^2);
keep = (rc >= 1.0);
xc = xc(keep);
yc = yc(keep);

c = 6;
keep = true(numel(xc), 1);
for j = 1:numel(xc)
    dx = xc(j) - xb;
    dy = yc(j) - yb;
    if min(dx.^2 + dy.^2) < (c * params.delta)^2
        keep(j) = false;
    end
end

xc = xc(keep);
yc = yc(keep);

if numel(xc) < nt_want
    error("Not enough valid target points. Increase ncand or L.");
end

xt = xc(1:nt_want);
yt = yc(1:nt_want);
nt = nt_want;

fprintf("\nUsing %d target points with r>=1 and dist>=%.3f.\n", ...
    nt, c * params.delta);

%% --------------------------------------------------
% reference scattered field at targets
% --------------------------------------------------
us_ref = zeros(nt, 1);
for j = 1:nt
    us_ref(j) = single_layer_extern( ...
        phiref, xt(j), yt(j), Cref, params.k);
end

% --------------------------------------------------
% compare each level to reference via target fields
% --------------------------------------------------
fprintf("\nComparing to reference level %d\n", NREF);
fprintf("------------------------------------------------------------------\n");
fprintf("lev   N        time(s)    maxdiff(ref)   log10(err)    rate\n");
fprintf("------------------------------------------------------------------\n");

errs = zeros(NREF, 1);
Ns   = zeros(NREF, 1);
err_prev = NaN;

for lev = 1:NREF
    Cc = Curves{lev};
    phic = Phis{lev};
    Ns(lev) = Cc.global.N;

    tic;
    us_c = zeros(nt, 1);
    for j = 1:nt
        us_c(j) = single_layer_extern( ...
            phic, xt(j), yt(j), Cc, params.k);
    end

    d = abs(us_c - us_ref) ./ us_ref;
    dmax = max(d);
    t_elapsed = toc;

    errs(lev) = dmax;

    if lev == 1
        rate_str = "   -";
    else
        rate = log(err_prev / dmax) / log(2);
        rate_str = sprintf("%6.2f", rate);
    end
    err_prev = dmax;

    fprintf("%3d  %6d   %8.3f   % .3e     %8.3f    %s\n", ...
        lev - 1, Ns(lev), t_elapsed, dmax, log10(dmax), rate_str);
end

fprintf("------------------------------------------------------------------\n");
fprintf("Done.\n\n");
%%
% --------------------------------------------------
% optional plot
% --------------------------------------------------
% figure;
% loglog(Ns, errs, 'o-');
% grid on;
% xlabel('N (coarse dofs)');
% ylabel('max_j |u_scat^N(x_j) - u_scat^{ref}(x_j)|');
% title('Self-convergence via scattered field at fixed targets');
% 
% % --------------------------------------------------
% % optional: visualize target points
% % --------------------------------------------------
% figure;
% plot_curve(Cref);
% hold on;
% plot(xt, yt, 'ro', 'MarkerFaceColor', 'r');
% axis equal;
% title('Reference curve and target points');
