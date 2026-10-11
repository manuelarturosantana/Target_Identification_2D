% diag_quadjet_plane5_visual.m
%
% Visual proof for the quadjet_plane5 wing/fuselage junction gap:
%   (1) a paper-style full-plane outline plot (same method as
%       Paper_Figure_Generator/plot_all_planes.m: raw per-patch line
%       plot, no fill) -- to show why the gap isn't obvious there.
%   (2) a zoomed inset at the junction point with the two curve
%       endpoints marked separately, to show the actual gap.
% Uses the exact same build_quadjet_params/quadjet_plane_geometry calls
% as data_gen_quadjet.m and our eigenfunction pipeline.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));

save_dir = fullfile(fileparts(mfilename('fullpath')), "figures");
if ~exist(save_dir, "dir")
    mkdir(save_dir);
end

gparams = build_quadjet_params('preset', 'plane5');
geom = quadjet_plane_geometry(gparams);

front2 = geom.pieces{2};
lwing_p1 = geom.pieces{4};
xfl = front2.x(front2.tf);
yfl = front2.y(front2.tf);
tstar = lwing_p1.t0;
xw = lwing_p1.x(tstar);
yw = lwing_p1.y(tstar);

fprintf("Fuselage endpoint (xfl,yfl)   = (%.6f, %.6f)\n", xfl, yfl);
fprintf("Wing-lead endpoint (xw,yw)    = (%.6f, %.6f)\n", xw, yw);
fprintf("Gap distance                  = %.6e\n", hypot(xw-xfl, yw-yfl));

fig = figure('Color', 'w', 'Position', [100 100 1400 650]);

% --- Panel 1: paper-style full view ---
subplot(1,2,1);
hold on;
for j = 1:numel(geom.pieces)
    P = geom.pieces{j};
    if isempty(P), continue; end
    t = linspace(-1, 1, 400);
    % pieces are parametrized on [t0,tf], not necessarily [-1,1] --
    % match plot_all_planes.m's *intent* (dense uniform sampling over
    % the piece's own domain) rather than its literal (buggy for
    % non-normalized domains) t range.
    t = linspace(P.t0, P.tf, 400);
    plot(P.x(t), P.y(t), 'k-', 'LineWidth', 1.5);
end
plot(xfl, yfl, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
axis equal; grid on;
title('quadjet\_plane5 -- full view (paper-style outline)');

% --- Panel 2: zoomed inset at the junction ---
subplot(1,2,2);
hold on;
t2 = linspace(front2.t0, front2.tf, 2000);
plot(front2.x(t2), front2.y(t2), 'b-', 'LineWidth', 2, 'DisplayName', 'front2 (fuselage)');
tw = linspace(lwing_p1.t0, lwing_p1.tf, 2000);
plot(lwing_p1.x(tw), lwing_p1.y(tw), 'r-', 'LineWidth', 2, 'DisplayName', 'lwing\_p1 (wing lead)');
plot(xfl, yfl, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'fuselage endpoint');
plot(xw, yw, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'wing-lead endpoint');
axis equal; grid on;
legend('Location', 'best');
pad = 0.15;
xlim([min(xfl,xw)-pad, max(xfl,xw)+pad]);
ylim([min(yfl,yw)-pad, max(yfl,yw)+pad]);
title(sprintf('Zoomed junction -- gap = %.4e', hypot(xw-xfl, yw-yfl)));

img_path = fullfile(save_dir, "diag_quadjet_plane5_junction_visual.png");
exportgraphics(fig, img_path, 'Resolution', 200);
fprintf("Saved: %s\n", img_path);
