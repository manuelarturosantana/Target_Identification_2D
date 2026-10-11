% replot_glrtvote_logscale.m
%
% Re-renders a saved dominant_eigenfunction_glrtvote_*.mat result with a
% signed-log (asinh) color transform instead of a linear/clipped scale,
% to fix the "field too strong inside the nacelle, leaking field looks
% too faint" problem: a hard linear cap can't help when the nacelle
% peak is orders of magnitude above the exterior leaked field, since
% even a small fraction of that peak still swamps the color scale.
%
%   transformed = asinh(Re(Efun) / eps_scale)
%
% asinh is linear near zero and behaves like a signed log far from zero,
% so it compresses the huge dynamic range between the resonant nacelle
% interior and the faint far-field leakage while preserving sign (and
% therefore the oscillatory/nodal structure of the eigenfunction) --
% unlike a magnitude-only dB scale, which would erase the sign changes.
%
% eps_scale is set from the median |field| magnitude (excluding exact
% zeros/NaNs), a robust "typical" scale unaffected by the nacelle's
% extreme outliers, so that most of the exterior field sits in the
% quasi-linear part of asinh and only the nacelle peak gets compressed.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

mat_path = "/scratch/vhojas/data_planes_2d/glrtvote_figures/dominant_eigenfunction_glrtvote_quadjet_plane1_ng1024_resolved.mat";
save_dir = "/scratch/vhojas/data_planes_2d/glrtvote_figures";

S = load(mat_path, "xv", "yv", "Efun", "C", "Solid", "target_pole", ...
    "Q_target", "resolved", "plane_name", "ng", "max_votes", "nrec", "true_max");

field = real(S.Efun);
finite_vals = abs(field(~isnan(field)));
finite_vals = finite_vals(finite_vals > 0);
eps_scale = median(finite_vals);

fprintf("true max = %.4e, median |field| (eps_scale) = %.4e, ratio = %.2f\n", ...
    S.true_max, eps_scale, S.true_max / eps_scale);

transformed = asinh(field / eps_scale);
cmax_t = max(abs(transformed(~isnan(transformed))));

fig = figure('Color', 'w', 'Visible', 'off');
imagesc(S.xv, S.yv, transformed);
set(gca, 'YDir', 'normal');
hold on;

has_solid = ~isempty(S.Solid);
if has_solid
    plot(S.Solid, 'FaceColor', [0.78 0.78 0.78], 'FaceAlpha', 1, ...
        'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 1);
end
for j = 1:numel(S.C.patches)
    Patch = S.C.patches{j};
    t = linspace(-1, 1, 100);
    plot(Patch.x(t), Patch.y(t), 'k-', 'LineWidth', 1);
end

axis equal tight;
colormap(parula(256));
if cmax_t > 0
    clim([-cmax_t, cmax_t]);
end
cb = colorbar;
cb.Label.String = sprintf('asinh(Re(eigenfunction) / %.2e)', eps_scale);
box on;
title(sprintf('%s eigenfunction (GLRT vote, signed-log scale), pole k=%.4f%+.4fi (Q~%.0f) [%s, ng=%d, %d/%d votes]', ...
    strrep(S.plane_name, '_', '\_'), real(S.target_pole), imag(S.target_pole), ...
    S.Q_target, ternary(S.resolved, "RESOLVED", "UNRESOLVED"), S.ng, S.max_votes, S.nrec));
hold off;

fname_base = sprintf("dominant_eigenfunction_glrtvote_%s_ng%d_%s_logscale", S.plane_name, S.ng, ...
    ternary(S.resolved, "resolved", "unresolved"));
img_path = fullfile(save_dir, fname_base+".png");
exportgraphics(fig, img_path, 'Resolution', 300);
fprintf("Saved: %s\n", img_path);

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
