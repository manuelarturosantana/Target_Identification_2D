% replot_dominant_eigenfunction.m
%
% Re-renders the dominant-eigenfunction PNGs from the already-saved
% figures/dominant_eigenfunction_*.mat files, without redoing the
% (expensive) field evaluation. Only the color-scale cap changes: instead
% of clipping at the clip_pct-th percentile of |field|, this clips at
% cap_frac * true_max.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

save_dir = fullfile(fileparts(mfilename('fullpath')), "figures");
cap_frac = 0.2;

mat_files = dir(fullfile(save_dir, "dominant_eigenfunction_*.mat"));
fprintf("Found %d saved result(s) in %s\n", numel(mat_files), save_dir);

for i = 1:numel(mat_files)
    mat_path = fullfile(mat_files(i).folder, mat_files(i).name);
    fprintf("\n[%d/%d] %s\n", i, numel(mat_files), mat_files(i).name);

    S = load(mat_path, "xv", "yv", "Efun", "C", "Solid", "target_pole", ...
        "Q_target", "resolved", "plane_name", "ng", "true_max");

    cmax = cap_frac * S.true_max;
    fprintf("  true max = %.4e, capping color scale at %.2f*max = %.4e\n", ...
        S.true_max, cap_frac, cmax);

    fig = figure('Color', 'w', 'Visible', 'off');
    imagesc(S.xv, S.yv, real(S.Efun));
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
    if cmax > 0
        clim([-cmax, cmax]);
    end
    cb = colorbar;
    cb.Label.String = 'Re(eigenfunction)';
    box on;
    title(sprintf('%s eigenfunction, pole k=%.4f%+.4fi (Q~%.0f) [%s, ng=%d, cap=%.0f%% max]', ...
        strrep(S.plane_name, '_', '\_'), real(S.target_pole), imag(S.target_pole), ...
        S.Q_target, ternary(S.resolved, "RESOLVED", "UNRESOLVED"), S.ng, 100*cap_frac));
    hold off;

    fname_base = sprintf("dominant_eigenfunction_%s_ng%d_%s", S.plane_name, S.ng, ...
        ternary(S.resolved, "resolved", "unresolved"));
    img_path = fullfile(save_dir, fname_base+".png");
    exportgraphics(fig, img_path, 'Resolution', 300);
    close(fig);
    fprintf("  Saved: %s\n", img_path);
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
