function cmap = diverging_colormap(n)
%DIVERGING_COLORMAP  Self-contained blue-white-red colormap for field plots.
% Avoids depending on toolbox-specific maps (e.g. redbluecmap).
%
% Local copy of rp2d_matlab/experiments/visual/diverging_colormap.m --
% duplicated here rather than pulled in from rp2d_matlab's experiments/
% folder (scratch/exploratory code, not part of its stable src/ API).

if nargin < 1
    n = 256;
end

anchors = [
    0.019, 0.188, 0.380   % dark blue
    0.259, 0.573, 0.776   % blue
    0.969, 0.969, 0.969   % near white
    0.839, 0.376, 0.302   % red
    0.404, 0.000, 0.121   % dark red
];

x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi);

end
