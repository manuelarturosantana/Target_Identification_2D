load("poles_convergence.mat")
poles = Poles{end};

params = build_params();
params.n = 12;
params.p_edge = 2;
params.p = 6;
params.delta = 0.05;
params.k = 65;
params.npolar = 80;
params.solver = 'direct';

% --------------------------------------------------
% User-specified frequency range for pole finding
% --------------------------------------------------
wlims = [63, 65];   % <-- edit these bounds
wimag = -0.2;
numw  = 400;

% --------------------------------------------------
% refinement levels
% --------------------------------------------------
NREF = 6;   % comparisons; finest is reference

% --------------------------------------------------
% base geometry
% --------------------------------------------------
gparams = build_twinjet_params();
geom = twinjet_plane_geometry(gparams);
Poles  = {};
all_svd_vals = {};

fprintf("Building nested curves + finding poles...\n");
for lev = 1:(NREF)

    % build LP object (needed for comp_poles)
    lp = RP2LP(params, geom, lev);
    N = lp.RP_Curve.global.N;

    % subdivide frequency range into sub-intervals (same rule as gen script)
    w_len       = wlims(2) - wlims(1);
    num_sub_ints = ceil(w_len / 2.0);
    jump        = w_len / num_sub_ints;

   
    svd_vals = [];
    for ii = 1:length(poles)
        S = svd(lp.bie_mat(poles(ii)));
        svd_vals(end+1) = S(end)
    end

    all_svd_vals{end+1} = svd_vals;

end