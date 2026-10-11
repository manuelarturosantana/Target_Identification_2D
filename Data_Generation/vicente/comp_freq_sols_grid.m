function f_sols = comp_freq_sols_grid(ps, lp, densities, Bkslow, pols, r_res)
%COMP_FREQ_SOLS_GRID  Vectorized analogue of TD_Solver_Code/comp_freq_sols.m.
%
% Same reasoning as comp_r_res_grid.m: replace the per-point call to the
% near-field evaluator with one vectorized call across all requested
% points per frequency. Math is identical to comp_freq_sols.m.
%
% Usage mirrors comp_freq_sols.m:
%   f_sols = comp_freq_sols_grid(ps,lp,densities,Bkslow)
%   f_sols = comp_freq_sols_grid(ps,lp,densities,Bkslow,pols,r_res)

npts = ps.num_spat_pts;
has_pols_res = (nargin >= 5);

f_sols = complex(zeros(npts, ps.numw, ps.numk));

for wind = 1:ps.numw
    paramsW = lp.params;
    paramsW.k = ps.ws(wind);

    solw = evaluate_helmholtz_potential( ...
        densities(:, wind), ps.xs, ps.ys, lp.RP_Curve, paramsW);

    if has_pols_res
        for pind = 1:length(pols)
            solw = solw-r_res(:, pind)./(ps.ws(wind)-pols(pind));
        end
    end

    for kind = 1:ps.numk
        f_sols(:, wind, kind) = Bkslow(kind, wind)*solw;
    end
end

end
