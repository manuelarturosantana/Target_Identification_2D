function r_res = comp_r_res_grid(ps, lp, pols, res_dens)
%COMP_R_RES_GRID  Vectorized analogue of TD_Solver_Code/comp_r_res.m.
%
% comp_r_res.m loops over every spatial point and calls the near-field
% evaluator once per point per contour sample, which is fine for a
% handful of far-field receivers but is far too slow for an image grid
% (thousands of points x n_res_w contour samples x num_pols). Since
% rp2d_matlab's evaluate_helmholtz_potential/single_layer_extern already
% accept whole point arrays (that's how plot_scattering_field.m uses
% them), this instead makes one vectorized call across every requested
% point per (pole, contour-point) pair. The residue formula itself is
% identical to comp_r_res.m -- only the evaluation strategy changes.
%
% Assumes ps.xs/ps.ys are already restricted to exterior points (the
% inside-the-solid test that comp_r_res.m does per point via ps.is_rp_curve
% is skipped here for that reason).
%
% Inputs / outputs: see comp_r_res.m.

npts = ps.num_spat_pts;
r_res = complex(zeros(npts, length(pols)));

eit = exp(1i*2*pi*(1:ps.n_res_w)/ps.n_res_w);
dt  = 2*pi/ps.n_res_w;
w = ps.cont_rad*eit(:);   % contour quadrature weights, n_res_w x 1

for pind = 1:length(pols)
    circ = pols(pind)+ps.cont_rad*exp(1i*2*pi*(1:ps.n_res_w)/ps.n_res_w);

    fvals_k = complex(zeros(npts, ps.n_res_w));
    for cind = 1:ps.n_res_w
        paramsC = lp.params;
        paramsC.k = circ(cind);
        fvals_k(:, cind) = evaluate_helmholtz_potential( ...
            res_dens(:, cind, pind), ps.xs, ps.ys, lp.RP_Curve, paramsC);
    end

    f2vals = 1./(circ-pols(pind));
    denom = 1i*dt*sum(f2vals(:).*w);
    numer = (1i*dt)*(fvals_k*w);

    r_res(:, pind) = numer/denom;
end

end
