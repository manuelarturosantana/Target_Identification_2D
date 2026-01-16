function [usol, scat_sol, r_res, pols,smoothint, f_sols, zeroint, num_inversions] = cts_denssv(ps,lp)
    % Version of the algorithm where the poles are computed using set valued aaa on the
    
    % F2 compute bkslow
    disp("Staring bk slow");
    Bkslow = comp_bkslow(ps);
    
    % Compute the rational approximant and the densities.
    disp("Starting recursive density/pol computation");
    [densities, Rs, pol_c,num_inversions] = ra_dense(ps, lp);
    
    disp("Starting res_dense_sv")
    [pols, res_dens] = res_dens_sv(ps, Rs, pol_c);
    
    % Compute the residues
    disp("starting residue computation");
    r_res = comp_r_res(ps,lp,pols,res_dens);
    
    % F4-5 Compute The frequency domain solutions.
    disp("starting freq sols computation")
    f_sols = comp_freq_sols(ps,lp, densities, Bkslow, pols, r_res);
    
    % T0-T5 Compute Time Solution
    disp("starting time sol computation")
    [usol, scat_sol, smoothint, zeroint] = comp_time_sol(ps,lp,f_sols,pols,r_res);
end
