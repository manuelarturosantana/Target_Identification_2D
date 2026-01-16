function [usol, scat_sol, r_res, pols, smoothint, f_sols, zeroint] = cts_wpols_LU(ps,lp, pols, L_rl, U_rl,L_pole,U_pole)
    % Version of the algorithm where the poles are computed using the adaptive search
    % method. Then residues are subtracted out in the inverse fourier transform.
    % Inputs:
    %   pols must be passed in.
    %   L_rl, U_rl is a cell array of the needed real frequency LUs
    %   L_pole, U_pole are cell arrays of cell arrays of the needed LUs for
    %   computing the density residue.

    disp("Starting with Res dens");
    res_dens = comp_res_dens(ps, lp, pols, L_pole, U_pole);
    
    % F2 compute bkslow
    disp("Starting Bkslow");
    Bkslow = comp_bkslow(ps);
   
    % F3 compute the densitites
    disp("Starting Densitites");
    densities = comp_dens(ps,lp, L_rl, U_rl);
   

    % Compute the residues
    disp("Starting r res");
    r_res = comp_r_res(ps,lp,pols,res_dens);
  
    
    % F4-5 Compute The frequency domain solutions.
    disp("starting freq sols")
    f_sols = comp_freq_sols(ps,lp, densities, Bkslow, pols, r_res);
    
    % T0-T5 Compute Time Solution
    disp("starting time sols")
    [usol, scat_sol,smoothint,zeroint] = comp_time_sol(ps,lp,f_sols,pols,r_res);
end