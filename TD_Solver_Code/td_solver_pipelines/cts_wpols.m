function [usol, scat_sol, r_res, pols, smoothint, f_sols, zeroint] = cts_wpols(ps,lp, pols)
    % Version of the algorithm where the poles are computed using the adaptive search
    % method. Then residues are subtracted out in the inverse fourier transform.
    % Inputs:
    %   the pols argument is optional, if it is not passed in poles are computed.
    if nargin == 2
        [pols, err] = comp_poles(ps,lp,true);
    end

    disp("Starting with Res dens");
    res_dens = comp_res_dens(ps, lp, pols);
    
    % F2 compute bkslow
    disp("Starting Bkslow");
    Bkslow = comp_bkslow(ps);
   
    % F3 compute the densitites
    disp("Starting Densitites");
    densities = comp_dens(ps,lp);
   

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