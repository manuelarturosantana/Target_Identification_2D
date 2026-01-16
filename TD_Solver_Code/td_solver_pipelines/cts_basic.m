function [usol, scat_sol, f_sols] = cts_basic(ps, lp)
    % Implementation of the basic algorithm from Anderson-Bruno-Lyons

    % F2
    Bkslow = comp_bkslow(ps);

    % F3 compute the densitites
    densities = comp_dens(ps,lp);
    
    % F4-5 Compute The frequency domain solutions.
    [f_sols] = comp_freq_sols(ps,lp, densities, Bkslow);
    
    % T0-T5 Compute Time Solution
    % We pass in empty for the last two so the parfor loop works.
    [usol, scat_sol] = comp_time_sol(ps,lp,f_sols);
end