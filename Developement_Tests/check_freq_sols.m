% Script to check that singularites really get subtracted out
% load("/scratch/msantana/plane_data_360/twinjet_plane3.mat")
load("/scratch/msantana/plane_data_360/twinjet_plane3.mat")
% load("/scratch/msantana/plane_data_360/quadjet_plane3.mat")

%%
rec_ind = 1;
figure(1)
clf
plot(psp.ws, real(f_sols_all(rec_ind,:)))

max(abs(r_res_all))

figure(2)
clf
plot(ts, real(uff_all(rec_ind,:)))
hold on

 ps2 = psp;
ps2.ws = pols; ps2.numw = length(pols);

bk_slow_pols = comp_bkslow(ps2,true);
% Multiply by the recentering term to compute the correct residue
bk_slow_pols = bk_slow_pols;
% We get the residue we need the contribution from B_k as it limits to
% the pole, and then add all the windows the windows together.
% Make it a column vector
bk_slow_pols = sum(bk_slow_pols,1);

plot(ts, real(-1i * bk_slow_pols * r_res_all(rec_ind) * exp(-1i *pols * (psp.ts - psp.t_0))))
hold off



%% This configuration doesn't seem to work.
% Script to check that singularites are real
% load("/scratch/msantana/plane_data_360/twinjet_plane3.mat")
% 
% %%
% figure(1)
% clf
% plot(real(f_sols_all(2,:)))
% 
% max(abs(r_res_all))
