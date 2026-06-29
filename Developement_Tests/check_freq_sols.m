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
plot(real(uff_all(rec_ind,:)))



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
