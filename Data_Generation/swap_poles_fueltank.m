
% Short script which swaps the poles from the twinjet to the twinjet with
% fuel tank. This allows the classification algorithm to use the computed
% twinjet poles for to classify the version with the fuel tank
clear
load('/scratch/msantana/plane_data_360/twinjet_plane1.mat')

pols2save = pols;

load('/scratch/msantana/plane_data_fuel_tank/twinjet1_with_fuel_tank.mat')

abs(pols2save - pols) / abs(pols)
pols = pols2save;
clear pols2save

save('/scratch/msantana/plane_data_fuel_tank/twinjet1_with_fuel_tank.mat')