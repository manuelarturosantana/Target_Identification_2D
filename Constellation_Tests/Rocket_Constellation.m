
box_lims = [40, 42, -0.2, 0];

N = ceil((box_lims(2) * 12) / 2);

rp = rec_params("xlims",box_lims(1:2),"ylims",box_lims(3:4),"ntb",150,"nlr",150,"m_rec",7);
figure(1)
clf

%% Compute_Single Circle Eigenvalues
curve = Build_Curve("missile cavity", N, 0.75);

xshift = [];
yshift = [];
draw_constellation = true;
poles1 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);


%%
xshift = [2.5];
yshift = [0];
poles2 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);

%%
% xshift = [2.5, 2.5];
% yshift = [0, 1];
xshift = [2.5, 2.2];
yshift = [0, 1];
[poles3, lp, mc] = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);
return
%%
for ii = 1:length(poles3)
    S = svd(lp.bie_mat(poles3(ii)));
    S(end)
end

%%
markers = {'.', 'o','x','+','s','d','^','<','>','p','h'};
figure(1)
clf
mc.draw()
axis equal

figure(2)
clf
hold on
plot(poles1,'*','Markersize',15)
plot(poles3,'s','Markersize',15)
legend("Single Scatterer Poles","Three Constellation Poles")
ax = gca;
ax.set("FontSize",20);
hold off

%%
figure(3)
clf
hold on
plot(poles1,'*','Markersize',15)
plot(poles3,'s','Markersize',15)
legend("Single Scatterer Poles","Three Constellation Poles")
ax = gca;
ax.set("FontSize",20);

xlim([41.95,42])
ylim([-0.08,-0.079])
hold off


