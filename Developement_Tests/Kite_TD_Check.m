% Simple sanity check for time domain scattering by the kite to make sure
% things look okay

curve = Kite(200,[],true);

figure(1)
clf
hold on
curve.draw()
% NF
% xs = [-1.2,1.2]; ys = [0,0];

% FF
xs = [-1,1]; ys = [0,0];

plot(xs,ys,'*')
hold off

[sigmas,wlims] = gauss_picker(5,2,1e-10);

lp = CombinedPotential(curve,-4);

ps = problem_data("xs",xs,"ys",ys,"wlims",wlims,"sigmas",sigmas,'Tlims',[0,25],'is_open_curve',false, ...
    'kappa',[-1,0],'numt',300,'numw',200,'w_0',5,'t_0',10,'is_far_field',true);

[usol, scat_sol, f_sols] = cts_basic(ps,lp);

%%
figure(2)
clf
tiledlayout(1,2)
nexttile
plot(ps.ts,usol(1,:))
nexttile
plot(ps.ts,usol(2,:))

%%
ps = problem_data("xlims",[-3,3],"ylims",[-3,3],"numx",200,"numy",200,"wlims",wlims,"sigmas",sigmas,'Tlims',[0,25],'is_open_curve',false, ...
    'kappa',[1,1],'numt',150,'numw',200,'w_0',5,'t_0',10);

[usol, scat_sol, f_sols] = cts_basic(ps,lp);

%% The corner is off in the plot by this reshape
to_plot = reshape(usol, 200, 200, 150);
% to_plot = reshape(scat_sol,200,200,150);
to_plot_max = max(abs(to_plot),[],'all');
tempxs  = linspace(ps.xlims(1),ps.xlims(2),ps.numx);
tempys  = linspace(ps.ylims(1),ps.ylims(2),ps.numy);

[X,Y] = meshgrid(tempxs,tempys); 
X = X'; Y = Y.';

for ii = 1:20
    figure(ii)
    clf
    % imagesc(real(to_plot(:,:,ii*5)))
    surf(X,Y, real(to_plot(:,:,ii*5)))
    % set(gca,'YDir','normal')
    colorbar
    caxis([0, 0.25])
    hold on
 
end

%%
figure(1)
clf
tempxs  = linspace(ps.xlims(1),ps.xlims(2),ps.numx);
tempys  = linspace(ps.ylims(1),ps.ylims(2),ps.numy);

[X,Y] = meshgrid(tempxs,tempys); 
X = X'; Y = Y.';

surf(X,Y, real(to_plot(:,:,1*5)))



