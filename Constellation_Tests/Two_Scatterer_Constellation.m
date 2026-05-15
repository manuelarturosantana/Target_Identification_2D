box_lims = [38, 40, -0.2, 0];

N = ceil((box_lims(2) * 12) / 2);

rp = rec_params("xlims",box_lims(1:2),"ylims",box_lims(3:4),"ntb",150,"nlr",150,"m_rec",7);
figure(1)
clf

%% Compute Single Rocket Evals
% curve_names = {"missile cavity"};
% curve_names = {"elliptic cavity"};
% aps = [0.75];
% xshift = [0];
% yshift = [0];
% draw_constellation = true;
% poles_rocket = compute_constellation_evals2(curve_names,  aps, xshift,yshift, N, ...
%     rp, draw_constellation);

%% Compute the Circle Evals
curve_names = {"circular cavity"};
aps = [2.5];
xshift = [0];
yshift = [0];
draw_constellation = true;
[poles_circ, lp_sing] = compute_constellation_evals2(curve_names,  aps, xshift,yshift, N, ...
    rp, draw_constellation);

%% Now together
curve_names = {"circular cavity", "circular cavity"};
% curve_names = {"circular cavity", "elliptic cavity"};
aps = [2.5, 2.5];
xshift = [0, 2.5];
yshift = [0, 0];
draw_constellation = true;
[poles_together, lp, mc] = compute_constellation_evals2(curve_names,  aps, xshift,yshift, N, ...
    rp, draw_constellation);

%%
figure(2)
clf
tiledlayout(1,2)
nexttile
mc.draw()
axis equal
legend("Scatterer 1", "Scatterer 2")
nexttile
hold on
plot(poles_circ,'o','markersize',15)
plot(poles_together,'kp','markersize',15)
hold off
legend("Scatterer 1 Poles", "Constellation Poles")
ax = gca; ax.set("Fontsize",20);
axis equal

%% Demonstrate the pole splitting
poles_circ = sort(poles_circ,'ComparisonMethod','real');
poles_together = sort(poles_together,'ComparisonMethod','real');
poles_circ(end)
pc1 = poles_together(end-6)
pc2 = poles_together(end-7)
[~,I] = min(abs(poles_circ - pc1));
psing = poles_circ(I)

%% Sanity check that the two poles on the end really are different
S = svd(lp.bie_mat(pc1)); S(end)
S = svd(lp.bie_mat(pc2)); S(end)
S = svd(lp.bie_mat((pc1 + pc2) / 2)); S(end)
abs(pc1 + pc2)

%% Compute the two eigen-densities
A1 = lp.bie_mat(pc1);
[~, ~, V] = svd(A1, 'econ');
dens1 = V(:, end); norm(A1 * dens1)

A2 = lp.bie_mat(pc2);
dens2 = comp_ns(A2); norm(A2 * dens2)

A3 = lp_sing.bie_mat(psing);
dens_sing = comp_ns(A3); norm(A3 * dens_sing)

%% Plot the two eigenfunctions
xlims = [-2, -2 + 2+ 2.5  + 2]; ylims = [-2,2];
[vals1, ~, ~] = scat_field_data(lp, mc, dens1, xlims, ylims, 40, 30);
[vals2, x, y] = scat_field_data(lp, mc, dens2, xlims, ylims, 40, 30);

[vals_sing,xs,ys] = scat_field_data(lp_sing, mc, dens_sing, [-2,2], [-2,2], 40, 30);
%%

figure(1)
clf
hold on
title(sprintf("Eigenfunction corresponding to pole at %.16g + %.16g", real(pc1), imag(pc1)))
[X,Y] = meshgrid(x,y);
pcolor(X,Y,abs(vals1.'))
shading flat
plot(mc.obs{1}.X,mc.obs{1}.Y,'k','linewidth',2);
plot(mc.obs{2}.X,mc.obs{2}.Y,'k','linewidth',2);
hold off


figure(2)
clf
hold on
title(sprintf("Eigenfunction corresponding to pole at %.16g + %.16g", real(pc2), imag(pc2)))
[X,Y] = meshgrid(x,y);
pcolor(X,Y,abs(vals2.'))
shading flat
plot(mc.obs{1}.X,mc.obs{1}.Y,'k','linewidth',2);
plot(mc.obs{2}.X,mc.obs{2}.Y,'k','linewidth',2);
hold off

%%
figure(3)
clf
hold on
title(sprintf("Eigenfunction corresponding to pole at %.16g + %.16g", real(psing), imag(psing)))
[X,Y] = meshgrid(xs,ys);
pcolor(X,Y,abs(vals_sing.'))
shading flat
plot(mc.obs{1}.X,mc.obs{1}.Y,'k','linewidth',2);
axis equal
hold off




%% Versions which allow one to do constellations with different curves in them
function [poles, lp, mc] = compute_constellation_evals2(...
    curve_names, aps, xshift, yshift, N,...
    rp, draw_constellation)
    % Function which builds a constellation and computes its poles.
    % 
%
%   WARNING: This assumes that each curve is okay having the same number of
%   unknowns
%   
    % Inputs:
    %   curve  : The base curve to build the constellation from
    %   xshift : For each curve shift in x for the
    %            location of the new curve.
    %   yshift : Same as xshift but for y
    %   rp     : Parameters for the aaa recursive calculation
    %   draw_constellation : If true plot the constellation before
    %   computing
    %
    %  Outputs:
    %       poles: The computed poles
    %       lp   : The layer potential object
    %       mc   : The multiple curves object

    [lp, mc] = setup_constellation2(curve_names, aps,  xshift, yshift, N, draw_constellation);

    N_total = length(xshift) * N;


    v = randn(N_total,1); u = rand(1,N_total);

    f = @(k) u * (lp.bie_mat(k) \ v);

    poles = aaa_recursive(f,rp);

end

function [lp, mc] = setup_constellation2(curve_names, aps, xshift, yshift, N, draw_constellation)
    % Function which builds a constellation and computes its poles.
    % 
    % Inputs:
    %   curve_names : Names of the base curves to build the constellation from
    %   aps    : The aperature sizes for the base curves
    %   xshift : For each curve beyond the first, shift in x for the
    %            location of the new curve.
    %   yshift : Same as xshift but for y
    %   All of the above need to be of the same length for this to work as
    %   expected.
    %  draw_constellation : If true plot the constellation 

    all_curves = {};
 
    for ii = 1:length(xshift)
        new_curve = Build_Curve(curve_names{ii},N,aps(ii));
        new_curve.X = new_curve.X  + xshift(ii);
        new_curve.Y = new_curve.Y  + yshift(ii);
        all_curves{end+1} = new_curve;
    end
    

    mc = MultipleCurves(false, all_curves{:});

    if draw_constellation
        figure; mc.draw(); axis equal; title('Constellation Geometry');
    end
    
    lp = OpenMultipleScattering(mc, false);
  

end

