% Script to compute the eigenvalues of the open circle, and compare to the
% eigenvalues of various constellations

box_lims = [48, 50, -0.2, 0];

N = (box_lims(2) * 12) / 2;

rp = rec_params("xlims",box_lims(1:2),"ylims",box_lims(3:4),"ntb",150,"nlr",150,"m_rec",7);
figure(1)
clf

%% Compute_Single Circle Eigenvalues
curve = Build_Curve("circular cavity", N, 0.75);

xshift = [];
yshift = [];
draw_constellation = true;
poles1 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);


%%
xshift = [2.5];
yshift = [0];
poles2 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);

%%
xshift = [2.5, 5.5];
yshift = [0, 0];
poles3 = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation);
return
%%
% Next verify everything is a actually a true pole via the SVD. 
% Think about what this is going to mean for GLRT etc.
% save the data in an easy to present form for Oscar and Vicente

figure(1)
clf
hold on
plot(poles1,'*')
plot(poles2,'^')
plot(poles3,'p')

save("Poles123.mat","poles3","poles2","poles1")
%% Functions that help
function poles = compute_constellation_evals(curve, xshift,yshift,rp, draw_constellation)
    % Function which builds a constellation and computes its poles.
    % 
    % Inputs:
    %   curve  : The base curve to build the constellation from
    %   xshift : For each curve beyond the first, shift in x for the
    %            location of the new curve.
    %   yshift : Same as xshift but for y
    %   rp     : Parameters for the aaa recursive calculation
    %   draw_constellation : If true plot the constellation before
    %   computing

    [lp, ~] = setup_constellation(curve, xshift, yshift, draw_constellation);


    N_total = (length(xshift) + 1) * curve.N;

    v = randn(N_total,1); u = rand(1,N_total);

    f = @(k) u * (lp.bie_mat(k) \ v);

    poles = aaa_recursive(f,rp);

end

function [lp, mc] = setup_constellation(curve, xshift, yshift, draw_constellation)
    % Function which builds a constellation and computes its poles.
    % 
    % Inputs:
    %   curve  : The base curve to build the constellation from
    %   xshift : For each curve beyond the first, shift in x for the
    %            location of the new curve.
    %   yshift : Same as xshift but for y
    %  draw_constellation : If true plot the constellation 

    all_curves = {};
    all_curves{end+1} = curve;

    for ii = 1:length(xshift)
        new_curve = curve;
        new_curve.X = new_curve.X  + xshift(ii);
        new_curve.Y = new_curve.Y  + yshift(ii);
        all_curves{end+1} = new_curve;
    end
    
    if isempty(xshift)
        lp = OpenSingleLayer(curve);
        mc = [];
        if draw_constellation
            figure; plot(curve.X, curve.Y); drawnow;
        end

    else


        mc = MultipleCurves(false, all_curves{:});
    
        if draw_constellation
            figure; mc.draw(); axis equal; title('Constellation Geometry');
        end
        
        lp = OpenMultipleScattering(mc, true);
    end

end

