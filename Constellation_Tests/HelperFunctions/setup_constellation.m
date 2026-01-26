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
