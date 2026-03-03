classdef RP2LP 
    % This class is to convert Vicente's rectangular polar API into my
    % layerpotential API so that I don't have to add a bunch of switch
    % statments into the time domain solver code.
    %
    % WARNING: This will only work with the time domain code for far field
    % evaluations and eigenvalue computations. 
    %
    % To do the near field the following changes would need to be made:
    % 1. Since other open curve code uses the Diffarcs Test_Distance
    %    function a corresponding version of that function would have to be
    %    implemented for the RP curves. (or just set to always return 1)
    % 2. For same reason the a version of the Scattered_Field diffarcs
    %    function would have to be changed for the near field evaluations.
    %

    properties
        curve; % struct containing curve.X, curve.Y which give the global 
               % curve values of x and y.
        RP_Curve; % curve data for the rectangular polar
        params;
    end

    methods
        function lp = RP2LP(params, geom, ppwl)
            % Inputs:
            %   params: Parameter struct from Vicentes code. Note the 
            %           k value which is in params should be the highest
            %           frequency k value which you will use this solver to
            %           solve.
            %
            %   geom  : Geometry struct from Vicentes code
            %   ppwl  : (optional) number of patches per wavelength,
            %           default 1.

           
            if (nargin < 3)
                ppwl = 1;
            end
            
            lp.params = params;

            RP_Curve = build_curve(geom, params);
            lp.RP_Curve = refine_curve_wavenumber(RP_Curve, 2*pi/params.k, ppwl);

            [xb, yb] = get_global_xy(lp.RP_Curve);
            curve.X = xb; curve.Y = yb;
            lp.curve = curve;          
        end

        function mat = bie_mat(lp, k)
            lp.params.k = k;
            mat = single_layer_mat(lp.RP_Curve, lp.params);
        end

        function val = eval_far_field(lp, k,density,xhat)
            % Inputs: 
            %   k      : The wave number
            %   density: The solution of the boundary integral equation
            %   xhat   : A 2 x n vector giving the coordinates of where to
            %            valuate the far field.
            
            val = single_layer_far_field(density, xhat(1,:).', xhat(2,:).', lp.RP_Curve, k);

        end

    end

end