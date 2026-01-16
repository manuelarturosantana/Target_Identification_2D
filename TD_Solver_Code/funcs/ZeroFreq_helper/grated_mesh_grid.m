function [intws,ccps,cs,hs] = grated_mesh_grid(w_c,npatch,N,q)
    % Function which created the grated mesh, and Chebyshev points in each grid.
    % See equation 40 in ABL.
    % Inputs:
    %   w_c    : The grid goes from 0 to w_c
    %   npatch : The number of patches in the grated mesh
    %   N      : The degree of the Chebyshev interpolation in each patch
    %   q      : Exponent controling the grating. q > N + 2. 
    %   
    % Outputs: 
    %   intws : The ws in each patch needed for integration. Note due to Clenshaw Curtis
    %           weights the patches share endpoints. Point run backwards going from [1,0)
    %   ccps  : The N+1 Clenshaw curtis points in [-1,1]
    %   cs    : Shift in affine transform to [-1,1] for each patch.
    %   hs    : Scaling in affine transform to [-1,1] for each patch.

    % Create the endpoints of the grated grid. Note there is one more endpoint than
    % patches
    p_endpts = w_c * ((1:npatch+1) / (npatch+1)).^q;

    %Step 2 Computation of needed frequences

    intws = []; % Omegas for integration for each patch execpt the first one
    cs    = [];   hs    = []; 
   
    % Closed-Closed Tchebyshev Mesh
    ccps = cos(pi *(0:N) / N);

    % Compute the needed frequencies for each patch
    for patind = 1:npatch
        c = (p_endpts(patind+1) + p_endpts(patind)) / 2;
        h = (p_endpts(patind+1) - p_endpts(patind)) / 2;
        ws = c + h * ccps;
        % Here we use the closed-closed mesh, so the points overlap. 
        % Because the Tchebyshev mesh if flipped (1 to -1) we don't take
        % the next patches point unless it is the end point.
        if patind ~= npatch
            ws = ws(2:end);
        end
        % Tchbyshev points are backwards
        intws = [ws, intws];
        % Store jacobians backwards too
        cs = [c;cs]; hs = [h;hs];
    end
end