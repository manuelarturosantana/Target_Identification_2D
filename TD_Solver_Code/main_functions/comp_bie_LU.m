function [L_cell, U_cell] = comp_bie_LU(ws,lp)
% Functions which computes an LU factorization at all the requested
% frequency points.
% Inputs:
%   ps  : The problem stucture
%   lp  : The layer potential object
%
% Output:
%   L_cell: Cell array corresponding the L matrices
%`  U_cell: Cell array corresponding the U matrices
    
    numw = length(ws);
    L_cell = cell(numw,1);
    U_cell = cell(numw,1);

    %   - For each omega solve the BIE with the incident plane wave. 
    parfor wind = 1:numw
        [L, U] = lu(lp.bie_mat(ws(wind)));
        L_cell{wind} = L; U_cell{wind} = U;
    end

end