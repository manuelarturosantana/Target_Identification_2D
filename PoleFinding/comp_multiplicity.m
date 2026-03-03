% Function which determines the multiplicity and compute the eigenvectors
% corresponding to a certain pol
% Inputs:
%   Fmat  : The matrix whose inverse corresponds to the resolvent.
%   tol   : Tolerance which to accept singular values as multiple
%           eigenvalues
%   nvs   : number of random vectors which to apply the resolvent to.
%           Should be greater than 1
%           Default 3
%   nvmax : Maximum number of random vectors to apply the resolvent to.
%           Default 5


function [mult, evec] = comp_multiplicity(Fmat, tol, nvs, nvmax)
    if nargin < 3
        nvs = 3;
    end

    if nargin < 4
        nvmax = 5;
    end
    if nvs == 1
        warning("comp_multiplicity: using 2 random vectors")
        nvs = nvs + 1;
    end

    N = size(Fmat,1);
    vs = rand(N,nvs);
    
    A = Fmat \ vs;
    
    total_vs = nvs;
    while total_vs <= nvmax
        
        if nargout == 2
            [U,S,~] = svd(A,'econ');
            S = diag(S);
            ind = si_diff(S, tol);
            if ~isnan(ind)
                mult = ind;
                evec  = U(:,1:ind);
                return
            end
   
        % When we only want the index we can save time by computing a 
        % smaller SVD
        else 
           S = svd(rand(total_vs,N) * A,'econ');
           ind = si_diff(S, tol);
            if ~isnan(ind)
                mult = ind;
                return
            end
        end

        A = [A, Fmat \ rand(N,1)];
        total_vs = total_vs + 1;
    end
end

%% Helper Functions
% Function which searches for the index where the singular values start to 
% decay.
% Inputs:
%   S   : A square matrix from an SVD
%   tol : The tolerance to check for the decay. Breaks if (S(ii,ii) / S(ii + 1,ii + 1)) < tol
% Outputs:
%   ind : The index the last large singular value.
function [ind] = si_diff(S, tol)

    if size(S,1) == 1
        ind = nan;
        return
    end

    for ii = 1:size(S,1) - 1
        if (S(ii+1) / S(ii)) < tol
            ind = ii;
            return
        end
    end

    ind = nan;
end














