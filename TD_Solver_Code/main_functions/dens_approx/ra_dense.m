function [densities, Rs, pol_c, num_inversions] = ra_dense(ps,lp)
    % Rational Approximation of the density function using set valued aaa, and a recursive
    % approach to the real line sampling.
    % Two notes for improvement:
    %   - At the moment we do not use the extra densities computed, or the rational approximants
    %   to make a refined ws grid to compute the solution from
    %   - In the case of zero frequency the recursive algorithm is only applied to the
    %     interval [w_c,1]
    % Inputs:
    %   ps  : The ubiquitous problem structure
    %   lp  : The layer potential object used for computing the densities.
    % Outputs: 
    %   densities: The computed densities at the origonally frequency points in ps. See 
    %              the long comment below on making this more efficient.
    %   Rs       : Cell array of Vector valued function which takes in w, and returns the rational approximation 
    %              to the density. Length depends on the number of recurssive calls.
    %   pol_c   : Cell array of poles from each subinterval in the recursive calls. 
    %   num_inversions: The total number of inversions for all steps within the adaptive algorithm.
    


    % Compute the densities
    densities = comp_dens(ps,lp);
    % Initially there is one inversion for each w
    num_inversions = ps.numw;

    % New random vector approximation.
    ut = rand(ps.nscal,ps.N,'like',1i);
    multscal = ut * densities;

    % The aaa_sv takes as input a matrix with function values in each row, and functions 
    % varying across the columns. Thus, because we treat each entry of the density as a
    % function value we transpose the densities.
    [~, pol_sv, ~, ~, zms, ~, wms, errvec] = aaa_sv(multscal.',ps.ws,'tol', ps.aaa_tol,'mmax',ps.mmax);

   
    % Case of non-recursive algorithm or no-recursive algorithm, needed.
    % Note errvec does not change size with cleanup.
    if ~ps.use_rec_alg || size(errvec,1) < ps.mmax
        inds = func_inds(ps.ws, zms);
        R = @(zz) reval_sv(zz,zms,densities(:,inds).',wms);

        % Filter poles use wimag, and wlims
        if ps.has_zero_freq
            pol_sv = inrectangle(pol_sv,0,ps.wlims(2),ps.wimag,0);
        else
            pol_sv = inrectangle(pol_sv,ps.wlims(1),ps.wlims(2),ps.wimag,0);
        end
        
        Rs     = {R};
        pol_c = {pol_sv};

    % Recursive Algorithm, split the densities left and right
    else
    % Initialize empty cell array of function handels, pol_c
    Rs = {}; pol_c = {};

    if ps.has_zero_freq
        % In this case we just do the algorithm from w_c to 1. We assume all the poles are
        % within this region
        % FC grid has the left endpoint included
        w_old = ps.ws(ps.gmend+1:end);
        % The Strange thing with the FC code is that it means the ws are a
        % column vector
        w_old = w_old.';
        densr = densities(:,ps.gmend+1:end);
    else
        % Here we add in the last end point, to make dividing in half easier.
        w_old = ps.wlims(2); w_old = [ps.ws, w_old];

        % For which end point we then need to compute the density.
        ps2 = ps; ps2.ws = ps.wlims(2); ps2.numw = 1;
        dens_end  = comp_dens(ps2,lp);
        densr = [densities, dens_end];
    end
    
    % The index of the middle point.
    mid_ind = ceil(length(w_old) / 2);

    % Left Side
    ws_left = w_old(1:mid_ind);
    dens_left = densr(:,1:mid_ind);

    [ws_left, dens_left, num_invl] = refine_ws(ps,lp,ws_left,dens_left);
                         
    [Rs, pol_c, num_lout] =  ra_dense_rec(ps,lp,Rs,pol_c,ws_left,dens_left);

    % Right Side
    ws_right   = w_old(mid_ind:end);
    dens_right = densr(:,mid_ind:end);

    [ws_right, dens_right, num_invr] = refine_ws(ps,lp,ws_right,dens_right);

    % note Rs, pol_c are updated in the function.
    [Rs, pol_c, num_rout] =  ra_dense_rec(ps,lp,Rs,pol_c,ws_right,dens_right);

    num_inversions = num_inversions + num_invl + num_invr + num_lout + num_rout;

    % Potential Feature Combine densities, ws and sort and set unique.
    % TODO: There is a subtelty here. What if one side recursed deeper than the other?
    % The strategy we will first implement is to not use the new densities and w values
    % just stick to the originonal amount.
    % To possible strategies to fix this:
    % (Preferred) Use the Rs to evaluate the densities at a refined equally spaced grid
    %             with grid refinement level based the accuracy of the rs.
    % Or try to figure out how to use only the computed densities by picking those densities
    % which correspond to only the least refined spot. This is difficult however because each
    % recursion may have a different refinement, leading to a need to chase the refinements.
    end
    
end

function [Rs, pol_c, num_inversions] =  ra_dense_rec(ps,lp,Rs,pol_c,ws,densities)
    % Recursive function to compute the densities
    % Inputs:
    %   ps  : The problem structure
    %   lp  : The layer potential object including the curve
    %   Rs  : The cell array containing the R function handels.
    %   pol_c: The cell array containing the poles
    %   ws  : The omega vlaue to compute the rational approximation at
    %   densities: The densities corresponding to the omega values.
    %   Outputs:
    %       Rs: The cell array of function handles to the rational approximations
    %       pol_c: The cell array of poles
    %       num_inversions: The number of inversions to compute the new densities.
    %   WARNING: This function assumes that floor(length(ws)) / 2 < ps.mmax 

    % Compute the scalarization
    disp("Starting New Recursion")
   
    ut = rand(ps.nscal,ps.N,'like',1i);
    multscal = ut * densities;

    [~, pol_sv, ~, ~, zms, ~, wms] = aaa_sv(multscal.',ws,'tol', ps.aaa_tol,'mmax',ps.mmax);

    inds = func_inds(ws, zms);
    R = @(zz) reval_sv(zz,zms,densities(:,inds).',wms);

    % Filter poles use wimag, and wlims
    pol_sv = inrectangle(pol_sv,ws(1),ws(end),ps.wimag,0);

    
    % Base case end
    if length(zms) < ps.mmax
        Rs{end+1}    = R;
        pol_c{end+1} = pol_sv;
        % Base case no new inversions were used.
        num_inversions = 0;
    % Otherwise split and go again.
    else    
        
        % Split ws and densities
        mid_ind = ceil(length(ws) / 2);

        % Left side
        ws_left = ws(1:mid_ind);
        dens_left = densities(:,1:mid_ind);
       
        [ws_left, dens_left, num_invl] = refine_ws(ps,lp,ws_left,dens_left);
                      
        [Rs, pol_c,num_lout] =  ra_dense_rec(ps,lp,Rs,pol_c,ws_left,dens_left);

        % Right Side
        ws_right   = ws(mid_ind:end);
        dens_right = densities(:,mid_ind:end);

        [ws_right, dens_right,num_invr] = refine_ws(ps,lp,ws_right,dens_right);

        [Rs, pol_c,num_rout] =  ra_dense_rec(ps,lp,Rs,pol_c,ws_right,dens_right);

        num_inversions = num_invl + num_invr + num_lout + num_rout;
    end
end

function [ws, densities, num_inversions] = refine_ws(ps,lp,ws,densities)
    % Function which takes even samples between ws, computes the new densities, and combines
    % returns it all together.
    % Inputs:
    %   ps  : The problem structure for computing the densities
    %   lp  : The Layer potential object for computing the densities.
    %   ws  : The previous omega values.
    %   densities: The previous densities
    % Outputs:
    %   ws       : An array of the new ws and the old ones, sorted.
    %   densities: A combined array of the new computed densities and the old ones.
    %   num_inversions: The number of inversions to compute the new densities.


    % First compute the new ws and their corresponding density
    ws_ref = (ws(1:end - 1) + ws(2:end)) / 2;    
    num_inversions = length(ws_ref);

    ps2 = ps; ps2.ws = ws_ref ; ps2.numw = length(ws_ref);
    dens_ref = comp_dens(ps2,lp);
    
    % Now align this with the new RHS
    ws = [ws, ws_ref]; densities = [densities, dens_ref];
    [ws,inds] = sort(ws);
    densities = densities(:,inds);
end


function inds = func_inds(z1,zj)
    % Function which gets the indicies of the support points from the function values.
    % z1 is all the valued passed into AAA, zj are the support points.
    % WARNING: ASSUMES NO DUPLICATE ENTRIES IN Z
    inds = [];
    for ii = 1:length(zj)
        [val,I] = min(abs(z1 - zj(ii)));
        if val ~= 0
            disp("You did what in a cup")
        end
        inds = [inds, I];
    end
end

function r = reval_sv(zz, zj, fj, wj)
    % Evaluate vector valued rational function in barycentric form.
    l = length(zz);
    zv = zz(:);                             % vectorize zz if necessary
    CC = 1./(zv-zj.');   % Cauchy matrix
    r = CC*(wj.*fj)./(CC*wj);

    % Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
    r(isinf(zv),:) = kron(ones(sum(isinf(zv)),1),sum(wj.*fj,1)./sum(wj));

    % Deal with NaN:
    ii = find(isnan(r));
    ii = [mod(ii(:),l),floor(ii(:)/(l+1))+1];
    ii(ii(:,1) == 0) = l;
    for jj = 1:size(ii,1)
        if ( isnan(zv(ii(jj,1))) || ~any(zv(ii(jj,1)) == zj) )
            % r(NaN) = NaN is fine.
            % The second case may happen if r(zv(ii)) = 0/0 at some point.
        else
            % Clean up values NaN = inf/inf at support points.
            % Find the corresponding node and set entry to correct value:
            r(ii(jj,1),ii(jj,2)) = fj(zv(ii(jj,1)) == zj,ii(jj,2));
        end
    end

    % Reshape to input format:
    % r = reshape(r, length(zz),size(fj,2));

end % End of REVAL()