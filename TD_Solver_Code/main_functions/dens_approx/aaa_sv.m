function [r, pol, res, zer, z, f, w, errvec] = aaa_sv(F, varargin)
    % CODE TAKEN FROM: https://academic.oup.com/imajna/article/42/2/1087/6139194
    %
    % Couple of Modifications to their code
    % Cleanup proceedure, following the origonal AAA proceedure
    % Enforced an edgecase that mmax = min(floor(M / 2),mmax), for small data sets.
    % UPDATE TO INCLUDE INITIAL PARAMETER, OR NO CLEAN UP IF NOT SUFFICIENTLY RESOVLED.
    
    
    %AAA_SV - Reworked version of the aaa function from Chebfun (http://www.chebfun.org/), 
    % see also: 'N AKATSUKASA , Y., S ÈTE , O. & T REFETHEN , L. N. (2016) 
    % The AAA algorithm for rational approximation. Preprint arXiv:1612.00337'
    %
    % Computes a AAA rational approximation where F can have multiple
    % outputs
    %   R = AAA(F, Z) computes the AAA rational approximant R (function handle) to
    %   data F on the set of sample points Z.  F may be given by its values at Z,
    %   or as a function handle or a chebfun.
    %
    %   [R, POL, RES, ZER] = AAA(F, Z) returns vectors of poles POL,
    %   residues RES, and zeros ZER of R.
    %
    %   [R, POL, RES, ZER, ZJ, FJ, WJ] = AAA(F, Z) also returns the vectors
    %   of support points ZJ, function values FJ, and weights WJ of the
    %   barycentric representation of R.
    %
    %   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC] = AAA(F, Z) also returns the
    %   vector of errors ||f-r||_infty in successive iteration steps of AAA.
    %
    %   R = AAA(F, Z, NAME, VALUE) sets the following parameters:
    %   - 'tol', TOL: relative tolerance (default TOL = 1e-13),
    %   - 'mmax', MMAX: maximal number of terms in the barycentric representation
    %       (default MMAX = 100).
    %   - 'dom', DOM: domain (default DOM = [-1, 1]). No effect if Z is provided.
    %
    %   One can also execute R = AAA(F), with no specification of a set Z.
    %   This is equivalent to defining Z = LINSPACE(DOM(1), DOM(2), LENGTH(F)) if F
    %   is a vector (by default DOM = [-1, 1]).
    %   If F is a function handle or a chebfun, AAA attempts to resolve F on its
    %   domain DOM.  By default, DOM = [-1, 1] for a function handle, and
    %   DOM = F.DOMAIN([1, END]) for a chebfun.
    %
    
    % Parse inputs:
    [F, Z, M, nF, dom, tol, mmax, needZ, mmax_flag, cleanup_tol, cleanup_flag] = ...
    parseInputs(F, varargin{:});
    
    if ( needZ )
        % Z was not provided.  Try to resolve F on its domain.
        [r, pol, res, zer, z, f, w, errvec] = ...
            aaa_autoZ(F, dom, tol, mmax, mmax_flag);
        return
    end
    
    % Scale the functions
    normF = max(abs(F),[],1);
    F = F./normF;
    
    Fclean = F;
    Zclean = Z;
    
    % Left scaling matrix:
    SF = spdiags(F, 0:-M:-M*(nF-1), M*nF, M);
    
    
    % Initialize values
    F = F(:);
    R = mean(F);
    errvec = zeros(mmax,1);
    z = zeros(mmax,1);
    f = zeros(mmax,nF);
    ind = zeros(mmax,nF);
    H = zeros(mmax,mmax-1);
    S = zeros(mmax,mmax-1);
    
    Q = zeros(M*nF,0);
    C = zeros(M,0);
    
    % AAA iteration:
    for m = 1:mmax
        [errvec(m),loc] = max(abs(F-R));               % Select next support point where error is largest
        if ( errvec(m) <= tol )
            m = m-1;
            break
        end
    
        loc = mod(loc,M);
        ind(m,:) =  loc + (M*(loc==0):M:(nF-1+(loc==0))*M);  % Get indices of the z_i
        z(m) = Z(ind(m,1));                           % Add interpolation point
        f(m,:) = F(ind(m,:));                         % Add function values
        C(:,end+1) = 1./(Z - z(m));                   % Get column of the Cauchy matrix.
        C(ind(1:m,1),m) = 0;                          % Set the selected elements to 0
    
        v = C(:,m)*f(m,:);                            % Compute the next vector of the basis.
        v = SF*C(:,m)-v(:);
    
        % Update H and S to compensate for the removal of the rows
        q = Q(ind(m,:),1:m-1);
        q = q*S(1:m-1,1:m-1);
        ee = eye(m-1,m-1)-q'*q;
        ee(1:size(ee,1)+1:end) = real(diag(ee));
        Si = chol(ee);
        H(1:m-1,1:m-1) = Si*H(1:m-1,1:m-1);
        S(1:m-1,1:m-1) = S(1:m-1,1:m-1)/Si;
        S(m,m) = 1;
        Q(ind(1:m,:),:) = 0;
    
        nv = norm(v);
        H(1:m-1,m) = Q'*v;
        H(1:m-1,m) = S(1:m-1,1:m-1)'*H(1:m-1,m);
        HH = S(1:m-1,1:m-1)*H(1:m-1,m);
        v = v-(Q*HH);
        H(m,m) = norm(v);
        % Reorthoganlization is necessary for higher precision
        it = 0;
        while (it < 3) && (H(m,m) < 1/sqrt(2)*nv)
            h_new = S(1:m-1,1:m-1)'*(Q'*v);
            v = v - Q*(S(1:m-1,1:m-1)*h_new);
            H(1:m-1,m) = H(1:m-1,m) + h_new;
            nv = H(m,m);
            H(m,m) = norm(v);
            it = it+1;
        end
        v = v/H(m,m);
    
        % Add v
        Q(:,end+1) = v;
    
        % Solve small least squares problem with H
        [~,~,V] = svd(H(1:m,1:m));
        w = V(:,end);
    
        % Get the rational approximation
        N = C*(w.*f(1:m,:));       % Numerator
        D = C*(w.*ones(m,nF));     % Denominator
        R = N(:)./D(:);
        R(ind(1:m,:)) = F(ind(1:m,:));
    end
    
    f = f(1:m,:).*normF;
    w = w(1:m);
    z = z(1:m);
    aaa_final_deg = m
    aaa_final_err = errvec(m)
    errvec = errvec(1:m).*normF;
    
    % Note: When M == 2, one weight is zero and r is constant.
    % To obtain a good approximation, interpolate in both sample points.
    if ( M == 2 )
        z = Z;
        f = F;
        w = [1; -1];       % Only pole at infinity.
        w = w/norm(w);   % Impose norm(w) = 1 for consistency.
        errvec(2) = 0;
    end
    
    % Remove support points with zero weight:
    I = find(w == 0);
    z(I) = [];
    w(I) = [];
    f(I,:) = [];
    
    % Construct function handle:
    r = @(zz) reval(zz, z, f, w);
    
    % Compute poles, residues and zeros:
    [pol, res, zer] = prz(r, z, f, w);
    
    % MY ADDITION: ADD CLEANUP STEP
    if cleanup_flag
        f = f ./ normF;
        [r, pol, res, zer, z, f, w] = ...
            cleanup(r, pol, res, zer, z, f, w, Zclean, Fclean, cleanup_tol);
        f = f .* normF;
    end
    
    end % of AAA()
    
    
    
    %% parse Inputs:
    
    function [F, Z, M, nF, dom, tol, mmax, needZ, mmax_flag, cleanup_tol, cleanup_flag] = ...
        parseInputs(F, varargin)
    % Input parsing for AAA.
    
    % Check if F is empty:
    if ( isempty(F) )
        error('No function given.')
    end
    
    % Sample points:
    if ( ~isempty(varargin) && isfloat(varargin{1}) )
        % Z is given.
        Z = varargin{1};
        if ( isempty(Z) )
            error('If sample set is provided, it must be nonempty.')
        end
        varargin(1) = [];
    end
    
    % Set defaults for other parameters:
    tol = 1e-13;        % Relative tolerance.
    mmax = 100;         % Maximum number of terms.
    
    if ~iscell(F), F = {F}; end
    
    % Domain:
    if ( isa(F{1}, 'chebfun') )
        dom = F{1}.domain([1, end]);
    else
        dom = [-1, 1];
    end
    mmax_flag = 0;
    cleanup_flag = 1; % Clean up automatically set to on.
    cleanup_set = 0;               % Checks if cleanup_tol manually specified
    % Check if parameters have been provided:
    while ( ~isempty(varargin) )
        if ( strncmpi(varargin{1}, 'tol', 3) )
            if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
                tol = varargin{2};
                if ( ~cleanup_set && tol > 0 ) % If not set, set cleanup_tol to tol
                    cleanup_tol = tol;
                end
            end
            varargin([1, 2]) = [];
    
        elseif ( strncmpi(varargin{1}, 'mmax', 4) )
            if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
                mmax = varargin{2};
                mmax_flag = 1;
            end
            varargin([1, 2]) = [];
    
        elseif ( strncmpi(varargin{1}, 'dom', 3) )
            if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 2]) )
                dom = varargin{2};
            end
            varargin([1, 2]) = [];
            if ( isa(F, 'chebfun') )
                if ( ~isequal(dom, F{1}.domain([1, end])) )
                    warning('CHEBFUN:aaa:dom', ...
                        ['Given domain does not match the domain of the chebfun.\n', ...
                        'Results may be inaccurate.'])
                end
            end 
        elseif ( strncmpi(varargin{1}, 'cleanuptol', 10) )
            if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
              cleanup_tol = varargin{2};
              cleanup_set = 1;
            end
            varargin([1, 2]) = [];
    
        elseif ( strncmpi(varargin{1}, 'cleanup', 7) )
            if ( strncmpi(varargin{2}, 'off', 3) || ( varargin{2} == 0 ) )
                cleanup_flag = 0;
            end
            varargin([1, 2]) = [];
    
        else
            error('CHEBFUN:aaa:UnknownArg', 'Argument unknown.')
        end
    end
    
    % Deal with Z and F:
    if ( ~exist('Z', 'var') && isfloat(F{1}) )
        % F is given as data values, pick same number of sample points:
        Z = linspace(dom(1), dom(2), size(F{1},1)).';
    end
    
    if ( exist('Z', 'var') )
        % Z is given:
        needZ = 0;
    
        % Work with column vector:
        Z = Z(:);
        M = length(Z);
    
        %%%% MY MODIFICATION 
        % This least squares problem is only well posed if M /2 > mmax
    
        mmax = min(floor(M / 2),mmax);
        %%%%%
    
        % Function values:
        for it = 1:length(F)
            if ( isa(F{it}, 'function_handle') || isa(F{it}, 'chebfun') )
                % Sample F on Z:
                F{it} = F{it}(Z);
            elseif ( isnumeric(F{it}) )
                % Work with column vector and check that it has correct length.
                if ( size(F{it},1) ~= M )
                    error('CHEBFUN:aaa:lengthFZ', ...
                        'Inputs F and Z must have the same length.')
                end
            else
                error('CHEBFUN:aaa:UnknownF', 'Input for F not recognized.')
            end
        end
        F = cell2mat(F);
        nF = size(F,2);
    else
        % Z was not given.  Set flag that Z needs to be determined.
        % Also set Z and M since they are needed as output.
        needZ = 1;
        Z = [];
        M = length(Z);
        nF = length(F);
    end
    end % End of PARSEINPUT().
    
    
    %% Evaluate rational function in barycentric form.
    
    function r = reval(zz, zj, fj, wj)
    % Evaluate rational function in barycentric form.
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
    
    end % End of REVAL().
    
    
    %% Compute poles, residues and zeros.
    
    function [pol, res, zer] = prz(r, zj, fj, wj)
    % Compute poles, residues, and zeros of rational function in barycentric form.
    m = length(wj);
    
    % Compute poles via generalized eigenvalue problem:
    B = eye(m+1);
    B(1,1) = 0;
    E = [0 wj.'; ones(m, 1) diag(zj)];
    pol = eig(E, B);
    % Remove zeros of denominator at infinity:
    pol = pol(~isinf(pol));
    
    % Compute residues via discretized Cauchy integral:
    dz = 1e-5*exp(2i*pi*(1:4)/4);
    
    pp = pol+dz;
    rvals = r(pp(:));
    res = zeros(length(pol),size(rvals,2));
    for it = 1:size(rvals,2)
        res(:,it) = reshape(rvals(:,it),[],4)*dz.'/4;
    end
    
    % Compute zeros via generalized eigenvalue problem:
    for it = 1:size(fj,2)
        E = [0 (wj.*fj(:,it)).'; ones(m, 1) diag(zj)];
        zer{it} = eig(E, B);
        % Remove zeros of numerator at infinity:
        zer{it} = zer{it}(~isinf(zer{it}));
    end
    end % End of PRZ().
    
    
    %% Automated choice of sample set
    
    function [r, pol, res, zer, zj, fj, wj, errvec] = ...
        aaa_autoZ(F, dom, tol, mmax, mmax_flag)
    %
    
    % Flag if function has been resolved:
    isResolved = 0;
    
    % Main loop:
    for n = 5:14
        % Sample points:
        % Next line enables us to do pretty well near poles
        Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
        [r, pol, res, zer, zj, fj, wj, errvec] = aaa_sv(F, Z, 'tol', tol, ...
            'mmax', mmax);
    
        % Test if rational approximant is accurate:
        if ~iscell(F), F = {F}; end
    
        r_val = r(Z);
        Zrefined = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), ...
                round(1.5 * (1 + 2^(n+1)))).';
        r_val_r = r(Zrefined);
    
        % Final check that the function is resolved, inspired by sampleTest().
        % Pseudo random sample points in [-1, 1]:
        xeval = [-0.357998918959666; 0.036785641195074];
        % Scale to dom:
        xeval = (dom(2) - dom(1))/2 * xeval + (dom(2) + dom(1))/2;
        r_val_f = r(xeval);
    
        isResolved = zeros(1,length(F));
    
        for it = 1:length(F)
            reltol = tol * norm(F{it}(Z), inf);
    
            err(1,1) = norm(F{it}(Z) - r_val(:,it), inf);
            err(2,1) = norm(F{it}(Zrefined) - r_val_r(:,it), inf);
    
            if ( all(err < reltol) )
                if ( norm(F{it}(xeval) - r_val_f(:,it), inf) < reltol )
                    isResolved(it) = 1;
                end
            end
        end
        if all(isResolved)
            break
        end
    end
    
    if ( ( prod(isResolved) == 0 ) && ~mmax_flag )
        warning('CHEBFUN:aaa:notResolved', ...
            'Function not resolved using %d pts.', length(Z))
    end
    
    end % End of AAA_AUTOZ().
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CLEANUP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [r, pol, res, zer, z, f, w] = ...
        cleanup(r, pol, res, zer, z, f, w, Z, F, cleanup_tol) 
        % Important inputs:
        %   pol : The computed poles of length num_pol
        %   res : The residue of size num_pol x num_random scale (ell)
        %   z   : The support points of size m
        %   f   : The function values of size m x ell
        %   Z   : All the initial Z points
        %   F   : All the functions values as a num(Z) x ell
    % Remove spurious pole-zero pairs.
    
    % Origonal AAA perscription.
    % Note mean returns a row vector containing the mean of each column
    if any(F)
        geometric_mean_of_absF = exp(mean(log(abs(F(F~=0)))));
    else
        geometric_mean_of_absF = 0;
    end
    
    
    Zdistances = NaN(size(pol));
    for j = 1:length(Zdistances)
        Zdistances(j) = min(abs(pol(j)-Z));
    end
    
    % abs(res)./Zdistances divides the residue of each function at each pol by the distance to
    % the pole, assuming that Zdistance is a column vector. 
    % Then we compare each residue of each pole of each function as in the scalar case
    % We then sum to see if all the residues are small;
    inds = sum(abs(res)./Zdistances < cleanup_tol * geometric_mean_of_absF,1);
    % If all the residues are small then the sum of the logicals evaulated to the number of 
    % of poles
    ii = find(inds == length(pol));
    ni = length(ii);
    if ( ni == 0 )
        return
    elseif ( ni == 1 )
        warning('AAA:Froissart','1 Froissart doublet');
    else
        warning('AAA:Froissart',[int2str(ni) ' Froissart doublets']);
    end
    
    % For each spurious pole find and remove closest support point:
    % and all function values corresponding to that
    for j = 1:ni
        azp = abs(z-pol(ii(j)));
        jj = find(azp == min(azp),1);
    
        % Remove support point(s):
        z(jj)  = []; 
        f(jj,:) = [];
    end
    
    % Remove support points z from sample set:
    for jj = 1:length(z)
        F(Z == z(jj),:) = [];
        Z(Z == z(jj)) = [];
    end
    m = length(z);
    M = length(Z);
    
    % A slightly slow implementation but we build the Loewner matrix in a loop.
    C = 1./(Z-z.'); % Cauchy Matrix
    A = [];
    for ii = 1:size(f,2)
        Fii = F(ii,:); Fii = Fii(:);
        SF = spdiags(Fii, 0, M, M);
        Sf = diag(f(ii,:));
                  % Cauchy matrix.
        A_i = SF*C - C*Sf;              % Loewner matrix for this function.
        A = [A; A_i];
    end
    
    % IF THIS IS SLOW, THEN CONSIDER DOING A * A WHICH HAS THE SAME SINGULAR VALUES. HOTCHMAN 2017
    % Solve least-squares problem to obtain weights:
    [~, ~, V] = svd(A, 0);
    w = V(:,m);
    
    
    % Build function handle and compute poles, residues and zeros:
    % r = @(zz) reval(zz, z, f, w);
    [pol, res, zer] = prz(z, f, w);
    
    end % End of CLEANUP.
    