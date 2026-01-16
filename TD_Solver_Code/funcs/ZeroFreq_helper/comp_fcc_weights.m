function weights = comp_fcc_weights(ts,hs,npatch,N,M, is_inverse)
    % Function which computes the Filon-Clenshaw Curtis weights for each patch in the graded mesh.
    % By default if the function is slowly oscillating enough (abs(ht) < 1/2) then the inverse
    % is computed via the usual clenshaw-curtis rule. The arbitrary choice of (1/2) follows
    % Filon-Clenshaw-Curtis rules for highly-oscillatory integrals with algebraic singularities and stationary points
    % by Dominguez, Graham, and Kim
    % Inputs: 
    %   ts         : The time to be evaluated at. Can be negative, (as in
    %                 the case of ts - sk in the windowing and recentering)
    %   hs         : The scaling of the patches to -1,1
    %   npatch     : The number of patches
    %   N          : The degree of the Tchebyshev interpolant used.
    %   M          : The value for the asymptotic expansion used to solve
    %                the tridaigonal system.
    %   is_inverse : If true compute weights for inverse fourier transform i.e. using e^(-iwt).
    %                If false use e^(iwt)
    %
    % Outputs:
    %   weights : A matrix of size npatch x numt x (N + 1) containing the weights for each patch,t, and polynomial.
    
    weights = zeros(npatch,length(ts),N+1);

    [~,cc_weights] = clenshaw_curtis(N); 

    for tind=1:length(ts)
        for patind = 1:npatch
            % The rescaling to [-1,1] rescales the time
            curr_time = hs(patind) * ts(tind);
            % If slowly enough oscillating just compute with clenshaw_curtis.
            if abs(curr_time) < 1/2
                weights(patind,tind,:) = cc_weights;
            else
                if is_inverse
                    curr_time = -curr_time;
                end
                weights(patind,tind,:) = fcc(curr_time,N,M);
            end
        end
    end
end




function weights = fcc(t,N,M)
    % Function which computes the weights for the Filon-Clenshaw-Curtis
    % integration. Notation for the most part follows: 
    % Stability and error estimates for Filon-Clenshaw-Curtis rules for 
    % highly-oscillatory integrals V. Domínguez I.G. Graham V.P. Smyshlyaev
    %
    % With the exception that we have switched t to k

    % Inputs:
    %   t   : The value in the complex exponential. Our convention is that the 
    %         sign of t determines the direction of the fourier transform.
    %         that is, weights or e^(its) are computed, t being positive or
    %         negative.
    %   N   : Degree of the Tchebysev interpolatn being used
    %   M   : Value for asymptotic expanion. The larger M is the more accurate
    %         the wieghts will be if t < N. However the a tridiagonal linear
    %         system of size M x M must be solved. M must satisfy M >=
    %         max(ceil(t)/2, N/2).
    % Outputs:
    %   weights: The t dependent integration weights. Column vector of size 1xN


    % Determine needed gammas based on cases in the paper.
    if N < abs(t)
        % Gammas from 0 to N are needed
        gamma_size = N+1;
    else
        % Gammas from 0 to 2M - 1 are needed
        gamma_size  = 2*M;
    end

    % Compute gammas
    gammas = zeros(gamma_size,1);
    gammas(1:2:end) = (2 * sin(t)) / t;      % even_terms (matlab is 1 indexed)
    gammas(2:2:end) = (2 * cos(t)) / (1i*t); % odd terms.

    % In this case we can just use recursive algorithm
    if N < abs(t)
        rhos    = rhok_rec(t,gammas,N);
        weights = fcc_weights(t,N,gammas,rhos);
    % Here use recursive algorithm and the tridiagonal phase
    else 
        % Asymptotic case
        % Use recursive algorithm for weights indexed less than n_0
        % Note the first index can always be computed stably, so we
        % make sure n_0 is at least 2.
        n_0 = max(ceil(abs(t)),2);
        rhos = rhok_rec(t,gammas,n_0-1);

        rhos2 = asym_tridiag(n_0,M,t,gammas,rhos(end));
        rhos  = [rhos;rhos2];

        weights = fcc_weights(t,N,gammas,rhos);
        return
    end
end

function rhos = rhok_rec(t, gammas, end_ind)
    % Function which uses the recurrence relation to compute rho_k
    % see equations 33a-33c
    % Inputs:
    %   t      : The exponential integration frequency.
    %   gammas: The gamma values
    %   end_ind: The index of the final value of rho to compute until >= 1
    % 
    % Outputs:
    %   rhos : The rho values from 1 to end_ind
    
    rhos = zeros(end_ind,1);
    % Note they index rho by 1 and gamma by 0 in the paper :(
    rhos(1) = gammas(1);
    % Edge case
    if end_ind == 1
        return;
    end
    rhos(2) = 2*gammas(2) - (2/(1i*t)) * gammas(1);

    for rho_ind = 3:end_ind
        rhos(rho_ind) = 2 * gammas(rho_ind) - ...
            ((2*(rho_ind - 1) )/(1i*t)) * rhos(rho_ind-1) + rhos(rho_ind-2);
    end
end

function weights = fcc_weights(t,N,gammas,rhos)
    % Function which uses the recurrance to compute the weights
    % See equation 34
    %   t      : The exponential integration frequency.
    %   N      : The number of weight values to compute
    %   gammas : The gamma values
    %   rhos   : The computed rhos

    % N+1 because weights go from 0:N
    weights = zeros(N+1,1);
    
    weights(1) = rhos(1);
    
    for wind = 2:(N+1)
        % wind is n + 1 due to one indexing, note there is no rho_0
        weights(wind) = gammas(wind) - ((wind - 1) / (1i * t))*rhos(wind - 1);
    end
end

function rhos = asym_tridiag(n_0,M,t,gammas,rho_n0m1)
    % Inputs:
    %   n0   : ceil(t). The value until which the recursive algorithm
    %         computed.
    %   M    : The parameter determining the size of the asymptotic
    %           expansion
    %  gamma : The computed values of gamma
    %  rho_n0m1 : The last computed term of rho

    % Create Am and bm following formula 25
    matsize = 2*M - n_0 ; 
    d = 2 * (n_0:(2*M-1)) / (1i*t);
    bin = [-ones(matsize,1),d(:), ones(matsize,1)];
    Am = spdiags(bin,-1:1,matsize,matsize);

    % Gamma is zero indexed hence we have to move up by 1
    bm = 2*gammas(n_0+1:2*M);
    bm(1) = bm(1) + rho_n0m1;

    rho_2Mt = comp_rho_asym(t,M);
    bm(end) = bm(end) - rho_2Mt;
    rhos = Am \ bm;

end

function rho_2Mt = comp_rho_asym(t,M)
    % Following the formulas in (36), (37) compute the asymptotic value of
    % rho

    p0 = 1/(2*M);
    p1 = t / (2*M)^3;
    p2 = 3*t^2 / (2*M)^5;
    p3 = ((15*t^2 - 4*M^2)*t) / (2*M)^7;
    p4 = ((105*t^2 - 60*M^2)*t^2) / (2*M)^9;
    p5 = ((945*t^4 - 840*t^2*M^2 + 16*M^4)*t) / (2*M)^11;

    rho_2Mt = sin(t) * (p0 - p2 + p4) + cos(t) * (p1 - p3 + p5);
    rho_2Mt = 2i * rho_2Mt;
end