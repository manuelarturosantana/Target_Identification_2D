function fcc_ints = gfcc_int(fvals,ts,ccps,hs,cs,weights,npatch, is_inverse)
    % Compute the integrals using e^(itw) of the graded filon-curtis-clenshaw rule
    % Inputs:
    %   fvals      : The function values at each point in each patch.
    %   ts         : The t values to compute the integrals at. Note these must match the 
    %                indicies of the weights computed for that time. Maybe in the future
    %                I implement a version where you can pass in the corresponding t indices.
    %   ccps       : The regular curtis clenshaw, or closed closed Tchebyshev mesh in [-1,1]
    %   hs, cs     : The linear scaling and shift in the change of varibles in each patch.
    %   weights    : The weights for each patch from comp_fcc_weights.
    %   is_inverse : (Optional) If true use e^(-itw) and not e^(itw). Note that this 
    %                parameter should match the input of comp_fcc_weights
    %
    % Outputs:
    %   fcc_ints  A vector of size numt containg the integrals for each t value

    if nargin <= 7
        is_inverse = true;
    end

    % The degree of the Chebyshev interpolating polynomial
    N = length(ccps) - 1;

    fcc_ints = zeros(1,length(ts));
    
    for tind = 1:length(ts)
        for patind = 1:npatch
            % Each patch as N+1 points, and execpt for the last points shares the
            % N+1 point with its neighbor. Thus patches are reall only N long,
            % then we grab the first point of the neighbor. Also don't forget 1
            % indexing
            fs = squeeze(fvals((patind-1)*N+1:patind*N + 1));
            fs = fs(:).'; % Enforce that fs is a row vector
        
            if is_inverse
                t = -ts(tind);
            else
                t= ts(tind);
            end
            % If t is small, just use Tchebyshev
            if abs(hs(patind) * ts(tind)) < 1/2
                % Here we integrate I_\tilde{k} eqn 2.20 in Dominguez, Gram, Kim
                % Jacobian is multiplied down below.
                fs = fs .* exp(1i * hs(patind)*t*ccps);
                wf = fs .* squeeze(weights(patind,tind,:)).';
            else
                % perform DCT-type 1. First we change from Matlab normalization to
                % the first and last term of the sum being multiplied by 1/2
                fs(1) = (1/sqrt(2)) * fs(1); fs(end) = (1/sqrt(2)) * fs(end);
                % Now we compute the Type 1 dct, changing square root normalization
                coeffs = sqrt(2 / N)*dct(fs,'Type',1);
                % Finally we change normalization on the first and last term
                coeffs(1) = sqrt(2) * coeffs(1); coeffs(end) = sqrt(2) * coeffs(end);
                
                
                % Sum with weights, taking into account special sum.
                wf = coeffs .* squeeze(weights(patind,tind,:)).';
                wf(1) = (1/2) * wf(1); wf(end) = (1/2) * wf(end);
        
            end
            % Multiply by Jacobian and add to the integral
            fcc_ints(tind) = fcc_ints(tind) + ...
            sum(wf) * hs(patind) * exp(1i * t * cs(patind));  
            
        end % For patches
    end % For num t
end