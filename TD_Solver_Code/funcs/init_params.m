function ps = init_params(ps)
% Initializes the following parameters from the hyperparameters in the
% problem structure.
% a_t  : The time domain function which depends on sigmas and w_0
% B_w  : The full fourier transform of a_t
% sk   : The centers of the time domain windows
% numk : The number of time domain windows.
% ts   : The times to compute the solutions at which are taken uniformly between Tlims.
% xs,ys: The spatial points r = (x,y) to evaluate the solution at.
% ws: The omega values.
%
%% Zero Frequency Parameters Initialized
% numw     : The total number of frequencies used from the graded mesh and fourier continuation
% hs, cs   : The shifts and scales in the affine transformation of each patch to [-1,1]. 
%            note these are stored started with the patch containing 1, and going backwards.
% ccps     : The curtis clenshaw points corrosponding to [-1,1]
% gmend    : The index of the last point in the graded mesh grid. The FC points start at gmend + 1.
%            If ps.do_win_zero == true, then there is overlap in the grid where the windows overlap (w_c, w_c/2);
% ps.win_zero : Function representing the window for the zero frequency component, if do_win_zero is true.
% ps.win_FC : Function representing the half windo for the integral away from zero, if do_win_zero is true.
%
%% Other parameters initialized
%  pou   : A function which is the appropriate partition of unity
%  a_t   : The time domain incident signal without r dependence
%  u_inc : The time domain incident signal with r dependence
%          Note: To implement more signals add them to this file and update the documentation in gen_params

% Now we initialize other parameters based on the inputs

    ps.ts   = linspace(ps.Tlims(1),ps.Tlims(2),ps.numt);
    ps.xs   = linspace(ps.xlims(1),ps.xlims(2),ps.numx);
    ps.ys   = linspace(ps.ylims(1),ps.ylims(2),ps.numy);
    
  
    if ps.wlims(1) <= 0 && ps.wlims(2) >= 0
        [wsGM, ccps,cs,hs] = grated_mesh_grid(ps.w_c,ps.npatch,ps.chebN,ps.q);

        if ps.do_win_zero
            % Note this is different because the FC_grid uses the endpoint, 
            % where as to avoid aliasing the normal DFT grid does not.
            % The window will begin from w_c /2
            wsFC = fc_grid(ps.w_c / 2,ps.wlims(2),ps.numwFC);
            ps.win_zero = @(w) erfc_win(w,ps.w_c);
            ps.win_FC   = @(w) erf_hwin(w,ps.w_c);
        else
            wsFC = fc_grid(ps.w_c,ps.wlims(2),ps.numwFC);
        end
        % In all functions it will be assumed that the graded mesh frequencies 
        % come before the FC frequencies. Also know the graded mesh points come backwards.
        ps.ws = [wsGM(:);wsFC(:)];
        ps.ccps = ccps;
        ps.cs   = cs; ps.hs = hs;
        ps.numw  = length(ps.ws);
        ps.gmend = length(wsGM);
        ps.has_zero_freq = true;
    else
        delta = (ps.wlims(1) + ps.wlims(2)) / 2;
        A     = (ps.wlims(2) - ps.wlims(1)) /2;
        % The number of frequency domain solutions
        if mod(ps.numw,2) ~= 0
            ps.numw = ps.numw + 1;
        end
        M  = ps.numw;
        ms = (-M/2:M/2-1);
        dx = (2*A) / M;
        ps.ws = ms * dx + delta;
        % We set to false in case the ps is reinitialized.
        ps.has_zero_freq = false;
    end

    % Initialize the partition of unity
    if strcmp(ps.pou_name,"bump")
        ps.pou = @(t,H) bump_win(t,H);
    elseif strcmp(ps.pou_name,"erfc")
        ps.pou = @(t,H) erfc_win(t,H);
    else
        error("init_params: Invalid name for partition of unity")
    end

    % Initialize the Gaussian Inputs
    if strcmp(ps.inc_field,"gaussian")
        % We shift the incident field by t_0
        ps.a_t  = @(t) at_gauss(t - ps.t_0,ps.sigmas,ps.w_0);
        ps.B_w  = @(w) exp(-(w - ps.w_0).^2 / ps.sigmas);
        ps.uinc = @(x,y,t) uinc_gauss(x,y,t - ps.t_0,ps.kappa,ps.w_0,ps.sigmas);

        % Because we use the exact Fourier transform there is not windowing needed.
        ps.numk = 1; 
        ps.sk = 0;
    elseif strcmp(ps.inc_field,"chirp")
        % Compute the window length for the chirp.
        center = (ps.csupp(2) + ps.csupp(1)) / 2;
        window_len = (ps.csupp(2) - ps.csupp(1)) / 2;
        ps.a_t = @(t) ps.pou(t - center,window_len).* at_chirp(t);
        ps.uinc = @(x,y,t) ps.a_t(t - [x,y] * ps.kappa(:));

        % The windows only depend on the support of a_k
        ps.sk   = 0:(3 * ps.H / 2):ps.csupp(2);
        ps.numk = length(ps.sk);
    else
        error("The following is not a valid incident field: " + ps.inc_field);
    end

    % An error check for wimag :)
    if ps.wimag >= 0
        error("PS InitParams: Wimag should be negative :)")
    end

end

function vals = at_gauss(t, sigmas, w_0)
    % This is the inverse fourier transform of exp(-(w - w_0)^2/sigmas)
    % sigmas is sigma squared
      vals = exp(-(1/4)*t .* (sigmas .* t + 4i * w_0));
      vals = vals / (sqrt(2) * sqrt(1/sigmas));
      vals = vals * (1/sqrt(2*pi));
end

function val = uinc_gauss(x,y,t,kappa, w_0, sigmas)
    % sigmas is sigma squared
      r = dot(kappa, [x;y]);
      val = exp(-0.25 *(r - t) * (r * sigmas - sigmas * t - 4i * w_0));
      val = val / (sqrt(2) * sqrt(1/sigmas));
      % This came from wolfram Alpha, which by default uses the symmetric
      % normalization, so we follow the paper to have the inverse have the
      % 2pi normalization.
      val = val * (1/sqrt(2*pi));
end

function vals = at_chirp(t)
    g = @(t) 4*t + 6 * cos(t / sqrt(12));
    vals = sin(g(t) + (1/4000)* g(t).^2);
end