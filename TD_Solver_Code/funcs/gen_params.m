function ps = gen_params(varargin)
% Function which generates the parameters for the time domain scattering problem
% Pass in as key, value pairs.
% N            : Number of points used for discretization of the scatterer default 160
% wlims        : The bounds on the values of omega to consider
% numw          : How many frequency points to use, must be even default 80

% Two Ways of choosing spatial evaluation points:
% Way 1:
%   Directly pass in xs and ys as 1D arrays representing pairs of points to evaluate the field at. 
% Way 2:
% Pass in the following:
%   xlims/ylims  : limit on x values to consider 
%   numx, numy   : How many x and y values to use in generating the solution 
%   This will create a grid of (x,y) pairs represented as a numy x numx
%   matrix which is then flattened column wise. So the final answer in the
%   x dimesion will have first dimension of size numy * numx and can be
%   reshaped to be back in the grid.
%   Note that ps.xs, and ps.ys will then contain the flattened mesh grids.
%


% Tlims        : Limits on which time values to consider. Default [0,10]
%                Note that the first window goes from [Tlims(1), H].
%                That is no windowing is implemented for large negative
%                times.
% numt         : How many time values to use. Default 50;
% kappa        : Direction of incident field. Default [0,1];
% H            : Window size for recentering. Default 10
% pou_name     : What function to use for the partition of unity. Options 
%                  - bump: The bump function used in Anderson-Bruno-Lyon. Exact partition of unity.
%                  - erfc: (default)  Pou based off of complementary error function. Exact only to machine percision
%                          but has faster decaying fourier modes, which may lead to better convergence.
%                  Note after parameters are initialized this will be a function not a string.
% is_far_field: If true, the spatial evaluation points are assumed to be
%               on the unit disk and the far field pattern is used. 
% is_open_curve: If true use diffarcs, syntax, otherwise use curve syntax.
% use_FC       : If true use and FC expansion for the integration step. DEPRECIATED: No more FC since window the subtraction
%
% Parameters for 2-D Zero Frequency Content
% These parameters are used automatically if 0 \in [wlims(1),wlims(2)].
% WARNING: The input time domain signal is assumed to be real so the inverse fourier transform
%           Can be computed using only frequencies [0, wlims(2)]. In this case wlims(1) is ignored, other than in checking forzero frequencycontent.
% w_c          : Use Filon-Clenshaw-Curtis rule from [0,w_c], and F_C otherwise. Default 1 (As recommended by Sabhrant Sachan)
% chebN        : Degree of Tchebyshev polynomial to use in each inverval of grated mesh, 
% q            : Exponent on the graded mesh. Must satisfy q > chebN + 2. Default 9.1;
%                 (so N+1 points are used) in each interval. Default 7.
% npatch       : Number of patches to use in the grade mesh. Default 4;
% M_asym       : Value determining size of tridiagonal system to solve when computing Filon-Clenshaw-Curtis Weights. Should be bigger than N/2Default 40
% numwFC       : Number of values to use in FC integration of the interval [w_c, wlims(2)]. Default 80
% do_win_zero  : If true create a window from [0,w_c], and another from [w_c,wlims(2)] so that
%                the singularity expansion can be evaluated on [w_c,wlims(2)], and the logarithmic 
%                integral from [0,w_c] it just computed for all time, as it decays slowly. Default false.
%                Thisis implemented by multiplying BKslow by the window.
% n_eint       : Number of terms to use in the exponential integral which arise from
%                adding the residue term. Default 3000
%               
%
% Incident Field Parameters: pass in as "inc_field". Other paremeters pass in as regular
%   gaussian   : Default field. Corresponds to Gaussian eqn (43) in Anderson-Bruno-Lyon.
%                Note if the gaussian field is used the full windowing strategy is not.
%       w_0          : Center of gaussian in frequency space Default 8.5
%       sigams       : Standard deviation of gaussian Default 0.5
%       t_0          : Shift parameter. If 0 corresponds to eqn 43, else is the fourier 
%                      transform shifted in time of (43) by t_0. Default 0
%   chirp      : The time domain chirp given in (44) of Anderson-Bruno_Lyon
%       csupp  : The chirp is windowed to be supported in [csupp(1),csupp(2)]. Default [1,20]
%
% Parameters for pole finding
% wimag          : Compute poles between wimag and 0. Since the recursive strip algorithm
%                  is implemented for pole finding will fail if large. Default 0.5
% p_numx, p_numy : Number of solutions to use in the x and y direction of the recursive strip. Default 50
% use_secant     : If true use the secant method to refine the pole computations. Default true.
% sec_its        : Number of iterations in secant method. Default 4
% sec_tol        : Tolerance on secant method. Default 1e-13;
%
% Parameters for set valued AAA based real line version of the algorithm
% mmax           : The maximum number of AAA iterations in set valued algorithm. Default 100.
%                  WARNING: set_valued AAA is currently implemented so that mmax is set to min(mmax, numw / 2);
% use_rec_alg    : If true use the recurive algorithm. Default true
% nscal          : Number of random vectors in the left scalarization. Default 20.
% aaa_tol        : Tolerance on set valued AAA algorithm. May need to be relaxed to 1e-10 Default 1e-13. 
%
% Parameters for computing the residues
% n_res_w  : The number of frequency domain solutions to use in each contour integral Default 10
% cont_rad : The radius of the circle around each pole. Default 1e-5
%


%% General Problem Parameters
ps.N               = 160;
ps.numw            = 80;
ps.xlims           = [];
ps.ylims           = [];
ps.numx            = 0;
ps.numy            = 0;
ps.xs              = [];
ps.ys              = [];
ps.Tlims           = [0,10];
ps.numt            = 50;
ps.kappa           = [0,1];
ps.H               = 10;
ps.pou_name        = "erfc";
ps.is_open_curve   = true;
ps.n_eint          = 3000;
ps.is_far_field    = false;

%% Zero Frequency Content Parameters
ps.w_c             = 1;
ps.chebN           = 7;
ps.q               = 9.1;
ps.npatch          = 4;
ps.M_asym          = 40;
ps.numwFC          = 80;
ps.do_win_zero     = false;

%% Incident Field Parameters
ps.inc_field       = "gaussian";
ps.w_0             = 8.5;
ps.sigmas          = 0.5;
ps.t_0             = 0;
ps.csupp           = [1,20];

%% Pole Finding Parameters
ps.wimag           = -0.5;
ps.p_numx          = 50; 
ps.p_numy          = 50;
ps.use_secant      = true;
ps.sec_its         = 4;
ps.sec_tol         = 1e-13;

%% Set Valued AAA Parameters
ps.mmax        = 100;
ps.use_rec_alg = true;
ps.nscal       = 20;
ps.aaa_tol     = 1e-13;

%% Residue Parameters
ps.n_res_w         = 10;
ps.cont_rad        = 1e-5;

if (mod(nargin, 2) ~= 0)
    error("Generate Hyper Params: odd number of arguments, pass in key value pairs")
end

validKeys = ...
    {'N','wlims','numw','xlims','ylims','numx','numy','xs','ys','Tlims','numt', ...
     'kappa','H','pou_name','nworkers','is_open_curve','n_eint', ... 
     'w_c','chebN','q','npatch','M_asym','numwFC','do_win_zero', ... % Zero Frequency
     'inc_field','w_0','sigmas','t_0','csupp', ... % Time domain
     'wimag','p_numx','p_numy','use_secant','sec_its','sec_tol', ... % pole finding
     'mmax','use_rec_alg','nscal','aaa_tol', ... % Residue and pole computation for aaa sv
     'n_res_w','cont_rad','is_far_field'}; % Residue computation 


ii = 1;
while ii < nargin
    key = varargin{ii};
    val = varargin{ii+1};

    if ~any(strcmp(key, validKeys))
        error("Error: '%s' is not a valid field.", key);
    end

    ps.(key) = val;

    if strcmp(key,'kappa')
        ps.kappa = ps.kappa / norm(ps.kappa);
    end

    ii = ii + 2;
end

%% Some sanity checks
    if ps.inc_field == "gaussian" && (ps.w_0 > ps.wlims(2) || ps.w_0 < ps.wlims(1))
        error("Gen Params: The center of the gaussian is not in wlims")
    end
end
