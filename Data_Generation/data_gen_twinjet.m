% paths
addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/chebfun"))

% Experiment params
w_0              = 50;             % The center frequency of the Gaussian
freq_gauss_width = 5;             % Gauss picker how far from w_0 to make the tolerance
gauss_picker_tol = 1e-10;         % At w_0 +/- freq_gauss_width the gaussian will be zero to this tolerance.
numw             = 100;           % Number points used in the integration interval
t_end            = 400;           % Number of points per wave length
t_offset         = 30;            % (Makes sure that the Gaussian starts at 0 :)
num_recievers_per_obs = 1;  
are_recivers_random = false;      % If true generate random angles
angle1 = pi / 3; 
angle2 = pi /3;                  % Angles in the far field to generate data
wimag  = -0.5;                    % Compute poles with imaginary part larger than this value.
save_data_dir = "/home/vhojas/Code/Target_Identification_2D"; % Make sure this ends with a /
% time-domain and AAA derived params
[sigmas, wlims] = gauss_picker(w_0, freq_gauss_width, gauss_picker_tol);
ws = linspace(wlims(1), wlims(2), 1000);
psp = problem_data(...
    'wlims', wlims, ...
    'sigmas', sigmas, ...
    'w_0', w_0, ...
    'xs', 0, ...
    'ys', 0, ...
    't_0', t_offset);
% Time interval to compute the solution in, and number of points in time.
lambda = (2 * pi) / wlims(2);
tlims = [0, t_end];
% Number of wavelengths in the time interval times 15 so we get 15 points 
% per wavelength.
numt = ceil((tlims(2) / lambda)  * 15);

% Geometry params
params = build_params( ...
    'k', w_0 + freq_gauss_width, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);
% Plane params and geometry
geom_params = build_plane_params();
geom = twinjet_plane_geometry(geom_params);

% all geoms (currently only 1 plane geometry)
all_geometries = {geom};
all_geometry_names = {"twinjet"};

% run experiment
start = tic;
for ii = 1:length(all_geometries)
    geom = all_geometries{ii};
    
    lp = RP2LP(params, geom);

    % total number of points in the curve
    N = lp.RP_Curve.global.N;
    
    % This is a rule of thumb for computing the complex resonances.
    % Divide into intervals of length 2. This should save some adaptive
    % steps in the recursive AAA algorithm.
    % Note for large frequency ranges this can be sped up by adapting the
    % discretization size for each subinterval, which we don't do here.
    w_len = wlims(2) - wlims(1);
    num_sub_ints = ceil(w_len / 2.0);
    jump = w_len / num_sub_ints;

    pols_wp = [];
    pstart = tic; 
    for jj = 1:num_sub_ints

        wlims_loc = [wlims(1) + (jj- 1) * jump, wlims(1) + jj * jump];
        % dummy problem data for creating the ws and computing the poles
        psp = problem_data('wlims',wlims_loc,'w_0', (wlims_loc(2) + wlims_loc(1))/2.0,...
            'numw',numw,'xs',0,'ys',0,'wimag',wimag,'N',N);
        
        % one time pole computation
        pols_wp_loc = comp_poles(psp,lp);
        pols_wp_loc = sort(pols_wp_loc,"comparisonMethod","real");
        pols_wp = [pols_wp; pols_wp_loc(:)];
    end

    pol_comp_time = toc(pstart)
    
    % Dummy problem for getting the wlimits for the LU decomp
    psp = problem_data('wlims',wlims,'w_0',w_0,'numw',numw,'xs',0,'ys',0,'wimag',wimag,'N',N);
   
    tic
    [L_rl,U_rl] = comp_bie_LU(psp.ws,lp);
    LU_time_real_line = toc

    tic
    [L_pole, U_pole] = comp_pole_LU(psp,lp,pols_wp);
    LU_time_pols = toc

    if are_recivers_random
        [xs, ys] = ff_points_random(angle1,angle2,num_recievers_per_obs);
    else
        [xs, ys] = ff_points_equally_spaced(angle1,angle2,num_recievers_per_obs);
    end
    
    uff_all     = [];
    r_res_all    = [];

    for rind = 1:num_recievers_per_obs
        x = xs(rind); y = ys(rind);
        psp = problem_data('wlims',wlims,'N',N,'Tlims',tlims,'kappa',[-x,-y],...
        'numt', numt,'numw',numw,'w_0',w_0,'t_0',t_offset,'sigmas',sigmas,'wimag',wimag,...
        'is_far_field',true,'xs',x,'ys',y);
        
        tic
        [uff, ~, r_res] = cts_wpols_LU(psp,lp, pols_wp, L_rl, U_rl,L_pole,U_pole);
        single_reciever_comp_time = toc
        uff = uff(:).'; r_res = r_res(:).';

        uff_all = [uff_all; uff];
        r_res_all = [r_res_all;r_res];
       
    end

    ts = psp.ts;
        
    curveX = lp.curve.X; curveY = lp.curve.Y;
    pols = pols_wp;
    % xs and ys are the reciever locations
    save(fullfile(save_data_dir, all_geometry_names{ii} + ".mat"), ...
        "r_res_all", "uff_all", "pols", "ts", "xs", "ys", ...
        "curveX", "curveY");

end
total_time = toc(start)
