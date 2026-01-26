addpath(genpath("/home/msantana/Progamming/MATLAB_PACKAGES/LayerPotentials"))
addpath(genpath("/home/msantana/Progamming/MATLAB_PACKAGES/DIFFARCS"))
addpath(genpath("/home/msantana/Progamming/Target_Identification_2D"))
addpath(genpath("/home/msantana/Progamming/MATLAB_PACKAGES/chebfun"))
%%
% Main parameters to play with
w_0              = 8;             % The center frequency of the Gaussian
freq_gauss_width = 1;             % Gauss picker how far from w_0 to make the tolerance
gauss_picker_tol = 1e-10;         % At w_0 +/- freq_gauss_width the gaussian will be zero to this tolerance.
numw             = 200;           % Number points used in the integration interval
t_end            = 200;            % Number of points per wave length
t_offset         = 10;            % (Makes sure that the Gaussian starts at 0 :)
num_recievers_per_obs = 1000;  
are_recivers_random = false;      % If true generate random angles
angle1 = 0; 
angle2 = 2 * pi;                  % Angles in the far field to generate data

% Specify the geometries
gap_sizes_ellipse = [0.5:0.5:2];
scat_name = [repmat("circular cavity",1,length(gap_sizes_ellipse)),...
             repmat("elliptic cavity",1,length(gap_sizes_ellipse))];
gap_sizes = [gap_sizes_ellipse, gap_sizes_ellipse];

xshifts = {[0,0], [5, 2.6], [2.7, 5], [2.5, 4.8], [10]};
yshifts = {[2.6, 5], [0, 2.4], [0, 0], [2.5, 4.8], [0]};

constellation_names = {"snowman","triangle","line","diagonal","two_far"};


%% Visual Sanity Check that you have picked good incident and time domain
% field parameters.
[sigmas,wlims] = gauss_picker(w_0,freq_gauss_width,gauss_picker_tol);
ws = linspace(wlims(1),wlims(2),1000);

psp = problem_data('wlims',wlims,'sigmas',sigmas,'w_0',w_0,'xs',0,'ys',0);
% 
% vals = [];
% for ii = 1:length(ws)
%     vals(ii) = psp.B_w(ws(ii));
% end
% % 
% figure(1)
% clf
% semilogy(ws,abs(vals))
% title("Frequency Domain Gaussian Decay")
% xlabel("ws")
% ylabel("Magnitude of Gaussian")
% 
% 
% ts = linspace(1,t_end,1000);
% vals = psp.a_t(ts - t_offset);
% figure(2)
% clf
% semilogy(ts,abs(vals))
% title("Time Domain Gaussian Decay")
% xlabel("ws")
% ylabel("Magnitude of Gaussian")

%% Parameters automatically calculated from the previous ones

% N is the Number of points for the discretization of the scatterer. A good
% rule of thumb in this context is 12 * wlims(2)
N = 12 * wlims(2);   

% Time interval to compute the solution in, and number of points in time.
lambda = (2 * pi) / wlims(2);
tlims = [0, t_end];
% Number of wavelengths in the time interval times 15 so we get 15 points 
% per wavelength.
numt = ceil((tlims(2) / lambda)  * 15);

%%
start = tic;
for ii = 1:length(gap_sizes)
% for ii = 1:1
    ii
    % Build the curve for all directions
    curve = Build_Curve(scat_name(ii), N, gap_sizes(ii));
    lp = OpenSingleLayer(curve);

    % dummy problem data for creating the ws and computing the poles
    psp = problem_data('wlims',wlims,'numw',numw,'xs',0,'ys',0,'wimag',-0.3,'N',N);
    
    % one time pole computation
    pstart = tic; 
    pols_wp = comp_poles(psp,lp);
    pol_comp_time = toc(pstart)
    pols_wp = sort(pols_wp,"comparisonMethod","real");

    for sind = 1:length(xshifts)

        [lp, mc] = setup_constellation(curve, xshifts{sind},yshifts{sind}, false);
        
        % poles for the constellation
        disp("Computing poles for scatterer " + num2str(ii) + " const " + constellation_names{sind});

        % Add one b/c of the base unshifted scattere.
        N_c = N * (length(xshifts{sind}) + 1);

        % A poor implementation nmae of is_open_curve b/c multiple
        % scattering should be has to have is open curve false. I should
        % unify the API
        psp = problem_data('wlims',wlims,'numw',numw,'xs',0,'ys',0,'wimag',-0.3,'N', N_c,'is_open_curve',false);
        pols_c = comp_poles(psp,lp);

        tic
        [L_rl,U_rl] = comp_bie_LU(psp.ws,lp);
        LU_time_real_line = toc
    
        tic
        [L_pole, U_pole] = comp_pole_LU(psp,lp, pols_c);
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
            psp = problem_data('wlims',wlims,'N',N_c,'Tlims',tlims,'kappa',[-x,-y],...
            'numt', numt,'numw',numw,'w_0',w_0,'t_0',t_offset,'sigmas',sigmas,'wimag',-0.3,...
            'is_far_field',true,'xs',x,'ys',y,'is_open_curve',false);
            
            tic
            [uff, ~, r_res] = cts_wpols_LU(psp,lp, pols_c, L_rl, U_rl,L_pole,U_pole);
            single_reciever_comp_time = toc
            uff = uff(:).'; r_res = r_res(:).';
    
            uff_all = [uff_all; uff];
            r_res_all = [r_res_all;r_res];
           
        end
    
        ts = psp.ts;
            
        curveX = mc.xs(1,:); curveY = mc.xs(2,:);
        pols_single = pols_wp;
        pols_constellation = pols_c;
        name = constellation_names{sind};

        % xs and ys are the reciever locations
        sname = strrep(scat_name(ii),' ','');
        save("/scratch/msantana/time_domain_data_2D/constellation_tests/"+sname + num2str(gap_sizes(ii)) + "cind" + num2str(sind) + ".mat",...
            "r_res_all","uff_all","pols_single","pols_constellation","ts","xs","ys","curveX","curveY","name")
    end

end
total_time = toc(start)