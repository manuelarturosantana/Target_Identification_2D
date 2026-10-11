addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
addpath(genpath("/home/vhojas/Code/chebfun"));

% Bound the auto-started parallel pool to what Slurm actually granted
% this job -- some internal function here (single_layer_mat's parfor)
% auto-starts a 'local' profile pool sized to ALL physical cores on the
% node if none exists yet, which is fine on a compute node under Slurm's
% cgroup cpuset but not something to leave unbounded.
num_workers = str2double(getenv("SLURM_CPUS_PER_TASK"));
if isnan(num_workers) || num_workers < 1
    num_workers = maxNumCompThreads;
end
if isempty(gcp("nocreate"))
    job_tag = getenv("SLURM_JOB_ID");
    if isempty(job_tag)
        job_tag = sprintf("pid%d", feature("getpid"));
    end
    pool_dir = fullfile(getenv("HOME"), ".matlab", "local_cluster_jobs", "per_run_" + job_tag);
    if ~exist(pool_dir, "dir")
        mkdir(pool_dir);
    end
    cluster = parcluster("local");
    cluster.JobStorageLocation = pool_dir;
    parpool(cluster, num_workers);
end

vhojas_dir = "/scratch/vhojas/data_planes_2d";
pole_data_file = fullfile(vhojas_dir, "twinjet_plane1.mat");

tilt_deg = 5;
kappa = -[cosd(tilt_deg), sind(tilt_deg)];

ppw = 1; npolar = 160; n_quad = 20; p_corner = 3; p_edge = 2; delta_bc = 0.05;
curve_k = 55;
num_bisections = 2;

gparams = build_twinjet_params('preset', 'plane1');
geom = twinjet_plane_geometry(gparams);

params = build_params('k', curve_k, 'formulation', 'sl', 'solver', 'direct', ...
    'delta', delta_bc, 'npolar', npolar, 'n', n_quad, 'p', p_corner, 'p_edge', p_edge);

C = build_curve(geom, params);
C = refine_curve_wavenumber(C, 2*pi/params.k, ppw);
for jj = 1:num_bisections
    C = bisect_all_patches(C);
    C = build_global_indexing(C);
    C = get_close_points(C);
end

lp = RP2LP(params, geom);
lp.RP_Curve = C;
[xb, yb] = get_global_xy(C);
lp.curve.X = xb; lp.curve.Y = yb; lp.curve.Flag = "vicente open";
bndpts = [xb, yb]';
N = lp.RP_Curve.global.N;
fprintf("N = %d\n", N);

p = 49.800512 - 0.014974i;   % refined pole 1 from the job log

fprintf("\n--- RCOND / min singular value scan vs distance from pole ---\n");
for rad = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    w = p + rad;   % just displace along real direction for a clean scan
    M = lp.bie_mat(w);
    r = rcond(M);
    fprintf("  radius=%.1e  w=%.6f%+.6fi   rcond=%.4e\n", rad, real(w), imag(w), r);
end

M_realp = lp.bie_mat(real(p));
fprintf("  real(p)=%.6f            rcond=%.4e\n", real(p), rcond(M_realp));

fprintf("\n--- Raw density magnitude vs contour radius ---\n");
n_res_w = 10;
rhs_real = -exp(1i*real(p)*(kappa*bndpts)).';
dens_real = lp.bie_mat(real(p)) \ rhs_real;
fprintf("  max|dens_real| = %.4e\n", max(abs(dens_real)));

for cont_rad = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    circ = p + cont_rad*exp(1i*2*pi*(1:n_res_w)/n_res_w);
    mx = 0;
    for cind = 1:n_res_w
        rhs = -exp(1i*circ(cind)*(kappa*bndpts)).';
        d = lp.bie_mat(circ(cind)) \ rhs;
        mx = max(mx, max(abs(d)));
    end
    fprintf("  cont_rad=%.1e  max|res_dens| = %.4e   ratio-to-real = %.4e\n", ...
        cont_rad, mx, mx/max(abs(dens_real)));
end
