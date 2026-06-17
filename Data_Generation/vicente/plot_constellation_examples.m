% plot_original_constellation_curves.m
clear
clc
close all

% cluster paths
% addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
% addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));
% addpath(genpath("/home/vhojas/chebfun"))

% laptop (windows) paths
addpath(genpath("C:\Users\vicen\Code\rp2d_matlab\src"));
addpath(genpath("C:\Users\vicen\Code\Target_Identification_2D"));

rng(1, "twister");

params = build_params( ...
    'k', 15, ...
    'formulation', 'sl', ...
    'solver', 'direct', ...
    'delta', 0.05, ...
    'npolar', 80, ...
    'n', 12, ...
    'p', 6, ...
    'p_edge', 2);

tags = ["circ_o050"; "ell_o075"];
aa = [1.0; 0.8];
bb = [1.0; 1.6];
openings = [0.50; 0.75];

families = [
    "same_direction"
    "facing_away_origin"
    "random_orientation"
    "facing_origin"
];

opts = struct();
opts.show_nodes = true;
opts.same_color = false;

same_angle = 3*pi/2;
spacing_factor = 2.6;

for iobj = 1:2
    for nobj = 1:5
        centers = make_centers(nobj, max(aa(iobj), bb(iobj)), spacing_factor);

        if nobj == 1
            families_here = "single";
        else
            families_here = families;
        end

        for family = families_here.'
            opening_dirs = make_angles(centers, family, same_angle);
            
            if tags(iobj) == "ell_o075"
                local_opening_dir = 3*pi/2;
                rotations = mod(opening_dirs - local_opening_dir, 2*pi);
            
                geom = elliptic_constellation_geometry( ...
                    aa(iobj), bb(iobj), openings(iobj), centers, rotations);
            else
                geom = cavity_constellation_geometry( ...
                    repmat(aa(iobj), nobj, 1), ...
                    repmat(bb(iobj), nobj, 1), ...
                    repmat(openings(iobj), nobj, 1), ...
                    centers, ...
                    opening_dirs);
            end

            lp = RP2LP(params, geom);
            plot_curve(lp.RP_Curve, opts);

            title(sprintf("%s, N%d, %s", ...
                tags(iobj), nobj, family), ...
                "Interpreter", "none");
        end
    end
end

function centers = make_centers(nobj, r, spacing_factor)

if nobj == 1
    centers = [0, 0];
    return;
end

R = spacing_factor*r/(2*sin(pi/nobj));
theta = linspace(0, 2*pi, nobj + 1).';
theta(end) = [];

centers = [R*cos(theta), R*sin(theta)];

end

function angles = make_angles(centers, family, same_angle)

nobj = size(centers, 1);

switch family
    case {"single", "same_direction"}
        angles = same_angle*ones(nobj, 1);

    case "facing_away_origin"
        angles = atan2(centers(:,2), centers(:,1));

    case "random_orientation"
        angles = 2*pi*rand(nobj, 1);

    case "facing_origin"
        angles = atan2(-centers(:,2), -centers(:,1));

    otherwise
        error("Unknown family: %s", family);
end

angles(vecnorm(centers, 2, 2) == 0) = same_angle;
angles = mod(angles, 2*pi);

end

function geom = elliptic_constellation_geometry(a, b, opening, centers, rotations)

geom = struct();
geom.is_closed = false;
geom.pieces = {};

for j = 1:size(centers, 1)
    gj = ellipticcavity(a, b, opening, rotations(j));

    for k = 1:numel(gj.pieces)
        geom.pieces{end + 1, 1} = translate_piece( ...
            gj.pieces{k}, centers(j, 1), centers(j, 2));
    end
end

end

function Q = translate_piece(P, cx, cy)

Q = P;
Q.x = @(t) P.x(t) + cx;
Q.y = @(t) P.y(t) + cy;

if isfield(P, "xj")
    Q.xj = P.xj + cx;
end

if isfield(P, "yj")
    Q.yj = P.yj + cy;
end

end