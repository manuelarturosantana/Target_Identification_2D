% diag_quadjet_plane5_junction.m
%
% Diagnoses the reported "holes" at the cabin/front-wing junction in
% quadjet_plane5. quadjet_plane_geometry.m's get_wing_params fits the
% left wing-lead ellipse via fsolve so it passes exactly through the
% fuselage front2 curve's endpoint (xfl,yfl), and the right wing is an
% exact Y-mirror -- so a genuine open gap there is unlikely for any
% preset. This checks the two remaining candidate causes instead, using
% the ACTUAL pieces returned by quadjet_plane_geometry.m (no
% reimplementation of its internals):
%   (a) the wing-lead curve crossing the fuselage curve a SECOND time
%       near the joint (a real self-intersection, not just failing to
%       meet at the declared endpoint), and
%   (b) a shallow/near-tangent join angle between the two curves, which
%       enclosed_polyshape.m's own docstring warns can fool its
%       half-edge angle-sort face tracer into carving out a spurious
%       face.
% Compares plane5 (b=0.10, the only preset that changes fuselage
% half-thickness) against plane1 (b=0.08, default) as a baseline.
%
% Piece order out of quadjet_plane_geometry.m (closed_front=false) is
% fixed for indices 1-9: front1, front2, back, lwing_p1, lwing_p2,
% lwing_p3, lwing_p4(trail), lwing_p5(trail), lwing_p6(trail) -- none of
% these ever return empty, unlike the engine pieces further down the
% list. So pieces{2}=front2 and pieces{4}=lwing_p1 reliably.

clear
clc

addpath(genpath("/home/vhojas/Code/rp2d_matlab/src"));
addpath(genpath("/home/vhojas/Code/Target_Identification_2D"));

for preset = ["plane1", "plane2", "plane3", "plane4", "plane5"]
    fprintf("\n================ quadjet %s ================\n", preset);
    gparams = build_quadjet_params('preset', char(preset));
    fprintf("b (fuselage half-thickness) = %.4f, front_angle = %.4f, closed_front = %d\n", ...
        gparams.b, gparams.front_angle, gparams.closed_front);

    geom = quadjet_plane_geometry(gparams);
    if gparams.closed_front
        % pieces = {front; back; lwing_p1; ...} -- only one front piece.
        frontA = geom.pieces{1};
        frontB = geom.pieces{1};
        lwing_p1 = geom.pieces{3};
        xfl = frontA.x(frontA.tf);
        yfl = frontA.y(frontA.tf);
    else
        % pieces = {front1; front2; back; lwing_p1; ...}
        frontA = geom.pieces{1};
        frontB = geom.pieces{2};
        lwing_p1 = geom.pieces{4};
        xfl = frontB.x(frontB.tf);
        yfl = frontB.y(frontB.tf);
    end
    front1 = frontA;
    front2 = frontB;
    tstar = lwing_p1.t0;

    % --- endpoint exactness (sanity check on the fsolve fit) ---
    dx = lwing_p1.x(tstar) - xfl;
    dy = lwing_p1.y(tstar) - yfl;
    mismatch = hypot(dx, dy);
    % enclosed_polyshape.m clusters endpoints into shared vertices using
    % tol = 1e-8*scale (scale ~ max coordinate magnitude, ~O(1) here) --
    % if the mismatch exceeds that, the two endpoints are NOT merged
    % into one vertex and the wing loop gets a spurious dangling end
    % right where it should meet the fuselage.
    approx_tol = 1e-8;
    fprintf("Endpoint mismatch |wing(tstar) - front2(tf)| = %.3e  (enclosed_polyshape's clustering tol ~ %.0e)\n", ...
        mismatch, approx_tol);
    if mismatch > approx_tol
        fprintf("  ** Exceeds clustering tolerance: this vertex will NOT be merged -- wing loop gets a dangling end here. **\n");
    end

    % --- (b) join angle between the two curves at the shared point ---
    v_front = [front2.xp(front2.tf), front2.yp(front2.tf)];
    v_wing  = [lwing_p1.xp(tstar), lwing_p1.yp(tstar)];
    v_front = v_front / norm(v_front);
    v_wing  = v_wing / norm(v_wing);
    join_angle_deg = acosd(max(-1, min(1, dot(v_front, v_wing))));
    fprintf("Tangent angle between front2 and wing-lead at joint = %.3f deg\n", join_angle_deg);
    fprintf("  (180 = smooth continuation, 0 = curves double back; near 0 or 180 is near-degenerate.)\n");

    % --- (a) does the wing-lead curve cross the fuselage curve again
    %     near the joint, beyond the shared endpoint? ---
    t_near = linspace(tstar, tstar + 0.15*(lwing_p1.tf - tstar), 2000);
    wx = lwing_p1.x(t_near); wy = lwing_p1.y(t_near);

    tf1 = linspace(front1.t0, front1.tf, 4000);
    tf2 = linspace(front2.t0, front2.tf, 4000);
    fx = [front1.x(tf1), front2.x(tf2)];
    fy = [front1.y(tf1), front2.y(tf2)];

    min_d = zeros(size(wx));
    for i = 1:numel(wx)
        min_d(i) = min(hypot(fx - wx(i), fy - wy(i)));
    end
    away_from_joint = t_near > tstar + 0.01*(lwing_p1.tf - tstar);
    second_min = min(min_d(away_from_joint));
    fprintf("Closest approach of wing-lead curve to fuselage curve, excluding the joint itself: %.3e\n", second_min);
    if second_min < 1e-3
        fprintf("  ** SUSPICIOUS: wing-lead curve comes very close to (or crosses) the fuselage curve again near the joint. **\n");
    else
        fprintf("  No second close approach found in the sampled neighborhood -- no evidence of a real self-intersection here.\n");
    end
end
