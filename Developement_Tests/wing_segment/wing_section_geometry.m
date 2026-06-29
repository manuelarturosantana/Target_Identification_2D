function geom = wing_section_geometry(gparams)

[pieces, data] = build_left_nacelle_wing_tip(gparams);

geom = struct();
geom.is_closed = false;
geom.pieces = pieces;
geom.data = data;

end

function [pieces, data] = build_left_nacelle_wing_tip(gparams)

if gparams.closed_front
    front = build_front(gparams);
    xfl = front.x(front.tf);
    yfl = front.y(front.tf);
else
    [~, front2] = build_front_split(gparams);
    xfl = front2.x(front2.tf);
    yfl = front2.y(front2.tf);
end

back = build_back(gparams);
xbl = back.x(back.t0);
ybl = back.y(back.t0);

[a, b, x0, y0, tstar] = get_wing_params(gparams, xfl, yfl);

lwing_p1 = build_left_wing_lead_1(gparams, a, b, x0, y0, tstar);

xtip = lwing_p1.x(gparams.wing_shape2 * pi);
ytip = lwing_p1.y(gparams.wing_shape2 * pi);

yeng  = lwing_p1.y(lwing_p1.tf);
xeng1 = xbl + (xtip - xbl) * (yeng - ybl) / (ytip - ybl);
xeng2 = lwing_p1.x(lwing_p1.tf);

[ae, be, x0e, y0e, t0e, tfe] = ...
    get_engine_params(gparams, xeng1, xeng2, yeng);

phi = gparams.wing_angle;
cl = cos(phi);
sl = sin(phi);

[t0e_top, tfe_top, t0t, t0l] = compute_engine_cuts( ...
    x0e, y0e, ae, be, gparams.engine_shape1, ...
    x0, y0, a, b, cl, sl, ...
    xbl, ybl, xeng1, yeng);

sep_lead = 0.0;
sep_trail = 0.0;

if isfield(gparams, "separation_lead")
    sep_lead = gparams.separation_lead;
end

if isfield(gparams, "separation_trail")
    sep_trail = gparams.separation_trail;
end

% Interpret separations as fractions of the remaining piece parameter width.
sep_lead = max(0.0, min(1.0, sep_lead));
sep_trail = max(0.0, min(1.0, sep_trail));

% Move the leading-edge start away from the engine, toward the wing tip.
t0l = t0l + sep_lead * (gparams.wing_shape2*pi - t0l);

% Move the trailing-edge start away from the engine, toward the frame.
t0t = t0t - sep_trail * (1 - t0t);

lwing_p2 = build_left_wing_lead_2(gparams, a, b, x0, y0, t0l);

x0t = xbl + t0t * (xeng1 - xbl);
y0t = ybl + t0t * (yeng - ybl);
lwing_p4 = straight_line(x0t, y0t, xtip, ytip);

leng_p1 = build_left_engine_p1(gparams, ae, be, x0e, y0e, t0e);
leng_p2 = build_left_engine_p2(gparams, ae, be, x0e, y0e, t0e, tfe);
leng_p3 = build_left_engine_p3(gparams, ae, be, x0e, y0e, tfe);
leng_p4 = build_left_engine_p4(gparams, ae, be, x0e, y0e, tfe_top);
leng_p5 = build_left_engine_p5(gparams, ae, be, x0e, y0e, t0e_top, ...
    tfe_top);
leng_p6 = build_left_engine_p6(gparams, ae, be, x0e, y0e, t0e_top);

% Closed wing/nacelle substructure.
lwing_p2 = mark_closed(lwing_p2);
lwing_p4 = mark_closed(lwing_p4);
leng_p5 = mark_closed(leng_p5);

% Open nacelle dangling pieces.
leng_p1 = mark_open(leng_p1);
leng_p3 = mark_open(leng_p3);
leng_p4 = mark_open(leng_p4);
leng_p6 = mark_open(leng_p6);
leng_p2 = mark_open(leng_p2);

% Change some normals orientations that need flipping.
lwing_p4 = flip_normals(lwing_p4);
leng_p5 = flip_normals(leng_p5);

pieces = {lwing_p2;
          lwing_p4;
          leng_p1;
          leng_p2;
          leng_p3;
          leng_p4;
          leng_p5;
          leng_p6};

pieces = pieces(~cellfun('isempty', pieces));

data = struct();
data.ae = ae;
data.be = be;
data.x0e = x0e;
data.y0e = y0e;
data.t0e = t0e;
data.tfe = tfe;
data.t0e_top = t0e_top;
data.tfe_top = tfe_top;
data.t0t = t0t;
data.t0l = t0l;
data.separation_lead = sep_lead;
data.separation_trail = sep_trail;
data.xeng1 = xeng1;
data.xeng2 = xeng2;
data.yeng = yeng;
data.xtip = xtip;
data.ytip = ytip;

end

function piece = flip_normals(piece)
%FLIP_NORMALS Reverse the orientation of stored normal vectors.

if isempty(piece)
    return
end

if ~isfield(piece, "nx") || ~isfield(piece, "ny")
    error("Piece does not have nx/ny fields.");
end

nx_old = piece.nx;
ny_old = piece.ny;

piece.nx = @(t) -nx_old(t);
piece.ny = @(t) -ny_old(t);

end

function piece = mark_closed(piece)

if isempty(piece)
    return
end

piece.is_closed = true;

piece.nx = @(t) piece.yp(t) ./ sqrt(piece.xp(t).^2 + piece.yp(t).^2);
piece.ny = @(t) -piece.xp(t) ./ sqrt(piece.xp(t).^2 + piece.yp(t).^2);

end

function piece = mark_open(piece)

if isempty(piece)
    return
end

piece.is_closed = false;

end

% FRAME
function piece = build_front(params)
% builds front of the plane until wing junction

% unpack params
a  = params.a;           % fuselage ellipse major axis (x-direction)
b  = params.b;           % fuselage ellipse minor axis (y-direction)
q  = params.q;           % deformation parameter to make point more pointy
wf = params.wing_frac;   % wing starts at wf fraction of parameter space

% build piece
piece = struct();
piece.t0 = -wf * pi;
piece.tf =  wf * pi;
piece.edge0 = true;
piece.edge1 = true;

% point deformation auxiliary functions
w  = @(t) (1 + cos(t))./2;
dw = @(t) -sin(t)./2;
m  = @(t) 1 - (q/(1+q)) * (1 - w(t)).^2;
dm = @(t) 2*(q/(1+q))*(1-w(t)) .* dw(t);

% piece parametrization
piece.x  = @(t) a * m(t) .* cos(t);
piece.y  = @(t) b * sin(t);
piece.xp = @(t) a * (dm(t) .* cos(t) - m(t) .* sin(t));
piece.yp = @(t) b * cos(t);

end

function [piece1,  piece2] = build_front_split(params)
% build front of the plane with gap in the middle

% unpack params
a  = params.a;           % fuselage ellipse major axis (x-direction)
b  = params.b;           % fuselage ellipse minor axis (y-direction)
q  = params.q;           % deformation parameter to make point more pointy
wf = params.wing_frac;   % wing starts at wf fraction of parameter space
fa = params.front_angle; % front opening angle

% build bottom piece
piece1 = struct();
piece1.t0 = -wf * pi;
piece1.tf = -fa/2;
piece1.edge0 = true;
piece1.edge1 = true;
% point deformation auxiliary functions
w  = @(t) (1 + cos(t))./2;
dw = @(t) -sin(t)./2;
m  = @(t) 1 - (q/(1+q)) * (1 - w(t)).^2;
dm = @(t) 2*(q/(1+q))*(1-w(t)) .* dw(t);
% piece parametrization
piece1.x  = @(t) a * m(t) .* cos(t);
piece1.y  = @(t) b * sin(t);
piece1.xp = @(t) a * (dm(t) .* cos(t) - m(t) .* sin(t));
piece1.yp = @(t) b * cos(t);

% build top piece
piece2 = struct();
piece2.t0 = fa/2;
piece2.tf = wf * pi;
piece2.edge0 = true;
piece2.edge1 = true;
% piece parametrization
piece2.x  = @(t) a * m(t) .* cos(t);
piece2.y  = @(t) b * sin(t);
piece2.xp = @(t) a * (dm(t) .* cos(t) - m(t) .* sin(t));
piece2.yp = @(t) b * cos(t);

end

function piece = build_back(params)
% builds back of the plane

% unpack params
a  = params.a;           % fuselage ellipse major axis (x-direction)
b  = params.b;           % fuselage ellipse minor axis (y-direction)
q  = params.q;           % deformation parameter to make point more pointy
wf = params.wing_frac;   % wing starts at wf fraction of parameter space
wd = params.wing_width;  % wing width in parameter space

% build piece
piece = struct();
piece.t0 = (wf + wd) * pi;
piece.tf = 2*pi - (wf + wd) * pi;
piece.edge0 = true;
piece.edge1 = true;

% point deformation auxiliary functions
w  = @(t) (1 + cos(t))./2;
dw = @(t) -sin(t)./2;
m  = @(t) 1 - (q/(1+q)) * (1 - w(t)).^2;
dm = @(t) 2*(q/(1+q))*(1-w(t)) .* dw(t);

% piece parametrization
piece.x  = @(t) a * m(t) .* cos(t);
piece.y  = @(t) b * sin(t);
piece.xp = @(t) a * (dm(t) .* cos(t) - m(t) .* sin(t));
piece.yp = @(t) b * cos(t);

end

% LEFT WING
function [a, b, x0, y0, tstar] = get_wing_params(params, xstar, ystar)
% get leading edge (ellipse section) parameters from global plane parameters

% unpack params
s1  = params.wing_shape1;
b   = params.wing_lead;
q   = params.q;
wl  = params.wing_length;
phi = params.wing_angle;
af  = params.a;
bf  = params.b;

% rotation params
c = cos(phi);
s = sin(phi);

% major axis
a = wl/(s1 + params.wing_shape2);

% find solution to nonlinear problem for intersection
f = @(u) wing_residual(u, a, b, c, s, xstar, ystar, q, af, bf, s1, phi);
opts = optimoptions("fsolve","Display","none");
u0 = [0.0; xstar - a*s1*cos(phi); ystar - a*s1*sin(phi); 0.2];
u = fsolve(f, u0, opts);

v = u(1);
tstar = pi / (1 + exp(-v));
x0 = u(2);
y0 = u(3);
% ttilde = u(4);    % unused outside residual

end

function r = wing_residual(u, a, b, c, s, xstar, ystar, ...
    q, af, bf, s1, phi)

v      = u(1);
x0     = u(2);
y0     = u(3);
ttilde = u(4);

% force tstar to lie in [0, pi]
tstar = pi / (1 + exp(-v));

% point deformation auxiliary functions
w  = @(t) (1 + cos(t))./2;
m  = @(t) 1 - (q/(1+q)) * (1 - w(t)).^2;
xf = @(t) af * m(t) .* cos(t);
yf = @(t) bf * sin(t);

p1 = a*cos(tstar) - (c*(xstar-x0) - s*(ystar-y0));
p2 = b*sin(tstar) - (s*(xstar-x0) + c*(ystar-y0));
p3 = x0 - a*s1*cos(phi) - xf(ttilde);
p4 = y0 - a*s1*sin(phi) - yf(ttilde);

r = [p1; p2; p3; p4];
end

function piece = build_left_wing_lead_1(params, a, b, x0, y0, tstar1)
% tstar1: intersection with the main frame
% tstar2: intersection with engine

% unpack params
phi  = params.wing_angle;
s2   = params.wing_shape2;
eloc = params.engine_loc;

% build piece
piece = struct();
piece.t0 = tstar1;
piece.tf = tstar1 + eloc * (s2*pi - tstar1);
piece.edge0 = 1;
piece.edge1 = 1;

% original ellipse parametrization
xe  = @(t) a * cos(t);
ye  = @(t) b * sin(t);
xep = @(t) -a * sin(t);
yep = @(t)  b * cos(t);

% rotation in wing angle
cph = cos(phi);
sph = sin(phi);

% geometry
piece.x  = @(t) (x0 + cph*xe(t) + sph*ye(t));
piece.y  = @(t) (y0 - sph*xe(t) + cph*ye(t));
piece.xp = @(t) cph*xep(t) + sph*yep(t);
piece.yp = @(t) -sph*xep(t) + cph*yep(t);

end

function piece = build_left_wing_lead_2(params, a, b, x0, y0, tstar2)
% leading edge from engine intersection to wing tip

phi = params.wing_angle;
s2  = params.wing_shape2;

piece = struct();
piece.t0 = tstar2;
piece.tf = s2 * pi;
piece.edge0 = 1;
piece.edge1 = 1;

% base ellipse
xe  = @(t) a*cos(t);
ye  = @(t) b*sin(t);
xep = @(t) -a*sin(t);
yep = @(t)  b*cos(t);

% rotation
cph = cos(phi);
sph = sin(phi);

piece.x  = @(t) x0 + cph*xe(t) + sph*ye(t);
piece.y  = @(t) y0 - sph*xe(t) + cph*ye(t);

piece.xp = @(t) cph*xep(t) + sph*yep(t);
piece.yp = @(t) -sph*xep(t) + cph*yep(t);
end

function piece = straight_line(x0, y0, xf, yf)
% straight line parametrization

% build piece
piece = struct();
piece.t0 = 0;
piece.tf = 1;
piece.edge0 = 1;
piece.edge1 = 1;

% straight line parametrization
piece.x  = @(t) x0 + t*(xf-x0);
piece.y  = @(t) y0 + t*(yf-y0);
piece.xp = @(t) (xf-x0)*ones(size(t));
piece.yp = @(t) (yf-y0)*ones(size(t));

end

% LEFT ENGINE
function r = engine_residual(u, a, b, c, x1, x2, y)

x0 = u(1);
y0 = u(2);
t0 = u(3);
tf = u(4);

p1 = a*cos(t0) + x0 - x1;
p2 = y0 + b*sin(t0) - c*sin(t0).^2 - y;
p3 = a*cos(tf) + x0 - x2;
p4 = y0 + b*sin(tf) - c*sin(tf).^2 - y;

r = [p1; p2; p3; p4];

end

function [a, b, x0, y0, t0, tf] = get_engine_params(params, x1, x2, y)
% get params for engine ellipse

% unpack params
a = params.engine_length;
b = params.engine_width;
c = params.engine_shape1;

% initial guess: center x at midpoint, center y roughly below y
x0g = 0.5 * (x1 + x2);
y0g = y + b;

% initial guess for angles from x positions (lower-half guess)
arg = (x2 - x1) / (2 * a);
arg = max(-1.0, min(1.0, arg));
theta = acos(arg);

t0g = pi + theta;
tfg = 2*pi - theta;

u0 = [x0g; y0g; t0g; tfg];

% solve
f = @(u) engine_residual(u, a, b, c, x1, x2, y);
opts = optimoptions("fsolve","Display","none", ...
    "FunctionTolerance", 1e-12, "StepTolerance", 1e-12);

u = fsolve(f, u0, opts);

% unpack
x0 = u(1);
y0 = u(2);
t0 = mod(u(3), 2*pi);
tf = mod(u(4), 2*pi);

% enforce that t0 corresponds to (x1,y) and tf to (x2,y)
x_at_t0 = x0 + a*cos(t0);
x_at_tf = x0 + a*cos(tf);

if abs(x_at_t0 - x1) > abs(x_at_tf - x1)
    tmp = t0; t0 = tf; tf = tmp;
end

end

function piece = build_left_engine_p1(params, a, b, x0, y0, t0)
% builds dangling piece of engine between trail intersection and engine
% intersection
%
% returns empty if the arc does not exist

% unpack params
s2 = params.engine_shape2;
c  = params.engine_shape1;

% trail-edge intersection parameter
tcut = pi + s2 * pi;

% check if piece exists
if t0 <= tcut
    piece = [];
    return
end

% build piece
piece = struct();
piece.t0 = tcut;
piece.tf = t0;
piece.edge0 = 1;
piece.edge1 = 1;

% engine parametrization
piece.x  = @(t) x0 + a*cos(t);
piece.y  = @(t) y0 + b*sin(t) - c*sin(t).^2;

% derivatives
piece.xp = @(t) -a*sin(t);
piece.yp = @(t) b*cos(t) - 2*c*cos(t).*sin(t);

end

function piece = build_left_engine_p2(params, a, b, x0, y0, t0, tf)
% builds inside-wing piece of engine between trail and lead intersections

% unpack params
c  = params.engine_shape1;

% build piece
piece = struct();
piece.t0 = t0;
piece.tf = tf;
piece.edge0 = 1;
piece.edge1 = 1;

% engine parametrization
piece.x  = @(t) x0 + a*cos(t);
piece.y  = @(t) y0 + b*sin(t) - c*sin(t).^2;

% derivatives
piece.xp = @(t) -a*sin(t);
piece.yp = @(t) b*cos(t) - 2*c*cos(t).*sin(t);

end

function piece = build_left_engine_p3(params, a, b, x0, y0, tf)
% builds dangling piece of engine between lead intersection and front opening
%
% returns empty if the arc does not exist

% unpack params
s3 = params.engine_shape3;   % opening at the front
c  = params.engine_shape1;

% front opening parameter
tcut = 2*pi - s3 * pi;

% check if piece exists
if tf >= tcut
    piece = [];
    return
end

% build piece
piece = struct();
piece.t0 = tf;
piece.tf = tcut;
piece.edge0 = 1;
piece.edge1 = 1;

% engine parametrization
piece.x  = @(t) x0 + a*cos(t);
piece.y  = @(t) y0 + b*sin(t) - c*sin(t).^2;

% derivatives
piece.xp = @(t) -a*sin(t);
piece.yp = @(t) b*cos(t) - 2*c*cos(t).*sin(t);

end

function piece = build_left_engine_p4(params, a, b, x0, y0, t_tr)
% top dangling piece near the back (analogue of p1), using the TOP arc
% and the TOP trail-intersection parameter t_tr.

s2 = params.engine_shape2;
c  = params.engine_shape1;

% top trailing cut is the mirror of (pi + s2*pi): 2*pi - (pi + s2*pi)
tcut = pi - s2 * pi;   % in (0, pi)

if t_tr > tcut
    piece = [];
    return
end

piece = struct();
piece.t0 = t_tr;
piece.tf = tcut;
piece.edge0 = true;
piece.edge1 = true;

piece.x  = @(t) x0 + a*cos(t);
piece.y  = @(t) y0 + b*sin(t) + c*sin(t).^2;

piece.xp = @(t) -a*sin(t);
piece.yp = @(t) b*cos(t) + 2*c*cos(t).*sin(t);
end

function piece = build_left_engine_p5(params, a, b, x0, y0, t_le, t_tr)
% top inside-wing piece between lead and trail intersections

c = params.engine_shape1;

% order on [0,pi]
t0 = t_le; tf = t_tr;
if tf < t0
    tmp = t0; t0 = tf; tf = tmp;
end

piece = struct();
piece.t0 = t0;
piece.tf = tf;
piece.edge0 = true;
piece.edge1 = true;

piece.x  = @(t) x0 + a*cos(t);
piece.y  = @(t) y0 + b*sin(t) + c*sin(t).^2;

piece.xp = @(t) -a*sin(t);
piece.yp = @(t) b*cos(t) + 2*c*cos(t).*sin(t);
end

function piece = build_left_engine_p6(params, a, b, x0, y0, t_le)
% top front dangling piece from opening cut to lead intersection

s3 = params.engine_shape3;
c  = params.engine_shape1;

tcut = s3 * pi;   % front opening on top arc

if t_le <= tcut
    piece = [];
    return
end

piece = struct();
piece.t0 = tcut;
piece.tf = t_le;
piece.edge0 = true;
piece.edge1 = true;

piece.x  = @(t) x0 + a*cos(t);
piece.y  = @(t) y0 + b*sin(t) + c*sin(t).^2;

piece.xp = @(t) -a*sin(t);
piece.yp = @(t) b*cos(t) + 2*c*cos(t).*sin(t);
end

function [t0e, tfe, t0t, t0l] = compute_engine_cuts( ...
    x0e, y0e, ae, be, ce, ...
    x0l, y0l, al, bl, cl, sl, ...
    xbl, ybl, xeng1, yeng)
% Solve for intersection of the TOP engine arc with:
%   (1) the wing leading ellipse (parameter t0l)
%   (2) the wing trailing segment from (xbl,ybl) to (xeng1,yeng) 
%       (parameter t0t)
%
% Unknowns u = [t0e; tfe; t0t; t0l]
%   t0e : engine parameter at lead intersection (top arc)
%   tfe : engine parameter at trail intersection (top arc)
%   t0t : trailing segment parameter in [0,1]
%   t0l : wing-lead ellipse parameter

% initial guess for t0e
t0e_g = 0.5*pi;
tfe_g = 0.5*pi;
% initial guess for wing lead parameter t0l
t0l_g = pi - 0.1;
% initial guess for trailing segment parameter t0t
t0t_g = 1.1;

u0 = [t0e_g; tfe_g; t0t_g; t0l_g];

f = @(u) engine_cut_res(u, x0e, y0e, ae, be, ce, ...
    x0l, y0l, al, bl, cl, sl, xbl, ybl, xeng1, yeng);

opts = optimoptions("fsolve", "Display", "none", ...
    "FunctionTolerance", 1e-12, "StepTolerance", 1e-12);

u = fsolve(f, u0, opts);

t0e = mod(u(1), 2*pi);
tfe = mod(u(2), 2*pi);
t0t = u(3);
t0l = mod(u(4), 2*pi);

end

function r = engine_cut_res(u, x0e, y0e, ae, be, ce, x0l, y0l, al, bl, ...
    cl, sl, xbl, ybl, xeng1, yeng)

t0e = u(1);     % t0 engine (coincides with t0 lead)
tfe = u(2);     % tf engine (coincides with t0 trail)
t0t = u(3);     % t0 trail (coincides with tf engine)
t0l = u(4);     % t0 lead (coincides with t0 engine)

xt = xbl + t0t*(xeng1-xbl);
yt = ybl + t0t*(yeng-ybl);

% equations
% p1: x lead with x engine
p1 = x0e + ae*cos(t0e) - (x0l + cl*al*cos(t0l) + sl*bl*sin(t0l));
% p2: y lead with y engine
p2 = (y0e + be*sin(t0e) + ce*sin(t0e).^2 - ...
    (y0l - sl*al*cos(t0l) + cl*bl*sin(t0l)));
% p3: x trail with x engine
p3 = x0e + ae * cos(tfe) - xt;
% p4: y trail with y engine
p4 = y0e + be*sin(tfe) + ce*sin(tfe).^2 - yt;

r = [p1; p2; p3; p4];

end

% RIGHT EVERYTHING (MIRROR IN Y)
function pieceR = mirror_piece_y(pieceL)
% Mirror a parametric piece across the x-axis: (x,y) -> (x,-y).

% guards against empty input piece
if isempty(pieceL)
    pieceR = [];
    return
end

pieceR = pieceL;

pieceR.y  = @(t) -pieceL.y(t);
pieceR.yp = @(t) -pieceL.yp(t);

% x/xp unchanged
pieceR.x  = pieceL.x;
pieceR.xp = pieceL.xp;
end
