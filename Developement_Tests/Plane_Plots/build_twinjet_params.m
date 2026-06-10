function gparams = build_twinjet_params(varargin)
%BUILD_PLANE_PARAMS  Plane geometry parameters with defaults + overrides.
%
% Usage:
%   gparams = build_plane_params();
%   gparams = build_plane_params('a', 0.45, 'q', 3.0);
%
% Unknown parameter names trigger an error.

gparams = struct();

% defaults
gparams.a = 1;      % fuselage half-length
gparams.b = 0.08;   % fuselage half-thickness
gparams.q = 10.0;   % nose sharpness (0 = ellipse)
gparams.wing_frac = 0.32;       % fraction where wing starts
gparams.wing_width = 0.05;      % wing width as fraction in param space
gparams.wing_angle = pi/3;      % wing angle
gparams.wing_length = 0.7;      % wing_length
gparams.wing_shape1 = 0.8;      % wing shape param 1, in (0, 1)
gparams.wing_shape2 = 0.9;      % wing shape param 2, in (0, 1)
gparams.wing_lead = 0.07;       % wing leading edge thickness
gparams.closed_engine = false;  % back of the engine closed flag
gparams.engine_loc = 0.2;       % engine location as fraction of wl
gparams.engine_length = 0.2;    % engine length (approx)
gparams.engine_width = 0.01;    % engine width (approx)
gparams.engine_shape1 = 0.07;   % engine curvature
gparams.engine_shape2 = 0.05;   % how closed engine is at the back
gparams.engine_shape3 = 0.3;    % how closed engine is at the front
gparams.closed_front = true;    % front is closed or not
gparams.front_angle = 0.4;      % front opening in radians

% validate number of parameters
if mod(nargin, 2) ~= 0
    error('Arguments must be name-value pairs.');
end

% extract preset
preset = '';
for i = 1:2:nargin
    if strcmp(varargin{i}, 'preset')
        preset = varargin{i+1};
    end
end

% apply preset params
switch lower(preset)

    case 'plane1'
        % same as defaults -> do nothing :)

    case 'plane2'
        % same as plane1 but with open front
        gparams.closed_front = false;
        gparams.front_angle = 0.4;

        gparams.engine_shape1 = 0.028;  % engine curvature
        gparams.engine_shape3 = 0.28;   % wider
    case 'plane3'
        % same as plane1 but with engine more open at the back
          gparams.engine_shape2 = 0.11; % was 0.2
        gparams.engine_width = 0.029;
    case 'plane4'
        % narrower engine with narrower openings and open front
        gparams.engine_shape1 = 0.06;   % narrower engine curvature
        gparams.engine_shape3 = 0.4;    % narrower front opening
        gparams.closed_front = false;
        gparams.front_angle = 0.48;
        gparams.q = 5;
    case 'plane5'
        % wider frame, open front, open narrow engines
        gparams.b = 0.10;
        gparams.wing_width = 0.1;
        gparams.closed_front = false;
        gparams.front_angle = 0.5;
        gparams.engine_loc = 0.2;
        gparams.engine_shape1 = 0.077;   % narrower engine curvature
        gparams.engine_shape2 = 0.13;   % wider back opening
        gparams.engine_shape3 = 0.18;    % narrower front opening
    case ''
        % no preset -> do nothing :)

    otherwise
        error("Unknown preset '%s'.", preset);
end

% override defaults with user-defined parameters
for i = 1:2:nargin
    % extract name-value pairs
    name  = varargin{i};
    value = varargin{i+1};
    % validate parameter is str
    if ~ischar(name) && ~isstring(name)
        error('Parameter names must be strings.');
    end
    name = char(name);
    % skip preset since it was already handled
    if strcmp(name, 'preset')
        continue;
    end
    % validate parameter in acceptable list
    if ~isfield(gparams, name)
        error("Unknown parameter '%s'.", name);
    end
    % set parameter
    gparams.(name) = value;
end

end